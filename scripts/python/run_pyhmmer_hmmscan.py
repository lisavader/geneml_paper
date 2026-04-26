#!/usr/bin/env python3

"""Run `hmmscan` with pyhmmer.

This script is a lightweight replacement for the HMMER `hmmscan` binary.
It keeps the common scoring and reporting options, but writes tabular output
using PyHMMER instead of spawning an external process.
"""

from __future__ import annotations

import argparse
import contextlib
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, cast
import sys

from pyhmmer import easel, hmmer, plan7  # type: ignore[import-not-found]  # pylint: disable=import-error


_WORKER_STATE: dict[str, object | None] = {
    "profiles": None,
    "alphabet": None,
    "pipeline_options": None,
    "cpu": None,
}


@dataclass
class Job:
    seqfile: str
    tblout: str | None
    domtblout: str | None
    pfamtblout: str | None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_pyhmmer_hmmscan.py",
        description="Scan protein sequences against an HMM database using pyhmmer.",
        allow_abbrev=False,
    )
    parser.add_argument("hmmfile", help="HMM database in HMMER format")
    parser.add_argument("seqfile", nargs="?", help="Protein sequence file in FASTA or another supported format (single-run mode)")
    parser.add_argument("--batch-tsv", help="Batch input list file with one sequence path per line")
    parser.add_argument("--output-dir", help="Output directory for batch mode")
    parser.add_argument("--batch-jobs", type=int, default=1, help="Number of batch worker processes (each loads DB once, then runs serially)")

    parser.add_argument("--tblout", help="Write per-target tabular output to this file")
    parser.add_argument("--domtblout", help="Write per-domain tabular output to this file")
    parser.add_argument("--pfamtblout", help="Write Pfam-style tabular output to this file")

    parser.add_argument("--cpu", type=int, default=1, help="Number of worker threads to use")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for the search pipeline")
    parser.add_argument("--E", type=float, default=10.0, help="Per-target reporting E-value threshold")
    parser.add_argument("--T", type=float, default=None, help="Per-target reporting bit score threshold")
    parser.add_argument("--domE", type=float, default=10.0, help="Per-domain reporting E-value threshold")
    parser.add_argument("--domT", type=float, default=None, help="Per-domain reporting bit score threshold")
    parser.add_argument("--incE", type=float, default=0.01, help="Per-target inclusion E-value threshold")
    parser.add_argument("--incT", type=float, default=None, help="Per-target inclusion bit score threshold")
    parser.add_argument("--incdomE", type=float, default=0.01, help="Per-domain inclusion E-value threshold")
    parser.add_argument("--incdomT", type=float, default=None, help="Per-domain inclusion bit score threshold")
    parser.add_argument("--Z", type=float, default=None, help="Effective number of targets searched")
    parser.add_argument("--domZ", type=float, default=None, help="Effective number of significant targets for domains")
    parser.add_argument("--F1", type=float, default=0.02, help="MSV filter threshold")
    parser.add_argument("--F2", type=float, default=0.001, help="Viterbi filter threshold")
    parser.add_argument("--F3", type=float, default=1e-05, help="Forward filter threshold")
    parser.add_argument("--cut_ga", action="store_true", help="Use gathering thresholds from the HMMs")
    parser.add_argument("--cut_nc", action="store_true", help="Use noise thresholds from the HMMs")
    parser.add_argument("--cut_tc", action="store_true", help="Use trusted thresholds from the HMMs")
    parser.add_argument("--max", action="store_true", help="Disable filters and use maximum sensitivity")
    parser.add_argument("--nobias", action="store_true", help="Disable the composition bias filter")
    parser.add_argument("--acc", action="store_true", help="Accepted for CLI compatibility; accession fields are written when available")
    parser.add_argument("--noali", action="store_true", help="Accepted for CLI compatibility; alignment display is not produced by this wrapper")
    parser.add_argument("--notextw", action="store_true", help="Accepted for CLI compatibility")
    parser.add_argument("--stream-profiles", action="store_true", help="Do not preload profiles into memory; stream from the HMM file (lower RAM, slower)")
    return parser


def open_output(path: str | None):
    if path in (None, "-"):
        return contextlib.nullcontext(sys.stdout.buffer)
    assert path is not None
    return open(path, "wb")


def make_pipeline_options(args: argparse.Namespace) -> dict[str, object]:
    options: dict[str, object] = {
        "E": args.E,
        "T": args.T,
        "domE": args.domE,
        "domT": args.domT,
        "incE": args.incE,
        "incT": args.incT,
        "incdomE": args.incdomE,
        "incdomT": args.incdomT,
        "Z": args.Z,
        "domZ": args.domZ,
        "F1": args.F1,
        "F2": args.F2,
        "F3": args.F3,
        "seed": args.seed,
        "bias_filter": not args.nobias,
    }
    if args.max:
        options.update({"bias_filter": False, "F1": 1.0, "F2": 1.0, "F3": 1.0})
    if args.cut_ga:
        options["bit_cutoffs"] = "gathering"
    elif args.cut_nc:
        options["bit_cutoffs"] = "noise"
    elif args.cut_tc:
        options["bit_cutoffs"] = "trusted"
    return options


def parse_batch_tsv(tsv_path: str) -> list[Job]:
    jobs: list[Job] = []
    with open(tsv_path, "r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                raise ValueError(f"Invalid --batch-tsv line {line_number}: expected one input path per line")
            jobs.append(Job(seqfile=line, tblout=None, domtblout=None, pfamtblout=None))
    if not jobs:
        raise ValueError("No jobs found in --batch-tsv")
    return jobs


def resolve_jobs(args: argparse.Namespace) -> list[Job]:
    if args.batch_tsv:
        if args.seqfile is not None:
            raise ValueError("Provide either seqfile or --batch-tsv, not both")
        if args.batch_jobs < 1:
            raise ValueError("--batch-jobs must be >= 1")
        if args.output_dir is None:
            raise ValueError("--output-dir is required when --batch-tsv is provided")
        if args.tblout or args.domtblout or args.pfamtblout:
            raise ValueError("Do not use --tblout/--domtblout/--pfamtblout with --batch-tsv")

        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        jobs = parse_batch_tsv(args.batch_tsv)
        resolved: list[Job] = []
        for job in jobs:
            stem = Path(job.seqfile).stem
            resolved.append(Job(seqfile=job.seqfile, tblout=None, domtblout=str(output_dir / f"{stem}.domtblout.tsv"), pfamtblout=None))
        return resolved

    if args.seqfile is None:
        raise ValueError("seqfile is required when --batch-tsv is not provided")
    if args.batch_jobs != 1:
        raise ValueError("--batch-jobs is only valid with --batch-tsv")
    if args.output_dir is not None:
        raise ValueError("--output-dir is only valid with --batch-tsv")

    return [
        Job(
            seqfile=args.seqfile,
            tblout=args.tblout,
            domtblout=args.domtblout,
            pfamtblout=args.pfamtblout,
        )
    ]


def load_profiles(hmm_path: str, stream_profiles: bool):
    if stream_profiles:
        hmm_file = plan7.HMMFile(hmm_path)
        return hmm_file, (hmm_file.optimized_profiles() if hmm_file.is_pressed() else hmm_file)

    with plan7.HMMFile(hmm_path) as hmm_file:
        profiles = list(hmm_file.optimized_profiles()) if hmm_file.is_pressed() else list(hmm_file)
    return contextlib.nullcontext(None), profiles


def run_one_job(
    job: Job,
    profiles: Iterable,
    alphabet: easel.Alphabet,
    pipeline_options: dict[str, object],
    cpu: int,
) -> None:
    output_specs = [
        (job.tblout, "targets"),
        (job.domtblout, "domains"),
        (job.pfamtblout, "pfam"),
    ]
    outputs: list[tuple[str | None, str]] = [(path, fmt) for path, fmt in output_specs if path is not None]
    if not outputs:
        outputs = [(None, "targets")]

    with contextlib.ExitStack() as stack:
        sequences = stack.enter_context(easel.SequenceFile(job.seqfile, digital=True, alphabet=alphabet))
        opened_outputs = []
        for path, fmt in outputs:
            handle = stack.enter_context(open_output(path))
            opened_outputs.append((handle, fmt))

        wrote_header: dict[int, bool] = {index: False for index in range(len(opened_outputs))}

        for hits in hmmer.hmmscan(
            sequences,
            profiles,
            cpus=cpu,
            backend="threading",
            **pipeline_options,
        ):
            for index, (handle, fmt) in enumerate(opened_outputs):
                hits.write(handle, format=fmt, header=not wrote_header[index])
                wrote_header[index] = True


def split_jobs(jobs: list[Job], parts: int) -> list[list[Job]]:
    groups: list[list[Job]] = [[] for _ in range(parts)]
    for index, job in enumerate(jobs):
        groups[index % parts].append(job)
    return [group for group in groups if group]


def init_batch_worker(hmmfile: str, pipeline_options: dict[str, object], cpu: int) -> None:
    _WORKER_STATE["alphabet"] = easel.Alphabet.amino()
    _WORKER_STATE["pipeline_options"] = pipeline_options
    _WORKER_STATE["cpu"] = cpu
    _, _WORKER_STATE["profiles"] = load_profiles(hmmfile, stream_profiles=False)


def run_batch_chunk(chunk: list[Job]) -> None:
    profiles = _WORKER_STATE["profiles"]
    alphabet = _WORKER_STATE["alphabet"]
    pipeline_options = _WORKER_STATE["pipeline_options"]
    cpu = _WORKER_STATE["cpu"]

    if profiles is None or alphabet is None or pipeline_options is None or cpu is None:
        raise RuntimeError("Batch worker was not initialized")

    profiles = cast(Iterable, profiles)
    alphabet = cast(easel.Alphabet, alphabet)
    pipeline_options = cast(dict[str, object], pipeline_options)
    cpu = cast(int, cpu)

    for job in chunk:
        run_one_job(job, profiles, alphabet, pipeline_options, cpu)


def run_hmmscan(args: argparse.Namespace) -> int:
    alphabet = easel.Alphabet.amino()
    pipeline_options = make_pipeline_options(args)
    jobs = resolve_jobs(args)

    if args.batch_tsv and args.batch_jobs > 1:
        if args.stream_profiles:
            raise ValueError("--stream-profiles is not supported with --batch-tsv and --batch-jobs > 1")
        chunks = split_jobs(jobs, args.batch_jobs)
        with ProcessPoolExecutor(
            max_workers=len(chunks),
            initializer=init_batch_worker,
            initargs=(args.hmmfile, pipeline_options, args.cpu),
        ) as pool:
            for future in [pool.submit(run_batch_chunk, chunk) for chunk in chunks]:
                future.result()
        return 0

    profile_ctx, profiles = load_profiles(args.hmmfile, stream_profiles=args.stream_profiles)
    with profile_ctx:
        for job in jobs:
            run_one_job(job, profiles, alphabet, pipeline_options, args.cpu)

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_hmmscan(args)


if __name__ == "__main__":
    raise SystemExit(main())
