import argparse
from typing import Set


def get_species_ids_from_species(clade_set: Set[str], orthodb_species: str
                         ) -> tuple[Set[str], Set[str]]:
    """
    Find species ids for given clade names from an OrthoDB species file.

    Args:
        clade_set (Set[str]): Set of clade names (case-insensitive).
        orthodb_species (str): Path to OrthoDB species file (.tsv).

    Returns:
        Tuple[Set[str], Set[str]]: Set of matching species ids and set of found clade names.

    Notes:
        - Assumes file is tab-separated, species_id in column 2, clade_name in column 3.
    """
    if not clade_set:
        return set(), set()

    species_ids: Set[str] = set()
    found: Set[str] = set()

    with open(orthodb_species, "r", encoding="utf-8") as stream:
        for line_num, line in enumerate(stream, 1):
            parts = line.strip().split("\t")
            try:
                species_id = parts[1]
                species_name = parts[2].lower()
                genus_name = species_name.split()[0]
            except IndexError as e:
                raise ValueError(f"Malformed line {line_num} in species file: "
                                 f"{line.strip()}") from e

            # Check if either full species or genus matches
            if species_name in clade_set:
                species_ids.add(species_id)
                found.add(species_name)
            if genus_name in clade_set:
                species_ids.add(species_id)
                found.add(genus_name)

    return species_ids, found


def get_species_ids_from_levels(clade_set: Set[str], orthodb_levels: str,
                                orthodb_level2species: str) -> tuple[Set[str], Set[str]]:
    """
    Find species ids for given clade names from an OrthoDB levels file.

    Args:
        clade_set (Set[str]): Set of clade names (case-insensitive).
        orthodb_levels (str): Path to OrthoDB levels file (.tsv).
        orthodb_level2species (str): Path to OrthoDB level to species mapping file (.tsv).

    Returns:
        Tuple[Set[str], Set[str]]: Set of matching species ids and set of found clade names.

    Notes:
        - Assumes file is tab-separated, clade_name in column 2, taxid in column 1.
    """
    if not clade_set:
        return set(), set()

    taxids: Set[str] = set()
    species_ids: Set[str] = set()
    found: Set[str] = set()

    with open(orthodb_levels, "r", encoding="utf-8") as stream:
        for line_num, line in enumerate(stream, 1):
            parts = line.strip().split("\t")
            try:
                clade_name = parts[1]
                taxid = parts[0]
            except IndexError as e:
                raise ValueError(f"Malformed line {line_num} in levels file: "
                                 f"{line.strip()}") from e
            if clade_name.lower() in clade_set:
                taxids.add(taxid)
                found.add(clade_name.lower())

    with open(orthodb_level2species, "r", encoding="utf-8") as stream:
        for line_num, line in enumerate(stream, 1):
            parts = line.strip().split("\t")
            try:
                species_id = parts[1]
                lineage = parts[3].strip("{}").split(",")
            except IndexError as e:
                raise ValueError(f"Malformed line {line_num} in level2species file: "
                                 f"{line.strip()}") from e
            if any(taxid in taxids for taxid in lineage):
                species_ids.add(species_id)

    return species_ids, found


def get_species_ids(clades: str, orthodb_species: str, orthodb_levels: str,
                    orthodb_level2species: str) -> Set[str]:
    """
    Get all species ids for a comma-separated list of clade names, searching both species and levels files.

    Args:
        clades (str): Comma-separated clade names (case-insensitive).
        orthodb_species (str): Path to OrthoDB species file (.tsv).
        orthodb_levels (str): Path to OrthoDB levels file (.tsv).
        orthodb_level2species (str): Path to OrthoDB level to species mapping file (.tsv).

    Returns:
        Set[str]: Set of matching species ids.
    """
    if not clades:
        return set()

    species_ids: Set[str] = set()
    clade_set = set(name.lower() for name in clades.split(","))

    species_ids, found = get_species_ids_from_species(clade_set, orthodb_species)
    not_found = clade_set - found

    if not_found:
        new_species_ids, new_found = get_species_ids_from_levels(not_found, orthodb_levels,
                                                                 orthodb_level2species)
        species_ids.update(new_species_ids)
        found.update(new_found)
        not_found = clade_set - found

    for missing in not_found:
        print(f"Warning: clade '{missing}' not found in species or levels file.")
    return species_ids


def write_selected_entries(selected_species_ids: Set[str], orthodb_aa: str, output: str) -> None:
    """
    Write FASTA entries for selected species ids to output file.

    Args:
        selected_species_ids (Set[str]): Set of species ids to include.
        orthodb_aa (str): Path to OrthoDB amino acid FASTA file.
        output (str): Path to output filtered FASTA file.

    Returns:
        None

    Notes:
        - Assumes species id is the second word in the header line (after '>').
    """
    with open(orthodb_aa, "r", encoding="utf-8") as stream_in, \
         open(output, "w", encoding="utf-8") as stream_out:
        write_entry = False
        for line_num, line in enumerate(stream_in, 1):
            if line.startswith(">"):
                try:
                    species_id = line.split()[1]
                except IndexError as e:
                    raise ValueError(f"Malformed FASTA header at line {line_num}: "
                                     f"{line.strip()}") from e
                if species_id in selected_species_ids:
                    write_entry = True
                    stream_out.write(line)
                else:
                    write_entry = False
            else:
                if write_entry:
                    stream_out.write(line)


def main(orthodb_aa: str, orthodb_species: str, orthodb_levels: str, orthodb_level2species: str,
         output: str, include: str, exclude: str) -> None:
    """
    Main function to filter OrthoDB FASTA by clade and species ids.

    Args:
        orthodb_aa (str): Path to OrthoDB amino acid FASTA file.
        orthodb_levels (str): Path to OrthoDB levels file (.tsv).
        orthodb_level2species (str): Path to OrthoDB level to species mapping file (.tsv).
        orthodb_species (str): Path to OrthoDB species file (.tsv).
        output (str): Path to output filtered FASTA file.
        include (str): Comma-separated clade names to include.
        exclude (str): Comma-separated clade names to exclude.

    Returns:
        None
    """
    include_species_ids = get_species_ids(include, orthodb_species,
                                          orthodb_levels, orthodb_level2species)
    exclude_species_ids = get_species_ids(exclude, orthodb_species,
                                          orthodb_levels, orthodb_level2species)
    selected_species_ids = include_species_ids - exclude_species_ids
    write_selected_entries(selected_species_ids, orthodb_aa, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("orthodb_aa", type=str,
                        help="Path to Orthodb amino acid file (.fasta)")
    parser.add_argument("orthodb_species", type=str,
                        help="Path to Orthodb species file (.tsv)")
    parser.add_argument("orthodb_levels", type=str,
                        help="Path to Orthodb levels file (.tsv)")
    parser.add_argument("orthodb_level2species", type=str,
                        help="Path to Orthodb level to species mapping file (.tsv)")
    parser.add_argument("-o", "--output", type=str,
                        help="Output filtered Orthodb amino acid file (.fasta)", required=True)
    parser.add_argument("--include", type=str,
                        help="Comma-separated list of clade names to include", required=True)
    parser.add_argument("--exclude", type=str,
                        help="Comma-separated list of clade names to exclude", default=None)
    args = parser.parse_args()
    main(**vars(args))
