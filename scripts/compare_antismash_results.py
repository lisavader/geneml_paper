import argparse

from antismash.common.secmet.features.cds_feature import CDSFeature
from antismash.common.secmet.features.pfam_domain import PFAMDomain
from antismash.common.secmet.features.protocluster import Protocluster
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.common.serialiser import AntismashResults

def get_overlapping(features1, features2):
    overlapping = []
    i = 0
    j = 0
    while i < len(features1) and j < len(features2):
        overlaps = features1[i].overlaps_with(features2[j])
        if overlaps:
            if isinstance(features1[i], Protocluster) and features1[i].product == features2[j].product:
                overlapping.extend([features1[i], features2[j]])
            else:
                overlapping.extend([features1[i], features2[j]])
        if features1[i].end < features2[j].end:
            i += 1
        else:
            j += 1
    return overlapping

def build_lines(file_names, features, record_name):
    lines = []
    overlapping = get_overlapping(*features)
    for file_name, feature_list in zip(file_names, features):
        for feature in feature_list:
            unique = feature not in overlapping
            if isinstance(feature, PFAMDomain):
                feature_name = "pfam"
                description = feature.identifier
            elif isinstance(feature, Protocluster):
                feature_name = "protocluster"
                description = feature.product
            elif isinstance(feature, CDSFeature):
                feature_name = str(feature.gene_function)
                description = ','.join([annotation.description for annotation in feature.gene_functions])
            line = '\t'.join([file_name, record_name, feature_name, str(feature.location.start), str(feature.location.end), description, str(unique)])
            lines.append(line)
    return '\n'.join(lines)+'\n'

def main(file1, file2, table_out):
    file_names = (file1, file2)
    as_results = tuple(AntismashResults.from_file(f) for f in file_names)
    with open(table_out, 'w') as s:
        for record_pair in zip(*(r.records for r in as_results)):
            record_names = set(r.name for r in record_pair)
            assert len(record_names) == 1, f"Record names don't match: {record_names}"
            record_name = record_names.pop()
            biosynth_cdses = []
            clusters = []
            pfams = []
            for record in record_pair:
                cdses = tuple(cds for cds in record.get_cds_features() if cds.gene_function != GeneFunction.OTHER)
                biosynth_cdses.append(cdses)
                clusters.append(record.get_protoclusters())
                pfams.append(record.get_pfam_domains())
            for results in biosynth_cdses, clusters, pfams:
                s.write(build_lines(file_names, results, record_name))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("table_out")
    args = parser.parse_args()
    main(**vars(args))
