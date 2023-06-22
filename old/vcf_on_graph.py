'Matches nodes of GFA to variants of VCF.'
from argparse import ArgumentParser, SUPPRESS
from Bio import Align
from seaborn import heatmap
from matplotlib import pyplot as plt


def vcf_parser(vcf_file: str) -> dict:
    """Extracts sequences and identifiers from a VCF file

    Args:
        vcf_file (str): path to a file

    Returns:
        dict: mapping id to sequence
    """
    vcf_datas: dict = {}
    with open(vcf_file, 'r', encoding='utf-8') as vcf_reader:
        for line in vcf_reader:
            if line[0] != '#':
                datas = line.split()
                vcf_datas[f"{datas[2]}_ref"] = datas[3]
                vcf_datas[f"{datas[2]}_alt"] = datas[4]
    return vcf_datas


def match_nodes_to_vcf(gfa_file: str, vcf_sequences: dict, target_score: float | None = None) -> dict:
    """Tries to do a perfect match between nodes and vcf sequences.

    Args:
        gfa_file (str): A path to a valid GFA-like file
        vcf_sequences (dict): A dict of names:seq issued from a VCF file
        target_score (float | None, optional): If not none, ties to solve a pairwise alignment between nodes and VCF, and uses value as score bound. Defaults to None.
        Setting 1.0 as score is equivalent to None, as it allows only sequences with 100% identity (but None should be provided, as it is faster)
    Returns:
        dict: mapping between variants names from VCF file and node names from GFA, in tuples with their scores.
    """
    if target_score is not None:
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0

    mapping: dict = {name: {} for name in vcf_sequences.keys()}
    with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
        for line in gfa_reader:
            contents = line.split()
            if contents[0] == 'S':
                for name, sequence in vcf_sequences.items():
                    if contents[2] == sequence:
                        mapping[name][contents[1]] = 1.0
                    elif target_score is not None:
                        score: float = (aligner.align(sequence, contents[2])[  # type: ignore
                                        0].score)/max(len(sequence), len(contents[2]))  # type: ignore
                        if score >= target_score:
                            mapping[name][contents[1]] = score
                        else:
                            mapping[name][contents[1]] = float('-inf')
                    else:
                        mapping[name][contents[1]] = float('-inf')

    return mapping


def vcf_heatmap(mapping: dict) -> None:
    """Plots the heatmap, mapping variants to the gfa by their score

    Args:
        mapping (dict): a dict of dicts
    """
    vcf_nodes: list = list(mapping.keys())
    gfa_nodes: list = sorted(list(
        set([gfax for gfanode in list(mapping.values()) for gfax in gfanode.keys()])), key=int)
    datas = [[mapping[vcf][gfa] if vcf in mapping and gfa in mapping[vcf]
              else float('-inf') for gfa in gfa_nodes] for vcf in vcf_nodes]
    graph = heatmap(datas, annot=True, xticklabels=gfa_nodes,  # type: ignore
                    yticklabels=vcf_nodes, vmin=0.0, vmax=1.0,  # type: ignore
                    cmap="Blues", linewidths=1, linecolor='white', square=True)
    graph.set_facecolor('lightsteelblue')
    plt.savefig("vcf_heatmap.png")
    plt.show()


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "gfa", type=str, help="Path to a gfa-like file.")
    parser.add_argument(
        "vcf", type=str, help="Path to a vcf-like file.")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Maps variants to graph.')
    parser.add_argument(
        "-s", "--score", help="Filters by percentage of identity", required=False, type=float, default=None)

    args = parser.parse_args()

    vcf_heatmap(match_nodes_to_vcf(args.gfa, vcf_parser(args.vcf), args.score))
