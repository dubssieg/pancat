"Isolates from tags from fasta file"
from argparse import ArgumentParser, SUPPRESS
from json import dump
from os import remove, path
from Bio import SeqIO


def export_mapping(paf_file: str, save: bool = False, threshold: int = 4000000) -> list:
    """Exports a list of dicts

    Args:
        paf_file (str): path to a paf-formated file
        save (bool, optional): if should save to disk. Defaults to False.

    Returns:
        list: dicts containing for each scaffold mapping between query and reference
    """
    mapping: list[dict] = sorted(
        [
            {
                'query_seq_name': l.split()[0],
                'ref_seq_name':l.split()[5],
                'sequence_length':int(l.split()[10])
            }
            for l in open(paf_file, "r", encoding="utf-8")
            if int(l.split()[10]) >= threshold
        ],
        key=lambda x: x['sequence_length']
    )[::-1]

    if save:
        dump(mapping, open(
            f"{paf_file.split('.')[0]}_chromosom_mapping.json", "w", encoding="utf-8"))
    return mapping


def isolate_scaffolds(fasta_file: str, out_file: str, paf_file: str, chromosom: str) -> None:
    """Isolate scaffolds based upon correspondances from reference to query

    Args:
        fasta_file (str): query file
        out_file (str): path for output
        paf_file (str): mapping of reference against query
        chromosom (str): chromosom identifier, name used on reference file
    """
    nb_seq, length = 0, 0
    if path.exists(out_file):
        print(f"Erasing {out_file}")
        remove(out_file)
    with open(out_file, 'a', encoding="utf-8") as handler:
        for fasta in SeqIO.parse(open(fasta_file, 'r', encoding="utf-8"), 'fasta'):
            if fasta.id in [x['query_seq_name'] for x in export_mapping(paf_file=paf_file, save=False) if x['ref_seq_name'] == chromosom]:
                print(f"Writing {fasta.id} to {out_file}...")
                SeqIO.write(fasta, handler, 'fasta')
                nb_seq += 1
                length += len(fasta.seq)
        print(
            f"Job done, {nb_seq} sequences were written for a total of {length} bases!")


def isolate_sequence(fasta_file: str, out_file: str, chromosom: str) -> None:
    """Isolate scaffolds based upon correspondances from reference to query

    Args:
        fasta_file (str): query file
        out_file (str): path for output
        paf_file (str): mapping of reference against query
        chromosom (str): chromosom identifier, name used on reference file
    """
    nb_seq, length = 0, 0
    if path.exists(out_file):
        print(f"Erasing {out_file}")
        remove(out_file)
    with open(out_file, 'w', encoding="utf-8") as handler:
        for fasta in SeqIO.parse(open(fasta_file, 'r', encoding="utf-8"), 'fasta'):
            if fasta.id == chromosom:
                print(f"Writing {fasta.id} to {out_file}...")
                SeqIO.write(fasta, handler, 'fasta')
                nb_seq += 1
                length += len(fasta.seq)
        print(
            f"Job done, {nb_seq} sequences were written for a total of {length} bases!")


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="fasta-like file")
    parser.add_argument(
        "out", type=str, help="fasta-like output")
    parser.add_argument(
        "paffile", type=str, help="paf-like file")
    parser.add_argument("-c",
                        "--chromosom", type=str, help="name of assembly on reference sequence", default=None)
    args = parser.parse_args()
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Extracts from a fasta-like file all sequences in a query assembly given a mapping to a reference and an identifier on reference.')

    if args.chromosom is None:
        export_mapping(paf_file=args.file, save=True)
    else:
        isolate_scaffolds(fasta_file=args.file, out_file=args.out,
                          paf_file=args.paffile, chromosom=args.chromosom)
