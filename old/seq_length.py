"Computes and print the length of each sequence in a fasta file; Also compares fasta files"
from argparse import ArgumentParser, SUPPRESS
from Bio import SeqIO


def get_assemblies_size(fasta_file: str) -> dict:
    """Show the size of each element inside a fasta file

    Args:
        fasta_file (str): a path to a fasta like file

    Returns:
        dict: a dict of dict mapping each seq-id to its number of bases
    """
    return {fasta_file: {fasta.id: len(fasta.seq) for fasta in SeqIO.parse(open(fasta_file, encoding="utf-8"), 'fasta')}}


def compare_fasta_sequences(fasta_files: list) -> list[bool]:
    """Checks if all fastas in a list have the same number of sequences, and if all sequences are equal.

    Args:
        fasta_files (list): a list of fasta-like files

    Returns:
        list[bool]: if sequence is the same across all given files.
    """
    fastas: list = []
    for fasta_file in fasta_files:
        fastas.append([(fasta.id, fasta.seq) for fasta in SeqIO.parse(
            open(fasta_file, encoding="utf-8"), 'fasta')])
    return [all([seq.strip() == x.strip() for fasta in fastas for id, x in fasta if fasta_id == id]) and all([len(fastas[0]) == len(fasta) for fasta in fastas[1:]]) for fasta_id, seq in fastas[0]]


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str,
                        help="fasta-like file", nargs='+')
    args = parser.parse_args()
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Calculates the length of all sequences in file.')

    for file in args.file:
        print(get_assemblies_size(fasta_file=file))

    if all(compare_fasta_sequences(args.file)):
        print('All sequences are identical')
    else:
        print('Sequences are not the same !')
