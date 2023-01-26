'Verifies to which sequences are mapped each node of a GFA.'
from itertools import chain
from re import finditer
from argparse import ArgumentParser, SUPPRESS
from Bio import SeqIO
import matplotlib.pyplot as plt


def get_nodes(input_file: str) -> dict:
    """Gets all nodes from a GFA-like file

    Args:
        input_file (str): a GFA-like formatted file with S-lines

    Returns:
        dict: mapping segment_name:segment_sequence
    """
    with open(input_file, 'r', encoding='utf-8') as gfa_reader:
        return {seq.split()[1]: seq.split()[2] for seq in gfa_reader if seq.split()[0] == 'S'}


def get_fastas(input_file: str) -> dict:
    """Gets all sequences from a FASTA-like file

    Args:
        input_file (str): a FASTA file

    Returns:
        dict: mapping id:sequence
    """
    return {
        fasta.id: fasta.seq for fasta in SeqIO.parse(open(input_file, encoding="utf-8"), 'fasta')
    }


def get_lengths(input_file: str) -> dict:
    """Gets all lengths of sequences from a FASTA-like file

    Args:
        input_file (str): a FASTA file

    Returns:
        dict: mapping id:sequence
    """
    return {
        fasta.id: len(fasta.seq) for fasta in SeqIO.parse(open(input_file, encoding="utf-8"), 'fasta')
    }


def mapper(gfa_file: str, fasta_file: str) -> dict:
    """Maps presence of node inside sequence

    Args:
        gfa_file (str): gfa input file
        fasta_file (str): fasta input file

    Returns:
        dict: id of sequences each node is in
    """
    return {
        name:
        [
            id for id, fasta in get_fastas(fasta_file).items() if node in fasta
        ]
        for name, node in get_nodes(gfa_file).items()
    }


def exact_mapper(gfa_file: str, fasta_file: str, threshold: int = 10) -> dict:
    """Maps presence of node inside sequence

    Args:
        gfa_file (str): gfa input file
        fasta_file (str): fasta input file
        threshold (int): defaults to 50, minimal length of node to plot.

    Returns:
        dict: id of sequences each node is in, and position of each occurence
    """
    return {
        name:
        {
            id: [
                # f'(?={node})'
                (m.start(), len(node)) for m in finditer(node, str(fasta))
            ]
            for id, fasta in get_fastas(fasta_file).items() if node in fasta
        }
        for name, node in get_nodes(gfa_file).items() if len(node) >= threshold
    }


def show_alignments(map: dict, sequences: dict, nodes_to_highlight: list = []) -> None:
    """_summary_

    Args:
        map (dict): _description_
        sequences (dict): _description_
    """

    fig, axs = plt.subplots(nrows=len(sequences), ncols=1, figsize=(
        14, 10), sharex=True, sharey=False)

    all_sequences: list = list(set(chain(*map.values())))
    all_nodes: list = list(map.keys())

    for seq_number, seq_name in enumerate(all_sequences):
        axs[seq_number].set_title(f"Alignment for {seq_name}")
        axs[seq_number].set_xlim(0, sequences[seq_name])
        axs[seq_number].set_ylim(-0.5, len(all_nodes)+1)
        axs[seq_number].plot([0, sequences[seq_name]], [
                             0, 0], linewidth=6, color='coral')
        for i, (name, mapping) in enumerate(map.items()):
            for seq_id, positions in mapping.items():
                if seq_id == seq_name:
                    datas = [x for x, _ in positions]  # positions
                    lengths = [y for _, y in positions]  # horizontal lengths
                    # vertical positions
                    offsets = [i+1 for _ in range(len(datas))]
                    if len(datas) > 0:
                        for i, data in enumerate(datas):
                            if name in nodes_to_highlight:
                                axs[seq_number].plot(
                                    [data, data+lengths[i]], [offsets[i], offsets[i]], linewidth=1.5, color='slateblue')
                            else:
                                axs[seq_number].plot(
                                    [data, data+lengths[i]], [offsets[i], offsets[i]], linewidth=1.5, color='coral')
        axs[seq_number].set_yticks(
            [i for i in range(len(all_nodes)+1)], labels=[seq_name]+all_nodes)

    plt.show()


def dotgrid_plot(map: dict, nodes_to_highlight: list = []) -> None:
    """Displays the alignments as a grid of dots

    Args:
        map (dict): mapping node-to-sequence
        nodes_to_highlight (list, optional): a list of nodes to display in orange. Defaults to [].
    """

    all_sequences: list = list(set(chain(*map.values())))
    all_nodes: list = list(map.keys())

    plt.xticks([i for i in range(len(all_sequences))], labels=all_sequences)
    plt.xlabel("Origin sequences")
    plt.ylabel("Nodes")
    plt.yticks([i for i in range(len(all_nodes))], labels=all_nodes)

    for i, node in enumerate(all_nodes):
        for j, seq in enumerate(all_sequences):
            if seq in map[node]:
                if node in nodes_to_highlight:
                    plt.scatter(x=j, y=i, c='darkorange')
                else:
                    plt.scatter(x=j, y=i, c='slateblue')
    plt.show()


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "gfa", type=str, help="gfa-like file")
    parser.add_argument(
        "fasta", type=str, help="fasta-like file")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Verifies to which sequences are mapped each node of a GFA.')
    args = parser.parse_args()

    fmap = mapper(args.gfa, args.fasta)
    print(fmap)
    dotgrid_plot(fmap, nodes_to_highlight=['13', '14', '15', '16'])

    show_alignments(exact_mapper(args.gfa, args.fasta),
                    get_lengths(args.fasta), nodes_to_highlight=['13', '14', '15', '16'])
