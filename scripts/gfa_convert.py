"Converts various GFA files"
from argparse import ArgumentParser, SUPPRESS
from re import sub


def rgfa_to_gfa(input_file: str, output_file: str, p_lines: bool = False, keep_tags: bool = False) -> None:
    """Converts rGFA (minigraph) to GFA1 files by adding header and P-lines
    This process is not lossless !!! We lost tags in-between.

    Args:
        input_file (str): a path to a rGFA file
        output_file (str): a path where file should be written
    """
    # Cleaning file
    edgelist: list = []
    number_of_nodes: int = 1
    with open(output_file, "w", encoding="utf-8") as gfa_writer:
        pass
    link_informations: list = []
    with open(output_file, "a", encoding="utf-8") as gfa_writer:
        # Header
        gfa_writer.write("H\tVN:Z:1.0")
        with open(input_file, "r", encoding="utf-8") as gfa_reader:
            for line in gfa_reader:
                datas: list = line.split()
                # Segment
                if datas[0] == 'S':
                    number_of_nodes += 1
                    if keep_tags:
                        gfa_writer.write(
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1])]+datas[2:]))
                    else:
                        gfa_writer.write(
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1]), datas[2]]))
                # Line
                elif datas[0] == 'L':
                    if keep_tags:
                        gfa_writer.write(
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1]), datas[2], sub('\D', '', datas[3])]+datas[4:]))
                    else:
                        gfa_writer.write(
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1]), datas[2], sub('\D', '', datas[3]), datas[4], datas[5]]))
                    # datas[5] == cigar
                    # datas[6][5:] == origin_sequence
                    link_informations.append(
                        (sub('\D', '', datas[1])+datas[2], sub('\D', '', datas[3])+datas[4], datas[5], datas[6][5:]))
                    edgelist.append((sub('\D', '', datas[1])+datas[2], [], sub(
                        '\D', '', datas[3])+datas[4], int(datas[6][5:])))
                # We don't really know linetype
                else:
                    gfa_writer.write('\n'+'\t'.join(datas))

        if p_lines:
            dnp: dict = create_p_lines(edgelist)
            # Writing P-lines
            for path_number, (origin, path) in enumerate(dnp.items()):
                gfa_writer.write(
                    f"\nP\t{path_number+number_of_nodes}\t{','.join(path)}\t*\tSR:i:{origin}")


def create_p_lines(edges: list[tuple]) -> dict:
    """Assuming edges is a four-part tuple list
    [(edge-a,middle|[],edge-b,sequence-of-edge)]
    computes the reconstruction of genomes.

    Args:
        edges (list[tuple]): edge description
    """
    # Edges in reference are shared.
    bank_of_edges: list[tuple] = [(start, middle, stop, seq)
                                  for (start, middle, stop, seq) in edges if seq == 0]

    # We iterate until convergence, which happens when no sequence could be merged
    while True:
        match iterate_edges(edges, bank_of_edges):
            case None:
                # We need to complete paths that are not part of the reference with the end and start of the reference sequence
                paths: dict = {str(seq): [start]+mid+[end]
                               for (start, mid, end, seq) in edges}
                refseq: list = paths['0']
                for id, path in paths.items():
                    if id != '0':
                        paths[id] = refseq[0:refseq.index(
                            paths[id][0])] + path + refseq[refseq.index(paths[id][-1])+1:]
                return paths
            case (new_edge, pos_to_edit, pos_to_suppress):
                edges[pos_to_edit] = new_edge
                if not pos_to_suppress > len(edges):
                    edges.pop(pos_to_suppress)


def iterate_edges(edges: list[tuple], bank_of_edges: list[tuple]) -> tuple | None:
    """Iterates over the edges and chains them

    Args:
        edges (list[tuple]): a list of all edges inside the graph
        bank_of_edges (list[tuple]): a list of all edges from the reference inside the graph

    Returns:
        tuple: (new_edge,position_to_edit,position_to_del)
    """
    for i, edge_a in enumerate(edges):
        for j, edge_b in enumerate(edges + bank_of_edges):
            if edge_a != edge_b:
                match (edge_a, edge_b):
                    case ((start, internal_left, bridge_a, seq_a), (bridge_b, internal_right, end, seq_b)) | ((bridge_a, internal_right, end, seq_a), (start, internal_left, bridge_b, seq_b)) if seq_a == seq_b and bridge_a == bridge_b:
                        # Same seq and chained, or same seq and chained and reversed
                        return (start, internal_left + [bridge_a] + internal_right, end, seq_a), i, j
                    case ((start, internal_left, bridge_a, 0), (bridge_b, internal_right, end, seq)) | ((start, internal_left, bridge_a, seq), (bridge_b, internal_right, end, 0)) | ((bridge_a, internal_right, end, 0), (start, internal_left, bridge_b, seq)) | ((bridge_a, internal_right, end, seq), (start, internal_left, bridge_b, 0)) if bridge_a == bridge_b:
                        # Different seq but chained to ref, or different seq but chained to ref and reversed
                        return (start, internal_left + [bridge_a] + internal_right, end, seq), i, j
                    case _:
                        # No edge can be merge
                        pass


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="rGFA file")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='This script aims to convert rGFA files issued from minigraph to GFA1 compatible format. It implies to rename nodes and add P-lines if asked for.')
    parser.add_argument(
        "-p", "--plines", help="Asks to calculate p-lines for graph", action='store_true')
    parser.add_argument(
        "-k", "--keep", help="Keeps rGFA-specific tags in output", action='store_true')
    args = parser.parse_args()

    rgfa_to_gfa(
        args.file,
        f"{args.file.split('.')[0]}_gfa1.gfa",
        p_lines=args.plines,
        keep_tags=args.keep)
