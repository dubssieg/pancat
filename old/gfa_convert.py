"Converts various GFA files"
from argparse import ArgumentParser, SUPPRESS
from re import sub


def rgfa_to_gfa(input_file: str, output_file: str,
                names_of_haplotypes: list, p_lines: bool = False, keep_tags: bool = False) -> None:
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
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1]),
                                            datas[2], sub('\D', '', datas[3])]+datas[4:]))
                    else:
                        gfa_writer.write(
                            '\n'+'\t'.join([datas[0], sub('\D', '', datas[1]),
                                            datas[2], sub('\D', '', datas[3]), datas[4], datas[5]]))
                    # datas[5] == cigar
                    # datas[6][5:] == origin_sequence
                    link_informations.append(
                        (sub('\D', '', datas[1])+datas[2],
                         sub('\D', '', datas[3])+datas[4], datas[5], datas[6][5:]))
                    edgelist.append((sub('\D', '', datas[1])+datas[2], [], sub(
                        '\D', '', datas[3])+datas[4], int(datas[6][5:])))
                # We don't really know linetype
                else:
                    gfa_writer.write('\n'+'\t'.join(datas))

        if p_lines:
            dnp: dict = create_p_lines(edgelist, names_of_haplotypes)
            # Writing P-lines
            for (origin, path) in dnp.items():
                gfa_writer.write(
                    f"\nP\t{names_of_haplotypes[int(origin)]}\t{','.join(path)}\t*\tSR:i:{origin}")


def create_p_lines(edges: list[tuple], tags: list) -> dict:
    """Assuming edges is a four-part tuple list
    [(edge-a,middle|[],edge-b,sequence-of-edge)]
    computes the reconstruction of genomes.

    Args:
        edges (list[tuple]): edge description
        tags (list): names of each haplotype
    """
    # Init a bank of edges
    bank_of_edges: dict = {int(x): []
                           for x in set([e for _, _, _, e in edges])}
    if len(bank_of_edges) != len(tags):
        raise ValueError(
            f"""Incorrect number of names given the number of haplotypes the file contains.
            In file : {len(bank_of_edges)} / in tags : {len(tags)}""")
    # We iterate until convergence, which happens when no sequence could be merged
    while edges:
        match iterate_edges(edges, bank_of_edges):
            case None:
                # We need to complete paths that are not part
                # of the reference with the end and start of the reference sequence
                paths: dict = {str(seq): [start]+mid+[end]
                               for (start, mid, end, seq) in bank_of_edges.values()}
                refseq: list = paths['0']
                for identifier, path in paths.items():
                    if identifier != '0':
                        paths[identifier] = refseq[0:refseq.index(
                            paths[identifier][0])] + path + refseq[refseq.index(paths[identifier][-1])+1:]
                return paths
            case (new_edge, pos_to_edit):
                bank_of_edges[pos_to_edit] = new_edge
    raise ValueError("Terminated without building correct paths.")


def iterate_edges(edges: list[tuple], bank_of_edges: dict[int, list]) -> tuple | None:
    """Iterates over the edges and chains them

    Args:
        edges (list[tuple]): a list of all edges inside the graph
        bank_of_edges (list[tuple]): a list of all newformed edges

    Returns:
        tuple: (new_edge,position_to_edit,position_to_del)
    """
    for edge_a in edges:
        for identifier, edge_b in bank_of_edges.items():
            match (edge_a, edge_b):
                case ((_, _, _, origin), []) if origin == identifier:
                    # Empty sequence
                    return edge_a, identifier
                case ((start, internal_left, bridge_a, seq_a), (bridge_b, internal_right, end, seq_b)) if seq_a == seq_b and bridge_a == bridge_b and end not in [start]+internal_left:
                    # Same seq and chained
                    return (start, internal_left + [bridge_a] + internal_right, end, seq_a), seq_a
                case ((bridge_a, internal_right, end, seq_a), (start, internal_left, bridge_b, seq_b)) if seq_a == seq_b and bridge_a == bridge_b and end not in [start]+internal_left:
                    # same seq and chained and reversed
                    return (start, internal_left + [bridge_a] + internal_right, end, seq_a), seq_a
                case ((start, internal_left, bridge_a, 0), (bridge_b, internal_right, end, seq)) if bridge_a == bridge_b and end not in [start]+internal_left:
                    # Different seq but chained to ref
                    return (start, internal_left + [bridge_a] + internal_right, end, seq), seq
                case ((start, internal_left, bridge_a, seq), (bridge_b, internal_right, end, 0)) if bridge_a == bridge_b and end not in [start]+internal_left:
                    # Different seq but chained to ref
                    return (start, internal_left + [bridge_a] + internal_right, end, seq), seq
                case ((bridge_a, internal_right, end, 0), (start, internal_left, bridge_b, seq)) if bridge_a == bridge_b and end not in [start]+internal_left:
                    # Different seq but chained to ref and reversed
                    return (start, internal_left + [bridge_a] + internal_right, end, seq), seq
                case ((bridge_a, internal_right, end, seq), (start, internal_left, bridge_b, 0)) if bridge_a == bridge_b and end not in [start]+internal_left:
                    # Different seq but chained to ref and reversed
                    return (start, internal_left + [bridge_a] + internal_right, end, seq), seq
                case _:
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
    parser.add_argument(
        "-n",
        "--haplotypes_names",
        help="Give one name per haplotype, ordered as the same order you did for your rGFA. First item of the list is the reference. Required if you ask for P-lines.",
        nargs='+',
        default=[]
    )
    args = parser.parse_args()

    rgfa_to_gfa(
        args.file,
        f"{args.file.split('.')[0]}_gfa1.gfa",
        names_of_haplotypes=args.haplotypes_names,
        p_lines=args.plines,
        keep_tags=args.keep)
