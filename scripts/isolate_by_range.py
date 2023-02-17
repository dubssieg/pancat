"Requires PO tag to work."
from argparse import ArgumentParser, SUPPRESS
from gfagraphs import Segment, Graph, Walk, Path, GfaStyle


def nodes_to_keep(list_of_nodes: list[Segment], start_point: int, end_point: int) -> set:
    """Given a range of positions, only keeps nodes that are embed within that range

    Args:
        list_of_nodes (list[Segment]): all nodes from the graph
        start_point (int): a starting point
        end_point (int): a ending point

    Returns:
        set: all nodes names that participates in the range
    """
    return set(node.datas["name"] for node in list_of_nodes if any(
        [start >= start_point and end <= end_point for _, (start, end, _) in node.datas['PO'].items()]))


def isolate(gfa_graph: str, output: str, start_point: int, end_point: int, gfa_style: str) -> None:
    """Given a GFA file, extracts subgraph within position range. Requires PO offset (other script).

    Args:
        gfa_graph (str): input graph
        output (str): output graph
        start_point (int): bp position
        end_point (int): bp position
        gfa_style (str): gfa subtype
    """
    graph: Graph = Graph(gfa_graph, gfa_style)
    extracted_nodes: set = nodes_to_keep(
        graph.segments, start_point, end_point)
    embed_paths: list[Walk | Path] = graph.get_path_list()
    del graph
    with open(output, 'w', encoding='utf-8') as gfa_writer:
        with open(gfa_graph, 'r', encoding='utf-8') as gfa_reader:
            for line in gfa_reader:
                if (x := line.split())[0] == 'S' and x[1] in extracted_nodes:
                    gfa_writer.write(line)
                elif x[0] == 'L' and x[1] in extracted_nodes and x[3] in extracted_nodes:
                    gfa_writer.write(line)
                elif x[0] == 'H':
                    if GfaStyle(gfa_style) != GfaStyle.RGFA:
                        gfa_writer.write(line)
    new_paths: dict = {path.datas["name"]: [
        (node, orient) for node, orient in path.datas["path"] if node in extracted_nodes] for i, path in enumerate(embed_paths)}
    match GfaStyle(gfa_style):
        case GfaStyle.GFA1:
            with open(output, 'a', encoding='utf-8') as gfa_writer:
                gfa_writer.writelines(
                    [f"P\t{key}\t{','.join([node+orientation.value for node,orientation in val])}\t*\n" for key, val in new_paths.items()])
        case GfaStyle.GFA1_1:
            with open(output, 'a', encoding='utf-8') as gfa_writer:
                gfa_writer.writelines(
                    [f"W\t{key}\t{i}\t{key}\t0\t0\t{''.join([orientation.value+node for node,orientation in val]).replace('+', '>').replace('-', '<')}\t*\n" for i, (key, val) in enumerate(new_paths.items())])
        case GfaStyle.RGFA:
            pass
        case _:
            raise NotImplementedError(
                "Functionnality currently not implemented.")


if __name__ == '__main__':
    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument("out", type=str, help="Output path (with extension)")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Extracts subgraphs from a GFA graph')
    parser.add_argument(
        "-g",
        "--gfa_version",
        help="Tells the GFA input style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
    )
    parser.add_argument(
        '-s',
        '--start',
        type=int,
        help='To specifiy a starting point (in bp) to create a subgraph',
    )
    parser.add_argument(
        '-e',
        '--end',
        type=int,
        help='To specifiy a end point (in bp) to create a subgraph',
    )
    args = parser.parse_args()

    isolate(args.file, args.out, args.start, args.end, args.gfa_version)
