"Adapting BubbleGun bfs for yielding paths/walks"
from argparse import ArgumentParser, SUPPRESS
from collections import deque
from BubbleGun.bfs import bfs
from BubbleGun.Graph import Graph
from gfatypes import GfaStyle, Record


def bfs_step(gfa_graph: Graph, starting_node: str, number_of_nodes: int) -> set:
    """Extracts nodes from GFA graph around a node position

    Args:
        input_file (str): a GFA file
        output_path (str): destination for subset of nodes
        starting_nodes (list): a node that should be the center of isolates
        number_of_nodes (int): number of nodes we want to isolate around each node position

    Returns:
        set : the nodes that has been extracted
    """
    return bfs(gfa_graph, starting_node, number_of_nodes)


def paths_step(origin_file: str, file_to_edit: str, extracted_nodes: set, gfa_version: str, gfa_output: str) -> None:
    """Aims to add subpaths (walks/paths) issued from the origin graph
    to the new file obtained by BubbleGun's bfs function

    Args:
        origin_file (str): the file we extracted nodes from
        file_to_edit (str): the output file for the bfs step
        extracted_nodes (set): the nodes that are emmbed in paths to search for
        gfa_version (str): identifier for gfa input version
        gfa_output (str): identifier for gfa output version
    """
    output_type: GfaStyle = GfaStyle(gfa_output)
    if GfaStyle(gfa_version) == GfaStyle.RGFA or output_type == GfaStyle.RGFA:
        raise NotImplementedError(
            "Nodes can be extracted, but paths could not be determined.")
    embed_paths: list[Record] = list()
    with open(file_to_edit, 'w', encoding='utf-8') as gfa_writer:
        with open(origin_file, 'r', encoding='utf-8') as gfa_reader:
            for line in gfa_reader:
                if (x := line.split())[0] == 'S' and x[1] in extracted_nodes:
                    gfa_writer.write(line)
                elif x[0] == 'L' and x[1] in extracted_nodes and x[3] in extracted_nodes:
                    gfa_writer.write(line)
                elif x[0] == 'P' or x[0] == 'W':
                    embed_paths += [Record(line, gfa_version)]
                elif x[0] == 'H':
                    gfa_writer.write(line)
    for i, path in enumerate(embed_paths):
        path.line.path = deque(path.line.path, maxlen=len(path.line.path))
        try:
            while not path.line.path[0][0] in extracted_nodes:
                _ = path.line.path.popleft()
            while not path.line.path[-1][0] in extracted_nodes:
                _ = path.line.path.pop()
        except IndexError:
            # No node inside path is in set, removing path
            embed_paths[i] = None
    match output_type:
        case GfaStyle.GFA1:
            with open(file_to_edit, 'a', encoding='utf-8') as gfa_writer:
                gfa_writer.writelines(
                    [f"P\t{path.line.name}\t{','.join([node+orientation.value for node,orientation in path.line.path])}\t*\n" for path in embed_paths if path is not None])
        case GfaStyle.GFA1_1:
            with open(file_to_edit, 'a', encoding='utf-8') as gfa_writer:
                gfa_writer.writelines(
                    [f"W\t{path.line.idf}\t{i}\t{path.line.name}\t0\t0\t{''.join([orientation.value+node for node,orientation in path.line.path]).replace('+', '>').replace('-', '<')}\t*\n" for i, path in enumerate(embed_paths) if path is not None])
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
        "-o",
        "--gfa_output",
        help="Tells the GFA output style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
    )
    parser.add_argument(
        '-s',
        '--start_node',
        type=str,
        help='To specifiy a starting node on reference to create a subgraph',
        nargs='+'
    )
    parser.add_argument("-c", "--count", type=int,
                        help="Number of nodes around each starting point")
    args = parser.parse_args()

    graph: Graph = Graph(args.file)
    for i, node in enumerate(args.start_node):
        output: str = f"{args.out.split('.')[0]}_{i}.gfa" if len(
            args.start_node) > 1 else args.out
        nodes: set = bfs_step(graph, node, args.count)
        print(f"Selected {len(nodes)} nodes.")
        paths_step(args.file, output, nodes,
                   args.gfa_version, args.gfa_output)
