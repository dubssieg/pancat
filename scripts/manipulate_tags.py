from argparse import ArgumentParser, SUPPRESS
from gfagraphs import Graph


if __name__ == "__main__":
    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Extracts tags from a GFA graph')
    parser.add_argument(
        "-g",
        "--gfa_version",
        help="Tells the GFA input style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
    )
    parser.add_argument(
        '-n',
        '--node',
        type=str,
        help='To specifiy a node name',
        nargs='+'
    )
    args = parser.parse_args()

    graph: Graph = Graph(args.file, args.gfa_version)
    nodes = {n.datas['name']: n for n in graph.segments}
    for anode in args.node:
        try:
            print(nodes[anode].datas["PO"])
        except KeyError:
            print(f"No node matching name {anode}.")
