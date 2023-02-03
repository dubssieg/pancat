"Creates a graph we can navigate in."
from argparse import ArgumentParser, SUPPRESS
from math import log10
from networkx import MultiDiGraph, isolates
from mycolorpy import colorlist
from pyvis.network import Network
from gfatypes import LineType, Record, GfaStyle


def html_graph(graph: MultiDiGraph, job_name: str) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        gfa_file (str): path to a rGFA file
        debug (bool, optional): plots less nodes in graph. Defaults to False.
        plines (bool, optional) : plots the P-lines as paths on the graph. Defaults to False.
    """

    graph_visualizer = Network(
        height='1000px', width='100%', directed=True)
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    try:
        graph_visualizer.show(f"{job_name}_graph.html")
    except FileNotFoundError:
        # Path indicated for file may not be correct regarding the lib but writes .html anyways, so ignore ^^
        pass

    with open(f"{job_name}_graph.html", "r", encoding="utf-8") as html_reader:
        outfile = html_reader.readlines()
        # <img src='{gfa_file.split('.')[0].split('/')[-1]}_cbar.png' align='center' rotate='90'>
    outfile[10] = f"<h1>Graph for <b>{job_name}</b></h1>"
    with open(f"{job_name}_graph.html", "w", encoding="utf-8") as html_writer:
        html_writer.writelines(outfile)


def init_graph(gfa_file: str, gfa_version: str, n_aligns: int) -> MultiDiGraph:
    """Initializes the graph without displaying it

    Args:
        gfa_file (str): GFA-like input file
        gfa_version (str): user-assumed GFA subversion
        n_aligns (int): number of distinct origin sequences

    Raises:
        ValueError: Occurs if graph specified format isn't correct given the file
        NitImplementedError : Occurs if the function is currently not impelemented yet

    Returns:
        MultiDiGraph: a graph representing the given pangenome
    """
    graph = MultiDiGraph()

    cmap: list = colorlist.gen_color_normalized(
        cmap="rainbow", data_arr=[i/n_aligns for i in range(n_aligns)])

    visited_paths: int = 0
    version: GfaStyle = GfaStyle(gfa_version)
    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            gfa_line: Record = Record(line, gfa_version)
            match (gfa_line.linetype, version):
                case LineType.SEGMENT, _:
                    graph.add_node(
                        gfa_line.line.name,
                        title=f"{gfa_line.line.length} bp.",
                        size=log10(gfa_line.line.length),
                        color='darkslateblue'
                    )
                case LineType.WALK, _:
                    if not gfa_line.line.idf == '_MINIGRAPH_':
                        for i in range(len(gfa_line.line.path)-1):
                            left_node, left_orient = gfa_line.line.path[i]
                            right_node, right_orient = gfa_line.line.path[i+1]
                            graph.add_edge(
                                left_node,
                                right_node,
                                title=gfa_line.line.name,
                                color=cmap[visited_paths],
                                label=f"{left_orient.value}/{right_orient.value}"
                            )
                        visited_paths += 1
                case LineType.PATH, _:
                    for i in range(len(gfa_line.line.path)-1):
                        left_node, left_orient = gfa_line.line.path[i]
                        right_node, right_orient = gfa_line.line.path[i+1]
                        graph.add_edge(
                            left_node,
                            right_node,
                            title=gfa_line.line.name,
                            color=cmap[visited_paths],
                            label=f"{left_orient.value}/{right_orient.value}"
                        )
                    visited_paths += 1
                case LineType.LINE, GfaStyle.RGFA:
                    graph.add_edge(
                        gfa_line.line.start,
                        gfa_line.line.end,
                        title=str(gfa_line.line.origin),
                        color=cmap[gfa_line.line.origin],
                        label=gfa_line.line.orientation
                    )
    graph.remove_nodes_from(list(isolates(graph)))
    return graph


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument("job_name", type=str,
                        help="Job identifier for output (ex : chr3_graph)")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Creates a representation of a pangenome graph.')
    parser.add_argument(
        "-n", "--number_alignments", help="Gives the number of origin sequences", required=True, type=int)
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    args = parser.parse_args()

    pangenome_graph: MultiDiGraph = init_graph(
        gfa_file=args.file,
        gfa_version=args.gfa_version,
        n_aligns=args.number_alignments
    )

    html_graph(
        pangenome_graph,
        args.job_name
    )
