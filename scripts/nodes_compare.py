'Compares nodes within a graph'
from itertools import combinations
from argparse import ArgumentParser, SUPPRESS
from re import sub
from networkx import MultiDiGraph, compose_all, add_path, isolates
from gfatypes import LineType, Record, GfaStyle
from pyvis.network import Network
import matplotlib.pyplot as plt


def are_nodes_matching(gfa_files: list) -> list:
    """Checks if some segments within graph are duplicates of another,
    and how nodes are equal in multiple graphs

    Args:
        gfa_files (list): paths to GFA-like files

    Returns:
        list: nodes grouped together if they are the same
    """
    nodes: dict = {}
    for gfa_file in gfa_files:
        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            nodes: dict = {
                **nodes,
                **{
                    gfa_file.split('/')[-1].split('.')[0]+'_'+sub('\D', '', seq.split()[1]):
                    seq.split()[2].strip()
                    for seq in gfa_reader
                    if seq.split()[0] == 'S'
                }
            }
    assoc: list = []
    stack: list = list(nodes.keys())
    while stack:
        pos: list = [i for i, x in enumerate(
            stack) if nodes[x] == nodes[stack[0]]]
        assoc.append([stack[x] for x in pos])
        for offset, position in enumerate(pos):
            stack.pop(position-offset)
    return assoc


def show_identity(gfa_files: list, gfa_versions: list, colors: list) -> MultiDiGraph:
    """Given a list of files, inits the networks for each graph,
    merges them, and display dotted links for identical nodes accross graph

    Args:
        gfa_files (list): a list of GFA-like files
        gfa_versions (list): user-assumed GFA subversions
        colors (list): one color per pangenome graph
    """
    if any(len(x) != len(y) for (x, y) in combinations(locals().values(), 2)):
        raise ValueError("Some arguments does not have the same length.")
    if 'rGFA' in gfa_versions:
        combined_view: MultiDiGraph = compose_all(
            [init_simple_graph(gfa, gfa_versions[i], colors[i]) for i, gfa in enumerate(gfa_files)])
    else:
        combined_view: MultiDiGraph = compose_all(
            [init_graph(gfa, gfa_versions[i], colors[i]) for i, gfa in enumerate(gfa_files)])
    identities: list = are_nodes_matching(gfa_files)
    for idies in identities:
        for (left, right) in combinations(idies, 2):
            combined_view.add_edge(
                left, right, style='dashed', color='forestgreen', arrows='')
    return combined_view


def init_simple_graph(gfa_file: str, gfa_version: str, color: str) -> MultiDiGraph:
    """Initializes the graph without displaying it

    Args:
        gfa_file (str): GFA-like input file
        gfa_version (str): user-assumed GFA subversion
        color (str): color to apply to all nodes and edges of graph

    Raises:
        ValueError: Occurs if graph specified format isn't correct given the file
        NitImplementedError : Occurs if the function is currently not impelemented yet

    Returns:
        MultiDiGraph: a graph representing the given pangenome
    """
    graph = MultiDiGraph()

    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            gfa_line: Record = Record(line, gfa_version)
            match gfa_line.linetype:
                case LineType.SEGMENT:
                    graph.add_node(
                        f"{gfa_file.split('/')[-1].split('.')[0]}_{gfa_line.line.name}",
                        title=f"{gfa_line.line.length} bp.",
                        color=color
                    )
                case LineType.LINE:
                    graph.add_edge(
                        f"{gfa_file.split('/')[-1].split('.')[0]}_{gfa_line.line.start}",
                        f"{gfa_file.split('/')[-1].split('.')[0]}_{gfa_line.line.end}",
                        color=color,
                        weight=4
                    )
    return graph


def get_palette(number_of_colors: int) -> list:
    """Returns a number_of_colors-sized palette, as a list,
    that one can access with colors[i].

    Args:
        number_of_colors (int): number of colors needed

    Returns:
        list: palette of colors
    """
    colormap = plt.cm.viridis  # type:ignore
    number_of_colors = min(colormap.N, number_of_colors)
    return [colormap(int(x*colormap.N/number_of_colors)) for x in range(number_of_colors)]


def init_graph(gfa_file: str, gfa_version: str, color: str, n_aligns: int = 3) -> MultiDiGraph:
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

    cmap: list = ['slateblue', 'darkslateblue', 'mediumslateblue']

    visited_paths: int = 0
    version: GfaStyle = GfaStyle(gfa_version)
    if version == GfaStyle.RGFA:
        raise ValueError(
            "This function is not rGFA compatible. Please convert first.")
    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            gfa_line: Record = Record(line, gfa_version)
            match gfa_line.linetype:
                case LineType.SEGMENT:
                    graph.add_node(
                        f"{gfa_file.split('/')[-1].split('.')[0]}_{gfa_line.line.name}",
                        title=f"{gfa_line.line.length} bp.",
                        color=color
                    )
                case LineType.WALK:
                    if not gfa_line.line.idf == '_MINIGRAPH_':
                        add_path(
                            graph,
                            [f"{gfa_file.split('/')[-1].split('.')[0]}_{node}" for (
                                node, _) in gfa_line.line.path],
                            title=gfa_line.line.name,
                            color=cmap[visited_paths],
                            weight=4
                        )
                        visited_paths += 1
                case LineType.PATH:
                    add_path(
                        graph,
                        [f"{gfa_file.split('/')[-1].split('.')[0]}_{node}" for (
                            node, _) in gfa_line.line.path],
                        title=gfa_line.line.name,
                        color=cmap[visited_paths],
                        weight=4
                    )
                    visited_paths += 1
    graph.remove_nodes_from(list(isolates(graph)))
    return graph


def display_graph(graph: MultiDiGraph) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
    """
    graph_visualizer = Network(
        height='1080px', width='100%', directed=True)
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    try:
        graph_visualizer.show("proximity_graph.html")
    except FileNotFoundError:
        pass


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to one or more gfa-like file(s).", nargs='+')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Evaluates identity of nodes within and across graphs')
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
    args = parser.parse_args()

    display_graph(show_identity(args.file, args.gfa_version,
                  ['rebeccapurple', 'crimson', 'orchid']))
