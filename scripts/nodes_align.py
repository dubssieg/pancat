'Compares nodes within a graph'
from itertools import combinations, chain
from argparse import ArgumentParser, SUPPRESS
from networkx import MultiDiGraph, compose_all, add_path, isolates
from gfatypes import LineType, Record, GfaStyle
from pyvis.network import Network
import matplotlib.pyplot as plt
from Bio import Align
from datetime import datetime


def node_aligner(node: str, nodes_to_align: list) -> list[float]:
    """Given a node and a list of nodes, compules all pairwise alignments and return their scores

    Args:
        node (str): node to compare sequences to
        nodes_to_align (list): all nodes to align

    Returns:
        list[float]: a list of scores
    """
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    scores: list[float] = [aligner.align(
        node, query)[0].score for query in nodes_to_align]
    max_score: float = max(scores)
    return [s/max_score for s in scores]


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
        graphs: list = [init_simple_graph(
            gfa, gfa_versions[i], colors[i]) for i, gfa in enumerate(gfa_files)]
    else:
        graphs: list = [init_graph(gfa, gfa_versions[i], colors[i])
                        for i, gfa in enumerate(gfa_files)]
    combined_view: MultiDiGraph = compose_all(graphs)
    graphs: list = [list(graph.nodes(data=True)) for graph in graphs]
    for idx, connex_component in enumerate(graphs):
        print(
            f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Started working on connex component #{idx}")
        all_other_nodes: list = list(
            chain.from_iterable(graphs[:idx] + graphs[idx+1:]))
        nodes_names: list = [n for n, _ in all_other_nodes]
        all_other_nodes = [data['seq'] for _, data in all_other_nodes]
        for node in connex_component:
            print(
                f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Aligning {node[0]} to graph...")
            scores: list = node_aligner(node[1]['seq'], all_other_nodes)
            for i, score in enumerate(scores):
                combined_view.add_edge(
                    node[0], nodes_names[i], color='forestgreen', arrows='', alpha=score)
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
                        seq=line.split()[2],
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
    print(f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Initializing graph for {gfa_file}")
    graph = MultiDiGraph()

    cmap: dict = {'SeqBt1': 'darkred',
                  'BtChar.1': 'indianred', 'BtChar.2': 'coral'}

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
                        seq=line.split()[2],
                        color=color
                    )
                case LineType.WALK:
                    if not gfa_line.line.idf == '_MINIGRAPH_':
                        add_path(
                            graph,
                            [f"{gfa_file.split('/')[-1].split('.')[0]}_{node}" for (
                                node, _) in gfa_line.line.path],
                            title=gfa_line.line.name,
                            color=cmap[gfa_line.line.name],
                            weight=4
                        )
                case LineType.PATH:
                    add_path(
                        graph,
                        [f"{gfa_file.split('/')[-1].split('.')[0]}_{node}" for (
                            node, _) in gfa_line.line.path],
                        title=gfa_line.line.name,
                        color=cmap[gfa_line.line.name],
                        weight=4
                    )
    graph.remove_nodes_from(list(isolates(graph)))
    print(f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Graph for {gfa_file} created!")
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

    print(f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Started nodes_align.py")
    display_graph(show_identity(args.file, args.gfa_version,
                  ['rebeccapurple', 'crimson', 'orchid']))
    print(
        f"[{datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}] Script ended sucessfully!")
