"Creates a graph we can navigate in."
from os import path, remove
from json import load
from networkx import MultiDiGraph, compose
from pyvis.network import Network
from pgGraphs import Graph
from tharospytools.path_tools import path_allocator
from pgGraphs import GFANetwork, Graph as pgGraph
from tharospytools.list_tools import grouper


def display_graph(graph: MultiDiGraph, colors_paths: dict[str, str], annotations: dict, output_path: str) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
        name (str): output name for graph render
        colors_paths (dict[str, str]): a set of colors to keep path colors consistent
    """
    output_path: str = path_allocator(
        output_path, particle='.html', default_name='graph')
    output_path_temp: str = path_allocator(
        output_path, particle='.tmp.html', default_name='graph')
    graph_visualizer = Network(
        height='1000px', width='100%', directed=True, select_menu=False, filter_menu=False, bgcolor='#ffffff')
    graph_visualizer.set_template_dir(path.dirname(__file__), 'template.html')
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    html = graph_visualizer.generate_html()
    legend: str = '\n'.join(
        [f"<li><span class='{key}'></span> <a href='#'>{key}</a></li>" for key in colors_paths.keys()])
    with open(output_path_temp, "w+", encoding='utf-8') as out:
        out.write(html)
    with open(output_path, "w", encoding="utf-8") as html_writer:
        with open(output_path_temp, "r", encoding="utf-8") as html_file:
            for line in html_file:
                if "<div class='sidenav'>" in line:
                    html_writer.write(
                        f"""{line}{''.join(["<a href='#' title=''>"+str(key)+" : <b>"+str(value)+"</b></a>" for key,value in annotations.items()])}\n<ul class='legend'>{legend}</ul>"""
                    )
                elif "/* your colors */" in line:
                    html_writer.write(''.join(
                        [".legend ."+key+" { background-color: "+val+"; }" for key, val in colors_paths.items()]))
                else:
                    html_writer.write(line)
    if path.exists(output_path_temp):
        remove(output_path_temp)


def compute_stats(
    graph: Graph,
    length_classes: tuple[list] = (
        [0, 1], [2, 10], [11, 50], [51, 200], [201, 500], [
            501, 1000], [1001, 10000], [10001, float('inf')]
    )
) -> dict:
    """Computes some basic metrics for the graph

    Args:
        graph (Graph): a gfagraphs Graph object

    Returns:
        dict: a container for metrics
    """

    stats: dict = {}
    stats["Number of segments"] = len(graph.segments)
    stats["Number of edges"] = len(graph.lines)

    # segment_sizes[size_class] = [number,cum_length]
    segment_sizes: dict = {x: [0, 0] for x in range(len(length_classes))}
    total_size: int = 0
    for seg_datas in graph.segments.values():
        for x in segment_sizes.keys():
            low_bound, high_bound = length_classes[x]
            if seg_datas['length'] >= low_bound and seg_datas['length'] <= high_bound:
                segment_sizes[x] = [segment_sizes[x]
                                    [0]+1, segment_sizes[x][1]+seg_datas['length']]
                total_size += seg_datas['length']
    for x, (number, cum_size) in segment_sizes.items():
        low_bound, high_bound = length_classes[x]
        stats[f'Number of {low_bound} bp - {high_bound} bp'] = f"{number} ({round((number/len(graph.segments))*100,ndigits=2)}%)"
        stats[f'Size of {low_bound} bp - {high_bound} bp'] = f"{cum_size} ({round((cum_size/total_size)*100,ndigits=2)}%)"
    stats["Total size of segments"] = total_size
    try:
        stats["Is graph acyclic"] = all([len(set([x for x, _ in path_datas["path"]])) == len(
            path_datas["path"]) for path_datas in graph.paths.values])
    except:
        pass
    return stats


def multigraph_viewer(
    file_A: str,
    file_B: str,
    file_editions: str,
    boundaries: list,
    output: str,
) -> None:
    """Renders a digraph alignment and explicits their connexions inbetween them.

    Args:
        file_A (str): first gfa file
        file_B (str): second gfa file
        file_editions (str): json file containing editions, must be graph-level
        boundaries (list): list of node class sizes
        output (str): output for html file
    """
    bounds = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )

    # Creating pgGraphs object
    gfa_graph_A: pgGraph = pgGraph(
        gfa_file=file_A,
        with_sequence=True
    )
    gfa_graph_B: pgGraph = pgGraph(
        gfa_file=file_B,
        with_sequence=True
    )

    # Computing NetworkX structure
    pangenome_graph_A: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph_A,
        node_prefix='A',
        node_size_classes=bounds
    )
    pangenome_graph_B: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph_B,
        node_prefix='B',
        node_size_classes=bounds
    )

    full_graph: MultiDiGraph = compose(
        pangenome_graph_A,
        pangenome_graph_B
    )

    editions: dict = load(open(file_editions, 'r', encoding='utf-8'))

    for path, edition in editions.items():
        current_counter_A: int = 0
        current_counter_B: int = 0
        node_index_A: int = 0
        node_index_B: int = 0
        for style in ['merges', 'splits']:
            for pos, [node] in edition[style]:
                while current_counter_A < pos:
                    current_counter_A += gfa_graph_A.segments[gfa_graph_A.paths[path]
                                                              ['path'][node_index_A][0]]['length']
                    node_index_A += 1
                while current_counter_B < pos:
                    current_counter_B += gfa_graph_B.segments[gfa_graph_B.paths[path]
                                                              ['path'][node_index_B][0]]['length']
                    node_index_B += 1
                if not full_graph.has_edge(
                    f"A_{gfa_graph_A.paths[path]['path'][node_index_A-1][0]}",
                    f"B_{gfa_graph_B.paths[path]['path'][node_index_B-1][0]}"
                ):
                    full_graph.add_edge(
                        f"A_{gfa_graph_A.paths[path]['path'][node_index_A-1][0]}",
                        f"B_{gfa_graph_B.paths[path]['path'][node_index_B-1][0]}",
                        arrows='',
                        alpha=.5,
                        color='red' if style == 'merges' else 'blue'
                    )
    stats_A = compute_stats(
        graph=gfa_graph_A,
        length_classes=tuple(bounds)
    )
    stats_B = compute_stats(
        graph=gfa_graph_B,
        length_classes=tuple(bounds)
    )

    # Computing stats and displaying
    display_graph(
        graph=full_graph,
        colors_paths={
            **gfa_graph_A.metadata['colors'], **gfa_graph_B.metadata['colors']},
        annotations={
            **{f'A_{key}': val for key, val in stats_A.items()},
            **{f'B_{key}': val for key, val in stats_B.items()}
        },
        output_path=output
    )


def graph_viewer(
    file: str,
    output: str,
    boundaries: list
) -> None:
    bounds = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )

    # Creating pgGraphs object
    gfa_graph: pgGraph = pgGraph(
        gfa_file=file,
        with_sequence=True
    )

    # Computing NetworkX structure
    pangenome_graph: MultiDiGraph = GFANetwork.compute_networkx(
        graph=gfa_graph,
        node_prefix=None,
        node_size_classes=bounds
    )

    # Computing stats and displaying
    stats = compute_stats(
        graph=gfa_graph,
        length_classes=tuple(bounds)
    )
    display_graph(
        graph=pangenome_graph,
        colors_paths=gfa_graph.metadata['colors'],
        annotations=stats,
        output_path=output
    )


def graph_stats(file: str, boundaries: list) -> None:
    bounds = grouper(
        [0] + [bound+x for bound in boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
    )

    # Creating pgGraphs object
    gfa_graph: pgGraph = pgGraph(
        gfa_file=file,
        with_sequence=True
    )

    # Computing stats and displaying
    graph_stats = compute_stats(
        graph=gfa_graph,
        length_classes=tuple(bounds)
    )
    for key, value in graph_stats.items():
        print(f"{key}: {value}")
