"Creates a graph we can navigate in."
from os import path, remove
from argparse import ArgumentParser, SUPPRESS
from networkx import MultiDiGraph
from pyvis.network import Network
from gfagraphs import Graph
from tharospytools import path_allocator


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


def compute_stats(graph: Graph, length_classes: tuple[list] = ([0, 1], [2, 10], [11, 50], [51, 200],
                                                               [201, 500], [501, 1000], [1001, 10000], [10001, float('inf')])) -> dict:
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
    for seg in graph.segments:
        size = len(seg.datas['seq'])
        for x in segment_sizes.keys():
            low_bound, high_bound = length_classes[x]
            if size >= low_bound and size <= high_bound:
                segment_sizes[x] = [segment_sizes[x]
                                    [0]+1, segment_sizes[x][1]+size]
                total_size += size
    for x, (number, cum_size) in segment_sizes.items():
        low_bound, high_bound = length_classes[x]
        stats[f'Number of {low_bound} bp - {high_bound} bp'] = f"{number} ({round((number/len(graph.segments))*100,ndigits=2)}%)"
        stats[f'Size of {low_bound} bp - {high_bound} bp'] = f"{cum_size} ({round((cum_size/total_size)*100,ndigits=2)}%)"
    stats["Total size of segments"] = total_size
    return stats


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Creates a representation of a pangenome graph.')
    parser.add_argument("output", type=str,
                        help="Output path for the html graph file.")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    args = parser.parse_args()

    pangenome_graph: MultiDiGraph = (pgraph := Graph(
        args.file, with_sequence=True)).compute_networkx()

    graph_stats = compute_stats(pgraph)

    display_graph(
        graph=pangenome_graph,
        colors_paths=pgraph.colors,
        annotations=graph_stats,
        output_path=args.output
    )
