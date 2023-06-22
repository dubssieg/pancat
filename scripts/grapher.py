"Creates a graph we can navigate in."
from os import path, remove
from argparse import ArgumentParser, SUPPRESS
from networkx import MultiDiGraph
from pyvis.network import Network
from gfagraphs import Graph


def display_graph(graph: MultiDiGraph, name: str, colors_paths: dict[str, str], annotations: dict, output_path: str) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
        name (str): output name for graph render
        colors_paths (dict[str, str]): a set of colors to keep path colors consistent
    """
    graph_visualizer = Network(
        height='1000px', width='100%', directed=True, select_menu=False, filter_menu=False, bgcolor='#ffffff')
    graph_visualizer.set_template_dir(path.dirname(__file__), 'template.html')
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    html = graph_visualizer.generate_html()
    legend: str = '\n'.join(
        [f"<li><span class='{key}'></span> <a href='#'>{key}</a></li>" for key in colors_paths.keys()])
    with open(path.join(output_path, f"{name}_tmp.html"), "w+", encoding='utf-8') as out:
        out.write(html)
    with open(path.join(output_path, f"{name}.html"), "w", encoding="utf-8") as html_writer:
        with open(path.join(output_path, f"{name}_tmp.html"), "r", encoding="utf-8") as html_file:
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
    if path.exists(f"{name}_tmp.html"):
        remove(f"{name}_tmp.html")


def compute_stats(graph: Graph) -> dict:
    """Computes some basic metrics for the graph

    Args:
        graph (Graph): a gfagraphs Graph object

    Returns:
        dict: a container for metrics
    """

    stats: dict = {}
    stats["Number of segments"] = len(graph.segments)
    stats["Number of edges"] = len(graph.lines)

    nb_seg_sub50: int = 0
    nb_seg_sup50: int = 0
    size_seg_sub50: int = 0
    size_seg_sup50: int = 0

    for seg in graph.segments:
        if (size := len(seg.datas['seq'])) < 50:
            nb_seg_sub50 += 1
            size_seg_sub50 += size
        else:
            nb_seg_sup50 += 1
            size_seg_sup50 += size

    stats["Number of segments >=50bp"] = nb_seg_sup50
    stats["Number of segments <50bp"] = nb_seg_sub50
    stats["Length of segments >=50bp"] = size_seg_sup50
    stats["Length of segments <50bp"] = size_seg_sub50

    return stats


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument("-j", "--job_name", type=str, default="pangenome_graph",
                        help="Specifies a job name (ex : chr3_graph)")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Creates a representation of a pangenome graph.')
    parser.add_argument("output", type=str,
                        help="Path where to output the graph")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    args = parser.parse_args()

    pangenome_graph: MultiDiGraph = (pgraph := Graph(
        args.file, args.gfa_version, with_sequence=True)).compute_networkx()

    graph_stats = compute_stats(pgraph)

    display_graph(
        graph=pangenome_graph,
        name=args.job_name,
        colors_paths=pgraph.colors,
        annotations=graph_stats,
        output_path=args.output
    )
