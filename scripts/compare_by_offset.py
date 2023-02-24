"Does position-based checks of segment status between graphs, following paths."
from argparse import ArgumentParser, SUPPRESS
from os import path, remove, getcwd
from typing import Generator
from itertools import combinations
from networkx import MultiDiGraph, compose_all
from pyvis.network import Network
from gfagraphs import Graph


def get_backbone(files: list, versions: list, with_sequences: bool = False) -> Generator:
    """Iterates through pairs of files, computes their graphs to extract nodes positions and paths

    Args:
        files (list): a list of files
        versions (list): a list of versions, ordered as the files

    Yields:
        Generator: tuple (paths,node data)
    """
    assert len(files) == len(versions), "Length of arguments does not match."
    for duet in combinations(zip(files, versions), 2):
        node_map: list = list()
        path_map: list = list()
        graph_map: list = list()
        colors: dict = dict()
        for file, ver in duet:
            graph: Graph = Graph(file, ver, with_sequences)
            node_map.append({node.datas["name"]: node.datas['PO']
                            for node in graph.segments})
            path_map.append(
                {path.datas["name"]: path.datas["path"] for path in graph.get_path_list()})
            graph_map.append(
                {(x := file.split('.')[0].split('/')[-1]): graph.compute_networkx(x)})
            colors: dict = {**graph.colors, **colors}
            del graph
        yield (path_map, node_map, graph_map, colors)


def compare_positions(paths: list[dict], nodes_datas: list[dict], graphs: list[dict], path_names: list | None = None, shifts: bool = True) -> tuple:
    """_summary_

    Args:
        paths (list): _description_
        nodes_datas (list): _description_

    Returns:
        dict: _description_
    """
    assert len(paths) == 2 or len(
        nodes_datas) == 2 or len(graphs) == 2, "Graph comparaison can only be made 2 at a time"
    path_of_first, path_of_second = paths
    if isinstance(path_names, list):
        path_of_first = {key: [(n, o) for (n, o) in val if n in nodes_datas[0].keys()] for key,
                         val in path_of_first.items() if key in path_names}
        path_of_second = {key: [(n, o) for (n, o) in val if n in nodes_datas[1].keys()] for key,
                          val in path_of_second.items() if key in path_names}
    assert len(set(path_of_first.keys())-set(path_of_second.keys())
               ) == 0, f"Your graphs does not use the same names for their paths. ({set(path_of_first.keys())-set(path_of_second.keys())} were unique)."
    nodes_of_a, nodes_of_b = nodes_datas
    combined_view: MultiDiGraph = compose_all(
        [g for graph in graphs for _, g in graph.items()])
    events: dict = {
        "Shift A -> B": {'number': 0, 'desc': 'Overlaps between A and B, where suffix of A is prefix of B'},
        "Shift B -> A": {'number': 0, 'desc': 'Overlaps between B and A, where suffix of B is prefix of A'},
        "Inclusion A -> B": {'number': 0, 'desc': 'Inclusions between A and B, where A is a subsequence of B'},
        "Inclusion B -> A": {'number': 0, 'desc': 'Inclusions between B and A, where B is a subsequence of A'},
        "Reverse Inclusion A -> B": {'number': 0, 'desc': 'Inclusions between A and B, where A is a subsequence of B, and B is reversed.'},
        "Reverse Inclusion B -> A": {'number': 0, 'desc': 'Inclusions between B and A, where B is a subsequence of A, and B is reversed.'},
        "Equivalence": {'number': 0, 'desc': 'Equivalences between A and B, where A is equal to B.'},
        "Reverse Equivalence": {'number': 0, 'desc': 'Equivalences between A and B, where A is equal to B, and B is reversed.'},
    }
    name_graph_a, name_graph_b = list(graphs[0].keys())[
        0], list(graphs[1].keys())[0]
    for name, path_a in path_of_first.items():
        path_b = path_of_second[name]
        a, b = 0, 0
        while a < len(path_a) and b < len(path_b):

            start_a, end_a, ori_a = nodes_of_a[path_a[a][0]][name]
            start_b, end_b, ori_b = nodes_of_b[path_b[b][0]][name]
            name_a = f"{name_graph_a}_{path_a[a][0]}"
            name_b = f"{name_graph_b}_{path_b[b][0]}"

            # print(f'{name_a}: {start_a}/{end_a} ; {name_b}: {start_b}/{end_b}')

            if start_a == start_b and end_a == end_b:
                # A & B are the same
                if ori_a == ori_b:
                    events["Equivalence"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='blueviolet', arrows='', label="E", weight=0.5)
                else:
                    events["Reverse Equivalence"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='blueviolet', arrows='', label="RE", weight=0.5)
                a += 1
                b += 1
            elif start_a >= start_b and end_a <= end_b:
                # Inclusion of A in B
                if ori_a == ori_b:
                    events["Inclusion A -> B"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkorange', label="I", weight=0.5)
                else:
                    events["Reverse Inclusion A -> B"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkorange', label="RI", weight=0.5)
                a += 1
            elif start_b >= start_a and end_b <= end_a:
                # Inclusion of B in A
                if ori_a == ori_b:
                    events["Inclusion B -> A"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkorange', label="I", weight=0.5)
                else:
                    events["Reverse Inclusion B -> A"]['number'] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkorange', label="RI", weight=0.5)
                b += 1
            elif shifts:
                if start_a < start_b and end_a < end_b and end_a != start_b:
                    # Shift of B after A
                    events["Shift A -> B"]['number'] += 1
                    a += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkgreen', label="S", weight=0.5)
                elif start_b < start_a and end_b < end_a and end_b != start_a:
                    # Shift of A after B
                    events["Shift B -> A"]['number'] += 1
                    b += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkgreen', label="S", weight=0.5)
                else:
                    if end_a == start_b:
                        a += 1
                    else:
                        b += 1
    return (events, combined_view)


def display_graph(graph: MultiDiGraph, name: str, paths: list, datas_events: dict, filenames: list[str], colors_paths: dict[str, str]) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
    """
    graph_visualizer = Network(
        height='1000px', width='100%', directed=True, select_menu=False, filter_menu=False, bgcolor='#ffffff')
    graph_visualizer.set_template_dir(path.dirname(__file__), 'template.html')
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    # graph_visualizer.show(f"{name}_tmp.html")
    html = graph_visualizer.generate_html()
    with open(f"{name}_tmp.html", "w+", encoding='utf-8') as out:
        out.write(html)
    # <img src='{gfa_file.split('.')[0].split('/')[-1]}_cbar.png' align='center' rotate='90'>
    # with open(f"{name}_tmp.html", "r", encoding="utf-8") as html_reader:
    sidebar: str = '\n'.join(['<a href=\'#\' title=\''+val['desc']+'\'>'+key +
                             ' : '+str(val['number'])+'</a>' for key, val in datas_events.items()])
    legend: str = '\n'.join(
        [f"<li><span class='{key}'></span> <a href='#'>{key}</a></li>" for key in colors_paths.keys()])
    with open(f"{name}.html", "w", encoding="utf-8") as html_writer:
        with open(f"{name}_tmp.html", "r", encoding="utf-8") as html_file:
            for line in html_file:
                if "<div class='sidenav'>" in line:
                    html_writer.write(
                        f"""{line}
<a href='#' title='Path(s) used for the comparaison between offsets'>Paths alignments : <b>{', '.join(paths)}</b></a>
<a href='#' title='First specified graph, considered as the default orientation'>Graph A : {filenames[0].split('/')[-1].split('.')[0]}</a>
<a href='#' title='Second specified graph, evaluated against the first'>Graph B : {filenames[1].split('/')[-1].split('.')[0]}</a>
{sidebar}\n<ul class='legend'>{legend}</ul>"""
                    )
                elif "/* your colors */" in line:
                    html_writer.write(''.join(
                        [".legend ."+key+" { background-color: "+val+"; }" for key, val in colors_paths.items()]))
                else:
                    html_writer.write(line)
    if path.exists(f"{name}_tmp.html"):
        remove(f"{name}_tmp.html")


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
    parser.add_argument("-j", "--job_name", type=str, required=True,
                        help="Job identifier for output (ex : chr3_graph)")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Does position-based checks of segment status between graphs, following paths.')
    parser.add_argument(
        "-p", "--paths", type=str, help="Path(s) to display.", nargs='+', default=None)
    args = parser.parse_args()

    # Handling input errors
    if len(args.file) < 2:
        parser.error("Please specify at least two GFA files as input.")
    if len(args.file) != len(args.gfa_version):
        parser.error(
            "Please match the number of args between files and gfa types.")

    for i, (path, nodes, graphs, colors) in enumerate(get_backbone(args.file, args.gfa_version)):
        datas, full_graph = compare_positions(path, nodes, graphs, args.paths)
        display_graph(full_graph, f"{args.job_name}_{i}", args.paths, datas, [
                      args.file[i], args.file[i+1]], colors)
