"Does position-based checks of segment status between graphs, following paths."
from argparse import ArgumentParser, SUPPRESS
from typing import Generator
from itertools import combinations
from networkx import MultiDiGraph, compose_all
from pyvis.network import Network
from gfagraphs import Graph


def get_backbone(files: list, versions: list) -> Generator:
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
        for file, ver in duet:
            graph: Graph = Graph(file, ver)
            node_map.append({node.datas["name"]: node.datas['PO']
                            for node in graph.segments})
            path_map.append(
                {path.datas["name"]: path.datas["path"] for path in graph.get_path_list()})
            graph_map.append(
                {(x := file.split('.')[0].split('/')[-1]): graph.compute_networkx(x)})
            del graph
        yield (path_map, node_map, graph_map)


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
        "Shift": 0,
        "Inclusion A -> B": 0,
        "Inclusion B -> A": 0,
        "Reverse Inclusion A -> B": 0,
        "Reverse Inclusion B -> A": 0,
        "Equivalence": 0,
        "Reverse Equivalence": 0
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

            if start_a == start_b and end_a == end_b:
                # A & B are the same
                if ori_a == ori_b:
                    events["Equivalence"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='blueviolet', arrows='', label="E", weight=0.5)
                else:
                    events["Reverse Equivalence"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='blueviolet', arrows='', label="RE", weight=0.5)
                a += 1
                b += 1
            elif start_a >= start_b and end_a <= end_b:
                # Inclusion of A in B
                if ori_a == ori_b:
                    events["Inclusion A -> B"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkorange', label="I", weight=0.5)
                else:
                    events["Reverse Inclusion A -> B"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkorange', label="RI", weight=0.5)
                a += 1
            elif start_b >= start_a and end_b <= end_a:
                # Inclusion of B in A
                if ori_a == ori_b:
                    events["Inclusion B -> A"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkorange', label="I", weight=0.5)
                else:
                    events["Reverse Inclusion B -> A"] += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkorange', label="RI", weight=0.5)
                b += 1
            elif shifts:
                if start_a >= start_b and end_a >= end_b:
                    # Shift of A after B
                    events["Shift"] += 1
                    b += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_b, name_a, color='darkgreen', label="S", weight=0.5)
                elif start_a <= start_b and end_a <= end_b:
                    # Shift of B after A
                    events["Shift"] += 1
                    a += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        combined_view.add_edge(
                            name_a, name_b, color='darkgreen', label="S", weight=0.5)
    return (events, combined_view)


def display_graph(graph: MultiDiGraph, name: str) -> None:
    """Creates a interactive .html file representing the given graph

    Args:
        graph (MultiDiGraph): a graph combining multiple pangenomes to highlight thier similarities
    """
    graph_visualizer = Network(
        height='1000px', width='100%', directed=True)
    graph_visualizer.toggle_physics(True)
    graph_visualizer.from_nx(graph)
    graph_visualizer.set_edge_smooth('dynamic')
    try:
        graph_visualizer.show(f"{name}.html")
    except FileNotFoundError:
        pass


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
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

    for i, (path, nodes, graphs) in enumerate(get_backbone(args.file, args.gfa_version)):
        datas, full_graph = compare_positions(path, nodes, graphs, args.paths)
        display_graph(full_graph, f"offset_{i}")
