"Does position-based checks of segment status between graphs, following paths."
from argparse import ArgumentParser, SUPPRESS
from os import path, remove
from enum import Enum
from typing import Generator, Iterable
from itertools import combinations
from networkx import MultiDiGraph, compose_all, isolates
from pyvis.network import Network
import matplotlib.pyplot as plt
from gfagraphs import Graph


class Event(Enum):
    "Describes local event types between nodes"
    SHIFT = 'S'
    INCLUSION = 'I'
    EQUIVALENCE = 'E'
    REVERSE_SHIFT = 'RS'
    REVERSE_INCLUSION = 'RI'
    REVERSE_EQUIVALENCE = 'RE'


def event_score(
    event_type: Event,
    value: int,
    scores: dict = {
        Event.SHIFT: -2,
        Event.INCLUSION: 0,
        Event.EQUIVALENCE: 1,
        Event.REVERSE_SHIFT: -3,
        Event.REVERSE_INCLUSION: -1,
        Event.REVERSE_EQUIVALENCE: 0
    }
) -> int:
    """Calculates the score for the current event

    Args:
        event_type (Event): description of event
        value (int): number of BP that are concerned
        scores (dict, optional): Describes the score for one bp per event.

    Returns:
        int: the score associated to the event
    """
    return scores[event_type]*value


def get_backbone(files: list, versions: list, with_sequences: bool = True) -> Generator:
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


def compare_positions(
        output_file: str,
        paths: list[dict],
        nodes_datas: list[dict],
        graphs: list[dict],
        reference: str,
        path_names: list | None = None,
) -> tuple:
    """Seeks for events within the alignement of two graphs.

    Args:
        output_file (str): _description_
        paths (list[dict]): _description_
        nodes_datas (list[dict]): _description_
        graphs (list[dict]): _description_
        reference (str): _description_
        path_names (list | None, optional): _description_. Defaults to None.

    Returns:
        tuple: (
            events : a list of all events,
            combined_view : a graph that packs the two graphs and their relations,
            score : a result of a score function,equivalence
            equivalences : a list of relations for the reference path
                --> needs to be refinded as "a list of sinks and sources"
            )
    """
    assert len(paths) == 2 or len(
        nodes_datas) == 2 or len(graphs) == 2, "Graph comparaison can only be made 2 at a time"
    path_of_first, path_of_second = paths
    if isinstance(path_names, list):
        # a path or a list of paths is specified by the user
        path_of_first = {key: [(n, o) for (n, o) in val if n in nodes_datas[0].keys()] for key,
                         val in path_of_first.items() if key in path_names}
        path_of_second = {key: [(n, o) for (n, o) in val if n in nodes_datas[1].keys()] for key,
                          val in path_of_second.items() if key in path_names}
    else:
        # no path specified, considering all paths
        path_of_first = {key: [(n, o) for (n, o) in val if n in nodes_datas[0].keys()] for key,
                         val in path_of_first.items()}
        path_of_second = {key: [(n, o) for (n, o) in val if n in nodes_datas[1].keys()] for key,
                          val in path_of_second.items()}
    assert len(set(path_of_first.keys())-set(path_of_second.keys())
               ) == 0, f"Your graphs does not use the same names for their paths. ({set(path_of_first.keys())-set(path_of_second.keys())} were unique)."

    nodes_of_a, nodes_of_b = nodes_datas
    combined_view: MultiDiGraph = compose_all(
        [g for graph in graphs for _, g in graph.items()])
    events: dict[str, dict[str, str | int]] = {
        "Shift A -> B": {'number': 0, 'desc': 'Overlaps between A and B, where a suffix of A is a prefix of B'},
        "Shift B -> A": {'number': 0, 'desc': 'Overlaps between B and A, where a suffix of B is a prefix of A'},
        "Inclusion A -> B": {'number': 0, 'desc': 'Inclusions between A and B, where A is a subsequence of B'},
        "Inclusion B -> A": {'number': 0, 'desc': 'Inclusions between B and A, where B is a subsequence of A'},
        "Reverse Inclusion A -> B": {'number': 0, 'desc': 'Inclusions between A and B, where A is a subsequence of B, and B is reversed.'},
        "Reverse Inclusion B -> A": {'number': 0, 'desc': 'Inclusions between B and A, where B is a subsequence of A, and B is reversed.'},
        "Equivalence": {'number': 0, 'desc': 'Equivalences between A and B, where A is equal to B.'},
        "Reverse Equivalence": {'number': 0, 'desc': 'Equivalences between A and B, where A is equal to B, and B is reversed.'},
    }
    name_graph_a, name_graph_b = list(graphs[0].keys())[
        0], list(graphs[1].keys())[0]
    align_score: int = 0
    equivalence_list: list = list()
    number_of_paths: int = len(set(path_of_first.keys()))
    with open(f"{output_file}.tsv", 'w', encoding='utf-8') as report:
        header: str = "PATH\tSOURCE\tSSTART\tSSTOP\tTARGET\tTSTART\tTSTOP\tEVENT\tREVERSED\tSCORE\tAMBIGUOUS\n"
        report.write(header)
        for name, path_a in path_of_first.items():
            path_b = path_of_second[name]
            a, b = 0, 0
            while a < len(path_a) and b < len(path_b):
                try:
                    start_a, end_a, ori_a = nodes_of_a[path_a[a][0]][name]
                    start_b, end_b, ori_b = nodes_of_b[path_b[b][0]][name]
                    ambiguous: bool = False
                except KeyError:
                    # Path follows reference (might be ambiguous)
                    start_a, end_a, ori_a = nodes_of_a[path_a[a][0]][reference]
                    start_b, end_b, ori_b = nodes_of_b[path_b[b][0]][reference]
                    ambiguous: bool = True
                name_a = f"{name_graph_a}_{path_a[a][0]}"
                name_b = f"{name_graph_b}_{path_b[b][0]}"

                if start_a == start_b and end_a == end_b:
                    if name == reference:
                        equivalence_list.append((
                            name_a,
                            name_b,
                            combined_view.degree[name_a] +
                            combined_view.degree[name_b] == 4 *
                            number_of_paths,
                            Event.EQUIVALENCE
                        ))
                    # A & B are the same
                    if ori_a == ori_b:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Equivalence"]['number'] += 1
                            combined_view.add_edge(
                                name_a, name_b, color='blueviolet', arrows='', label="E", weight=0.5, title=(score := event_score(Event.EQUIVALENCE, abs(end_a-start_a))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tE\t{False}\t{score}\t{ambiguous}\n")
                    else:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Reverse Equivalence"]['number'] += 1
                            combined_view.add_edge(
                                name_a, name_b, color='blueviolet', arrows='', label="RE", weight=0.5, title=(score := event_score(Event.REVERSE_EQUIVALENCE, abs(end_a-start_a))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tE\t{True}\t{score}\t{ambiguous}\n")
                    a += 1
                    b += 1

                elif start_a >= start_b and end_a <= end_b:
                    if name == reference:
                        equivalence_list.append((
                            name_a,
                            name_b,
                            combined_view.degree[name_a] +
                            combined_view.degree[name_b] == 4 *
                            number_of_paths,
                            Event.INCLUSION
                        ))
                    # Inclusion of A in B
                    if ori_a == ori_b:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Inclusion A -> B"]['number'] += 1
                            combined_view.add_edge(
                                name_a, name_b, color='darkorange', label="I", weight=0.5, title=(score := event_score(Event.INCLUSION, abs(end_a-start_a))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tI\t{False}\t{score}\t{ambiguous}\n")
                    else:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Reverse Inclusion A -> B"]['number'] += 1
                            combined_view.add_edge(
                                name_a, name_b, color='darkorange', label="RI", weight=0.5, title=(score := event_score(Event.REVERSE_INCLUSION, abs(end_a-start_a))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tI\t{True}\t{score}\t{ambiguous}\n")
                    a += 1
                elif start_b >= start_a and end_b <= end_a:
                    if name == reference:
                        equivalence_list.append((
                            name_a,
                            name_b,
                            combined_view.degree[name_a] +
                            combined_view.degree[name_b] == 4 *
                            number_of_paths,
                            Event.INCLUSION
                        ))
                    # Inclusion of B in A
                    if ori_a == ori_b:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Inclusion B -> A"]['number'] += 1
                            combined_view.add_edge(
                                name_b, name_a, color='darkorange', label="I", weight=0.5, title=(score := event_score(Event.INCLUSION, abs(end_b-start_b))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tI\t{False}\t{score}\t{ambiguous}\n")
                    else:
                        if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                            events["Reverse Inclusion B -> A"]['number'] += 1
                            combined_view.add_edge(
                                name_b, name_a, color='darkorange', label="RI", weight=0.5, title=(score := event_score(Event.REVERSE_INCLUSION, abs(end_b-start_b))))
                            report.write(
                                f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tI\t{True}\t{score}\t{ambiguous}\n")
                    b += 1
                elif start_a < start_b and end_a < end_b and end_a != start_b and start_b < end_a:
                    if name == reference:
                        equivalence_list.append((
                            name_a,
                            name_b,
                            combined_view.degree[name_a] +
                            combined_view.degree[name_b] == 4 *
                            number_of_paths,
                            Event.SHIFT
                        ))
                    # Shift of B after A
                    a += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        events["Shift A -> B"]['number'] += 1
                        combined_view.add_edge(
                            name_a, name_b, color='darkgreen', label="O", weight=0.5, title=(score := event_score(Event.SHIFT, abs(end_a-start_b))))
                        report.write(
                            f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tS\t{False}\t{score}\t{ambiguous}\n")
                elif start_b < start_a and end_b < end_a and end_b != start_a and start_a < end_b:
                    if name == reference:
                        equivalence_list.append((
                            name_a,
                            name_b,
                            combined_view.degree[name_a] +
                            combined_view.degree[name_b] == 4 *
                            number_of_paths,
                            Event.SHIFT
                        ))
                    # Shift of A after B
                    b += 1
                    if not combined_view.has_edge(name_a, name_b) and not combined_view.has_edge(name_b, name_a):
                        events["Shift B -> A"]['number'] += 1
                        combined_view.add_edge(
                            name_b, name_a, color='darkgreen', label="O", weight=0.5, title=(score := event_score(Event.SHIFT, abs(end_b-start_a))))
                        report.write(
                            f"{name}\t{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\tS\t{False}\t{score}\t{ambiguous}\n")
                else:
                    score = 0  # no particular event
                    if end_a <= start_b:
                        # segment is before, or ends at the position of next to compare, we shift the comparison
                        a += 1
                    else:
                        # same, but for b
                        b += 1
                align_score += score
        report.write(f"Total align score: {align_score}\n")
        for event, data in events.items():
            report.write(f"## {event} = {data['number']} ({data['desc']})\n")
    return (events, combined_view, align_score, equivalence_list)


def extract_components(
        combined_graph: MultiDiGraph,
        equivalence_list: list
) -> MultiDiGraph:
    """Given a networkx graph, keeps zones that are preceded and succeded by a "equivalence" relation
    that implies all genomes to go through the equivalence

    Args:
        combined_graph (MultiDiGraph): a graph comparison

    Returns:
        MultiDiGraph: a set of connex components
    """
    sources_sinks: list = list()
    bounds_A: list = list()
    bounds_B: list = list()
    merge_source_sinks: list = list()
    idx: int = 0
    node_A_prefix = '_'.join(equivalence_list[0][0].split('_')[:-1])
    node_B_prefix = '_'.join(equivalence_list[0][1].split('_')[:-1])
    source_nodes: tuple = (equivalence_list[0][0], equivalence_list[0][1])
    in_component: bool = equivalence_list[0][3] != Event.EQUIVALENCE
    for event in equivalence_list:
        # status = True if all genomes go through node_a and node_b
        node_A, node_B, status, event_style = event
        if event_style == Event.EQUIVALENCE and status:
            if in_component:
                # we found a sink
                in_component = False
                idx += 1
                sources_sinks.extend(
                    [node_A, node_B, source_nodes[0], source_nodes[1]])
                sink_nodes: tuple = (node_A, node_B)
                upper_bound_A, lower_bound_A = int(node_A.split(
                    '_')[-1]), int(source_nodes[0].split('_')[-1])
                bounds_A.extend(boundedA := [f"{node_A_prefix}_{num_node}" for num_node in range(
                    lower_bound_A, upper_bound_A+1)])
                upper_bound_B, lower_bound_B = int(node_B.split(
                    '_')[-1]), int(source_nodes[1].split('_')[-1])
                bounds_B.extend(boundedB := [f"{node_B_prefix}_{num_node}" for num_node in range(
                    lower_bound_B, upper_bound_B+1)])
                merge_source_sinks.append((source_nodes, sink_nodes))

                draw_component(combined_graph, boundedA, boundedB, idx, list(set(item for sublist in [
                               combined_graph.nodes[node]['offsets'].keys() for node in boundedA+boundedB] for item in sublist)))
            else:
                # we found a potential source
                source_nodes: tuple = (node_A, node_B)
        elif event_style != Event.EQUIVALENCE:
            in_component = True
            sources_sinks.extend([node_A, node_B])
        elif event_style == Event.EQUIVALENCE and not status:
            sources_sinks.extend([node_A, node_B])
    nodes_to_remove: list = [
        node for node in combined_graph.nodes if node not in bounds_A and node not in bounds_B]
    combined_graph.remove_nodes_from(nodes_to_remove)
    # Cleaning graph
    deg = combined_graph.degree()
    to_remove = [n for (n, degree) in deg if degree == 1]
    combined_graph.remove_nodes_from(to_remove)
    combined_graph.remove_nodes_from(list(isolates(combined_graph)))
    for i in range(len(merge_source_sinks[:-1])):
        # to disable double sinks/sources from appearing
        if not merge_source_sinks[i][0] == merge_source_sinks[i+1][0]:
            sink_A, sink_B = merge_source_sinks[i][1]
            next_source_A, next_source_B = merge_source_sinks[i+1][0]
            combined_graph.add_edge(sink_A, next_source_A,
                                    color='grey', weight=10, arrows='', alpha=0.1)
            combined_graph.add_edge(sink_B, next_source_B,
                                    color='grey', weight=10, arrows='', alpha=0.1)
    return combined_graph


def draw_component(
        graph: MultiDiGraph,
        nodes_of_A: list,
        nodes_of_B: list,
        idx: int,
        alignments: list
) -> None:
    """Draws on a matplotlib plot each of the alignement for each of the pair of paths

    Args:
        graph (MultiDiGraph): diGraph from networkx
        nodes_of_A (list): a list of the nodes of the graph A
        nodes_of_B (list): a list of the nodes of the graph B
        idx (int): index of current component
        alignments (list) a list of all the paths within the subgraph
    """
    fig, axs = plt.subplots(len(alignments), figsize=(20, len(alignments)*2))
    for i, align in enumerate(alignments):
        axs[i].set_title(
            f"{align} on zone {idx}. (top: {'_'.join(nodes_of_A[0].split('_')[:-1])}, bottom: {'_'.join(nodes_of_B[0].split('_')[:-1])})")
        axs[i].axis('off')
        axs[i].set_ylim(-4, 4)

    for i, subgraph in [(0.5, nodes_of_A), (-0.5, nodes_of_B)]:
        for node_name in subgraph:
            offsets: dict = graph.nodes[node_name]['offsets']
            for wayline in alignments:
                if wayline in offsets.keys():
                    start, end, _ = offsets[wayline]
                    axs[alignments.index(wayline)].text(
                        (start+end)/2, i*3, node_name.split('_')[-1], fontsize=5, va='center', ha='center', rotation=90)
                    axs[alignments.index(wayline)].plot([start, end], [i, i])

    plt.savefig(f"component_{idx}.png", bbox_inches='tight')


def display_graph(graph: MultiDiGraph, name: str, paths: list, datas_events: dict, filenames: list[str], colors_paths: dict[str, str], score: int) -> None:
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
    html = graph_visualizer.generate_html()
    with open(f"{name}_tmp.html", "w+", encoding='utf-8') as out:
        out.write(html)
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
<a href='#' title='Path(s) used for the comparaison between offsets'>Paths alignments : <b>{', '.join(paths) if isinstance(paths,Iterable) else 'paired mode'}</b></a>
<a href='#' title='First specified graph, considered as the default orientation'>Graph A : {filenames[0].split('/')[-1].split('.')[0]}</a>
<a href='#' title='Second specified graph, evaluated against the first'>Graph B : {filenames[1].split('/')[-1].split('.')[0]}</a>
<a href='#' title='Total alignment score between the two graphs'>Score : {score}</a>
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
    parser.add_argument("-r", "--reference", type=str, required=True,
                        help="Path to refer to if position is ambiguous.")
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
        datas, full_graph, score, equivalence_list = compare_positions(
            f"{args.job_name}_{i}", path, nodes, graphs, args.reference, path_names=args.paths)
        display_graph(full_graph, f"{args.job_name}_{i}", args.paths, datas, [
                      args.file[i], args.file[i+1]], colors, score)
        purged_graph: MultiDiGraph = extract_components(
            full_graph, equivalence_list)
        display_graph(purged_graph, f"{args.job_name}_{i}_extracted", args.paths, datas, [
                      args.file[i], args.file[i+1]], colors, score)
