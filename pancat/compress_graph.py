"""
This WIP algorithm is designed to losselesly compress graph in order to reduce the number of nodes, by:
    - compressiong 1-sized substitutions bubbles
    - tba
"""
from pgGraphs import Graph
from pgGraphs.abstractions import Orientation
from itertools import pairwise
from collections import Counter
from tharospytools.bio_tools import revcomp


def load_graph(gfa_file: str) -> Graph:
    """Inits a gfagraphs object with all the supplementary info computed.

    Args:
        gfa_file (str): path to a gfa file to reduce

    Returns:
        Graph: a graph with pre-computed info
    """
    graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True,
        low_memory=False,
        regexp=".*",
    )

    # Adding fully connected source and sinks to help with bubble detection
    graph.add_node(
        name='source',
        sequence=''
    )
    graph.add_node(
        name='sink',
        sequence=''
    )

    for path_data in graph.paths.values():
        graph.add_edge(
            'source',
            '+',
            path_data['path'][0][0],
            path_data['path'][0][1]
        )
        graph.add_edge(
            path_data['path'][-1][0],
            path_data['path'][-1][1],
            'sink',
            '+'
        )
        path_data['path'] = [
            ('source', Orientation('+')),
            *path_data['path'],
            ('sink', Orientation('+'))
        ]

    graph.compute_neighbors()
    graph.sequence_offsets()

    return graph


def get_substitutor(graph: Graph, successors_set: set, base_orientation: str, orientation_view: dict) -> str:
    """Returs the appropriate substitutor for the sequence of nucleotides to be abstracted

    Args:
        graph (Graph): a gfa graph
        successors_set (set): the nodes that shall be replaced
        base_orientation (str): orientation of the source node
        orientation_view (dict): a orientation mapping that collects info across all paths

    Returns:
        str: the substitutor string
    """
    x = set([graph.segments[successor]['seq'] if orientation_view[successor] == base_orientation else revcomp(
        graph.segments[successor]['seq']) for successor in successors_set])
    if x == {'G', 'A'}:
        substitutor: str = 'R'
    elif x == {'C', 'T'}:
        substitutor: str = 'Y'
    elif x == {'G', 'T'}:
        substitutor: str = 'K'
    elif x == {'C', 'A'}:
        substitutor: str = 'M'
    elif x == {'G', 'C'}:
        substitutor: str = 'S'
    elif x == {'T', 'A'}:
        substitutor: str = 'W'
    elif x == {'G', 'T', 'C'}:
        substitutor: str = 'B'
    elif x == {'G', 'T', 'A'}:
        substitutor: str = 'D'
    elif x == {'A', 'T', 'C'}:
        substitutor: str = 'H'
    elif x == {'G', 'A', 'C'}:
        substitutor: str = 'V'
    elif x == {'A', 'T', 'C', 'G'}:
        substitutor: str = 'N'
    else:
        try:
            (substitutor,) = x
        except ValueError:
            print(f"Error at x={x}.")
    return substitutor


def get_removable_bubbles(graph: Graph, orientation_view: dict) -> list[tuple[str, set[str], str]]:
    """Returns the list of bubbles that can be removed (that are 1-sized substitution bubbles)

    Args:
        graph (Graph): a loaded gfa graph
        orientation_view (dict): a mapping between nodes and absolute reading direction (+/-/?)

    Returns:
        list[tuple[str,set[str],str]]: a list of tuples, each tuple is a bubble, and each tuple contains source, sink and nodes inside
    """
    return [
        (
            source,
            node_data['successors'],
            list(sink)[0]
        ) for source, node_data in graph.segments.items() if all(
            [
                graph.segments[successor]['length'] ==
                1 for successor in node_data['successors']
            ]
        )
        and len(node_data['successors']) > 1
        and len(
            (sink := set().union(*[
                graph.segments[successor]['successors'] for successor in node_data['successors']
            ]
            ))
        ) == 1
        and len(
            set().union(
                *[
                    graph.segments[successor]['predecessors']
                    for successor in node_data['successors']
                ]
            )
        ) == 1
        and not orientation_view[source] == '?'
        and not orientation_view[list(sink)[0]] == '?'
        and not any(
            [
                orientation_view[successor] ==
                '?' for successor in node_data['successors']
            ]
        )
        and source not in ['source', 'sink']
        and sink not in ['source', 'sink']
        and set(graph.segments[source]['PO'].keys()) == set(graph.segments[list(sink)[0]]['PO'].keys()) == set().union(
            *[
                set(x) for x in [
                    list(graph.segments[successor]['PO'].keys()) for successor in node_data['successors']
                ]
            ]
        )
    ]


def get_orientation_view(graph: Graph) -> dict:
    """Returns a dict which maps nodes to their absolute direction.
    + means all paths read this node in forward
    - means all paths are in reverse
    ? means a bit of both

    Args:
        graph (Graph): _description_

    Returns:
        dict: _description_
    """
    orientation_view: dict[str, str] = dict()

    for path_data in graph.paths.values():
        for node, orientation in path_data['path']:
            orientation_view[node] = orientation.value if graph.segments.get(
                'orientation_view', orientation.value) == orientation.value else '?'

    return orientation_view


def compress_graph(gfa_file: str, gfa_output: str, minimized: bool = False) -> None:
    """Main call for graph compression

    Args:
        gfa_file (str): _description_
        gfa_output (str): _description_
        minimized (bool, optional): _description_. Defaults to False.
    """
    removed: dict[str, tuple] = dict()
    print(f"Loading file {gfa_file}")
    graph: Graph = load_graph(
        gfa_file=gfa_file
    )
    orientation_view: dict = get_orientation_view(
        graph=graph
    )
    removable_nodes: list[tuple[str, set[str], str]] = sorted(
        get_removable_bubbles(
            graph=graph,
            orientation_view=orientation_view
        ),
        key=lambda x: int(x[0])
    )

    # Safety for bubble chains
    mappings_sources: dict[str, str] = dict()

    # Safety for multi-bubbles branching at the same node
    counts_duplicates: list = [
        x for x, v in Counter(
            [
                s for s, _, _ in removable_nodes
            ]).items() if v > 1
    ] + [
        x for x, v in Counter(
            [
                s for _, _, s in removable_nodes
            ]).items() if v > 1
    ]
    # elements_in_bubbles: set = set().union(*[s for _, s, _ in removable_nodes])

    # Safe storing of sequences
    store_sequences: dict = dict()
    # Safety for keeping path information
    for seg_name, seg_data in graph.segments.items():
        store_sequences[seg_name] = seg_data['seq']
        seg_data['alternate_paths'] = dict()

    total_node_number: int = len(graph.segments)

    removed_node_count: int = 0

    for source, successors, sink in removable_nodes:
        if not source in counts_duplicates and not sink in counts_duplicates:
            real_source: str = mappings_sources.get(source, source)
            while real_source not in graph.segments:
                real_source = mappings_sources[real_source]

            graph.segments[real_source]['seq'] = graph.segments[real_source]['seq'] + get_substitutor(
                graph=graph,
                successors_set=successors,
                base_orientation=orientation_view[source],
                orientation_view=orientation_view
            ) + (
                graph.segments[sink]['seq'] if orientation_view[sink] == orientation_view[source] else revcomp(graph.segments[sink]['seq']))
            graph.segments[real_source]['length'] = len(
                graph.segments[real_source]['seq'])

            for path_name in graph.segments[real_source]['PO'].keys():
                try:
                    graph.segments[real_source]['alternate_paths'][path_name] = graph.segments[real_source]['alternate_paths'].get(
                        path_name,
                        [
                            (
                                real_source,
                                store_sequences[real_source],
                                orientation_view[real_source]
                            )
                        ]
                    ) + [
                        (
                            subpath_node,
                            store_sequences[subpath_node],
                            orientation_view[subpath_node]
                        ) for subpath_node in [
                            [
                                potential for potential in successors if path_name in graph.segments[potential]['PO'].keys()
                            ][0],
                            sink
                        ]
                    ]
                except IndexError:
                    # Edgecase where we reach the end of a path at the start of a bubble
                    pass

            for node in successors.union(set([sink])):
                mappings_sources[sink] = real_source
                try:
                    del graph.segments[node]
                    removed[node] = (source, sink)
                    removed_node_count += 1
                except KeyError:
                    pass

    print(f"Removed {removed_node_count} nodes form the graph ({round(removed_node_count/total_node_number*100,ndigits=2)}%)")

    # Cleaning the graph
    del graph.segments['source']
    del graph.segments['sink']

    for path_data in graph.paths.values():
        path_data['path'] = [
            (node, ori) for node, ori in path_data['path'] if node in graph.segments]
    """
    mark_for_destruction: list[tuple] = list()
    for doves in graph.lines.keys():
        if any([dove not in graph.segments or dove in ['source', 'sink'] for dove in doves]):
            mark_for_destruction.append(doves)
    for marked in mark_for_destruction:
        del graph.lines[marked]
    """

    for seg_data in graph.segments.values():
        del seg_data['predecessors']
        del seg_data['successors']

    # Cleaning edges (radical, but should work in theory)
    graph.lines = {}

    added_edges_count: int = 0
    # Validate all edges are according to paths
    for path_data in graph.paths.values():
        for (xi, yi), (xj, yj) in pairwise(path_data['path']):
            if (xi, xj) not in graph.lines:
                added_edges_count += 1
                graph.lines[(xi, xj)] = {
                    "start": xi,
                    "end": xj,
                    "orientation": f"{yi.value}/{yj.value}"
                }
    print(f"Added {added_edges_count} edges to the graph.")

    for seg_data in graph.segments.values():
        if seg_data['alternate_paths'] == {}:
            del seg_data['alternate_paths']

    if not minimized:
        graph.sequence_offsets(recalculate=True)
    graph.save_graph(gfa_output, minimal=minimized)
