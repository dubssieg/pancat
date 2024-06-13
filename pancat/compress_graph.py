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
    node_strs: list[str] = list(
        set(
            [
                graph.segments[successor]['seq'] if orientation_view[successor] ==
                base_orientation else revcomp(graph.segments[successor]['seq']) for successor in successors_set
            ]
        )
    )
    substitutor: str = ''
    for i in range(len(node_strs[0])):
        x: set[str] = set([strs[i] for strs in node_strs])
        if x == {'G', 'A'}:
            substitutor += 'R'
        elif x == {'C', 'T'}:
            substitutor += 'Y'
        elif x == {'G', 'T'}:
            substitutor += 'K'
        elif x == {'C', 'A'}:
            substitutor += 'M'
        elif x == {'G', 'C'}:
            substitutor += 'S'
        elif x == {'T', 'A'}:
            substitutor += 'W'
        elif x == {'G', 'T', 'C'}:
            substitutor += 'B'
        elif x == {'G', 'T', 'A'}:
            substitutor += 'D'
        elif x == {'A', 'T', 'C'}:
            substitutor += 'H'
        elif x == {'G', 'A', 'C'}:
            substitutor += 'V'
        elif x == {'A', 'T', 'C', 'G'}:
            substitutor += 'N'
        else:
            try:
                (a,) = x
                substitutor += a
            except ValueError:
                print(f"Error at x={x}.")
    return substitutor


def get_removable_bubbles(graph: Graph, orientation_view: dict, max_len_to_collapse: int) -> list[tuple[str, set[str], str]]:
    """Returns the list of bubbles that can be removed (that are 1-sized substitution bubbles)

    Args:
        graph (Graph): a loaded gfa graph
        orientation_view (dict): a mapping between nodes and absolute reading direction (+/-/?)

        /!\
            Notion of predecessors/successors might need to be refined as we only take into account one reading direction with this method.
            As of now it is not a problem as we only have + orientations and ? orientaions, because mixed can't be resolved the same way
            We will also have to take into account a better way to represent edges' orientation within such bubbles.
            Some edges are missing in some graphs (orientations showing in graph paths but not in edges)

    Returns:
        list[tuple[str,set[str],str]]: a list of tuples, each tuple is a bubble, and each tuple contains source, sink and nodes inside

    """
    return [
        (
            source,
            node_data['successors'],
            list(sink)[0]
        ) for source, node_data in graph.segments.items() if len(
            set(
                [
                    graph.segments[successor]['length'] for successor in node_data['successors']
                ]
            )
        ) == 1
        and all([graph.segments[successor]['length'] <= max_len_to_collapse for successor in node_data['successors']])
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
        and orientation_view.get(source, '?') != '?'
        and orientation_view.get(list(sink)[0], '?') != '?'
        and all(
            [
                orientation_view.get(successor, '?') !=
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
        graph (Graph): a gfagraphs object

    Returns:
        dict: a descriptor for nodes orientation for nodes that are in paths
    """
    orientation_view: dict[str, str] = dict()

    for path_data in graph.paths.values():
        for node, orientation in path_data['path']:
            orientation_view[node] = orientation.value if orientation_view.get(
                node, orientation.value) == orientation.value else '?'

    return orientation_view


def compress_graph(gfa_file: str, gfa_output: str, minimized: bool = False, max_len_to_collapse: int = float('inf')) -> None:
    """Main call for graph compression

    Args:
        gfa_file (str): a gfa standard formatted file
        gfa_output (str): a container for the output, in gfa format
        minimized (bool, optional): if the graph should contain the minimum number of informations to be a gfa file. Defaults to False.
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
            orientation_view=orientation_view,
            max_len_to_collapse=max_len_to_collapse,
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

    for seg_data in graph.segments.values():
        del seg_data['predecessors']
        del seg_data['successors']

    # Cleaning edges (radical, but should work in theory)
    graph.lines = {}

    # Validate all edges are according to paths
    for path_data in graph.paths.values():
        for (xi, yi), (xj, yj) in pairwise(path_data['path']):
            graph.add_edge(
                source=xi,
                ori_source=yi,
                sink=xj,
                ori_sink=yj,
            )

    for seg_data in graph.segments.values():
        if seg_data['alternate_paths'] == {}:
            del seg_data['alternate_paths']

    if not minimized:
        graph.sequence_offsets(recalculate=True)
    graph.save_graph(gfa_output, minimal=minimized)
