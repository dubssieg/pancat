from networkx import all_simple_paths
from argparse import ArgumentParser
from gfagraphs import Graph, GFANetwork, GFAParser
from json import dump
from typing import Any
from networkx import DiGraph


def load_graph(graph_path: str) -> Graph:
    """Load a graph in memory and computes its offsets

    Parameters
    ----------
    graph_path : str
        path to a valid gfa file

    Returns
    -------
    Graph
        annotated gfagraphs object
    """
    gfa: Graph = Graph(
        gfa_file=graph_path,
        with_sequence=True,
        with_reverse_edges=False,
        low_memory=False,
        regexp='.*',
    )
    gfa.compute_neighbors()
    gfa.sequence_offsets()
    return gfa


def find_anchors(gfa: Graph, interval: tuple[int, int], reference: str) -> tuple[str, str]:
    """Finds the pair of nodes that corresponds to start and end positions

    Parameters
    ----------
    gfa : Graph
        loaded and anotated gfa graph
    positions : tuple[int,int]
        a couple of (start,stop) positions
    path_name : str
        name of the path for positions

    Returns
    -------
    tuple[str,str]
        a tuple of node names
    """
    if interval[0] > interval[1]:
        interval = (interval[1], interval[0])
    if interval[0] < 0:
        interval = (0, interval[1])
    if interval[1] > (ln := sum(gfa.segments[x]['length'] for x, _ in gfa.paths[reference]['path'])) or interval[0] > interval[1]:
        interval = (interval[0], ln)

    encountered_nodes: dict[str, int] = {}
    start: bool | str = False
    stop: bool | str = False
    for node, _ in gfa.paths[reference]['path']:
        if gfa.segments[node]['PO'][reference][encountered_nodes.get(node, 0)][1] >= interval[0] and not start:
            start = node
            stac: int = encountered_nodes.get(node, 0)
        if gfa.segments[node]['PO'][reference][encountered_nodes.get(node, 0)][1] >= interval[1] and start and not stop:
            stop = node
            stoc: int = encountered_nodes.get(node, 0)
        encountered_nodes[node] = encountered_nodes.get(node, 0) + 1
    return ((start, stac), (stop, stoc))


def extract_subgraph(gfa_path: str, x: int, y: int, sequence: str, output: str) -> None:
    """Gathers the subgraph from two points on a reference

    Parameters
    ----------
    gfa_path : str
        path to file
    x : int
        start point
    y : int
        end point
    sequence : str
        reference genome for positions
    output : str
        path to output file
    """
    gfa_graph: Graph = load_graph(gfa_path)

    ((source, stac), (sink, stoc)) = find_anchors(
        gfa=gfa_graph,
        interval=(x, y),
        reference=sequence,
    )

    path_set: set[str] = set(gfa_graph.segments[source]['PO'].keys()).intersection(
        set(gfa_graph.segments[sink]['PO'].keys()))

    # nodes_between_set: set[str] = {node for path in all_simple_paths(nx_graph, source=source, target=sink) for node in path}
    selected_nodes_lists: list[str] = [(pl := [x for x, _ in gfa_graph.paths[path_name]['path']])[
        pl.index(source):pl.index(sink)+1] for path_name in path_set]

    nodes_between_set: set[str] = {
        node for nl in selected_nodes_lists for node in nl
    }

    # créer fonction d'extraction à partir de noeuds dans la lib
    GFAParser.save_subgraph(
        graph=gfa_graph,
        output_path=output,
        nodes=nodes_between_set,
        force_format=False,
        minimal_graph=True
    )
