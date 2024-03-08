"""
pancat correct is about to iterate through the graph, and fix small mistakes made by graph creation tools.
    - (TBI) correct missing edges in the graph according to paths taht are described in the file
"""
from itertools import pairwise
from pgGraphs import Graph


def correct_edges(gfa_graph: Graph) -> None:
    """Adds edges to the graph according to described paths

    Args:
        gfa_graph (Graph): _description_
    """
    # Removing old edges
    old_count: int = len(gfa_graph.lines)
    gfa_graph.lines = {}

    # Adding back edges in graph, based on paths
    for path_data in gfa_graph.paths.values():
        for (node_a, ori_a), (node_b, ori_b) in pairwise(path_data['path']):
            gfa_graph.add_edge(node_a, ori_a, node_b, ori_b)

    print(
        f"Corrected graph has {len(gfa_graph.lines)} edges whereas old graph had {old_count} edges."
    )
    # Old system for edges but still huge difference. Should be investigated!


def correct_graph(gfa_file: str, gfa_output: str) -> None:

    # Loading graph in memory
    graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True,
        low_memory=False,
        regexp=".*",
    )
    # Correcting edges
    correct_edges(graph)
    graph.save_graph(gfa_output)
