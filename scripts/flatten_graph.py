from gfagraphs import Graph, Segment


def load_graph(gfa_file: str, gfa_version: str, with_sequences: bool = False) -> Graph:
    """Loads a graph from a file

    Args:
        gfa_file (str): path to file
        gfa_version (str): gfa subformat
        with_sequences (bool, optional): if sequences should be loaded. Defaults to False.

    Returns:
        Graph: a gfagraphs graph object
    """
    return Graph(gfa_file, gfa_version, with_sequence=with_sequences)


def merge_nodes(nodes_to_merge: list[Segment]) -> Segment:
    """Given a list of nodes, edits the lexicographically first node 

    Args:
        nodes_to_merge (list[Segment]): _description_

    Returns:
        Segment: _description_
    """
