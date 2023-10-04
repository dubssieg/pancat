'Extract sequences given a graph'
from typing import Generator
from tharospytools import revcomp
from gfagraphs import Graph


def reconstruct_paths(gfa_file: str, gfa_version: str) -> Generator:
    """Given a GFA file with paths, follows those paths to reconstruct the linear sequences

    Args:
        gfa_file (str): gfa graph file, with paths
        gfa_version (str): the version of the file

    Raises:
        ValueError: if no path in graph

    Yields:
        Generator: mapping between path name and sequence
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_version,
        with_sequence=True
    )
    if len(graph_paths := gfa_graph.get_path_list()) > 0:
        yield from {path.datas["name"]: [
            gfa_graph.get_segment(x).datas['seq'] if vect == '+' else revcomp(
                gfa_graph.get_segment(x).datas['seq']
            ) for x, vect in path
        ] for path in graph_paths}
    else:
        raise ValueError(
            "Graph has no embed paths. Could'nt determine what to reconstruct!"
        )
