"Computes edition between pangenome graphs"
from json import dump
from pgGraphs import Graph


def perform_edition(
        gfa_A: str,
        gfa_B: str,
        output_path: str,
        selection: list[str] | bool = True,
) -> tuple:
    """
    In this function, we do calculate the distance between G1 and G2, by trying to modify G2 into G1.
    Note the two graphs can be freely swapped, we just need to invert scores for reciprocal events and operations

    Args:
        gfa_A (str): a path to a gfa file
        gfa_B (str): a path to another gfa file
        output_path (str): path where to store editions
        selection (list[str] | bool, optional):
            - if true, compute edition on intersection of paths at path-level.
            - if false, compute edition at graph-level.
            - if a list[str], compute edition on specified paths
            Defaults to True.

    Returns:
        tuple: results from the edition
    """
    graph_A: Graph = Graph(gfa_file=gfa_A, with_sequence=True)
    graph_B: Graph = Graph(gfa_file=gfa_B, with_sequence=True)

    results: dict = dict()
    if selection:
        # We compute the intersection of paths in both graphs
        path_intersect: set[str] = set(
            graph_A.paths.keys()).intersection(set(graph_B.paths.keys()))
        if isinstance(selection, list):
            # We perform edition on selected paths, if all paths are in both graph
            if not all([x in path_intersect for x in selection]):
                raise ValueError()
            results = path_level_edition(graph_A, graph_B, set(selection))
        else:
            # We perform edition on shared paths, hoping the best for non-common paths \o/
            # (Best practice is to validate before if all paths are shared)
            results = path_level_edition(graph_A, graph_B, path_intersect)
    else:
        results = graph_level_edition()
        raise NotImplementedError(
            "Needs to be implemented, placeholder function"
        )
    dump(results, open(output_path, 'w', encoding='utf-8'))


def path_level_edition(graph_A: Graph, graph_B: Graph, selected_paths: set[str]) -> dict:
    """Compute edition, path by path, between the two graphs.
    The graph_A will be used as reference

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph
        selected_paths (set[str]): the paths where the edition needs to be computed

    Returns:
        dict: results of edition
    """
    edition_results: dict = dict()

    # Iterating on each pair of paths
    for path_name in selected_paths:

        i: int = 0  # counter of segmentations on graph_A
        j: int = 0  # counter of segmentations on graph_B

        merges: set[int] = set()  # set for merges
        splits: set[int] = set()  # set for merges

        pos_A: int = 0  # Absolute position in BP on A
        pos_B: int = 0  # Absolute position in BP on B

        global_pos: int = 0  # Position across both genomes

        # Iterating until we did not go through both segmentations
        while i < len(graph_A.paths[path_name]['path']) and j < len(graph_B.paths[path_name]['path']):
            # Currently evaluated nodes
            current_node_A: str = graph_A.paths[path_name]['path'][i][0]
            current_node_B: str = graph_B.paths[path_name]['path'][j][0]

            # We compute the next closest breakpoint
            global_pos = min(
                pos_A+graph_A.segments[current_node_A]['length'],
                pos_B+graph_B.segments[current_node_B]['length']
            )

            # We added the interval to current positions
            match (global_pos-pos_A == graph_A.segments[current_node_A]['length'], global_pos-pos_B == graph_B.segments[current_node_B]['length']):
                case (True, True):
                    # Iterating on both, no edition needed
                    pos_A += graph_A.segments[current_node_A]['length']
                    pos_B += graph_B.segments[current_node_B]['length']
                    i += 1
                    j += 1
                case (True, False):
                    # Iterating on top, split required
                    splits.add(global_pos)
                    pos_A += graph_A.segments[current_node_A]['length']
                    i += 1
                case (False, True):
                    # Iterating on bottom, merge required
                    merges.add(global_pos)
                    pos_B += graph_B.segments[current_node_B]['length']
                    j += 1
                case (False, False):
                    raise ValueError()

        edition_results[path_name] = {
            'merges': sorted(list(merges)),
            'splits': sorted(list(splits))
        }

    return edition_results


def graph_level_edition(graph_A: Graph, graph_B: Graph) -> dict:
    """Compute edition, at graph level, between the two graphs.
    The graph_A will be used as reference

    Args:
        graph_A (Graph): pangenome graph
        graph_B (Graph): pangenome graph

    Returns:
        dict: results of edition
    """
    edition_results: dict = dict()

    return edition_results
