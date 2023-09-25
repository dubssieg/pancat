from typing import Generator
from gfagraphs import Graph


def grouper(iterable, n=2, m=1):
    """Collect data into overlapping fixed-length chunks or blocks"""
    return [iterable[i:i+n] for i in range(0, len(iterable)-1, n-m)]


def common_members(elements: list[set]):
    first_path, other_paths = elements[0], elements[1:]
    return sorted(list(first_path.intersection(*other_paths)))


def bubble_caller(gfa_graph: Graph) -> list[dict]:
    """Calls out the bubbles in the graph.
    A bubble can be defined as having a starting and an ending node
    with a in and out node with degree equal to the number of paths
    for superbubble level we don't have to watch the order, as 

    Args:
        gfa_file (str): path to a gfa-like file

    Returns:
        list[dict]: a list of mappings between paths names and the subchain in the bubble
                    one element per bubble
    """
    gfa_paths: list = gfa_graph.get_path_list()

    all_sets = {
        path.datas['name']:
            [
                node_name for node_name, _ in path.datas['path']
        ]
        for path in gfa_paths
    }

    bubbles_endpoints: list = sorted(common_members(
        list(
            set(x) for x in all_sets.values()
        )
    ), key=int)
    bubbles: list[dict] = [{}
                           for _ in range(len(bubbles_endpoints)-1)]
    for path in gfa_paths:
        # Computing endpoint positions in list for each path
        endpoints_indexes: list = grouper(
            [
                all_sets[
                    path.datas['name']
                ].index(
                    endpoint
                ) for endpoint in bubbles_endpoints
            ],
            2
        )
        print(endpoints_indexes)
        # Getting bubble chains
        for i, (start, end) in enumerate(endpoints_indexes):
            bubbles[i][path.datas['name']
                       ] = all_sets[path.datas['name']][start:end+1]
    return bubbles


def call_variants(gfa_file: str, gfa_type: str, reference_name: str) -> Generator:
    """Given a GFA file and a path name, calls rank 1 variants against it

    Args:
        gfa_file (str): path to a gfa file
        gfa_type (str): subformat
        reference_name (str): a path name in the gfa file
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_type,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
    print(bubbles)
    for bubble in bubbles:
        yield {path_name: ''.join([gfa_graph.get_segment(node=node).datas['seq'] for node in path_chain]) for path_name, path_chain in bubble.items()}


"""
def flatten_graph(gfa_file: str, gfa_type: str) -> None:
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_type,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
"""
