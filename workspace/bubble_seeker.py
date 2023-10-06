from typing import Generator
from gfagraphs import Graph, Walk, Path
from tharospytools import path_allocator, flatten


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

    # loop until convergence
    convergence: int = 0
    while (len(bubbles_endpoints) != convergence):
        convergence: int = len(bubbles_endpoints)
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
            print(bubbles_endpoints)
            # Getting bubble chains
            for i, (start, end) in enumerate(endpoints_indexes):
                bubbles[i][path.datas['name']
                           ] = all_sets[path.datas['name']][start:end+1]

            # Decompressing all paths
            embed_nodes: set = set(flatten(
                [chain_content[1:-1] for bubble in bubbles for _, chain_content in bubble.items()]))
            # Search for endpoints that are in the set
            # if in the middle of the chain we notice a endpoint, THIS IS NO ENDPOINT and we need to clear it
            bubbles_endpoints = [
                endpoint for endpoint in bubbles_endpoints if endpoint not in embed_nodes]
            # If we need to suppress endpoints, we will have different length, so we will loop
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
    for bubble in bubbles:
        yield {
            path_name: ''.join(
                [gfa_graph.get_segment(node=node).datas['seq']
                 for node in path_chain]
            ) for path_name, path_chain in bubble.items()
        }


"""
def flatten_graph(gfa_file: str, gfa_type: str) -> None:
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_type,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
"""


def flattenable_bubbles(bubbles: list[dict]) -> Generator:
    """Only returns superbubbles given a list of bubbles

    Args:
        bubbles (list[dict]): a list of mixed bubbles

    Yields:
        list[dict]: a list of superbubbles
    """
    yield from [bubble for bubble in bubbles if any([len(chain) > 3 for chain in bubble.values()])]


def linearize_bubbles(gfa_file: str, gfa_type: str, output: str) -> Generator:
    """Given a GFA file, flattens the bubbles

    Args:
        gfa_file (str): path to a gfa file
        gfa_type (str): subformat
        output (str): output file path
    """
    output_path: str = path_allocator(
        output, particle=".gfa", default_name="graph")
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_type,
        with_sequence=True)
    bubbles: list[dict] = bubble_caller(gfa_graph=gfa_graph)
    # Init return graph
    output_graph: Graph(
        gfa_file=None,
        gfa_type=gfa_type,
        with_sequence=True
    )
    copied_nodes: set = set()
    # For each bubble, we compute new nodes
    for bubble in bubbles:
        for path_name, node_chain in bubble.items():
            pass


if __name__ == "__main__":
    all_bubbles: list[dict] = bubble_caller(Graph(gfa_file="/home/sidubois/Workspace/Notes/graph_cactus_sandra_extract.gfa",
                                                  gfa_type='GFA1.1', with_sequence=True))
    for bubble in all_bubbles:
        print(bubble)


def chain_readway(graph: Graph, chain: list, path_name: str) -> list[str]:
    """Extracts from paths the information about the way we must read.
    Beware: this does not function with current position and will stop at the first match
    meaning that if 

    Args:
        graph (Graph): a gfa graph from gfagraphs library
        chain (list): the subchain of the bubble we seeking for
        path_name (str): name of the current path

    Returns:
        list[str]: orientations segments must be read
    """
    wayline: Walk | Path = graph.get_path(path_name)
    chain_length: int = len(chain)
    for ind in (i for i, (e, _) in enumerate(wayline.datas['path']) if e == chain[0]):
        if [x for x, _ in wayline.datas['path'][ind:ind+chain_length]] == chain:
            return [z for _, z in wayline.datas['path'][ind:ind+chain_length]]

# [(1, 3), (8, 10)]

# WARNING bubble seeker needs to take into account the direction it reads nodes :(

# bubbles must be expanded > search for pattern in the bubble after to get a longer chain. if subchain > expand, else go next
# it means we verify if bubbles are strictly independant
# source, chain, sink = x[0],x[1:-1],x[-1]
# if x_next[1] in chain:
# merge_index = chain.index(x_next[1])


def expand_bubbles(bubbles: list[dict]) -> list[dict]:
    current_bubble: dict = bubbles[0]
