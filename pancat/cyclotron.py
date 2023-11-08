from gfagraphs import Graph
from tharospytools.bio_tools import revcomp


def get_graph_cycles(gfa_file: str, gfa_ver: str, mode: str = 'names') -> dict:
    """Extracts the graph cycles 

    Args:
        gfa_file (str): path to the gfa file
        gfa_ver (str): version of the gfa file
        chain_length (bool, optional): tells if we want to return cycle length instead of node names. Defaults to False.

    Returns:
        dict: a mapping between path names and the loops
    """

    graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_ver,
        with_sequence=True
    )
    cycles: dict = detect_cycles(graph)
    if mode == 'lengths':
        return {name: [sum([len(graph.mapping[node]) for node, _ in chain]) for chain in chains] for name, chains in cycles.items()}
    elif mode == 'sequences':
        return {name: [[graph.mapping[node] if vect.value == '+' else revcomp(graph.mapping[node]) for node, vect in chain] for chain in chains] for name, chains in cycles.items()}
    return {name: [[(node, vect.value) for node, vect in chain] for chain in chains] for name, chains in cycles.items()}


def detect_cycles_old(gfa_graph: Graph) -> dict[str, list]:
    """Detects the cycles in a graph, returning nodes names

    Args:
        gfa_graph (Graph): a loaded graph from the gfagraph lib

    Returns:
        dict[str,list]: a list of lists per path, featiring nodes names of each loop in graph
    """
    all_paths: list = gfa_graph.get_path_list()
    chains: dict[str, list[list[str]]] = {
        path.datas['name']: list() for path in all_paths}
    for path in all_paths:
        mappings: dict[str, list[int]] = {}
        for i, (node, _) in enumerate(path.datas['path']):
            if node not in mappings:
                mappings[node] = [i]
            else:
                chains[path.datas['name']].append(
                    path.datas['path'][mappings[node][-1]:i])
                for nd, _ in path.datas['path'][mappings[node][-1]:i]:
                    if nd in mappings:
                        del mappings[nd]
                mappings[node] = [i]
    return chains


def detect_cycles(gfa_graph: Graph) -> dict[str, list]:
    all_paths: list = gfa_graph.get_path_list()
    chains: dict[str, list[list[str]]] = {
        path.datas['name']: list() for path in all_paths}
    # We iterate over paths, resetting all counters
    for path in all_paths:
        mappings: dict[str, list[int]] = {
            nd.datas['name']: 0 for nd in gfa_graph.segments}
        iterations: list = [0 for _ in range(len(path.datas['path']))]
        # We construct the iteration list
        for i, (node, _) in enumerate(path.datas['path']):
            itr = mappings[node] + 1
            iterations[i] = itr
            mappings[node] = itr
        # Now we dissassemble the list
        for x in range(max(iterations), 1, -1):
            indices = [n for n, lvl in enumerate(iterations) if lvl == x]
        # We correct values in order to say we performed analysis on those loops


def cyclotron(gfa_graph: Graph) -> dict[str, list]:
    all_paths: list = gfa_graph.get_path_list()
    chains: dict[str, list[list[str]]] = {
        path.datas['name']: list() for path in all_paths}
    # We iterate over paths, resetting all counters
    for path in all_paths:
        x: int = 0
        dt: list = path.datas['path']
        mappings: dict[str, list[int]] = {nd: list() for nd, _ in dt}
        while x < len(dt):
            node, ori = dt[x]
            if node in mappings:
                pass
            mappings[node] = list(x)
