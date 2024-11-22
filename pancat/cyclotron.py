from pgGraphs import Graph


def revcomp(string: str, compl: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}) -> str:
    """Tries to compute the reverse complement of a sequence

    Args:
        string (str): original character set
        compl (dict, optional): dict of correspondances. Defaults to {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.

    Raises:
        IndexError: Happens if revcomp encounters a char that is not in the dict

    Returns:
        str: the reverse-complemented string
    """
    try:
        return ''.join([compl[s] for s in string][::-1])
    except IndexError as exc:
        raise IndexError(
            "Complementarity does not include all chars in sequence.") from exc


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
        with_sequence=True
    )
    cycles: dict = detect_cycles(graph)
    if mode == 'lengths':
        return {name: [sum([len(graph.segments[node]['length']) for node, _ in chain]) for chain in chains] for name, chains in cycles.items()}
    elif mode == 'sequences':
        return {name: [[graph.segments[node]['seq'] if vect.value == '+' else revcomp(graph.segments[node]['seq']) for node, vect in chain] for chain in chains] for name, chains in cycles.items()}
    return {name: [[(node, vect.value) for node, vect in chain] for chain in chains] for name, chains in cycles.items()}


def detect_cycles_old(gfa_graph: Graph) -> dict[str, list]:
    """Detects the cycles in a graph, returning nodes names

    Args:
        gfa_graph (Graph): a loaded graph from the gfagraph lib

    Returns:
        dict[str,list]: a list of lists per path, featiring nodes names of each loop in graph
    """
    chains: dict[str, list[list[str]]] = {
        path: list() for path in gfa_graph.paths.keys()}
    for path_name, path_datas in gfa_graph.paths.items():
        mappings: dict[str, list[int]] = {}
        for i, (node, _) in enumerate(path_datas['path']):
            if node not in mappings:
                mappings[node] = [i]
            else:
                chains[path_name].append(
                    path_datas['path'][mappings[node][-1]:i])
                for nd, _ in path_datas['path'][mappings[node][-1]:i]:
                    if nd in mappings:
                        del mappings[nd]
                mappings[node] = [i]
    return chains


def detect_cycles(gfa_graph: Graph) -> dict[str, list]:
    chains: dict[str, list[list[str]]] = {
        path: list() for path in gfa_graph.paths.keys()}
    # We iterate over paths, resetting all counters
    for path_name, path_datas in gfa_graph.paths.items():
        mappings: dict[str, list[int]] = {
            nd: 0 for nd in gfa_graph.segments.keys()}
        iterations: list = [0 for _ in range(len(path_datas['path']))]
        # We construct the iteration list
        for i, (node, _) in enumerate(path_datas['path']):
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
