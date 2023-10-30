from gfagraphs import Graph


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
        return {name: [sum([len(graph.mapping[node]) for node in chain]) for chain in chains] for name, chains in cycles.items()}
    elif mode == 'sequences':
        return {name: [[graph.mapping[node] for node in chain] for chain in chains] for name, chains in cycles.items()}
    return cycles


def detect_cycles(gfa_graph: Graph) -> dict[str, list]:
    """Detects the cycles in a graph, returning nodes names

    Args:
        gfa_graph (Graph): a loaded graph from the gfagraph lib

    Returns:
        dict[str,list]: a list of lists per path, featiring nodes names of each loop in graph
    """
    all_paths: list = gfa_graph.get_path_list()
    aliases: dict = {path.datas['name']: [
        node_name for node_name, _ in path.datas['path']] for path in all_paths}
    chains: dict[str, list[list[str]]] = {
        path.datas['name']: list() for path in all_paths}
    for path in all_paths:
        mappings: dict[str, list[int]] = {}
        for i, (node, _) in enumerate(path.datas['path']):
            if node not in mappings:
                mappings[node] = [i]
            else:
                chains[path.datas['name']].append(
                    aliases[path.datas['name']][mappings[node][-1]:i])
                for node in aliases[path.datas['name']][mappings[node][-1]:i]:
                    del mappings[node]
    return chains
