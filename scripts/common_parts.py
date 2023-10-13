from gfagraphs import Graph


def get_common_nodes(node_sequences: dict[str, set], number_of_haplotypes: int) -> dict[str, list]:
    """Return the nodes shared between *number_of_haplotypes* haplotypes

    Args:
        node_sequences (dict[str,set]): a mapping between haplotype name and nodes in this haplotype
        number_of_haplotypes (int): the cutoff

    Returns:
        dict[str,list]: mapping between node names and the haplotypes names
    """
    node_bank: dict[str, list] = {}
    for haplotype, sequence in node_sequences:
        for node in sequence:
            if node not in node_bank.keys():
                node_bank[node] = [haplotype]
            else:
                node_bank[node].append(haplotype)
    return {node: haplotypes for node, haplotypes in node_bank.items() if len(haplotypes) == number_of_haplotypes}


if __name__ == "__main__":
    gfa_graph: str = "path/to/gfa"
    gfa_ver: str = "GFA1.1"
    number_of_haplotypes: int = 2

    # Load the graph in memory
    graph: Graph = Graph(gfa_graph, gfa_ver, with_sequence=True)

    # Get nodes in each path
    seqs: dict[str, set] = {
        path.datas["name"]: set(
            [
                node for node, _ in path.datas["path"]
            ]
        ) for path in graph.get_path_list()
    }

    # Get shared nodes between haplotypes
    shared_nodes = get_common_nodes(seqs, number_of_haplotypes)
