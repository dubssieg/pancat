"Extract connex components where events are occuring."
from pathlib import Path
from gfagraphs import Graph

# We will figure event style as True/False => True will be equivalence


def get_edit_bubbles(graph_1_path: str, graph_1_version: str, graph_2_path: str, graph_2_version: str, reference: str = "heifer_diploid"):
    """From two GFA graphs, function must return a list of pairs of graphs
    that are the edit components of each graph
    """
    # Loading the two graphs
    graph_1_gfagraphs: Graph = Graph(
        graph_1_path, graph_1_version, with_sequence=True)
    graph_2_gfagraphs: Graph = Graph(
        graph_2_path, graph_2_version, with_sequence=True)

    # Extracting each node offset
    node_map_1_offsets: dict = {
        node.datas["name"]: node.datas['PO'] for node in graph_1_gfagraphs.segments}
    node_map_2_offsets: dict = {
        node.datas["name"]: node.datas['PO'] for node in graph_2_gfagraphs.segments}

    # Extracting paths we will follow
    paths_1_genomes: dict = {walk.datas['name']: walk.datas["path"]
                             for walk in graph_1_gfagraphs.get_path_list()}
    paths_2_genomes: dict = {walk.datas['name']: walk.datas["path"]
                             for walk in graph_2_gfagraphs.get_path_list()}

    path_of_first = {key: [(n, o) for (n, o) in val if n in node_map_1_offsets.keys()] for key,
                     val in paths_1_genomes.items()}
    path_of_second = {key: [(n, o) for (n, o) in val if n in node_map_2_offsets.keys()] for key,
                      val in paths_2_genomes.items()}

    name_graph_a, name_graph_b = Path(
        graph_1_path).stem, Path(graph_2_path).stem

    equivalence_list: list = list()

    for name, path_a in path_of_first.items():
        path_b = path_of_second[name]
        pointer_1, pointer_2 = 0, 0
        while pointer_1 < len(path_a) and pointer_2 < len(path_b):

            try:
                start_a, end_a, _ = node_map_1_offsets[path_a[pointer_1][0]][name]
                start_b, end_b, _ = node_map_2_offsets[path_b[pointer_2][0]][name]
            except KeyError:
                # Path follows reference (might be ambiguous)
                start_a, end_a, _ = node_map_1_offsets[path_a[pointer_1]
                                                       [0]][reference]
                start_b, end_b, _ = node_map_2_offsets[path_b[pointer_2]
                                                       [0]][reference]

            name_a = f"{name_graph_a}_{path_a[pointer_1][0]}"
            name_b = f"{name_graph_b}_{path_b[pointer_2][0]}"

            if start_a == start_b and end_a == end_b:
                equivalence_list.append((
                    name_a,
                    name_b,
                    name,
                    True
                ))
                pointer_1 += 1
                pointer_2 += 1

            elif start_a >= start_b and end_a <= end_b:
                equivalence_list.append((
                    name_a,
                    name_b,
                    name,
                    False
                ))
                pointer_1 += 1
            elif start_b >= start_a and end_b <= end_a:
                equivalence_list.append((
                    name_a,
                    name_b,
                    name,
                    False
                ))
                pointer_2 += 1
            elif start_a < start_b and end_a < end_b and end_a != start_b and start_b < end_a:

                equivalence_list.append((
                    name_a,
                    name_b,
                    name,
                    False
                ))
                pointer_1 += 1
            elif start_b < start_a and end_b < end_a and end_b != start_a and start_a < end_b:
                equivalence_list.append((
                    name_a,
                    name_b,
                    name,
                    False
                ))
                pointer_2 += 1

            else:
                if end_a <= start_b:
                    pointer_1 += 1
                else:
                    pointer_2 += 1

    return equivalence_list


def extract_bubble(
        equivalence_list: list,
    reference: str = "heifer_diploid"
) -> list:
    """Given a equivalence series, find souces and sinks of edition bubbles

    Args:
        equivalence_list (list): equivalence series

    Returns:
        list: list of tuples ((source_1,source_2),(sink_1,sink_2))
    """
    sources_sinks: list = list()
    idx: int = 0
    source_nodes: tuple = (equivalence_list[0][0], equivalence_list[0][1])
    in_component: bool = not equivalence_list[0][3]
    for event in equivalence_list:
        # status = True if all genomes go through node_a and node_b
        node_1, node_2, status, event_style = event
        if event_style and status == reference:
            if in_component:
                # we found a sink
                in_component = False
                idx += 1
                sink_nodes: tuple = (node_1, node_2)
                sources_sinks.append((source_nodes, sink_nodes))
            else:
                # we found a potential source
                source_nodes: tuple = (node_1, node_2)
        elif not event_style:
            in_component = True

    return sources_sinks


def compute_extraction(sources_sinks: list, ratio_to_add: int = 2) -> list:
    """Given each set of sources and sinks
    compute middle node and number of nodes to extract 
    returns a list of tuple ((node_graph_1_center,nb_nodes),(node_graph_2_center,nb_nodes))
    we add a ratio as BFS step is not the best :( (this ratio is cursed btw)
    """
    extraction_data: list = list()
    for ((source_1, source_2), (sink_1, sink_2)) in sources_sinks:
        nb_nodes_1 = int(sink_1.split('_')[-1]) - int(source_1.split('_')[-1])
        nb_nodes_2 = int(sink_2.split('_')[-1]) - int(source_2.split('_')[-1])
        center_node_1 = int(source_1.split('_')[-1]) + nb_nodes_1//2
        center_node_2 = int(source_2.split('_')[-1]) + nb_nodes_2//2
        extraction_data.append(
            ((center_node_1, nb_nodes_1+nb_nodes_1//ratio_to_add), (center_node_2, nb_nodes_2+nb_nodes_2//ratio_to_add)))
    return extraction_data


if __name__ == "__main__":

    graph_1 = "/udd/sidubois/Documents/Code/data/FIRST_1M/cactus_default.gfa"
    graph_2 = "/udd/sidubois/Documents/Code/data/FIRST_1M/abberant.gfa"
    ver = "GFA1.1"
    equivalences = get_edit_bubbles(graph_1, ver, graph_2, ver)
    source_sink = extract_bubble(equivalences)
    extraction_points = compute_extraction(source_sink)

    print(extraction_points)
