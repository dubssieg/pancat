"Requires PO tag to work."
from json import dumps
from gfagraphs import Graph
from workspace.offset_in_gfa import calculate_sequence_offsets
from workspace.find_bubbles import common_members
from tharospytools import path_allocator


def range_isolate(gfa_file: str, gfa_ver: str, output: str, reference_name: str, start: int, stop: int) -> None:
    # Load graph in memory
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        gfa_type=gfa_ver,
        with_sequence=True
    )
    # Computing offsets
    embed_paths: list = gfa_graph.get_path_list()
    nodes_information: dict = {
        node.datas["name"]: len(node.datas["seq"]) for node in gfa_graph.segments
    }
    sequence_offsets, _ = calculate_sequence_offsets(
        nodes_information,
        embed_paths
    )

    # Getting sources and sinks
    all_sets = {
        path.datas['name']:
            [
                node_name for node_name, _ in path.datas['path']
        ]
        for path in embed_paths
    }
    bubbles_endpoints: list = common_members(
        list(
            set(x) for x in all_sets.values()
        )
    )
    # We filter true sources and sinks
    mapping: dict = {x: sequence_offsets[x][reference_name][0] for x in bubbles_endpoints if all([len(
        sequence_offsets[x][name]) == 1 for name in all_sets.keys()])}

    source, penalty_source = False, float('inf')
    sink, penalty_sink = False, float('inf')
    for node_name, (offset_x, offset_y, _) in mapping.items():
        if min(penalty_source, abs(start-offset_x)) == abs(start-offset_x):
            penalty_source = abs(start-offset_x)
            source = node_name
        elif min(penalty_sink, abs(offset_y-stop)) == abs(offset_y-stop):
            penalty_sink = abs(offset_y-stop)
            sink = node_name

    # We check we have valid source and sink
    if not (source and sink):
        raise ValueError(
            "Please select other boundaries, as it was not possible to find appropriate source and sink.")

    node_set: set = set()
    # Filtering walks
    for walk in embed_paths:
        walk_name: str = walk.datas['name']
        walk_start_index: int = 0
        walk_stop_index: int = 0

        for i, (node, _) in enumerate(walk.datas['path']):
            if node == source:
                walk_start_index = i
                walk.datas['start_offset'] = sequence_offsets[node][walk_name][0][0]
            elif node == sink:
                walk_stop_index = i
                walk.datas['stop_offset'] = sequence_offsets[node][walk_name][0][1]
        walk.datas['path'] = walk.datas['path'][walk_start_index:walk_stop_index]
        # Retrieve nodes of walk, and adding them to the set
        node_set |= set([node for node, _ in walk.datas['path']])

    # Filtering nodes
    gfa_graph.segments = [
        seg for seg in gfa_graph.segments if seg.datas['name'] in node_set]
    # Adding coordinates
    for seg in gfa_graph.segments:
        seg.datas['PO'] = dumps(
            sequence_offsets[seg.datas['name']], indent=0, separators=(',', ':'))

    # Filtering edges
    gfa_graph.lines = [line for line in gfa_graph.lines if line.datas['start']
                       in node_set and line.datas['end'] in node_set]

    # Allocate output path and write the file
    gfa_graph.save_graph(
        output_path=path_allocator(
            path_to_validate=output,
            particle='.gfa',
            default_name=f'extracted_graph_{start}_{stop}'
        ),
        output_format=gfa_graph.version
    )
