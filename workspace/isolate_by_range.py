"Requires PO tag to work."
from argparse import ArgumentParser, SUPPRESS
from collections import deque
from gfagraphs import Segment, Graph, Walk, Path, GfaStyle
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
    sequence_offsets, walks_offsets = calculate_sequence_offsets(
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
    mapping: dict = {}
    # This way does not work. we need to check positions relative to reference. See formulas in paper.

    node_set: set = set()
    # Filtering walks
    for walk_name, walk_sequence in walks_offsets.items():
        path_to_edit: Walk | Path = gfa_graph.get_path(walk_name)
        walk_start_index: int = 0
        walk_stop_index: int = 0
        for i, (x, y) in enumerate(walk_sequence):
            if x <= start and y >= start:
                walk_start_index = i
                path_to_edit.datas['start_offset'] = x
            elif x <= stop and y >= stop:
                walk_stop_index = i
                path_to_edit.datas['stop_offset'] = y
        path_to_edit.datas['path'] = path_to_edit.datas['path'][walk_start_index, walk_stop_index]
        # Retrieve nodes of walk, and adding them to the set
        node_set.add([node for node, _ in path_to_edit.datas['path']])

    # Filtering nodes
    gfa_graph.segments = [
        seg for seg in gfa_graph.segments if seg.datas['name'] in node_set]
    # Adding coordinates
    for seg in gfa_graph.segments:
        seg.datas['PO'] = sequence_offsets[seg.datas['name']]

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
