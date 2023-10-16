"""
Here, we aim to extend the current GFA tag format by adding tags
that do respect the GFA naming convention.
A JSON string, PO (Path Offset) positions, relative to paths.
Hence, PO:J:{'w1':(334,335,'+'),'w2':(245,247,'-')} tells that the walk/path w1
contains the sequence starting at position 334 and ending at position 335,
and the walk/path w2 contains the sequence starting at the offset 245 (ending 247),
and that the sequences are reversed one to each other.
Note that any non-referenced walk in this field means that the node
is not inside the given walk.
"""
from argparse import ArgumentParser, SUPPRESS
from json import dumps
from gfagraphs import Graph


def calculate_sequence_offsets(node_data: dict, walks: list) -> tuple[dict[str, dict[str, list[tuple]]], dict[str, list[tuple[int, int]]]]:
    """Given a set of paths, calculates the offsets within each path

    Args:
        node_data (list[tuple]): series of tuple (name,length) for each node
        walks (list[Record]): the GFA walks

    Returns:
        tuple[dict[str, dict[str, list[tuple]]],dict[str,list[tuple[int,int]]]]: a mapping node_name:{walk:offset_in_walk} for each walk, and for each node + a mapping path name : coords of each node
    """
    sequence_offsets: dict[str, dict[str, list[tuple[int, int, str]]]] = {
        node_name: {} for node_name in node_data.keys()}

    walk_offsets: dict = dict()
    for walk in walks:
        walk_name: str = walk.datas["name"]
        offsets_of_items: list = list()
        start_offset: int = int(
            walk.datas['start_offset']) if 'start_offset' in walk.datas.keys() else 0
        for node, vect in walk.datas["path"]:
            offsets_of_items.append(
                (start_offset, start_offset+node_data[node]))
            if walk_name in sequence_offsets[node]:
                # We already encountered the node in this path
                sequence_offsets[node][walk_name].append(
                    (start_offset, start_offset+node_data[node], vect.value))
            else:
                # First time we encounter this node for this path
                sequence_offsets[node][walk_name] = [
                    (start_offset, start_offset+node_data[node], vect.value)]
            start_offset += node_data[node]
        walk_offsets[walk_name] = offsets_of_items

    return sequence_offsets, walk_offsets


def add_offsets_to_gfa(gfa_file: str, output_file: str, gfa_version: str) -> None:
    """Given a GFA file, calculates the offsets within paths and stores it as complementary tags

    Args:
        gfa_file (str): input gfa file
        output_file (str): output gfa file
        gfa_version (str): the user-assumed gfa subformat
    """
    gfa_graph: Graph = Graph(gfa_file, gfa_version, with_sequence=True)
    embed_paths: list = gfa_graph.get_path_list()
    nodes_information: dict = {node.datas["name"]: len(
        node.datas["seq"]) for node in gfa_graph.segments}

    nodes_offsets: list = calculate_sequence_offsets(
        nodes_information,
        embed_paths
    )[0]

    with open(output_file, 'w', encoding='utf-8') as gfa_writer:
        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            for line in gfa_reader:
                if line.startswith('S'):
                    node_name: str = line.split()[1]
                    gfa_writer.write(
                        f"{line.strip()}\tPO:J:{dumps(nodes_offsets[node_name])}\n")
                else:
                    gfa_writer.write(line)


if __name__ == '__main__':
    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument("out", type=str, help="Output path (with extension)")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Adds path offsets information to GFA')
    parser.add_argument(
        "-g",
        "--gfa_version",
        help="Tells the GFA input style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
    )
    args = parser.parse_args()

    add_offsets_to_gfa(args.file, args.out, args.gfa_version)
