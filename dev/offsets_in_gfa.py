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


def calculate_sequence_offsets(graph: Graph) -> tuple[dict[str, dict[str, list[tuple]]], dict[str, list[tuple[int, int]]]]:
    """Given a set of paths, calculates the offsets within each path

    Args:
        graph (Graph): a pangenome graph loaded in pgGraph format.

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


def add_offsets_to_gfa_old(gfa_file: str, output_file: str) -> None:
    """Given a GFA file, calculates the offsets within paths and stores it as complementary tags

    Args:
        gfa_file (str): input gfa file
        output_file (str): output gfa file
        gfa_version (str): the user-assumed gfa subformat
    """
    gfa_graph: Graph = Graph(gfa_file, with_sequence=True)

    nodes_offsets: list = calculate_sequence_offsets(
        graph=gfa_graph
    )[0]

    with open(output_file, 'w', encoding='utf-8') as gfa_writer:
        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            for line in gfa_reader:
                if line.startswith('S'):
                    node_name: str = line.split()[1]
                    gfa_writer.write(
                        f"{line.strip()}\tPO:J:{dumps(nodes_offsets[node_name],indent=0, separators=(',', ':'))}\n")
                else:
                    gfa_writer.write(line)
