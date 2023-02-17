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
from gfagraphs import Record, GfaStyle


def calculate_sequence_offsets(node_data: list[tuple], walks: list[Record]) -> list[dict]:
    """Given a set of paths, calculates the offsets within each path

    Args:
        node_data (list[tuple]): series of tuple (name,length) for each node
        walks (list[Record]): the GFA walks

    Returns:
        list[dict]: a mapping {walk:offset_in_walk} for each walk, and for each node
    """
    offsets: list[int] = [0 for _ in range(len(walks))]
    orientations: list[str] = ['' for _ in range(len(walks))]
    sequence_offsets: list[dict] = list()
    walks_nodes: list[list] = [
        [node for node, _ in walk.datas["path"]] for walk in walks]
    walks_orientations: list[list] = [
        [ori for _, ori in walk.datas["path"]] for walk in walks]
    for (name, length) in node_data:
        inside_walks: list[int] = list()
        for i, walk in enumerate(walks_nodes):
            try:
                pos: int = walk.index(name)
                inside_walks.append(i)
                offsets[i] += length
                orientations[i] = str(walks_orientations[i][pos].value)
            except ValueError:
                pass
        sequence_offsets.append(
            {walks[x].datas["name"]: (offsets[x]-length, offsets[x], orientations[x])
             for x in inside_walks}
        )

    return sequence_offsets


def add_offsets_to_gfa(gfa_file: str, output_file: str, gfa_version: str) -> None:
    """Given a GFA file, calculates the offsets within paths and stores it as complementary tags

    Args:
        gfa_file (str): input gfa file
        output_file (str): output gfa file
        gfa_version (str): the user-assumed gfa subformat
    """
    if (GfaStyle(gfa_version)) == GfaStyle.RGFA:
        raise NotImplementedError(
            "Nodes can be extracted, but paths could not be determined.")
    embed_paths: list[Record] = list()
    nodes_information: list[tuple] = list()
    with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
        for line in gfa_reader:
            if (x := line.split())[0] == 'S':
                nodes_information.append((x[1], len(x[2])))
            elif x[0] == 'P' or x[0] == 'W':
                embed_paths += [Record(line, gfa_version)]
    nodes_offsets: list = calculate_sequence_offsets(
        nodes_information,
        embed_paths
    )
    collected_nodes: int = 0
    with open(output_file, 'w', encoding='utf-8') as gfa_writer:
        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            for line in gfa_reader:
                if line.startswith('S'):
                    gfa_writer.write(
                        f"{line.strip()}\tPO:J:{dumps(nodes_offsets[collected_nodes])}\n")
                    collected_nodes += 1
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
