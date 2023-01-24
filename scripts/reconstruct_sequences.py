'Extract sequences given a graph'
from argparse import ArgumentParser, SUPPRESS
from typing import Generator
from gfatypes import LineType, Record, GfaStyle


def grab_paths(gfa_file: str, gfa_version) -> list[Record]:
    """_summary_

    Args:
        gfa_file (str): _description_
        gfa_version (_type_): _description_

    Raises:
        NotImplementedError: _description_

    Returns:
        list[Record]: _description_
    """
    paths: list[Record] = []
    version: GfaStyle = GfaStyle(gfa_version)
    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            gfa_line: Record = Record(line, gfa_version)
            match (gfa_line.linetype, version):
                case LineType.WALK, _:
                    if not gfa_line.line.idf == '_MINIGRAPH_':  # type:ignore
                        paths.append(gfa_line)
                case LineType.PATH, _:
                    paths.append(gfa_line)
                case LineType.LINE, GfaStyle.RGFA:
                    raise NotImplementedError('rGFA currently not supported')
                case _:
                    pass
    return paths


def node_range(
    paths_to_follow: list[Record],
        node_start: str | None = None,
        node_end: str | None = None
) -> list[Record]:
    """Assuming first record is the reference
    edits the list of records to keep only nodes in subpath

    Args:
        paths_to_follow (list[Record]): _description_
        node_start (str): _description_
        node_end (str): _description_

    Raises:
        ValueError: _description_

    Returns:
        list[Record]: _description_
    """
    if node_start is None and node_end is None:
        return paths_to_follow
    path_seq: list = [seq for seq,
                      _ in paths_to_follow[0].line.path]  # type:ignore
    if node_start is None:
        node_start = path_seq[0]
    if node_end is None:
        node_end = path_seq[-1]
    try:
        paths_to_follow[0].line.path = paths_to_follow[0].line.path[path_seq.index(  # type:ignore
            node_start):path_seq.index(node_end)]
        path_seq: list = [seq for seq,
                          _ in paths_to_follow[0].line.path]  # type:ignore
    except ValueError as exc:
        raise ValueError(
            'Query nodes are not in reference sequence. Please specify nodes within it.') from exc

    for path in paths_to_follow[1:]:
        path_query: list = [seq for seq,
                            _ in path.line.path]  # type:ignore
        start_index: int = get_node_position(path_seq, path_query)
        stop_index: int = get_node_position(path_seq, path_query[::-1])
        path.line.path = path.line.path[start_index:stop_index]  # type:ignore
    return paths_to_follow


def get_node_position(path_seq: list, path_query: list) -> int:
    """Gets the first overlap of two lists

    Args:
        path_seq (list): _description_
        path_query (list): _description_

    Returns:
        int: _description_
    """
    for node in path_seq:
        try:
            return path_query.index(node)
        except ValueError:
            pass
    return -1


def reconstruct(gfa_file: str, gfa_version: str, paths_to_follow: list[Record]) -> Generator:
    """_summary_

    Args:
        gfa_file (str): _description_
        gfa_version (str): _description_
        paths_to_follow (list[Record]): _description_

    Yields:
        Generator: _description_
    """
    for path in paths_to_follow:
        path_seq: list = [seq for seq, _ in path.line.path]  # type:ignore
        reconstruct_seq: list = ['' for _ in range(len(path_seq))]
        with open(gfa_file, "r", encoding="utf-8") as reader:
            for line in reader:
                gfa_line: Record = Record(line, gfa_version)
                if isinstance(gfa_line.line, Record.Segment):
                    try:
                        idx: int = path_seq.index(
                            gfa_line.line.name)  # type:ignore
                        reconstruct_seq[idx] = line.split()[2]
                    except ValueError:
                        pass
            yield reconstruct_seq


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument("out", type=str, help="Output path (with extension)")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Reconstruct linear sequences from a GFA graph')
    parser.add_argument(
        "-g",
        "--gfa_version",
        help="Tells the GFA input style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
    )
    parser.add_argument(
        '--start',
        type=str,
        help='To specifiy a starting node on reference to create a subgraph',
        default=None
    )
    parser.add_argument(
        '--stop',
        type=str,
        help='To specifiy a ending node on reference to create a subgraph',
        default=None
    )
    args = parser.parse_args()

    followed_paths: list = node_range(grab_paths(
        args.file, args.gfa_version), args.start, args.stop)

    with open(args.out, "w", encoding="utf-8") as writer:
        for i, sequence in enumerate(reconstruct(args.file, args.gfa_version, followed_paths)):
            writer.write(
                f"> {followed_paths[i].line.name}\n{''.join(sequence)}\n")  # type:ignore
