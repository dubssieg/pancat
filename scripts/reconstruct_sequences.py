'Extract sequences given a graph'
from argparse import ArgumentParser, SUPPRESS
from typing import Generator
from gfagraphs import Record, GfaStyle, Walk, Path, Line, Segment


def grab_paths(gfa_file: str, gfa_version: str, reference: str) -> list[Record]:
    """From a given gfa-like file, grabs the path of each haplotype

    Args:
        gfa_file (str): a path to a gfa-like file
        gfa_version (str): a GFA-format identifier

    Raises:
        NotImplementedError: in case of rGFA which is currently not supportef

    Returns:
        list[Record]: list of paths, one for each haplotype
    """
    paths: list[Record] = []
    version: GfaStyle = GfaStyle(gfa_version)
    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            gfa_line: Record = Record(line, gfa_version)
            if isinstance(gfa_line, Walk):
                if not gfa_line.datas['id'] == '_MINIGRAPH_':
                    if gfa_line.datas['id'] == reference:
                        paths = [gfa_line] + paths
                    else:
                        paths.append(gfa_line)
            elif isinstance(gfa_line, Path):
                if gfa_line.datas['name'] == reference:
                    paths = [gfa_line] + paths
                else:
                    paths.append(gfa_line)
            elif isinstance(gfa_line, Line) and version == GfaStyle.RGFA:
                raise NotImplementedError('rGFA currently not supported')
            else:
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
        paths_to_follow (list[Record]): a list of paths
        node_start (str): optional, defaults to none ; node we cut graph from
        node_end (str): optional, defaults to none ; node we cut graph to

    Raises:
        ValueError: If taget nodes aren't in the backbone of the graph

    Returns:
        list[Record]: the subpaths within the nodes
    """
    if node_start is None and node_end is None:
        return paths_to_follow
    path_seq: list = [seq for seq,
                      _ in paths_to_follow[0].datas['path']]
    if node_start is None:
        node_start = path_seq[0]
    if node_end is None:
        node_end = path_seq[-1]
    try:
        paths_to_follow[0].datas['path'] = paths_to_follow[0].datas['path'][path_seq.index(
            node_start):path_seq.index(node_end)]
        path_seq: list = [seq for seq,
                          _ in paths_to_follow[0].datas['path']]
    except ValueError as exc:
        raise ValueError(
            'Query nodes are not in reference sequence. Please specify nodes within it.') from exc

    for path in paths_to_follow[1:]:
        path_query: list = [seq for seq,
                            _ in path.datas['path']]
        start_index: int = get_node_position(path_seq, path_query)
        stop_index: int = get_node_position(path_seq, path_query[::-1])
        path.datas['path'] = path.datas['path'][start_index:stop_index]
    return paths_to_follow


def get_node_position(path_seq: list, path_query: list) -> int:
    """Gets the first overlap of two lists

    Args:
        path_seq (list): a list of nodes
        path_query (list): another one

    Returns:
        int: position on query of the first overlap with ref
    """
    for node in path_seq:
        try:
            return path_query.index(node)
        except ValueError:
            pass
    return -1


def reconstruct(gfa_file: str, gfa_version: str, paths_to_follow: list[Record]) -> Generator:
    """Given a list of paths and a file, recreates sequence of the current path and yields it

    Args:
        gfa_file (str): a file to parse
        gfa_version (str): GFA subversion of the parse
        paths_to_follow (list[Record]): paths for the haplotypes inside the graph

    Yields:
        Generator: a list of ordered subsequences that are the DNA sequences the path goes through
    """
    for path in paths_to_follow:
        path_seq: list = [seq for seq, _ in path.datas['path']]
        reconstruct_seq: list = ['' for _ in range(len(path_seq))]
        with open(gfa_file, "r", encoding="utf-8") as reader:
            for line in reader:
                gfa_line: Record = Record(line, gfa_version)
                if isinstance(gfa_line, Segment):
                    try:
                        idx: int = path_seq.index(
                            gfa_line.datas['name'])
                        reconstruct_seq[idx] = line.split()[2]
                    except ValueError:
                        pass
            yield reconstruct_seq


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="Path to a gfa-like file")
    parser.add_argument(
        "out", type=str, help="Output path (without extension)")
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
        "-r",
        "--reference",
        help="Tells the reference sequence we seek start and stop into",
        required=True,
        type=str
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
    parser.add_argument(
        "-s", "--split", help="Tells to split in different files", action='store_true')
    args = parser.parse_args()

    followed_paths: list = node_range(grab_paths(
        args.file, args.gfa_version, args.reference), args.start, args.stop)

    if args.split:
        for i, sequence in enumerate(reconstruct(args.file, args.gfa_version, followed_paths)):
            with open(f"{args.out}_{i}.fasta", "w", encoding="utf-8") as writer:
                writer.write(
                    f"> {followed_paths[i].datas['name']}\n{''.join(sequence)}\n")
    else:
        with open(f"{args.out}.fasta", "w", encoding="utf-8") as writer:
            for i, sequence in enumerate(reconstruct(args.file, args.gfa_version, followed_paths)):
                writer.write(
                    f"> {followed_paths[i].datas['name']}\n{''.join(sequence)}\n")
