"Inspect graphs."
from typing import Generator
from argparse import ArgumentParser, SUPPRESS
from gfatypes import LineType, Record


def extract_lines(gfa_file: str, line_type: LineType, gfa_version: str) -> Generator:
    """Extract one type of lines form a GFA.

    Args:
        gfa_file (str): path to a valid GFA-like file
        line_type (LineType): type for lines we want to extract
        gfa_version (str): GFA file subversion

    Returns:
        list[Record]: all lines matching LineType

    Yields:
        Iterator[Record]: iterates over lines matching the desired LineType
    """
    with open(gfa_file, "r", encoding="utf-8") as reader:
        for line in reader:
            if line[0] == line_type.value:
                yield Record(line, gfa_version)


def grab_paths(gfa_file: str, gfa_version: str) -> list[Record]:
    """Grabs all paths

    Args:
        gfa_file (str): _description_
        gfa_version (str): _description_

    Raises:
        ExceptionGroup: _description_

    Returns:
        list[Record]: _description_
    """
    try:
        return [rec for rec in extract_lines(gfa_file, LineType.WALK, gfa_version)]
    except (ValueError, NotImplementedError) as ex1:
        try:
            return [rec for rec in extract_lines(gfa_file, LineType.PATH, gfa_version)]
        except (ValueError, NotImplementedError) as ex2:
            raise ExceptionGroup('Invalid file format', [ex1, ex2])


def nodepairs(paths: list[Record]) -> Generator:
    """If used on multiple contigs, must merge paths to be used

    Args:
        paths (list[Record]): _description_

    Yields:
        Generator: _description_
    """
    for i, element in enumerate(paths[0].line.path[:-1]):  # type:ignore
        try:
            indexes = [paths[i].line.path.index(  # type:ignore
                element) for i in range(1, len(paths))]
            if all([paths[0].line.path[i+1] == paths[p+1].line.path[idx+1] for p, idx in enumerate(indexes)]):
                yield element, paths[0].line.path[i+1]
        except (ValueError, IndexError):
            pass


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to one or more gfa-like file(s).")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Evaluates inconstistencies in topology.')
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    args = parser.parse_args()

    i: int = 0
    for _ in nodepairs(grab_paths(args.file, args.gfa_version)):
        i += 1
    print(i)
