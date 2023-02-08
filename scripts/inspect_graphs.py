"Inspect graphs."
from typing import Generator
from argparse import ArgumentParser, SUPPRESS
from gfatypes.gfatypes import LineType, Record


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


def count_snps(input_file: str, paths: dict) -> dict:
    """Counts the SNP of all paths in a gfa file

    Args:
        input_file (str): path to a GFA-like file

    Returns:
        Counter: length distribution of sequences
    """
    count_snp = {name: 0 for name in paths.keys()}
    with open(input_file, 'r', encoding='utf-8') as gfa_reader:
        for i, seq in enumerate(gfa_reader):
            datas = seq.split()
            if datas[0] == 'S' and len(datas[2]) == 1:
                for name in count_snp.keys():
                    try:
                        paths[name].remove(datas[1])
                        count_snp[name] += 1
                    except:
                        pass

    return count_snp


if __name__ == '__main__':

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path to a gfa-like file.")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Evaluates inconstistencies in topology.')
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    args = parser.parse_args()

    paths = grab_paths(args.file, args.gfa_version)
    dict_paths = {path.line.name: [node for node,
                                   _ in path.line.path] for path in paths}
    print(len(dict_paths))
    counts = count_snps(args.file, dict_paths)
    print(counts)

    with open('report.txt', 'w', encoding='utf-8') as report:
        i: int = 0
        for _ in nodepairs(paths):
            i += 1
        print(f"Total number of splitted sequences : {i}", file=report)
        print(f"SNP counts per walk :\n{counts}")
