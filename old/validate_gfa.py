"This scripts aims to validate one or multiple GFA files, by checking correspondance between path names, and their type"
from argparse import ArgumentParser, SUPPRESS
from collections import Counter
from itertools import chain
from gfagraphs import Graph


def load_graphs(files: list[str]) -> dict[str, Graph]:
    """Loads all graphs topologies

    Args:
        files (list[str]): a set of files to evaluate

    Returns:
        dict[str, Graph]: mappng filename -> Graph object
    """
    return {file: Graph(file, 'unknown') for file in files}


def evaluate(graphs: dict[str, Graph]) -> list[dict]:
    """Given a dict of graphs, extracts comparaison information between graphs

    Args:
        graphs (dict[str, Graph]): graph topologies

    Returns:
        list(dict): informations about each file
    """
    graph_infos: list[dict] = [{} for _ in range(len(graphs))]
    for i, (file_path, graph) in enumerate(graphs.items()):
        graph_infos[i]['filename'] = file_path.split('/')[-1]
        graph_infos[i]['version'] = (ver := graph.assert_format().value)
        if ver not in ['RGFA', 'unknown', 'GFA2']:
            graph_infos[i]['embed_paths'] = set([path.datas['name']
                                                 for path in graph.get_path_list()])
    return graph_infos


def assert_compatibility(infos: list[dict]) -> bool | None:

    if any([dct['version'] in ['RGFA', 'unknown', 'GFA2'] for dct in infos]):
        return None
    freqs: Counter = Counter(chain.from_iterable(
        [dct['embed_paths']for dct in infos]))
    return len({idx for idx in freqs if freqs[idx] == 1}) == 0


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='This scripts aims to validate one or multiple GFA files, by checking correspondance between path names, and their type.')
    args = parser.parse_args()

    if len(graph_datas := evaluate(load_graphs(args.file))) > 1:
        if (status := assert_compatibility(graph_datas)) is None:
            print(
                "Beware, some files are in formats that does not work well with pangraphs")
        elif status:
            print("Files are compatible")
        else:
            print("Files have conflicting path names. Please fix.")
