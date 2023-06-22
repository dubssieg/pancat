"This scripts aims to detect loops inside pangenome graphs, by checking on each path if a node is crossed multiple times or not"
from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from gfagraphs import Graph
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def assess_loop(gfa_graph: str, gfa_ver: str) -> None:
    """Plots nodes implicated in loops in the graph

    Args:
        gfa_graph (str): a path to a gfa file
        gfa_ver (str): a gfa sub-version
    """

    datas: dict = dict()
    graph_to_inspect: Graph = Graph(gfa_file=gfa_graph, gfa_type=gfa_ver)
    all_paths: list = graph_to_inspect.get_path_list()

    for path in all_paths:
        visits: dict = dict()
        for i, (node, _) in enumerate(path.datas['path']):
            visits[node] = visits.get(node, []) + [i]
        x_data: list = [len(value)
                        for value in visits.values() if len(value) > 1]
        y_data: list = [(value[-1]-value[0])
                        for value in visits.values() if len(value) > 1]
        datas = {'x_data': datas.get('x_data', []) + x_data, 'y_data': datas.get('y_data', []) +
                 y_data, 'path_name': datas.get('path_name', []) + [path.datas['name'] for _ in x_data]}
        print(
            f"{path.datas['name']} => {Counter([len(value) for value in visits.values()])}")

    df = pd.DataFrame.from_dict(datas)
    plt.figure(figsize=(10, 8))
    sns.jointplot(data=df, x='x_data', y='y_data', hue='path_name', kind='kde')
    plt.show()


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path to a gfa-like file.")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input styles", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Computes the number of times a node is crossed')
    args = parser.parse_args()

    assess_loop(args.file, args.gfa_version)
