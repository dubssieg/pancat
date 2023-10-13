import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
from os import listdir, path
from collections import Counter
from workspace.edit_distance import perform_edition

if __name__ == "__main__":

    graph_folder: str = "/home/sidubois/Workspace/Notes/graphs_to_compare/graphs"
    output_folder: str = "/home/sidubois/Workspace/Notes/graphs_to_compare/output"
    categories: str = ['String', 'SuperString', 'SubString', 'SubSuffix', 'SuperPrefix',
                       'SubPrefix', 'SuperSuffix', 'OverlapSuffixPrefix', 'OverlapPrefixSuffix']
    all_graphs: list = listdir(graph_folder)
    # fig, axs = plt.subplots(ncols=len(categories))

    table: dict[str, list[list[int]]] = {category: [
        [0 for _ in range(len(all_graphs))] for _ in range(len(all_graphs))] for category in categories}

    for i, graph_1 in enumerate(all_graphs):
        for j, graph_2 in enumerate(all_graphs):
            results: Counter = perform_edition([path.join(graph_folder, graph_1), path.join(
                graph_folder, graph_2)], ['GFA1.1', 'GFA1.1'], output_folder, do_edition=False)
            for category in categories:
                table[category][i][j] = results[category]

    for i, (category, datas) in enumerate(table.items()):
        # axs[i].set_title(category)
        fig, axs = plt.figure(figsize=(6, 8))
        plt.title(category)
        ax = sns.heatmap(datas, linewidth=0.5, annot=True,
                         cbar=False, ax=axs)  # ax=axs[i]
        plt.savefig(
            f"/home/sidubois/Workspace/Notes/graphs_to_compare/output/graph_{category}.png")
        plt.clf()
