"Plots graph statistics"
from re import sub
from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from statsmodels.stats.weightstats import DescrStatsW
from networkx import MultiDiGraph, is_isolate
import matplotlib.pyplot as plt


def size_ratio(counts: list, threshold: int) -> float:
    """Computes the ratio of nodes below or equal a threshold / nodes above a threshold
    1 means all nodes are above threshold

    Args:
        counts (list): a list of tuples containing (length:number)
        threshold (int): a limit between nodes considered as 'small' and 'big'

    Returns:
        float: the propotion of nodes with length beneath threshold
    """
    num, denom = .0, .0
    for size, count in counts:
        if size <= threshold:
            num += count
        else:
            denom += count
    return num/(denom+num)


def get_palette(number_of_colors: int) -> list:
    """Returns a number_of_colors-sized palette, as a list,
    that one can access with colors[i].

    Args:
        number_of_colors (int): number of colors needed

    Returns:
        list: palette of colors
    """
    colormap = plt.cm.viridis  # type:ignore
    number_of_colors = min(colormap.N, number_of_colors)
    return [colormap(int(x*colormap.N/number_of_colors)) for x in range(number_of_colors)]


def plot_ratio(counts: list, names: list, x_min: int, x_max: int, x_step: int = 1) -> None:
    """Given a set of counts of segments frequencies,
    displays the ratio of length and their number givent their length.

    Args:
        counts (list): a list of tuples containing (length:number)
        names (list): names of files to annotate the graph
        x_min (int): length lower bound for graph
        x_max (int): length upper bound for graph
        x_step (int, optional): Step between two ratios calculation. Defaults to 1.
    """
    mapcolors: list = get_palette(len(counts))
    _, ax1 = plt.subplots(figsize=(20, 4), sharex=True)
    ax2 = ax1.twinx()
    ratio_min: float = 1.0
    for i, count in enumerate(counts):
        ratios: list = [size_ratio(count, x)
                        for x in range(x_min, x_max, x_step)]
        ax1.plot([i for i in range(x_min, x_max, x_step)], ratios,
                 label=names[i], color=mapcolors[i])
        ax2.plot([l for l, _ in count], [c for _, c in count],
                 alpha=0.5, color=mapcolors[i], linestyle=':')
        ratio_min = min(ratio_min, min(ratios))
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax1.legend(loc='upper center', bbox_to_anchor=(
        0.5, -0.05), ncol=len(counts))
    ax1.set_ylim(ratio_min, 1)
    ax2.set_ylim(1, max([c for count in counts for (_, c) in count]))
    plt.xlim(x_min, x_max)
    plt.savefig('test.png', bbox_inches='tight')
    plt.show()


def lonely_nodes(graph: MultiDiGraph) -> list:
    """Those are some of the saddest nodes. :(

    Args:
        graph (MultiDiGraph): _description_

    Returns:
        list: _description_
    """
    return sorted(
        [(k, v) for k, v in Counter(
            [
                int(sub('\D', '', datas['title'])) for node, datas in graph.nodes(data=True)
                if is_isolate(graph, node)
            ]
        ).items()]
    )


def neighboured_nodes(graph: MultiDiGraph) -> list:
    """Those are not lonely.

    Args:
        graph (MultiDiGraph): _description_

    Returns:
        list: _description_
    """
    return sorted(
        [(k, v) for k, v in Counter(
            [
                int(sub('\D', '', datas['title'])) for node, datas in graph.nodes(data=True)
                if not is_isolate(graph, node)
            ]
        ).items()]
    )


def parse_gfa(input_file: str) -> list:
    """Counts the length distribution of all segments in a gfa file

    Args:
        input_file (str): path to a GFA-like file

    Returns:
        Counter: length distribution of sequences
    """
    with open(input_file, 'r', encoding='utf-8') as gfa_reader:
        return sorted([(k, v) for k, v in Counter([len(seq.split()[2]) for seq in gfa_reader if seq.split()[0] == 'S']).items()])


def plot_distribution(counts: list[list[tuple]], lonely: list[list[tuple]] | None, graph_names: list) -> None:
    """Given a Counter of length distribution, plots ths distribution

    Args:
        counts (list): list of sorted Counter, length distribution of segments of graph
    """
    fig, axs = plt.subplots(nrows=len(counts), ncols=1,
                            sharex=True, figsize=(10, 12))
    for i in range(len(counts)):
        x: list = [k for (k, _) in counts[i]]
        y: list = [v for (_, v) in counts[i]]
        w_stats: DescrStatsW = DescrStatsW(x, weights=y, ddof=0)
        if lonely is not None:
            a: list = [k for (k, _) in lonely[i]]
            b: list = [v for (_, v) in lonely[i]]
        try:
            axs[i].set_title(
                f"Distribution for {graph_names[i]}, |V|={sum(y)}, $\mu$={round(w_stats.mean,ndigits=2)}, $\sigma$={round(w_stats.std,ndigits=2)}")
            axs[i].plot(x, y)
            axs[i].plot(a, b)
            axs[i].set_yscale('log')
        except TypeError:
            axs.set_title(
                f"Distribution for {graph_names[i]}, |V|={sum(y)}, $\mu$={round(w_stats.mean,ndigits=2)}, $\sigma$={round(w_stats.std,ndigits=2)}")
            axs.plot(x, y)
            axs.plot(a, b)
            axs.set_yscale('log')
    try:
        for ax in axs.flat:
            ax.set(xlabel='Sequence size', ylabel='Number of occurences')
            ax.label_outer()
    except AttributeError:
        axs.set(xlabel='Sequence size', ylabel='Number of occurences')
    plt.xscale('log')
    fig.savefig('test.png', bbox_inches='tight')
    plt.show()


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument("file", type=str, help="gfa-like file", nargs='+')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Plot distribution of node length across graph')
    parser.add_argument(
        "-g",
        "--gfa_version",
        help="Tells the GFA input style",
        required=True,
        choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'],
        nargs='+'
    )
    args = parser.parse_args()

    file_names: list = [filepath.split('.')[0].split('/')[-1] for filepath in args.file] if isinstance(
        args.file, list) else [args.file.split('.')[0].split('/')[-1]]

    lengths: list = [parse_gfa(name) for name in args.file]
    print(lengths)

    plot_ratio(lengths, file_names, 1, 10000, 1)
