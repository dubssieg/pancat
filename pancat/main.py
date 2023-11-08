#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from os.path import exists
from os import stat
from json import dump
from pathlib import Path
from rich import print
from networkx import MultiDiGraph
from gfagraphs import Graph as pGraph, supplementary_datas
from workspace.isolate_by_range import range_isolate
from workspace.offset_in_gfa import add_offsets_to_gfa
from workspace.grapher import compute_stats, display_graph
from workspace.reconstruct_sequences import reconstruct_paths, graph_against_fasta
from workspace.edit_distance import perform_edition
from workspace.find_bubbles import linearize_bubbles
from workspace.cyclotron import get_graph_cycles
from tharospytools.path_tools import path_allocator

from rich.traceback import install
install(show_locals=True)

parser: ArgumentParser = ArgumentParser(
    description='GFA manipulation tools.', add_help=True)
subparsers = parser.add_subparsers(
    help='Available subcommands', dest="subcommands")

parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'

## Subparser for isolate_by_range ##

parser_isolate: ArgumentParser = subparsers.add_parser(
    'isolate',
    help="Isolates a subgraph within a graph.\n"
    "Relies on position in base pairs, requires the PO tag (built by pangraphs offset). In order to output a correct graph, you should provide a graph that has paths or walks to describe nodes chaining, and your range should be valid (meaning, base positions must be in the graph)."
)

parser_isolate.add_argument("file", type=str, help="Path to a gfa-like file")
parser_isolate.add_argument(
    "out", type=str, help="Output path (with extension)")
parser_isolate.add_argument(
    '-s',
    '--start',
    type=int,
    help='To specifiy a starting point (in bp) to create a subgraph',
)
parser_isolate.add_argument(
    '-e',
    '--end',
    type=int,
    help='To specifiy a end point (in bp) to create a subgraph',
)
parser_isolate.add_argument(
    '-r',
    '--reference',
    type=str,
    help='To specifiy the path to follow',
)

## Subparser for offset_in_gfa ##

parser_offset: ArgumentParser = subparsers.add_parser(
    'offset', help="Add path offsets to a graph\n"
    "Adds a JSON string, PO (Path Offset) positions, relative to paths. Hence, PO:J:{'w1':[(334,335,'+')],'w2':[(245,247,'-'),(336,338,'-')]} tells that the walk/path w1 contains the sequence starting at position 334 and ending at position 335, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247) and crosses it a second time between position 336 and 338, and that the sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk.")

parser_offset.add_argument("file", type=str, help="Path to a gfa-like file")
parser_offset.add_argument(
    "out", type=str, help="Output path (with extension)")

## Subparser for linearize ##

parser_linear: ArgumentParser = subparsers.add_parser(
    'linearize', help="Tries to simplify the graph structure, compressing superbubble chains in simple nodes.")

parser_linear.add_argument("file", type=str, help="Path to a gfa-like file")
parser_linear.add_argument("output", type=str,
                           help="Output path for the gfa graph simplified file.")

## Subparser for grapher ##

parser_grapher: ArgumentParser = subparsers.add_parser(
    'grapher', help="Creates a html view of the graph.\n"
    "Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pangraphs isolate or pangraphs neighborhood for instance), then computes the vizualisation on the selected part.")

parser_grapher.add_argument("file", type=str, help="Path to a gfa-like file")
parser_grapher.add_argument("output", type=str,
                            help="Output path for the html graph file.")
parser_grapher.add_argument(
    "-b", "--boundaries", type=int, help="One or a list of ints to use as boundaries for display (ex : -b 50 2000 will set 3 colors : one for nodes in range 0-50bp, one for nodes in range 51-2000 bp and one for nodes in range 2001-inf bp).", nargs='+', default=[50])

## Subparser for stats ##

parser_stats: ArgumentParser = subparsers.add_parser(
    'stats', help="Retrieves basic stats on a pangenome graph.")

parser_stats.add_argument("file", type=str, help="Path to a gfa-like file")
parser_stats.add_argument(
    "-b", "--boundaries", type=int, help="One or a list of ints to use as boundaries for display (ex : -b 50 2000 will set 3 colors : one for nodes in range 0-50bp, one for nodes in range 51-2000 bp and one for nodes in range 2001-inf bp).", nargs='+', default=[50])

## Subparser for assert_complete_pangenome ##

parser_complete: ArgumentParser = subparsers.add_parser(
    'complete', help="Asserts if the graph is a completet pangenome graph"
)

parser_complete.add_argument(
    "file", type=str, help="Path to a gfa-like file")
parser_complete.add_argument(
    "pipeline", type=str, help="Tab-separated mapping between path names and path to files")


## Subparser for reconstruct_sequences ##

parser_reconstruct: ArgumentParser = subparsers.add_parser(
    'reconstruct', help="Reconstruct linear sequences from a GFA graph.\n"
    "Given a GFA file with paths, reconstruct a linear sequence for each haplotype between two offsets."
)

parser_reconstruct.add_argument(
    "file", type=str, help="Path to a gfa-like file")
parser_reconstruct.add_argument(
    "out", type=str, help="Output path (without extension)")
parser_reconstruct.add_argument(
    "-s", "--split", help="Tells to split in different files", action='store_true')
parser_reconstruct.add_argument(
    "--selection", type=str, help="Name(s) for the paths you want to reconstruct.", nargs='*', default=None)


## Subparser for cyclotron ##

parser_cycles: ArgumentParser = subparsers.add_parser(
    'cycles', help="Extracts all cycles from the graph file, path by path.\n"
    "It can return those as lists of nodes names, sequences, or lengths of chains."
)

parser_cycles.add_argument(
    "file", type=str, help="Path to a gfa-like file")
parser_cycles.add_argument(
    "output", type=str, help="Path to folder or json file")
parser_cycles.add_argument(
    "-m",
    "--mode",
    help="Tells the return mode",
    required=True,
    choices=['names', 'sequences', 'lengths']
)

## Subparser for edit_distance ##

parser_edit: ArgumentParser = subparsers.add_parser(
    'edit', help="Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another.")

parser_edit.add_argument(
    "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
parser_edit.add_argument(
    "-o", "--output_folder", required=True, type=str, help="Path to a folder for results.")
parser_edit.add_argument(
    "-p", "--perform_edition", help="Asks to perform edition on graph and outputs it.", action='store_true')
parser_edit.add_argument(
    "--selection", type=str, help="Name(s) for the paths you want to reconstruct.", nargs='*', default=None)


#######################################

args = parser.parse_args()


def get_gfa_subtype(gfa_file_path: str | list[str]) -> str | list[str]:
    """Given a file, or more, returns the gfa subtypes, and raises error if file is invalid or does not exists

    Args:
        gfa_file_path (str | list[str]): one or more file paths

    Returns:
        str | list[str]: a gfa subtype descriptor per input file
    """
    styles: list[str] = list()
    if isinstance(gfa_file_path, str):
        gfa_file_path = [gfa_file_path]
    for gfa_file in gfa_file_path:
        # Checking if path exists
        if not exists(gfa_file):
            raise OSError(
                "Specified file does not exists. Please check provided path."
            )
        # Checking if file descriptor is valid
        if not gfa_file.endswith('.gfa'):
            raise IOError(
                "File descriptor is invalid. Please check format, this lib is designed to work with Graphical Fragment Assembly (GFA) files."
            )
        # Checking if file is not empty
        if stat(gfa_file).st_size == 0:
            raise IOError(
                "File is empty."
            )
        with open(gfa_file, 'r', encoding='utf-8') as gfa_reader:
            header: str = gfa_reader.readline()
            if header[0] != 'H':
                styles.append('rGFA')
            else:
                try:
                    version_number: str = supplementary_datas(
                        header.strip('\n').split('\t'), 1
                    )["VN"]
                    if version_number == '1.0':
                        styles.append('GFA1')
                    elif version_number == '1.1':
                        styles.append('GFA1.1')
                    elif version_number == '1.2':
                        styles.append('GFA1.2')
                    elif version_number == '2.0':
                        styles.append('GFA2')
                    else:
                        styles.append('unknown')
                except KeyError:
                    styles.append('rGFA')
    if len(styles) == 1:
        return styles[0]
    return styles


def main() -> None:
    "Main call for subprograms"
    if len(argv) == 1:
        print(
            "[dark_orange]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit()
    gfa_version_info: str | list[str] = get_gfa_subtype(args.file)
    print(
        f"Working on file(s) {args.file} of respective detected formats {gfa_version_info}")
    if args.subcommands == 'isolate':
        range_isolate(
            gfa_file=args.file,
            gfa_ver=gfa_version_info,
            output=args.out,
            start=args.start,
            stop=args.end,
            reference_name=args.reference
        )
    elif args.subcommands == 'cycles':
        dump(get_graph_cycles(args.file, gfa_version_info, args.mode),
             open(path_allocator(args.output, particle='.json', default_name=f"cycles_{Path(args.file).stem}"), mode='w', encoding='utf-8'), indent=4)
    elif args.subcommands == 'complete':
        graph_against_fasta(args.file, gfa_version_info, args.pipeline)
    elif args.subcommands == 'offset':
        add_offsets_to_gfa(args.file, args.out, gfa_version_info)
    elif args.subcommands == 'stats':
        pgraph: pGraph = pGraph(
            args.file, gfa_version_info, with_sequence=True)
        bounds: list = []
        boundaries = [
            0] + [bound+x for bound in args.boundaries for x in [0, 1]] + [float('inf')]
        for i in range(0, len(boundaries), 2):
            x = i
            bounds.append([boundaries[x], boundaries[x+1]])

        graph_stats = compute_stats(pgraph, length_classes=tuple(bounds))

        for key, value in graph_stats.items():
            print(f"{key}: {value}")
    elif args.subcommands == 'grapher':
        bounds: list = []
        boundaries = [
            0] + [bound+x for bound in args.boundaries for x in [0, 1]] + [float('inf')]
        for i in range(0, len(boundaries), 2):
            x = i
            bounds.append([boundaries[x], boundaries[x+1]])

        pangenome_graph: MultiDiGraph = (pgraph := pGraph(
            args.file, gfa_version_info, with_sequence=True)).compute_networkx(node_size_classes=tuple(bounds))

        graph_stats = compute_stats(pgraph, length_classes=tuple(bounds))
        display_graph(
            graph=pangenome_graph,
            colors_paths=pgraph.colors,
            annotations=graph_stats,
            output_path=args.output
        )
    elif args.subcommands == 'reconstruct':
        sequences: dict = reconstruct_paths(
            args.file, gfa_version_info, args.selection)
        if args.split:
            for i, (label, sequence) in enumerate(sequences.items()):
                with open(f"{args.out}_{i}.fasta", "w", encoding="utf-8") as writer:
                    writer.write(
                        f">{label}\n{''.join(sequence)}\n"
                    )
        else:
            with open(f"{args.out}.fasta", "w", encoding="utf-8") as writer:
                for label, sequence in sequences.items():
                    writer.write(
                        f">{label}\n{''.join(sequence)}\n"
                    )
    elif args.subcommands == 'edit':
        perform_edition(args.file, gfa_version_info,
                        args.output_folder, args.perform_edition, args.selection)
    elif args.subcommands == 'linearize':
        linearize_bubbles(
            gfa_file=args.file,
            gfa_type=gfa_version_info,
            output=args.output
        )
    else:
        print(
            "[dark_orange]Unknown command. Please use the help to see available commands.")
        exit(1)
