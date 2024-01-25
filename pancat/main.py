#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from json import dump
from pathlib import Path
from rich import print
from networkx import MultiDiGraph
from pgGraphs import Graph as pgGraph, GFANetwork
from tharospytools.path_tools import path_allocator
from tharospytools.list_tools import grouper
# TODO Done
from pancat.offset_in_gfa import add_offsets_to_gfa
from pancat.grapher import compute_stats, display_graph
from pancat.reconstruct_sequences import reconstruct_paths, graph_against_fasta
# TODO Testing
from pancat.edit_distance import perform_edition
# TODO Rebuilding
from pancat.find_bubbles import linearize_bubbles
from pancat.cyclotron import get_graph_cycles

from rich.traceback import install
install(show_locals=True)

parser: ArgumentParser = ArgumentParser(
    description='GFA manipulation tools.', add_help=True)
subparsers = parser.add_subparsers(
    help='Available subcommands', dest="subcommands")

parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'

## Subparser for isolate_by_range ##

"""
parser_isolate: ArgumentParser = subparsers.add_parser('isolate',help="Isolates a subgraph within a graph.\n""Relies on position in base pairs, requires the PO tag (built by pangraphs offset). In order to output a correct graph, you should provide a graph that has paths or walks to describe nodes chaining, and your range should be valid (meaning, base positions must be in the graph).")
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
"""

## Subparser for offset_in_gfa ##

parser_offset: ArgumentParser = subparsers.add_parser(
    'offset', help="Add path offsets to a graph\n"
    "Adds a JSON string, PO (Path Offset) positions, relative to paths. Hence, PO:J:{'w1':[(334,335,'+')],'w2':[(245,247,'-'),(336,338,'-')]} tells that the walk/path w1 contains the sequence starting at position 334 and ending at position 335, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247) and crosses it a second time between position 336 and 338, and that the sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk.")

parser_offset.add_argument("file", type=str, help="Path to a gfa-like file")
parser_offset.add_argument(
    "out", type=str, help="Output path (with extension)")

## Subparser for linearize ##

# parser_linear: ArgumentParser = subparsers.add_parser('linearize', help="Tries to simplify the graph structure, compressing superbubble chains in simple nodes.")
# parser_linear.add_argument("file", type=str, help="Path to a gfa-like file")
# parser_linear.add_argument("output", type=str,help="Output path for the gfa graph simplified file.")

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
    "--selection", type=str, help="Name(s) for the paths you want to reconstruct.", nargs='*', default=True)


## Subparser for cyclotron ##

# parser_cycles: ArgumentParser = subparsers.add_parser('cycles', help="Extracts all cycles from the graph file, path by path.\n""It can return those as lists of nodes names, sequences, or lengths of chains.")

# parser_cycles.add_argument("file", type=str, help="Path to a gfa-like file")
# parser_cycles.add_argument("output", type=str, help="Path to folder or json file")
# parser_cycles.add_argument("-m","--mode",help="Tells the return mode",required=True,choices=['names', 'sequences', 'lengths'])

## Subparser for edit_distance ##

parser_edit: ArgumentParser = subparsers.add_parser(
    'edit', help="Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another.")

parser_edit.add_argument(
    "graph_A", type=str, help="Path to a GFA-like file.")
parser_edit.add_argument(
    "graph_B", type=str, help="Path to a GFA-like file.")
parser_edit.add_argument(
    "-o", "--output_path", required=True, type=str, help="Path to a .json output for results.")
parser_edit.add_argument(
    "-p", "--pattern", type=str, help="Regexp to filer on for path names. Dealults to no filtering.")
parser_edit.add_argument(
    "-g", "--graph_level", help="Asks to perform edition computation at graph level.", action='store_true', default=False)
parser_edit.add_argument(
    "-c", "--cores", help="Number of cores for computing edition", type=int, default=1)
parser_edit.add_argument(
    "-s", "--selection", type=str, help="Name(s) for the paths you want to compute edition on.", nargs='*', default=True)
parser_edit.add_argument(
    "-t", "--trace_memory", help="Print to log file memory usage of data structures.", action='store_true', default=False)

## Subparser for unfold ##

parser_edit: ArgumentParser = subparsers.add_parser(
    'unfold', help="Break cycles of a graph with P-lines or W-lines.")

parser_edit.add_argument(
    "graph", type=str, help="Path to a >=GFA1-like file.")
parser_edit.add_argument(
    "-o", "--output_path", required=True, type=str, help="Path to a .gfa output for new graph.")


#######################################

args = parser.parse_args()


def main() -> None:
    "Main call for subprograms"
    print("PANCAT initialisation sucessful!")
    if len(argv) == 1:
        print(
            "[dark_orange]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit()

    ##############################################################################
    #                     AVAILABLE & WORKING COMMANDS                           #
    ##############################################################################

    if args.subcommands == 'offset':
        "This command adds PO offsets to the graph, and then saves it."
        add_offsets_to_gfa(
            gfa_file=args.file,
            output_file=args.out
        )

    elif args.subcommands == 'complete':
        "This command assess that the graph is a complete pangenome graph."
        graph_against_fasta(
            gfa_graph=args.file,
            pipeline_txt=args.pipeline
        )

    elif args.subcommands == 'grapher':
        "This command aims to render a graph as a PyVis network displayed in html file"
        bounds = grouper(
            [0] + [bound+x for bound in args.boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
        )

        # Creating pgGraphs object
        gfa_graph: pgGraph = pgGraph(
            gfa_file=args.file,
            with_sequence=True
        )

        # Computing NetworkX structure
        pangenome_graph: MultiDiGraph = GFANetwork.compute_networkx(
            graph=gfa_graph,
            node_prefix=None,
            node_size_classes=bounds
        )

        # Computing stats and displaying
        graph_stats = compute_stats(
            graph=gfa_graph,
            length_classes=tuple(bounds)
        )
        display_graph(
            graph=pangenome_graph,
            colors_paths=gfa_graph.metadata['colors'],
            annotations=graph_stats,
            output_path=args.output
        )

    elif args.subcommands == 'stats':
        "This command aims to extract stats from a graph for basic analysis"
        bounds = grouper(
            [0] + [bound+x for bound in args.boundaries for x in [0, 1]] + [float('inf')], n=2, m=1
        )

        # Creating pgGraphs object
        gfa_graph: pgGraph = pgGraph(
            gfa_file=args.file,
            with_sequence=True
        )

        # Computing stats and displaying
        graph_stats = compute_stats(
            graph=gfa_graph,
            length_classes=tuple(bounds)
        )
        for key, value in graph_stats.items():
            print(f"{key}: {value}")

    elif args.subcommands == 'reconstruct':
        "This command reconstruct sequences from the graph"
        sequences: dict = reconstruct_paths(
            gfa_file=args.file,
            selected_paths=args.selection
        )
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

    ##############################################################################
    #                        REQUIRES IN-DEPTH TESING                            #
    ##############################################################################

    elif args.subcommands == 'unfold':
        graph_to_unfold: pgGraph = pgGraph(
            gfa_file=args.graph,
            with_sequence=True
        )
        graph_to_unfold.unfold()
        graph_to_unfold.save_graph(
            output_file=path_allocator(
                path_to_validate=args.output_path,
                particle='.gfa',
                default_name='unfolded_graph',
                always_yes=True
            )
        )

    elif args.subcommands == 'edit':
        perform_edition(
            gfa_A=args.graph_A,
            gfa_B=args.graph_B,
            output_path=args.output_path,
            graph_level=args.graph_level,
            selection=args.selection,
            cores=args.cores,
            regular_expression=args.pattern,
            trace_memory=args.trace_memory,
        )

    ##############################################################################
    #                      BUGGY & DE-ACTIVATED FOR NOW                          #
    ##############################################################################

    elif args.subcommands == 'isolate':
        """
        range_isolate(
            gfa_file=args.file,
            output=args.out,
            start=args.start,
            stop=args.end,
            reference_name=args.reference)
        """
        raise NotImplementedError()

    elif args.subcommands == 'linearize':
        "TODO: better detection of cycles (buggy for now) / already done by GFAffix in part"
        linearize_bubbles(
            gfa_file=args.file,
            output=args.output
        )

    elif args.subcommands == 'cycles':
        "TODO: better detection of cycles (buggy for now)"
        dump(
            get_graph_cycles(
                args.file,
                args.mode
            ),
            open(
                path_allocator(
                    args.output,
                    particle='.json',
                    default_name=f"cycles_{Path(args.file).stem}"
                ),
                mode='w',
                encoding='utf-8'
            ),
            indent=4
        )

    else:
        print(
            "[dark_orange]Unknown command. Please use the help to see available commands.")
        exit(1)
