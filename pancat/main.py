#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from json import dump
from os.path import exists
from pathlib import Path
from rich import print
from pgGraphs import Graph as pgGraph
from tharospytools.path_tools import path_allocator
from pancat.constants import *
# TODO Done
from pancat.offset_in_gfa import add_offsets_to_gfa
from pancat.grapher import graph_stats, graph_viewer, multigraph_viewer
from pancat.reconstruct_sequences import reconstruct_paths, graph_against_fasta, graph_against_multifasta
# TODO Testing
from pancat.edit_distance import perform_edition
from pancat.compress_graph import compress_graph
from pancat.correct import correct_graph
from pancat.isolate_subgraph import extract_subgraph
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

parser_isolate: ArgumentParser = subparsers.add_parser(
    'isolate',
    help=HELP_COMMAND_ISOLATE,
)
parser_isolate.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_isolate.add_argument(
    "out",
    type=str,
    help=HELP_OUTPUT_FILE_GFA,
)
parser_isolate.add_argument(
    '-r',
    '--reference',
    type=str,
    help=HELP_INPUT_REFERENCE,
)
parser_isolate.add_argument(
    '-s',
    '--start',
    type=int,
    help=HELP_INPUT_START_POSITION,
)
parser_isolate.add_argument(
    '-e',
    '--end',
    type=int,
    help=HELP_INPUT_END_POSITION,
)

## Subparser for offset_in_gfa ##

parser_offset: ArgumentParser = subparsers.add_parser(
    'offset',
    help=HELP_COMMAND_OFFSET,
)
parser_offset.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_offset.add_argument(
    "out",
    type=str,
    help=HELP_OUTPUT_FILE_GFA,
)

## Subparser for grapher ##

parser_grapher: ArgumentParser = subparsers.add_parser(
    'grapher',
    help=HELP_COMMAND_GRAPHER,
)
parser_grapher.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_GENERAL,
)
parser_grapher.add_argument(
    "output",
    type=str,
    help=HELP_OUTPUT_FILE_HTML,
)
parser_grapher.add_argument(
    "-b",
    "--boundaries",
    type=int,
    help=HELP_PARAM_BOUNDARIES,
    nargs='+',
    default=[50],
)

## Subparser for multigrapher ##

parser_multigrapher: ArgumentParser = subparsers.add_parser(
    'multigrapher',
    help=HELP_COMMAND_MULTIGRAPHER,
)
parser_multigrapher.add_argument(
    "file_A",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_multigrapher.add_argument(
    "file_B",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_multigrapher.add_argument(
    "editions",
    type=str,
    help=HELP_INPUT_EDITION_FILE,
)
parser_multigrapher.add_argument(
    "output",
    type=str,
    help=HELP_OUTPUT_FILE_HTML,
)
parser_multigrapher.add_argument(
    "-b",
    "--boundaries",
    type=int,
    help=HELP_PARAM_BOUNDARIES,
    nargs='+',
    default=[50],
)

## Subparser for stats ##

parser_stats: ArgumentParser = subparsers.add_parser(
    'stats',
    help=HELP_COMMAND_STATS,
)
parser_stats.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_GENERAL,
)
parser_stats.add_argument(
    "-b",
    "--boundaries",
    type=int,
    help=HELP_PARAM_BOUNDARIES,
    nargs='+',
    default=[50],
)

## Subparser for assert_complete_pangenome ##

parser_complete: ArgumentParser = subparsers.add_parser(
    'complete',
    help=HELP_COMMAND_COMPLETE,
)
parser_complete.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_complete.add_argument(
    "pipeline",
    type=str,
    help=HELP_INPUT_PIPELINE,
)
parser_complete.add_argument(
    "-m",
    "--multifasta_mode",
    help=HELP_PARAM_MULTIFASTA,
    action='store_true',
    default=False,
)

## Subparser for reconstruct_sequences ##

parser_reconstruct: ArgumentParser = subparsers.add_parser(
    'reconstruct',
    help=HELP_COMMAND_RECONSTRUCT,
)
parser_reconstruct.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_reconstruct.add_argument(
    "out",
    type=str,
    help=HELP_OUTPUT_FOLDER,
)
parser_reconstruct.add_argument(
    "-s",
    "--split",
    help=HELP_PARAM_SPLIT,
    action='store_true',
)
parser_reconstruct.add_argument(
    "--selection",
    type=str,
    help=HELP_PARAM_PATHSELECT,
    nargs='*',
    default=True,
)

## Subparser for edit_distance ##

parser_edit: ArgumentParser = subparsers.add_parser(
    'edit',
    help=HELP_COMMAND_EDIT,
)
parser_edit.add_argument(
    "graph_A",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_edit.add_argument(
    "graph_B",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_edit.add_argument(
    "-o",
    "--output_path",
    required=True,
    type=str,
    help=HELP_OUTPUT_JSON_LOG,
)
parser_edit.add_argument(
    "-p",
    "--pattern",
    type=str,
    help=HELP_PARAM_REGEXP,
)
parser_edit.add_argument(
    "-g",
    "--graph_level",
    help=HELP_PARAM_GRAPH_LEVEL,
    action='store_true',
    default=False,
)
parser_edit.add_argument(
    "-c",
    "--cores",
    help=HELP_PARAM_THREADS,
    type=int,
    default=1,
)
parser_edit.add_argument(
    "-s",
    "--selection",
    type=str,
    help=HELP_PARAM_PATHSELECT,
    nargs='*',
    default=True,
)
parser_edit.add_argument(
    "-t",
    "--trace_memory",
    help=HELP_PARAM_VERBOSE,
    action='store_true',
    default=False,
)


## Subparser for compress_graph ##

parser_compress: ArgumentParser = subparsers.add_parser(
    'compress',
    help=HELP_COMMAND_COMPRESS,
)
parser_compress.add_argument(
    "file",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_compress.add_argument(
    "-o",
    "--output_file",
    required=True,
    type=str,
    help=HELP_OUTPUT_FILE_GFA,
)
parser_compress.add_argument(
    "-m",
    "--minimize",
    help=HELP_PARAM_MINIMIZE,
    action='store_true',
    default=False,
)
parser_compress.add_argument(
    "-l",
    "--length",
    help=HELP_PARAM_LENGTH,
    type=int,
    default=float('inf'),
)

## Subparser for unfold ##

parser_edit: ArgumentParser = subparsers.add_parser(
    'unfold',
    help=HELP_COMMAND_UNFOLD,
)
parser_edit.add_argument(
    "graph",
    type=str,
    help=HELP_INPUT_FILE_WITH_PATHS,
)
parser_edit.add_argument(
    "-o",
    "--output_path",
    required=True,
    type=str,
    help=HELP_OUTPUT_FILE_GFA,
)


#######################################

args = parser.parse_args()


def main() -> None:
    "Main call for subprograms"
    if len(argv) == 1:
        print(
            "[dark_orange]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands."
        )
        exit()

    for identifier, syspath in [(key, path) for key, path in args.__dict__.items() if key in ['file', 'graph_A', 'graph_B', 'file_A', 'file_B']]:
        if not exists(syspath):
            raise RuntimeError(
                f"Specified path '{syspath}' for argument '{identifier}' does not exists."
            )

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
        if args.multifasta_mode:
            graph_against_multifasta(
                gfa_graph=args.file,
                pipeline_txt=args.pipeline
            )
        else:
            graph_against_fasta(
                gfa_graph=args.file,
                pipeline_txt=args.pipeline
            )

    elif args.subcommands == 'grapher':
        "This command aims to render a graph as a PyVis network displayed in html file"
        graph_viewer(
            file=args.file,
            output=args.output,
            boundaries=args.boundaries
        )

    elif args.subcommands == 'multigrapher':
        "This command aims to render the alignment of two graphs"
        multigraph_viewer(
            file_A=args.file_A,
            file_B=args.file_B,
            file_editions=args.editions,
            boundaries=args.boundaries,
            output=args.output,
        )

    elif args.subcommands == 'stats':
        "This command aims to extract stats from a graph for basic analysis"
        graph_stats(
            file=args.file,
            boundaries=args.boundaries
        )

    elif args.subcommands == 'reconstruct':
        "This command reconstruct sequences from the graph"
        reconstruct_paths(
            gfa_file=args.file,
            folder=args.out,
            selected_paths=args.selection,
            split=args.split,
        )

    ##############################################################################
    #                        REQUIRES IN-DEPTH TESING                            #
    ##############################################################################

    elif args.subcommands == 'compress':
        compress_graph(
            args.file,
            args.output_file,
            minimized=args.minimize,
            max_len_to_collapse=args.length,
        )

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
        extract_subgraph(
            gfa_path=args.file,
            x=args.start,
            y=args.end,
            sequence=args.reference.upper(),
            output=args.out
        )

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
