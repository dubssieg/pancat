#!/usr/bin/env python3
from argparse import ArgumentParser
from sys import argv
from networkx import MultiDiGraph
from BubbleGun.Graph import Graph as bGraph
from gfagraphs import Graph as pGraph
from tharospytools import get_palette
from scripts.isolate_by_range import isolate
from scripts.offset_in_gfa import add_offsets_to_gfa
from scripts.paths_bubblegun_bfs import bfs_step, paths_step
from scripts.parse_genomes import isolate_scaffolds
from scripts.grapher import html_graph
from scripts.levenshtein_distance import show_identity, display_graph
from scripts.compare_by_offset import display_graph as compare_display_graph, get_backbone, compare_positions
from scripts.gfa_convert import rgfa_to_gfa
from scripts.length_ratios import parse_gfa, plot_ratio
from scripts.vcf_on_graph import vcf_heatmap, match_nodes_to_vcf, vcf_parser
from scripts.reconstruct_sequences import reconstruct, node_range, grab_paths
from scripts.map_graph_nodes import mapper, dotgrid_plot, show_alignments, exact_mapper, get_lengths
from scripts.create_vcf import render_vcf, get_graph_structure
from rich import print

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
    "-g",
    "--gfa_version",
    help="Tells the GFA input style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
)
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
    "Adds a JSON string, PO (Path Offset) positions, relative to paths. Hence, PO:J:{'w1':(334,335,'+'),'w2':(245,247,'-')} tells that the walk/path w1 contains the sequence starting at position 334 and ending at position 335, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247), and that the sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk.")

parser_offset.add_argument("file", type=str, help="Path to a gfa-like file")
parser_offset.add_argument(
    "out", type=str, help="Output path (with extension)")
parser_offset.add_argument(
    "-g",
    "--gfa_version",
    help="Tells the GFA input style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
)

## Subparser for paths_bubblegun_bfs ##

parser_neighborhood: ArgumentParser = subparsers.add_parser(
    'neighborhood', help="Extracts subgraph given a starting node\n"
    "Given a node and a number of neighbors, attempts to extract paths, nodes and edges around the target node. Beware : if you select a node at one of the ends of the graph, you may stuck yourself in a infinite loop.")

parser_neighborhood.add_argument(
    "file", type=str, help="Path to a gfa-like file")
parser_neighborhood.add_argument(
    "out", type=str, help="Output path (with extension)")
parser_neighborhood.add_argument(
    "-g",
    "--gfa_version",
    help="Tells the GFA input style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
)
parser_neighborhood.add_argument(
    "-o",
    "--gfa_output",
    help="Tells the GFA output style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
)
parser_neighborhood.add_argument(
    '-s',
    '--start_node',
    type=str,
    help='To specifiy a starting node on reference to create a subgraph',
    nargs='+'
)
parser_neighborhood.add_argument("-c", "--count", type=int,
                                 help="Number of nodes around each starting point")

## Subparser for parse_genomes ##


parser_scaffold: ArgumentParser = subparsers.add_parser(
    'scaffold', help="Cuts fasta file to isolate chromosoms/scaffolds from PAF file.\n"
    "Extracts from a fasta-like file all sequences in a query assembly given a mapping to a reference and an identifier on reference.")
parser_scaffold.add_argument(
    "file", type=str, help="fasta-like file")
parser_scaffold.add_argument(
    "out", type=str, help="fasta-like output")
parser_scaffold.add_argument(
    "paffile", type=str, help="paf-like file")
parser_scaffold.add_argument(
    "chromosom", type=str, help="name of assembly on reference sequence")

## Subparser for grapher ##

parser_grapher: ArgumentParser = subparsers.add_parser(
    'grapher', help="Creates a html view of the graph.\n"
    "Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pangraphs isolate or pangraphs neighborhood for instance), then computes the vizualisation on the selected part.")

parser_grapher.add_argument("file", type=str, help="Path to a gfa-like file")
parser_grapher.add_argument("job_name", type=str,
                            help="Job identifier for output (ex : chr3_graph)")
parser_grapher.add_argument(
    "-n", "--number_alignments", help="Gives the number of origin sequences", required=True, type=int)
parser_grapher.add_argument(
    "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])

## Subparser for levenshtein_distance ##

parser_levenshtein: ArgumentParser = subparsers.add_parser(
    'levenshtein', help="Evaluates identity of nodes within and across graphs.\n"
    "Given multiple graphs, aims to compute Levenshtein distance between nodes in order to evaluate identity. This implementation is quite slow, and we moved to a quicker and more informative comparaison method (pangraphs compare) but we keep it here as it's not exactly the same information."
)

parser_levenshtein.add_argument(
    "file", type=str, help="Path(s) to one or more gfa-like file(s).", nargs='+')
parser_levenshtein.add_argument(
    "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
parser_levenshtein.add_argument("-s", "--score", help="Filters by percentage of identity",
                                required=False, type=float, default=None)
parser_levenshtein.add_argument(
    "-b", "--backbone", help="Asks to hide paths within graphs.", action='store_true')

## Subparser for compare_by_offset ##

parser_compare: ArgumentParser = subparsers.add_parser(
    'compare', help="Does position-based checks of segment status between graphs, following paths.\n"
    "For each path, tries to evaluate, based on position, the existence of shifts, inclusions and equivalences between graphs using the same set of coordinates."
)

parser_compare.add_argument(
    "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
parser_compare.add_argument("-j", "--job_name", type=str, required=True,
                            help="Job identifier for output (ex : chr3_graph)")
parser_compare.add_argument("-r", "--reference", type=str, required=True,
                            help="Path to refer to if position is ambiguous.")
parser_compare.add_argument(
    "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
parser_compare.add_argument(
    "-p", "--paths", type=str, help="Path(s) to display.", nargs='+', default=None)
parser_compare.add_argument(
    "-s", "--with_sequences", help="Asks to show full sequences in node info.", action='store_true')


## Subparser for gfa_convert ##

parser_convert: ArgumentParser = subparsers.add_parser(
    'convert', help="(Experimental) Attempts to convert rGFA to GFA1.\n"
    "Converts rGFA files issued from minigraph to GFA1 compatible format. It implies to rename nodes and add P-lines if asked for. As the P-lines could not be precisely defined from the rGFA, and at is does not realign sequences on the graph, any ambiguous path will go through reference sequence."
)

parser_convert.add_argument("file", type=str, help="rGFA file")

parser_convert.add_argument(
    "-p", "--plines", help="Asks to calculate p-lines for graph", action='store_true')
parser_convert.add_argument(
    "-k", "--keep", help="Keeps rGFA-specific tags in output", action='store_true')
parser_convert.add_argument(
    "-n",
    "--haplotypes_names",
    help="Give one name per haplotype, ordered as the same order you did for your rGFA. First item of the list is the reference. Required if you ask for P-lines.",
    nargs='+',
    default=[]
)

## Subparser for length_ratios ##

parser_length: ArgumentParser = subparsers.add_parser(
    'length', help="Plot distribution of nodes lengths across graph.\n"
    ""
)

parser_length.add_argument("file", type=str, help="gfa-like file", nargs='+')
parser_length.add_argument("-x", "--xmax", type=int,
                           help="maximum length size to plot. recommanded : 10000.", required=True)
parser_length.add_argument(
    "-g",
    "--gfa_version",
    help="Tells the GFA input style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'],
    nargs='+'
)

## Subparser for length_ratios ##

parser_vcfmatch: ArgumentParser = subparsers.add_parser(
    'vcfmatch', help="Maps variants to graph.\n"
    "Given a VCF file and a GFA graph, evaluates as a heatmap alignement score based on Levenshtein distance.")

parser_vcfmatch.add_argument(
    "gfa", type=str, help="Path to a gfa-like file.")
parser_vcfmatch.add_argument(
    "vcf", type=str, help="Path to a vcf-like file.")
parser_vcfmatch.add_argument(
    "-s", "--score", help="Filters by percentage of identity", required=False, type=float, default=None)

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
    "-g",
    "--gfa_version",
    help="Tells the GFA input style",
    required=True,
    choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2']
)
parser_reconstruct.add_argument(
    "-r",
    "--reference",
    help="Tells the reference sequence we seek start and stop into",
    required=True,
    type=str
)
parser_reconstruct.add_argument(
    '--start',
    type=str,
    help='To specifiy a starting node on reference to create a subgraph',
    default=None
)
parser_reconstruct.add_argument(
    '--stop',
    type=str,
    help='To specifiy a ending node on reference to create a subgraph',
    default=None
)
parser_reconstruct.add_argument(
    "-s", "--split", help="Tells to split in different files", action='store_true')

## Subparser for map_graph_nodes ##

parser_align: ArgumentParser = subparsers.add_parser(
    'align', help="Verifies to which sequences are mapped each node of a GFA, and where. Two figures are produced : a dotgrid displaying a haplotype/segment mapping, and a alignment where segments are matched back on the linear genome."
)
parser_align.add_argument(
    "gfa", type=str, help="gfa-like file")
parser_align.add_argument(
    "fasta", type=str, help="fasta-like file")
parser_align.add_argument(
    "-n", "--highlight", type=str, help="nodes to highlight in another color", nargs='+')

## Subparser for create_vcf ##

parser_vcf: ArgumentParser = subparsers.add_parser(
    'vcf', help="(Experimental) Aims to extract variants from a GFA graph. Calls bubbles inside the graph and proceeds to check if node is in reference or not. If not, interprets it as a variant, and seeks a sequence it refers too. Requires PO offset (pangraphs offset) on your GFA input file.")
parser_vcf.add_argument("file", type=str, help="Input GFA file")
parser_vcf.add_argument("output", type=str,
                        help="Output path for VCF (with extension)")
parser_vcf.add_argument(
    "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'])

parser_vcf.add_argument(
    "-c", "--chromosom", help="Name of the chromosom the graph is from", required=True, type=str)
parser_vcf.add_argument(
    "-r", "--reference", help="Name of the reference path inside graph", required=True, type=str)


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
    if args.subcommands == 'isolate':
        isolate(args.file, args.out, args.start, args.end,
                args.gfa_version, args.reference)
    elif args.subcommands == 'offset':
        add_offsets_to_gfa(args.file, args.out, args.gfa_version)
    elif args.subcommands == 'neighborhood':
        graph: bGraph = bGraph(args.file)
        for i, node in enumerate(args.start_node):
            output: str = f"{args.out.split('.')[0]}_{i}.gfa" if len(
                args.start_node) > 1 else args.out
            nodes: set = bfs_step(graph, node, args.count)
            paths_step(args.file, output, nodes,
                       args.gfa_version, args.gfa_output)
    elif args.subcommands == 'scaffold':
        isolate_scaffolds(fasta_file=args.file, out_file=args.out,
                          paf_file=args.paffile, chromosom=args.chromosom)
    elif args.subcommands == 'grapher':
        pangenome_graph: MultiDiGraph = pGraph(
            args.file,
            args.gfa_version
        ).compute_networkx()

        html_graph(
            pangenome_graph,
            args.job_name
        )
    elif args.subcommands == 'levenshtein':
        display_graph(
            show_identity(
                args.file,
                args.gfa_version,
                get_palette(len(args.file)),
                args.score,
                args.backbone
            )
        )
    elif args.subcommands == 'compare':
        if len(args.file) < 2:
            parser.error("Please specify at least two GFA files as input.")
        if len(args.file) != len(args.gfa_version):
            parser.error(
                "Please match the number of args between files and gfa types.")

        for i, (path, nodes, graphs, colors) in enumerate(get_backbone(args.file, args.gfa_version, args.with_sequences)):
            datas, full_graph = compare_positions(f"{args.job_name}_{i}" if i > 0 else f"{args.job_name}",
                                                  path, nodes, graphs, args.reference, path_names=args.paths, shifts=True)  # type: ignore
            compare_display_graph(
                full_graph, f"{args.job_name}_{i}" if i > 0 else f"{args.job_name}", args.paths, datas, [args.file[i], args.file[i+1]], colors)
    elif args.subcommands == 'convert':
        rgfa_to_gfa(
            args.file,
            f"{args.file.split('.')[0]}_gfa1.gfa",
            names_of_haplotypes=args.haplotypes_names,
            p_lines=args.plines,
            keep_tags=args.keep)
    elif args.subcommands == 'length':
        file_names: list = [filepath.split('.')[0].split('/')[-1] for filepath in args.file] if isinstance(
            args.file, list) else [args.file.split('.')[0].split('/')[-1]]
        lengths: list = [parse_gfa(name) for name in args.file]
        plot_ratio(lengths, file_names, 1, args.xmax, 1)
    elif args.subcommands == 'vcfmatch':
        vcf_heatmap(match_nodes_to_vcf(
            args.gfa, vcf_parser(args.vcf), args.score))
    elif args.subcommands == 'reconstruct':
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
    elif args.subcommands == 'align':
        fmap = mapper(args.gfa, args.fasta)
        dotgrid_plot(fmap, nodes_to_highlight=args.highlight)

        show_alignments(exact_mapper(args.gfa, args.fasta),
                        get_lengths(args.fasta), nodes_to_highlight=args.highlight)
    elif args.subcommands == 'vcf':
        render_vcf(args.output, get_graph_structure(
            args.file, args.gfa_version, args.reference, args.chromosom))
    else:
        print(
            "[dark_orange]Unknown command. Please use the help to see available commands.")
        exit(1)
