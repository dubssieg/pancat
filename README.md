[![](https://img.shields.io/badge/Python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/Documentation-unfinished-orange.svg)]()
[![](https://img.shields.io/badge/build-passing-green.svg)]()

# PANGRAPHS - GFA visualisation and exploration

Does pangenome graphs visualisation (sort of) and many scripts for performing various actions on such graphs.

## Installation

Requires **python >=3.10**.

You can install pangraphs by running `python setup.py install`, which will enlable command-line usage of the program.

For a manual install, note that dependancies are in `requirements.txt`.
Install all with `pip install -r requirements.txt`.
A conda env is provided in the `env.sh` script.

## Quick start : provided commands

This tool is a collection of small scripts. Not every function or script is accessible through the front-end `pangraphs`, but this front-end showcase what the tools can do.

```text
usage: pangraphs [-h] {isolate,offset,neighborhood,scaffold,grapher,levenshtein,compare,convert,length,vcfmatch,reconstruct,align} ...

GFA manipulation tools.

Subcommands:
  {isolate,offset,neighborhood,scaffold,grapher,levenshtein,compare,convert,length,vcfmatch,reconstruct,align}
                        Available subcommands
    isolate             Isolates a subgraph within a graph. Relies on position in base pairs, requires the PO tag (built by pangraphs offset). In order to output a correct graph, you should provide
                        a graph that has paths or walks to describe nodes chaining, and your range should be valid (meaning, base positions must be in the graph).
    offset              Add path offsets to a graph Adds a JSON string, PO (Path Offset) positions, relative to paths. Hence, PO:J:{'w1':(334,335,'+'),'w2':(245,247,'-')} tells that the walk/path w1
                        contains the sequence starting at position 334 and ending at position 335, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247), and that the
                        sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk.
    neighborhood        Extracts subgraph given a starting node Given a node and a number of neighbors, attempts to extract paths, nodes and edges around the target node. Beware : if you select a
                        node at one of the ends of the graph, you may stuck yourself in a infinite loop.
    scaffold            Cuts fasta file to isolate chromosoms/scaffolds from PAF file. Extracts from a fasta-like file all sequences in a query assembly given a mapping to a reference and an
                        identifier on reference.
    grapher             Creates a html view of the graph. Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pangraphs isolate
                        or pangraphs neighborhood for instance), then computes the vizualisation on the selected part.
    levenshtein         Evaluates identity of nodes within and across graphs. Given multiple graphs, aims to compute Levenshtein distance between nodes in order to evaluate identity. This
                        implementation is quite slow, and we moved to a quicker and more informative comparaison method (pangraphs compare) but we keep it here as it's not exactly the same
                        information.
    compare             Does position-based checks of segment status between graphs, following paths. For each path, tries to evaluate, based on position, the existence of shifts, inclusions and
                        equivalences between graphs using the same set of coordinates.
    convert             (Experimental) Attempts to convert rGFA to GFA1. Converts rGFA files issued from minigraph to GFA1 compatible format. It implies to rename nodes and add P-lines if asked for.
                        As the P-lines could not be precisely defined from the rGFA, and at is does not realign sequences on the graph, any ambiguous path will go through reference sequence.
    length              Plot distribution of nodes lengths across graph.
    vcfmatch            Maps variants to graph. Given a VCF file and a GFA graph, evaluates as a heatmap alignement score based on Levenshtein distance.
    reconstruct         Reconstruct linear sequences from a GFA graph. Given a GFA file with paths, reconstruct a linear sequence for each haplotype between two offsets.
    align               Verifies to which sequences are mapped each node of a GFA, and where. Two figures are produced : a dotgrid displaying a haplotype/segment mapping, and a alignment where
                        segments are matched back on the linear genome.
```
