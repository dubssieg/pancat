[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/python-3.12-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![https://tharos-ux.github.io/pangenome-notes/](https://img.shields.io/badge/docs-unfinished-orange.svg)]()

# PANCAT - PANgenome Comparison and Anlaysis Toolkit

> [!WARNING]\
> A paper is in preparation about this work. If you consider to use this tool, please contact the author for attribution.

Implementations of many functions for performing various actions on GFA-like graphs in a command-line tool, such as extracting or offseting a pangenome graph.
Is capable of comparing graphs topology between graphs that happen to contain the same set of sequences. Does pangenome graphs visualisation with interactive html files.
Uses the [gfagraphs library](https://pypi.org/project/gfagraphs/) to load and manipulate pangenome graphs.

![](https://media.discordapp.net/attachments/874430800802754623/1180182798968033280/graph_big.png)

> [!NOTE]\
> Want to contribute? Feel free to open a PR on an issue about a missing, buggy or incomplete feature! **Please do bug reports in the issue tracker!**.

## Installation

Requires **python $\geq$ 3.10**.

Installation can be made with the following command line, and updates may be run using `just` (requires [just](https://github.com/casey/just))

```bash
git clone https://github.com/Tharos-ux/pancat.git
cd pancat
pip install -r requirements.txt --upgrade
python -m pip install . --quiet
```

## Troubleshooting

> [!WARNING]\
> This tool is under heavy devlopment, and so it's [associated library](https://github.com/Tharos-ux/gfagraphs). I advise to update `pip install gfagraphs --upgrade` every now and then, when you update the tool. Any issue to this project is more than welcome, as I could not test all usecases! Feel free to [open one here](https://github.com/Tharos-ux/pancat/issues) if any problems occurs. **Please do bug reports in the issue tracker as well**.

## Quick start : provided commands

This program is a collection of tools.
Not every function or script is accessible through the front-end `pancat`, but this front-end showcase what the tools can do.
Other tools are in the `scripts` folder. 

Are available in `pancat`:

+ **compare** - Aligns shared genomes in graphs and compute the minimal set of events that explains how we go from one graph to another. A fastest RUST implementation is available [here](https://github.com/Tharos-ux/rs-pancat-compare).
+ **isolate** - Isolates a subgraph within a graph in-between two positions, according to a path or walk in the graph. Extraction only considers other paths and walks that cross both positions.
+ **offset** - Add path offsets to a graph. Adds a JSON string, PO (Path Offset) positions, relative to paths in the GFA file. Hence, PO:J:{'w1':[(334,336,'+')],'w2':[(245,247,'-'),(336,338,'-')]} tells that the walk/path w1 contains the sequence starting at position 334 and ending at position 336, and the walk/path w2 contains the sequence starting at the offset 245 (ending 247) and crosses it a second time between position 336 and 338, and that the sequences are reversed one to each other. Note that any non-referenced walk in this field means that the node is not inside the given walk.
+ **complete** - Asserts if the graph is a complete pangenome graph. A complete pangenome graph is a sequence graph where all the genomes used to build the graph are exactly contained in the graph, and that have paths.
+ **grapher** - Creates a html view of the graph. Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pancat isolate for instance), then computes the vizualisation on the selected part.
+ **multigrapher** - Creates a html view of the alignment of two graphs. Huge graphs may take long time to display, or might be messy. Advice would be to extract parts you want to display (with pancat isolate for instance), then computes the vizualisation on the selected part.
+ **stats** - Retrieves basic stats on a pangenome graph.
+ **reconstruct** - Reconstruct linear sequences from a GFA graph. Given a GFA file with paths, reconstruct a linear sequence for each haplotype between two offsets.
+ **unfold** - [WIP] Break cycles of a graph with paths or walks. Aims to linearize the graph by creating new nodes and connecting them in the graph.
+ **compress** - [WIP] Does a compression of the graph, merging simple non-conflicting substitution bubbles.

At any point, one can use `pancat -h` to get the manpage of the tool and its available subcommands, and also can display the help for a specific command (let's say **compare**) with `pancat compare -h`.
