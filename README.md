[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/python-3.12-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![https://tharos-ux.github.io/pangenome-notes/](https://img.shields.io/badge/docs-in progress-orange.svg)]()

# PANCAT - PANgenome Comparison and Anlaysis Toolkit

Implementations of many functions for performing various actions on GFA-like graphs in a command-line tool, such as extracting or offseting a pangenome graph.
Is capable of comparing graphs topology between graphs that happen to contain the same set of sequences. Does pangenome graphs visualisation with interactive html files.
Uses the [gfagraphs library](https://pypi.org/project/gfagraphs/) to load and manipulate pangenome graphs.
Details about implementation can be [found here](https://hal.science/hal-04213245) (in french only, sorry).

![](https://media.discordapp.net/attachments/874430800802754623/1180182798968033280/graph_big.png)

> [!NOTE]\
> Want to contribute? Feel free to open a PR on an issue about a missing, buggy or incomplete feature!

## Installation

Requires **python $\geq$ 3.10**.

Installation can be made with the following command line, and updates may be run using `just` (requires [just](https://github.com/casey/just))

```bash
git clone https://github.com/Tharos-ux/pancat.git
cd pancat
pip install -r requirements.txt --upgrade
python -m pip install . --quiet
```

## Quick start : provided commands

This program is a collection of tools. Not every function or script is accessible through the front-end `pancat`, but this front-end showcase what the tools can do.
Other tools are in the `scripts` folder. 

Are available through `pancat`:

- **offset** adds relative position information as a tag in GFA file
- **grapher** creates interactive graph representation from a GFA file
- **stats** gathers basic stats from the input GFA 
- **complete** assesses if the graph is a complete pangenome graph (all genomes fully embedded in the graph)
- **reconstruct** recreates the linear sequences from the graph
- **edit** computes a edit distance between variation graphs

Were available before (and will be back soon):
- **isolate** extracts a subgraph from positions in the paths
- **neigborhood** extracts a subgraph from a set of nodes around a node
- **cycles** detect and (optionnally) linearizes all loops in graph

## Render interactive html view

With this command, you can create a html interactive view of your graph, with sequence in the nodes (S-lines) and nodes connected by edges (L-lines). If additional information is given (as such as W-lines or P-lines), supplementary edges will be drawn in order to show the path that the genomes follows in the graph.

```bash
pancat grapher [-h] [-b BOUNDARIES [BOUNDARIES ...]] file output

positional arguments:
  file                  Path to a gfa-like file
  output                Output path for the html graph file.

options:
  -h, --help            show this help message and exit
  -b BOUNDARIES [BOUNDARIES ...], --boundaries BOUNDARIES [BOUNDARIES ...]
                        One or a list of ints to use as boundaries for display (ex : -b 50 2000 will set 3 colors : one for nodes in range 0-50bp, one for nodes in range 51-2000 bp
                        and one for nodes in range 2001-inf bp).
```

When using this command, please only work with graphs with under 10k nodes. To do so, you may flatten the graph or extract subgraphs (using for instance **pancat neighborhood** or **pancat isolate**).

The `-b`/`--boundaries` option lets you choose size classes to differentiate. They will have a different color, and their number will be computed separately.

The `output` argument may be : a path to a folder (existing or not) or a path to a file (with .HTML extension or not).

## Compute stats on your graph

With this command, you can output basic stats on your graph.

```bash
pancat stats [-h] [-b BOUNDARIES [BOUNDARIES ...]] file

positional arguments:
  file                  Path to a gfa-like file

options:
  -h, --help            show this help message and exit
  -b BOUNDARIES [BOUNDARIES ...], --boundaries BOUNDARIES [BOUNDARIES ...]
                        One or a list of ints to use as boundaries for display (ex : -b 50 2000 will set 3 colors : one for nodes in range 0-50bp, one for nodes in range 51-2000 bp
                        and one for nodes in range 2001-inf bp).
```

This program displays stats in command-line (stdout). You may pipe it to a file if you want to use it on a cluster. (pancat stats graph.gfa > out.txt)

The `-b`/`--boundaries` option lets you choose size classes to differentiate. Their number will be computed separately.

## Extract sequences from the graph

With this command, you can reconstruct linear sequences from the graph.

```bash
pancat reconstruct [-h] -r REFERENCE [--start START] [--stop STOP] [-s] file out

positional arguments:
  file                  Path to a gfa-like file
  out                   Output path (without extension)

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Tells the reference sequence we seek start and stop into
  --start START         To specifiy a starting node on reference to create a subgraph
  --stop STOP           To specifiy a ending node on reference to create a subgraph
  -s, --split           Tells to split in different files
```

For this function, the `-r`/`--reference` option is needed only if you specify starting and ending points.

## Adding coordinate system

With this command, you ca add a JSON GFA-compatible string to each S-line of the graph (each node). This field will contain starting position, ending position and orientation, for each path in the graph.

```bash
pancat offset [-h] file out

positional arguments:
  file        Path to a gfa-like file
  out         Output path (with extension)

options:
  -h, --help  show this help message and exit
```

## Compute edition between graphs

In order to compare two graphs, they need to :
+ have at least some shared paths
+ the reconstruction of those shared paths must yield the same sequences

If those criteria are met, you may compare your graphs.

```bash
pancat edit [-h] -o OUTPUT_PATH [-g] [--selection [SELECTION ...]] graph_A graph_B

positional arguments:
  graph_A               Path to a GFA-like file.
  graph_B               Path to a GFA-like file.

options:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to a .json output for results.
  -g, --graph_level     Asks to perform edition computation at graph level.
  --selection [SELECTION ...]
                        Name(s) for the paths you want to reconstruct.
```