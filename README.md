[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![](https://img.shields.io/badge/wiki-nonexistent-red.svg)]()
[![](https://img.shields.io/badge/comments-done-green.svg)]()
[![](https://img.shields.io/badge/build-passing-green.svg)]()

# PANGRAPHS - GFA visualisation and exploration

Implementations of many functions for performing various actions on GFA-like graphs in a command-line tool, such as extracting or offseting a pangenome graph. Is capable of comparing graphs topology between graphs that happen to contain the same set of sequences. Does pangenome graphs visualisation with interactive html files.
Uses the [gfagraphs library](https://pypi.org/project/gfagraphs/) to load and manipulate pangenome graphs.

## Installation

Requires **python >=3.10**.

```bash
git clone git@github.com:Tharos-ux/pangraphs.git
cd pangraphs
python -m pip install . --quiet
```

## Quick start : provided commands

This program is a collection of tools. Not every function or script is accessible through the front-end `pangraphs`, but this front-end showcase what the tools can do.
Other tools are in the `scripts` folder. 

Are available through `pangraphs`:

- **grapher** creates interactive graph representation from a GFA file
- **reconstruct** recreates the linear sequences from the graph
- **offset** adds relative position information as a tag in GFA file
- **isolate** extracts a subgraph from positions in the paths
- **neigborhood** extracts a subgraph from a set of nodes around a node
- **edit** computes a edit distance between variation graphs