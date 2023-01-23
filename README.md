# pangraphs
Does pangenome graphs (sort of)

## Installation

Requires python >=3.10
Dependancies are in `requirements.txt`. Install all with `pip install -r requirements.txt`.
Alternatively, a conda env is provided in the `env.sh` script.

## To create a pangenome graph

Graphs are interactive .html files.

```bash
usage: grapher.py [-h] -n NUMBER_ALIGNMENTS -g {rGFA,GFA1,GFA1.1,GFA1.2,GFA2} file job_name

positional arguments:
  file                  Path to a gfa-like file
  job_name              Job identifier for output (ex : chr3_graph)

options:
  -h, --help            Creates a representation of a pangenome graph.
  -n NUMBER_ALIGNMENTS, --number_alignments NUMBER_ALIGNMENTS
                        Please put the number of origin sequences (n# of genomes/scaffolds depending usecase)
  -g {rGFA,GFA1,GFA1.1,GFA1.2,GFA2}, --gfa_version {rGFA,GFA1,GFA1.1,GFA1.2,GFA2}
                        Tells the GFA input style
```

### Example usage : 

```bash
python grapher.py toy_examples/pggb_graph.gfa toy_examples/pggb_view -n 3 -g GFA1

python grapher.py toy_examples/cactus_graph.gfa toy_examples/cactus_view -n 3 -g GFA1.1
```