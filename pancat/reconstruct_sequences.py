'Extract sequences given a graph'
from Bio import SeqIO
from gfagraphs import Graph
from os.path import join
from pathlib import Path


def reconstruct_paths(gfa_file: str, folder: str, selected_paths: list | bool = False, split: bool = True) -> dict:
    """Given a GFA file with paths, follows those paths to reconstruct the linear sequences

    Args:
        gfa_file (str): gfa graph file, with paths
        gfa_version (str): the version of the file

    Raises:
        ValueError: if no path in graph

    Returns:
        dict: mapping between path name and sequence
    """
    gfa_graph: Graph = Graph(
        gfa_file=gfa_file,
        with_sequence=True
    )
    if len(gfa_graph.paths) > 0:
        genomes: dict = gfa_graph.reconstruct_sequences()
        selected_paths: list = (selected_paths, sorted(list(gfa_graph.paths.keys())))[
            isinstance(selected_paths, bool)]
        for label in selected_paths:
            with open(join(folder, f'{label}.fasta') if split else join(folder, f'{Path(gfa_file).stem}.fasta'), 'w' if split else 'a', encoding="utf-8") as writer:
                writer.write(
                    f">{label}\n{''.join(genomes[label])}\n"
                )
                del genomes[label]
    else:
        raise ValueError(
            "Graph has no embed paths. Could'nt determine what to reconstruct!"
        )


def graph_against_fasta(gfa_graph: str, pipeline_txt: str) -> bool:
    with open(pipeline_txt, 'r', encoding='utf-8') as pipeline:
        lines = pipeline.readlines()

        # We load the graph
        gfa_graph: Graph = Graph(
            gfa_file=gfa_graph,
            with_sequence=True
        )
        if len(gfa_graph.paths) > 0:
            genomes: dict = {
                x.lower(): y for x, y in gfa_graph.reconstruct_sequences().items()
            }

        print(f"Graph paths: \t{', '.join(genomes.keys())}")
        print(f"Sequence headers: \t{', '.join(genomes.keys())}\n")

        complete_pangenome_graph: list[bool] = [
            False for _ in range(len(lines))]
        for y, line in enumerate(lines):
            # We assume the pipeline file is in the cactus format
            seq_name, seq_path = line.strip().split(
                '\t')[0].lower(), line.strip().split('\t')[1]
            if seq_name in genomes:
                # Do comparison
                # We load the file
                with open(seq_path, 'r', encoding='utf-8') as reader:
                    for record in SeqIO.parse(reader, "fasta"):
                        seq1: str = record.seq.lower()
                        seq2: str = ''.join(genomes[seq_name]).lower()
                        if len(seq1) != len(seq2):
                            print(
                                f"[{seq_name}] Length of sequence with header {record.id} does not match the path {seq_name} from a size-perspective (resp. {len(seq1)} vs {len(seq2)}).")
                        else:
                            densities: list = [0 for _ in range(101)]
                            for i in range(len(seq1)):
                                if seq1[i] != seq2[i]:
                                    densities[round((i/len(seq1))*100)] += 1
                            if sum(densities) != 0:
                                print(
                                    f"[{seq_name}] Sequence {record.id} and path have the same length ({len(seq1)}). {round((sum(densities)/len(seq1))*100,2)}% of base pairs are not shared.")
                            else:
                                print(
                                    f"[{seq_name}] Sequence {record.id} and path represents the same sequence.")
                                complete_pangenome_graph[y] = True
            else:
                print(
                    f"[{seq_name}] Can't find {seq_name} in the paths of the graph.")
        if all(complete_pangenome_graph):
            print(f"\n{gfa_graph} is a complete pangenome graph")
        else:
            print(f"\n{gfa_graph} is not a complete pangenome graph")
    return all(complete_pangenome_graph)


def graph_against_multifasta(gfa_graph: str, pipeline_txt: str) -> bool:

    with open(pipeline_txt, 'r', encoding='utf-8') as pipeline:
        lines = pipeline.readlines()

        # We load the graph
        gfa_graph: Graph = Graph(
            gfa_file=gfa_graph,
            with_sequence=True
        )
        if len(gfa_graph.paths) > 0:
            genomes: dict = {
                x.lower(): y for x, y in gfa_graph.reconstruct_sequences().items()
            }

        print(f"Graph paths: \t{', '.join(genomes.keys())}")
        print(f"Sequence headers: \t{', '.join(genomes.keys())}\n")

        complete_pangenome_graph: dict[str, bool] = {
            path_name: False for path_name in gfa_graph.paths.keys()
        }
        for y, line in enumerate(lines):
            # We assume the pipeline file is in the cactus format
            seq_name, seq_path = line.strip().split(
                '\t')[0].lower(), line.strip().split('\t')[1]

            with open(seq_path, 'r', encoding='utf-8') as reader:
                for record in SeqIO.parse(reader, "fasta"):
                    if candidate := gfa_graph.compare_pathnames_to_string(record.id):
                        seq1: str = record.seq.lower()
                        seq2: str = ''.join(genomes[candidate]).lower()
                        if len(seq1) != len(seq2):
                            print(
                                f"[{candidate}] Length of sequence with header {record.id} does not match the path {candidate} from a size-perspective (resp. {len(seq1)} vs {len(seq2)}).")
                        else:
                            densities: list = [0 for _ in range(101)]
                            for i in range(len(seq1)):
                                if seq1[i] != seq2[i]:
                                    densities[round((i/len(seq1))*100)] += 1
                            if sum(densities) != 0:
                                print(
                                    f"[{candidate}] Sequence {record.id} and path have the same length ({len(seq1)}). {round((sum(densities)/len(seq1))*100,2)}% of base pairs are not shared.")
                            else:
                                print(
                                    f"[{candidate}] Sequence {record.id} and path represents the same sequence.")
                                complete_pangenome_graph[y] = True
                    else:
                        print(
                            f"[{record.id}] Can't find {record.id} in the paths of the graph."
                        )
        for path, status in complete_pangenome_graph:
            if not status:
                print(
                    f"[{path}] Missing {path} in fasta files.")
        if all(complete_pangenome_graph.values()):
            print(f"\n{gfa_graph} is a complete pangenome graph")
        else:
            print(f"\n{gfa_graph} is not a complete pangenome graph")
    return all(complete_pangenome_graph.values())
