"""
Here, we aim to extend the current GFA tag format by adding tags
that do respect the GFA naming convention.
A JSON string, PO (Path Offset) positions, relative to paths.
Hence, PO:J:{'w1':[(334,335,'+')],'w2':[(245,247,'-')]} tells that the walk/path w1
contains the sequence starting at position 334 and ending at position 335,
and the walk/path w2 contains the sequence starting at the offset 245 (ending 247),
and that the sequences are reversed one to each other.
Note that any non-referenced walk in this field means that the node
is not inside the given walk.
"""
from pgGraphs import Graph


def add_offsets_to_gfa(gfa_file: str, output_file: str) -> None:
    """Given a GFA file, calculates the offsets within paths and stores it as complementary tags

    Args:
        gfa_file (str): input gfa file
        output_file (str): output gfa file
        gfa_version (str): the user-assumed gfa subformat
    """
    gfa_graph: Graph = Graph(gfa_file, with_sequence=True)
    gfa_graph.sequence_offsets()
    gfa_graph.save_graph(output_file=output_file)
