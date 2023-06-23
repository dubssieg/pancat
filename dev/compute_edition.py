"Modelizes edition between paths and computes if edition is acceptable inside a graph"
from itertools import accumulate
from argparse import ArgumentParser, SUPPRESS
from os import path
from enum import Enum
from indexed import IndexedOrderedDict
from gfagraphs import Graph
from tharospytools import futures_collector


class EditEvent(Enum):
    "Relations between a couple of nodes"
    STR = "String"
    SUR = "SubString"
    SPR = "SuperString"
    SUP = "SubPrefix"
    SPP = "SuperPrefix"
    SUS = "SubSuffix"
    SPS = "SuperSuffix"
    OPS = "OverlapPrefixSuffix"
    OSP = "OverlapSuffixPrefix"


def edition_between_paths(reference_path: IndexedOrderedDict, query_path: IndexedOrderedDict) -> IndexedOrderedDict:
    """Analyses each relation between nodes of the paths, from left to right,
    and edits the dictionnary of editions containing a mapping ref_node:editions

    Raises:
        ValueError: if encountered event is invalid (bad index comparison, should not happend, but here in case)
    """
    edition: IndexedOrderedDict = IndexedOrderedDict()

    # Index for iterating in the paths
    idx_reference: int = 0
    idx_query: int = 0
    # Offsets to keep track of position inside alignment
    offset_reference: int = 0
    offset_query: int = 0
    # Counters for editions and relations
    required_merges: int = 0
    required_splits: int = 0

    # Computing endpos for each node
    reference_endpos: list = [0] + \
        list(accumulate(len(x) for x in reference_path.values()))
    query_endpos: list = [0] + \
        list(accumulate(len(x) for x in query_path.values()))

    while idx_reference < len(reference_path):

        # Init edit list for unbound
        edit_list: list = list()

        # Its only one event while we did not change both nodes
        current_pos: int = idx_reference
        # While we're on current node according to offset, we append edition to node
        while current_pos == idx_reference:

            reference_start: int = reference_endpos[idx_reference]
            reference_end: int = reference_endpos[idx_reference+1]

            query_start: int = query_endpos[idx_query]
            query_end: int = query_endpos[idx_query+1]

            if reference_end == query_end and reference_start == query_start:
                # No edition to be made, its an equivalence
                edit_list.append(EditEvent.STR)
                offset_reference += reference_end - reference_start
                offset_query += query_end - query_start

            elif reference_start > query_start and reference_end < query_end:
                # Reference is a subpart of other node
                edit_list.append(EditEvent.SPR)
                offset_reference += reference_end - reference_start
                offset_query += reference_end - reference_start
                required_splits += 1

            elif reference_start < query_start and reference_end > query_end:
                # Other node is included in reference
                edit_list.append(EditEvent.SUR)
                offset_reference += query_end - query_start
                offset_query += query_end - query_start

            elif reference_start == query_start and reference_end > query_end:
                edit_list.append(EditEvent.SUP)
                offset_reference += query_end - query_start
                offset_query += query_end - query_start

            elif reference_start == query_start and reference_end < query_end:
                edit_list.append(EditEvent.SPP)
                offset_reference += reference_end - reference_start
                offset_query += reference_end - reference_start
                required_splits += 1

            elif reference_start < query_start and reference_end == query_end:
                # Amogus
                edit_list.append(EditEvent.SUS)
                offset_reference += query_end - query_start
                offset_query += query_end - query_start

            elif reference_start > query_start and reference_end == query_end:
                edit_list.append(EditEvent.SPS)
                offset_reference += reference_end - reference_start
                offset_query += reference_end - reference_start

            elif reference_start < query_start and reference_end < query_end:
                edit_list.append(EditEvent.OPS)
                offset_reference += reference_end - query_start
                offset_query += reference_end - query_start
                required_splits += 1

            elif reference_start > query_start and reference_end > query_end:
                edit_list.append(EditEvent.OSP)
                offset_reference += query_end - reference_start
                offset_query += query_end - reference_start
            else:
                raise ValueError(
                    f"Unknown event ; Nodes pos={idx_reference}/{idx_query}, offset ref : {offset_reference}, offset query : {offset_query}\nRef start/end : {reference_start}/{reference_end}, Query start/end : {query_start}/{query_end}")

            # We add what we did treat to next potential event
            if offset_query >= query_endpos[idx_query+1]:
                idx_query += 1
            if offset_reference >= reference_endpos[idx_reference+1]:
                idx_reference += 1

        edition[reference_path.keys()[current_pos]] = edit_list
        required_merges += max(0, len(edit_list)-1)
        edit_list: list = list()
    return edition


def perform_edition(
        files: list,
        gfa_versions: list,
        output_folder: str,
) -> None:
    """Main loop for edition ; iterates over couples of files
    and performs edition on each couple if possible

    Args:
        files (list): a list of gfa files (in the same order as versions)
        gfa_versions (list): a list of gfa versions (in the same order as files)
        output_folder (str): a folder where to store results

    Raises:
        ValueError: files and gfa versions argument error
    """
    # In this function, we do calculate the distance between G1 and G2, by trying to modify G2 into G1.
    # Note the two graphs can be freely swapped, we just need to invert scores

    # Unpacking arguments
    file_1, file_2 = files
    version_1, version_2 = gfa_versions

    # Loading graphs
    graph1: Graph = Graph(file_1, version_1, True)
    graph2: Graph = Graph(file_2, version_2, True)

    paths_of_g1: dict = {walk.datas["name"]: walk.datas["path"]
                         for walk in graph1.get_path_list()}
    paths_of_g2: dict = {walk.datas["name"]: walk.datas["path"]
                         for walk in graph2.get_path_list()}

    union_of_paths: list[str] = list(set(
        list(paths_of_g1.keys())+list(paths_of_g2.keys())))

    for label in union_of_paths:

        pog1 = IndexedOrderedDict()
        pog2 = IndexedOrderedDict()

        for node, _ in paths_of_g1[label]:
            pog1[node] = graph1.get_segment(node).datas['seq']

        paths_of_g1[label] = pog1

        for node, _ in paths_of_g2[label]:
            pog2[node] = graph2.get_segment(node).datas['seq']

        paths_of_g2[label] = pog2

        if ''.join(pog1.values()) != ''.join(pog2.values()):
            raise ValueError("Paths does not represent the same genomes")

    results: list[IndexedOrderedDict] = futures_collector(
        edition_between_paths, [(paths_of_g1[label], paths_of_g2[label]) for label in union_of_paths])

    required_merges: int = 0
    required_splits: int = 0
    with open(path.join(output_folder, f"out.log"), "w", encoding='utf-8') as output:
        output.write(f"PATH_NAME\tNODE_ID\tRELATIONS\n")
        for i, result in enumerate(results):
            for node_id, relations in result.items():
                output.write(
                    f"{union_of_paths[i]}\t{node_id}\t{[rel.name for rel in relations]}\n")
                required_merges += len(relations)-1
                for relation in relations:
                    if relation == EditEvent.OPS or relation == EditEvent.SPP or relation == EditEvent.SPR:
                        required_splits += 1
        output.write("-"*60)
        output.write(
            f"\nTotal required merges: {required_merges}\nTotal required splits: {required_splits}")


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Paths to TWO gfa-like files.", nargs=2)
    parser.add_argument(
        "-o", "--output_folder", required=True, type=str, help="Path to a folder for results.")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the TWO corresponding GFA input styles", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs=2)
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another.')
    args = parser.parse_args()

    perform_edition(
        args.file,
        args.gfa_version,
        args.output_folder,
    )
