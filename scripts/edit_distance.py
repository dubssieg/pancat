"Modelizes edition between paths and computes if edition is acceptable inside a graph"
from itertools import combinations, accumulate
from argparse import ArgumentParser, SUPPRESS
from indexed import IndexedOrderedDict
from gfagraphs import Graph
from tharospytools import revcomp


class EditEvent():
    "Modelizes a elementary operation between a couple of nodes"

    def __init__(self, seq_a: str, seq_b: str, offset: int, start: int, end: int) -> None:
        """_summary_

        Args:
            seq_a (str): sequence of reference
            seq_b (str): sequence of query
            offset (int): the offset inside the query sequence
            start (int): the start position inside the segment (0=first base)
            end (int): the end position inside the segment
        """
        self.event = type(self).__name__
        self.revcomp: bool = revcomp(seq_b) == seq_a
        self.positions: tuple = (start, end)
        self.offsets: tuple = (offset, offset+(end-start))

    def __str__(self) -> str:
        return f"{self.event}(revcomp={self.revcomp}, positions={self.positions})"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, __o: object) -> bool:
        return isinstance(__o, type(self)) and self.positions == __o.positions and self.revcomp == __o.revcomp


class String(EditEvent):
    "Situation of two equal nodes"


class Overlap_Prefix_Prefix(EditEvent):
    "Particular kind of inclusion"


class Overlap_Suffix_Suffix(EditEvent):
    "Particular case of inclusion"


class Overlap_Prefix_Suffix(EditEvent):
    "Prefix/Suffix means edit = Split"


class Overlap_Suffix_Prefix(EditEvent):
    "Prefix/suffix means edit = split"


class Substring(EditEvent):
    "Affix succession means that edit = merge"


class Superstring(EditEvent):
    "Cut on both ends"


class PathEdit:
    "Defines edition between two paths"

    def __init__(self, path_name: str, path_a: IndexedOrderedDict, path_b: IndexedOrderedDict):
        """Creates a PathEdit object

        Args:
            path_name (str): a user-friendly name for the PathEdit object (eg. name of the haplotype)
            path_a (IndexedOrderedDict): a dict mapping node names to sequences for the first path
                ! (keys:values) pairs MUST BE ORDERED as they are traversed by the path!
            path_b (IndexedOrderedDict): same format, path we compare the first one to

        Raises:
            ValueError: if reconstructed sequences are different between the two paths
        """
        self.alignment_name: str = path_name
        # Storing paths (association of ordered segment names and segment sequences)
        self.reference_path: IndexedOrderedDict = path_a
        self.query_path: IndexedOrderedDict = path_b
        # Basic evaluation of path comparison
        self.same_segment_count: bool = len(path_a) == len(path_b)
        # Can be revcomp or same seq, or revcomp of anysubseq
        self.can_be_aligned: bool = ''.join(path_a.values()) == ''.join(
            path_b.values()) or ''.join(path_a.values()) == revcomp(''.join(path_b.values()))
        if not self.can_be_aligned:
            raise ValueError(
                "Could not compute edition between the paths as they don't describe the same sequences.")
        self.edition: IndexedOrderedDict = IndexedOrderedDict()
        self.compute_edition()

    def compute_edition(self):
        "We give a list of strings that are the segments in the order the path traverses them"
        # Clearing previous edition calculations
        self.edition: IndexedOrderedDict = IndexedOrderedDict()
        # Index for iterating in the paths
        idx_reference: int = 0
        idx_query: int = 0
        # Offsets to keep track of position inside alignment
        offset_reference: int = 0
        offset_query: int = 0

        # Computing endpos for each node
        reference_endpos: list = [0] + \
            list(accumulate(len(x) for x in self.reference_path.values()))
        query_endpos: list = [0] + \
            list(accumulate(len(x) for x in self.query_path.values()))

        while idx_reference < len(self.reference_path):
            # Current targeted nodes
            seq_reference: str = self.reference_path.values()[idx_reference]
            seq_query: str = self.query_path.values()[idx_query]

            # Init edit list for unbound
            edit_list: list = list()

            # Its only one event while we did not change both nodes
            current_pos: int = idx_reference
            while current_pos == idx_reference:

                reference_start: int = reference_endpos[idx_reference]
                reference_end: int = reference_endpos[idx_reference+1]

                query_start: int = query_endpos[idx_query]
                query_end: int = query_endpos[idx_query+1]

                if reference_end == query_end and reference_start == query_start:
                    # Segments are the same
                    # No edition to be made, its an equivalence
                    # (Its also a specific case of inclusion)
                    edit_list.append(
                        String(seq_reference, seq_query, offset_query, 0, query_end-query_start))
                    offset_reference += reference_end - reference_start
                    offset_query += query_end - query_start

                elif reference_start > query_start and reference_end < query_end:
                    # Reference is a subpart of other node
                    edit_list.append(Superstring(
                        seq_reference, seq_query, offset_query, reference_start - query_start, reference_end - reference_start))
                    offset_reference += reference_end - reference_start
                    offset_query += reference_end - reference_start

                elif reference_start < query_start and reference_end > query_end:
                    # Other node is included in reference
                    edit_list.append(
                        Substring(seq_reference, seq_query, offset_query, 0, query_end-query_start))
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start == query_start and reference_end > query_end:
                    edit_list.append(
                        Overlap_Prefix_Prefix(seq_reference, seq_query, offset_query, 0, query_end-query_start))
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start < query_start and reference_end == query_end:
                    edit_list.append(
                        Overlap_Suffix_Suffix(seq_reference, seq_query, offset_query, 0, query_end-query_start))
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start < query_start and reference_end < query_end:
                    edit_list.append(
                        Overlap_Prefix_Suffix(seq_reference, seq_query, offset_query, 0, reference_end-query_start))
                    offset_reference += reference_end - query_start
                    offset_query += reference_end - query_start

                else:
                    edit_list.append(Overlap_Suffix_Prefix(
                        seq_reference, seq_query, offset_query, reference_start-offset_query, query_end-query_start))
                    offset_reference += query_end - reference_start
                    offset_query += query_end - reference_start

                # We add what we did treat to next potential event

                if offset_query >= query_endpos[idx_query+1]:
                    idx_query += 1
                if offset_reference >= reference_endpos[idx_reference+1]:
                    idx_reference += 1

            self.edition[self.reference_path.keys()[current_pos]] = edit_list
            edit_list: list = list()


def evaluate_consensus(*paths_editions: PathEdit) -> bool:
    """Given a series of PathEdit, evaluates if they are equal.

    Args:
        paths (PathEdit): any number of PathEdit objects

    Returns:
        bool: if paths edits are non ambiguous
    """
    for path_A, path_B in combinations(paths_editions, 2):
        shared_nodes = set(path_A.edition.keys()).intersection(
            set(path_B.edition.keys()))
        if not (path_A.can_be_aligned and path_B.can_be_aligned and all([path_A.edition[x] == path_B.edition[x] for x in shared_nodes])):
            return False
    return True

################################################ EDIT ######################################


class Edition():
    pass


class Merge(Edition):
    "This operation means all its child elements must be merged, left to right"

    def __init__(self, *edit_events: Edition) -> None:
        self.events: list = list(edit_events)

    def __str__(self) -> str:
        return f"{type(self).__name__}({', '.join(self.events)})"

    def __repr__(self) -> str:
        return self.__str__()


class Split(Edition):
    "This operation means a node should be splitted in two"


class Keep(Edition):
    "This operation means that all child elements must not be changed"


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Does position-based checks of segment edition between graphs, by following paths.')
    args = parser.parse_args()

    # Handling input errors
    if len(args.file) < 2:
        parser.error("Please specify at least two GFA files as input.")
    if len(args.file) != len(args.gfa_version):
        parser.error(
            "Please match the number of args between files and gfa types.")

    for (f1, v1), (f2, v2) in combinations(zip(args.file, args.gfa_version), 2):
        all_dipaths: list[PathEdit] = list()

        graph1: Graph = Graph(f1, v1, True)
        graph2: Graph = Graph(f2, v2, True)

        paths = set(path.datas['name'] for path in graph1.get_path_list()).intersection(
            set(path.datas['name'] for path in graph2.get_path_list()))

        for path in paths:
            path_in_g1 = graph1.get_path(path)
            path_in_g2 = graph2.get_path(path)

            pog1 = IndexedOrderedDict()
            pog2 = IndexedOrderedDict()

            for node, _ in path_in_g1.datas['path']:
                pog1[node] = graph1.get_segment(node).datas['seq']

            for node, _ in path_in_g2.datas['path']:
                pog2[node] = graph2.get_segment(node).datas['seq']

            all_dipaths.append(PathEdit(path, pog1, pog2))

        del graph1, graph2
        with open("out.log", "w", encoding='utf-8') as output:
            for dipath in all_dipaths:
                for key, value in dipath.edition.items():
                    if not isinstance(value[0], String):
                        output.write(
                            f"{dipath.alignment_name}\t{key}\t{value}\n")

        print(evaluate_consensus(*all_dipaths))
