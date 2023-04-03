"Modelizes edition between paths and computes if edition is acceptable inside a graph"
from itertools import combinations, accumulate, chain
from argparse import ArgumentParser, SUPPRESS
from os import path
from pathlib import Path
from copy import deepcopy
from indexed import IndexedOrderedDict
from rich import print
from gfagraphs import Graph
from tharospytools import revcomp


class EditEvent():
    "Modelizes a elementary operation between a couple of nodes"

    def __init__(
            self,
            ref_name: str,
            target_name: str,
            seq_a: str,
            seq_b: str,
            offset_target: int,
            start_target: int,
            end_target: int,
            offset_ref: int,
            start_ref: int,
            end_ref: int,
            need_split: bool
    ) -> None:
        """Inits a EditEvent object

        Args:
            seq_a (str): sequence of reference
            seq_b (str): sequence of query
            offset (int): the offset inside the query sequence
            start (int): the start position inside the segment (0=first base)
            end (int): the end position inside the segment
        """
        self.reference_node = ref_name
        self.target_node = target_name
        self.event = type(self).__name__
        self.revcomp: bool = revcomp(seq_b) == seq_a
        self.positions_target: tuple = (start_target, end_target)
        self.offsets_target: tuple = (
            offset_target, offset_target+(end_target-start_target))
        self.positions_ref: tuple = (start_ref, end_ref)
        self.offsets_ref: tuple = (offset_ref, offset_ref+(end_ref-start_ref))
        self.need_split: bool = need_split
        self.relative_pos: tuple = (0, end_ref-start_ref)
        self.length: int = end_ref-start_ref
        self.length_alt: int = end_target-start_target

    def __str__(self) -> str:
        """Ovverides the default str

        Returns:
            str: a description of the EditEvent
        """
        return f"{self.event}(revcomp={self.revcomp}, positions={self.offsets_target}, refnode={self.reference_node}, qnode={self.target_node})"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, __o: object) -> bool:
        return isinstance(__o, type(self)) and self.positions_target == __o.positions_target and self.revcomp == __o.revcomp


class String(EditEvent):
    "Situation of two equal nodes"


class SubPrefix(EditEvent):
    "Particular kind of inclusion"


class SubSuffix(EditEvent):
    "Particular case of inclusion"


class SuperPrefix(EditEvent):
    "Particular kind of split"


class SuperSuffix(EditEvent):
    "Particular case of split"


class OverlapPrefixSuffix(EditEvent):
    "Prefix/Suffix means edit = Split"


class OverlapSuffixPrefix(EditEvent):
    "Prefix/suffix means edit = split"


class SubString(EditEvent):
    "Affix succession means that edit = merge"


class SuperString(EditEvent):
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
                f"Could not compute edition between the paths {self.alignment_name} as they don't describe the same sequences.\n\n\t{''.join(path_a.values())}\n\t{''.join(path_b.values())}")
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

        print(f"\n\n\t\t>>>{self.alignment_name}<<<\n\n")

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
                        String(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            False
                        )
                    )
                    offset_reference += reference_end - reference_start
                    offset_query += query_end - query_start

                elif reference_start > query_start and reference_end < query_end:
                    # Reference is a subpart of other node
                    edit_list.append(
                        SuperString(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            reference_start - query_start,
                            reference_end - reference_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            True
                        )
                    )
                    offset_reference += reference_end - reference_start
                    offset_query += reference_end - reference_start

                elif reference_start < query_start and reference_end > query_end:
                    # Other node is included in reference
                    edit_list.append(
                        SubString(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            False
                        )
                    )
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start == query_start and reference_end > query_end:
                    edit_list.append(
                        SubPrefix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            False
                        )
                    )
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start == query_start and reference_end < query_end:
                    edit_list.append(
                        SuperPrefix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            True
                        )
                    )
                    offset_reference += reference_end - reference_start
                    offset_query += reference_end - reference_start

                elif reference_start < query_start and reference_end == query_end:
                    edit_list.append(
                        SubSuffix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            False
                        )
                    )
                    offset_reference += query_end - query_start
                    offset_query += query_end - query_start

                elif reference_start > query_start and reference_end == query_end:
                    edit_list.append(
                        SuperSuffix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            True
                        )
                    )
                    offset_reference += reference_end - reference_start
                    offset_query += reference_end - reference_start

                elif reference_start < query_start and reference_end < query_end:
                    edit_list.append(
                        OverlapPrefixSuffix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            0,
                            reference_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            True
                        )
                    )
                    offset_reference += reference_end - query_start
                    offset_query += reference_end - query_start

                elif reference_start > query_start and reference_end > query_end:
                    edit_list.append(
                        OverlapSuffixPrefix(
                            self.reference_path.keys()[idx_reference],
                            self.query_path.keys()[idx_query],
                            seq_reference,
                            seq_query,
                            offset_query,
                            reference_start-offset_query,
                            query_end-query_start,
                            offset_reference,
                            reference_start,
                            reference_end,
                            True
                        )
                    )
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

                print(
                    f"IdxRef={idx_reference}, Oref={offset_reference} || IdxQry={idx_query}, Oqry={offset_query}")

            self.edition[self.reference_path.keys()[current_pos]] = edit_list
            print(edit_list)
            edit_list: list = list()


def evaluate_consensus(*paths_editions: PathEdit) -> dict:
    """Given a series of PathEdit, evaluates if they are equal.

    Args:
        paths (PathEdit): any number of PathEdit objects

    Returns:
        bool: if paths edits are non ambiguous
    """
    for path_A, path_B in list(combinations(paths_editions, 2)):
        shared_nodes = set(path_A.edition.keys()).intersection(
            set(path_B.edition.keys()))
        if not (path_A.can_be_aligned and path_B.can_be_aligned and all([path_A.edition[x] == path_B.edition[x] for x in shared_nodes])):
            raise ValueError("No consensus could be made between paths.")
    all_edits: dict = dict()
    for pe in paths_editions:
        all_edits.update(pe.edition)
    return all_edits

################################################ EDIT ######################################


def edit_graph(
        graph_to_edit: Graph,
        edition_to_perform: dict[str, list[EditEvent]]
) -> Graph:
    """Given a graph and the consensus of a multipath edit, performs the required edition

    Args:
        graph_to_edit (Graph): a graph to be edited
        edition_to_perform (dict): a series of events to compute on the graph, per node of the reference

    Returns:
        Graph: edited graph object
    """
    # Edited graph will handle the splits of nodes
    edited_graph: Graph = Graph()
    # Output graph will handle merges and final output
    output_graph: Graph = Graph()

    mapping: dict = {}
    # We iterate over nodes. If we have to split, we split only in two, with the name being 'temp',<original_node_name>
    for nedit, qedit in edition_to_perform.items():
        print(f"{nedit} |-> {qedit}")

        # Deepcopy to block edits from propagating
        edited_graph.segments += deepcopy(
            [graph_to_edit.get_segment(segname.target_node) for segname in qedit])
        edited_graph.lines += deepcopy(
            list(chain.from_iterable([graph_to_edit.get_edges(segname.target_node) for segname in qedit])))

        for edition_event in qedit:
            if type(edition_event) in [SuperPrefix, OverlapPrefixSuffix, SuperString]:
                # Left border of node
                edited_graph.split_segments(
                    edition_event.target_node,
                    nodes_names := ['-1', edition_event.target_node],
                    splits := [
                        (
                            0,
                            edition_event.length
                        ),
                        (
                            edition_event.length,
                            edition_event.length_alt
                        )
                    ]
                )
                edition_event.target_node = '-1'
                print(vars(edition_event))
                print(f"Splitting at {splits} for resp. {nodes_names}")
            # Récupérer le segment splitté et le coller à nouveau si on le rerencontre en merge ?
            # Faire en sorte que le split renvoie simplement du graphe ... pb des aretes et chemins, du renommage peut-être plus efficace

        targets: list[str] = [x.target_node for x in qedit]
        output_graph.segments += deepcopy([edited_graph.get_segment(segname)
                                          for segname in targets])
        output_graph.lines += deepcopy(
            list(chain.from_iterable([edited_graph.get_edges(segname) for segname in targets])))
        if len(targets) > 1:
            # We merge only if we have multiple nodes to consider
            print(
                f"Nodes in output graph : {[n.datas['name'] for n in output_graph.segments]}")
            print(f"Merging {targets} in r{nedit}")
            output_graph.merge_segments(*targets, merge_name=f"r{nedit}")
            mapping[f"r{nedit}"] = nedit
        else:
            output_graph.rename_node(targets[0], f"r{targets[0]}")
            mapping[f"r{targets[0]}"] = nedit
    # TODO fix fonction dans la lib qui ne cherche que par le nombre et non par le nom complet qui empêche donc la fonction de réussir sa tâche
    for old_name, new_name in mapping.items():
        output_graph.rename_node(old_name, new_name)

    return output_graph


def perform_edition(files: list, gfa_versions: list, output_folder: str, do_edition: bool = True) -> None:
    """Main call

    Args:
        args (Namespace): arguments for the call
    """

    # Handling input errors
    if len(files) < 2:
        raise ValueError("Please specify at least two GFA files as input.")
    if len(files) != len(gfa_versions):
        raise ValueError(
            "Please match the number of args between files and gfa types.")

    # Iterating over couples of files
    for (f1, v1), (f2, v2) in combinations(zip(files, gfa_versions), 2):
        all_dipaths: list[PathEdit] = list()

        graph1: Graph = Graph(f1, v1, True)
        graph2: Graph = Graph(f2, v2, True)

        paths = set(path.datas['name'] for path in graph1.get_path_list()).intersection(
            set(path.datas['name'] for path in graph2.get_path_list()))

        # Iterating over paths
        for dpath in paths:
            path_in_g1 = graph1.get_path(dpath)
            path_in_g2 = graph2.get_path(dpath)

            pog1 = IndexedOrderedDict()
            pog2 = IndexedOrderedDict()

            for node, _ in path_in_g1.datas['path']:
                pog1[node] = graph1.get_segment(node).datas['seq']

            for node, _ in path_in_g2.datas['path']:
                pog2[node] = graph2.get_segment(node).datas['seq']

            all_dipaths.append(PathEdit(dpath, pog1, pog2))

        with open("out.log", "w", encoding='utf-8') as output:
            for dipath in all_dipaths:
                for key, value in dipath.edition.items():
                    if not isinstance(value[0], String):
                        output.write(
                            f"{dipath.alignment_name}\t{key}\t{value}\n")

        # We evaluate we can perform the edition befor actually doing it
        edit_moves: dict = evaluate_consensus(*all_dipaths)

        if do_edition:

            edited: Graph = edit_graph(graph2, edit_moves)

            edited.save_graph(path.join(
                output_folder, f"edited_graph_{Path(f1).stem}_{Path(f2).stem}.gfa"))

        del graph1, graph2


if __name__ == "__main__":

    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "file", type=str, help="Path(s) to two or more gfa-like file(s).", nargs='+')
    parser.add_argument(
        "-o", "--output_folder", required=True, type=str, help="Path to a folder for results.")
    parser.add_argument(
        "-g", "--gfa_version", help="Tells the GFA input style", required=True, choices=['rGFA', 'GFA1', 'GFA1.1', 'GFA1.2', 'GFA2'], nargs='+')
    parser.add_argument(
        "-p", "--perform_edition", help="Asks to perform edition on graph and outputs it.", action='store_true')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another.')
    args = parser.parse_args()

    perform_edition(args.file, args.gfa_version,
                    args.output_folder, args.perform_edition)
