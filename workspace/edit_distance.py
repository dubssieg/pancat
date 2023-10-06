"Modelizes edition between paths and computes if edition is acceptable inside a graph"
from itertools import combinations, accumulate, chain
from argparse import ArgumentParser, SUPPRESS
from collections import Counter
from os import path
from pathlib import Path
from copy import deepcopy
from indexed import IndexedOrderedDict
from gfagraphs import Graph
from tharospytools import revcomp
from workspace.grapher import display_graph


class EditEvent():
    """Modelizes a elementary operation between a couple of nodes."""

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
        return f"{self.event}({round(self.affine_score(),ndigits=1)})"  # (revcomp={self.revcomp}, positions={self.offsets_target}, refnode={self.reference_node}, qnode={self.target_node})

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, __o: object) -> bool:
        return isinstance(__o, type(self)) and self.positions_target == __o.positions_target and self.revcomp == __o.revcomp

    def affine_score(self) -> int:
        "Basic implementation of score"
        return 0


class String(EditEvent):
    "Situation of two equal nodes"

    def affine_score(
        self,
            event_score: float = 0,
            shared_score: float = 0.1,
            diverging_score: float = 0
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + self.length*shared_score + 0*diverging_score


class SubPrefix(EditEvent):
    "Particular kind of inclusion"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class SubSuffix(EditEvent):
    "Particular case of inclusion"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class SuperPrefix(EditEvent):
    "Particular kind of split"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class SuperSuffix(EditEvent):
    "Particular case of split"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class OverlapPrefixSuffix(EditEvent):
    "Prefix/Suffix means edit = Split"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.2
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class OverlapSuffixPrefix(EditEvent):
    "Prefix/suffix means edit = split"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.2
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class SubString(EditEvent):
    "Affix succession means that edit = merge"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


class SuperString(EditEvent):
    "Cut on both ends"

    def affine_score(
        self,
            event_score: float = -1,
            shared_score: float = 0.1,
            diverging_score: float = -0.1
    ) -> float:
        """Computes the score for edit event
        Idea is that the number of shared bases account for n*shared_score,
        the number of diverging bases account for n*diverging_score,
        and the event has a base score named event_score.

        Args:
            event_score (float, optional): Base score for event. Defaults to -1.
            shared_score (float, optional): Additional score for each shared nucleotide. Defaults to 0.1.
            diverging_score (float, optional): Additional score for each not shared nucleotide. Defaults to -0.1.

        Returns:
            float: score for the event
        """
        return event_score + min(self.length, self.length_alt)*shared_score + (max(self.length, self.length_alt)-min(self.length, self.length_alt))*diverging_score


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
        # Validation of coordinates before starting computing edition
        if not len(''.join(path_a.values())) == len(''.join(path_b.values())):
            raise ValueError(
                f"Path does not have the same length, probably not defining the same sequence ({len(''.join(path_a.values()))} vs. {len(''.join(path_b.values()))})."
            )
        self.compute_edition()
        self.counts = self.get_counts()

    def get_counts(self) -> Counter:
        """Counts each event type in pathedit

        Returns:
            Counter: the counts of each type
        """
        counter: Counter = Counter()
        for edit in self.edition.values():
            counter += Counter([type(e).__name__ for e in edit])
        return counter

    def compute_edition(self):
        """Analyses each relation between nodes of the paths, from left to right,
        and edits the dictionnary of editions containing a mapping ref_node:editions

        Raises:
            ValueError: if encountered event is invalid (bad index comparison, should not happend, but here in case)
        """
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
            # While we're on current node according to offset, we append edition to node
            while current_pos == idx_reference:

                reference_start: int = reference_endpos[idx_reference]
                reference_end: int = reference_endpos[idx_reference+1]

                query_start: int = query_endpos[idx_query]
                query_end: int = query_endpos[idx_query+1]

                if reference_end == query_end and reference_start == query_start:
                    # No edition to be made, it's an equivalence
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

            self.edition[self.reference_path.keys()[current_pos]] = edit_list
            edit_list: list = list()


def edit_by_paths(
        graph_to_edit: Graph,
        edit_on_paths: list[PathEdit],
        out_path: str,
        plot_graph: bool = False,
) -> tuple[Graph, int, int]:
    """Edition loop over graph.
    Path-aware edition. Limits collision between successive edits on multiple paths.

    Args:
        graph_to_edit (Graph): the graph to be edited
        edit_on_paths (list[PathEdit]): for each path of the graph, the events to apply per node
        plot_graph (bool, optional): if input and output graph should be plotted. Defaults to False.

    Returns:
        tuple[Graph, int, int]: the edited graph, the number of elementary merges done, the number of elementary splits done
    """
    # Init empty return graph
    full_graph: Graph = Graph()

    graph_path: dict[str, Graph] = {
        edit.alignment_name: Graph() for edit in edit_on_paths}

    number_of_splits: int = 0
    number_of_merges: int = 0
    score: float = 0

    # We iterate on each path
    for pathedit in edit_on_paths:
        # Init values for alignment edition
        alignment: str = pathedit.alignment_name

        # We copy the path we want to follow
        graph_path[alignment].walks = [
            deepcopy(w) for w in graph_to_edit.walks if w.datas['name'] == pathedit.alignment_name]
        graph_path[alignment].paths = [
            deepcopy(p) for p in graph_to_edit.paths if p.datas['name'] == pathedit.alignment_name]

        mapping: dict = {}
        negative_edits: int = 0
        # We iterate over nodes. If we have to split, we split only in two, with the name being 'temp',<original_node_name>
        for nedit, qedit in pathedit.edition.items():

            # Deepcopy to block edits from propagating
            graph_path[alignment].segments += deepcopy(
                [graph_to_edit.get_segment(segname.target_node) for segname in qedit if segname.target_node not in [s.datas['name'] for s in graph_path[alignment].segments]])

            for edition_event in qedit:
                score += edition_event.affine_score()
                if type(edition_event) in [SuperPrefix, OverlapPrefixSuffix, SuperString]:
                    negative_edits -= 1
                    # We want alternate positions for not erasing nodes in specific cases
                    # I do deserve hell for this
                    nodes_names: list = [f'{negative_edits}', edition_event.target_node][::-(
                        type(edition_event) not in [OverlapPrefixSuffix, SuperString]) or 1]
                    # Left border of node
                    future_split = graph_to_edit.get_segment(
                        edition_event.target_node)
                    # Splitting in two (elementary split)
                    graph_path[alignment].split_segments(
                        edition_event.target_node,
                        nodes_names,
                        [
                            (
                                0,
                                min(edition_event.length,
                                    edition_event.length_alt)
                            ),
                            (
                                min(edition_event.length,
                                    edition_event.length_alt),
                                future_split.datas['length']
                            )
                        ]
                    )
                    number_of_splits += 1
                    # Setting temporary name as target
                    edition_event.target_node = f'{negative_edits}'

        # Estimating targets
        targets_of_target: dict = {}
        cumulative_len: int = 0
        for nedit, qedit in pathedit.edition.items():
            number_of_edits = len(qedit)
            names_in_path: list = [
                x for p in graph_path[alignment].get_path_list() for (x, _) in p.datas['path']]
            targets_of_target[nedit] = names_in_path[cumulative_len:cumulative_len+number_of_edits]
            cumulative_len += number_of_edits

        for nedit, targets in targets_of_target.items():
            if len(targets) > 1:
                # We merge only if we have multiple nodes to consider
                graph_path[alignment].merge_segments(
                    *targets, merge_name=f"r{nedit}")
                mapping[f"r{nedit}"] = nedit
                # Adding number of required elementary merges for this merge
                number_of_merges += len(targets)-1
            else:
                # Otherwise we simply rename the node
                # Does not count as any elementary operation
                if not '-' in targets[0]:
                    graph_path[alignment].rename_node(
                        targets[0], f"r{targets[0]}")
                    mapping[f"r{targets[0]}"] = nedit

        for old_name, new_name in mapping.items():
            graph_path[alignment].rename_node(old_name, new_name)

    # Merging time in full graph!
    # We first merge walks (addition of walks of each graph)
    full_graph.walks = list(
        chain(*[alg.walks for alg in graph_path.values()]))
    full_graph.paths = list(
        chain(*[alg.paths for alg in graph_path.values()]))

    # Then we select segments, keeping memory of which we already copied and only by selecting them in path (edited in semigraph)
    copied_segments: list = list()
    for output_path in list(chain(full_graph.get_path_list())):
        following: str = output_path.datas['name']
        for seg, _ in output_path.datas['path']:
            if seg not in copied_segments:
                full_graph.segments.append(
                    deepcopy(graph_path[following].get_segment(seg)))
                copied_segments.append(seg)

    # We plot the output graph if required
    if plot_graph:
        total_counts: Counter = Counter()
        for dipath in edit_on_paths:
            total_counts += dipath.counts

        display_graph(
            full_graph.compute_networkx(),
            "output",
            full_graph.colors,
            IndexedOrderedDict({
                'Score': round(score, ndigits=1),
                'Merges': number_of_merges,
                'Splits': number_of_splits,
                'Total edits': number_of_splits+number_of_merges,
                **total_counts
            }),
            output_path=out_path
        )

    return full_graph, number_of_merges, number_of_splits


def perform_edition(
        files: list,
        gfa_versions: list,
        output_folder: str,
        do_edition: bool = True,
        do_plot: bool = False
) -> None:
    """Main loop for edition ; iterates over couples of files
    and performs edition on each couple if possible

    Args:
        files (list): a list of gfa files (in the same order as versions)
        gfa_versions (list): a list of gfa versions (in the same order as files)
        output_folder (str): a foler where to store results
        do_edition (bool, optional): if editions should be made. Defaults to True.

    Raises:
        ValueError: files and gfa versions argument error
    """
    # In this function, we do calculate the distance between G1 and G2, by trying to modify G2 into G1.
    # Note the two graphs can be freely swapped, we just need to invert scores

    # Handling input errors
    if len(files) < 2:
        raise ValueError("Please specify at least two GFA files as input.")
    if len(files) != len(gfa_versions):
        raise ValueError(
            "Please match the number of args between files and gfa types.")

    # Iterating over couples of files
    for i, ((file_1, version_1), (file_2, version_2)) in enumerate(combinations(zip(files, gfa_versions), 2)):
        all_dipaths: list[PathEdit] = list()

        graph1: Graph = Graph(file_1, version_1, True)
        graph2: Graph = Graph(file_2, version_2, True)

        # We perform edition on shared paths, hoping the best for non-common paths \o/
        # (Best practice is to validate before if all paths are shared)
        paths = set(path.datas['name'] for path in graph1.get_path_list()).intersection(
            set(path.datas['name'] for path in graph2.get_path_list()))

        # Iterating over paths
        for dpath in paths:
            print(f"Working on {dpath}...")
            path_in_g1 = graph1.get_path(dpath)
            path_in_g2 = graph2.get_path(dpath)

            pog1 = IndexedOrderedDict()
            pog2 = IndexedOrderedDict()

            for node, vect in path_in_g1.datas['path']:
                pog1[node] = graph1.get_segment(
                    node).datas['seq'] if vect == '+' else revcomp(graph1.get_segment(node).datas['seq'])

            for node, vect in path_in_g2.datas['path']:
                pog2[node] = graph2.get_segment(
                    node).datas['seq'] if vect == '+' else revcomp(graph2.get_segment(node).datas['seq'])

            # We append new path edit to stack
            all_dipaths.append(PathEdit(dpath, pog1, pog2))

        total_counts: Counter = Counter()
        for dipath in all_dipaths:
            total_counts += dipath.counts

        with open(path.join(output_folder, f"out_{Path(file_1).stem}_against_{Path(file_2).stem}_{i}.log"), "w", encoding='utf-8') as output:
            # Write edited stuff in log
            relations_counter: Counter = Counter()
            required_merges: int = 0
            for dipath in all_dipaths:
                for key, value in dipath.edition.items():
                    required_merges += max(0, len(value)-1)
                    relations_counter.update([type(val) for val in value])
                    # Variable value is a list of EditEvent
                    edits: str = ','.join(
                        [f"{val},{val.offsets_ref},{val.offsets_target}" for val in value])
                    output.write(
                        f"{dipath.alignment_name}\t{key}\t{edits}\n")
            output.write(f"Required merges : {required_merges} ")
            output.write(
                f"Required splits : {relations_counter[OverlapPrefixSuffix]+relations_counter[SuperPrefix]+relations_counter[SuperString]} ")
            output.write(str(relations_counter))

            # If we ask to compute edition, we do it
            if do_edition:

                edited, _, _ = edit_by_paths(
                    graph2, all_dipaths, output_folder, do_plot)

                edited.save_graph(path.join(
                    output_folder, f"edited_graph_{Path(file_1).stem}_{Path(file_2).stem}.gfa"))

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
    parser.add_argument(
        "-l", "--do_plot", help="Asks to plot the graph, with metrics.", action='store_true')
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Aligns graphs and tires to compute the minimal set of events that explains how we go from one graph to another.')
    args = parser.parse_args()

    perform_edition(
        args.file,
        args.gfa_version,
        args.output_folder,
        args.perform_edition,
        args.do_plot
    )
