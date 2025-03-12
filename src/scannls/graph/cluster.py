"""cluster NLpaths."""

from itertools import combinations

import networkx as nx
from loguru import logger
from networkx import connected_components

from scannls.base import MicroHomology, NovelInsertion

from .basic_graph import NLPath, Node
from .merge_condition import MergeCondition


class Ruler:
    """Calculate the similarity distance between two series.

    using longer one as the reference
    """

    def __init__(self, prune_threshold: int) -> None:
        """Initialize Ruler.

        :param logger: logger
        """
        self.prune_threshold = prune_threshold

    def __repr__(self) -> str:
        """Represent Ruler."""
        return f"{self.__class__.__name__}()"

    def __call__(self, nlpath_a: NLPath, nlpath_b: NLPath) -> float:
        """Call Ruler to calculate the distance between two series.

        :param series_a: series a
        :param series_b: series b
        :return: distance

        :Example:

        >>> ruler = Ruler()
        >>> ruler(series_a, series_b)
        0.5
        """
        if len(nlpath_a) < len(nlpath_b):
            nlpath_a, nlpath_b = nlpath_b, nlpath_a

        head_node_a = nlpath_a[0]
        head_node_b = nlpath_b[0]

        tail_node_a = nlpath_a[-1]
        tail_node_b = nlpath_b[-1]

        middle_nodes_a: list[Node] = nlpath_a[1:-1]  # type: ignore
        middle_nodes_b: list[Node] = nlpath_b[1:-1]  # type: ignore
        merge_condition: MergeCondition = MergeCondition(self.prune_threshold)

        if (
            merge_condition.head2head4distance(head_node_a, head_node_b)
            or merge_condition.tail2tail4distance(tail_node_a, tail_node_b)
            or merge_condition.head2head(head_node_a, head_node_b)
            or merge_condition.tail2tail(tail_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_b, tail_node_a)
            or _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
                middle_nodes_a,
                middle_nodes_b,
            )
        ):
            return 0.0

        for middle_node_b in middle_nodes_b:
            if merge_condition.head2mid(
                head_node_a,
                middle_node_b,
            ) or merge_condition.tail2mid(tail_node_a, middle_node_b):
                return 0.0

        for middle_node_a in middle_nodes_a:
            if merge_condition.head2mid(
                head_node_b,
                middle_node_a,
            ) or merge_condition.tail2mid(tail_node_b, middle_node_a):
                return 0.0

        return 1.0


def create_merge_key_for_node(node: Node):
    return f"{node.chrom}_{node.introns!s}"


def create_merge_key_for_nlpath(nodes_key: list[str]):
    return "_".join(nodes_key)


def create_sort_key_for_node(node: Node):
    return middle_node_signature(node)


def create_sort_key_for_nlpath(nlpath: NLPath):
    return (
        len(nlpath),
        *[create_sort_key_for_node(node) for node in nlpath],
    )

def obtain_edge_info_signature_for_nlpath(nlpath: NLPath):
    """Obtain edge info signature.
       for every hop
       blunt end: 2
       microhomology: 1
       microinsertion: 0
    """
    edge_info_signature = 0
    for event_id, current_node in enumerate(nlpath.nodes[:-1], 1):
        current_edge = nlpath.next_edge(current_node, event_id - 1)
        if current_edge.insertion_info:
            if isinstance(current_edge.insertion_info[1], NovelInsertion):
                edge_info_signature += 0
            elif isinstance(current_edge.insertion_info[1], MicroHomology):
                edge_info_signature += 1
        else:
            # blunt end
            edge_info_signature += 2
    return edge_info_signature


def create_sort_key_by_merge_factor(nlpath: NLPath):
    return (
        len(nlpath),
        nlpath.merge_factor,
        obtain_edge_info_signature_for_nlpath(nlpath),
        *[create_sort_key_for_node(node) for node in nlpath],
    )


def sort_cluster(cluster, key=create_sort_key_for_nlpath, *, reverse=False):
    return sorted(
        cluster,
        key=key,
        reverse=reverse,
    )


def middle_node_signature(node: Node) -> str:
    """Middle node signature using chrom, exons and strand."""
    exons_string = (f"{x[0]}-{x[1]}" for x in node.exons)  # type: ignore
    return f"{node.chrom}:{';'.join(exons_string)};{node.strand}"


def _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
    node_list1: list[Node],
    node_list2: list[Node],
) -> bool:
    """Check if two node list of middle nodes have shared node or not.
    :param node_list1:  node_list1
    :param node_list2:  node_list2
    :return:  True if two node list have shared node, otherwise False.
    """
    if len(node_list1) == 0 and len(node_list2) == 0:
        return False

    shared_middle_nodes = set(map(middle_node_signature, node_list1)) & set(
        map(middle_node_signature, node_list2),
    )

    return len(shared_middle_nodes) > 0


def merge_nlpath(path1: NLPath, path2: NLPath, start_index: int):
    logger.debug(f"merge {path1=}")
    logger.debug(f"merge {path2=}")
    for idx, (updated_node, current_node) in enumerate(
        zip(path1[start_index : start_index + len(path2)], path2),  # type: ignore
    ):
        updated_node.merge(current_node)

        # update edge data
        if (
            node1_edge := path1.next_edge(
                nodes=updated_node,
                nodes_idx=start_index + idx,
            )
        ) is not None and (node2_edge := path2.next_edge(nodes=current_node, nodes_idx=idx)) is not None:
            node1_edge.merge(node2_edge, current_node.strand, path2[idx + 1].strand, merge_insertion_info=False)


def merge_same_len_node_list(
    path1: NLPath,
    path2: NLPath,
    start_index: int,
    threshold: int,
) -> bool:
    """seires1 is equal than series2 and series1 merge series2.

    orignial s1: [ ] - [ ] - [ ] - [ ]
    s2:                [ ] - [ ] - [ ]
    """

    merge_condition = MergeCondition(threshold)

    for idx, (node1, node2) in enumerate(
        zip(path1[start_index : start_index + len(path2)], path2),  # type: ignore
    ):
        node1_edge = path1.next_edge(nodes=node1, nodes_idx=start_index + idx)
        node2_edge = path2.next_edge(nodes=node2, nodes_idx=idx)

        same_edge = True
        if node2_edge is not None and node1_edge is not None:
            same_edge = node1_edge.insertion_info == node2_edge.insertion_info and node1_edge.merged(
                node2_edge,
                compared_break_point=False,
            )

        if not same_edge:
            return False

        if node1.self_identity.is_head() and node2.self_identity.is_head():
            if not merge_condition.head2head(node1, node2):
                return False

        elif node1.self_identity.is_mid() and node2.self_identity.is_head():
            if not merge_condition.mid2head(node1, node2):
                return False

        elif node1.self_identity.is_mid() and node2.self_identity.is_mid():
            if not merge_condition.mid2mid(node1, node2):
                return False

        elif node1.self_identity.is_mid() and node2.self_identity.is_tail():
            if not merge_condition.mid2tail(node1, node2):
                return False

        elif node1.self_identity.is_tail() and node2.self_identity.is_tail():
            if not merge_condition.tail2tail(node1, node2):
                return False
        else:
            msg = f"invalid node identity node1.self_identity={node1.self_identity!r} node2.self_identity={node2.self_identity!r}"
            raise ValueError(
                msg,
            )

    return True


class ClusterFinder:
    """Find cliques in a graph based on series level.

     which will help to construct splice graph base on nodes level in the future.

    :param intact_series_list: list of intact series for a bam file of one sample
    :param logger: logger
    :param threshold: threshold to determine whether two series are connected

    .. note::
        :function: `networkx.algorithms.components.connected.connected_components` is used to find clusters.

    .. example::

    >>> from loguru import  logger
    >>> clique_finder = CliqueFinder([], logger)
    >>> clique_finder.find_clique()
    """

    def __init__(
        self,
        intact_nlpaths: list[NLPath],
        prune_threshold: int,
        threshold: float = 0.2,
    ) -> None:
        """Initialize CliqueFinder."""
        self.ruler = Ruler(prune_threshold)
        self.intact_nlpaths = intact_nlpaths
        self.intact_nlpaths_len = len(intact_nlpaths)
        self.threshold = threshold
        self.distance_dict: dict[tuple[int, int], float] = {}
        self._graph = nx.Graph()

    def _calculate_distance(self, x: int, y: int) -> float:
        """Calculate distance between two series. If distance has been calculated before.

        return True and distance value. Otherwise, calculate distance and return False and
        distance value.

        :param x: nlpath x
        :param y: nlpath y
        :return: is_calculated, distance value
        """
        return self.ruler(self.intact_nlpaths[x], self.intact_nlpaths[y])

    def _add_edge_between_two_nlpath(self, x: int, y: int) -> None:
        if self._calculate_distance(x, y) < self.threshold:
            self._graph.add_edge(x, y)
            if not self.intact_nlpaths[x].is_in_graph:
                self.intact_nlpaths[x].is_in_graph = True
            if not self.intact_nlpaths[y].is_in_graph:
                self.intact_nlpaths[y].is_in_graph = True

    def _create_graph_for_nlpath(self) -> None:
        """Create graph for all series in intact_series_list.

        add edge between two series in terms of the distance value

        :return: None
        """
        last_x = 0
        ind_x, ind_y = 0, 0

        for ind_x, ind_y in combinations(range(self.intact_nlpaths_len), 2):
            self._add_edge_between_two_nlpath(ind_x, ind_y)
            if last_x != ind_x:
                if not self.intact_nlpaths[last_x].is_in_graph:
                    self._graph.add_node(last_x)
                last_x = ind_x

        # solve last two node
        if not self.intact_nlpaths[ind_x].is_in_graph:
            self._graph.add_node(last_x)

        if not self.intact_nlpaths[ind_y].is_in_graph:
            self._graph.add_node(ind_y)

    def find_cluster_index(self):
        """Find cluster in graph with help of :func:`networkx.algorithms.components.connected.connected_components`.

        :return:  every clique in graph as a iterator (List[int])
        """
        self._create_graph_for_nlpath()
        yield from connected_components(self._graph)

    @staticmethod
    def creat_merge_indexs(cluster) -> dict[int, list[str]]:
        result = {}

        for series in cluster:
            nodes_key = []
            for node in series:
                nodes_key.append(create_merge_key_for_node(node))

            result[series.id] = nodes_key

        return result

    @staticmethod
    def check_if_two_nlpath_merge(
        path1: NLPath,
        path2: NLPath,
        merge_keys: dict[int, list[str]],
        threadhold: int,
    ) -> bool:
        nlpath_2_nodes_key = "".join(merge_keys[path2.id])

        for start_index in range(len(path1) - len(path2) + 1):
            nlpath_1_nodes_key = "".join(
                merge_keys[path1.id][start_index : start_index + len(path2)],
            )

            if nlpath_1_nodes_key == nlpath_2_nodes_key and merge_same_len_node_list(
                path1,
                path2,
                start_index,
                threadhold,
            ):
                merge_nlpath(path1, path2, start_index)  # type: ignore
                path1.polish_edges()
                path1.merge_factor += 1
                return True

        return False

    @staticmethod
    def check_merge(path1: NLPath, path2: NLPath, merge_keys: dict[int, list[str]], threshold: float):
        """Check if nlpath1 can merge nlpath2."""

        if len(merge_keys[path1.id]) >= len(merge_keys[path2.id]):
            return ClusterFinder.check_if_two_nlpath_merge(path1, path2, merge_keys, threshold)

        msg = "nlpath1 is shorter than nlpath2"
        raise ValueError(msg)

    def merge_cluster(self):
        for cluster_index in self.find_cluster_index():
            nlpaths = []
            for i in cluster_index:
                current_nlpath = self.intact_nlpaths[i]
                nlpaths.append(current_nlpath)

            sorted_nlpaths = sort_cluster(nlpaths)
            logger.debug(f"{len(sorted_nlpaths)=} nlpaths for merge: {sorted_nlpaths}")

            merge_keys = self.creat_merge_indexs(sorted_nlpaths)
            new_cluster = ClusterFinder._merge_cluster(sorted_nlpaths, merge_keys, self.ruler.prune_threshold)
            yield sort_cluster(
                new_cluster,
                key=create_sort_key_by_merge_factor,  # type: ignore
                reverse=True,
            )

    @staticmethod
    def _merge_cluster(nlpaths, merge_keys, threshold: int):
        result = []
        removed_nlpaths = set()  # Use a set to keep track of removed paths for efficiency

        while nlpaths:
            selected_nlpath = nlpaths.pop()
            # Iterate over a copy of the list to safely modify the original list
            for current_nlpath in list(nlpaths):
                if current_nlpath in removed_nlpaths:
                    continue  # Skip this one as it's already marked for removal

                # Check if the selected and current nlpaths should be merged
                if ClusterFinder.check_merge(selected_nlpath, current_nlpath, merge_keys, threshold):
                    removed_nlpaths.add(current_nlpath)

            # Remove all marked nlpaths from the main list after checking
            nlpaths = [nlpath for nlpath in nlpaths if nlpath not in removed_nlpaths]
            result.append(selected_nlpath)
        return result
