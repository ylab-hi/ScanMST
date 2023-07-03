"""Find cliques in a graph.

@Filename:    cluster.py
@license:     MIT Licence
@Time:        1/19/22 7:59 PM
"""
from itertools import combinations
from typing import Any

import networkx as nx
from loguru import logger
from networkx import connected_components

from ..graph import NLPath
from ..graph import Node
from ..utils import timeit
from .mergeCondition import MergeCondition


class Ruler:
    """Calculate the similarity distance between two series.

    using longer one as the reference
    """

    def __init__(self, prune_threshold: int = 10) -> None:
        """Initialize Ruler.

        :param logger: logger
        """
        self.prune_threshold = prune_threshold

    def __repr__(self):
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
        connection = False
        merge_condition: MergeCondition = MergeCondition(self.prune_threshold)

        if (
            merge_condition.head2head(head_node_a, head_node_b)
            or merge_condition.tail2tail(tail_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_b, tail_node_a)
            or _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
                middle_nodes_a, middle_nodes_b
            )
        ):
            connection = True

        if not connection:
            for _middle_node_b in middle_nodes_b:
                if merge_condition.head2mid(head_node_a, _middle_node_b):
                    connection = True
                    break
                if merge_condition.tail2mid(tail_node_a, _middle_node_b):
                    connection = True
                    break

            if not connection:
                for _middle_node_a in middle_nodes_a:
                    if merge_condition.head2mid(head_node_a, _middle_node_a):
                        connection = True
                        break
                    if merge_condition.tail2mid(tail_node_b, _middle_node_a):
                        connection = True
                        break

        return 0.0 if connection else 1.0


def create_merge_key_for_node(node: Node):
    introns = node.introns
    introns_key = "-".join([f"{i}-{j}" for i, j in introns]) if introns else "None"
    chrom = node.chrom

    return f"{chrom}_{introns_key}"


def create_merge_key_for_nlpath(nodes_key: list[str]):
    return "_".join(nodes_key)


def create_sort_key_for_node(node: Node):
    return middle_node_signature(node)


def creat_sort_key_for_nlpath(nlpath: NLPath):
    return (
        len(nlpath),
        *[create_sort_key_for_node(node) for node in nlpath],
    )


def create_sort_key_by_merge_factor(nlpath: NLPath):
    return (
        len(nlpath),
        nlpath.merge_factor,
        *[create_sort_key_for_node(node) for node in nlpath],
    )


def sort_cluster(cluster, key=lambda x: creat_sort_key_for_nlpath(x), reverse=False):
    return sorted(
        cluster,
        key=key,
        reverse=reverse,
    )


def middle_node_signature(node: Node) -> str:
    """Middle node signature using chrom, exons and strand."""
    chrom = node.chrom
    exons = node.exons
    strand = node.strand
    exons_string = map(lambda x: f"{x[0]}-{x[1]}", exons)  # type: ignore

    return f"{chrom}:{';'.join(exons_string)};{strand}"


def _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
    node_list1: list[Node], node_list2: list[Node]
) -> bool:
    """Check if two node list of middle nodes have shared node or not.
    :param node_list1:  node_list1
    :param node_list2:  node_list2
    :return:  True if two node list have shared node, otherwise False
    """
    shared_middle_nodes = set(map(middle_node_signature, node_list1)) & set(
        map(middle_node_signature, node_list2)
    )

    return len(shared_middle_nodes) > 0


def merge_nlpath(path1: NLPath, path2: NLPath, start_index: int):
    for idx, (updated_node, current_node) in enumerate(
        zip(path1[start_index : start_index + len(path2)], path2)  # type: ignore
    ):
        # update exon coordinates
        updated_node.ref_start = min(  # type: ignore
            updated_node.exons[0][0], current_node.exons[0][0]  # type: ignore
        )
        updated_node.exons[0] = updated_node.ref_start, updated_node.exons[0][1]  # type: ignore

        updated_node.ref_end = max(  # type: ignore
            updated_node.exons[-1][1], current_node.exons[-1][1]  # type: ignore
        )
        updated_node.exons[-1] = updated_node.exons[-1][0], updated_node.ref_end  # type: ignore

        # update edge data
        node1_edge = path1.get_edge(nodes=updated_node, nodes_idx=start_index + idx)
        node2_edge = path2.get_edge(nodes=current_node, nodes_idx=idx)

        if node2_edge is not None and node1_edge is not None:
            node1_edge.sr += node2_edge.sr
            # WARN: Update break point in covering way <Yangyang Li>
            node1_edge.break_point1 = node2_edge.break_point1
            node1_edge.break_point2 = node2_edge.break_point2


def merge_same_len_node_list(
    path1: NLPath, path2: NLPath, start_index: int, threashold: int
) -> bool:
    """seires1 is equal than series2 and series1 merge series2.

    orignial s1: [ ] - [ ] - [ ] - [ ]
    s2:                [ ] - [ ] - [ ]
    """
    logger.debug(f"merge: nlpath1:{path1} nlpath2:{path2}")
    assert len(path1) == len(path2)

    merge_condition = MergeCondition(threashold)

    flag = True

    for idx, (node1, node2) in enumerate(
        zip(path1[start_index : start_index + len(path2)], path2)  # type: ignore
    ):
        node1_edge = path1.get_edge(nodes=node1, nodes_idx=start_index + idx)
        node2_edge = path2.get_edge(nodes=node2, nodes_idx=idx)

        same_edge = True
        if node2_edge is not None and node1_edge is not None:
            same_edge = node1_edge.is_merged(
                node2_edge, compared_break_point=False, break_point_threshold=threashold
            )

        if node1.self_identity.is_head() and node2.self_identity.is_head():
            if not (same_edge and merge_condition.head2head(node1, node2)):
                flag = False
                break

        elif node1.self_identity.is_mid() and node2.self_identity.is_head():
            if not (same_edge and merge_condition.mid2head(node1, node2)):
                flag = False
                break

        elif node1.self_identity.is_mid() and node2.self_identity.is_mid():
            if not (same_edge and merge_condition.mid2mid(node1, node2)):
                flag = False
                break

        elif node1.self_identity.is_mid() and node2.self_identity.is_tail():
            if not merge_condition.mid2tail(node1, node2):
                flag = False
                break

        elif node1.self_identity.is_tail() and node2.self_identity.is_tail():
            if not merge_condition.tail2tail(node1, node2):
                flag = False
                break

        else:
            raise ValueError(
                f"invalid node identity {node1.self_identity=} {node2.self_identity=}"
            )

    return flag


class ClusterFinder:
    """Find cliques in a graph based on series level.

     which will help to construct splice graph base on nodes level in the future.

    :param intact_series_list: list of intact series for a bam file of one sample
    :param logger: logger
    :param threshold: threshold to determine whether two series are connected

    .. note::
        :function: `networkx.algorithms.components.connected.connected_components` is used to find clusters.

    :Example:

    >>> from loguru import  logger
    >>> clique_finder = CliqueFinder([], logger)
    >>> clique_finder.find_clique()
    """

    def __init__(
        self,
        intact_nlpaths: list[NLPath],
        threshold: float = 0.2,
    ) -> None:
        """Initialize CliqueFinder."""
        self.ruler = Ruler()
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

    @timeit
    def find_cluster(self) -> Any:
        """Find clique in graph with help of :func:`networkx.algorithms.components.connected.connected_components`.

        :return:  every clique in graph as a iterator (List[Series])

        :Example:

        >>> from loguru import logger
        >>> clique_finder = CliqueFinder([], logger)
        >>> cliques = clique_finder.find_clique()
        >>> for clique in cliques:
        ...     for series_list in clique:
        ...         assert isinstance(series_list, Series)
        """
        self._create_graph_for_nlpath()

        for clique_index in connected_components(self._graph):
            yield (self.intact_nlpaths[i] for i in clique_index)

    def find_cluster_index(self):
        """Find clique in graph with help of :func:`networkx.algorithms.components.connected.connected_components`.

        :return:  every clique in graph as a iterator (List[int])
        """
        self._create_graph_for_nlpath()
        yield from connected_components(self._graph)

    def creat_merge_indexs(self, cluster) -> dict[int, list[str]]:
        result = {}

        for series in cluster:
            nodes_key = []
            for node in series:
                nodes_key.append(create_merge_key_for_node(node))

            result[series.id] = nodes_key

        return result

    @staticmethod
    def check_if_two_nlpath_merge(
        path1: NLPath, path2: NLPath, merge_keys: dict[int, list[str]]
    ):
        series_2_nodes_key = "".join(merge_keys[path2.id])

        for start_index in range(0, len(path1) - len(path2) + 1):
            series_1_nodes_key = "".join(
                merge_keys[path1.id][start_index : start_index + len(path2)]
            )

            if series_1_nodes_key == series_2_nodes_key:
                if merge_same_len_node_list(path1, path2, start_index, 1):  # type: ignore
                    merge_nlpath(path1, path2, start_index)  # type: ignore
                    path1.merge_factor += 1
                    return True

        return False

    @staticmethod
    def check_merge(path1: NLPath, path2: NLPath, merge_keys: dict[int, list[str]]):
        """Check if series1 can merge series2."""
        if len(merge_keys[path1.id]) == len(merge_keys[path2.id]):
            # reduce duplication
            return False
        elif len(merge_keys[path1.id]) > len(merge_keys[path2.id]):
            return ClusterFinder.check_if_two_nlpath_merge(path1, path2, merge_keys)
        else:
            raise ValueError("series1 is shorter than series2")

    def merge_cluster(self):
        for cluster_index in self.find_cluster_index():
            nlpaths = []
            for i in cluster_index:
                current_nlpath = self.intact_nlpaths[i]
                current_nlpath.id = i
                nlpaths.append(current_nlpath)

            sorted_nlpaths = sort_cluster(nlpaths)
            logger.debug(f"sorted_series:{len(sorted_nlpaths)} {sorted_nlpaths}")

            merge_keys = self.creat_merge_indexs(sorted_nlpaths)
            new_cluster = []
            ClusterFinder._merge_cluster(sorted_nlpaths, new_cluster, merge_keys)
            yield sort_cluster(
                new_cluster,
                key=lambda x: create_sort_key_by_merge_factor(x),  # type: ignore
                reverse=True,
            )

    @staticmethod
    def _merge_cluster(nlpaths, result, merge_keys):
        while nlpaths:
            slected_nlpath = nlpaths.pop()

            for current_nlpath in nlpaths:
                if ClusterFinder.check_merge(
                    slected_nlpath, current_nlpath, merge_keys
                ):
                    nlpaths.remove(current_nlpath)

            result.append(slected_nlpath)
        return result
