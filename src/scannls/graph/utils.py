from collections import defaultdict
from pathlib import Path

import HTSeq
import intervaltree

from .nlgraph import NLGraph, Node


def get_neighbors(node: Node) -> list[Node]:
    return node.predecessors + node.successors


def is_weakly_connected(graph: NLGraph) -> bool:
    """Check if the graph is weakly connected."""
    visited = set()
    start_node = list(graph.get_start_nodes())

    if len(start_node) == 0:
        return False

    queue = [start_node[0]]

    while len(queue) > 0:
        current_node = queue.pop()
        visited.add(current_node)

        for neighbor in get_neighbors(current_node):
            if neighbor not in visited:
                queue.append(neighbor)

    return len(visited) == len(graph.nodes)


def build_interval_tree(
    annotation_source: Path,
) -> dict[str, intervaltree.IntervalTree]:
    """Build interval tree for all genes in terms of chromosome.

    Build interval tree for all genes in terms of chromosome from a gff file.
    """

    interval_trees = defaultdict(intervaltree.IntervalTree)

    for feature in HTSeq.GFF_Reader(annotation_source):
        if feature.type == "gene":
            interval_trees[feature.iv.chrom].add(
                intervaltree.Interval(
                    feature.iv.start,
                    feature.iv.end,
                    feature.name,
                ),
            )

    return interval_trees


def annotate_node(node: Node, interval_trees: dict[str, intervaltree.IntervalTree]):
    """Annotate a node with a string."""
    interval_tree = interval_trees[node.chrom]
    if interval_tree.is_empty():
        return

    genes = interval_tree.overlap(node.ref_start, node.ref_end)

    if len(genes) == 0:
        return

    node.gene_names = [gene.data for gene in genes]


def annotate_graph(graph: NLGraph, node: Node, annotation_source: Path):
    """Annotate a node with a string.

    .. note::
         1. build interval tree for all genes in terms of chromosome
         2. for each node, search the specific interval tree for genes that overlap based on chromosome
         3. annotate the node with the gene name
         4. annotate the graph with the high frequent gene name
    """
    from collections import Counter

    interval_trees = build_interval_tree(annotation_source)

    for node in graph:
        annotate_node(node, interval_trees)

    gene_names = Counter(gene for node in graph for gene in node.gene_names)

    return gene_names.most_common(1)[0][0]
