# !/usr/bin/env python
"""Plot Graphs.

@Filename:    plotGraph.py
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        1/28/22 8:46 PM
"""
from typing import Any
from typing import Dict
from typing import Union

import networkx as nx
from matplotlib import pyplot as plt

from .basicClass import Insertion
from .basicClass import Node

NodeType = Union[Node, Insertion]


def get_label_from_node(node: NodeType) -> str:
    """Get label from node."""
    query_name = node.query_name.split(",")
    return f"{node.chrom}_{node.ref_start}_{node.ref_end}:{len(query_name)}{node.is_start_node()}"


def plot_graph_helper(
    start_node: NodeType, path, nx_graph: nx.Graph, labels: Dict[NodeType, str]
) -> None:
    """Plot graph helper."""
    if not start_node or start_node in path:
        return
    else:
        if successors := start_node.successors:
            for successor in successors:
                labels.update({successor: get_label_from_node(successor)})
                nx_graph.add_edge(
                    get_label_from_node(start_node), get_label_from_node(successor)
                )
                plot_graph_helper(successor, path + [start_node], nx_graph, labels)
        else:
            # successor be [] or None
            plot_graph_helper(
                successors, path + [start_node], nx_graph, labels  # type: ignore
            )


def plot_graph(graph: Any, figure_name: str, is_matplotlib=True) -> None:
    """Plot graph."""
    g = nx.DiGraph()
    labels = {}
    for start_node in graph.get_start_nodes():
        labels.update({start_node: get_label_from_node(start_node)})
        plot_graph_helper(start_node, [], g, labels)  # type: ignore
    options = {
        "font_size": 10,
        "node_size": 100,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 2,
        "width": 2,
    }
    if is_matplotlib:
        fig, ax = plt.subplots(figsize=(20, 20))
        nx.draw_networkx(g, **options, ax=ax)
        ax.set_title(f"Node number: {len(list(graph))}")
        plt.axis("off")
        fig.tight_layout()
        plt.savefig(f"graph_{figure_name}.png")
    else:
        from pyvis.network import Network

        nt = Network(height="750px", directed=True, width="100%")
        nt.from_nx(g)
        nt.save_graph(f"graph_{figure_name}.html")
