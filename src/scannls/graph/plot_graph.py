"""Plot Graphs.

@Filename:    plotGraph.py
@Author:      YangyangLi
@Time:        1/28/22 8:46 PM
"""
from pathlib import Path
from typing import Any

import networkx as nx

from . import Node


def get_label_from_node(node: Node) -> str:
    """Get label from node."""
    return (
        f"{node.chrom}_{node.ref_start}_{node.ref_end}_{node.sr}_{node.is_start_node()}"
    )


def plot_graph_helper(
    start_node: Node, path, nx_graph: nx.Graph, labels: dict[Node, str]
) -> None:
    """Plot graph helper."""
    if not start_node or start_node in path:
        return

    if successors := start_node.successors:
        for successor in successors:
            labels.update({successor: get_label_from_node(successor)})
            nx_graph.add_edge(
                get_label_from_node(start_node), get_label_from_node(successor)
            )
            plot_graph_helper(successor, [*path, start_node], nx_graph, labels)
    else:
        # successor be [] or None
        plot_graph_helper(
            successors, [*path, start_node], nx_graph, labels  # type: ignore
        )


def export_graph(graph: Any, file_name: Path) -> None:
    """Export graph."""
    g = nx.DiGraph()
    labels = {}
    for start_node in graph.get_start_nodes():
        if not start_node.successors:
            g.add_node(get_label_from_node(start_node))
        else:
            labels.update({start_node: get_label_from_node(start_node)})
            plot_graph_helper(start_node, [], g, labels)  # type: ignore

    nx.write_adjlist(g, file_name)


def plot_graph(graph: Any, figure_name: str, *, is_matplotlib=True) -> None:
    """Plot graph."""
    g = nx.DiGraph()
    labels = {}
    for start_node in graph.get_start_nodes():
        if not start_node.successors:
            g.add_node(get_label_from_node(start_node))
        else:
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
        from matplotlib import pyplot as plt  # type: ignore

        fig, ax = plt.subplots(figsize=(20, 20))
        nx.draw_networkx(g, **options, ax=ax)
        ax.set_title(f"Node number: {len(list(graph))}")
        plt.axis("off")
        fig.tight_layout()
        plt.savefig(f"graph_{figure_name}.png")
    else:
        from pyvis.network import Network  # type: ignore

        nt = Network(height="750px", directed=True, width="100%")
        nt.from_nx(g)
        nt.save_graph(f"graph_{figure_name}.html")
