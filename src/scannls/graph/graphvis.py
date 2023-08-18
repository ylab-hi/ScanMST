"""Plot Graphs.

@Author:      YangyangLi
@Time:        1/28/22 8:46 PM
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from . import Node

GRAPH_LINK_DATA = {"link": "edges", "source": "from", "target": "to"}


def read_graph(file_name: str | Path):
    if isinstance(file_name, str):
        file_name = Path(file_name)

    if not file_name.exists():
        msg = f"{file_name} not exists."
        raise FileNotFoundError(msg)

    with file_name.open() as f:
        data = json.load(f)
        return nx.node_link_graph(data, **GRAPH_LINK_DATA)


def plot_graph(graph, figure_name: str, *, is_matplotlib=True) -> None:
    """Plot graph."""
    g = create_nxgraph(graph)

    if is_matplotlib:
        visualize_graph_via_matplot(g, figure_name)
    else:
        visualize_graph_via_pyvis(g, figure_name)


def output_graph(
    graph,
    file_name: str | Path,
):
    if isinstance(file_name, str):
        file_name = Path(file_name)

    g = create_nxgraph(graph)

    data = nx.node_link_data(g, **GRAPH_LINK_DATA)
    with Path(f"{file_name}.json").open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=4)


def get_label_from_node(node: Node) -> str:
    """Get label from node."""
    head_node = "H" if node.is_start_node() else "T"
    return f"{node.chrom}_{node.ref_start}_{node.ref_end}_{head_node}{node.strand}"


def create_nxgraph(graph) -> nx.DiGraph:
    g = nx.DiGraph()
    labels = {}

    for start_node in graph.get_start_nodes():
        if not start_node.successors:
            g.add_node(get_label_from_node(start_node))
        else:
            labels.update({start_node: get_label_from_node(start_node)})
            traverse_graph(start_node, [start_node], g, graph, labels)  # type: ignore

    return g


def traverse_graph(
    start_node: Node,
    path,
    nx_graph: nx.Graph,
    graph,
    labels: dict[Node, str],
) -> None:
    """Plot graph helper."""
    if not start_node:
        return

    if successors := start_node.successors:
        for successor in successors:
            for edge in graph.get_possible_edges(path, start_node, successor, 1):
                labels.update({successor: get_label_from_node(successor)})
                nx_graph.add_edge(
                    get_label_from_node(start_node),
                    get_label_from_node(successor),
                    weight=edge.sr,
                )
                traverse_graph(
                    successor,
                    [*path, edge, successor],
                    nx_graph,
                    graph,
                    labels,
                )
    else:
        # successor be [] or None
        traverse_graph(successors, [*path], nx_graph, graph, labels)


# https://networkx.org/documentation/latest/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
# https://networkx.org/documentation/latest/reference/drawing.html


def visualize_graph_via_matplot(graph, figure_name: str) -> None:
    from matplotlib import pyplot as plt  # type: ignore

    fig, ax = plt.subplots(figsize=(15, 15))

    options = {
        "font_size": 10,
        "node_size": 1000,
        "node_color": ["red" if "H" in n else "white" for n in graph],
        "edgecolors": "black",
        "linewidths": 2,
        "width": 3,
    }

    pos = nx.spring_layout(graph, seed=42)

    nx.draw_networkx(graph, pos=pos, arrows=True, **options)

    edge_labels = nx.get_edge_attributes(graph, "weight")
    nx.draw_networkx_edge_labels(graph, pos, edge_labels)

    ax.set_title(f"Node number: {len(list(graph))}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"graph_{figure_name}.png")


def visualize_graph_via_pyvis(graph, figure_name: str) -> None:
    from pyvis.network import Network  # type: ignore

    nt = Network(height="750px", directed=True, width="100%")
    nt.from_nx(graph)
    nt.save_graph(f"graph_{figure_name}.html")
