"""Plot Graphs.

@Author:      YangyangLi
@Time:        1/28/22 8:46 PM
"""
import json
from pathlib import Path

import networkx as nx

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


def get_label_from_node(node: Node) -> str:
    """Get label from node."""
    head_node = "S" if node.is_start_node() else "N"
    return f"{node.chrom}_{node.ref_start}_{node.ref_end}_{node.sr}_{node.strand}_{head_node}"


def output_graph(
    graph,
    file_name: str | Path,
):
    if isinstance(file_name, str):
        file_name = Path(file_name)

    g = nx.DiGraph()
    labels = {}

    for start_node in graph.get_start_nodes():
        if not start_node.successors:
            g.add_node(get_label_from_node(start_node))
        else:
            labels.update({start_node: get_label_from_node(start_node)})
            plot_graph_helper(start_node, [], g, graph)  # type: ignore

    data = nx.node_link_data(g, **GRAPH_LINK_DATA)
    with open(f"{file_name}.json", "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=4)


def plot_graph_helper(
    start_node: Node,
    path,
    nx_graph: nx.Graph,
    graph,
    labels: dict[Node, str],
) -> None:
    """Plot graph helper."""
    if not start_node or start_node in path:
        return

    if successors := start_node.successors:
        for successor in successors:
            for edge in graph.get_possible_edges(path, start_node, successor):
                labels.update({successor: get_label_from_node(successor)})
                nx_graph.add_edge(
                    get_label_from_node(start_node),
                    get_label_from_node(successor),
                )
                plot_graph_helper(
                    successor,
                    [*path, start_node, edge],
                    nx_graph,
                    graph,
                    labels,
                )
    else:
        # successor be [] or None
        plot_graph_helper(successors, [*path, start_node], nx_graph, graph, labels)


def plot_graph(graph, figure_name: str, *, is_matplotlib=True) -> None:
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
