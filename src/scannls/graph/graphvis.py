"""Plot Graphs.

@Author:      YangyangLi
@Time:        1/28/22 8:46 PM
"""
from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar

import networkx as nx

if TYPE_CHECKING:
    from . import NLGraph, Node


def default_visitors(graph: NLGraph, figure_name: str, support_reads: int) -> GraphVis:
    return GraphVis.from_visitors(
        graph,
        [MatplotlibVisualizeGraph(figure_name, support_reads)],
    )


class GraphVis:
    GRAPH_LINK_DATA: ClassVar = {
        "link": "edges",
        "source": "from",
        "target": "to",
    }

    def __init__(self, nlgraph: NLGraph):
        self.nlgraph = nlgraph
        self.visitors: list[GraphVisitor] = []

    @classmethod
    def from_visitors(cls, nlgraph: NLGraph, visitors: list[GraphVisitor]) -> GraphVis:
        """Create GraphVis from visitors."""
        graph_vis = cls(nlgraph)
        graph_vis.visitors = visitors
        return graph_vis

    @staticmethod
    def load(file_name: str | Path):
        if isinstance(file_name, str):
            file_name = Path(file_name)

        if not file_name.exists():
            msg = f"{file_name} not exists."
            raise FileNotFoundError(msg)

        with file_name.open() as f:
            data = json.load(f)
            return nx.node_link_graph(data, **GraphVis.GRAPH_LINK_DATA)

    def register(self, visitor: GraphVisitor) -> None:
        """Register visitor."""
        self.visitors.append(visitor)

    def visualize(self, *, visitors: list[GraphVisitor] | None = None) -> None:
        """Plot graph."""
        g = self.create_nxgraph()
        if visitors is None:
            visitors = []

        for visitor in self.visitors + visitors:
            visitor.visit(g)

    @staticmethod
    def get_label_from_node(node: Node) -> str:
        """Get label from node."""
        head_node = "H" if node.is_start_node() else "T"
        return f"{node.chrom}_{node.ref_start}_{node.ref_end}_{head_node}{node.strand}"

    def create_nxgraph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        labels = {}

        for start_node in self.nlgraph.get_start_nodes():
            if not start_node.successors:
                g.add_node(self.get_label_from_node(start_node))
            else:
                labels.update({start_node: self.get_label_from_node(start_node)})
                self._traverse_graph(start_node, [start_node], g, self.nlgraph, labels)  # type: ignore

        return g

    def _traverse_graph(
        self,
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
                    labels.update({successor: self.get_label_from_node(successor)})
                    nx_graph.add_edge(
                        self.get_label_from_node(start_node),
                        self.get_label_from_node(successor),
                        weight=edge.sr,
                    )
                    self._traverse_graph(
                        successor,
                        [*path, edge, successor],
                        nx_graph,
                        graph,
                        labels,
                    )
        else:
            # successor be [] or None
            self._traverse_graph(successors, [*path], nx_graph, graph, labels)


# https://networkx.org/documentation/latest/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
# https://networkx.org/documentation/latest/reference/drawing.html


class GraphVisitor(ABC):
    @abstractmethod
    def visit(self, graph: nx.DiGraph) -> None:
        """Visit graph."""


class MatplotlibVisualizeGraph(GraphVisitor):
    def __init__(self, figure_name: str, support_reads: int):
        self.figure_name = figure_name
        self.support_reads = support_reads

    def visit(self, graph: nx.DiGraph) -> None:
        """Visit graph."""
        visualize_graph_via_matplot(graph, self.figure_name, self.support_reads)


class PyvisVisualizeGraph(GraphVisitor):
    def __init__(self, figure_name: str):
        self.figure_name = figure_name

    def visit(self, graph: nx.DiGraph) -> None:
        """Visit graph."""
        visualize_graph_via_pyvis(graph, self.figure_name)


class GraphExporter(GraphVisitor):
    def __init__(self, file_name: str | Path, graph_link_data):
        if isinstance(file_name, str):
            self.file_name = Path(file_name)

        self.graph_link_data = graph_link_data

    def visit(self, graph: nx.DiGraph):
        data = nx.node_link_data(graph, **self.graph_link_data)
        with Path(f"{self.file_name}.json").open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=4)


def visualize_graph_via_matplot(
    graph: nx.DiGraph,
    figure_name: str,
    support_reads: int,
) -> None:
    from matplotlib import pyplot as plt  # type: ignore

    fig, ax = plt.subplots(figsize=(15, 15))

    pos = nx.spring_layout(graph, seed=42)
    edge_labels = nx.get_edge_attributes(graph, "weight")
    options = {
        "font_size": 10,
        "node_size": 1000,
        "node_color": ["red" if "H" in n else "white" for n in graph],
        "edgecolors": "black",
        "edge_color": [
            "red" if weight >= support_reads else "black"
            for weight in edge_labels.values()
        ],
        "linewidths": 2,
        "width": 3,
    }

    nx.draw_networkx(graph, pos=pos, arrows=True, **options)

    nx.draw_networkx_edge_labels(graph, pos, edge_labels)

    ax.set_title(f"Node number: {len(list(graph))}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"graph_{figure_name}.png")


def visualize_graph_via_pyvis(graph: nx.DiGraph, figure_name: str | Path) -> None:
    from pyvis.network import Network  # type: ignore

    nt = Network(height="750px", directed=True, width="100%")
    nt.from_nx(graph)
    nt.save_graph(f"graph_{figure_name}.html")
