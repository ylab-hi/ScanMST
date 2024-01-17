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
from loguru import logger
from matplotlib import pyplot as plt  # type: ignore

if TYPE_CHECKING:
    from . import Edge, NLGraph, Node


def default_visitors(graph: NLGraph, figure_name: str, support_reads: int = 1) -> GraphVis:
    return GraphVis.from_visitors(
        graph,
        [
            GraphCytoscapeExporter(figure_name),
        ],
        support_reads,
    )


class GraphVis:
    GRAPH_LINK_DATA: ClassVar[dict[str, str]] = {
        "link": "edges",
        "source": "source",
        "target": "target",
    }

    def __init__(self, nlgraph: NLGraph, visitors: list[GraphVisitor] | None = None, min_supprt_reads=1):
        self.nlgraph = nlgraph
        self.visitors: list[GraphVisitor] = [] if visitors is None else visitors
        self.min_supprt_reads = min_supprt_reads

    @classmethod
    def from_visitors(cls, nlgraph: NLGraph, visitors: list[GraphVisitor], min_supprt_reads) -> GraphVis:
        """Create GraphVis from visitors."""
        return cls(nlgraph, visitors, min_supprt_reads)

    @staticmethod
    def load(file_name: str | Path):
        if isinstance(file_name, str):
            file_name = Path(file_name)

        if not file_name.exists():
            msg = f"{file_name} not exists."
            raise FileNotFoundError(msg)

        with file_name.open() as f:
            data = json.load(f)
            return nx.node_link_graph(data, **GraphVis.GRAPH_LINK_DATA)  # type: ignore

    @staticmethod
    def load_cytoscape(file_name: str | Path) -> nx.Graph:
        if isinstance(file_name, str):
            file_name = Path(file_name)
        if not file_name.exists():
            msg = f"{file_name} not exists."
            raise FileNotFoundError(msg)
        with file_name.open() as f:
            data = json.load(f)
            return nx.cytoscape_graph(data)

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
        return f"{node.chrom}_{node.ref_start}_{node.ref_end}_{head_node}"

    @staticmethod
    def add_node_to_graph(node: Node, graph: nx.Graph) -> None:
        """Add node to graph."""
        node_label = GraphVis.get_label_from_node(node)
        graph.add_node(
            node_label,
            chrom=node.chrom,
            ref_start=node.ref_start,
            ref_end=node.ref_end,
            strand=str(node.strand),
            is_head=node.is_start_node(),
            trace_id=node.trace_id,
        )

    @staticmethod
    def add_edge_to_graph(
        node1: Node,
        node2: Node,
        edge: Edge,
        graph: nx.Graph,
    ) -> None:
        """Add edge to graph."""
        edge_label = f"{edge.variation_type}_{edge.insertion_info}_{edge.sr}"
        node1_label = GraphVis.get_label_from_node(node1)
        node2_label = GraphVis.get_label_from_node(node2)

        if graph.has_edge(node1_label, node2_label):
            current_edge_label = [i["label"] for i in graph[node1_label][node2_label].values()]
            if edge_label not in current_edge_label:
                graph.add_edge(
                    node1_label,
                    node2_label,
                    label=edge_label,
                    weight=edge.sr,
                    read_ids=edge.read_ids,
                )
        else:
            graph.add_edge(
                node1_label,
                node2_label,
                label=edge_label,
                weight=edge.sr,
                read_ids=edge.read_ids,
            )

    def create_nxgraph(self) -> nx.DiGraph:
        # https://networkx.org/documentation/stable/reference/classes/multidigraph.html

        g = nx.MultiDiGraph()

        try:
            for start_node in self.nlgraph.get_start_nodes():
                self._traverse_graph(start_node, [start_node], g, self.nlgraph)  # type: ignore
        except RecursionError:
            logger.error("RecursionError: maximum recursion depth exceeded when export graph")

        return g

    def _traverse_graph(
        self,
        start_node: Node,
        path,
        nx_graph: nx.Graph,
        graph,
    ) -> None:
        """Plot graph helper."""
        if not start_node:
            return

        self.add_node_to_graph(start_node, nx_graph)

        if successors := start_node.successors:
            for successor in successors:
                for edge in graph.get_possible_edges(
                    path,
                    start_node,
                    successor,
                    self.min_supprt_reads,  # minimal support_reads,
                    filter_edges=False,
                ):
                    self.add_node_to_graph(successor, nx_graph)
                    self.add_edge_to_graph(start_node, successor, edge, nx_graph)
                    self._traverse_graph(
                        successor,
                        [*path, edge, successor],
                        nx_graph,
                        graph,
                    )


# https://networkx.org/documentation/latest/auto_examples/drawing/plot_weighted_graph.html#sphx-glr-auto-examples-drawing-plot-weighted-graph-py
# https://networkx.org/documentation/latest/reference/drawing.html


class GraphVisitor(ABC):
    @abstractmethod
    def visit(self, graph: nx.DiGraph, *args, **kwargs) -> None:
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


class GraphJsonExporter(GraphVisitor):
    def __init__(self, file_name: str | Path, graph_link_data):
        if isinstance(file_name, str):
            self.file_name = Path(file_name)

        self.graph_link_data = graph_link_data

    def visit(self, graph: nx.DiGraph):
        data = nx.node_link_data(graph, **self.graph_link_data)
        with Path(f"graph_{self.file_name}.json").open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=4)


class GraphCytoscapeExporter(GraphVisitor):
    def __init__(self, file_name: str | Path):
        if isinstance(file_name, str):
            self.file_name = Path(file_name)

    def visit(self, graph: nx.DiGraph):
        data = nx.cytoscape_data(graph)
        with Path(f"{self.file_name}_cy.json").open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=4)


def _cal_figure_size(nodes_size: int):
    unit = 5
    max_size = 25
    size = min(nodes_size * unit, max_size)
    return (size, size)


def draw_multiedge_labels(graph, pos):
    edge_labels = {}
    for u, v, _key, data in graph.edges(keys=True, data=True):
        if (u, v) not in edge_labels:
            edge_labels[(u, v)] = []
        edge_labels[(u, v)].append(str(data["label"]))

    for (u, v), labels in edge_labels.items():
        label = "\n".join(labels)
        x_pos = (pos[u][0] + pos[v][0]) / 2
        y_pos = (pos[u][1] + pos[v][1]) / 2
        plt.text(x_pos, y_pos, label, horizontalalignment="center")


def visualize_graph_via_matplot(
    graph: nx.DiGraph,
    figure_name: str,
    support_reads: int,
) -> None:
    node_numbers = len(list(graph))

    _fig, ax = plt.subplots(figsize=_cal_figure_size(node_numbers))

    pos = nx.spring_layout(graph, seed=42)

    edge_weight = nx.get_edge_attributes(graph, "weight")

    options = {
        "font_size": 13,
        "node_size": 1500,
        "node_color": ["red" if "H" in n else "white" for n in graph],
        "edgecolors": "black",
        "edge_color": ["green" if weight >= support_reads else "black" for weight in edge_weight.values()],
        "linewidths": 2,
        "width": 3,
        "connectionstyle": "arc3, rad = 0.1",
    }

    nx.draw_networkx(graph, pos=pos, arrows=True, **options)
    draw_multiedge_labels(graph, pos)

    ax.set_title(f"Node number: {node_numbers}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f"graph_{figure_name}.png")
    plt.close()


def visualize_graph_via_pyvis(graph: nx.Graph, figure_name: str | Path) -> None:
    from pyvis.network import Network  # type: ignore

    nt = Network(height="750px", directed=True, width="100%")
    nt.from_nx(graph)
    nt.save_graph(f"graph_{figure_name}.html")


def visualize_graph_via_d3(graph: nx.Graph, figure_name: str | Path) -> None:
    # https://d3blocks.github.io/d3blocks/pages/html/d3graph.html
    # https://docs.bokeh.org/en/latest/docs/examples/topics/graph/from_networkx.html
    # https://github.com/d3blocks/d3blocks
    # https://docs.bokeh.org/en/latest/docs/examples/topics/graph/node_and_edge_attributes.html
    # https://observablehq.com/@d3/force-directed-graph/2?intent=fork
    # https://observablehq.com/@d3/gallery?utm_source=d3js-org&utm_medium=nav&utm_campaign=try-observable
    raise NotImplementedError
