"""Visualize Graphs."""

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


def default_visitors(
    graph: NLGraph,
    figure_name: str,
    support_reads: int = 1,
    possible_paths: dict[str, list[str]] | None = None,
) -> GraphVis:
    return GraphVis.from_visitors(
        graph,
        [
            GraphCytoscapeExporter(figure_name),
        ],
        support_reads,
        possible_paths,
    )


def get_label_from_node(node: Node) -> str:
    """Get label from node."""
    return node.id


def add_node_to_nxgraph(node: Node, graph: nx.Graph) -> None:
    """Add node to graph."""
    node_label = get_label_from_node(node)
    idendities = []
    for rid, identity in node.identities.items():
        if identity.is_head():
            identity_str = "SO"
        elif identity.is_tail():
            identity_str = "SI"
        else:
            identity_str = "IN"
        idendities.append(f"{rid}:{identity_str}")

    graph.add_node(
        node_label,
        chrom=node.chrom,
        ref_start=node.ref_start,
        ref_end=node.ref_end,
        strand=str(node.strand),
        is_head=node.is_start_node(),
        exons=str(node.exons),
        reads=",".join(idendities),
        ptc=node.ptc,
        ptf=node.ptf,
    )


def add_edge_to_nxgraph(
    node1: Node,
    node2: Node,
    edge: Edge,
    graph: nx.Graph,
) -> None:
    """Add edge to graph."""
    edge_label = edge.id
    node1_label = get_label_from_node(node1)
    node2_label = get_label_from_node(node2)
    link_type = edge.nclt_link_type
    breakpoints = f"{edge.break_point1.chrom},{edge.break_point2.chrom},{edge.break_point1.pos},{edge.break_point2.pos},{link_type}"

    link_attributes = edge.link_attributes(rescue_sr=False)
    mode1 = edge.mode1.to_str()
    mode2 = edge.mode2.to_str()

    insertion = "" if edge.insertion_info is None else f"{edge.insertion_info[1]}"

    if graph.has_edge(node1_label, node2_label):
        current_edge_label = [i["id"] for i in graph[node1_label][node2_label].values()]
        if edge_label not in current_edge_label:
            logger.warning(f"vis: multiple edges between {node1_label} and {node2_label}")
            graph.add_edge(
                node1_label,
                node2_label,
                id=edge.id,
                weight=edge.sr,
                read_ids=edge.read_ids,
                gene1=edge.gene1,
                gene2=edge.gene2,
                breakpoints=breakpoints,
                insertion_info=insertion,
                mode1=mode1,
                mode2=mode2,
                **link_attributes._asdict(),
            )
    else:
        graph.add_edge(
            node1_label,
            node2_label,
            id=edge.id,
            weight=edge.sr,
            read_ids=edge.read_ids,
            gene1=edge.gene1,
            gene2=edge.gene2,
            breakpoints=breakpoints,
            insertion_info=insertion,
            mode1=mode1,
            mode2=mode2,
            **link_attributes._asdict(),
        )


def create_nxgraph(
    nlgraph,
    min_support_reads: int = 1,
    possible_paths: dict[str, list[str]] | None = None,
) -> nx.DiGraph:
    # https://networkx.org/documentation/stable/reference/classes/multidigraph.html

    g = nx.MultiDiGraph() if possible_paths is None else nx.MultiDiGraph(possible_paths=possible_paths)

    try:
        for start_node in nlgraph.get_start_nodes():
            _create_nxgraph(start_node, [start_node], g, nlgraph, min_support_reads)  # type: ignore
    except RecursionError:
        logger.error("RecursionError: maximum recursion depth exceeded when export graph")
    return g


def _create_nxgraph(
    start_node: Node,
    path,
    nx_graph: nx.Graph,
    graph,
    min_support_reads,
) -> None:
    """Plot graph helper."""
    if not start_node:
        return

    add_node_to_nxgraph(start_node, nx_graph)

    if successors := start_node.successors:
        for successor in successors:
            for idx, edge in enumerate(
                graph.get_possible_edges(
                    path,
                    start_node,
                    successor,
                    min_support_reads,  # minimal support_reads,
                    filter_edges=False,
                )
            ):
                if idx > 0:
                    logger.warning("Vis: multiple edges between {} and {}", start_node, successor)

                add_node_to_nxgraph(successor, nx_graph)
                add_edge_to_nxgraph(start_node, successor, edge, nx_graph)
                _create_nxgraph(
                    successor,
                    [*path, edge, successor],
                    nx_graph,
                    graph,
                    min_support_reads,
                )


class GraphVis:
    GRAPH_LINK_DATA: ClassVar[dict[str, str]] = {
        "link": "edges",
        "source": "source",
        "target": "target",
    }

    def __init__(
        self,
        nlgraph: NLGraph,
        visitors: list[GraphVisitor],
        min_support_reads: int = 1,
        possible_paths: dict[str, list[str]] | None = None,
    ):
        self.nlgraph = nlgraph
        self.visitors: list[GraphVisitor] = [] if visitors is None else visitors
        self.min_support_reads = min_support_reads
        self.possible_paths = possible_paths

    @classmethod
    def from_visitors(
        cls,
        nlgraph: NLGraph,
        visitors: list[GraphVisitor],
        min_support_reads: int,
        possible_paths: dict[str, list[str]] | None = None,
    ) -> GraphVis:
        """Create GraphVis from visitors."""
        return cls(nlgraph, visitors, min_support_reads, possible_paths)

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
        g = create_nxgraph(self.nlgraph, self.min_support_reads, self.possible_paths)

        if visitors is None:
            visitors = []

        for visitor in self.visitors + visitors:
            visitor.visit(g)


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
            json.dump(data, f, ensure_ascii=False, indent=2)


class TSGraphExporter(GraphVisitor):
    def __init__(self, file_name: str | Path):
        if isinstance(file_name, str):
            self.file_name = Path(file_name)

    def visit(self, graph: nx.DiGraph):
        """Export graph to TSG format.

        # Header information
        H	TSG	1.0
        H	reference	GRCh38
        # Nodes
        N	n1	chr1:+:1000-1200,1500-1700	read1:SO,read2:SO	ACGTACGT
        N	n2	chr1:+:2000-2200	read4:SO,read5:SO	TGCATGCA
        N	n3	chr1:+:2500-2700	read1:IN,read2:IN,read3:IN,read4:IN	CTGACTGA
        N	n4	chr1:+:2500-2700	read1:SI,read2:SI	CTGACTGA
        N	n5	chr1:+:2500-2700	read3:SI,read4:SI	CTGACTGA
        # Edges
        E	e1	n1	n3	chr1,chr1,1700,2000,splice
        E	e2	n3	n4	chr1,chr1,1700,2000,splice
        E	e3	n2	n3	chr1,chr1,2200,2500,splice
        E	e4	n3	n5	chr1,chr1,1700,2500,splice
        # Chains (building the graph)
        C	chain1	n1	e1	n3	e2	n4
        C	chain2	n2	e3	n3  e4  n5
        # Paths (traversals through the constructed graph)
        P	transcript1	n1+	e1+	n3+	e2+	n4+
        P	transcript2	n2+	e3+	n3+ e4+ n5+
        # Sets (grouping elements)
        U	exon_set	n1	n2	n3
        # Attributes (metadata)
        A	N	n1	expression:f:10.5
        A	O	transcript1	tpm:f:8.2
        A	O	transcript2	tpm:f:3.7
        """
        # write header
        with Path(f"{self.file_name}.tsg").open("w", encoding="utf-8") as f:
            f.write("H\tTSG\t1.0\n")
            f.write("H\treference\tGRCh38\n")

            # write nodes
            for node in graph.nodes(data=True):
                f.write(f"N\t{node[0]}\t{node[1]['chrom']}:{node[1]['strand']!s}:{node[1]['exons'][1:-1]!s}\t{node[1]['reads']}\n")

            # write edges
            for edge in graph.edges(data=True):
                f.write(f"E\t{edge[2]['id']}\t{edge[0]}\t{edge[1]}\t{edge[2]['breakpoints']}\n")

            # write node attributes sr
            for node in graph.nodes(data=True):
                f.write(f"A\tN\t{node[0]}\tptc:i:{node[1]['ptc']}\n")
                f.write(f"A\tN\t{node[0]}\tptf:f:{node[1]['ptf']}\n")

            # write edge attributes sr
            for edge in graph.edges(data=True):
                f.write(f"A\tE\t{edge[2]['id']}\tsr:i:{edge[2]['weight']}\n")

                # write insertion info if exists
                if "insertion_info" in edge[2]:
                    insertion_info = edge[2]["insertion_info"]

                    if insertion_info != "":
                        insertion_type, insertion_seq = insertion_info.split("(")
                        if insertion_type == "NovelInsertion":
                            seq = insertion_seq.strip(")").split(":")[0]
                            f.write(f"A\tE\t{edge[2]['id']}\tnovel_insertion:Z:{seq}\n")
                        elif insertion_type == "MicroHomology":
                            seq = insertion_seq.strip(")")
                            f.write(f"A\tE\t{edge[2]['id']}\tmicrohomology:Z:{seq}\n")


def _cal_figure_size(nodes_size: int):
    unit = 5
    max_size = 25
    size = min(nodes_size * unit, max_size)
    return (size, size)


def draw_multiedge_labels(graph, pos):
    edge_labels = {}
    for u, v, _key, data in graph.edges(keys=True, data=True):
        if (u, v) not in edge_labels:
            edge_labels[u, v] = []
        edge_labels[u, v].append(str(data["label"]))

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
