from .nlgraph import NLGraph, Node


def get_neighbors(node: Node) -> list[Node]:
    return node.predecessors + node.successors


def is_weakly_connected(graph: NLGraph):
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
