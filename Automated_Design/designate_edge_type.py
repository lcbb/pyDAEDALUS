from networkx.algorithms.mst import prim_mst

def designate_edge_type(full_graph):
    """
    Calculate minimum spanning tree, and label the links in the full graph that are also in the min spanning tree as such:
        1 is non-spanning tree edge: DX edge with 1 scaffold crossover
        2 is spanning tree edge: DX edge with 0 scaffold crossovers
    """
    tree = prim_mst(full_graph)

    full_tree = tree #networkx undirected networks already implicitly can be considered a symmetrical matrix

    tree_edges = full_tree.edges()
    for edge in full_graph.edges():
        if edge in tree_edges:
            type = 2
        else:
            type = 1
        i = edge[0]
        j = edge[1]
        full_graph[i][j]['type'] = type


    return full_graph
