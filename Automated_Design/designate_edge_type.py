from networkx.algorithms.mst import prim_mst

def designate_edge_type(full_graph):
    """
    ...
    1 is non-spanning tree edge: DX edge with 1 scaffold crossover
    2 is spanning tree edge: DX edge with 0 scaffold crossovers
    """
    tree = prim_mst(full_graph)

    full_tree = tree #networkx undirected networks already implicitly can be considered a symmetrical matrix

    #TODO: is `edge_type_mat = full_tree + full_graph` a no-op?  if not, what's the equivalent in this representation?
    tree_edges = full_tree.edges()
    for edge in full_graph.edges():
        print edge
        if edge in tree_edges:
            type = 2
        else:
            type = 1
        i = edge[0]
        j = edge[1]
        full_graph[i][j]['type'] = type


    return full_graph
