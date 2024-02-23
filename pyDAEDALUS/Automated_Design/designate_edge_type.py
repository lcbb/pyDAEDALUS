from networkx.algorithms import tree

def designate_edge_type(full_graph):
    """
    Calculate minimum spanning tree, and label the links in the full graph
    that are also in the min spanning tree as such:
        1 is non-spanning tree edge: DX edge with 1 scaffold crossover
        2 is spanning tree edge: DX edge with 0 scaffold crossovers
    """
    mst = tree.minimum_spanning_tree(full_graph, algorithm="prim")
    edgelist = mst.edges
    
    for edge in full_graph.edges():
        if edge in edgelist:
            edgetype = 2
        else:
            edgetype = 1
        i = edge[0]
        j = edge[1]
        full_graph[i][j]['type'] = edgetype
    
    return full_graph