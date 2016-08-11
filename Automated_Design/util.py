import networkx as nx


def generate_graph(num_vert, edges, edge_length_vec=None):
    G = nx.Graph()
    for n in range(num_vert):
        G.add_node(n)
    for edge, length in zip(edges, edge_length_vec):
        if edge_length_vec:
            G.add_edge(edge[0], edge[1], length=length)
        else:
            G.add_edge(edge[0], edge[1])
        # TODO: verfiy the length property is not used as 'weight' in spanning tree alg.

    full_graph = G
    return full_graph