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

    # TODO: here's the matlab code... ambiguous if the above is enough or not until I get further into the rest of the code
    # graph = sparse(edges(:, 1), edges(:, 2), ones(num_edges, 1), num_vert, num_vert)
    # edge_length_mat = sparse(edges(:, 1), edges(:, 2), edge_length_vec, num_vert, num_vert)

    # Convert lower triangle graphs to full symmetric graphs
    # full_graph = graph + graph'
    # edge_length_mat_full = edge_length_mat + edge_length_mat'
    full_graph = G
    return full_graph