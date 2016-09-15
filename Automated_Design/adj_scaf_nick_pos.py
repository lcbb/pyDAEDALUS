import numpy as np

from Automated_Design.util import intersect_lists


def get_scaf_nick_pos(edges, route_real, edge_length_vec):
    first_node_in_route = route_real[0]
    second_node_in_route = route_real[1]
    min_node = min(first_node_in_route, second_node_in_route)
    max_node = max(first_node_in_route, second_node_in_route)

    edges = np.array(edges)
    from_nodes = np.where(edges[:, 0] == max_node)[0]
    to_nodes = np.where(edges[:, 1] == min_node)[0]

    first_edge = intersect_lists(from_nodes, to_nodes)[0]  #TODO: assured to always have 1 result?

    first_edge_length = edge_length_vec[first_edge]
    if first_edge_length < 42:  # == 31
        scaf_nick_pos = 16
    else:
        scaf_nick_pos = 19

    return scaf_nick_pos


def adj_scaf_nick_pos(scaf_to_edge, scaf_nick_pos, num_bases):

    scaf_to_edge_adj = None

    return scaf_to_edge_adj
