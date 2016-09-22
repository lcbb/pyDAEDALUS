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
    # Adjust the scaffold nick position so that there are scaf_nick_pos bases
    # upstream of scaffold nick.
    # Inputs: scaf_to_edge = Ex2 cell array, where each row corresponds to one
    #            edge, 1st column is duplex from low ID to high ID vertex,
    #            2nd column is from high to low. Each element is a row vector
    #            containing the scaffold base IDs in order on that duplex.
    #         scaf_nick_pos = number of bases upstream of scaffold nick
    #         num_bases = total number of scaffold bases
    # Outputs: scaf_to_edge_adj = scaf_to_edge with adjusted numbering
    ###########################################################################
    # by Sakul Ratanalert, MIT, Bathe Lab, 2016
    #
    # Copyright 2016. Massachusetts Institute of Technology. Rights Reserved. 
    # M.I.T. hereby makes following copyrightable material available to the
    # public under GNU General Public License, version 2 (GPL-2.0). A copy of
    # this license is available at https://opensource.org/licenses/GPL-2.0
    ###########################################################################

    # # Initialize output vector
    scaf_to_edge_adj = scaf_to_edge   # row vec

    # # Adjust scaf_to_edge
    import numpy as np
    num_edges, num_sides, dontcare = np.array(scaf_to_edge).shape # num_sides should be 2

    for edge_ID in range(num_edges):
        for sides_ID in range(num_sides):
            scaf_to_edge_adj[edge_ID][sides_ID] = adjust(scaf_to_edge[edge_ID][sides_ID], scaf_nick_pos, num_bases)


    return scaf_to_edge_adj

def adjust(old_num_vec, scaf_nick_pos, num_bases):
    # Adjusts numbered vector by shifting by scaf_nick_pos bases, with
    # wraparound
    # Inputs: old_num_vec = numbered vector to shift
    #         scaf_nick_pos = number of bases to shift by
    #         num_bases = total number of bases in scaffold
    # Output: new_num_vec = adjusted numbered vector

    # Initialize
    new_num_vec = old_num_vec

    old_num_vec_len = len(old_num_vec)
    for i in range(old_num_vec_len):
        new_num = (int(old_num_vec[i]) - int(scaf_nick_pos)) % int(num_bases)

        if new_num == num_bases: # if is now last base
            new_num_vec[i] = num_bases  # set to last base
        else:
            new_num_vec[i] = new_num

    return new_num_vec
