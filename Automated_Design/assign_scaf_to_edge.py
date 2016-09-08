def find(iterable, val):
    return [i for i in range(len(iterable)) if iterable[i] == val]

def intersect(a, b):
    return sorted(list(set(a).intersection(set(b))))

def assign_scaf_to_edge(edges, num_edges, edge_type_mat, edge_bgn_vec,
                        edge_fin_vec, edge_type_vec):
    """
    # Assign enumerated scaffold bases to edges. Create vectors for each duplex
    # in each edge, identify to which scaffold base each edge base corresponds
    # Inputs: edges = Ex2 matrix where each row corresponds to one edge,
    #           denoting the vertices being connected. 1st column > 2nd column
    #         num_edges = number of edges, E
    #         edge_length_mat_full = VxV sparse matrix of edge lengths
    #         edge_bgn_vec = row vector of scaff nt IDs at which edge begins
    #         edge_fin_vec = row vector of scaff nt IDs at which edge finishes
    #         edge_type_vec = row vector of edge types, corresponding to
    #                          edge_length_mat_full
    #   2 is spanning tree edge: DX edge with 0 scaffold crossovers
    #  -3 is half of a non-spanning tree edge, connecting to vertex at 3' end
    #  -5 is half of a non-spanning tree edge, connecting to vertex at 5' end
    # Outputs: scaf_to_edge = Ex2 cell array, where each row corresponds to one
    #            edge, 1st column is duplex from low ID to high ID vertex,
    #            2nd column is from high to low. Each element is a row vector
    #            containing the scaffold base IDs in order on that duplex.
    ###########################################################################
    # by Sakul Ratanalert, MIT, Bathe Lab, 2016
    #
    # Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    # M.I.T. hereby makes following copyrightable material available to the
    # public under GNU General Public License, version 2 (GPL-2.0). A copy of
    # this license is available at https://opensource.org/licenses/GPL-2.0
    ###########################################################################
    """

    scaf_to_edge = []
    for edge_ID in range(num_edges):  #TODO: convert to `for edge in edges:`
        # first column low to high, second column high to low
        row = [None, None]
        for high_to_low in [1, 2]:  # for each duplex direction on edge
            col = 2 - high_to_low

            if high_to_low == 1:  # high to low 5' to 3'
                edge_bgn = edges[edge_ID][0]
                edge_fin = edges[edge_ID][1]
            else:  # low_to_high 5' to 3'
                edge_bgn = edges[edge_ID][1]
                edge_fin = edges[edge_ID][0]

            edge_type = edge_type_mat[edge_bgn][edge_fin]['type']
            if edge_type == 2:  # tree edge  #TODO: extract into constant
                bases = intersect(find(edge_bgn_vec, edge_bgn),
                                  find(edge_fin_vec, edge_fin))
            else:  # non-tree edge
                bases_all = intersect(find(edge_bgn_vec, edge_bgn),
                                      find(edge_fin_vec, edge_fin))
                bases_5 = intersect(bases_all, find(edge_type_vec, -5))
                bases_3 = intersect(bases_all, find(edge_type_vec, -3))

                bases = bases_5 + bases_3

            row[col] = bases

        scaf_to_edge.append(row)

    return scaf_to_edge