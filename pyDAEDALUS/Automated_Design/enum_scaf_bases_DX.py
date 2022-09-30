import math


def enum_scaf_bases_DX(route_real, route_vals, edge_length_mat_full, Aform=False):
    """
    Enumerates scaffold bases for N-arm DX tile-based cages
    Inputs: route_real = row vector of vertices listed in visitation order
                (only real vertex IDs)
             route_vals = row vector of edge types, where the value at index
                j in route_vals is the edge type of the edge between the
                vertices route_real[j:j+1), wrapping around at end
            edge_length_mat_full = VxV sparse matrix of edge lengths
    Outputs: edge_bgn_vec = row vector of scaff nt IDs at which edge begins
             edge_fin_vec = row vector of scaff nt IDs at which edge finishes
             edge_type_vec = row vector of edge types, corresponding to
                             edge_length_mat_full
      2 is spanning tree edge: DX edge with 0 scaffold crossovers
     -3 is half of a non-spanning tree edge, connecting to vertex at 3' end
     -5 is half of a non-spanning tree edge, connecting to vertex at 5' end
    ##########################################################################
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https://opensource.org/licenses/GPL-2.0
    ##########################################################################
    """
    #
    # If even: EdgeLength/2-8
    #          EdgeLength/2-3
    # Elif odd: ceil(EdgeLength/2)+2
    #           floor(EdgeLength/2)-2
    #
    # # Initialize output variables and length of route
    len_route = len(route_real)  # length of route
    edge_bgn_vec = []
    edge_fin_vec = []
    edge_type_vec = []

    for route_ID in range(len_route):  # for each vertex in the route
        edge_bgn = route_real[route_ID]  # current edge

        # # Find edge_fin, the other end of the edge
        if route_ID == len_route - 1:  # if last vertex...
            edge_fin = route_real[0]  # use first vertex in route.
        else:  # else use next vertex
            edge_fin = route_real[route_ID + 1]

        # # Find edge_before, the vertex visited prior to edge_bgn
        if route_ID == 0:  # if first vertex...
            edge_before = route_real[-1]  # use last vertex in route.
        else:  # else use previous vertex
            edge_before = route_real[route_ID - 1]

        # Get edge type: 2 (no scafold crossover) or -1 (half of an edge
        # with 1 scaffold crossover)
        route_val = route_vals[route_ID]

        len_edge = edge_length_mat_full[edge_bgn][edge_fin]['length']

        if route_val == 2:  # tree edge, edge with no scaff crossovers
            len_edge = int(len_edge)  # len_edge should be int-able by now
            edge_bgn_vec = edge_bgn_vec + ([edge_bgn] * len_edge)
            edge_fin_vec = edge_fin_vec + ([edge_fin] * len_edge)
            edge_type_vec = edge_type_vec + ([route_vals[route_ID]] * len_edge)

        elif route_val == -1:  # half non-tree edge
            if Aform: # H-form
                # # Is it an even or odd multiple of 11? (Is the length even/odd?)
                quotient = round(len_edge / 11)
                quot_is_odd = quotient % 2  # 1 if odd, 0 if even multiple

                # # Adjust len_edge according to parity of len and quotient
                if quot_is_odd:  # quotient is odd, scaff cross in middle +/-2
                    len_edge_adjusted = len_edge

                    # if returning half-strand, 3' end, 1 shorter
                    if edge_before == edge_fin:
                        route_val = -3  # -3 = returning half-strand
                        len_half = math.ceil(len_edge_adjusted / 2.0) # changed
                    else:  # outgoing half-strand, 5' end, 1 longer
                        route_val = -5  # -5 = outgoing half-strand
                        len_half = math.floor(len_edge_adjusted / 2.0) # changed

                else:  # quotient is even, scaffold 3/8 bp from middle
                        # len is even, outgoing + returning = len_edge +/-11
                        # if returning half-strand, 3' end, 1 shorter:
                    if edge_before == edge_fin:
                        route_val = -3  # -3 = returning half-strand

                        if edge_bgn > edge_fin:  # shorter half
                            len_edge_adjusted = len_edge - 10
                        else:  # longer half
                            len_edge_adjusted = len_edge + 10

                        len_half = math.floor(len_edge_adjusted / 2.0) # changed

                    else:  # outgoing half-strand, 5' end, 1 longer
                        route_val = -5  # -5 = outgoing half-strand

                        if edge_bgn < edge_fin:  # shorter half
                            len_edge_adjusted = len_edge - 10
                        else:  # longer half
                            len_edge_adjusted = len_edge + 10

                        len_half = math.ceil(len_edge_adjusted / 2.0)

            else: # B-form
                # # Is it an even or odd multiple of 10.5? Is the length even/odd?
                quotient = round(len_edge / 10.5)
                len_is_odd = len_edge % 2  # 1 if odd, 0 if even length
                quot_is_odd = quotient % 2  # 1 if odd, 0 if even multiple

                # # Adjust len_edge according to parity of len and quotient
                if quot_is_odd:  # quotient is odd, scaff cross in middle
                    if len_is_odd:  # outgoing + returning = len_edge
                        len_edge_adjusted = len_edge

                        # if returning half-strand, 3' end, 1 shorter
                        if edge_before == edge_fin:
                            route_val = -3  # -3 = returning half-strand
                            len_half = math.floor(len_edge_adjusted / 2.0)
                        else:  # outgoing half-strand, 5' end, 1 longer
                            route_val = -5  # -5 = outgoing half-strand
                            len_half = math.ceil(len_edge_adjusted / 2.0)

                    else:  # len is even, outgoing + returning = len_edge +/-1
                        # if returning half-strand, 3' end, 1 shorter
                        if edge_before == edge_fin:
                            route_val = -3  # -3 = returning half-strand

                            if edge_bgn > edge_fin:  # shorter half-tile
                                len_edge_adjusted = len_edge - 1
                            else:  # longer half
                                len_edge_adjusted = len_edge + 1

                            len_half = math.floor(len_edge_adjusted / 2.0)

                        else:  # outgoing half-strand, 5' end, 1 longer
                            route_val = -5  # -5 = outgoing half-strand

                            if edge_bgn < edge_fin:  # shorter half-tile
                                len_edge_adjusted = len_edge - 1
                            else:  # longer half-tile
                                len_edge_adjusted = len_edge + 1

                            len_half = math.ceil(len_edge_adjusted / 2.0)

                else:  # quotient is even, scaffold 5/6 bp from middle
                    if len_is_odd:  # outgoing + returning = len_edge +/-10
                        # if returning half-strand, 3' end, 1 shorter
                        if edge_before == edge_fin:
                            route_val = -3  # -3 = returning half-strand

                            if edge_bgn > edge_fin:  # shorter half
                                len_edge_adjusted = len_edge - 10
                            else:  # longer half
                                len_edge_adjusted = len_edge + 10

                            len_half = math.floor(len_edge_adjusted / 2.0)

                        else:  # outgoing half-strand, 5' end, 1 longer
                            route_val = -5  # -5 = outgoing half-strand

                            if edge_bgn < edge_fin:  # shorter half
                                len_edge_adjusted = len_edge - 10
                            else:  # longer half
                                len_edge_adjusted = len_edge + 10

                            len_half = math.ceil(len_edge_adjusted / 2.0)

                    else:  # len is even, outgoing + returning = len_edge +/-11
                        # if returning half-strand, 3' end, 1 shorter:
                        if edge_before == edge_fin:
                            route_val = -3  # -3 = returning half-strand

                            if edge_bgn > edge_fin:  # shorter half
                                len_edge_adjusted = len_edge - 11
                            else:  # longer half
                                len_edge_adjusted = len_edge + 11

                            len_half = math.floor(len_edge_adjusted / 2.0)

                        else:  # outgoing half-strand, 5' end, 1 longer
                            route_val = -5  # -5 = outgoing half-strand

                            if edge_bgn < edge_fin:  # shorter half
                                len_edge_adjusted = len_edge - 11
                            else:  # longer half
                                len_edge_adjusted = len_edge + 11

                            len_half = math.ceil(len_edge_adjusted / 2.0)

            len_half = int(len_half)  # len_half should be int-able by now

            # # Add to edge_{bgn,fin,type}_vec
            edge_bgn_vec = edge_bgn_vec + ([edge_bgn] * len_half)
            edge_fin_vec = edge_fin_vec + ([edge_fin] * len_half)
            edge_type_vec = edge_type_vec + ([route_val] * len_half)

    return [edge_bgn_vec, edge_fin_vec, edge_type_vec]
