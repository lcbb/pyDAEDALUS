def arrange_neighbors(edge_type_mat_allNodes, pseudo_vert, vert_ID, neighbors,
                      vert_to_face, face_assign):
    """
    Identifies which neighboring vertices share faces (which neighbors are
    adjacent and adds N nodes at the vertex to route the scaffold around the
    vertex vert_ID accordingly, where N is the degree of the vertex vert_ID.
    Inputs: edge_type_mat_wHalfs = VxV sparse matrix where
     -1 is half of a non-spanning tree edge (one side of scaffold crossover)
      2 is spanning tree edge: DX edge with 0 scaffold crossovers
            pseudo_vert = row vector where value j at index i indicate that
               vertex i corresponds to vertex j, one of the V real vertices
            vert_ID = ID of the central vertex that neighbors surround
            neighbors = Vx1 column vector of vertices neighboring vert_ID
               0 indicates not a neighbor, nonzero value indicates neighbor
            vert_to_face = Vx1 cell array, each row has a row vector listing
                the face IDs the particular vertex belongs to
            face_assign = vector that stores face_ID of pseudo-vertices
    Outputs: edge_type_mat_allNodes = VxV sparse matrix where
     -1 is half of a non-spanning tree edge (one side of scaffold crossover)
      2 is spanning tree edge: DX edge with 0 scaffold crossovers
           pseudo_vert = updated to include new pseudo-vertices
           face_assign = vector that stores face_ID of pseudo-vertices
    ##########################################################################
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https://opensource.org/licenses/GPL-2.0
    ##########################################################################
    """

    n_rows, dontcare, n_vals = zip(*neighbors)  # n for neighbors

    # # Identify degree of vertex vert_ID, number of vertices connected to
    # # central vertex. This includes pseudo-vertices, so degree may be > N
    degree = len(n_rows)

    # # Identify neighboring pairs and add a new node between each pair
    list_a = range(degree - 1)
    for row_ID1 in list_a:

        # # Find neighbor 1 info
        neighbor_1 = n_rows[row_ID1]  # extract ID of neighbor
        val_1 = n_vals[row_ID1]  # extract edge type connecting to neighbor
        face_1 = face_assign[neighbor_1]  # designated face, if assigned

        list_b = range(row_ID1 + 1, degree)
        for row_ID2 in list_b:

            # # Find neighbor info
            neighbor_2 = n_rows[row_ID2]  # extract ID of neighbor
            val_2 = n_vals[row_ID2]  # extract edge type connecting to neighbor
            face_2 = face_assign[neighbor_2]  # designated face, if assigned

            # # If these potential neighbors do not share the same real vertex,
            # # i.e. are not the same vertex with different pseudo-vert names
            # # And if this pair hasn't already been checked
            if pseudo_vert[neighbor_1] != pseudo_vert[neighbor_2]:
                # # If these potential neighbors share a face with vert_ID and
                # # each other, then they are neighbors
                v2f_vert_ID = vert_to_face[vert_ID]
                v2f_n1 = vert_to_face[pseudo_vert[neighbor_1]]
                v2f_n2 = vert_to_face[pseudo_vert[neighbor_2]]

                shared_neighbor = set(v2f_vert_ID)\
                    .intersection(set(v2f_n1))\
                    .intersection(set(v2f_n2))
                # TODO: there has to be a way to clean up at
                # least the for-for-if above...

                if shared_neighbor:
                    # # Record the face this new node belongs to

                    face_3 = list(shared_neighbor)[0]

                    # # Check if connecting to another pseudo node if faces
                    # # match, otherwise don't add node
                    #  if face_1 == face_2 == face_3, but
                    # allowing for empty face_1 or face_2:
                    if (face_1 == None or face_1 == face_3) and \
                            (face_2 == None or face_2 == face_3):  # noqa: E711

                        # # Create a new node and store in pseudo_vert
                        new_node_id = len(pseudo_vert)  # create a new node
                        edge_type_mat_allNodes.add_node(new_node_id)
                        pseudo_vert.append(vert_ID)  # store in pseudo_vert

                        # Incorporate new node into edge_type_mat_allNodes
                        edge_type_mat_allNodes.add_edge(
                            neighbor_1, new_node_id, attr_dict=val_1)
                        edge_type_mat_allNodes.add_edge(
                            new_node_id, neighbor_1, attr_dict=val_1)
                        edge_type_mat_allNodes.add_edge(
                            new_node_id, neighbor_2, attr_dict=val_2)
                        edge_type_mat_allNodes.add_edge(
                            neighbor_2, new_node_id, attr_dict=val_2)

                        # # Store in face_assign
                        face_assign.append(face_3)

    # Remove/disconnect vert_ID from edge_type_mat_allNodes, the real vertex
    # is no longer needed to route, having been replaced by pseudo-vertices
    edge_type_mat_allNodes.remove_edges_from(
        edge_type_mat_allNodes.in_edges(vert_ID))
    edge_type_mat_allNodes.remove_edges_from(
        edge_type_mat_allNodes.out_edges(vert_ID))

    # TODO: replace pseudo_vert with a property on the node `original_id`?
    return edge_type_mat_allNodes, pseudo_vert, face_assign
