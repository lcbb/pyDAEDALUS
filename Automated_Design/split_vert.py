from Automated_Design.arrange_neighbors import arrange_neighbors


def split_vert(edge_type_mat_wHalfs, pseudo_vert, num_vert, vert_to_face):
    # Splits vertices into N nodes, where N is the degree of the vertex
    # Inputs: edge_type_mat_wHalfs = VxV sparse matrix where
    #  -1 is half of a non-spanning tree edge (one side of scaffold crossover)
    #   2 is spanning tree edge: DX edge with 0 scaffold crossovers
    #         pseudo_vert = row vector where value j at index i indicate that
    #            vertex i corresponds to vertex j, one of the V real vertices
    #         num_vert = number of vertices, V
    #         coordinates = Vx3 matrix of spatial coordinates of vertices
    #         vert_to_face = Vx1 cell array, each row has a row vector listing
    #             the face IDs the particular vertex belongs to
    # Outputs: edge_type_mat_allNodes = sparse matrix where
    #  -1 is half of a non-spanning tree edge (one side of scaffold crossover)
    #   2 is spanning tree edge: DX edge with 0 scaffold crossovers
    #        pseudo_vert = updated to include new pseudo-vertices
    ###########################################################################
    # by Sakul Ratanalert, MIT, Bathe Lab, 2016
    #
    # Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    # M.I.T. hereby makes following copyrightable material available to the
    # public under GNU General Public License, version 2 (GPL-2.0). A copy of
    # this license is available at https://opensource.org/licenses/GPL-2.0
    ###########################################################################

    # # Initialize edge_type_mat_allNodes. It will be augmented from _wHalfs
    edge_type_mat_allNodes = edge_type_mat_wHalfs.copy()

    # # Initialize face_storage matrix for pseudo nodes added in this function
    face_assign = [None] * len(pseudo_vert)

    for vert_ID in range(num_vert):  # for each real vertex
        neighbors = edge_type_mat_allNodes.in_edges(vert_ID, data=True)

        # # arrange_neighbors identifies which neighboring vertices share faces
        # # and adds nodes at the vertex to route the scaffold around the
        # # vertex vert_ID accordingly
        edge_type_mat_allNodes, pseudo_vert, face_assign = arrange_neighbors(
            edge_type_mat_allNodes, pseudo_vert, vert_ID, neighbors,
            vert_to_face, face_assign)
    return edge_type_mat_allNodes, pseudo_vert
