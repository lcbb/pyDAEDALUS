import numpy as np
from matplotlib import pyplot as plt

from Automated_Design.contsants import REDPURPLE, VERMILLION, WHITE, SKYBLUE
from Automated_Design.util import generate_graph


def gen_schlegel(edges, coordinates, faces, edge_type_mat=None, schlegel_filename=None):
    """
    Plots vertices and edges in Schlegel diagram, tracing scaffold path in
    colored lines
    Inputs: edges = Ex2 matrix where each row corresponds to one edge,
              denoting the vertices being connected. 1st column > 2nd column
              E = number of edges
            coordinates = Vx3 matrix of spatial coordinates of vertices,
              V = number of vertices
            faces = Fx2 cell matrix, where F is the number of faces.        #TODO: Polish up these 'faces' docs
                The first column details how many vertices the face has
                The second column details the vertex IDs of the face
            edge_type_mat = sparse matrix where
      1 is non-spanning tree edge: DX edge with 1 scaffold crossover
      2 is spanning tree edge: DX edge with 0 scaffold crossovers
    Output: figure of vertices and edges in Schlegel diagram and coordinates
    of projected vertices (xycoord, Vx2 matrix).  Saved into folder
    `results_filename` if defined, else it's shown on screen
    ##########################################################################
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https:##opensource.org#licenses#GPL-2.0
    ##########################################################################
    """

    # Choose Biggest Face:
    big_face = max(faces, key=len)

    # Initialize:
    num_face_vert = len(big_face)  # for best results choose most vertices per face
    num_vert = len(coordinates)  # number of vertices
    num_edges = len(edges)  # number of edges

    # If edge_type_mat isn't specified, generate it.
    if edge_type_mat is None:
        edge_type_mat = generate_graph(num_vert=num_vert, edges=edges,)

    # Initialize xycoord
    xycoord = np.zeros(shape=(num_vert, 2))
    num_repeat = 10 * num_vert  # number of iterations to calculate xycoord

    # Set big face on unit circle
    angle = 2.*np.pi/num_face_vert

    for i in range(num_face_vert):  # for each vertex in big face
        xycoord[big_face[i]] = [np.cos(angle * i), np.sin(angle * i)]


    # # Calculate positions of other vertices through iterative process
    for repeat_ID in range(num_repeat):

        for vert_ID in range(num_vert): # For each vertex,
            if vert_ID not in big_face:  # that's not in big_face,
                # move the xycoord of the given vertex to the centerpont of all neighboring nodes.

                # # Find vert_ID
                find_vert_row, find_vert_col = np.where(edges == vert_ID)

                # # Find neighbors
                neighbors = []
                for row, col in zip(find_vert_row, find_vert_col):
                    find_vert_col_nbr = 1 - col
                    neighbors.append(edges[row, find_vert_col_nbr])
                neighbors = np.array(neighbors, dtype=np.int)

                # # Get neighbor coords
                nbr_coord = xycoord[neighbors, :]

                # # # Find midpoint
                num_nbr = len(find_vert_col)
                midpt = np.array([sum(nbr_coord[:, 0]), sum(nbr_coord[:, 1])]) / num_nbr


                # # Store as new xycoord for vert_ID
                xycoord[vert_ID, :] = midpt

    f = plt.figure(2, figsize=(8, 8))
    f.clf()
    plt.xlim((-1.2, 1.2))
    plt.ylim((-1.2, 1.2))
    #TODO: `axis equal` resolvede by figsize?
    plt.axis('off')

    # TODO:  `for thickfirst = 1:2` needed?  Yes.  Plot blue lines then pink line.
    for edge in edges:
        # Get end vertices
        edge_bgn, edge_fin = edge

        # Get end vertex coordinates
        x_bgn = xycoord[edge_bgn, 0]
        y_bgn = xycoord[edge_bgn, 1]

        x_fin = xycoord[edge_fin, 0]
        y_fin = xycoord[edge_fin, 1]

        # Create vectors from beginning vertex to end
        x_vec = [x_bgn, x_fin]
        y_vec = [y_bgn, y_fin]

        # TODO rename `edge_type_mat` to something that involves `graph`
        if edge_type_mat.edge[edge_bgn][edge_fin]['type'] == 2:  #if this is a tree edge:
            plt.plot(x_vec, y_vec, linestyle='-', color=REDPURPLE, linewidth=8)
        else:  #if this is a non-tree edge:
            plt.plot(x_vec, y_vec, linestyle='-', color=SKYBLUE, linewidth=5)

    # Plot nodes:
    plt.plot(xycoord[:, 0], xycoord[:, 1], 'o',
             color=VERMILLION, markeredgecolor=VERMILLION, markersize=15)

    # Plot node labels:
    for i, node_xy in enumerate(xycoord):
        x, y = node_xy
        plt.text(x, y, str(i), color=WHITE, fontname='serif',
                 horizontalalignment='center', verticalalignment='center')

    if schlegel_filename:
        plt.savefig(schlegel_filename, bbox_inches='tight')
    else:
        #Should I make showing an option?  Either way, should I also save the earlier distribution plots?
        plt.show()

    return xycoord
