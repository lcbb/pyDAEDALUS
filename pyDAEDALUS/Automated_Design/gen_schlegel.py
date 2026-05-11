import numpy as np
from matplotlib import pyplot as plt

from Automated_Design.constants import VERMILLION, REDPURPLE, SKYBLUE, WHITE


def create_2d_mapping(edges, coordinates, faces):
    # Choose Biggest Face:
    big_face = max(faces, key=len)

    # Initialize:
    num_face_vert = len(big_face)
    num_vert = len(coordinates)

    # Initialize xycoord
    xycoord = np.zeros(shape=(num_vert, 2))
    num_repeat = 10 * num_vert  # number of iterations to calculate xycoord

    # Set big face on unit circle
    angle = 2. * np.pi / num_face_vert

    for i in range(num_face_vert):  # for each vertex in big face
        xycoord[big_face[i]] = [np.cos(angle * i), np.sin(angle * i)]

    # Calculate positions of other vertices through iterative process
    for repeat_ID in range(num_repeat):

        for vert_ID in range(num_vert):  # For each vertex,
            if vert_ID not in big_face:  # that's not in big_face,
                # move the xycoord of the given vertex to the centerpont
                # of all neighboring nodes.

                # # Find vert_ID
                find_vert_row, find_vert_col = np.where(edges == vert_ID)

                # # Find neighbors
                neighbors = []
                for row, col in zip(find_vert_row, find_vert_col):
                    find_vert_col_nbr = 1 - col
                    neighbors.append(edges[row, find_vert_col_nbr])
                neighbors = np.array(neighbors, dtype=int)

                # # Get neighbor coords
                nbr_coord = xycoord[neighbors, :]

                # # # Find midpoint
                num_nbr = len(find_vert_col)
                midpt = np.array(
                    [sum(nbr_coord[:, 0]), sum(nbr_coord[:, 1])]
                ) / num_nbr

                # # Store as new xycoord for vert_ID
                xycoord[vert_ID, :] = midpt
    return xycoord


def plot_schlegel(edges, edge_type_graph, xycoord,
                  schlegel_filename):  # pragma: no cover
    f = plt.figure(2, figsize=(8, 8))
    f.clf()
    plt.xlim((-1.2, 1.2))
    plt.ylim((-1.2, 1.2))
    plt.axis('off')

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

        # Only check the edge type if `edge_type_graph` is not None.  If it is
        # None, default to coloring edges the same as if all edges weren't a
        # part of the spanning tree.
        if edge_type_graph is not None \
                and (edge_type_graph.edges[edge_bgn,edge_fin]['type'] == 1):
            plt.plot(x_vec, y_vec, linestyle='-', color=SKYBLUE, linewidth=5)
        else:  # if this is a non-tree edge:
            plt.plot(x_vec, y_vec, linestyle='-', color=REDPURPLE, linewidth=8)

    # Plot nodes:
    plt.plot(xycoord[:, 0], xycoord[:, 1], 'o',
             color=VERMILLION, markeredgecolor=VERMILLION, markersize=15)

    # Plot node labels:
    for i, node_xy in enumerate(xycoord):
        x, y = node_xy
        plt.text(x, y, str(i), color=WHITE, fontname='serif',
                 horizontalalignment='center', verticalalignment='center')

    plt.savefig(schlegel_filename, bbox_inches='tight')


def gen_schlegel(edges, coordinates, faces, schlegel_filename,
                 edge_type_graph=None):    # pragma: no cover
    """
    Plots vertices and edges in Schlegel diagram, tracing scaffold path in
    colored lines.

    The three dimensional structure is mapped onto a 2d surface.  This new
    surface is visualized, saved, and optionally displayed to screen.  Edges
    that are a part of the previously calculated spanning tree (that is, the
    scaffold path) as marked in `edge_type_graph` are highlighted a different
    color.

    Parameters
    ----------
    edges : numpy.ndarray
        Ex2 matrix where each row corresponds to one edge, denoting the
        vertices being connected. 1st column > 2nd column; E = number of edges
    coordinates : numpy.ndarray
        Vx3 matrix of spatial coordinates of vertices; V = number of vertices.
    faces : list
        Fx2 cell matrix, where F is the number of faces.  The first column
        details how many vertices the face has.  The second column details the
        vertex IDs of the face.
    schlegel_filename : str
        Filename pointing at where to save the generated plot.
    edge_type_graph : networkx.classes.graph.Graph
        A graph the links have the 'type' property.  That property can have one
        of two values:
            - 1 is non-spanning tree edge: DX edge with 1 scaffold crossover
            - 2 is spanning tree edge: DX edge with 0 scaffold crossovers

    Returns
    -------
    xycoords
        The 2d mapping used to plot the figure.
    """

    xycoord = create_2d_mapping(edges, coordinates, faces)

    plot_schlegel(edges, edge_type_graph, xycoord, schlegel_filename)

    return xycoord
