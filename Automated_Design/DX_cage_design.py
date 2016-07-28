import networkx as nx
import numpy as np

from Automated_Design.constants import SCAF_SEQ
from Automated_Design.util import generate_graph
from designate_edge_type import designate_edge_type
from gen_schlegel import gen_schlegel
from gen_vert_to_face import gen_vert_to_face


def DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name):
    """
    Creates scaffold routing and staple placement of a DX-based DNA origami
    nano cage.
    Inputs: coordinates = Vx3 matrix of spatial coordinates of vertices,
              V = number of vertices
            edges = Ex2 matrix where each row corresponds to one edge,
              denoting the vertices being connected. 1st column > 2nd column
            faces = F cell matrix, where F is the number of faces.
              The second column details the vertex IDs of the face                #TODO: polish up these 'faces' docs
            edge_length_vec = column vector of edge lengths
            file_name = string to name structure
            staple_name = string to name staples to order (can be same as
              file_name if length is not an issue)
            singleXOs = 1 if single crossover vertex staples should be used,
              0 if double crossover vertex staples should be used.
            scaf_seq = string containing the sequence of the scaffold. If
              using default, input [].
            scaf_name = string containing the name of the scaffold. If using
              default, input [].
    Outputs: dnaInfo = Matlab file containing all spatial and routing
               information.
             routeInfo = Matlab file containing key variables, for debugging
             seq_..._.txt = text file to visualize each edge's sequences and
               nick/crossover information
             staples_..._.xlsx = Excel spreadsheet containing staple seqs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https://opensource.org/licenses/GPL-2.0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """


    # Determine the minimum length scaffold fragment to use.
    len_scaf_used = 2 * sum(edge_length_vec)  # length of scaffold used

    # Determine the default scaffold sequence to use.
    if not scaf_seq:  # if a scaffold sequence was not provided
        if len_scaf_used <= 7249:  # default to using M13 sequence
            # Set scaffold sequence as full M13 scaffold, 7249 nucleotides
            scaf_seq = SCAF_SEQ  # from NEB
            scaf_name = 'full_M13'  # scaffold name
        else:  # generate a random scaffold sequence
            scaf_seq = ''
            for i in range(2 * sum(edge_length_vec)):
                roll = np.random.uniform()
                if roll < 0.25:
                    chosen_letter = 'a'
                elif roll < 0.50:
                    chosen_letter = 't'
                elif roll < 0.75:
                    chosen_letter = 'g'
                else:
                    chosen_letter = 'c'
                scaf_seq += chosen_letter

    scaf_name = 'randomscaf'  # scaffold name
    print(scaf_name)  # display the scaffold name
    print(len(scaf_seq))  # display the length of the scaffold

    # Count number of vertices and edges
    num_vert = len(coordinates)
    num_edges = len(edges)

    ### 1. Generate sparse matrix of connectivities and vertex-face indexing ###
    # Create sparse matrices of connectivities and edge lengths
    full_graph = generate_graph(num_vert, edges, edge_length_vec)

    # Identify presence of every vertex in every face
    vert_to_face = gen_vert_to_face(num_vert, faces)
    # print('\nvert_to_face\n')
    # print(vert_to_face)


    ## 2. Generate spanning tree ##############################################
    # # Designate edges as type 1 or 2:
    # # Type 1: Non-spanning tree, i.e. 1 scaffold crossover in DX cage
    # # Type 2: Spanning tree edges, i.e. 0 scaffold crossovers in DX cage
    edge_type_mat = designate_edge_type(full_graph)
    print("graph with spanning tree data:")
    for edge in edge_type_mat.edges():
        i = edge[0]
        j = edge[1]
        print(i, j, edge_type_mat[i][j])

    gen_schlegel(edges, coordinates, faces, edge_type_mat, 2)


    #TODO: plot circles at each vertex

    #TODO: label each vertex

    full_file_name = None
    return full_file_name
