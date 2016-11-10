import networkx as nx
import numpy as np
from os import path

from Automated_Design.adj_scaf_nick_pos import adj_scaf_nick_pos, \
    get_scaf_nick_pos
from Automated_Design.assign_scaf_to_edge import assign_scaf_to_edge
from Automated_Design.assign_staples_wChoices import assign_staples_wChoices
from Automated_Design.constants import SCAF_SEQ
from Automated_Design.dna_info import DnaInfo
from Automated_Design.enum_scaf_bases_DX import enum_scaf_bases_DX
from Automated_Design.gen_stap_seq import gen_stap_seq
from Automated_Design.set_routing_direction import set_routing_direction
from Automated_Design.split_edge import split_edge
from Automated_Design.split_vert import split_vert
from Automated_Design.util import generate_graph
from designate_edge_type import designate_edge_type
from gen_schlegel import gen_schlegel
from gen_vert_to_face import gen_vert_to_face


def DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name, results_foldername=None):
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
    scaf_seq_default_used = False
    if not scaf_seq:  # if a scaffold sequence was not provided
        if len_scaf_used <= 7249:  # default to using M13 sequence
            scaf_seq_default_used = True
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

    # Count number of vertices and edges
    num_vert = len(coordinates)
    num_edges = len(edges)

    ### 1. Generate sparse matrix of connectivities and vertex-face indexing ###
    # Create sparse matrices of connectivities and edge lengths
    file_name_without_containing_folder = path.split(file_name)[1]
    full_graph = generate_graph(num_vert, edges, edge_length_vec)
    full_graph.name = file_name_without_containing_folder
    nx.write_gml(full_graph,
                 'tests/sample_files/' + full_graph.name + '_initial_graph.gml')

    # Identify presence of every vertex in every face
    vert_to_face = gen_vert_to_face(num_vert, faces)


    ## 2. Generate spanning tree ##############################################
    # # Designate edges as type 1 or 2:
    # # Type 1: Non-spanning tree, i.e. 1 scaffold crossover in DX cage
    # # Type 2: Spanning tree edges, i.e. 0 scaffold crossovers in DX cage
    edge_type_mat = designate_edge_type(full_graph)
    graph_with_spanning_tree_marked = edge_type_mat   #TODO: this rename
    if results_foldername:
        schlegel_filename = file_name_without_containing_folder + '_schlegel.png'
        full_schlegel_filename = path.join(results_foldername, schlegel_filename)
    else:
        full_schlegel_filename = None
    gen_schlegel(edges, coordinates, faces, edge_type_mat=edge_type_mat, schlegel_filename=full_schlegel_filename)
    edge_type_mat = edge_type_mat.to_directed()  # MST in networkx requires an undirected graph?  Later code requires directed?

    ## 3. Add nodes to edges ##################################################
    # Add two nodes to each nontree edge to implement scaffold crossovers
    edge_type_mat_wHalfs, pseudo_vert = split_edge(edge_type_mat, num_vert)
    graph_with_edges_split = edge_type_mat_wHalfs  #TODO: this rename

    ## 4. Add nodes to vertices ###############################################
    # # Split each vertex into N nodes, where N is degree of vertex
    edge_type_mat_allNodes, pseudo_vert = split_vert(
        edge_type_mat_wHalfs, pseudo_vert, num_vert, vert_to_face)

    ## 5. Set direction of routing ############################################
    [route_real, route_vals] = set_routing_direction(
        edge_type_mat_allNodes, num_vert, pseudo_vert, faces, vert_to_face)

    ## 6. Enumerate scaffold bases ############################################
    edge_length_mat_full = full_graph  #TODO: Did I save the ege lengths onto this one, too?  If not, need to propogate edge lengths to this point
    [edge_bgn_vec, edge_fin_vec, edge_type_vec] = enum_scaf_bases_DX(
        route_real, route_vals, edge_length_mat_full)

    num_bases = len(edge_type_vec)

    ## 7. Assign enumerated scaffold bases to edges ###########################
    scaf_to_edge = assign_scaf_to_edge(edges, num_edges, edge_type_mat,
                                       edge_bgn_vec, edge_fin_vec,
                                       edge_type_vec)

    ## 8. Adjust scaffold nick position #######################################
    scaf_nick_pos = get_scaf_nick_pos(edges, route_real, edge_length_vec)
    scaf_to_edge_adj = adj_scaf_nick_pos(scaf_to_edge, scaf_nick_pos,
                                         num_bases)
    scaf_to_edge = scaf_to_edge_adj

    ## 9. Add staples #########################################################
    staples = assign_staples_wChoices(edges, num_edges, edge_type_mat,
                                      scaf_to_edge, num_bases, num_vert,
                                      singleXOs)

    ## 10. Assign sequence to staples #########################################
    if scaf_seq:  # if a scaffold sequence has been input
        # TODO: won't this always be true, since even if scaf_seq started as
        # `[]`, it would have still been defined above?
        [stap_seq, stap_seq_list, stap_list,
         named_stap_seq_list] = gen_stap_seq(staples, num_edges, scaf_seq,
                                             staple_name, scaf_name,
                                             len_scaf_used)


    full_file_name = None
    return full_file_name
