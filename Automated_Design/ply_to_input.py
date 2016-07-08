from os import path
import numpy as np

def ply_to_input(fname_no_ply, min_len_nt):
    """
    Converts PLY file into design variables for DX_cage_design input
    Inputs: fname_no_ply = string containing name of PLY file without '.ply'
            min_len_nt = the number of nucleotides long the smallest edge
               will have. Each edge must be a multiple of 10.5 bp, min 31 bp.
    Outputs: coordinates = Vx3 matrix of spatial coordinates of vertices,
               V = number of vertices
             edges = Ex2 matrix where each row corresponds to one edge,
               denoting the vertices being connected. 1st column > 2nd column
             faces = Fx2 cell matrix, where F is the number of faces.
               The first column details how many vertices the face has
               The second column details the vertex IDs of the face
             edge_length_vec = column vector of edge lengths
             file_name = string to name structure
             staple_name = string to name staples to order (can be same as
               file_name if length is not an issue)
             singleXOs = 1 if single crossover vertex staples should be used,
               0 if double crossover vertex staples should be used.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https://opensource.org/licenses/GPL-2.0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    """


    file_name = fname_no_ply + '_' + str(min_len_nt)  #TODO: rename to disambiguate from fname?

    fname = fname_no_ply + '.ply'
    assert path.isfile(fname)
    fid = open(fname)

    def extract_number_from_keyword_in_ply_file(filestream, keyword):
        # TODO: make regular expression!!  And funcitonize it?
        # TODO: verify this filestream doesn't restart
        keyword_length = len(keyword)
        temp = ''
        while temp[:keyword_length] != keyword:
            temp = fid.readline()

        full_line_with_keyword = temp
        number = int(full_line_with_keyword[keyword_length:])

        return number

    num_vert = extract_number_from_keyword_in_ply_file(fid, 'element vertex ')
    num_faces = extract_number_from_keyword_in_ply_file(fid, 'element face ')


    coordinates = np.zeros(dtype=np.float64, shape=(num_vert, 3))  #double

    #`int`s:
    # edges
    # faces  # (he used cell structure)

    [edges, faces, edge_length_vec, file_name, staple_name, singleXOs] = [None, None, None, None, None, None]

    return [coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs]
