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
        # TODO: make regular expression?
        line = ''
        while keyword not in line:
            line = filestream.readline()

        number = int(line.strip().split()[-1]) # grab last whitespace-delimited bit from line and convert to int
        return number

    num_vert = extract_number_from_keyword_in_ply_file(fid, 'element vertex')
    num_faces = extract_number_from_keyword_in_ply_file(fid, 'element face')

    def extract_coordinates_from_file(filestream, number_of_vertices):
        # read lines up to and including the one containing 'end header':
        temp = ''
        while 'end_header' not in temp:  # not strict equality, because
            temp = filestream.readline()

        coordinates_as_list = []
        for i in range(number_of_vertices):
            line = filestream.readline()
            line_as_list = line.split()
            coords_on_this_line = map(float, line_as_list)
            coordinates_as_list.append(coords_on_this_line)

        #TODO: make this into an array?
        # coordinates_as_array = np.array(coordinates_as_list, dtype=np.float64)
        return coordinates_as_list

    coordinates = extract_coordinates_from_file(fid, num_vert)

    def extract_faces_from_file(filestream, number_of_faces):
        # TODO: verify this filestream doesn't restart

        faces_as_list = []
        for face_id in range(number_of_faces):
            line = filestream.readline()
            line_as_list_of_ints = map(int, line.strip().split())
            number_of_vertices = line_as_list_of_ints[0]
            vertices = line_as_list_of_ints[1:]

            assert number_of_vertices == len(vertices)
            faces_as_list.append([number_of_vertices, vertices])

        #TODO: make this into an array?
        return faces_as_list

    faces = extract_faces_from_file(fid, num_faces)

    #`int`s:
    # edges
    # faces  # (he used cell structure)

    [edges, edge_length_vec, file_name, staple_name, singleXOs] = [None, None, None, None, None]

    return [coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs]
