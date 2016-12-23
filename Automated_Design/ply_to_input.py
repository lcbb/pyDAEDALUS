from os import path
import numpy as np

from Automated_Design.plotters import plot_edge_length_distributions


def ply_as_filename_to_input(fname_no_ply, min_len_nt=31):
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
    fname = fname_no_ply + '.ply'
    assert path.isfile(fname)
    f = open(fname)
    return ply_to_input(fname_no_ply, f, min_len_nt)


# TODO: clean up the almost-redundant `fname_no_ply` and `f`:
def ply_to_input(fname_no_ply, f, min_len_nt):

    file_name = fname_no_ply + '_' + str(min_len_nt)

    def extract_number_from_keyword_in_ply_file(filestream, keyword):
        # TODO: make regular expression?

        # Eat (read through) all lines up to and including line with `keyword`:
        line = ''
        while keyword not in line:
            line = filestream.readline()

        # grab last whitespace-delimited bit from line and convert to int:
        number = int(line.strip().split()[-1])

        return number

    num_vert = extract_number_from_keyword_in_ply_file(f, 'element vertex')
    num_faces = extract_number_from_keyword_in_ply_file(f, 'element face')

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

        return coordinates_as_list

    coordinates = extract_coordinates_from_file(f, num_vert)

    def extract_faces_from_file(filestream, number_of_faces):

        faces_as_list = []
        for face_id in range(number_of_faces):
            line = filestream.readline()
            line_as_list_of_ints = map(int, line.strip().split())
            number_of_vertices = line_as_list_of_ints[0]
            vertices = line_as_list_of_ints[1:]

            assert number_of_vertices == len(vertices)
            faces_as_list.append(vertices)

        return faces_as_list

    faces = extract_faces_from_file(f, num_faces)

    def remove_unused_vertices(coordinates, faces, number_of_vertices):
        # Determine if you need to clean the vertex indices:
        used_face_ids_as_list = []
        for vertices in faces:
            used_face_ids_as_list += vertices
        unique_used_faces_as_set = set(used_face_ids_as_list)
        # sort to enforce consistency:
        unique_used_faces_as_list = sorted(list(unique_used_faces_as_set))
        cleaning_needed = len(unique_used_faces_as_set) < number_of_vertices

        # Then clean if needed:
        if cleaning_needed:
            # Remove unused row from coordinates data.
            new_coordinates = []
            for i, row in enumerate(coordinates):
                if i in unique_used_faces_as_set:
                    new_coordinates.append(row)
            coordinates = new_coordinates

            new_faces = []
            for current_vertices in faces:
                new_vertices = []
                for vertex in current_vertices:
                    new_vertex = unique_used_faces_as_list.index(vertex)
                    new_vertices.append(new_vertex)
                new_faces.append(new_vertices)
            faces = new_faces

        return coordinates, faces

    coordinates, faces = remove_unused_vertices(coordinates, faces, num_vert)

    def get_edges_from_faces(faces):
        edges = []
        for vertices in faces:
            curr_face = list(vertices)  # force python to make make copy
            curr_face += [curr_face[0]]
            for i in range(len(curr_face)-1):
                if curr_face[i + 1] > curr_face[i]:
                    edges.append((curr_face[i + 1], curr_face[i]))
        return np.array(edges)

    edges = get_edges_from_faces(faces)

    def get_edge_lengths(edges, coordinates):
        edge_length_vec = []
        for edge in edges:
            beginning, end = edge
            length = np.linalg.norm(np.array(coordinates[beginning]) -
                                    np.array(coordinates[end]))
            edge_length_vec.append(length)
        return edge_length_vec

    edge_length_PLY = get_edge_lengths(edges, coordinates)

    def get_scaled_and_rounded_edge_lengths(edge_lengths, min_len_nt):
        min_edge_PLY = min(edge_lengths)
        scale = float(min_len_nt)/min_edge_PLY
        scale_edge_length_PLY = np.rint(scale * np.array(edge_lengths))
        rounded = np.rint(scale_edge_length_PLY / 10.5)

        rounded_edge_length_PLY = []
        for edge_ID in range(len(edge_lengths)):
            rounded_times_ten_point_five = rounded[edge_ID]*10.5
            remainder = rounded[edge_ID] % 2
            if remainder == 0:
                final_length = rounded_times_ten_point_five
            elif scale_edge_length_PLY[edge_ID] > rounded_times_ten_point_five:
                final_length = rounded_times_ten_point_five + 0.5
            else:
                final_length = rounded_times_ten_point_five - 0.5
            rounded_edge_length_PLY.append(int(final_length))
        return scale_edge_length_PLY, rounded_edge_length_PLY

    scale_edge_length_PLY, rounded_edge_length_PLY = \
        get_scaled_and_rounded_edge_lengths(edge_length_PLY, min_len_nt)
    edge_length_vec = rounded_edge_length_PLY

    # Other parameters
    staple_name = file_name  # set short name as file name by default

    if min_len_nt < 42:
        singleXOs = 0
    else:
        singleXOs = 1

    plot_edge_length_distributions(fname_no_ply,
                                   scale_edge_length_PLY,
                                   rounded_edge_length_PLY)

    coordinates = np.array(coordinates)
    # faces = np.array(faces)

    return [coordinates, edges, faces, edge_length_vec, file_name,
            staple_name, singleXOs]
