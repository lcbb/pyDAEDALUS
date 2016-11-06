from os import path
import numpy as np


def ply_as_filename_to_input(fname_no_ply, min_len_nt=31, results_foldername=None):
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
    return ply_to_input(fname_no_ply, f, min_len_nt, results_foldername)

#TODO: clean up the almost-redundant `fname_no_ply` and `f` being passed in here.
def ply_to_input(fname_no_ply, f, min_len_nt, results_foldername):

    file_name = fname_no_ply + '_' + str(min_len_nt)

    def extract_number_from_keyword_in_ply_file(filestream, keyword):
        # TODO: make regular expression?
        line = ''
        while keyword not in line:
            line = filestream.readline()

        number = int(line.strip().split()[-1]) # grab last whitespace-delimited bit from line and convert to int
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

        #TODO: make this into an array?
        # coordinates_as_array = np.array(coordinates_as_list, dtype=np.float64)
        return coordinates_as_list

    coordinates = extract_coordinates_from_file(f, num_vert)

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

    faces = extract_faces_from_file(f, num_faces)

    def remove_unused_vertices(coordinates, faces, number_of_vertices):
        # Determine if you need to clean the vertex indices:
        used_face_ids_as_list = []
        for n, vertices in faces:
            used_face_ids_as_list += vertices
        unique_used_faces_as_set = set(used_face_ids_as_list)
        unique_used_faces_as_list = sorted(list(unique_used_faces_as_set))  # sort to enforce consistency
        cleaning_needed = len(unique_used_faces_as_set) < number_of_vertices

        # Then clean if needed:
        if cleaning_needed:
            #TODO: raise warning if cleaning is needed / used?
            # Remove unused row from coordinates data.
            new_coordinates = []
            for i, row in enumerate(coordinates):
                if i in unique_used_faces_as_set:
                    new_coordinates.append(row)
            coordinates = new_coordinates

            new_faces = []
            for face in faces:
                number_of_vertices = face[0]
                current_vertices = face[1]
                new_vertices = []
                for vertex in current_vertices:
                    new_vertex = unique_used_faces_as_list.index(vertex)
                    new_vertices.append(new_vertex)
                new_faces.append([number_of_vertices, new_vertices])
            faces = new_faces

        return coordinates, faces

    coordinates, faces = remove_unused_vertices(coordinates, faces, num_vert)

    def get_edges_from_faces(faces):
        edges = []
        for number_of_vertices, vertices in faces:
            curr_face = list(vertices)  #force python to make copy rather than create reference
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
            length = np.linalg.norm(np.array(coordinates[beginning]) - np.array(coordinates[end]))
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

    scale_edge_length_PLY, rounded_edge_length_PLY = get_scaled_and_rounded_edge_lengths(edge_length_PLY, min_len_nt)
    edge_length_vec = rounded_edge_length_PLY

    # Other parameters
    staple_name = file_name  # set short name as file name by default

    if min_len_nt < 42:
        singleXOs = 0
    else:
        singleXOs = 1

    def plot_edge_length_distributions(scale_edge_length_PLY, rounded_edge_length_PLY, results_foldername):
        from matplotlib import pyplot as plt
        min_len_nt = min(rounded_edge_length_PLY)  #TODO: This is right, right?  It's a little cleaner to not have to import it if the info is already in a variable we're passing in.
        bins_for_hist = range(min_len_nt, max(rounded_edge_length_PLY) + 3, 1)
        bins_for_plotting = [x - 0.25 for x in bins_for_hist[:-1]]


        fig_31 = plt.figure(31, figsize=(8, 6))
        fig_31.clf()
        y31 = np.histogram(scale_edge_length_PLY, bins=bins_for_hist)[0]
        plt.bar(bins_for_plotting, y31, width=0.5)
        plt.xlim((min_len_nt - 0.5, max(rounded_edge_length_PLY) + 1.5))
        plt.ylim((0, max(y31)+1))
        plt.title('Minimum edge length {} bp'.format(min_len_nt));
        plt.xlabel('Edge length (bp)')
        plt.ylabel('Number of edges')
        plt.xticks(bins_for_hist)


        fig_32 = plt.figure(32, figsize=(8, 6))
        fig_32.clf()
        y32 = np.histogram(rounded_edge_length_PLY, bins=bins_for_hist)[0]
        bins_a = [x - 0.125 for x in bins_for_plotting]
        bins_b = [x + 0.38 for x in bins_for_plotting]
        plt.bar(bins_a, y31, width=0.25, color='b')
        plt.bar(bins_b, y32, width=0.25, color='r')

        plt.xlim((min_len_nt - 0.5, max(rounded_edge_length_PLY) + 1.5))
        plt.ylim((0, max(max(y31)+1, 1)))
        plt.title('Edges rounded to nearest 10.5 bp')
        plt.xlabel('Edge length (bp)')
        plt.ylabel('Number of edges')
        plt.xticks(bins_for_hist)

        if results_foldername:
            shape_name = path.basename(path.normpath(fname_no_ply))
            shape_name_with_len = shape_name + '_{}_'.format(min_len_nt)
            base_filename = path.join(results_foldername, shape_name_with_len)
            fig_31.savefig(base_filename+'min_edge_length_dist.png', bbox_inches='tight')
            fig_32.savefig(base_filename+'edges_rounded_to_10_5.png', bbox_inches='tight')
        else:
            plt.show()

    plot_edge_length_distributions(scale_edge_length_PLY, rounded_edge_length_PLY, results_foldername)



    return [coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs]
