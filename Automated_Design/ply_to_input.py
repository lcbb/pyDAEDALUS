from os import path
import numpy as np

def ply_to_input(fname_no_ply, min_len_nt):

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
