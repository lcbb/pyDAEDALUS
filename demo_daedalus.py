from os import path
from Automated_Design.ply_to_input import ply_as_filename_to_input
from Automated_Design.DX_cage_design import DX_cage_design

def run_demo(fname_no_ply, min_len_nt=52, results_foldername=None):
    [coordinates, edges, faces_with_count, edge_length_vec, file_name, staple_name, singleXOs] = ply_as_filename_to_input(fname_no_ply, min_len_nt)
    faces = [face for count, face in faces_with_count]

    scaf_seq = [] # Using default scaffold sequence
    scaf_name = [] # Using default scaffold name
    full_file_name = DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name, results_foldername=results_foldername)



if __name__ == '__main__':
    fname_no_ply = path.join('PLY_Files', '05_icosahedron') # No '.ply' extension required
    fname_no_ply = path.join('PLY_Files', '01_tetrahedron')
    min_len_nt = 52
    run_demo(fname_no_ply, min_len_nt=min_len_nt)