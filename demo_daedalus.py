from os import path
from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.DX_cage_design import DX_cage_design
fname_no_ply = path.join('PLY_Files', '05_icosahedron_with_unused_vertex') # No '.ply' extension required
min_len_nt = 52 # minimum edge length. Input [] for default minimum 31.
[coordinates, edges, faces_with_count, edge_length_vec, file_name, staple_name, singleXOs] = ply_to_input(fname_no_ply, min_len_nt)
faces = [face for count, face in faces_with_count]

# print("\ncoordinates\n")
# print(coordinates)
#
# print("\nfaces\n")
# print(faces)
#
# print("\nedges\n")
# print(edges)
#
# print("\nedges lengths\n")
# print(edge_length_vec)


scaf_seq = [] # Using default scaffold sequence
scaf_name = [] # Using default scaffold name
full_file_name = DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name)




