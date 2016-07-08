from os import path
from Automated_Design.ply_to_input import ply_to_input

fname_no_ply = path.join('PLY_Files', '05_icosahedron_with_unused_vertex') # No '.ply' extension required
min_len_nt = 52 # minimum edge length. Input [] for default minimum 31.
[coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs] = ply_to_input(fname_no_ply, min_len_nt)

print("\ncoordinates\n")
print(coordinates)
print("\nfaces\n")
print(faces)

