from StringIO import StringIO

from mock import patch

from Automated_Design.ply_to_input import ply_to_input


def open_string_as_file(string):
    file = StringIO()
    file.write(string)
    file.seek(0)
    return file


# Sample files as strings:
ply_file_01_tetrahedron = """
ply
format ascii 1.0
element vertex 4
property float32 x
property float32 y
property float32 z
element face 4
property list uint8 int32 vertex_indices
end_header
0.000000 0.000000 0.612372
-0.288675 -0.500000 -0.204124
-0.288675 0.500000 -0.204124
0.577350 0.000000 -0.204124
3 0 2 1
3 0 1 3
3 0 3 2
3 1 2 3
"""

ply_file_05_icosahedron = """ply
format ascii 1.0
element vertex 12
property float32 x
property float32 y
property float32 z
element face 20
property list uint8 int32 vertex_indices
end_header
0.000000 0.000000 1.176000
1.051000 0.000000 0.526000
0.324000 1.000000 0.525000
-0.851000 0.618000 0.526000
-0.851000 -0.618000 0.526000
0.325000 -1.000000 0.526000
0.851000 0.618000 -0.526000
0.851000 -0.618000 -0.526000
-0.325000 1.000000 -0.526000
-1.051000 0.000000 -0.526000
-0.325000 -1.000000 -0.526000
0.000000 0.000000 -1.176000
3 0 1 2
3 0 2 3
3 0 3 4
3 0 4 5
3 0 5 1
3 1 5 7
3 1 7 6
3 1 6 2
3 2 6 8
3 2 8 3
3 3 8 9
3 3 9 4
3 4 9 10
3 4 10 5
3 5 10 7
3 6 7 11
3 6 11 8
3 7 10 11
3 8 11 9
3 9 11 10
"""


@patch('Automated_Design.ply_to_input'
       '.extract_file_reader_and_shape_name_from_input_filename')
def test_ply_input_for_01_tetrahedron(extract_file_and_reader_mock):
    f = open_string_as_file(ply_file_01_tetrahedron)
    shape_name = '01_tetrahedron'
    extract_file_and_reader_mock.return_value = [f, shape_name]

    coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
        singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
    assert len(coordinates) == 4
    assert len(faces) == 4


@patch('Automated_Design.ply_to_input'
       '.extract_file_reader_and_shape_name_from_input_filename')
def test_ply_input_for_05_icosahedron(extract_file_and_reader_mock):
    f = open_string_as_file(ply_file_05_icosahedron)
    shape_name = '05_icosahedron'
    extract_file_and_reader_mock.return_value = [f, shape_name]

    coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
        singleXOs = ply_to_input("05_tetrahedron", '', min_len_nt=52)
    assert len(coordinates) == 12
    assert len(faces) == 20
