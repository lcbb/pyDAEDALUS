from Automated_Design.ply_to_input import ply_to_input
from tests.sample_data import ply_file_01_tetrahedron, ply_file_05_icosahedron
from tests.utils import open_string_as_file


def test_ply_input_for_01_tetrahedron():
    f = open_string_as_file(ply_file_01_tetrahedron)
    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_to_input(
            "01_tetrahedron", f, min_len_nt=52, results_foldername=False)
    assert len(coordinates) == 4
    assert len(faces) == 4


def test_ply_input_for_05_icosahedron():
    f = open_string_as_file(ply_file_05_icosahedron)
    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_to_input(
            "05_tetrahedron", f, min_len_nt=52, results_foldername=False)
    assert len(coordinates) == 12
    assert len(faces) == 20
