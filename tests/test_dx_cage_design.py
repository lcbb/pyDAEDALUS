from StringIO import StringIO
from unittest import TestCase

from Automated_Design.ply_to_input import ply_to_input
from tests.sample_data import ply_file_01_tetrahedron, ply_file_05_icosahedron
from tests.utils import open_string_as_file


class TestReadingPlyFile(TestCase):
    def test_01_tetrahedron(self):
        f = open_string_as_file(ply_file_01_tetrahedron)
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs = ply_to_input("01_tetrahedron", f, min_len_nt=52)
        self.assertEqual(len(coordinates), 4)
        self.assertEqual(len(faces), 4)
        #TODO: improve these assertions

    def test_05_icosahedron(self):
        f = open_string_as_file(ply_file_05_icosahedron)
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs = ply_to_input("05_tetrahedron", f, min_len_nt=52)
        self.assertEqual(len(coordinates), 12)
        self.assertEqual(len(faces), 20)
        # TODO: improve these assertions

class TestSplitEdge(TestCase):
    def test_ending_node_count(self):
        self.fail()

