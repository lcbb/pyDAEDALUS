from unittest import TestCase
import networkx as nx
from os import path

import numpy as np
import scipy.io as sio
from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.split_vert import split_vert
from Automated_Design.split_edge import split_edge
from tests.sample_data import ply_file_01_tetrahedron, ply_file_05_icosahedron
from tests.utils import open_string_as_file


class TestPlyImport(TestCase):
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


def load_mat_file(filename):
    full_data = sio.loadmat(path.join('tests', 'targets', filename))
    # full_data.keys()
    # >>> ['thing', '__version__', '__header__', '__globals__']
    # where `thing` is whatever var was originally saved.
    key = (key for key in full_data.keys() if key[0] is not '_').next()
    data = full_data[key]
    return data
    #TODO: declare all these files currently used as 1_tetra...

def load_graph_from_mat(filename):
    graph_as_sparse_matrix = load_mat_file(filename)
    graph = nx.from_scipy_sparse_matrix(graph_as_sparse_matrix,
                                        create_using=nx.DiGraph(),
                                        edge_attribute='type')
    return graph


def load_pseudonodes_from_mat(filename):
    pseudonodes = load_mat_file(filename)
    pseudonodes = pseudonodes - 1
    return list(pseudonodes.flatten())  #convert shape(N,1) to shape(N)


def load_vert_to_face_from_mat(filename):
    """
        The raw read starts of as a, (N,1)-size length array of each face where
    a given face is a (1,m)-size array inside of a (1,)-shape array.
        I have it return as a list of lists.
    """
    vert_to_face = load_mat_file(filename) - 1
    formatted_vert_to_face = [list(face[0].flatten()) for face in vert_to_face]
    return formatted_vert_to_face


class TestRouting(TestCase):
    # maxDiff = None

    target_1_vert_to_face = load_vert_to_face_from_mat('1_vert_to_face.mat')

    target_2_edge_type_mat = load_graph_from_mat('2_edge_type_mat.mat')
    target_3_edge_type_mat_wHalfs = load_graph_from_mat('3_edge_type_mat_wHalfs.mat')
    target_3_pseudo_vert = load_pseudonodes_from_mat('3_pseudo_vert.mat')
    target_4_edge_type_mat_allNodes = load_graph_from_mat('4_edge_type_mat_allNodes.mat')
    target_4_pseudo_vert = load_pseudonodes_from_mat('4_pseudo_vert.mat')

    def test_split_edge(self):
        edge_type_mat = self.target_2_edge_type_mat
        num_vert = edge_type_mat.size()

        actual_edge_type_mat_wHalfs, actual_pseudo_vert = split_edge(edge_type_mat, num_vert)

        target_edge_type_mat_wHalfs = self.target_3_edge_type_mat_wHalfs
        target_pseudo_vert = self.target_3_pseudo_vert
        self.assertEqual(actual_pseudo_vert, target_pseudo_vert)
        self.assertEqual(actual_edge_type_mat_wHalfs.nodes(),
                         target_edge_type_mat_wHalfs.nodes())
        self.assertEqual(actual_edge_type_mat_wHalfs.edges(),
                         target_edge_type_mat_wHalfs.edges())

    def test_split_vert(self):
        edge_type_mat_wHalfs = self.target_3_edge_type_mat_wHalfs
        pseudo_vert = self.target_3_pseudo_vert
        num_vert = len(np.unique(pseudo_vert))
        vert_to_face = self.target_1_vert_to_face
        actual_edge_type_mat_allNodes, actual_pseudo_vert = split_vert(edge_type_mat_wHalfs, pseudo_vert, num_vert, vert_to_face)

        target_edge_type_allNodes = self.target_4_edge_type_mat_allNodes
        target_pseudo_vert = self.target_4_pseudo_vert
        self.assertEqual(actual_pseudo_vert, target_pseudo_vert)
        self.assertEqual(actual_edge_type_mat_allNodes.nodes(),
                         target_edge_type_allNodes.nodes())
        self.assertEqual(actual_edge_type_mat_allNodes.edges(),
                         target_edge_type_allNodes.edges())


