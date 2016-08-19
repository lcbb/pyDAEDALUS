from unittest import TestCase
import networkx as nx
from os import path

import numpy as np
import scipy.io as sio

from Automated_Design.adj_scaf_nick_pos import adj_scaf_nick_pos
from Automated_Design.adj_scaf_nick_pos import get_scaf_nick_pos
from Automated_Design.assign_scaf_to_edge import assign_scaf_to_edge
from Automated_Design.assign_staples_wChoices import assign_staples_wChoices
from Automated_Design.enum_scaf_bases_DX import enum_scaf_bases_DX
from Automated_Design.set_routing_direction import set_routing_direction
from Automated_Design.gen_stap_seq import gen_stap_seq
from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.split_vert import split_vert
from Automated_Design.split_edge import split_edge
from Automated_Design.toCanDo import toCanDo
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


class TestIntegrationsUsing01Tetrahedron(TestCase):
    """
        Walk through whole chain of processing using the 01_tetrahedron as the
    example comparing to what the actual output from the matlab version were.

        All contained tests follow the pattern:
            - Import original state.
            - Run function to test, catching actual outputs
            - Import target state.
            - Make assertions comparing target to actual.
    """
    # maxDiff = None

    #TODO: move all these test files into an `01_tetrahedron
    target_1_vert_to_face = load_vert_to_face_from_mat('1_vert_to_face.mat')

    target_2_edge_type_mat = load_graph_from_mat('2_edge_type_mat.mat')

    target_3_edge_type_mat_wHalfs = load_graph_from_mat('3_edge_type_mat_wHalfs.mat')
    target_3_pseudo_vert = load_pseudonodes_from_mat('3_pseudo_vert.mat')

    target_4_edge_type_mat_allNodes = load_graph_from_mat('4_edge_type_mat_allNodes.mat')
    target_4_pseudo_vert = load_pseudonodes_from_mat('4_pseudo_vert.mat')

    #TODO: All the following `load`s probably need further parsing out of raw state
    target_0_edges = load_mat_file('0_edges.mat')
    target_0_singleXOs = load_mat_file('0_singleXOs.mat')

    target_1_edge_length_mat_full = load_mat_file('1_edge_length_mat_full.mat')

    target_5_route_real = load_mat_file('5_route_real.mat')
    target_5_route_vals = load_mat_file('5_route_vals.mat')

    target_6_edge_bgn_vec = load_mat_file('6_edge_bgn_vec.mat')
    target_6_edge_fin_vec = load_mat_file('6_edge_fin_vec.mat')
    target_6_edge_type_vec = load_mat_file('6_edge_type_vec.mat')

    target_7_scaf_to_edge = load_mat_file('7_scaf_to_edge.mat')

    target_8_scaf_nick_pos = load_mat_file('8_scaf_nick_pos.mat')
    target_8_scaf_to_edge = load_mat_file('8_scaf_to_edge.mat')

    target_9_staples = load_mat_file('9_staples.mat')


    # 2
    def test_generate_spanning_tree(self):
        self.fail("Write me.")

    # 3
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

    # 4
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

    # 5
    def test_set_routing_direction(self):
        faces = None  #TODO:  import this `0_faces.mat`
        vert_to_face = self.target_1_vert_to_face
        edge_type_mat_allNodes = self.target_4_edge_type_mat_allNodes
        pseudo_vert = self.target_4_pseudo_vert
        num_vert = len(np.unique(pseudo_vert))

        actual_route_real, actual_route_vals = set_routing_direction(
            edge_type_mat_allNodes, num_vert, pseudo_vert, faces, vert_to_face
        )

        target_route_real = self.target_5_route_real
        target_route_vals = self.target_5_route_vals
        self.fail("Write these assertions...")

    # 6
    def test_enum_scaf_bases_DX(self):
        edge_length_mat_full = self.target_1_edge_length_mat_full
        route_real = self.target_5_route_real
        route_vals = self.target_5_route_vals

        actual_edge_bgn_vec, actual_edge_fin_vec, actual_edge_type_vec = \
            enum_scaf_bases_DX(route_real, route_vals, edge_length_mat_full)

        target_edge_bgn_vec = self.target_6_edge_bgn_vec
        target_edge_fin_vec = self.target_6_edge_fin_vec
        target_edge_type_vec = self.target_6_edge_type_vec
        self.fail("Write these assertions...")

    # 7
    def test_assign_scaf_to_edge(self):
        edges = self.target_0_edges
        num_edges = len(edges)
        edge_length_mat_full = self.target_1_edge_length_mat_full
        edge_bgn_vec = self.target_6_edge_bgn_vec
        edge_fin_vec = self.target_6_edge_fin_vec
        edge_type_vec = self.target_6_edge_type_vec

        actual_scaf_to_edge = assign_scaf_to_edge(
            edges, num_edges, edge_length_mat_full,
            edge_bgn_vec, edge_fin_vec, edge_type_vec)

        target_edge_type_vec = self.target_6_edge_type_vec
        target_scaf_to_edge = self.target_7_scaf_to_edge
        self.fail("Write these assertions...")

    # 8
    def test_get_scaf_nick_pos(self):
        edges = self.target_0_edges
        route_real = self.target_5_route_real

        actual_scaf_nick_pos = get_scaf_nick_pos(edges, route_real)

        target_scaf_nick_pos = self.target_8_scaf_nick_pos
        self.fail("Write these assertions...")

    def test_adj_scaf_nick_pos(self):
        num_bases = len(self.target_6_edge_type_vec)
        scaf_to_edge = self.target_7_scaf_to_edge
        scaf_nick_pos = self.target_8_scaf_nick_pos

        scaf_to_edge_adj = adj_scaf_nick_pos(scaf_to_edge, scaf_nick_pos, num_bases)

        target_scaf_to_edge = self.target_8_scaf_to_edge
        self.fail("Write these assertions...")

    # 9
    def test_assign_staples_wChoices(self):
        singleXOs = self.target_0_singleXOs
        edges = self.target_0_edges
        num_edges = len(edges)
        edge_type_mat = self.target_2_edge_type_mat
        scaf_to_edge = self.target_8_scaf_to_edge
        num_bases = len(self.target_6_edge_type_vec)
        num_vert = edge_type_mat.size()

        staples = assign_staples_wChoices(
            edges, num_edges, edge_type_mat, scaf_to_edge,
            num_bases, num_vert, singleXOs)

        target_staples = self.target_9_staples
        self.fail("Write these assertions...")

    # 10
    def test_gen_stap_seq(self):
        staples, num_edges, scaf_seq, staple_name, scaf_name, len_scaf_used = [None] * 6
        stuff = gen_stap_seq(staples, num_edges, scaf_seq, staple_name, scaf_name, len_scaf_used)
        self.assertEqual(stuff, "TODO: bring in real targets.  Write Assertions")

    # 11
    def test_toCanDo(self):
        scaf_to_edge, scaf_seq, stap_list, stap_seq_list, coordinates, edges, edge_length_vec, faces, vert_to_face, fig = [None] * 10
        stuff = toCanDo(scaf_to_edge, scaf_seq, stap_list, stap_seq_list, coordinates,
                        edges, edge_length_vec, faces, vert_to_face, fig)

        self.assertEqual(stuff, "TODO: Bring in real targets.  How to properly test this function?")