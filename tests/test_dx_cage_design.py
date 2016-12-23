from StringIO import StringIO
from unittest import TestCase, expectedFailure

import networkx as nx
import numpy as np
from mock.mock import patch
from networkx.algorithms.isomorphism.isomorph import is_isomorphic

from Automated_Design.adj_scaf_nick_pos import adj_scaf_nick_pos
from Automated_Design.adj_scaf_nick_pos import get_scaf_nick_pos
from Automated_Design.assign_scaf_to_edge import assign_scaf_to_edge
from Automated_Design.assign_staples_wChoices import assign_staples_wChoices
from Automated_Design.designate_edge_type import designate_edge_type
from Automated_Design.enum_scaf_bases_DX import enum_scaf_bases_DX
from Automated_Design.gen_stap_seq import gen_stap_seq
from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.set_routing_direction import set_routing_direction
from Automated_Design.split_edge import split_edge
from Automated_Design.split_vert import split_vert
from Automated_Design.dna_info import DnaInfo, calc_buff, d, wDX, gen_FE_norms
from tests.sample_data import ply_file_01_tetrahedron, ply_file_05_icosahedron
from tests.target_parsers import load_mat_file, load_graph_from_mat, \
    load_pseudonodes_from_mat, load_1d_list_from_mat, load_single_value, \
    load_edges_from_mat, load_faces_from_mat, load_vert_to_face_from_mat, \
    load_scaf_to_edge_from_mat, load_staples_from_mat, load_stap_seq_file, \
    load_stap_seq_list_file, load_stap_list_file, \
    load_named_stap_seq_list_file, load_dna_info
from tests.utils import open_string_as_file


@patch('Automated_Design.plotters.plot_edge_length_distributions')
def test_ply_input_for_01_tetrahedron(self, savefig_mock):
    f = open_string_as_file(ply_file_01_tetrahedron)
    coordinates, edges, faces, edge_length_vec, file_name, \
    staple_name, singleXOs = ply_to_input(
        "01_tetrahedron", f, min_len_nt=52)
    self.assertEqual(len(coordinates), 4)
    self.assertEqual(len(faces), 4)
    # TODO: mock fig itself and assert about what's `fig_mock.plot`ed


@patch('Automated_Design.plotters.plot_edge_length_distributions')
def test_ply_input_for_05_icosahedron(self, savefig_mock):
    f = open_string_as_file(ply_file_05_icosahedron)
    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_to_input(
            "05_tetrahedron", f, min_len_nt=52)
    self.assertEqual(len(coordinates), 12)
    self.assertEqual(len(faces), 20)
    # TODO: mock fig itself and assert about what's `fig_mock.plot`ed


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

    target_0_edges = load_edges_from_mat('0_edges.mat')
    target_0_faces = load_faces_from_mat('0_faces.mat')
    target_0_singleXOs = load_single_value('0_singleXOs.mat')
    target_0_edge_length_vec = load_1d_list_from_mat('0_edge_length_vec.mat')
    target_0_scaf_seq = load_mat_file("0_scaf_seq.mat")[0]

    target_0_staple_name = load_mat_file("0_staple_name.mat")[0]
    target_0_scaf_name = load_mat_file("0_scaf_name.mat")[0]
    target_0_coordinates = load_mat_file("0_coordinates.mat")  # TODO: Format

    target_1_vert_to_face = load_vert_to_face_from_mat('1_vert_to_face.mat')
    target_1_edge_length_mat_full = load_graph_from_mat(
        '1_edge_length_mat_full.mat',
        edge_attribute='length',
        graph_type=nx.Graph())
    target_1_full_graph = load_graph_from_mat(
        '1_full_graph.mat', graph_type=nx.Graph())

    target_2_edge_type_mat = load_graph_from_mat('2_edge_type_mat.mat',
                                                 graph_type=nx.Graph())
    target_3_edge_type_mat_wHalfs = load_graph_from_mat(
        '3_edge_type_mat_wHalfs.mat')
    target_3_pseudo_vert = load_pseudonodes_from_mat('3_pseudo_vert.mat')

    target_4_edge_type_mat_allNodes = load_graph_from_mat(
        '4_edge_type_mat_allNodes.mat')
    target_4_pseudo_vert = load_pseudonodes_from_mat('4_pseudo_vert.mat')

    target_5_route_real = load_1d_list_from_mat('5_route_real.mat',
                                                are_1_indexed_nodes=True)
    target_5_route_vals = load_1d_list_from_mat('5_route_vals.mat')

    target_6_edge_bgn_vec = list(
        load_mat_file('6_edge_bgn_vec.mat').flatten() - 1)
    target_6_edge_fin_vec = list(
        load_mat_file('6_edge_fin_vec.mat').flatten() - 1)
    target_6_edge_type_vec = list(
        load_mat_file('6_edge_type_vec.mat').flatten())

    target_7_scaf_to_edge = load_scaf_to_edge_from_mat('7_scaf_to_edge.mat')

    target_8_scaf_to_edge = load_scaf_to_edge_from_mat('8_scaf_to_edge.mat')
    target_8_scaf_nick_pos = load_single_value('8_scaf_nick_pos.mat')

    target_9_staples = load_staples_from_mat('9_staples.mat')

    target_10_stap_seq = load_stap_seq_file("10_stap_seq.mat")
    target_10_stap_seq_list = load_stap_seq_list_file("10_stap_seq_list.mat")
    target_10_stap_list = load_stap_list_file("10_stap_list.mat")
    target_10_named_stap_seq_list = load_named_stap_seq_list_file(
        "10_named_stap_seq_list.mat")

    target_11_buff_nt = load_1d_list_from_mat('11_buff_nt.mat')
    target_11_edge_norms = load_mat_file('11_edge_norms.mat')  # TODO: Format
    target_11_dna_info = load_dna_info("11_dna_info.mat")

    # 1:  I'm using a networkx.Graph rather than sparse matrix.  No direct
    # assertion possible, though we still could do assertions on the
    # properties.

    # 2
    def test_generate_spanning_tree(self):
        full_graph = self.target_1_full_graph

        actual_edge_type_mat = designate_edge_type(full_graph)

        target_edge_type_mat = self.target_2_edge_type_mat
        self.assertEqual(actual_edge_type_mat.nodes(),
                         target_edge_type_mat.nodes())
        self.assertEqual(len(actual_edge_type_mat.edges()),
                         len(target_edge_type_mat.edges()))
        self.assertEqual(actual_edge_type_mat.edges(data=True),
                         target_edge_type_mat.edges(data=True))
        # TODO: assert graphs are isomorphic instead of directly equal?

    # 3
    def test_split_edge(self):
        edge_type_mat = self.target_2_edge_type_mat.to_directed()
        # TODO:  Figure out where `.to_directed` transition needs to happen!
        num_vert = len(edge_type_mat.nodes())

        actual_edge_type_mat_wHalfs, actual_pseudo_vert = split_edge(
            edge_type_mat, num_vert)

        target_edge_type_mat_wHalfs = self.target_3_edge_type_mat_wHalfs
        target_pseudo_vert = self.target_3_pseudo_vert
        self.assertEqual(actual_pseudo_vert, target_pseudo_vert)
        self.assertEqual(actual_edge_type_mat_wHalfs.nodes(),
                         target_edge_type_mat_wHalfs.nodes())
        self.assertEqual(set(actual_edge_type_mat_wHalfs.edges()),
                         set(target_edge_type_mat_wHalfs.edges()))

    # 4
    def test_split_vert(self):
        edge_type_mat_wHalfs = self.target_3_edge_type_mat_wHalfs
        pseudo_vert = self.target_3_pseudo_vert
        num_vert = len(np.unique(pseudo_vert))
        vert_to_face = self.target_1_vert_to_face

        actual_edge_type_mat_allNodes, actual_pseudo_vert = split_vert(
            edge_type_mat_wHalfs, pseudo_vert, num_vert, vert_to_face)

        target_edge_type_mat_allNodes = self.target_4_edge_type_mat_allNodes
        target_pseudo_vert = self.target_4_pseudo_vert
        self.assertEqual(actual_pseudo_vert, target_pseudo_vert)
        self.assertTrue(is_isomorphic(actual_edge_type_mat_allNodes,
                                      target_edge_type_mat_allNodes))

        # TODO: Graphs are isomorphic but not directly equal. ... okay?

    # 5
    def test_set_routing_direction(self):
        faces = self.target_0_faces
        vert_to_face = self.target_1_vert_to_face
        edge_type_mat_allNodes = self.target_4_edge_type_mat_allNodes
        pseudo_vert = self.target_4_pseudo_vert
        num_vert = len(np.unique(pseudo_vert))

        actual_route_real, actual_route_vals = set_routing_direction(
            edge_type_mat_allNodes, num_vert, pseudo_vert, faces, vert_to_face
        )

        target_route_real = self.target_5_route_real
        target_route_vals = self.target_5_route_vals
        self.assertEqual(actual_route_real, target_route_real)
        self.assertEqual(actual_route_vals, target_route_vals)

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
        self.assertEqual(actual_edge_bgn_vec, target_edge_bgn_vec)
        self.assertEqual(actual_edge_fin_vec, target_edge_fin_vec)
        self.assertEqual(actual_edge_type_vec, target_edge_type_vec)

    # 7
    def test_assign_scaf_to_edge(self):
        edges = self.target_0_edges
        num_edges = len(edges)
        edge_type_mat = self.target_2_edge_type_mat
        edge_bgn_vec = self.target_6_edge_bgn_vec
        edge_fin_vec = self.target_6_edge_fin_vec
        edge_type_vec = self.target_6_edge_type_vec

        actual_scaf_to_edge = assign_scaf_to_edge(
            edges, num_edges, edge_type_mat,
            edge_bgn_vec, edge_fin_vec, edge_type_vec)

        target_scaf_to_edge = self.target_7_scaf_to_edge
        self.assertEqual(actual_scaf_to_edge, target_scaf_to_edge)

    # 8
    def test_get_scaf_nick_pos(self):
        edges = self.target_0_edges
        route_real = self.target_5_route_real
        edge_length_vec = self.target_0_edge_length_vec

        actual_scaf_nick_pos = get_scaf_nick_pos(
            edges, route_real, edge_length_vec)

        target_scaf_nick_pos = self.target_8_scaf_nick_pos
        self.assertEqual(actual_scaf_nick_pos, target_scaf_nick_pos)

    def test_adj_scaf_nick_pos(self):
        num_bases = len(self.target_6_edge_type_vec)
        scaf_to_edge = self.target_7_scaf_to_edge
        scaf_nick_pos = self.target_8_scaf_nick_pos

        actual_scaf_to_edge_adj = adj_scaf_nick_pos(
            scaf_to_edge, scaf_nick_pos, num_bases)

        target_scaf_to_edge = self.target_8_scaf_to_edge
        self.assertEqual(actual_scaf_to_edge_adj, target_scaf_to_edge)

    # 9
    def test_assign_staples_wChoices(self):
        singleXOs = self.target_0_singleXOs
        edges = self.target_0_edges
        num_edges = len(edges)
        edge_type_mat = self.target_2_edge_type_mat.to_directed()
        scaf_to_edge = self.target_8_scaf_to_edge
        num_bases = len(self.target_6_edge_type_vec)
        num_vert = len(edge_type_mat.nodes())

        actual_staples = assign_staples_wChoices(
            edges, num_edges, edge_type_mat, scaf_to_edge,
            num_bases, num_vert, singleXOs)

        target_staples = self.target_9_staples
        self.assertEqual(actual_staples, target_staples)

    # 10
    def test_gen_stap_seq(self):
        staples = self.target_9_staples
        num_edges = len(self.target_0_edges)
        len_scaf_used = 2*sum(self.target_0_edge_length_vec)
        scaf_seq = self.target_0_scaf_seq
        staple_name = self.target_0_staple_name
        scaf_name = self.target_0_scaf_name
        actual_stap_seq, actual_stap_seq_list, \
            actual_stap_list, actual_named_stap_seq_list \
            = gen_stap_seq(staples, num_edges, scaf_seq,
                           staple_name, scaf_name, len_scaf_used)

        target_stap_seq = self.target_10_stap_seq
        target_stap_seq_list = self.target_10_stap_seq_list
        target_stap_list = self.target_10_stap_list
        target_named_stap_seq_list = self.target_10_named_stap_seq_list

        self.assertEqual(actual_stap_seq, target_stap_seq)
        self.assertEqual(actual_stap_seq_list, target_stap_seq_list)
        self.assertEqual(actual_stap_list, target_stap_list)
        self.assertEqual(actual_named_stap_seq_list,
                         target_named_stap_seq_list)

    # 11
    def test_calc_buff(self):
        faces = self.target_0_faces
        num_vert = len(self.target_2_edge_type_mat.nodes())
        coordinates = self.target_0_coordinates

        actual_buff_nt = calc_buff(faces, num_vert, coordinates, d, wDX)
        target_buff_nt = self.target_11_buff_nt

        self.assertEqual(actual_buff_nt, target_buff_nt)

    def test_gen_FE_norms(self):
        faces = self.target_0_faces
        edges = self.target_0_edges
        vert_to_face = self.target_1_vert_to_face
        coordinates = self.target_0_coordinates

        actual_edge_norms = gen_FE_norms(
            coordinates, faces, edges, vert_to_face)
        target_edge_norms = self.target_11_edge_norms

        self.assertTrue(np.allclose(actual_edge_norms, target_edge_norms))

    def test_DnaInfo_init(self):
        scaf_to_edge = self.target_8_scaf_to_edge
        scaf_seq = self.target_0_scaf_seq
        stap_list = self.target_10_stap_list
        stap_seq_list = self.target_10_stap_seq_list
        coordinates = self.target_0_coordinates
        edges = self.target_0_edges
        edge_length_vec = self.target_0_edge_length_vec
        faces = self.target_0_faces
        vert_to_face = self.target_1_vert_to_face

        actual_dnaInfo = DnaInfo(scaf_to_edge, scaf_seq, stap_list,
                                 stap_seq_list, coordinates, edges,
                                 edge_length_vec, faces, vert_to_face)
        target_dnaInfo = self.target_11_dna_info

        # assertions for dnaGenom:
        actual_dNode = actual_dnaInfo.dnaGeom.dNode
        target_dNode = target_dnaInfo['dnaGeom']['dNode']
        self.assertTrue(np.allclose(actual_dNode, target_dNode))

        actual_id_nt = actual_dnaInfo.dnaGeom.id_nt
        target_id_nt = target_dnaInfo['dnaGeom']['id_nt']
        self.assertTrue(np.allclose(actual_id_nt, target_id_nt))

        actual_triad = actual_dnaInfo.dnaGeom.triad
        target_triad = target_dnaInfo['dnaGeom']['triad']
        self.assertTrue(np.allclose(actual_triad, target_triad))

        # assertions for dnaInfo:
        actual_dnaTop = actual_dnaInfo.dnaTop
        target_dnaTop = target_dnaInfo['dnaTop']
        # for row in target_dnaTop:
        #     print row
        self.assertEqual(len(actual_dnaTop), len(target_dnaTop))
        for actual, target in zip(actual_dnaTop, target_dnaTop):
            # TODO: write __eq__ opposite into DnaTop to clean this up?
            self.assertEqual(str(actual), str(target))

    @expectedFailure
    def test_dna_info_print_to_cando(self):
        # scaf_to_edge = self.target_8_scaf_to_edge
        # scaf_seq = self.target_0_scaf_seq
        # stap_list = self.target_10_stap_list
        # stap_seq_list = self.target_10_stap_seq_list
        # coordinates = self.target_0_coordinates
        # edges = self.target_0_edges
        # edge_length_vec = self.target_0_edge_length_vec
        # faces = self.target_0_faces
        # vert_to_face = self.target_1_vert_to_face

        # dnaInfo = DnaInfo(scaf_to_edge, scaf_seq, stap_list,
        #                   stap_seq_list, coordinates, edges,
        #                   edge_length_vec, faces, vert_to_face)
        # dnaInfo.save_dna_info_to_cando_file("THE_CANDO_THING.txt")

        # TODO: mock out file writer to be string writer
        self.fail("Write these assertions!")
        self.fail("And be sure to 1-index the file!")

    @expectedFailure
    def test_seq_to_text(self):
        # scaf_to_edge = self.target_8_scaf_to_edge
        # scaf_seq = self.target_0_scaf_seq
        # stap_list = self.target_10_stap_list
        # stap_seq_list = self.target_10_stap_seq_list
        # coordinates = self.target_0_coordinates
        # edges = self.target_0_edges
        # edge_length_vec = self.target_0_edge_length_vec
        # faces = self.target_0_faces
        # vert_to_face = self.target_1_vert_to_face
        #
        # dnaInfo = DnaInfo(scaf_to_edge, scaf_seq, stap_list,
        #                   stap_seq_list, coordinates, edges,
        #                   edge_length_vec, faces, vert_to_face)
        #
        # seqtoText(scaf_to_edge, edges, dnaInfo, "test_DONTEXIST.txt",
        #           "buttercup", 0, "threebuttercups.txt")

        self.fail("write me!")

    def print_a_thing(self):
        print("Test String To Catch")

    # 13?
    @patch('sys.stdout', new_callable=StringIO)
    def test_printing_with_mocks(self, mock_stdout):
        self.print_a_thing()

        # raise Exception("..{}..".format(mock_stdout.getvalue()))
        assert mock_stdout.getvalue() == 'Test String To Catch\n'
