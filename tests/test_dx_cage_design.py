from unittest import TestCase, expectedFailure
import networkx as nx
from os import path

import numpy as np
import scipy.io as sio

from Automated_Design.adj_scaf_nick_pos import adj_scaf_nick_pos
from Automated_Design.adj_scaf_nick_pos import get_scaf_nick_pos
from Automated_Design.assign_scaf_to_edge import assign_scaf_to_edge
from Automated_Design.assign_staples_wChoices import assign_staples_wChoices
from Automated_Design.designate_edge_type import designate_edge_type
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
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs = ply_to_input("01_tetrahedron", f, min_len_nt=52, results_foldername=None)
        self.assertEqual(len(coordinates), 4)
        self.assertEqual(len(faces), 4)
        #TODO: improve these assertions

    def test_05_icosahedron(self):
        f = open_string_as_file(ply_file_05_icosahedron)
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs = ply_to_input("05_tetrahedron", f, min_len_nt=52, results_foldername=None)
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

def load_graph_from_mat(filename, edge_attribute='type', graph_type=nx.DiGraph()):
    graph_as_sparse_matrix = load_mat_file(filename)
    graph = nx.from_scipy_sparse_matrix(graph_as_sparse_matrix,
                                        create_using=graph_type,
                                        edge_attribute=edge_attribute)
    #intify 'type' if that's the property used (defaults to float)
    if edge_attribute == 'type':
        for i, j, attribute in graph.edges(data=True):
            graph[i][j]['type'] = int(attribute['type'])
    return graph.copy()


def load_pseudonodes_from_mat(filename):
    pseudonodes = load_mat_file(filename)
    pseudonodes = pseudonodes - 1
    return list(pseudonodes.flatten())  #convert shape(N,1) to shape(N)


def load_1d_list_from_mat(filename, are_1_indexed_nodes = False):
    data = load_mat_file(filename)
    if are_1_indexed_nodes:  # convert to 0 indexed nodes
        data = data - 1
    return list(data.flatten())

def load_single_value(filename):
    data = load_mat_file(filename)
    val = data[0][0]
    return val

def load_edges_from_mat(filename):
    data = load_mat_file(filename)
    edges = []
    for row in data:
        zero_indexed_row = [x-1 for x in row]
        edges.append(zero_indexed_row)
    return edges

def load_faces_from_mat(filename):
    data = load_mat_file(filename)
    data = data - 1
    faces = []
    for row in data:
        n, node_data = row
        faces.append(list(node_data.flatten()))
    return faces


def load_vert_to_face_from_mat(filename):
    """
        The raw read starts of as a, (N,1)-size length array of each face where
    a given face is a (1,m)-size array inside of a (1,)-shape array.
        I have it return as a list of lists.
    """
    vert_to_face = load_mat_file(filename) - 1
    formatted_vert_to_face = [list(face[0].flatten()) for face in vert_to_face]
    return formatted_vert_to_face


def load_scaf_to_edge_from_mat(filename):
    """
    target structure:
    final_data = [
                   [ [0,1,2,...], [56, 57, 58, ...] ],
                   [ ..., ... ]
                   ...
                 ]
    """
    data = load_mat_file(filename)
    data = data - 1  # convert to 0 index

    final_data = []
    for a, b in data:
        a = a.flatten()  # each cell is shape=(N,1).  Needs to be shape=(N)
        b = b.flatten()
        final_data.append([list(a), list(b)])
    return final_data


def load_staples_from_mat(filename):
    data = load_mat_file(filename)
    data = np.delete(data, 4, axis=1)
    data = np.delete(data, 1, axis=1)

    listized_data = []
    for row in data:
        cleaned_row = []
        for subrow in row:
            if subrow.size > 0:
                ids = [x-1 if x > 0 else None for x in subrow[0]]
                cleaned_row.append(ids)
            else:
                cleaned_row.append([])
        listized_data.append(cleaned_row)
    return listized_data


def load_stap_seq_file(filename):
    data = load_mat_file(filename)
    data = np.delete(data, 4, axis=1)
    data = np.delete(data, 1, axis=1)

    stap_seq = []
    for row in data:
        cleaned_row = []
        for cell in row:
            if cell.size > 0:
                cleaned_row.append(cell[0])
            else:
                cleaned_row.append(u'')
        stap_seq.append(cleaned_row)
    return stap_seq


def load_stap_seq_list_file(filename):
    data = load_mat_file(filename)
    stap_seq_list = []
    for row in data:
        cell = row[0]
        if cell.size > 0:
            stap_seq_list.append(cell[0])
        else:
            stap_seq_list.append(None)
    return stap_seq_list


def load_stap_list_file(filename):
    data = load_mat_file(filename)
    stap_list = []
    for row in data:
        subrow = row[0]
        if subrow.size > 0:
            # TODO: how to de-dupe work from load_staples_from_mat
            ids = [x - 1 if x > 0 else None for x in subrow[0]]
            stap_list.append(ids)
        else:
            stap_list.append(None)
    return stap_list


def load_named_stap_seq_list_file(filename):

    # These next two function could be so much cleaner with regex...
    def de_1_index_name_a(name):
        parts = name.split('-', 2)
        parts = parts[0].rsplit('_', 1) + parts[1:]
        decrimented_name = '{}_{}-{}-{}'.format(
            parts[0], str(int(parts[1])-1), str(int(parts[2])-1), parts[3])
        return decrimented_name

    def de_1_index_name_b(name):
        parts = name.split('(')
        parts = [parts[0] ] + parts[1].split('-')
        parts = parts[:2] + parts[2].split(')')
        decrimented_name = '{}({}-{}){}'.format(
            parts[0], str(int(parts[1])-1), str(int(parts[2])-1), parts[3])
        return decrimented_name

    data = load_mat_file(filename)
    named_stap_seq_list = []
    for row in data:
        cleaned_row = []
        for cell in row:
            if cell.size > 0:
                cleaned_row.append(cell[0])
            else:
                cleaned_row.append(None)
        named_stap_seq_list.append(cleaned_row)

    #take 1 indexing in names down to 0 indexing
    splitter = named_stap_seq_list.index([None, None])
    for row in named_stap_seq_list[:splitter]:
        row[0] = de_1_index_name_a(row[0])

    for row in named_stap_seq_list[splitter+1:]:
        row[0] = de_1_index_name_b(row[0])

    return named_stap_seq_list


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

    target_1_vert_to_face = load_vert_to_face_from_mat('1_vert_to_face.mat')
    target_1_edge_length_mat_full = load_graph_from_mat(
        '1_edge_length_mat_full.mat',
        edge_attribute='length',
        graph_type=nx.Graph())
    target_1_full_graph = load_graph_from_mat(
        '1_full_graph.mat', graph_type=nx.Graph())

    target_2_edge_type_mat = load_graph_from_mat('2_edge_type_mat.mat',
                                                 graph_type=nx.Graph())
    target_3_edge_type_mat_wHalfs = load_graph_from_mat('3_edge_type_mat_wHalfs.mat')
    target_3_pseudo_vert = load_pseudonodes_from_mat('3_pseudo_vert.mat')

    target_4_edge_type_mat_allNodes = load_graph_from_mat('4_edge_type_mat_allNodes.mat')
    target_4_pseudo_vert = load_pseudonodes_from_mat('4_pseudo_vert.mat')

    target_5_route_real = load_1d_list_from_mat('5_route_real.mat', are_1_indexed_nodes = True)
    target_5_route_vals = load_1d_list_from_mat('5_route_vals.mat')

    target_6_edge_bgn_vec = list(load_mat_file('6_edge_bgn_vec.mat').flatten()-1)
    target_6_edge_fin_vec = list(load_mat_file('6_edge_fin_vec.mat').flatten()-1)
    target_6_edge_type_vec = list(load_mat_file('6_edge_type_vec.mat').flatten())

    target_7_scaf_to_edge = load_scaf_to_edge_from_mat('7_scaf_to_edge.mat')

    target_8_scaf_to_edge = load_scaf_to_edge_from_mat('8_scaf_to_edge.mat')
    target_8_scaf_nick_pos = load_single_value('8_scaf_nick_pos.mat')

    target_9_staples = load_staples_from_mat('9_staples.mat')

    target_10_stap_seq = load_stap_seq_file("10_stap_seq.mat")
    target_10_stap_seq_list = load_stap_seq_list_file("10_stap_seq_list.mat")
    target_10_stap_list = load_stap_list_file("10_stap_list.mat")
    target_10_named_stap_seq_list = load_named_stap_seq_list_file("10_named_stap_seq_list.mat")


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
        #TODO: assert graphs are isomorphic instead of directly equal?

    # 3
    def test_split_edge(self):
        edge_type_mat = self.target_2_edge_type_mat.to_directed()  #TODO:  Figure out where this transition needs to happen, and move it out of this test!!!
        num_vert = len(edge_type_mat.nodes())

        actual_edge_type_mat_wHalfs, actual_pseudo_vert = split_edge(edge_type_mat, num_vert)

        target_edge_type_mat_wHalfs = self.target_3_edge_type_mat_wHalfs
        target_pseudo_vert = self.target_3_pseudo_vert
        self.assertEqual(actual_pseudo_vert, target_pseudo_vert)
        self.assertEqual(actual_edge_type_mat_wHalfs.nodes(),
                         target_edge_type_mat_wHalfs.nodes())
        self.assertEqual(actual_edge_type_mat_wHalfs.edges(),
                         target_edge_type_mat_wHalfs.edges())

    # 4
    @expectedFailure
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

        actual_scaf_to_edge_adj = adj_scaf_nick_pos(scaf_to_edge, scaf_nick_pos, num_bases)

        target_scaf_to_edge = self.target_8_scaf_to_edge
        self.assertEqual(actual_scaf_to_edge_adj, target_scaf_to_edge)

    # 9
    def test_assign_staples_wChoices(self):
        singleXOs = self.target_0_singleXOs
        edges = self.target_0_edges
        num_edges = len(edges)
        edge_type_mat = self.target_2_edge_type_mat.to_directed()  #TODO:  Same as other time I call `.to_directed`!
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
        len_scaf_used = 2*sum(self.target_0_edge_length_vec)  # TODO: it seems like this variable will be undefined if random scaf seq is used!
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
        self.assertEqual(actual_named_stap_seq_list, target_named_stap_seq_list)

    # 11
    @expectedFailure
    def test_toCanDo(self):
        scaf_to_edge, scaf_seq, stap_list, stap_seq_list, coordinates, edges, edge_length_vec, faces, vert_to_face, fig = [None] * 10
        stuff = toCanDo(scaf_to_edge, scaf_seq, stap_list, stap_seq_list, coordinates,
                        edges, edge_length_vec, faces, vert_to_face, fig)

        self.assertEqual(stuff, "TODO: Bring in real targets.  How to properly test this function?")
