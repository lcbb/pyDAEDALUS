"""
All of the functions in this file are used to read in the matlab-exported
code and convert it into a format that can be directly compared to the
python functions within this directory.
"""
from os import path

import networkx as nx
import numpy as np
from scipy import io as sio

from Automated_Design.dna_info import DnaTop

TARGETS_FOLDERNAME = 'Targets01Tetrahedron'


def load_mat_file(filename):
    full_data = sio.loadmat(path.join('tests', TARGETS_FOLDERNAME, filename))
    # full_data.keys()
    # >>> ['thing', '__version__', '__header__', '__globals__']
    # where `thing` is whatever var was originally saved.
    key = (key for key in full_data.keys() if key[0] is not '_').next()
    data = full_data[key]
    return data
    # TODO: declare all these files currently used as 1_tetra...


def load_graph_from_mat(filename, edge_attribute='type',
                        graph_type=nx.DiGraph()):
    graph_as_sparse_matrix = load_mat_file(filename)
    graph = nx.from_scipy_sparse_matrix(graph_as_sparse_matrix,
                                        create_using=graph_type,
                                        edge_attribute=edge_attribute)
    # intify 'type' if that's the property used (defaults to float)
    if edge_attribute == 'type':
        for i, j, attribute in graph.edges(data=True):
            graph[i][j]['type'] = int(attribute['type'])
    return graph.copy()


def load_pseudonodes_from_mat(filename):
    pseudonodes = load_mat_file(filename)
    pseudonodes = pseudonodes - 1
    return list(pseudonodes.flatten())  # convert shape(N,1) to shape(N)


def load_1d_list_from_mat(filename, are_1_indexed_nodes=False):
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
        parts = [parts[0]] + parts[1].split('-')
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

    # take 1 indexing in names down to 0 indexing
    splitter = named_stap_seq_list.index([None, None])
    for row in named_stap_seq_list[:splitter]:
        row[0] = de_1_index_name_a(row[0])

    for row in named_stap_seq_list[splitter+1:]:
        row[0] = de_1_index_name_b(row[0])

    return named_stap_seq_list


def load_dna_info(filename):
    """
    There has to be a much better way to do this.  Namely, using the info you'd
    get out of `thing.dtype` every step of the way (I don't know if order is
    consistent across .mat targets if regenerated).  But it's probably not
    worth figuring out, since these tests should have initial states hard coded
    into the tests rather than depending on the from-matlab dumps.
    """
    full_data = sio.loadmat(path.join('tests', TARGETS_FOLDERNAME, filename))
    data = full_data['dnaInfo'][0][0]
    rawDnaTop = data[0][0]
    dnaTop = []
    for row in rawDnaTop:
        this_id = row[0][0][0] - 1
        this_up = row[1][0][0]
        this_up = this_up - 1 if this_up > 0 else -1
        this_down = row[2][0][0]
        this_down = this_down - 1 if this_down > 0 else -1
        this_across = row[3][0][0]
        this_across = this_across - 1 if this_across > 0 else -1
        this_seq = row[4][0][0]
        dnaTop.append(
            DnaTop(this_id, this_up, this_down, this_across, this_seq)
        )

    dnaGenom = data[1][0][0]
    dNode = dnaGenom[0]
    triad = dnaGenom[1]
    id_nt = dnaGenom[2] - 1

    dna_info = {'dnaTop': dnaTop,
                'dnaGeom': {
                    'dNode': dNode,
                    'triad': triad,
                    'id_nt': id_nt,
                }}
    return dna_info
