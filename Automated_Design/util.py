import networkx as nx
import numpy as np


def generate_graph(num_vert, edges, edge_length_vec=None):
    G = nx.Graph()
    for n in range(num_vert):
        G.add_node(n)
    for edge, length in zip(edges, edge_length_vec):
        if edge_length_vec:
            G.add_edge(edge[0], edge[1], length=length)
        else:
            G.add_edge(edge[0], edge[1])
        # TODO: verfiy the length property is not used as 'weight' in spanning tree alg.

    full_graph = G
    return full_graph


def intersect_lists(a, b):
    #in case `a` and `b` are np.arrays.
    if type(a) == np.ndarray:
        a = list(a.flatten())
    if type(b) == np.ndarray:
        b = list(b.flatten())
    #TODO: rewrite in a way that preverves order they're seen in?

    thing = set(a).intersection(set(b))
    return sorted(list(thing))


def find(iterable, val):
    return [i for i in range(len(iterable)) if iterable[i] == val]