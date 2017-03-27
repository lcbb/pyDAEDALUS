import numpy as np

from Automated_Design.util import generate_graph, intersect_lists, find


class TestGenerateGraph:
    def test_number_of_vertices(self):
        num_vert = 5
        edges = []
        edge_length_vec = []
        G = generate_graph(num_vert, edges, edge_length_vec)
        assert num_vert == len(G.nodes())

    def test_edge_length_is_accesble_on_resulting_graph(self):
        num_vert = 3
        edges = [[0, 1],
                 [1, 2]]
        edge_length_vec = [3, 5]
        G = generate_graph(num_vert, edges, edge_length_vec)
        assert G[1][2]['length'] == 5.


class TestIntersectLists:
    def test_no_match(self):
        a = [1]
        b = [2, 3]
        y = intersect_lists(a, b)
        assert y == []

    def test_partial_match(self):
        a = [1, 2]
        b = [2, 3]
        y = intersect_lists(a, b)
        assert y == [2]

    def test_full_match(self):
        a = [1, 2, 3]
        b = [1, 2, 3]
        y = intersect_lists(a, b)
        assert y == [1, 2, 3]

    def test_input_can_be_numpy_array(self):
        a = np.array([1, 2, 3])
        b = np.array([2, 3, 4, 5])
        y = intersect_lists(a, b)
        assert y == [2, 3]


class TestFind:
    def test_nothing_to_find(self):
        iterable = [1, 3, 3, 5, 6]
        val = 2
        y = find(iterable, val)
        assert y == []

    def test_only_a_few_hits(self):
        iterable = [1, 3, 3, 5, 6]
        val = 3
        y = find(iterable, val)
        assert y == [1, 2]

    def test_everything_is_a_match(self):
        iterable = [3, 3, 3, 3]
        val = 3
        y = find(iterable, val)
        assert y == [0, 1, 2, 3]

    def test_input_can_be_numpy_array(self):
        iterable = np.array([[1, 2], [3, 4], [4, 5]])
        val = np.array([4, 5])
        y = find(iterable, val)
        assert y == [2]
