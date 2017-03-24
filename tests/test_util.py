from Automated_Design.util import generate_graph


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
