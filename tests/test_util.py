from Automated_Design.util import generate_graph


class TestGenerateGraph:
    def test_number_of_vertices(self):
        num_vert = 5
        edges = []
        edge_length_vec = []
        G = generate_graph(num_vert, edges, edge_length_vec)
        assert num_vert == len(G.nodes())
