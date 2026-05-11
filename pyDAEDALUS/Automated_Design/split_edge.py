def split_edge(graph_with_spanning_tree):
    """
    Add two nodes for each nontree edge to implement scaffold crossovers.
    Added nodes (pseudo-vertices) have a reference vertex (one of the V real
    vertices) to maintain relative spatial coordinates among all
    vertices.

    Parameters
    ----------
    graph_with_spanning_tree : networkx.classes.digraph.DiGraph
        sparse matrix where
            1 is non-spanning tree edge: DX edge with 1 scaffold crossover
            2 is spanning tree edge: DX edge with 0 scaffold crossovers
    Returns
    -------
    graph_with_edges_split
        VxV sparse matrix where
            -1 is half of a non-spanning tree edge (one side of scaffold
            crossover)
            2 is spanning tree edge: DX edge with 0 scaffold crossovers
    pseudo_vert
        row vector where value j at index i indicate that vertex i corresponds
        to vertex j, one of the V real vertices
    """
    # # Initialize output variables. These will be augmented from input vars
    graph_with_edges_split = graph_with_spanning_tree.copy()

    # Start with known nodes.  Pseudo nodes added to graph as needed.
    pseudo_vert = list(graph_with_edges_split)
    # ^updated from below line; new networkx does not allow appending 
    # to .NodeView object.
    # pseudo_vert = graph_with_edges_split.nodes() 
    
    non_tree_edges = [edge for edge in graph_with_spanning_tree.edges(data=True)
                      if edge[2]['type'] == 1]

    # # ID non-tree edges (value 1), replace with 2 half-edges (value = -1)
    # edges = graph_with_spanning_tree.edges()
    # for i in range(num_vert):
    #     for j in range(num_vert):
    #         if (i, j) in edges:
    #             properties = graph_with_spanning_tree[i][j]
    #             if properties['type'] == 1:
    # !! replace the next two lines with the above 6 if you want an exact
    # match with the matlab code.  (else, order of pseudonodes not assured)
    for edge in non_tree_edges:
        i, j, properties = edge
        # `i` is the node that stays. `j` is cut off and replaced by pseudo
        # vert.  Since this digraph was built from a graph, a link from 2, 3
        # will also have a link from 3, 2.  Meaning in the 2, 3 case, you'll
        # keep 2 and replace 3.  And in the 3, 2 case, you'll keep 3 and
        # replace 2.

        index_of_new_vert = len(pseudo_vert)
        properties['type'] = -1

        pseudo_vert.append(j)  # vert_list[new_node_id] = original_node_id
        graph_with_edges_split.remove_edge(i, j)

        graph_with_edges_split.add_node(index_of_new_vert)
        graph_with_edges_split.add_edges_from([(i, index_of_new_vert,
                                      properties)])
        graph_with_edges_split.add_edges_from([(index_of_new_vert, i,
                                      properties)])

    return graph_with_edges_split, pseudo_vert
