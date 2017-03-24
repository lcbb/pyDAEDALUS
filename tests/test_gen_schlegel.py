import numpy as np

from Automated_Design.gen_schlegel import create_2d_mapping


def coords_are_almost_equal(a, b, threshold=0.001):
    a_x, a_y = a
    b_x, b_y = b
    x_almost_equal = abs(a_x - b_x) < threshold
    y_almost_equal = abs(a_y - b_y) < threshold
    return x_almost_equal and y_almost_equal


def test_create_2d_mapping_using_tetrahedron():
    edges = np.array([[2, 0],
                      [1, 0],
                      [3, 1],
                      [3, 0],
                      [2, 1],
                      [3, 2]])
    coordinates = np.array([[0., 0., 0.612372],
                            [-0.288675, -0.5, -0.204124],
                            [-0.288675, 0.5, -0.204124],
                            [0.57735, 0., -0.204124]])
    faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    xycoords = create_2d_mapping(edges, coordinates, faces)

    # Ideally, these assertions would assure that three of the four nodes
    # are arranged in an equilateral triangle scaled to sit on top of the
    # unit circle with the fourth node exactly in the middle of the other
    # three.  But, since I'm not clever enough to know how to write those
    # assertions exactly, I'm starting with asserting it results in both a
    # possible solution and the solution I know it will reach given how
    # the function currently works.

    m = 3.**(1./2.)/2

    # The fourth node happened to be picked as center:
    assert coords_are_almost_equal([0, 0], xycoords[3])
    # The first three nodes all surround the fourth node.
    # The first node happened to be picked to be directly above center:
    assert coords_are_almost_equal([1, 0], xycoords[0])
    # The second is bottom left:
    assert coords_are_almost_equal([-.5, -m], xycoords[1])
    # The third is bottom right:
    assert coords_are_almost_equal([-.5, m], xycoords[2])
