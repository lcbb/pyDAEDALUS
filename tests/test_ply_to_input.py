from StringIO import StringIO
from os import path

from mock import patch, MagicMock, mock_open

from Automated_Design.ply_to_input import ply_to_input, \
    extract_file_reader_and_shape_name_from_input_filename


def get_first_call_first_param(mock):
    return mock.call_args[0][0]


class TestExtractFileReaderAndShapeNameFromInputFilename:
    def setup_method(self, method):
        # Always mock out `open()` so you don't actually hit the file system.
        # You also need to mock out `path.is_file()` to both have it not
        # really check for a file and to be able to declare its intended
        # resolution.

        self.open_patcher = patch('__builtin__.open', mock_open())
        self.open_mock = self.open_patcher.start()

        self.is_filename_patcher = patch(
            'os.path.isfile', new=MagicMock(return_value=True))
        self.is_filename_mock = self.is_filename_patcher.start()

    def teardown_method(self, method):
        self.open_patcher.stop()
        self.is_filename_patcher.stop()

    # Test the part that optionally adds the .ply extension:
    def test_filename_with_ply_works(self):
        f, shape_name = extract_file_reader_and_shape_name_from_input_filename(
            'some_file')

        # Assert the file being opened has exactly one .ply at the end
        filename = get_first_call_first_param(self.is_filename_mock)
        assert filename[-4:] == '.ply'

        # Assert the shape name has no .ply
        assert shape_name[-4:] != '.ply'

    def test_filename_without_dot_ply_works(self):
        f, shape_name = extract_file_reader_and_shape_name_from_input_filename(
            'some_file.ply')

        # Assert the file being opened has exactly one .ply at the end
        filename = get_first_call_first_param(self.is_filename_mock)
        assert filename[-4:] == '.ply'

        # Assert the shape name has no .ply
        assert shape_name[-4:] != '.ply'

    # Test the part that takes the `.ply`-less filename and also removes the
    # preceding folder structure to give just the shape name.
    def test_shape_name_is_extracted_from_filename_within_folders(self):
        filename = path.join('folder1', 'folder2', 'some_shape.ply')

        f, shape_name = extract_file_reader_and_shape_name_from_input_filename(
            filename)

        assert shape_name == 'some_shape'


def open_string_as_file(string):
    file = StringIO()
    file.write(string)
    file.seek(0)
    return file


class TestPlyImportOnTetrahedron:
    @classmethod
    def setup_class(cls):
        cls.ply_file_01_tetrahedron = """
            ply
            format ascii 1.0
            element vertex 4
            property float32 x
            property float32 y
            property float32 z
            element face 4
            property list uint8 int32 vertex_indices
            end_header
            0.000000 0.000000 0.612372
            -0.288675 -0.500000 -0.204124
            -0.288675 0.500000 -0.204124
            0.577350 0.000000 -0.204124
            3 0 2 1
            3 0 1 3
            3 0 3 2
            3 1 2 3
        """
        cls.shape_name = '01_tetrahedron'

    def setup_method(self, method):
        # You need to re-initialize the 'file reader' at the beginning
        # of every function, else one test will read the file and the next
        # test will begin reading the file at the end, find no data, and hang
        # indefinitely.
        f = open_string_as_file(self.ply_file_01_tetrahedron)

        self.extract_file_and_reader_patcher = patch(
            'Automated_Design.ply_to_input'
            '.extract_file_reader_and_shape_name_from_input_filename',
            new=MagicMock(return_value=[f, self.shape_name]))
        self.extract_file_and_reader_mock = \
            self.extract_file_and_reader_patcher.start()

    def teardown_method(self, method):
        self.extract_file_and_reader_patcher.stop()

    def test_list_lengths_are_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
        assert len(coordinates) == 4
        assert len(faces) == 4

    def test_arbitrary_coordinate_is_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
        assert list(coordinates[1]) == [-0.288675, -0.500000, -0.204124]

    def test_arbitrary_face_is_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
        assert faces[3] == [1, 2, 3]

    # singleXOs tests:
    def test_single_xos_are_not_used_when_min_len_nt_is_short(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=41)
        assert singleXOs == 0

    def test_single_xos_used_when_min_len_nt_is_long(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=43)
        assert singleXOs == 1

    # Plotting tests:
    # On these two tests, you mock out the plotting function so only the mock
    # and not the plotting function itself is called.  That way, the sample
    # foldername defined in the second tests doesn't have to exist and you
    # don't waste time re-generating plots every time you run tests.
    @patch('Automated_Design.ply_to_input.plot_edge_length_distributions')
    def test_plot_function_is_called_if_results_foldername_set(
            self, plotter_mock):
        ply_to_input("01_tetrahedron", results_foldername='')
        assert plotter_mock.call_count == 0

    @patch('Automated_Design.ply_to_input.plot_edge_length_distributions')
    def test_plot_function_is_not_called_if_results_foldername_not_set(
                self, plotter_mock):
        ply_to_input("01_tetrahedron", results_foldername='i_am_defined')

        assert plotter_mock.call_count == 1


class TestPlyImportOnTetrahedronWithUnusedVertices:
    @classmethod
    def setup_class(cls):
        # This is a tetrahedron with the 4th vertex becoming the 5th and
        # the added 4th vertex is left unused in the faces given.
        cls.ply_file_01_tetrahedron = """
            ply
            format ascii 1.0
            element vertex 5
            property float32 x
            property float32 y
            property float32 z
            element face 4
            property list uint8 int32 vertex_indices
            end_header
            0.000000 0.000000 0.612372
            -0.288675 -0.500000 -0.204124
            -0.288675 0.500000 -0.204124
            1137. 113.7 11.37
            0.577350 0.000000 -0.204124
            3 0 2 1
            3 0 1 4
            3 0 4 2
            3 1 2 4
        """
        cls.shape_name = '01_tetrahedron'

    def setup_method(self, method):
        f = open_string_as_file(self.ply_file_01_tetrahedron)

        self.extract_file_and_reader_patcher = patch(
            'Automated_Design.ply_to_input'
            '.extract_file_reader_and_shape_name_from_input_filename',
            new=MagicMock(return_value=[f, self.shape_name]))
        self.extract_file_and_reader_mock = \
            self.extract_file_and_reader_patcher.start()

    def teardown_method(self, method):
        self.extract_file_and_reader_patcher.stop()

    def test_unused_vertex_is_removed(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)

        # Note there are 5 coordinates in the starting file:
        assert len(coordinates) == 4

        assert len(faces) == 4  # Shouldn't have changed

        # Started as coordinates[4]:
        assert list(coordinates[3]) == [0.577350, 0.000000, -0.204124]

        # Started as [1, 2, 4]
        assert faces[3] == [1, 2, 3]


class TestPlyImportOnIcosahedron:
    """
    This is here to make sure a larger file still generally works, though
    more detailed testing was already handled above.
    """
    @classmethod
    def setup_class(cls):
        cls.ply_file_05_icosahedron = """
            ply
            format ascii 1.0
            element vertex 12
            property float32 x
            property float32 y
            property float32 z
            element face 20
            property list uint8 int32 vertex_indices
            end_header
            0.000000 0.000000 1.176000
            1.051000 0.000000 0.526000
            0.324000 1.000000 0.525000
            -0.851000 0.618000 0.526000
            -0.851000 -0.618000 0.526000
            0.325000 -1.000000 0.526000
            0.851000 0.618000 -0.526000
            0.851000 -0.618000 -0.526000
            -0.325000 1.000000 -0.526000
            -1.051000 0.000000 -0.526000
            -0.325000 -1.000000 -0.526000
            0.000000 0.000000 -1.176000
            3 0 1 2
            3 0 2 3
            3 0 3 4
            3 0 4 5
            3 0 5 1
            3 1 5 7
            3 1 7 6
            3 1 6 2
            3 2 6 8
            3 2 8 3
            3 3 8 9
            3 3 9 4
            3 4 9 10
            3 4 10 5
            3 5 10 7
            3 6 7 11
            3 6 11 8
            3 7 10 11
            3 8 11 9
            3 9 11 10
            """
        cls.shape_name = '05_icosahedron'

    def setup_method(self, method):
        # You need to re-initialize the 'file reader' at the beginning
        # of every function, else one test will read the file and the next
        # test will begin reading the file at the end, find no data, and hang
        # indefinitely.
        f = open_string_as_file(self.ply_file_05_icosahedron)

        self.extract_file_and_reader_mock = patch(
            'Automated_Design.ply_to_input'
            '.extract_file_reader_and_shape_name_from_input_filename',
            new=MagicMock(return_value=[f, self.shape_name]))
        self.extract_file_and_reader_mock.start()

    def teardown_method(self, method):
        self.extract_file_and_reader_mock.stop()

    def test_list_lengths_are_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("05_icosahedron", '', min_len_nt=52)
        assert len(coordinates) == 12
        assert len(faces) == 20

    def test_arbitrary_coordinate_is_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
        assert list(coordinates[-5]) == [0.851000, -0.618000, -0.526000]

    def test_arbitrary_face_is_correct(self):
        coordinates, edges, faces, edge_length_vec, file_name, staple_name, \
            singleXOs = ply_to_input("01_tetrahedron", '', min_len_nt=52)
        assert faces[11] == [3, 9, 4]
