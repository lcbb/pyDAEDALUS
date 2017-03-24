import math
from copy import deepcopy

import mpmath
import numpy as np

from Automated_Design.constants import VERMILLION, REDPURPLE, WHITE, BLU, ORANG
from Automated_Design.util import intersect_lists


# parameters
d = 3.4  # angstroms, distance between two nucleotides
IHD = 20  # angstroms, inter-helical distance
r = 10  # angstroms, radius of double helix
wDX = IHD + 2 * r  # angstroms, width of DX tile


# Running with defining a class rather than a dict to:
# 1) Give you the dot-notation for subparts.  e.g., `.triad`
# 2) Keep the option open for moving some of the later calculations
#    into the DNAGenom class.
class DNAGenom(object):
    def __init__(self, n_bp):
        self.dNode = np.zeros(shape=(n_bp, 3))
        self.triad = np.zeros(shape=(3, 3, n_bp))
        self.id_nt = np.zeros(shape=(n_bp, 2))


class DnaTop(object):
    def __init__(self, ident, up, down, across, seq):
        self.id = ident
        self.up = up
        self.down = down
        self.across = across
        self.seq = seq

    def __str__(self):
        s = "<DnaTop Object. id={}, up={}, down={}, across={}, seq={}>".format(
            self.id, self.up, self.down, self.across, self.seq
        )
        return s

    def __repr__(self):
        return self.__str__()


def calc_buff(faces, num_vert, coordinates, d, wDX):
    """
    Calculates buffer distance between edge of DX tile and center of vertex.

    Parameters
    ----------
    faces : list
        Fx2 cell matrix, where F is the number of faces.
        The first column details how many vertices the face has
        The second column details the vertex IDs of the face
    num_vert : int
        number of vertices, V
    coordinates : numpy.ndarray
        Vx3 matrix of spatial coordinates of vertices,
        V = number of vertices
    d : float
        distance between two nucleotides, in angstroms
    wDX : int
        width of DX tile, in angstroms

    Returns
    -------
    buff_nt
        distance, in nucleotides, between edge of DX arm and vertex coordinate
    """

    vert_angles = [[] for i in range(num_vert)]

    # For each face
    for curr_face in faces:
        face_coords = coordinates[curr_face]

        # for each vert
        for face_vert_ID in range(len(curr_face)):

            center_vert = curr_face[face_vert_ID]
            center_coords = face_coords[face_vert_ID]

            if face_vert_ID == 0:  # if first
                prev_coords = face_coords[-1]
            else:
                prev_coords = face_coords[face_vert_ID - 1]

            if face_vert_ID == len(curr_face) - 1:  # if last
                next_coords = face_coords[0]
            else:
                next_coords = face_coords[face_vert_ID + 1]

            a_vec = np.array(prev_coords - center_coords)
            b_vec = np.array(next_coords - center_coords)

            angle = math.atan2(np.linalg.norm(np.cross(a_vec, b_vec)),
                               np.dot(a_vec, b_vec))

            vert_angles[center_vert] += [angle]

    buff_nt = [0] * num_vert

    # For creating no discontinuities
    for vert_ID in range(num_vert):
        sum_of_angles = sum(vert_angles[vert_ID])
        theta_max = max(vert_angles[vert_ID])
        theta_P_max = theta_max * (2. * np.pi / sum_of_angles)
        r = float(wDX / 2. * mpmath.cot(theta_P_max / 2.))

        if sum_of_angles <= 2. * np.pi:
            s = r * (2. * np.pi / sum_of_angles)
        else:
            s = r

        buff_nt[vert_ID] = int(np.ceil(s / d))

    return buff_nt


def gen_FE_norms(coordinates, faces, edges, vert_to_face):
    """
    Generate vectors normal to Faces and Edges.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Vx3 array of spatial coordinates of vertices,
        V = number of vertices
    faces : list
        Fx2 list, where F is the number of faces.
        The first column details how many vertices the face has.
        The second column details the vertex IDs of the face.
    edges : numpy.ndarray
        Ex2 array where each row corresponds to one edge,
        denoting the vertices being connected. 1st column > 2nd column
    vert_to_face : list
        V-long list, each row containing a list of the face IDs the
        particular vertex belongs to.

    Returns
    -------
    face_norms
        Fx3 matrix containing outward normal for each face.
    edge_norms
        Ex3 matrix containing outward normal for each edge.
    """

    # Initialize
    coordinates = np.array(coordinates)
    face_norms = np.zeros(shape=(len(faces), 3))
    edge_norms = np.zeros(shape=(len(edges), 3))

    for face_ID in range(len(faces)):
        curr_face = np.array(faces[face_ID])  # current face
        curr_coords = coordinates[curr_face[0:3]]  # coordinates of face verts

        a_vec = curr_coords[0] - curr_coords[1]
        b_vec = curr_coords[2] - curr_coords[1]

        normal_vec = np.cross(b_vec, a_vec)  # cross-product for normal

        # normalize
        face_norms[face_ID] = normal_vec / np.linalg.norm(normal_vec)

    for edge_ID in range(len(edges)):
        edge_bgn = edges[edge_ID][0]
        edge_fin = edges[edge_ID][1]

        bordering_faces = intersect_lists(vert_to_face[edge_bgn],
                                          vert_to_face[edge_fin])

        normal_vec = face_norms[bordering_faces[0]] + \
            face_norms[bordering_faces[1]]

        # normalize
        edge_norms[edge_ID] = normal_vec / np.linalg.norm(normal_vec)

    return edge_norms


def tensor_product(a, b):
    """
    Assume a 1 dimensional vector of length n you want to mulply with
    itself to get an (n,n) shaped matrix.
    """
    # a `-1` in numpy array length tells numpy to 'make this dimension be
    # however long it has to be to fit the data into the rest of the supplied
    # dimensions.
    m = a.reshape(1, -1) * b.reshape(-1, 1)
    return m


class DnaInfo(object):
    def __init__(self, scaf_to_edge, scaf_seq, stap_list, stap_seq_list,
                 coordinates, edges, edge_length_vec, faces, vert_to_face):
        """
        Generate file to sent to CanDo for rendering and finite element
        simulation.
        Inputs: scaf_to_edge = Ex2 cell array, where each row corresponds to
                   one edge, 1st column is duplex from low ID to high ID
                   vertex, 2nd column is from high to low. Each element is a
                   row vector containing the scaffold base IDs in order on that
                   duplex.
                scaf_seq = string of scaffold, length may be longer than
                   required. First nucleotide will be placed at index 1
                stap_list = columnar cell array, each cell contains string of
                   scaffold indices that staple bases correspond to
                stap_seq_list = same as stap_list but with sequence instead of
                   index information
                coordinates = Vx3 matrix of spatial coordinates of vertices
                edges = Ex2 matrix where each row corresponds to one edge,
                   denoting the vertices being connected. 1st column > 2nd
                   column
                num_edges = number of edges, E
                edge_length_vec = row vector of edge types, corresponding to
                   edge_length_mat_full
                fig = figure number for scaffold routing display
        Output: dnaInfo = structure that describes sequence, topology, and
        geometry of a programmed DNA assembly. Two substuctures:
                   dnaTop(n) = structure with info of nucleotide n. Fields:
                      id = identification number
                      up = id of nucleotide upstream (towards 5') (-1 if N/A)
                      down = id of nucleotide downstream (towards 3') (-1 if
                      N/A)
                      across = id of base paired nucleotide (-1 if N/A)
                      seq = identity of base (char A,T,C,G)
                   dnaGeom = structure with geometric info of base pairs.
                       Arrays:
                      dNode = n_bp x 3 matrix of spatial location
                      triad = 3 x 3 x n_bp matrix of spatial orientation
                      id_nt = n_bp x 2 matrix of nucleotides that compose base
                         pair
                      Note that single stranded nucleotides are ignored here.
        Note: Helicity of the nucleotides has not been implemented. The
        orientation is correct at each scaffold and staple crossover, but not
        necessarily in between. CanDo pre-mechanical model calculates the
        correct orientation of these nucleotides.
        (see CanDo for more detailed information)
        #######################################################################
        by Sakul Ratanalert, MIT, Bathe Lab, 2016

        Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
        M.I.T. hereby makes following copyrightable material available to the
        public under GNU General Public License, version 2 (GPL-2.0). A copy of
        this license is available at https://opensource.org/licenses/GPL-2.0
        #######################################################################
        """
        stap_seq_list = deepcopy(stap_seq_list)
        stap_list = deepcopy(stap_list)

        # Initialize numbers of vertices, edges, staples, base paired
        # nucleotides:
        num_vert = len(coordinates)
        num_edges = len(edge_length_vec)
        num_stap = len(stap_list)
        n_bp = 2*sum(edge_length_vec)  # number of bp, not length(scaf_seq)
        # because could be inputting seq larger than necessary

        # TODO: Delete this?  Not used elsewhere.
        # Calculate length of each staple
        '''
        len_stap = zeros(num_stap,1)
        for stap_ID in range(len(stap_list)):
            len_stap(stap_ID) = length(stap_list{stap_ID})
        '''

        # Initialize dnaGeom:
        self.dnaGeom = DNAGenom(n_bp=n_bp)

        # Note: preferred nucleotide always on scaffold strand

        # Generate dnaGeom first:

        # Generate buffers at vertices and unit normal vectors at edges
        buff_nt = calc_buff(faces, num_vert, coordinates, d, wDX)
        edge_norms = gen_FE_norms(coordinates, faces, edges, vert_to_face)

        # Scale coordinates. Structure is centered at origin, assume relative
        # lengths already correct.
        edge_to_scale = edges[0]  # arbitrary choose first edge to scale
        vert_1 = edge_to_scale[0]  # first vertex of edge
        vert_2 = edge_to_scale[1]  # second vertex of edge
        edge_length = edge_length_vec[0]  # length of vertex, base pairs
        coord_1 = coordinates[vert_1]  # coordinates of first vertex
        coord_2 = coordinates[vert_2]  # coordinates of second vertex
        current_distance = np.linalg.norm(coord_2 - coord_1)
        desired_distance = d*(buff_nt[vert_1] + edge_length + buff_nt[vert_2])
        scale = desired_distance/current_distance  # scaling factor
        coordinates *= scale  # scale coordinates
        self.scaled_coordinates = coordinates

        # Edges
        for edge_ID in range(num_edges):  # for each edge
            for low_to_high in [0, 1]:  # for each direction per edge
                if low_to_high == 0:
                    edge_bgn = edges[edge_ID][1]
                    edge_fin = edges[edge_ID][0]
                else:
                    edge_bgn = edges[edge_ID][0]
                    edge_fin = edges[edge_ID][1]

                # Obtain coordinates of end vertices
                coord_bgn = coordinates[edge_bgn]
                coord_fin = coordinates[edge_fin]

                z_axis = (coord_fin - coord_bgn)  # row vector
                z_axis = z_axis/np.linalg.norm(z_axis)  # normalize

                scaf_part = scaf_to_edge[edge_ID][low_to_high]  # scaff bases
                len_edge = len(scaf_part)  # length of edge, base pairs

                for bp_ID in range(len_edge):  # for each base pair in edge
                    numerator = buff_nt[edge_bgn] + bp_ID
                    denominator = (buff_nt[edge_bgn] + len_edge +
                                   buff_nt[edge_fin])
                    relative_location = float(numerator) / denominator

                    axis_pos = coord_bgn*(1-relative_location) + \
                        coord_fin*relative_location  # 3D position on edge axis
                    edge_norm = edge_norms[edge_ID]

                    y_axis = np.cross(edge_norm, z_axis)  # row vector
                    y_axis /= np.linalg.norm(y_axis)  # normalize

                    nt_ID = scaf_part[bp_ID]  # nucleotide id
                    axis_pos = axis_pos + y_axis*IHD/2  # shift to correct side

                    # Search for scaffold crossover, 4 base pairs per crossover
                    scaf_nick_pos = 0

                    if bp_ID > 0 and bp_ID < (len_edge-1):
                        nt_down = scaf_part[bp_ID + 1]

                        # set larger int type in scaf test data to negate
                        # the need for the following `int()`:
                        five_prime_side_of_scaff_nick = \
                            abs(int(nt_down) - nt_ID) > 1
                        not_scaffold_nick = nt_ID > 1 and nt_ID < n_bp
                        if five_prime_side_of_scaff_nick and not_scaffold_nick:
                            scaf_nick_pos = bp_ID

                    x_axis = np.cross(y_axis, z_axis)  # calculate x-axis

                    # Store in dnaGeom
                    self.dnaGeom.dNode[nt_ID] = axis_pos
                    self.dnaGeom.triad[:, 0, nt_ID] = x_axis
                    self.dnaGeom.triad[:, 1, nt_ID] = y_axis
                    self.dnaGeom.triad[:, 2, nt_ID] = z_axis
                    self.dnaGeom.id_nt[nt_ID] = [nt_ID, nt_ID + n_bp]

                # tensor product of z, row vector
                zXz = tensor_product(z_axis, z_axis)
                assert zXz.shape == (3, 3)

                zX = np.array([
                        [0,          -z_axis[2],  z_axis[1]],
                        [z_axis[2],           0, -z_axis[0]],
                        [-z_axis[1],  z_axis[0],          0]
                    ])

                if scaf_nick_pos > 0:
                    for bp_ID in range(scaf_nick_pos):
                        num_nt = scaf_nick_pos
                        start_nt = 0
                        nt_ID = scaf_part[bp_ID]  # nucleotide id

                        turnAngle = self.get_turn_angle_for_pos_scaf_nick_pos(
                            num_nt)

                        rot_matrix = self.get_rot_matrix(bp_ID, start_nt,
                                                         turnAngle, zX, zXz)

                        prev_triad = self.dnaGeom.triad[:, :, nt_ID]
                        y_axis = prev_triad[:, 1]

                        y_axis = np.inner(rot_matrix, y_axis)

                        x_axis = np.cross(y_axis, z_axis)  # calculate x-axis

                        # Store in dnaGeom
                        self.dnaGeom.triad[:, 0, nt_ID] = x_axis
                        self.dnaGeom.triad[:, 1, nt_ID] = y_axis
                        self.dnaGeom.triad[:, 2, nt_ID] = z_axis

                    for bp_ID in range(scaf_nick_pos, len_edge):
                        num_nt = len_edge-scaf_nick_pos
                        start_nt = scaf_nick_pos
                        nt_ID = scaf_part[bp_ID]  # nucleotide id

                        turnAngle = self.get_turn_angle_for_pos_scaf_nick_pos(
                            num_nt)

                        rot_matrix = self.get_rot_matrix(bp_ID, start_nt,
                                                         turnAngle, zX, zXz)

                        prev_triad = self.dnaGeom.triad[:, :, nt_ID]
                        y_axis = prev_triad[:, 1]

                        y_axis = np.inner(rot_matrix, y_axis)

                        x_axis = np.cross(y_axis, z_axis)  # calculate x-axis

                        # Store in dnaGeom
                        self.dnaGeom.triad[:, 0, nt_ID] = x_axis
                        self.dnaGeom.triad[:, 1, nt_ID] = y_axis
                        self.dnaGeom.triad[:, 2, nt_ID] = z_axis

                else:
                    for bp_ID in range(len_edge):
                        num_nt = len_edge
                        start_nt = 0
                        nt_ID = scaf_part[bp_ID]  # nucleotide id

                        turnAngle = self.get_turn_angle_for_zero_scaf_nick_pos(
                            num_nt)  # scalar

                        rot_matrix = self.get_rot_matrix(bp_ID, start_nt,
                                                         turnAngle, zX, zXz)

                        prev_triad = self.dnaGeom.triad[:, :, nt_ID]
                        y_axis = prev_triad[:, 1]

                        y_axis = np.inner(rot_matrix, y_axis)

                        x_axis = np.cross(y_axis, z_axis)  # calculate x-axis

                        # Store in dnaGeom
                        self.dnaGeom.triad[:, 0, nt_ID] = x_axis
                        self.dnaGeom.triad[:, 1, nt_ID] = y_axis
                        self.dnaGeom.triad[:, 2, nt_ID] = z_axis

        # Initialize dnaTop
        self.dnaTop = []

        for scaf_seq_ID in range(n_bp):
            # Set id in dnaTop
            this_id = scaf_seq_ID

            # Set up in dnaTop
            if scaf_seq_ID == 0:
                this_up = -1
            else:
                this_up = scaf_seq_ID-1

            # Set down in dnaTop
            if scaf_seq_ID == n_bp - 1:
                this_down = -1
            else:
                this_down = scaf_seq_ID+1

            # Set across in dnaTop
            this_across = scaf_seq_ID + n_bp

            # Set seq in dnaTop
            this_seq = scaf_seq[scaf_seq_ID].upper()

            self.dnaTop.append(
                DnaTop(this_id, this_up, this_down, this_across, this_seq))

        # Staples:
        polyT_ID = 2*n_bp  # initialize polyT counter at end of all paired nts

        for stap_ID in range(num_stap):  # for each staple
            stap = stap_list[stap_ID]

            # IDs: Change polyT from 0 to nonzero and staple to staple + n_bp
            for stap_base in range(len(stap)):
                if stap[stap_base] is None:  # part of polyT
                    stap[stap_base] = polyT_ID  # change from 0 to polyT

                    while not polyT_ID < len(self.dnaTop):
                        self.dnaTop.append(
                            DnaTop(None, None, None, None, None))

                    self.dnaTop[polyT_ID].across = -1  # no across
                    polyT_ID += 1  # increment
                else:
                    # new staple base ID:
                    stap[stap_base] = stap[stap_base] + n_bp

            for stap_base in range(len(stap)):  # for each base in staple

                stap_seq_ID = stap[stap_base]  # ID of staple base
                while not stap_seq_ID < len(self.dnaTop):
                    self.dnaTop.append(DnaTop(None, None, None, None, None))

                # Set id in dnaTop
                self.dnaTop[stap_seq_ID].id = stap_seq_ID

                # Set up in dnaTop
                if stap_base == 0:  # if 5' end
                    self.dnaTop[stap_seq_ID].up = -1
                else:
                    self.dnaTop[stap_seq_ID].up = stap[stap_base-1]

                # Set down in dnaTop
                if stap_base == len(stap)-1:  # if 3' end
                    self.dnaTop[stap_seq_ID].down = -1
                else:
                    self.dnaTop[stap_seq_ID].down = stap[stap_base+1]

                # Set across in dnaTop
                if self.dnaTop[stap_seq_ID].across is None:  # not in polyT
                    self.dnaTop[stap_seq_ID].across = stap_seq_ID - n_bp

                # Set upper in dnaTop
                self.dnaTop[stap_seq_ID].seq = \
                    stap_seq_list[stap_ID][stap_base].upper()

    @staticmethod
    def get_turn_angle_for_zero_scaf_nick_pos(num_nt):
        return 2 * np.pi * int(round(num_nt / 10.5)) / (num_nt + 1)

    @staticmethod
    def get_rot_matrix(bp_ID, start_nt, turnAngle, zX, zXz):
        rot_matrix = np.cos(turnAngle * (bp_ID - start_nt + 0.5)) * np.eye(3) \
                     + np.sin(turnAngle * (bp_ID - start_nt + 0.5)) * zX \
                     + (1 - np.cos(turnAngle * (bp_ID - start_nt + 0.5))) * zXz
        return rot_matrix

    @staticmethod
    def get_turn_angle_for_pos_scaf_nick_pos(num_nt):
        return 2 * np.pi * (int(np.floor(num_nt / 10.5)) + 0.5) / (num_nt + 1)

    def plot_3d_model(self, filename):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        fig = plt.figure(3, figsize=(12, 12))
        fig.clf()
        ax = fig.add_subplot(111, projection='3d')

        # Plot numbered circles at each vertex
        first_d = [loc[0] for loc in self.scaled_coordinates]
        second_d = [loc[1] for loc in self.scaled_coordinates]
        third_d = [loc[2] for loc in self.scaled_coordinates]
        ax.scatter(first_d, second_d, third_d,
                   marker='o',  # marker shape
                   color=VERMILLION,  # marker color
                   edgecolor=VERMILLION,
                   s=200)  # marker size

        # Plot positions of each nucleotide as a line
        first_d = [dNode[0] for dNode in self.dnaGeom.dNode]
        second_d = [dNode[1] for dNode in self.dnaGeom.dNode]
        third_d = [dNode[2] for dNode in self.dnaGeom.dNode]
        ax.plot(first_d, second_d, third_d,
                '-',
                color=BLU,
                linewidth=5)

        first_d = [loc[0] for loc in self.scaled_coordinates]
        second_d = [loc[1] for loc in self.scaled_coordinates]
        third_d = [loc[2] for loc in self.scaled_coordinates]
        for i, (first, second, third) in enumerate(zip(
                first_d, second_d, third_d)):
            ax.text(first, second, third,
                    str(i), color=WHITE,
                    horizontalalignment='center', verticalalignment='center')

        # Plot 5' end (orange square)
        ax.scatter([self.dnaGeom.dNode[0][0]],
                   [self.dnaGeom.dNode[0][1]],
                   [self.dnaGeom.dNode[0][2]],
                   marker='s',
                   color=ORANG,
                   edgecolor=ORANG,
                   s=500)
        # Plot 3' end (purple circle)
        ax.scatter([self.dnaGeom.dNode[-1][0]],
                   [self.dnaGeom.dNode[-1][1]],
                   [self.dnaGeom.dNode[-1][2]],
                   marker='o',
                   color=REDPURPLE,
                   edgecolor=REDPURPLE,
                   s=400)  # 'LineWidth',4,

        ax.set_xlabel('Angstroms', fontdict={'size': 16})
        ax.set_ylabel('Angstroms', fontdict={'size': 16})
        ax.set_zlabel('Angstroms', fontdict={'size': 16})

        plt.savefig(filename, bbox_inches='tight', pad_inches=0)

    def save_dna_info_to_cando_file(self, cando_filename):
        """
        Converts dnaInfo into a CanDo text file format.
        See http://cando-dna-origami.org/cndo-file-converter/ for file format

        WARNING!  While all of the internal bits here are 0-indexed, the
        output file needs to be 1-indexed.
        """

        def handle_1_indexing(val):
            one_indexed = True
            if one_indexed:
                val = int(val) + 1 if val >= 0 else val
            return val

        fid = open(cando_filename, 'w')

        # File header
        fid.write(
            '"CanDo (.cndo) file format version 1.0, Keyao Pan, Laboratory '
            'for Computational Biology and Biophysics, Massachusetts '
            'Institute of Technology, November 2015"\n')
        fid.write('\n')

        # dnaInfo.dnaTop
        fid.write('dnaTop,id,up,down,across,seq\n')
        for i, topology in enumerate(self.dnaTop):
            fid.write('{},{},{},{},{},{}\n'.format(
                handle_1_indexing(i),
                handle_1_indexing(topology.id),
                handle_1_indexing(topology.up),
                handle_1_indexing(topology.down),
                handle_1_indexing(topology.across),
                topology.seq))
        fid.write('\n')

        # dnaInfo.dnaGeom.dNode
        fid.write('dNode,"e0(1)","e0(2)","e0(3)"\n')
        for i, node in enumerate(self.dnaGeom.dNode):
            fid.write('{},{},{},{}\n'.format(
                handle_1_indexing(i),
                node[0],
                node[1],
                node[2]))
        fid.write('\n')

        # dnaInfo.dnaGeom.triad
        fid.write('triad,"e1(1)","e1(2)","e1(3)","e2(1)",'
                  '"e2(2)","e2(3)","e3(1)","e3(2)","e3(3)"\n')
        for i in range(self.dnaGeom.triad.shape[-1]):
            triad = self.dnaGeom.triad[:, :, i]
            fid.write('{},{},{},{},{},{},{},{},{},{}\n'.format(
                handle_1_indexing(i),
                -triad[0, 0],
                -triad[1, 0],
                -triad[2, 0],
                triad[0, 1],
                triad[1, 1],
                triad[2, 1],
                -triad[0, 2],
                -triad[1, 2],
                -triad[2, 2]))
        fid.write('\n')

        # dnaInfo.dnaGeom.id_nt
        fid.write('id_nt,id1,id2\n')
        for i, id_nt in enumerate(self.dnaGeom.id_nt):
            fid.write('{},{},{}\n'.format(
                handle_1_indexing(i),
                handle_1_indexing(id_nt[0]),
                handle_1_indexing(id_nt[1])))

        fid.close()
