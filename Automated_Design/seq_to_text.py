from datetime import datetime


def seqtoText(scaf_to_edge, edges, dnaInfo, file_name, scaf_name, singleXOs,
              full_file_name):
    """
    Generate text file to visualize each edge's sequences and nick/crossover
    information
    Inputs: scaf_to_edge = Ex2 cell array, where each row corresponds to one
               edge, 1st column is duplex from low ID to high ID vertex,
               2nd column is from high to low. Each element is a row vector
               containing the scaffold base IDs in order on that duplex.
            edges = Ex2 matrix where each row corresponds to one edge,
               denoting the vertices being connected. 1st column > 2nd column
            dnaInfo = MATLAB file containing all spatial and routing info
            file_name = name of structure
            scaf_name = name of scaffold
            singleXOs = 1 if single crossover vertex staples should be used,
               0 if double crossover vertex staples should be used.
            full_file_name = file_name with more information
    Output: Text file to visualize each edge's sequences and nick/crossover
              information
    ##########################################################################
    by Sakul Ratanalert, MIT, Bathe Lab, 2016

    Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    M.I.T. hereby makes following copyrightable material available to the
    public under GNU General Public License, version 2 (GPL-2.0). A copy of
    this license is available at https://opensource.org/licenses/GPL-2.0
    ##########################################################################
    """

    # TODO: Add bit about 0-indexing!!

    fid = open(full_file_name, 'w')

    num_edges = len(edges)

    fid.write('Structure: {}\r\n'.format(file_name))
    fid.write('Scaffold: {}\r\n'.format(scaf_name))

    if singleXOs > 0:
        fid.write('Vertex staples: Single crossovers\r\n')
    else:
        fid.write('Vertex staples: Double crossovers\r\n')

    fid.write('Produced: \r\n'.format(str(datetime.now())))
    fid.write(
        'Key:    outer strands = scaffold (sequence and index shown)\r\n')
    fid.write('        inner strands = staples\r\n')
    fid.write('                    | = nucleotide\r\n')
    fid.write('                    . = 5'' end\r\n')
    fid.write('               v or ^ = 3'' end\r\n')
    fid.write('   +>>+   or   +<<+   = staple crossover\r\n')
    fid.write(' >> >> >> or << << << = scaffold crossover\r\n')
    fid.write('Not shown: polyT loops in staple strands between edges\r\n\r\n')

    for edge_ID in range(num_edges):

        edge_bgn = edges[edge_ID][1]
        edge_fin = edges[edge_ID][0]

        scaf_right = scaf_to_edge[edge_ID][0]
        scaf_left = scaf_to_edge[edge_ID][1]

        fid.write('Vertex {} (top) to {} (bottom), Edge {}\r\n\r\n'.format(
            edge_bgn, edge_fin, edge_ID))
        fid.write("         3'5' 3'5' \r\n")
        fid.write('              | | {} {}\r\n'.format(
            dnaInfo.dnaTop[scaf_right[0]].seq, scaf_right[0]
        ))

        # Cut off first of scaf_right and invert scaf_left
        scaf_right = scaf_right[1:-1]
        scaf_left = scaf_left[-1::-1]

        for nt_ID in range(len(scaf_right)):
            # Left nucleotides
            scaf_left_ID = scaf_left[nt_ID]
            stap_left_ID = dnaInfo.dnaTop[scaf_left_ID].across

            # Right nucleotides
            scaf_right_ID = scaf_right[nt_ID]
            stap_right_ID = dnaInfo.dnaTop[scaf_right_ID].across

            #    text pos '01234567'
            middle = list('| |  | |')

            # Left scaffold nick or crossover?
            if dnaInfo.dnaTop[scaf_left_ID].up == -1:
                middle[0] = '.'
            elif dnaInfo.dnaTop[scaf_left_ID].down == -1:
                middle[0] = '^'
            elif dnaInfo.dnaTop[scaf_left_ID].up == scaf_right_ID:
                middle[0] = '<'
                middle[1] = '<'
                middle[3] = '<'
            elif dnaInfo.dnaTop[scaf_left_ID].down == scaf_right_ID:
                middle[0] = '>'
                middle[1] = '>'
                middle[3] = '>'

            # Right scaffold nick or crossover?
            if dnaInfo.dnaTop[scaf_right_ID].up == -1:
                middle[7] = '.'
            elif dnaInfo.dnaTop[scaf_right_ID].down == -1:
                middle[7] = 'v'
            elif dnaInfo.dnaTop[scaf_right_ID].up == scaf_left_ID:
                middle[4] = '>'
                middle[6] = '>'
                middle[7] = '>'
            elif dnaInfo.dnaTop[scaf_right_ID].down == scaf_left_ID:
                middle[4] = '<'
                middle[6] = '<'
                middle[7] = '<'

            # Left staple nick or crossover?
            if dnaInfo.dnaTop[stap_left_ID].up == -1:
                middle[2] = '.'
            elif dnaInfo.dnaTop[stap_left_ID].down == -1:
                middle[2] = 'v'
            elif dnaInfo.dnaTop[stap_left_ID].up == stap_right_ID:
                middle[2:3+1] = '+<'
            elif dnaInfo.dnaTop[stap_left_ID].down == stap_right_ID:
                middle[2:3+1] = '+>'

            # Right staple nick or crossover?
            if dnaInfo.dnaTop[stap_right_ID].up == -1:
                middle[5] = '.'
            elif dnaInfo.dnaTop[stap_right_ID].down == -1:
                middle[5] = '^'
            elif dnaInfo.dnaTop[stap_right_ID].up == stap_left_ID:
                middle[4:5+1] = '>+'
            elif dnaInfo.dnaTop[stap_right_ID].down == stap_left_ID:
                middle[4:5+1] = '<+'

            middle = ''.join(middle)

            # Left scaffold nucleotide
            fid.write('{:>6} {} {} {} {:<6}\r\n'.format(
                scaf_left_ID, dnaInfo.dnaTop[scaf_left_ID].seq,
                middle,
                dnaInfo.dnaTop[scaf_right_ID].seq, scaf_right_ID
            ))  # TODO: to give an exact match, `{:<6}` -> `{}`

        fid.write('{:>6} {} | |\r\n'.format(
            scaf_left[-1], dnaInfo.dnaTop[scaf_left[-1]].seq))
        fid.write("         5'3' 5'3' \r\n\r\n\r\n")

    fid.close()
