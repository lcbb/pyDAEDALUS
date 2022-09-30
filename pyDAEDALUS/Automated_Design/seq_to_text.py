from datetime import datetime


def seqtoText(scaf_to_edge, edges, dnaInfo, file_name, scaf_name, singleXOs,
              full_file_name, Aform=False):  # pragma: no cover
    """
    Generate text file to visualize each edge's sequences and nick/crossover
    information.  The result can be physically printed and used to visualize
    each edge's sequences and nick/crossover information.

    Parameters
    ----------
    scaf_to_edge : list
        Ex2 list, where each row corresponds to one edge, 1st column is duplex
        from low ID to high ID vertex, 2nd column is from high to low. Each
        element is a row vector containing the scaffold base IDs in order on
        that duplex.
    edges : numpy.ndarray
        Ex2 array where each row corresponds to one edge, denoting the vertices
        being connected. 1st column > 2nd column
    dnaInfo : Automated_Design.dna_info.DnaInfo
        Structure containing all spatial and routing info
    file_name : str
        name of structure
    scaf_name : str
        name of scaffold
    singleXOs : int
        1 if single crossover vertex staples should be used,
        0 if double crossover vertex staples should be used.
    full_file_name : str
        file_name with more information

    Returns
    -------
    None
    """

    # TODO: Add bit about 0-indexing!!

    fid = open(full_file_name, 'w')

    num_edges = len(edges)

    fid.write('Structure: {}\n'.format(file_name))
    fid.write('Scaffold: {}\n'.format(scaf_name))

    if singleXOs > 0:
        fid.write('Vertex staples: Single crossovers\n')
    else:
        fid.write('Vertex staples: Double crossovers\n')

    fid.write('Produced: \n'.format(str(datetime.now())))
    fid.write(
        'Key:    outer strands = scaffold (sequence and index shown)\n')
    fid.write('        inner strands = staples\n')
    fid.write('                    | = nucleotide\n')
    fid.write('                    . = 5'' end\n')
    fid.write('               v or ^ = 3'' end\n')
    fid.write('   +>>+   or   +<<+   = staple crossover\n')
    fid.write(' >> >> >> or << << << = scaffold crossover\n')
    fid.write('Not shown: polyT loops in staple strands between edges\n\n')

    for edge_ID in range(num_edges):

        edge_bgn = edges[edge_ID][1]
        edge_fin = edges[edge_ID][0]

        scaf_right = scaf_to_edge[edge_ID][0]
        scaf_left = scaf_to_edge[edge_ID][1]

        fid.write('Vertex {} (top) to {} (bottom), Edge {}\n\n'.format(
            edge_bgn, edge_fin, edge_ID))
        fid.write("         3'5' 3'5' \n")
        if Aform: # A-form
			# Invert scaf_left
            scaf_right = scaf_right[0:]
            scaf_left = scaf_left[-1::-1]
        else: # B-form
            fid.write('              | | {} {}\n'.format(
				dnaInfo.dnaTop[scaf_right[0]].seq, scaf_right[0]
			))

			# Cut off first of scaf_right and invert scaf_left
            scaf_right = scaf_right[1:]
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
            fid.write('{:>6} {} {} {} {:<6}\n'.format(
                scaf_left_ID, dnaInfo.dnaTop[scaf_left_ID].seq,
                middle,
                dnaInfo.dnaTop[scaf_right_ID].seq, scaf_right_ID
            ))  # TODO: to give an exact match, `{:<6}` -> `{}`

        if not Aform:
            fid.write('{:>6} {} | |\n'.format(
                scaf_left[-1], dnaInfo.dnaTop[scaf_left[-1]].seq))

        fid.write("         5'3' 5'3' \n\n\r\n")

    fid.close()
