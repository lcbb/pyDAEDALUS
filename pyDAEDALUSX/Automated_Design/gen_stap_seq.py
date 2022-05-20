import numpy as np


def gen_stap_seq(staples, scaf_seq, staple_name, scaf_name,
                 len_scaf_used):
    """
    Generate staple sequences from scaffold sequence.

    Parameters
    ----------
    staples : list
        List with E rows, each row contains another list.  Some sub-lists may
        be empty as fragments are combined into full staples.  Columns 0-5
        contain vertex staples, while 6+ contain edge staples.
    scaf_seq : str
        String of scaffold, length may be longer than
       required. First nucleotide will be placed at index 1
    short_name : str
        string to name staples, may be same as file_name
    scaf_name : str
        name of scaffold (string)
    len_scaf_used : int
        length of scaffold used, = 2*sum(edge_length_vec)

    Returns
    -------
    stap_seq
        List same dimensions as staples, each cell contains string of
        staple's sequence
    stap_seq_list
        stap_seq arranged in a 1-column list
    stap_list
        staples arranged in the same 1-column list
    named_stap_seq_list
        List with second column identical to stap_seq_list, first column
        contains string of staple name
    """

    # # Initialize outputs
    num_edges = len(staples)

    stap_seq = []
    stap_seq_list = []
    stap_list = []
    named_stap_seq_list = []

    for edge_ID in range(num_edges):
        staps_for_this_edge = []
        for stap_ID in range(len(staples[edge_ID])):
            stap = staples[edge_ID][stap_ID]  # obtain staple index information
            seq = u''  # initialize string
            for nt_ID in stap:
                # # because A binds to T, and G binds to C...
                if (nt_ID is None) or (scaf_seq[nt_ID] == 'A'):
                    seq += 'T'
                elif scaf_seq[nt_ID] == 'T':
                    seq += 'A'
                elif scaf_seq[nt_ID] == 'G':
                    seq += 'C'
                elif scaf_seq[nt_ID] == 'C':
                    seq += 'G'
                elif scaf_seq[nt_ID] == 'U': # in case user inputs RNA seq
                    seq += 'A'
                else:
                    # Error -- unrecognized nucleotide in scaffold sequence
                    assert False, 'Unrecognized nucleotide in scaffold sequence'
                
            staps_for_this_edge.append(seq)

            if seq:  # if seq exists (may have been blank)
                stap_seq_list.append(seq)  # store seq in stap_seq_list

                # # Store seq in named_stap_seq_list. First number is Edge ID,
                # # Second is 5' end of staple. Third indicates vertex or edge
                # # staple.
                if stap_ID <= 3:  # Vertex staple
                    named_stap_seq_a = staple_name + '_{}-{}-V'.format(
                        edge_ID, stap[0])
                else:  # Edge staple
                    named_stap_seq_a = staple_name + '_{}-{}-E'.format(
                        edge_ID, stap[0])
                named_stap_seq_b = seq
                named_stap_seq_list.append(
                    [named_stap_seq_a, named_stap_seq_b])

                stap_list.append(stap)  # store stap in stap_list
        stap_seq.append(staps_for_this_edge)

    named_stap_seq_list.append([None, None])

    # # Add sequence of scaffold at end
    block_size = 10000
    if len_scaf_used > block_size:  # now for readabiliy, was for excel limits
        num_mults = int(np.floor(float(len_scaf_used) / block_size))
        # leftover = len_scaf_used % 10000
        for i in range(num_mults):
            start_i = i * block_size
            end_i = (i + 1) * block_size
            named_stap_seq_a = scaf_name + '_({}-{})'.format(start_i, end_i-1)
            named_stap_seq_b = scaf_seq[start_i:end_i]
            named_stap_seq_list.append([named_stap_seq_a, named_stap_seq_b])

        start_i = num_mults * block_size
        end_i = len_scaf_used
        named_stap_seq_a = scaf_name + '_({}-{})'.format(start_i, end_i-1)
        named_stap_seq_b = scaf_seq[start_i:end_i]
        named_stap_seq_list.append([named_stap_seq_a, named_stap_seq_b])

    else:
        named_stap_seq_a = scaf_name + '_(0-{})'.format(len_scaf_used-1)
        named_stap_seq_b = scaf_seq[:len_scaf_used]
        named_stap_seq_list.append([named_stap_seq_a, named_stap_seq_b])

    return [stap_seq, stap_seq_list, stap_list, named_stap_seq_list]
