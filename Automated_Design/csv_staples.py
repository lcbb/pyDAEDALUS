from os import path

from Automated_Design.constants import RESULTS_FOLDERNAME


def csv_staples(full_file_name, named_stap_seq_list):
    """
    # Convert staple list to .csv file
    # Inputs: full_file_name = string containing full name of file
    #         named_stap_seq_list = cell array with second column identical to
    #             stap_seq_list, first column contains string of staple name
    # Outputs: CSV file of named_stap_seq_list
    ###########################################################################
    # by Sakul Ratanalert, MIT, Bathe Lab, 2016
    #
    # Copyright 2016. Massachusetts Institute of Technology. Rights Reserved.
    # M.I.T. hereby makes following copyrightable material available to the
    # public under GNU General Public License, version 2 (GPL-2.0). A copy of
    # this license is available at https://opensource.org/licenses/GPL-2.0
    ###########################################################################
    """

    full_filename = path.join(RESULTS_FOLDERNAME,
                              'staples_' + full_file_name + '.csv')
    fid = open(full_filename, 'w')

    num_stap = len(named_stap_seq_list)

    for stap_ID in range(num_stap):
        fid.write("{},{}\r\n".format(
            named_stap_seq_list[stap_ID][0], named_stap_seq_list[stap_ID][1]))

    fid.close()
