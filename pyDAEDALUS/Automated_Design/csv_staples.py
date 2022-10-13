from os import path


def csv_staples(full_file_name, named_stap_seq_list,
                results_foldername):  # pragma: no cover
    """
    Convert staple list to .csv file.  File is saved into `results_foldername`.

    Parameters
    ----------
    full_file_name : str
        String containing full name of file
    named_stap_seq_list : list
        List with second column identical to stap_seq_list, first column
        contains string of staple name.
    results_foldername : str
        Foldername to save this file into.

    Returns
    -------
    None
    """

    full_filename = path.join(results_foldername,
                              'staples_' + full_file_name + '.csv')
    fid = open(full_filename, 'w')

    num_stap = len(named_stap_seq_list)

    for stap_ID in range(num_stap):
        fid.write("{},{}\n".format(
            named_stap_seq_list[stap_ID][0], named_stap_seq_list[stap_ID][1]))

    fid.close()
