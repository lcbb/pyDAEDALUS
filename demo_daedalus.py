import shutil
from os import listdir, path, makedirs

import click
from tqdm import tqdm

from Automated_Design.ply_to_input import ply_to_input
from Automated_Design.DX_cage_design import DX_cage_design


def grab_all_ply_filenames_from_directory(directory):
    filenames = sorted(f for f in listdir(directory)
                       if path.isfile(path.join(directory, f)))
    filenames = [path.join(directory, filename) for filename in filenames]
    return filenames


def create_directory(directory, reset=True):
    if reset:
        if path.exists(directory):
            shutil.rmtree(directory)
    if not path.exists(directory):
        makedirs(directory)


@click.command()
@click.option('--input_filename', default=None,
              help='Ply filename to read from without the .ply.  Please only'
                   'supply a filename or a folder name')
@click.option('--input_foldername', default=None,
              help='Foldername containing PLY files to be rendered.  Please '
                   'only supply a folder name or a filename.')
@click.option('--results_foldername', default='Results',
              help='Foldername to save all results.')
@click.option('--reset_results_folder', default=False, is_flag=True,
              help="""Set to True if you want to clear past data (if any) from
              your specified results folder""")
@click.option('--min_len_nt', default=52,
              help='min length nt')
@click.option('--display_plots', default=False, is_flag=True,
              help="""Set to True if you want to also display plots to
              screen.""")
@click.option('--suppress_console_output', default=False, is_flag=True,
              help="""Set to False is you want to disable logging to console
              for a single file.  This option is ignored and console output is
              always suppressed when converting a batch.
              """)
def run_demo_from_command_line(input_filename, input_foldername,
                               results_foldername, reset_results_folder,
                               min_len_nt, display_plots,
                               suppress_console_output):

    # make sure input parameters make sense:
    if input_filename and input_foldername:
        # Both an input filename and an input foldername were specified.
        raise Exception("Please only supply a foldername or a filename to "
                        "specify your input. Supplying both is ambiguous.")
    if (input_filename is None) and (input_foldername is None):
        # Neither an input filename or an input foldername were specified.
        raise Exception("Please supply a foldername or a filename to specify"
                        "your input.")

    # Prepare system:
    print_to_console = not suppress_console_output
    # Create `Results` directory if it doesn't already exist.  Optionally
    # erase data that was in it.
    create_directory(results_foldername, reset=reset_results_folder)

    # Either run as batch or run single file.
    if input_foldername:
        # When running as batch, don't print to console (it gets messy).
        print_to_console = False
        run_batch(input_foldername, min_len_nt, display_plots,
                  print_to_console, results_foldername)
    elif input_filename:
        run_single_file(input_filename, min_len_nt, results_foldername,
                        display_plots, print_to_console)

    else:
        # It was already asserted only one of (foldername, filename) exist.
        raise Exception("Program should never reach here.  Please file a bug "
                        "report including a copy of the paramters you "
                        "entered.")


def run_single_file(input_filename, min_len_nt, results_foldername,
                    display_plots=False, print_to_console=True):

    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_to_input(
            input_filename, results_foldername, min_len_nt)

    scaf_seq = []  # Using default scaffold sequence
    scaf_name = []  # Using default scaffold name
    full_file_name = DX_cage_design(  # noqa: F841
        coordinates, edges, faces, edge_length_vec, file_name,
        staple_name, singleXOs, scaf_seq, scaf_name, results_foldername,
        print_to_console=print_to_console)

    if display_plots:
        from matplotlib import pyplot as plt
        plt.show()


def run_batch(input_foldername, min_len_nt, display_plots, print_to_console,
              results_foldername):
    ply_filenames = grab_all_ply_filenames_from_directory(input_foldername)

    assert len(ply_filenames) > 0, "No Ply files found in given foldername."

    for input_filename in tqdm([f for f in ply_filenames if f[-4:] == '.ply']):
        run_single_file(input_filename, min_len_nt, results_foldername,
                        display_plots=display_plots,
                        print_to_console=print_to_console)


if __name__ == '__main__':
    run_demo_from_command_line()
