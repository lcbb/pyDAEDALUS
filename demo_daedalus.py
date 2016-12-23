from os import path, makedirs

import click

from Automated_Design.constants import RESULTS_FOLDERNAME
from Automated_Design.ply_to_input import ply_as_filename_to_input
from Automated_Design.DX_cage_design import DX_cage_design

# No '.ply' extension required
default_fname_no_ply = path.join('PLY_Files', '05_icosahedron')


@click.command()
@click.option('--fname_no_ply', default=default_fname_no_ply,
              help='Ply filename to read from without the .ply')
@click.option('--min_len_nt', default=52,
              help='min length nt')
@click.option('--display_plots', default=False, is_flag=True,
              help="""Set to True if you want to also display plots to
              screen.""")
@click.option('--suppress_console_output', default=False, is_flag=True,
              help="""Set to False is you want to disable logging to console.
              """)
def run_demo_from_command_line(fname_no_ply, min_len_nt,
                               display_plots, suppress_console_output):
    print_to_console = not suppress_console_output
    run_demo(fname_no_ply, min_len_nt, display_plots, print_to_console)


def run_demo(fname_no_ply, min_len_nt,
             display_plots=False, print_to_console=True):

    # Create `Results` directory if it doesn't already exist
    if not path.exists(RESULTS_FOLDERNAME):
        makedirs(RESULTS_FOLDERNAME)

    coordinates, edges, faces, edge_length_vec, file_name, \
        staple_name, singleXOs = ply_as_filename_to_input(
            fname_no_ply, min_len_nt)

    scaf_seq = []  # Using default scaffold sequence
    scaf_name = []  # Using default scaffold name
    full_file_name = DX_cage_design(  # noqa: F841
        coordinates, edges, faces, edge_length_vec, file_name,
        staple_name, singleXOs, scaf_seq, scaf_name,
        print_to_console=print_to_console)

    if display_plots:
        from matplotlib import pyplot as plt
        plt.show()


if __name__ == '__main__':
    run_demo_from_command_line()
