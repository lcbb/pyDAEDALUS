from os import path

import click

from Automated_Design.ply_to_input import ply_as_filename_to_input
from Automated_Design.DX_cage_design import DX_cage_design


default_fname_no_ply = path.join(
    'PLY_Files', '05_icosahedron')  # No '.ply' extension required


@click.command()
@click.option('--fname_no_ply', default=default_fname_no_ply,
              help='Ply filename to read from without the .ply')
@click.option('--min_len_nt', default=52,
              help='min length nt')
@click.option('--display_plots', default=False, is_flag=True,
              help="""Set to True if you want to .""")
def run_demo_from_command_line(fname_no_ply, min_len_nt, display_plots):
    run_demo(fname_no_ply, min_len_nt, display_plots)


def run_demo(fname_no_ply, min_len_nt, display_plots=False):

    [coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs] = ply_as_filename_to_input(fname_no_ply, min_len_nt)

    scaf_seq = [] # Using default scaffold sequence
    scaf_name = [] # Using default scaffold name
    full_file_name = DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name)

    if display_plots:
        from matplotlib import pyplot as plt
        plt.show()

if __name__ == '__main__':
    run_demo_from_command_line()
