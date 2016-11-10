from os import path

import click

from Automated_Design.ply_to_input import ply_as_filename_to_input
from Automated_Design.DX_cage_design import DX_cage_design


default_fname_no_ply = path.join('PLY_Files',
                         '05_icosahedron')  # No '.ply' extension required
default_fname_no_ply = path.join('PLY_Files', '01_tetrahedron')


@click.command()
@click.option('--fname_no_ply', default=default_fname_no_ply,
              help='Ply filename to read from without the .ply')
@click.option('--min_len_nt', default=52,
              help='min length nt')
@click.option('--results_foldername', default=None,
              help="""Define this if you want to save plots to a results folder
                      rather than display them immediately to screen.""")
def run_demo_from_command_line(fname_no_ply, min_len_nt, results_foldername):
    run_demo(fname_no_ply, min_len_nt, results_foldername)


def run_demo(fname_no_ply, min_len_nt, results_foldername):

    if results_foldername:
        if not path.isdir(results_foldername):
            raise click.ClickException(
                "`--results_foldername` has to already exist and be a folder.")

    [coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs] = ply_as_filename_to_input(fname_no_ply, min_len_nt, results_foldername=results_foldername)

    scaf_seq = [] # Using default scaffold sequence
    scaf_name = [] # Using default scaffold name
    full_file_name = DX_cage_design(coordinates, edges, faces, edge_length_vec, file_name, staple_name, singleXOs, scaf_seq, scaf_name, results_foldername=results_foldername)


if __name__ == '__main__':
    run_demo_from_command_line()
