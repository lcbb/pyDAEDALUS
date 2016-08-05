import shutil
from os import listdir, path, mkdir

from tqdm import tqdm

from demo_daedalus import run_demo


def grab_all_ply_filenames_from_directory(directory):
    filenames = sorted(f for f in listdir(directory) if path.isfile(path.join(directory, f)))
    filenames = [path.join(directory, filename) for filename in filenames]
    return filenames

def create_empty_directory(directory):
    if path.exists(directory):
        shutil.rmtree(directory)
    mkdir(directory)

if __name__ == '__main__':
    ply_directory = 'PLY_Files'
    out_directory = 'Drawings'
    ply_filenames = grab_all_ply_filenames_from_directory(ply_directory)

    print ply_filenames[:10]

    create_empty_directory(out_directory)

    for filename in tqdm([f for f in ply_filenames if f[-4:] == '.ply']):
        fname_no_ply = filename[:-4]
        run_demo(fname_no_ply, results_foldername=out_directory)
