A port of DAEDALUS from Matlab to Python.


# Installation

1. clone this repo
1. (optional) Set up a virtualenv.
1. `pip -r requirements` in this project's root.

# Usage

`python demo_daedalus.py` to run with everything set to default.
`python demo_daedalus.py --help` to see what options are available, including their default values.
`python demo_daedalus.py --results_foldername=Results` to run program and save all plots to `Results` folder rather than print them to screen.  Note foldername is relative to this project's root and has to already exist.
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron` to run program on `05_icosahedron` file within folder `PLY_files`.  Note filename unless specified relative to project root and the format has to match your OS (this example uses Linux filename structure).
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron --results_foldername=Results` to both run on given shape and save plots into specified folder.


Saving to Results includes shape name and `min_len_nt` in all output files so you can run multiple shapes and not have filename collisions.

# Tests
Currently using Python's built in `unittest` package.  All tests import into `test.py`, so calling the
 * Simply run the tests file to run all tests: `python test.py`
 * To run a test class, specify the test class name as the first arg to `test.py`: `python test.py TestIntegrationsUsing01Tetrahedron`
 * To run a single test, additionally specify the function name like this: `python test.py TestIntegrationsUsing01Tetrahedron.test_split_edge`
