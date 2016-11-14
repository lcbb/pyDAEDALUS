A port of DAEDALUS from Matlab to Python.


# Installation

1. clone this repo
1. (optional) Set up a virtualenv.
1. `pip -r requirements` in this project's root.

# Usage

`python demo_daedalus.py` to run with everything set to default.
`python demo_daedalus.py --help` to see what options are available, including their default values.
`python demo_daedalus.py --display_plots` to run everything as default, and also display plots to your screen (plots will still be saved to disk as well).
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron` to run program on `05_icosahedron` file within folder `PLY_files`.  Note filename unless specified is relative to project root and the format has to match your OS (this example uses Linux filename structure).
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron --display_plots` to both run on given shape and to display plots to your screen.


All files saved into Results includes both shape name and `min_len_nt` in the filename to avoid undesired filename collisions.

# Tests
Currently using Python's built in `unittest` package.  All tests import into `test.py`, so calling the
 * Simply run the tests file to run all tests: `python test.py`
 * To run a test class, specify the test class name as the first arg to `test.py`: `python test.py TestIntegrationsUsing01Tetrahedron`
 * To run a single test, additionally specify the function name like this: `python test.py TestIntegrationsUsing01Tetrahedron.test_split_edge`
