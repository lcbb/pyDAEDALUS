[Daedalus](http://daedalus-dna-origami-portal.org/) - software to render nearly any 3D geometry as a scaffolded DNA origami nanoparticle.

Your 3D models can be input using the Polygon File Format (.ply).

[![Build Status](https://travis-ci.org/TheGrimmScientist/pyDAEDALUS.svg?branch=master)](https://travis-ci.org/TheGrimmScientist/pyDAEDALUS) [![Coverage Status](https://coveralls.io/repos/github/TheGrimmScientist/pyDAEDALUS/badge.svg?branch=master)](https://coveralls.io/github/TheGrimmScientist/pyDAEDALUS?branch=master) [![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://raw.githubusercontent.com/lcbb/pyDAEDALUS/master/LICENSE.txt)

# Installation

1. Clone this repo
1. (optional) Set up a virtualenv.
1. `pip install -r requirements` in this project's root.
1. (optional) Also run `pip install -r dev_requirements` to give you the libraries used for testing / code development.

# Usage
At the minimum when calling `demo_daedalus.py`, you need to tell it what ply file to run.  For example, running pyDAEDALUS on the included ply file `05_icosahedron.ply` can be done with the command:
`python demo_daedalus.py --input_filename=PLY_Files/05_icosahedron.ply`
 Note including the `.ply` extension is optional.  So the following command would be equivalent: `python demo_daedalus.py --input_filename=PLY_Files/05_icosahedron`

You can also run a batch of ply files at a time by specifying an input foldername rather than a single ply file:
`python demo_daedalus.py --input_foldername=PLY_Files`

Adding the flag `--display_plots` asks the program, in addition to saving the plots to your output folder, to print them to screen.  For example:
`python demo_daedalus.py --input_filename=PLY_Files/37_enneagonal_trapezohedron.ply --display_plots`

Additional options include:
 * `--results_foldername=yourfoldername` where 'yourfoldername' is the name of the folder you want to save this command's results in.  Folder will be created if it doesn't already exist.
 * `--suppress_console_output` (flag) if you want to not have things printed to your console at runtime. (always suppressed when running a batch)
 * `--reset_results_folder` (flag) if you want to erase all data within your output folder before running the current command.

You can find the complete list of arguments with descriptions with `python demo_daedalus.py --help`.


All files saved into Results includes both shape name and `min_len_nt` in the filename to avoid undesired filename collisions.

# Tests
Be sure to have installed dev_requirements.txt if you plan on running tests.  Or at the very least, install `pytest`.
Testing uses [`py.test`](http://docs.pytest.org/en/latest/usage.html).  Ways to run tests manually include:
 * Simply run `make test`.  This leverages `Makefile` to run tests with the right syntax.
 * Run all tests by directly calling pytest: `py.test`
 * Run all tests and drop into a [`bpython`](https://bpython-interpreter.org/)-enabled debugger at failure: `py.test --bpdb`
 * To run a test class or stand-alone function: `py.test tests/test_ply_to_input.py::test_ply_input_for_05_icosahedron`
 * To run a single function within a test class: `py.test tests/test_ply_to_input.py::TestIntegrationsUsing01Tetrahedron::test_generate_spanning_tree`
 
## Coverage

 Testing coverage done with `pytest-cov` package.  Manually run coverage report with:
 * `py.test --cov=Automated_Design tests/` to run coverage on folder `Automated_Design` using all tests found in `tests` folder.
 * `py.test --cov-report html --cov-report term --cov=Automated_Design tests/` to run coverage report as above, but with also generating a line-by-line report formatted as html (viewable in your browser by opening the created files)

The only code deliberately excluded from tests are plotting functions (e.g., `ply_to_input.plot_edge_length_distributions`) and functions that simply write raw internal data to disk (e.g., `csv_staples`).  All such functions are marked to not be included in the coverage report by adding the comment `  # pragma: no cover` on the same line as the function definition.

## Makefile

For convenience, `Makefile` is used to add several command shortcuts to do common testing operations:
* `make` to run both linting and tests (linting is only ran if tests pass)
* `make test` to run only tests.
* `make lint` to run only linting.
* `make coverage` to run coverage report.  The report is returned as both top level statistics to terminal and a line-by-line report saved into folder `htmlcov` as html files.