[Daedalus](http://daedalus-dna-origami-portal.org/) - software to render nearly any 3D geometry as a scaffolded DNA origami nanoparticle.

Your 3D models can be input using the Polygon File Format (.ply).

[![Build Status](https://travis-ci.com/TheGrimmScientist/Daedalus.svg?token=ygEd7xu7GnyJQrcup1bE&branch=master)](https://travis-ci.com/TheGrimmScientist/Daedalus) [![Coverage Status](https://coveralls.io/repos/github/TheGrimmScientist/Daedalus/badge.svg?t=xRpjvC)](https://coveralls.io/github/TheGrimmScientist/Daedalus)

# Installation

1. Clone this repo
1. (optional) Set up a virtualenv.
1. `pip -r requirements` in this project's root.

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
Testing uses [`py.test`](http://docs.pytest.org/en/latest/usage.html).  Ways to run tests manually include:
 * Simply run `make test`.  This leverages `Makefile` to run tests with the right syntax.
 * Run all tests by directly calling pytest: `py.test`
 * Run all tests and drop into a [`bpython`](https://bpython-interpreter.org/)-enabled debugger at failure: `py.test --bpdb`
 * To run a test class or stand-alone function: `py.test tests/test_dx_cage_design.py::test_ply_input_for_05_icosahedron`
 * To run a single function within a test class: `py.test tests/test_dx_cage_design.py::TestIntegrationsUsing01Tetrahedron::test_generate_spanning_tree`
