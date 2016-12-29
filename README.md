[Daedalus](http://daedalus-dna-origami-portal.org/) - software to render nearly any 3D geometry as a scaffolded DNA origami nanoparticle.

Your 3D models can be input using the Polygon File Format (.ply).

[![Build Status](https://travis-ci.com/TheGrimmScientist/Daedalus.svg?token=ygEd7xu7GnyJQrcup1bE&branch=master)](https://travis-ci.com/TheGrimmScientist/Daedalus) [![Coverage Status](https://coveralls.io/repos/github/TheGrimmScientist/Daedalus/badge.svg?t=xRpjvC)](https://coveralls.io/github/TheGrimmScientist/Daedalus)

# Installation

1. Clone this repo
1. (optional) Set up a virtualenv.
1. `pip -r requirements` in this project's root.

# Usage

`python demo_daedalus.py` to run with everything set to default.
`python demo_daedalus.py --help` to see what options are available, including their default values.
`python demo_daedalus.py --display_plots` to run everything as default, and also display plots to your screen (plots will still be saved to disk as well).
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron` to run program on `05_icosahedron` file within folder `PLY_files`.  Note filename unless specified is relative to project root and the format has to match your OS (this example uses Linux filename structure).
`python demo_daedalus.py --fname_no_ply=PLY_Files/05_icosahedron --display_plots` to both run on given shape and to display plots to your screen.

Additionally, pass in the `--suppress_console_output` if you want to not have things printed to your console at runtime.

All files saved into Results includes both shape name and `min_len_nt` in the filename to avoid undesired filename collisions.

# Tests
Testing uses [`py.test`](http://docs.pytest.org/en/latest/usage.html).  Ways to run tests manually include:
 * Simply run `make test`.  This leverages `Makefile` to run tests with the right syntax.
 * Run all tests by directly calling pytest: `py.test`
 * Run all tests and drop into a [`bpython`](https://bpython-interpreter.org/)-enabled debugger at failure: `py.test --bpdb`
 * To run a test class or stand-alone function: `py.test tests/test_dx_cage_design.py::test_ply_input_for_05_icosahedron`
 * To run a single function within a test class: `py.test tests/test_dx_cage_design.py::TestIntegrationsUsing01Tetrahedron::test_generate_spanning_tree`
