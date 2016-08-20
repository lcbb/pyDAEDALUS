A port of DAEDALUS from Matlab to Python.


# Installation

1. clone this repo
1. (optional) Set up a virtualenv.
1. `pip -r requirements` in this project's root.

# Tests
Currently using Python's built in `unittest` package.  All tests import into `test.py`, so calling the
 * Simply run the tests file to run all tests: `python test.py`
 * To run a test class, specify the test class name as the first arg to `test.py`: `python test.py TestIntegrationsUsing01Tetrahedron`
 * To run a single test, additionally specify the function name like this: `python test.py TestIntegrationsUsing01Tetrahedron.test_split_edge`
