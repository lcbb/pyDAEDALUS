all: test lint
.PHONY: all test lint coverage

test:
	py.test

lint:
	flake8

coverage:
	py.test --cov-report html --cov-report term --cov=Automated_Design tests/

bpdb:
	py.test --bpdb