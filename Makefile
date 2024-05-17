.PHONY: help clean clean-pyc clean-build lint test test-all coverage docs release sdist

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "lint - check style with flake8 and format with black"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "sdist - package"

clean: clean-build clean-pyc

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

lint:
	black vienna test

test:
	pytest

test-all:
	tox

coverage:
	coverage run --source vienna -m pytest
	coverage report -m
	coverage html
	open htmlcov/index.html || xdg-open htmlcov/index.html || start htmlcov/index.html

docs:
	rm -f docs/vienna.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ vienna
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html || xdg-open docs/_build/html/index.html || start docs/_build/html/index.html

release: clean
	python setup.py sdist upload
	python setup.py bdist_wheel upload

sdist: clean
	python setup.py sdist
	python setup.py bdist_wheel upload
	ls -l dist
