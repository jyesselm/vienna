.PHONY: help clean clean-pyc clean-build lint test test-all coverage docs release sdist install-dev

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
	@echo "install-dev - install development dependencies"

install-dev:
	pip install -e ".[dev]"

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
	flake8 vienna test --max-line-length=88 --extend-ignore=E203,W503

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
	pip install -q build twine
	python -m build
	twine upload dist/*

sdist: clean
	pip install -q build
	python -m build
	ls -l dist
