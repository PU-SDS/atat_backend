# Usage:
# make          # setup packages and install required dependencies
# make setup    # install poetry and pre-commit hooks
# make update   # upgrade poetry dependencies
# make upgrade  # upgrade pre-commit hooks
# make clean    # remove ALL Clean out cached pre-commit files
# make tidy     # run pre-commit hooks
# make test		# run all tests using PyTest

ifeq ($(OS),Windows_NT)
	SHELL := cmd
	POETRY_HOME := $(APPDATA)\Python\Scripts
	PATH := $(POETRY_HOME)\bin:$(PATH)
	PYTHON := python
	export PATH
	export POETRY_HOME
else
	SHELL := /bin/bash
	POETRY_HOME := $(HOME)/.poetry
	PATH := $(POETRY_HOME)/bin:$(PATH)
	PYTHON := python3
	export PATH
	export POETRY_HOME
endif

.PHONY: all setup install update upgrade clean tidy
all: setup install resources update tidy

setup:
	@curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | ${PYTHON} - -p
	@docker-compose up -d

install:
	@poetry install
	@poetry run pre-commit install -f

resources:
	@poetry run setup

update:
	@poetry self update
	@poetry update

upgrade:
	@poetry run pre-commit autoupdate

clean:
	@poetry cache clear pypi --all
	@poetry run pre-commit clean
	@poetry run pre-commit gc

tidy:
	@poetry run pre-commit run --all-files

test:
	@pytest
