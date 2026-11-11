# rs25-gone-fluorine: Space Shuttle Main Engine (RS-25) performance model
#
#   make setup      create .venv and install requirements.txt
#   make            run the model (SSMEMain.py)
#   make test       run the regression checks
#   make figures    redraw the README graphics in docs/
#   make clean      remove __pycache__
#   make distclean  also remove .venv

PYTHON ?= python3
VENV   ?= .venv
PY     := $(abspath $(VENV))/bin/python

.DEFAULT_GOAL := run
.PHONY: help setup run test figures clean distclean

help:
	@sed -n '3,8p' Makefile | sed 's/^# \{0,1\}//'

$(PY):
	$(PYTHON) -m venv $(VENV)
	$(PY) -m pip install --upgrade pip
	$(PY) -m pip install -r requirements.txt

setup: $(PY)

run: $(PY)
	$(PY) SSMEMain.py

test: $(PY)
	$(PY) tests/check_regression.py

figures: $(PY)
	$(PY) docs/make_figures.py

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

distclean: clean
	rm -rf $(VENV)
