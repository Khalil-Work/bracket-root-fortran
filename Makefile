# Bracket Root - native Fortran safeguarded scalar solvers
#
#   make            build the test binary (default)
#   make test       build and run the assertion test suite
#   make run        build and run the comparison driver -> build/csv/*.csv
#   make results    run the driver, then Python figures/GIF/summary
#   make clean      remove build artifacts

FC = gfortran
FFLAGS = -std=f2018 -Wall -Wextra -O2 -fall-intrinsics

OBJDIR = build

ifeq ($(OS),Windows_NT)
  EXE := .exe
  MKDIR = cmd /c "if not exist build mkdir build"
  MKCSV = cmd /c "if not exist build\csv mkdir build\csv"
  CLEAN = cmd /c "if exist build rmdir /s /q build"
else
  EXE :=
  MKDIR = mkdir -p build
  MKCSV = mkdir -p build/csv
  CLEAN = rm -rf build
endif

PYTHON ?= python

.PHONY: all test run results clean
.DEFAULT_GOAL := all

all: $(OBJDIR)/test_roots$(EXE)

$(OBJDIR):
	@$(MKDIR)

$(OBJDIR)/csv:
	@$(MKCSV)

$(OBJDIR)/roots.o: src/roots.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -J$(OBJDIR) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/test_roots.o: tests/test_roots.f90 $(OBJDIR)/roots.o | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/test_roots$(EXE): $(OBJDIR)/test_roots.o $(OBJDIR)/roots.o
	$(FC) $(FFLAGS) -o $@ $(OBJDIR)/test_roots.o $(OBJDIR)/roots.o

test: all
	$(OBJDIR)/test_roots$(EXE)

run: $(OBJDIR)/compare$(EXE) | $(OBJDIR)/csv
	$(OBJDIR)/compare$(EXE)

$(OBJDIR)/compare.o: examples/compare.f90 $(OBJDIR)/roots.o | $(OBJDIR)
	$(FC) $(FFLAGS) -I$(OBJDIR) -c $< -o $@

$(OBJDIR)/compare$(EXE): $(OBJDIR)/compare.o $(OBJDIR)/roots.o
	$(FC) $(FFLAGS) -o $@ $(OBJDIR)/compare.o $(OBJDIR)/roots.o

results: run
	$(PYTHON) tools/make_results.py

clean:
	@$(CLEAN)
