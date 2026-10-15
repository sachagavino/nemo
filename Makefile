#==============================================================================
# NEMO dynamical gas-grain code -- skeleton extracted from NMGC-2.0
#
# Build:   make            (default, -O2)
#          make OPT=-O0    (reference-equivalence runs, matches nmgc-2.0's -O0)
#          make debug      (bounds checking, warnings, backtrace)
#          make clean
#
# The binary is written to bin/nmgc, objects and .mod files to build/.
#==============================================================================

FC      := gfortran
SRCDIR  := src
ODEDIR  := $(SRCDIR)/odepack
BUILD   := build
BINDIR  := bin
EXE     := $(BINDIR)/nmgc

# STAGE 4 note: nmgc-2.0 shipped with FFLAGS empty, i.e. -O0. The default here
# is -O2. Use OPT=-O0 to reproduce the reference timings bit for bit.
OPT     ?= -O2
FFLAGS  := $(OPT) -J$(BUILD)
DBGFLAGS := -O0 -g -fbacktrace -fcheck=bounds -Wall -Wextra -J$(BUILD)

# The old ODEPACK sources have type mismatches that gfortran >= 10 rejects.
GFORTRAN_MAJOR := $(shell $(FC) -dumpversion | cut -d. -f1)
ifeq ($(shell [ $(GFORTRAN_MAJOR) -ge 10 ] && echo yes),yes)
    MISMATCH := -fallow-argument-mismatch
else
    MISMATCH :=
endif
ODEFLAGS := $(OPT) -w $(MISMATCH) -J$(BUILD)

# Compilation order = module dependency order. Keep it explicit: with ~10 files
# a generated dependency graph buys nothing and hides the layering.
MODULES := \
  $(SRCDIR)/iso_fortran_env.f90 \
  $(SRCDIR)/numerical_types.f90 \
  $(SRCDIR)/utilities.f90 \
  $(SRCDIR)/global_variables.f90 \
  $(SRCDIR)/shielding.f90 \
  $(SRCDIR)/environment.f90 \
  $(SRCDIR)/dust/dust_grid.f90 \
  $(SRCDIR)/input_output.f90 \
  $(SRCDIR)/ode_solver.f90 \
  $(SRCDIR)/dust/dustevolution.f90 \
  $(SRCDIR)/gasgrain.f90 \
  $(SRCDIR)/outputs.f90

MAIN    := $(SRCDIR)/main.f90
ODESRC  := $(ODEDIR)/opkda1.f90 $(ODEDIR)/opkda2.f90 $(ODEDIR)/opkdmain.f90

MODOBJ  := $(patsubst $(SRCDIR)/%.f90,$(BUILD)/%.o,$(MODULES))
MAINOBJ := $(BUILD)/main.o
ODEOBJ  := $(patsubst $(ODEDIR)/%.f90,$(BUILD)/%.o,$(ODESRC))

.PHONY: all debug clean
all: $(EXE)

debug: FFLAGS := $(DBGFLAGS)
debug: clean $(EXE)

$(EXE): $(MODOBJ) $(MAINOBJ) $(ODEOBJ) | $(BINDIR)
	$(FC) $(FFLAGS) $(MODOBJ) $(MAINOBJ) $(ODEOBJ) -o $@

# Modules are compiled strictly in the order given by MODULES: each one is
# forced to wait for the previous ones through the .mod files in $(BUILD).
# $(dir $@) lets objects live in nested build/ subdirs (e.g. build/dust/).
$(BUILD)/%.o: $(SRCDIR)/%.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD)/%.o: $(ODEDIR)/%.f90 | $(BUILD)
	@mkdir -p $(dir $@)
	$(FC) -c $(ODEFLAGS) $< -o $@

$(BUILD) $(BINDIR):
	mkdir -p $@

# --- explicit module ordering (module M must exist before its users compile)
$(BUILD)/utilities.o:        $(BUILD)/iso_fortran_env.o $(BUILD)/numerical_types.o
$(BUILD)/global_variables.o: $(BUILD)/utilities.o
$(BUILD)/shielding.o:        $(BUILD)/numerical_types.o
$(BUILD)/environment.o:      $(BUILD)/global_variables.o
$(BUILD)/dust/dust_grid.o:   $(BUILD)/global_variables.o
$(BUILD)/input_output.o:     $(BUILD)/global_variables.o $(BUILD)/dust/dust_grid.o $(BUILD)/environment.o
$(BUILD)/ode_solver.o:       $(BUILD)/global_variables.o $(BUILD)/shielding.o
$(BUILD)/gasgrain.o:         $(BUILD)/input_output.o $(BUILD)/ode_solver.o $(BUILD)/environment.o
$(BUILD)/outputs.o:          $(BUILD)/gasgrain.o
$(BUILD)/main.o:             $(BUILD)/gasgrain.o $(BUILD)/outputs.o

clean:
	@rm -rf $(BUILD) $(BINDIR)
	@echo "build/ and bin/ removed."
