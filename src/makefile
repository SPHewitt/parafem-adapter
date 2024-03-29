# -----------------
# Compiler
# -----------------
#FC=mpif90
#FCFLAGS=-Wall -g 
include ../build/$(MACHINE).inc


# -----------------
# Directories
# -----------------
SRCDIR=../src
OBJDIR=../obj
MODDIR=../mod
BINDIR=../bin

# -----------------
# Libraries
# -----------------
PRECICEROOT=${PRECICE_ROOT}
PRECICELIB=$(PRECICEROOT)/build/last

PARAFEMROOT=${PARAFEM_HOME}
PARAFEMLIB=$(PARAFEMROOT)/lib
PARAFEMMOD=$(PARAFEMROOT)/include/mpi

# Executable name
EXE=parafem

SRC= \
	$(SRCDIR)/SolverInterfaceF2003.f90 \
    $(SRCDIR)/utils.f90 \
    $(SRCDIR)/run_nonlinear.f90 \
	$(SRCDIR)/run_nonlinear_bd.f90 \
    $(SRCDIR)/input_precice.f90 \
    $(SRCDIR)/parafem.f90


OBJ:=$(SRC:$(SRCDIR)%.f90=$(OBJDIR)%.o)

# -----------------
# Rules
# -----------------
#
all: check-env release

release: check-env
	@echo $(SRC)
	@echo $(OBJ)
	$(MAKE) $(EXE) 

debug: check-env
	@echo $(SRC)
	@echo $(OBJ)
	$(MAKE) $(EXE) 

$(EXE): $(OBJ)
	@echo Compiling main
	$(FC) $(FCFLAGS) -I$(PARAFEMMOD) -I$(MODDIR) -o $(BINDIR)/$@ $^ -L$(PRECICELIB) -lprecice -L$(PARAFEMLIB) -lParaFEM_mpi.5.0.3 -larpack_linuxdesktop

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@echo Compiling "$@"
	$(FC) $(FCFLAGS) -I$(PARAFEMMOD) -o $@ -c $< -J$(MODDIR) 

.PHONY: check-env
check-env:
ifndef PARAFEM_HOME
	$(error PARAFEM_HOME is undefined)
endif

.PHONY: clean
clean:
	@echo cleaning
	rm -rf ../obj/*.o
	rm -rf ../bin/*
	rm -rf ../mod/*
	rm -rf ../*.log

