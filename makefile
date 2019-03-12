# -----------------
# Compiler
# -----------------
FC=/usr/local/bin/mpif90
FCFLAGS=
FFLAGS=

# -----------------
# Directories
# -----------------
SRCDIR=src
OBJDIR=obj
MODDIR=mod
BINDIR=bin

# -----------------
# Libraries
# -----------------
PRECICEROOT=${PRECICE_ROOT}
PRECICELIB=$(PRECICEROOT)/build/last

PARAFEMROOT=${PARAFEM_HOME}
PARAFEMLIB=$(PARAFEMROOT)/lib
PARAFEMMOD=$(PARAFEMROOT)/include/mpi

# Executable and main file
EXE=parafem

SRC= \
	$(SRCDIR)/SolverInterfaceF2003.f90 \
    $(SRCDIR)/parafeml.f90 \
    $(SRCDIR)/parafemutils.f90 \
    $(SRCDIR)/parafemnl.f90

OBJ:=$(SRC:$(SRCDIR)%.f90=$(OBJDIR)%.o)

# -----------------
# Rules
# -----------------
all: check-env
	@echo $(SRC)
	@echo $(OBJ)
	$(MAKE) $(EXE) 

$(EXE): $(OBJ)
	@echo Compiling main
	$(FC) $(FCFLAGS) -I$(PARAFEMMOD) $(EXE).f90 -o $(BINDIR)/$@ $^ -L$(PRECICELIB) -lprecice -L$(PARAFEMLIB) -lParaFEM_mpi.5.0.3 -larpack_linuxdesktop

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@echo Compiling "$@"
	$(FC) $(FCFLAGS) -I$(PARAFEMMOD) -o $@ -c $< -J$(MODDIR) -L$(PARAFEMLIB) -lParaFEM_mpi.5.0.3 -larpack_linuxdesktop

.PHONY: check-env
check-env:
ifndef PARAFEM_HOME
	$(error PARAFEM_HOME is undefined)
endif

.PHONY: clean
clean:
	@echo cleaning
	rm -rf obj/*.o
	rm -rf main
	rm -rf *.log

.PHONY: test
test:
	./bin/parafem precice-config.xml SolverTwo MeshTwo 
