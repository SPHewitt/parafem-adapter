# -----------------
# Compiler
# -----------------
FC=mpif90
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

PARAFEMROOT=/home/mclsssh5/paraFem/parafem
PARAFEMLIB=$(PARAFEMROOT)/lib
PARAFEMMOD=$(PARAFEMROOT)/include/mpi

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
all:
	@echo $(SRC)
	@echo $(OBJ)
	$(MAKE) $(EXE) 

$(EXE): $(OBJ)
	@echo Compiling main
	$(FC) $(FCFLAGS) -I$(PARAFEMMOD) $(EXE).f90 -o $(BINDIR)/$@ $^ -L$(PRECICELIB) -lprecice -L$(PARAFEMLIB) -lParaFEM_mpi.5.0.3 -larpack_linuxdesktop

#$(OBJ): $(SRC)
#	@echo Compiling Fortran Bindings
#	$(FC) $(FFLAGS) -c $< -o $@ -J$(MODDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@echo Compiling "$@"
	$(FC) $(CFLAGS) -I$(PARAFEMMOD) -o $@ -c $< -J$(MODDIR) -L$(PARAFEMLIB) -lParaFEM_mpi.5.0.3 -larpack_linuxdesktop

clean:
	@echo cleaning
	rm -rf obj/*.o
	rm -rf main
	rm -rf *.log
run:
	./bin/parafem precice-config.xml SolverTwo MeshTwo 
