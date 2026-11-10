F90=gfortran
# FFLAGS=-O0 -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
# FFLAGS=-O2 -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
FFLAGS=-O3 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall
EXE=exner_1d

all : $(EXE)

$(EXE) : mod_parametres.o mod_qte_physiques.o mod_fonctions.o mod_legendre.o mod_interpolation.o mod_limiteur.o mod_resolution.o mod_resolution_ODE.o exner.o
	$(F90) -o $@ $^

%.o : %.f90
	$(F90) -c $(FFLAGS) $<

clean :
	rm -f *.mod *.o $(EXE) ./SORTIE/*.dat