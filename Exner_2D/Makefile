F90=gfortran
# FFLAGS=-O0 -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all
# FFLAGS=-O2 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -fopenmp
# FFLAGS=-O3 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -fopenmp
FFLAGS=-O3 -march=native -ffast-math -c -cpp -ffree-line-length-0 -fstack-arrays -fopenmp
EXE=exner_2d
LIBS=-llapack -lblas -fopenmp

all : $(EXE)

$(EXE) : mod_parametres.o mod_qte_physiques.o mod_fonctions.o mod_types_derives.o mod_tri_maillage.o mod_maillage.o mod_legendre.o mod_base.o mod_cl.o mod_interpolation.o mod_limiteur.o mod_resolution.o exner_2d.o
	$(F90) -o $@ $^ $(LIBS)

%.o : %.f90
	$(F90) $(FFLAGS) $<

clean :
	rm -f *.mod *.o $(EXE) ./SORTIE/*.dat