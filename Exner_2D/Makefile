F90=gfortran
# FFLAGS=-O0 -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -JSRC
# FFLAGS=-O2 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -fopenmp -JSRC
# FFLAGS=-O3 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -fopenmp -JSRC
FFLAGS=-O3 -march=native -ffast-math -c -cpp -ffree-line-length-0 -fstack-arrays -fopenmp -JSRC
EXE=exner_2d
LIBS=-llapack -lblas -fopenmp

# Dossier contenant les sources
SRC_DIR = SRC

# Liste des fichiers objets localisés dans le dossier SRC/
OBJS = $(SRC_DIR)/mod_parametres.o \
       $(SRC_DIR)/mod_qte_physiques.o \
       $(SRC_DIR)/mod_fonctions.o \
       $(SRC_DIR)/mod_types_derives.o \
       $(SRC_DIR)/mod_tri_maillage.o \
       $(SRC_DIR)/mod_maillage.o \
       $(SRC_DIR)/mod_legendre.o \
       $(SRC_DIR)/mod_base.o \
       $(SRC_DIR)/mod_cl.o \
       $(SRC_DIR)/mod_interpolation.o \
       $(SRC_DIR)/mod_limiteur.o \
       $(SRC_DIR)/mod_resolution.o \
       $(SRC_DIR)/exner_2d.o

all : $(EXE)

# Édition de liens pour générer l'exécutable à la racine
$(EXE) : $(OBJS)
	@echo "Création de l'exécutable : $@"
	$(F90) -o $@ $^ $(LIBS)

# Règle générique de compilation pour tout fichier .f90 situé dans SRC/
$(SRC_DIR)/%.o : $(SRC_DIR)/%.f90
	@echo "Compilation de $<..."
	$(F90) $(FFLAGS) $< -o $@

clean :
	rm -f $(SRC_DIR)/*.mod $(SRC_DIR)/*.o $(EXE) ./SORTIE/*.dat
	@echo "Nettoyage effectué."

.PHONY : all clean