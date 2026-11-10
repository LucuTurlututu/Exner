F90=gfortran
# FFLAGS=-O0 -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -JSRC
# FFLAGS=-O2 -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -JSRC
FFLAGS=-O3 -c -cpp -g -fcheck=all -ffree-line-length-0 -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -JSRC
EXE=exner_1d

# Dossier contenant les sources
SRC_DIR = SRC

# Liste des fichiers objets localisés dans le dossier SRC/
OBJS = $(SRC_DIR)/mod_parametres.o \
       $(SRC_DIR)/mod_qte_physiques.o \
       $(SRC_DIR)/mod_fonctions.o \
       $(SRC_DIR)/mod_legendre.o \
       $(SRC_DIR)/mod_interpolation.o \
       $(SRC_DIR)/mod_limiteur.o \
       $(SRC_DIR)/mod_resolution.o \
       $(SRC_DIR)/mod_resolution_ODE.o \
       $(SRC_DIR)/exner.o

all : $(EXE)

# Édition de liens pour générer l'exécutable à la racine
$(EXE) : $(OBJS)
	@echo "Création de l'exécutable : $@"
	$(F90) -o $@ $^

# Règle générique de compilation pour tout fichier .f90 situé dans SRC/
$(SRC_DIR)/%.o : $(SRC_DIR)/%.f90
	@echo "Compilation de $<..."
	$(F90) $(FFLAGS) $< -o $@

clean :
	rm -f $(SRC_DIR)/*.mod $(SRC_DIR)/*.o $(EXE) ./SORTIE/*.dat
	@echo "Nettoyage effectué."

.PHONY : all clean