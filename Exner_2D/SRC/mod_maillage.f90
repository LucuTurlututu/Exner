module mod_maillage

    use mod_parametres
    use mod_types_derives
    use mod_tri_maillage

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_maillage :
!   - maillage
!   - lecture_maillage
!   - extraire_nom_maillage
!   - connectivite
!   - aire
!   - cl_arete

! Module inspire du module mod_maillage.f90 cree par Luc Mieussens :
! https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK2/VF/TP/mod_maillage.f90
! ==================================================================================================== !

    contains

        subroutine make_Mesh(fichier, Mesh)
! ==================================================================================================== !
! Subroutine recuperant toutes les donnes de construction du maillage et renvoyant les donnes
! essentielles

! Entree de la subroutine :
!   - fichier : Chemin vers le fichier contenant le maillage

! Entree/Sortie de la subroutine
!   - Mesh : Type derive contenant l'ensemble des sorties de la subroutine maillage

! Sorties de la subroutine :
!   - N_elements : Nombre de mailles
!   - N_aretes : Nombre d'aretes
!   - N_faces_Omega_i(1:N_elements) : Tableau contenant le nombre de faces de la maille Omega_i
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i
!   - Aire_Omega_i(1:N_elements) : Tableau contenant l'aire de toutes les mailles Omega_i
!   - Long_arete(1:N_aretes) : Tableau contenant la longueur de l'arete a
!   - Diam_Omega_i(1:N_elements) : Tableau contenant le diametre de la maille Omega_i
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille Omega_r,
! i.e. la maile adjacente a la maille Omega_i par la face F
!   - Faces_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r
!   - N_Omega_r_P(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F dans le cas de conditions periodiques
!   - Faces_Omega_r_P(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r, i.e. F_r dans le cas de conditions periodiques
! De plus, F_r est negatif si les faces F et F_r ont le meme sens de parcours et positif sinon
!   - n_F(1:2, 1:N_max_sommets, 1:N_elements) : Tableau contenant les composantes (n_x, n_y) de la 
! normale unitaire sortante de la face F de la maille Omega_i
!   - CL_Arete_bord(1:N_aretes) : Tableau contenant des lettres referant a la condition aux limites
! appliquee sur le bord :
!       "D" : Condition aux limites de Dirichlet
!       "N" : Conditions aux limites de Neumann
!       "P" : Coniditions aux limites periodiques
!       "O" : L'arete a n'est pas une arete de bord

! N.B. : Si Omega_i a (N_max_sommets - 1) sommets, Sommets_Omega_i(N_max_sommets, i) = 0 et ainsi
! de suite
! De meme, si Omega_i a (N_max_sommets - 1) aretes, Aretes_Omega_i(N_max_sommets, i) = 0 et ainsi
! de suite
! ==================================================================================================== !

! Entree de la subroutine
            character(len = *), intent(in)                              :: fichier

! Entree/Sortie de la subroutine
            type(Maillage), intent(inout)                               :: Mesh

! Variables locales
            integer                                                     :: N_noeuds, a, s, i, a_s_i, F
            integer, dimension(:, :), allocatable                       :: Sommets_Arete_a, Trig, Trig_F

            real(kind = pr)                                             :: Peri_i, Scal
            real(kind = pr), dimension(2)                               :: Coord_s_1, Coord_s_2, G_Omega_i, Milieu_arete_a

! Lecture du maillage
            call lecture_maillage(fichier, Mesh%N_faces_Omega_i, Mesh%N_sommets_Omega_i, Mesh%Sommets_Omega_i,              &
        &                                                                                                 Mesh%Coord_Sommets)

            Mesh%N_elements = SIZE(Mesh%Sommets_Omega_i, 2) ; N_noeuds = SIZE(Mesh%Coord_Sommets, 2)
            Mesh%N_max_sommets = SIZE(Mesh%Sommets_Omega_i, 1)

            print *, "Nombre de sommets : ", N_noeuds
            print *, "Nombre de mailles : ", Mesh%N_elements

! Connectivite
            call connectivite(Mesh%N_sommets_Omega_i, Mesh%Sommets_Omega_i, Mesh%Coord_Sommets, Sommets_Arete_a,            &
        &                     Mesh%Aretes_Omega_i, Trig, Trig_F, Mesh%N_Omega_r, Mesh%Faces_Omega_r)

            Mesh%N_aretes = SIZE(Sommets_Arete_a, 2)

! Aire des mailles
            call aire(Mesh%N_sommets_Omega_i, Mesh%Sommets_Omega_i, Mesh%Coord_Sommets, Mesh%Aire_Omega_i)

! Longueur des aretes
            allocate(Mesh%Long_arete(1:Mesh%N_aretes))
            
            do a = 1, Mesh%N_aretes
                Coord_s_1 = Mesh%Coord_Sommets(:, Sommets_Arete_a(1, a))
                Coord_s_2 = Mesh%Coord_Sommets(:, Sommets_Arete_a(2, a))

                Mesh%Long_arete(a) = SQRT(SUM((Coord_s_2 - Coord_s_1)**2))
            end do

! Calcul des normales, du diametre et du perimetre de Omega_i
            allocate(Mesh%Diam_Omega_i(1:Mesh%N_elements))
            allocate(Mesh%n_F(1:2, 1:Mesh%N_max_sommets, 1:Mesh%N_elements))

            Mesh%n_F = 0._pr
            do i = 1, Mesh%N_elements
                Peri_i = 0._pr

! Barycentre de la maille Omega_i
                G_Omega_i = 0._pr
                do s = 1, Mesh%N_sommets_Omega_i(i)
                    G_Omega_i = G_Omega_i + Mesh%Coord_Sommets(:, Mesh%Sommets_Omega_i(s, i))
                end do
                G_Omega_i = G_Omega_i/Mesh%N_sommets_Omega_i(i)

                do F = 1, Mesh%N_sommets_Omega_i(i)
                    a_s_i = Mesh%Aretes_Omega_i(F, i)
                    
                    Peri_i = Peri_i + Mesh%Long_arete(a_s_i)

                    Coord_s_1 = Mesh%Coord_Sommets(:, Sommets_Arete_a(1, a_s_i))
                    Coord_s_2 = Mesh%Coord_Sommets(:, Sommets_Arete_a(2, a_s_i))

! Calcul de la normale de la face F
                    Mesh%n_F(1, F, i) = (Coord_s_2(2) - Coord_s_1(2))/Mesh%Long_arete(a_s_i)
                    Mesh%n_F(2, F, i) = -(Coord_s_2(1) - Coord_s_1(1))/Mesh%Long_arete(a_s_i)

                    Milieu_arete_a = 1._pr/2*(Coord_s_1 + Coord_s_2)

                    Scal = DOT_PRODUCT((Milieu_arete_a - G_Omega_i), Mesh%n_F(:, F, i))

! Si le produit scalaire est negatif, alors la normale est entrante
                    if (Scal < 0._pr) then
                        Mesh%n_F(:, F, i) = -Mesh%n_F(:, F, i)
                    end if
                end do

                ! Diametre de la maille i
                Mesh%Diam_Omega_i(i) = Mesh%Aire_Omega_i(i)/Peri_i
            end do

! Conditions aux limites sur les aretes de bord
            call cl_arete(Mesh%N_sommets_Omega_i, Mesh%Coord_Sommets, Sommets_Arete_a, Mesh%Aretes_Omega_i,                 &
        &                 Mesh%Faces_Omega_r, Mesh%CL_Arete_bord)

! Remplissage de N_Omega_r_P et Faces_Omega_r_P
            if (ANY(Mesh%CL_Arete_Bord == "P")) then
                call fill_N_Omega_r_P_Faces_Omega_r_P(Mesh%N_sommets_Omega_i, Mesh%Sommets_Omega_i, Mesh%Coord_Sommets,     &
            &           Sommets_Arete_a, Mesh%N_aretes, Trig, Trig_F, Mesh%N_Omega_r, Mesh%N_Omega_r_P, Mesh%Faces_Omega_r_P)
            end if

            deallocate(Sommets_Arete_a, Trig, Trig_F)

        end subroutine make_Mesh

        subroutine lecture_maillage(fichier, N_faces_Omega_i, N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets)
! ==================================================================================================== !
! Subroutine lecture_maillage recuperant les informations directement contenues dans le fichier
! de maillage

! Entree de la subroutine :
!   - fichier : Chemin vers le fichier contenant le maillage

! Sorties de la subroutine :
!   - N_faces_Omega_i(1:N_elements) : Tableau contenant le nombre de faces de la maille Omega_i
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i
! ==================================================================================================== !

! Entree de la subroutine
            character(len = *), intent(in)                              :: fichier

! Sorties de la subroutine
            integer, dimension(:), allocatable, intent(out)             :: N_faces_Omega_i, N_sommets_Omega_i
            integer, dimension(:, :), allocatable, intent(out)          :: Sommets_Omega_i

            real(kind = pr), dimension(:, :), allocatable, intent(out)  :: Coord_Sommets

! Variables locales
            character(len = 4)                                          :: Extension
            character(len = 255)                                        :: Nom_maillage

            integer                                                     :: i, N_noeuds, N_max_sommets, N_elements
            integer                                                     :: len_Nom_maillage, N_sommets, id_noeud
            integer                                                     :: trash_1, trash_2, trash_3, trash_4, type_maille
            integer                                                     :: max_type_maille

            open(unit = 10, file = TRIM(ADJUSTL(fichier)))

! Recupere le nom du maillage sans le chemin
            call extraire_nom_maillage(fichier, Nom_maillage)

! Recupere l'extension de Nom_maillage
            len_Nom_maillage = LEN(TRIM(ADJUSTL(Nom_maillage)))
            extension = Nom_maillage(len_Nom_maillage-3:len_Nom_maillage)

            if (extension .EQ. "typ2") then
! Lit la premiere ligne et le nombre de sommets
                read(10, *)
                read(10, *) N_noeuds

                allocate(Coord_Sommets(1:2, 1:N_noeuds))

                do i = 1, N_noeuds
                    read(10, *) Coord_Sommets(1:2, i)
                end do

! Lit "cells" puis lit le nombre de mailles
                read(10, *)
                read(10, *) N_elements

! 1e lecture pour avoir N_max_sommets
                N_max_sommets = 0
                do i = 1, N_elements
                    read(10, *) N_sommets

                    if (N_sommets > N_max_sommets) N_max_sommets = N_sommets
                end do

                do i = 1, N_elements
                    backspace(10)
                end do

                allocate(N_faces_Omega_i(1:N_elements), N_sommets_Omega_i(1:N_elements))
                allocate(Sommets_Omega_i(N_max_sommets, 1:N_elements))

                Sommets_Omega_i = 0

                do i = 1, N_elements
! Lit le nombre de sommets de la maille puis leur numero
                    read(10, *) N_sommets_Omega_i(i), Sommets_Omega_i(1:N_sommets_Omega_i(i), i)

                    N_faces_Omega_i(i) = N_sommets_Omega_i(i)
                end do

            else if (extension .EQ. ".msh") then
! Lit les huit premieres lignes puis le nombre de sommets
                do i = 1, 8
                    read(10, *)
                end do
                read(10, *) N_noeuds

                allocate(Coord_Sommets(1:2, 1:N_noeuds))

                do i = 1, N_noeuds
                    read(10, *) id_noeud, Coord_Sommets(1, id_noeud), Coord_Sommets(2, id_noeud)
                end do

                Long = MAXVAL(Coord_Sommets(1, :)) - MINVAL(Coord_Sommets(1, :))
                Larg = MAXVAL(Coord_Sommets(2, :)) - MINVAL(Coord_Sommets(2, :))

! Lit "$EndNodes", "$Elements" puis le nombre d'elements
                read(10, *)
                read(10, *)
                read(10, *) N_elements

! 1e lecture pour avoir le nombre maximal de sommets
                N_max_sommets = 0 ; max_type_maille = 0
                do i = 1, N_elements
                    read (10, *) trash_1, type_maille

                    if (type_maille > max_type_maille) max_type_maille = type_maille
                end do

                if (max_type_maille == 2) then
                    N_max_sommets = 3
                else if (max_type_maille == 3) then
                    N_max_sommets = 4
                else
                    print *, "lecture_mailage : Nombre maximal de sommets non pris en compte"
                    stop
                end if

                do i = 1, N_elements
                    backspace(10)
                end do

                allocate(N_faces_Omega_i(1:N_elements), N_sommets_Omega_i(1:N_elements))
                allocate(Sommets_Omega_i(N_max_sommets, 1:N_elements))

                Sommets_Omega_i = 0

                do i = 1, N_elements
! Lit le nombre de sommets de la maille puis leur numero
                    read(10, *) trash_1, type_maille, trash_2, trash_3, trash_4, Sommets_Omega_i(1:type_maille+1, i)

                    N_sommets_Omega_i(i) = type_maille + 1
                    N_faces_Omega_i(i) = N_sommets_Omega_i(i)
                end do
            end if

            close(10)

        end subroutine lecture_maillage

        subroutine extraire_nom_maillage(fichier, Nom_maillage)
! ==================================================================================================== !
! Subroutine recuperant le nom du maillage sans le chemin associe vers le dossier contenant le
! fichier

! Entree de la subroutine :
!   - fichier : Chemin vers le fichier contenant le maillage

! Sortie de la subroutine
!   - Nom_maillage : Nom du fichier contenant le maillage
! ==================================================================================================== !

! Entree de la subroutine
            character(len = *), intent(in)                              :: fichier

! Sortie de la subroutine
            character(len = *), intent(out)                             :: Nom_maillage

! Variables locales
            integer                                                     :: i, Slash_chemin

            Nom_maillage = " "

! On cherche la position du dernier /
            Slash_chemin = 0
            do i = LEN(fichier), 1, -1
                if (fichier(i:i) == "/") then
                    Slash_chemin = i
                    exit
                end if
            end do

! Si un / est trouve, alors on extrait l portion apres le dernier /
            if (Slash_chemin > 0) then
                Nom_maillage = fichier(Slash_chemin+1:)
! Si aucun / n'est trouve, alors on a deja Nom_maillage
            else
                Nom_maillage = fichier
            end if

        end subroutine extraire_nom_maillage

        subroutine connectivite(N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets, Sommets_Arete_a,     &
        &                       Aretes_Omega_i, Trig, Trig_F, N_Omega_r, Faces_Omega_r)
! ==================================================================================================== !
! Subroutine generant toutes les informations de connectivite

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i

! Sorties de la subroutine :
!   - Sommets_Arete_a(1:2, 1:N_aretes) : Tableau contenant les deux sommets relies par l'arete a
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i
!   - Trig(1:2, 1:N_aretes) : Tableau contenant le numero des deux mailles adjacentes par l'arete a
!   - Trig_F(1:2, 1:N_aretes) : Tableau contenant le numero des faces 1 <= F, F_r <= N_max_sommets
! reliant les mailles Omega_i et Omega_r
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F
!   - Faces_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r
! De plus, F_r est negatif si les faces F et F_r ont le meme sens de parcours et positif sinon
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i

            real(kind = pr), dimension(:, :), intent(in)                :: Coord_Sommets

! Sorties de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: Sommets_Arete_a, Aretes_Omega_i
            integer, dimension(:, :), allocatable, intent(out)          :: Trig, Trig_F, N_Omega_r, Faces_Omega_r

! Variables locales
            integer                                                     :: N_aretes
            integer, dimension(:, :), allocatable                       :: Flag

            call fill_Flag(N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets, N_aretes, Flag)

            call fill_Sommets_Arete_a(N_sommets_Omega_i, Sommets_Omega_i, N_aretes, Flag, Sommets_Arete_a)

            call fill_Aretes_Omega_i(N_sommets_Omega_i, Sommets_Omega_i, Flag, Aretes_Omega_i)

            call fill_Trig_Trig_F(N_sommets_Omega_i, Sommets_Omega_i, Flag, N_aretes, Trig, Trig_F)

            call fill_N_Omega_r_Faces_Omega_r(N_sommets_Omega_i, Sommets_Omega_i, Flag, Trig, Trig_F, N_Omega_r, Faces_Omega_r)

            deallocate(Flag)

        end subroutine connectivite

        subroutine aire(N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets, Aire_Omega_i)
! ==================================================================================================== !
! Subroutine calculant l'aire de la maille Omega_i

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i

! Sortie de la subroutine :
!   - Aire_Omega_i(1:N_elements) : Tableau contenant l'aire de toutes les mailles Omega_i
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i

            real(kind = pr), dimension(:, :), intent(in)                :: Coord_Sommets

! Sortie de la subroutine
            real(kind = pr), dimension(:), allocatable, intent(out)     :: Aire_Omega_i

! Variables locales
            integer                                                     :: i, s, s_k, s_l, N_elements

            real(kind = pr)                                             :: Somme_Aire
            real(kind = pr), dimension(2)                               :: Coord_s_k, Coord_s_l

            N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(Aire_Omega_i(1:N_elements))

            do i = 1, N_elements
                Somme_Aire = 0._pr
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)

                    Coord_s_k = Coord_Sommets(:, s_k) ; Coord_s_l = Coord_Sommets(:, s_l)

                    Somme_Aire = Somme_Aire + 1._pr/2*(Coord_s_k(1)*Coord_s_l(2) - Coord_s_l(1)*Coord_s_k(2))
                end do
                Aire_Omega_i(i) = Somme_Aire
            end do

        end subroutine aire

        subroutine cl_arete(N_sommets_Omega_i, Coord_Sommets, Sommets_Arete_a, Aretes_Omega_i,          &
        &                   Faces_Omega_r, CL_Arete_bord)
! ==================================================================================================== !
! Subroutine affectant une condition aux limites ou non a chaque arete

! Entrees de la subroutine :
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i
!   - Sommets_Arete_a(1:2, 1:N_aretes) : Tableau contenant les deux sommets relies par l'arete a
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i
!   - Faces_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r, i.e. F_r

! Sortie de la subroutine :
!   - CL_Arete_bord(1:N_aretes) : Tableau contenant des lettres referant a la condition aux limites
! appliquee sur le bord :
!       "D" : Condition aux limites de Dirichlet
!       "N" : Conditions aux limites de Neumann
!       "P" : Coniditions aux limites periodiques
!       "O" : L'arete a n'est pas une arete de bord
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Arete_a, Aretes_Omega_i, Faces_Omega_r

            real(kind = pr), dimension(:, :), intent(in)                :: Coord_Sommets

! Sortie de la subroutine
            character(len = 1), dimension(:), allocatable, intent(out)  :: CL_Arete_bord

! Variables locales
            integer                                                     :: a, F, N_aretes, N_elements, i, s_k, s_l

            real(kind = pr)                                             :: Tol
            real(kind = pr), dimension(2)                               :: Coord_s_k, Coord_s_l

            N_aretes = SIZE(Sommets_Arete_a, 2) ; N_elements = SIZE(Aretes_Omega_i, 2)

            allocate(CL_Arete_bord(1:N_aretes))

            Tol = 1.e-13_pr
            CL_Arete_bord = "O"

            do i = 1, N_elements
                do F = 1, N_sommets_Omega_i(i)
                    if (Faces_Omega_r(F, i) == 0) then
                        a = Aretes_Omega_i(F, i)

                        s_k = Sommets_Arete_a(1, a) ; s_l = Sommets_Arete_a(2, a)
                        Coord_s_k = Coord_Sommets(:, s_k) ; Coord_s_l = Coord_Sommets(:, s_l)
                    
! On teste les differentes cotes de la maille
                        if (ABS(Coord_s_k(1)) <= Tol .AND. ABS(Coord_s_l(1)) <= Tol) then ! Bord Gauche
                            select case (cas)

                            case (1, 2, 3, 4, 5, 6, 10)
                                CL_Arete_bord(a) = "P"

                            case (7, 8, 9)
                                CL_Arete_bord(a) = "D"
                            
                            end select
                        else if (ABS(Coord_s_k(1)) >= Long-Tol .AND. ABS(Coord_s_l(1)) >= Long-Tol) then ! Bord Droit
                            select case (cas)

                            case (1, 2, 3, 4, 5, 6, 10)
                                CL_Arete_bord(a) = "P"

                            case (7, 8, 9)
                                CL_Arete_bord(a) = "D"
                            
                            end select
                        else if (ABS(Coord_s_k(2)) <= Tol .AND. ABS(Coord_s_l(2)) <= Tol) then ! Bord Bas
                            select case (cas)

                            case (1, 2, 3, 4, 5, 6, 10)
                                CL_Arete_bord(a) = "P"

                            case (7, 8, 9)
                                CL_Arete_bord(a) = "N"
                            
                            end select
                        else if (ABS(Coord_s_k(2)) >= Larg-Tol .AND. ABS(Coord_s_l(2)) >= Larg-Tol) then ! Bord haut
                            select case (cas)

                            case (1, 2, 3, 4, 5, 6, 10)
                                CL_Arete_bord(a) = "P"

                            case (7, 8, 9)
                                CL_Arete_bord(a) = "N"
                            
                            end select
                        end if
                    end if
                end do
            end do

        end subroutine cl_arete

end module mod_maillage