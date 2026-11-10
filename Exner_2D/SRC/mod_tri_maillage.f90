module mod_tri_maillage

    use mod_parametres

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_tri_maillage :
!   - get_S_i_1_back
!   - fill_Flag
!   - fill_Sommets_Arete_a
!   - fill_Aretes_Omega_i_Trig

! Module inspire du module mod_maillage.f90 cree par Luc Mieussens :
! https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK2/VF/TP/mod_maillage.f90
! ==================================================================================================== !

    contains

        subroutine get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)
! ==================================================================================================== !
! Subroutine renvoyant les deux sommets successifs s_k et s_l relies par l'arete a
! N.B. Elle retourne Sommets_Omega_i(1, i) si s = N_sommets_Omega_i(i)

! Entrees de la subroutine
!   - i : Entier caracterisant la maille Omega_i
!   - s : s-eme sommet de la maille Omega_i
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i

! Sorties de la subroutine
!   - s_k : Sommets_Omega_i(s, i)
!   - s_l : Sommet suivant le sommet s_k
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i, s
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i

! Sorties de la subroutine
            integer, intent(out)                                        :: s_k, s_l

            s_k = Sommets_Omega_i(MOD(s - 1, N_sommets_Omega_i(i)) + 1, i)
            s_l = Sommets_Omega_i(MOD(s, N_sommets_Omega_i(i)) + 1, i)

        end subroutine get_S_i_1_back

        subroutine fill_Flag(N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets, N_aretes, Flag)
! ==================================================================================================== !
! Subroutine renvoyant le tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant les
! sommets s_k et s_l en fin de boucle

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i

! Sorties de la subroutine :
!   - N_aretes : Nombre d'aretes du maillage
!   - Flag(1:N_noeuds, 1:N_noeuds) : Tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant 
! les sommets s_k et s_l
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i

            real(kind = pr), dimension(:, :), intent(in)                :: Coord_Sommets

! Sorties de la subroutine
            integer                                                     :: N_aretes
            integer, dimension(:, :), allocatable, intent(out)          :: Flag

! Variables locales
            integer                                                     :: i, s, s_k, s_l, N_elements, N_noeuds

! Nombre de mailles
            N_elements = SIZE(Sommets_Omega_i, 2)

! Nombre de sommets
            N_noeuds = SIZE(Coord_Sommets, 2)

            allocate(Flag(1:N_noeuds, 1:N_noeuds))

! Remplissage de Flag
            N_aretes = 0 ; Flag = 0

            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)

! Marque (ou pose un drapeau !) l'arete entre les sommets s_k et s_l si elle n'est pas deja marquee
                    if (Flag(s_k, s_l) == 0 .OR. Flag(s_l, s_k) == 0) then
                        N_aretes = N_aretes + 1
                        Flag(s_k, s_l) = N_aretes ; Flag(s_l, s_k) = N_aretes
                    end if
                end do
            end do

        end subroutine fill_Flag

        subroutine fill_Sommets_Arete_a(N_sommets_Omega_i, Sommets_Omega_i, N_aretes, Flag, Sommets_Arete_a)
! ==================================================================================================== !
! Subroutine renvoyant le tableau Sommets_Arete_a content les deux sommets relies par l'arete a

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - N_aretes : Nombre d'aretes du maillage
!   - Flag(1:N_noeuds, 1:N_noeuds) : Tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant 
! les sommets s_k et s_l

! Sortie de la subroutine :
!   - Sommets_Arete_a(1:2, 1:N_aretes) : Tableau contenant les deux sommets relies par l'arete a
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: N_aretes
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i, Flag

! Sortie de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: Sommets_Arete_a

! Variables locales
            integer                                                     :: i, s, s_k, s_l, N_elements

            N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(Sommets_Arete_a(1:2, 1:N_aretes))

            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)

                    Sommets_Arete_a(1, Flag(s_k, s_l)) = s_k ; Sommets_Arete_a(2, Flag(s_k, s_l)) = s_l
                end do
            end do

        end subroutine fill_Sommets_Arete_a

        subroutine fill_Aretes_Omega_i(N_sommets_Omega_i, Sommets_Omega_i, Flag, Aretes_Omega_i)
! ==================================================================================================== !
! Subroutine renvoyant les tableaux Aretes_Omega_i et Trig

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Flag(1:N_noeuds, 1:N_noeuds) : Tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant 
! les sommets s_k et s_l

! Sortie de la subroutine :
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i, Flag

! Sortie de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: Aretes_Omega_i

! Variables locales
            integer                                                     :: i, s, s_k, s_l, N_max_sommets, N_elements

            N_max_sommets = SIZE(Sommets_Omega_i, 1) ; N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(Aretes_Omega_i(1:N_max_sommets, 1:N_elements))

! Remplissage de Aretes_Omega_i
            Aretes_Omega_i = 0

            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)

! Recupere le numero associe a l'arete liant les sommets s_k et s_l
                    Aretes_Omega_i(s, i) = Flag(s_k, s_l)
                end do
            end do

        end subroutine fill_Aretes_Omega_i

        subroutine fill_Trig_Trig_F(N_sommets_Omega_i, Sommets_Omega_i, Flag, N_aretes, Trig, Trig_F)
! ==================================================================================================== !
! Subroutine remplissant les tableaux Trig et Trig_F

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Flag(1:N_noeuds, 1:N_noeuds) : Tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant 
! les sommets s_k et s_l
!   - N_aretes : Nombre d'aretes du maillage

! Sorties de la subroutine :
!   - Trig(1:2, 1:N_aretes) : Tableau contenant le numero des deux mailles adjacentes par l'arete a
!   - Trig_F(1:2, 1:N_aretes) : Tableau contenant le numero des faces 1 <= F, F_r <= N_max_sommets
! reliant les mailles Omega_i et Omega_r
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: N_aretes
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i, Flag
            
! Sorties de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: Trig, Trig_F

! Variables locales
            integer                                                     :: N_elements, i, s, s_k, s_l, a

            N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(Trig(1:2, 1:N_aretes), Trig_F(1:2, 1:N_aretes))
            Trig = 0 ; Trig_F = 0

! Remplissage de Trig et Trig_F
            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)
                    a = Flag(s_k, s_l)

                    if (Trig(1, a) == 0) then
                        ! Première rencontre : on enregistre la maille et la face
                        Trig(1, a) = i
                        Trig_F(1, a) = s
                    else if (Trig(1, a) /= i) then
                        ! Deuxième rencontre : on enregistre le voisin
                        Trig(2, a) = i
                        Trig_F(2, a) = s
                    end if
                end do
            end do

        end subroutine fill_Trig_Trig_F

        subroutine fill_N_Omega_r_Faces_Omega_r(N_sommets_Omega_i, Sommets_Omega_i, Flag, Trig, Trig_F, &
    &                                           N_Omega_r, Faces_Omega_r)
! ==================================================================================================== !
! Subroutine remplissant les tableaux N_Omega_r et Faces_Omega_r

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Flag(1:N_noeuds, 1:N_noeuds) : Tableau Flag tel que Flag(s_k, s_l) contiennent l'arete a liant 
! les sommets s_k et s_l
!   - N_aretes : Nombre d'aretes du maillage

! Sorties de la subroutine :
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F
!   - Faces_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r, i.e. F_r
! ==================================================================================================== !

! Entrees de la subroutine
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i, Flag, Trig, Trig_F

! Sorties de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: N_Omega_r, Faces_Omega_r

! Variables locales
            integer                                                     :: i, s, s_k, s_l, a, N_max_sommets, N_elements

            N_max_sommets = SIZE(Sommets_Omega_i, 1) ; N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(N_Omega_r(1:N_max_sommets, 1:N_elements), Faces_Omega_r(1:N_max_sommets, 1:N_elements))
            N_Omega_r = 0 ; Faces_Omega_r = 0

! Remplissage de N_Omega_r et Faces_Omega_r
            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)
                    a = Flag(s_k, s_l)

                    if (Trig(1, a) == i) then
                        N_Omega_r(s, i) = Trig(2, a)
                        Faces_Omega_r(s, i) = Trig_F(2, a)
                    else if (Trig(2, a) == i) then
                        N_Omega_r(s, i) = Trig(1, a)
                        Faces_Omega_r(s, i) = Trig_F(1, a)
                    end if
                end do
            end do

        end subroutine fill_N_Omega_r_Faces_Omega_r

        subroutine fill_N_Omega_r_P_Faces_Omega_r_P(N_sommets_Omega_i, Sommets_Omega_i, Coord_Sommets, Sommets_Arete_a,     &
    &   N_aretes, Trig, Trig_F, N_Omega_r, N_Omega_r_P, Faces_Omega_r_P)
! ==================================================================================================== !
! Subroutine construisant les tableaux N_Omega_r_P et Faces_Omega_r_P

! Entrees de la subroutine :
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i
!   - Sommets_Arete_a(1:2, 1:N_aretes) : Tableau contenant les deux sommets relies par l'arete a
!   - N_aretes : Nombre d'aretes du maillage
!   - Trig(1:2, 1:N_aretes) : Tableau contenant le numero des deux mailles adjacentes par l'arete a
!   - Trig_F(1:2, 1:N_aretes) : Tableau contenant le numero des faces 1 <= F, F_r <= N_max_sommets
! reliant les mailles Omega_i et Omega_r
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F

! Sorties de la subroutine :
!   - N_Omega_r_P(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F dans le cas de conditions periodiques
!   - Faces_Omega_r_P(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r, i.e. F_r dans le cas de conditions periodiques
! De plus, F_r est negatif si les faces F et F_r ont le meme sens de parcours et positif sinon
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: N_aretes
            integer, dimension(:), intent(in)                           :: N_sommets_Omega_i
            integer, dimension(:, :), intent(in)                        :: Sommets_Omega_i, Sommets_Arete_a, Trig, Trig_F
            integer, dimension(:, :), intent(in)                        :: N_Omega_r

            real(kind = pr), dimension(:, :), intent(in)                :: Coord_Sommets

! Sorties de la subroutine
            integer, dimension(:, :), allocatable, intent(out)          :: N_Omega_r_P, Faces_Omega_r_P

! Variables locales
            integer                                                     :: N_max_sommets, N_elements, i, s, s_k, s_l, s_n, s_m, a

            real(kind = pr), dimension(2)                               :: Coord_Centre_F, Coord_Centre_P, Coord_Candidat

            N_max_sommets = SIZE(Sommets_Omega_i, 1) ; N_elements = SIZE(Sommets_Omega_i, 2)

            allocate(N_Omega_r_P(1:N_max_sommets, 1:N_elements))
            allocate(Faces_Omega_r_P(1:N_max_sommets, 1:N_elements))

            N_Omega_r_P = 0 ; Faces_Omega_r_P = 0

            do i = 1, N_elements
                do s = 1, N_sommets_Omega_i(i)
                    if (N_Omega_r(s, i) == 0) then
                        call get_S_i_1_back(i, s, N_sommets_Omega_i, Sommets_Omega_i, s_k, s_l)

                        Coord_Centre_F = 1._pr/2*(Coord_Sommets(:, s_k) + Coord_Sommets(:, s_l))

                        Coord_Centre_P = Coord_Centre_F
                        if (ABS(Coord_Centre_F(1)) < eps_P) Coord_Centre_P(1) = Coord_Centre_F(1) + Long
                        if (ABS(Coord_Centre_F(1) - Long) < eps_P) Coord_Centre_P(1) = Coord_Centre_F(1) - Long

                        if (ABS(Coord_Centre_F(2)) < eps_P) Coord_Centre_P(2) = Coord_Centre_F(2) + Larg
                        if (ABS(Coord_Centre_F(2) - Larg) < eps_P) Coord_Centre_P(2) = Coord_Centre_F(2) - Larg

                        do a = 1, N_aretes
                            if (Trig(2, a) == 0) then
                                s_n = Sommets_Arete_a(1, a) ; s_m = Sommets_Arete_a(2, a)

                                Coord_Candidat = 1._pr/2*(Coord_Sommets(:, s_n) + Coord_Sommets(:, s_m))

                                if (MAXVAL(ABS(Coord_Centre_P - Coord_Candidat)) < eps_P) then
                                    N_Omega_r_P(s, i) = Trig(1, a)

                                    if (DOT_PRODUCT(Coord_Sommets(:, s_l) - Coord_Sommets(:, s_k),                          &
                                &                   Coord_Sommets(:, s_m) - Coord_Sommets(:, s_n)) > 0._pr) then
! Les faces F et F_r sont orientees dans le meme sens
                                        Faces_Omega_r_P(s, i) = -Trig_F(1, a)
                                    else
                                        Faces_Omega_r_P(s, i) = Trig_F(1, a)
                                    end if
                                end if
                            end if
                        end do
                    end if
                end do
            end do

        end subroutine fill_N_Omega_r_P_Faces_Omega_r_P

end module mod_tri_maillage