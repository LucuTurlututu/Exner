module mod_base

    use mod_parametres
    use mod_tri_maillage
    use mod_maillage
    use mod_legendre
    use mod_types_derives

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_base :
!   - hat_phi_i
!   - hat_phi_F
!   - gamma_F
!   - make_Quad_ref
!   - make_Base_ref
!   - make_Base_phy
! ==================================================================================================== !

    contains

        function hat_phi_i(Mesh, i, xi, eta) result(hat_phi_i_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant l'application envoyant l'element de reference sur la maille physique

! Entrees de la fonction :
!   - Mesh : Type derive contenant les informations de maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i

!   - i : Entier referant a la maille Omega_i
!   - xi : Coordonnee x des elements de refrence
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction :
!   - hat_phi_i_xi_eta(1:2) : Vecteur (x, y) des coordonnes physiques
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: xi, eta

            type(Maillage), intent(in)                                  :: Mesh

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: hat_phi_i_xi_eta

! Variables locales
            integer                                                     :: s_k, s_l, s_n, s_m

            if (Mesh%N_sommets_Omega_i(i) == 3) then
                s_k = Mesh%Sommets_Omega_i(1, i) ; s_l = Mesh%Sommets_Omega_i(2, i)
                s_n = Mesh%Sommets_Omega_i(3, i)

                hat_phi_i_xi_eta = (1._pr - xi - eta)*Mesh%Coord_Sommets(:, s_k)                        &
            &                    + xi*Mesh%Coord_Sommets(:, s_l) + eta*Mesh%Coord_Sommets(:, s_n)

            else if (Mesh%N_sommets_Omega_i(i) == 4) then
                s_k = Mesh%Sommets_Omega_i(1, i) ; s_l = Mesh%Sommets_Omega_i(2, i)
                s_n = Mesh%Sommets_Omega_i(3, i) ; s_m = Mesh%Sommets_Omega_i(4, i)

                hat_phi_i_xi_eta = 1._pr/4*(1 - xi)*(1 - eta)*Mesh%Coord_Sommets(:, s_k)                &
            &                    + 1._pr/4*(1 + xi)*(1 - eta)*Mesh%Coord_Sommets(:, s_l)                &
            &                    + 1._pr/4*(1 + xi)*(1 + eta)*Mesh%Coord_Sommets(:, s_n)                &
            &                    + 1._pr/4*(1 - xi)*(1 + eta)*Mesh%Coord_Sommets(:, s_m)

            else
                print *, "hat_phi_i : Forme non prise en charge"
                return
            end if

        end function hat_phi_i

        function hat_phi_F(Mesh, i, F, xi) result(hat_phi_F_xi)
! ==================================================================================================== !
! Fonction renvoyant l'application envoyant la face de reference sur la fave physique

! Entrees de la fonction :
!   - Mesh : Type derive contenant les informations de maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i
!   - xi : Coordonnee x des elements de refrence

! Sortie de la fonction :
!   - hat_phi_F_xi(1:2) : Vecteur (x, y) d'un point appartenant a une face
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i, F

            real(kind = pr), intent(in)                                 :: xi

            type(Maillage), intent(in)                                  :: Mesh

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: hat_phi_F_xi

! Variables locales
            integer                                                     :: s_k, s_l

            call get_S_i_1_back(i, F, Mesh%N_sommets_Omega_i, Mesh%Sommets_Omega_i, s_k, s_l)

            hat_phi_F_xi = (1._pr - xi)/2*Mesh%Coord_Sommets(:, s_k)                                    &
        &                                                     + (1._pr + xi)/2*Mesh%Coord_Sommets(:, s_l)

        end function hat_phi_F

        function gamma_F(N_sommets, F, xi) result(gamma_F_xi)
! ==================================================================================================== !
! Fonction renvoyant l'application gamma_F etant la trace de la face F

! Entrees de la fonction :
!   - N_sommets : Nombre de sommets de la maille Omega_i
!   - F : Entier compris en 1 et N_max_sommets tel que F = 1 refere a la face 1 de la
! Omega_i, et ainsi de suite
!   - xi : Coordonnee 1D des elements de refrence

! Sortie de la fonction :
!   - gamma_F_xi : Trace de la face F de la maille Omega_i
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: N_sommets, F
            
            real(kind = pr), intent(in)                                 :: xi

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: gamma_F_xi

! Variables locales
            real(kind = pr), dimension(2)                               :: Coord_Ref_1, Coord_Ref_2

            Coord_Ref_1 = 0._pr ; Coord_Ref_2 = 0._pr

            if (N_sommets == 3) then
                if (F == 1) then
                    Coord_Ref_1 = [0._pr, 0._pr] ; Coord_Ref_2 = [1._pr, 0._pr]
                else if (F == 2) then
                    Coord_Ref_1 = [1._pr, 0._pr] ; Coord_Ref_2 = [0._pr, 1._pr]
                else if (F == 3) then
                    Coord_Ref_1 = [0._pr, 1._pr] ; Coord_Ref_2 = [0._pr, 0._pr]
                end if

            else if (N_sommets == 4) then
                if (F == 1) then
                    Coord_Ref_1 = [-1._pr, -1._pr] ; Coord_Ref_2 = [ 1._pr, -1._pr]
                else if (F == 2) then
                    Coord_Ref_1 = [ 1._pr, -1._pr] ; Coord_Ref_2 = [ 1._pr,  1._pr]
                else if (F == 3) then
                    Coord_Ref_1 = [ 1._pr,  1._pr] ; Coord_Ref_2 = [-1._pr,  1._pr]
                else if (F == 4) then
                    Coord_Ref_1 = [-1._pr,  1._pr] ; Coord_Ref_2 = [-1._pr, -1._pr]
                end if

            else
                print *, "gamma_F : Cas non pris en compte"
                return
            end if

            gamma_F_xi = (1._pr - xi)/2*Coord_Ref_1 + (1._pr + xi)/2*Coord_Ref_2

        end function gamma_F

        subroutine make_Quad_ref(racines_Legendre, w_j_Legendre, Quad_ref)
! ==================================================================================================== !
! Subroutine renvoyant les coordonnees des points de Gauss et leur poids associe sur l'element
! de reference

! Entrees de la subroutine :
!   - racines_Legendre(1:N_G) : Tableau contenant les N_G racines du polynome de Legendre 1D de
! de degre N_G
!   - w_j_Legendre(1:N_G) : Tableau contenant les poids de quadrature associes aux points de
! Gauss-Legendre

! Entree/Sortie de la subroutine :
!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - N_G_Face : Nombre de points de Gauss pour la face de reference
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - Coord_Ref_Face(1:2, 1:N_G_Face) : Coordonnees des points de Gauss pour la face de
! reference
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T
!   - w_j_Ref_Face(1:N_G_Face) : Poids d'interpolation associes aux points de Gauss pour la
! face de reference
! ==================================================================================================== !

! Entrees de la subroutine
            real(kind = pr), dimension(:), intent(in)                   :: racines_Legendre, w_j_Legendre

! Entree/Sortie de la subroutine
            type(Quadrature_hat_Omega), intent(inout)                   :: Quad_ref

! Variables locales
            integer                                                     :: i, j, N_point, N_Dunavant

            real(kind = pr), dimension(:), allocatable                  :: w_j_Dunavant
            real(kind = pr), dimension(:, :), allocatable               :: racines_Dunavant

            call quadrature_Dunavant(N_Dunavant, racines_Dunavant, w_j_Dunavant)

! Remplissages des tableaux a l'interieur des mailles
            Quad_ref%N_G_Quad = N_G**2
            Quad_ref%N_G_Tri = N_Dunavant

            allocate(Quad_ref%Coord_Ref_Quad(1:2, 1:Quad_ref%N_G_Quad), Quad_ref%w_j_Ref_Quad(1:Quad_ref%N_G_Quad))
            allocate(Quad_ref%Coord_Ref_Tri(1:2, 1:Quad_ref%N_G_Tri), Quad_ref%w_j_Ref_Tri(1:Quad_ref%N_G_Tri))

! Remplissage des tableaux pour hat_Omega = Q
            N_point = 0
            do j = 1, N_G
                do i = 1, N_G
                    N_point = N_point + 1

                    Quad_ref%Coord_Ref_Quad(1, N_point) = racines_Legendre(i)    ! xi
                    Quad_ref%Coord_Ref_Quad(2, N_point) = racines_Legendre(j)    ! eta

                    Quad_ref%w_j_Ref_Quad(N_point) = w_j_Legendre(i)*w_j_Legendre(j)
                end do
            end do

! Remplissage des tableaux pour hat_Omega = T
            Quad_ref%Coord_Ref_Tri = racines_Dunavant

            Quad_ref%w_j_Ref_Tri = w_j_Dunavant

            deallocate(w_j_Dunavant, racines_Dunavant)

! Remplissage des tableaux pour les faces
            Quad_ref%N_G_Face = N_G

            allocate(Quad_ref%Coord_Ref_Face(1:Quad_ref%N_G_Face), Quad_ref%w_j_Ref_Face(1:Quad_ref%N_G_Face))

            Quad_ref%Coord_Ref_Face = racines_Legendre
            Quad_ref%w_j_Ref_Face = w_j_Legendre

        end subroutine make_Quad_ref

        subroutine make_Base_ref(Mesh, Quad_ref, Base_ref)
! ==================================================================================================== !
! Subroutine construisant le vecteur de base et sa derivee sur l'element de reference Omega_i

! Entrees de la subroutine :
!   - Mesh : Type derive contenant les informations de maillage
!   - N_elements : Nombre de mailles
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - N_G_Face : Nombre de points de Gauss pour la face de reference
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - Coord_Ref_Face(1:2, 1:N_G_Face) : Coordonnees des points de Gauss pour la face de
! reference

! Entree/Sortie de la subroutine :
!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Quad : Nombre de DDL pour une maille quadrangulaire
!   - N_ddl_i_Tri : Nombre de DDL pour une maille triangulaire
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre
!   - Phi_F_Quad(0:N_ddl_i-1, 1:N_G_Face, 1:4) : Vecteur contenant la base de Legendre evaluee aux
! points de Gauss-Legendre des 4 faces F
!   - Phi_F_Tri(0:N_ddl_i-1, 1:N_G_Face, 1:3) : Vecteur contenant la base de Legendre evaluee aux
! points de Gauss-Legendre des 3 faces F
!   - Grad_xi_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la derivee par rapport a xi de la base
! de Legendre evaluee aux points de Gauss-Legendre
!   - Grad_eta_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la derivee par rapport a eta de la base
! de Legendre evaluee aux points de Gauss-Legendre
!   - Grad_xi_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la derivee par rapport a xi de la base
! de Dubiner evaluee aux points de Gauss-Legendre
!   - Grad_eta_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la derivee par rapport a eta de la base
! de Dubiner evaluee aux points de Gauss-Legendre
! ==================================================================================================== !

! Entrees de la subroutine
            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref

! Entree/Sortie de la subroutine
            type(Base_hat_Omega), intent(inout)                         :: Base_ref

! Variables locales
            integer                                                     :: N_point, N_base, i, k, m, F, j

            real(kind = pr)                                             :: xi, eta, c_k, c_m
            real(kind = pr), dimension(2)                               :: gamma_F_xi

            Base_ref%N_ddl_i_Quad = (p+1)**2
            Base_ref%N_ddl_i_Tri = ((p+1)*(p+2))/2

            allocate(Base_ref%N_ddl_i(1:Mesh%N_elements))

            Base_ref%N_ddl_i = 0

            do i = 1, Mesh%N_elements
                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    Base_ref%N_ddl_i(i) = Base_ref%N_ddl_i_Tri
                else if (Mesh%N_Sommets_Omega_i(i) == 4) then
                    Base_ref%N_ddl_i(i) = Base_ref%N_ddl_i_Quad
                end if
            end do

! Phi_i(0, :) represente la moyenne dans la maille, ie le degre polynomial 0

! Mailles rectangulaires
            allocate(Base_ref%Phi_i_Quad(0:Base_ref%N_ddl_i_Quad-1, 1:Quad_ref%N_G_Quad))
            allocate(Base_ref%Grad_xi_Quad(0:Base_ref%N_ddl_i_Quad-1, 1:Quad_ref%N_G_Quad))
            allocate(Base_ref%Grad_eta_Quad(0:Base_ref%N_ddl_i_Quad-1, 1:Quad_ref%N_G_Quad))

! Mailles triangulaires
            allocate(Base_ref%Phi_i_Tri(0:Base_ref%N_ddl_i_Tri-1, 1:Quad_ref%N_G_Tri))
            allocate(Base_ref%Grad_xi_Tri(0:Base_ref%N_ddl_i_Tri-1, 1:Quad_ref%N_G_Tri))
            allocate(Base_ref%Grad_eta_Tri(0:Base_ref%N_ddl_i_Tri-1, 1:Quad_ref%N_G_Tri))
        
! Interfaces des mailles
            allocate(Base_ref%Phi_F_Tri(0:Base_ref%N_ddl_i_Tri-1, 1:Quad_ref%N_G_Face, 1:3))
            allocate(Base_ref%Phi_F_Quad(0:Base_ref%N_ddl_i_Quad-1, 1:Quad_ref%N_G_Face, 1:4))

! Base de Legendre (hat_Omega = Q)
            do N_point = 1, Quad_ref%N_G_Quad
                xi = Quad_ref%Coord_Ref_Quad(1, N_point)
                eta = Quad_ref%Coord_Ref_Quad(2, N_point)

                N_base = -1
                do k = 0, p
                    do m = 0, p
                        N_base = N_base + 1
                        c_k = SQRT(real(2*k+1)/2._pr) ; c_m = SQRT(real(2*m+1)/2._pr)

                        Base_ref%Phi_i_Quad(N_base, N_point) = Legendre_2D(k, m, xi, eta)

                        Base_ref%Grad_xi_Quad(N_base, N_point) = c_k*Legendre_xi_1D(k, xi)*c_m*Legendre_1D(m, eta)
                        Base_ref%Grad_eta_Quad(N_base, N_point) = c_k*Legendre_1D(k, xi)*c_m*Legendre_xi_1D(m, eta)
                    end do
                end do
            end do

! Base de Dubiner (hat_Omega = T)
            do N_point = 1, Quad_ref%N_G_Tri
                xi = Quad_ref%Coord_Ref_Tri(1, N_point)
                eta = Quad_ref%Coord_Ref_Tri(2, N_point)

                N_base = -1
                do k = 0, p
                    do m = 0, p-k
                        N_base = N_base + 1

                        Base_ref%Phi_i_Tri(N_base, N_point) = Dubiner_2D(k, m, xi, eta)

                        Base_ref%Grad_xi_Tri(N_base, N_point) = Dubiner_xi_2D(k, m, xi, eta)
                        Base_ref%Grad_eta_Tri(N_base, N_point) = Dubiner_eta_2D(k, m, xi, eta)
                    end do
                end do
            end do

! Base de Legendre (hat_F)
            do F = 1, 4
                do j = 1, Quad_ref%N_G_Face
                    gamma_F_xi = gamma_F(4, F, Quad_ref%Coord_Ref_Face(j))

                    N_base = -1
                    do k = 0, p
                        do m = 0, p
                            N_base = N_base + 1
                            Base_ref%Phi_F_Quad(N_base, j, F) = Legendre_2D(k, m, gamma_F_xi(1), gamma_F_xi(2))
                        end do
                    end do
                end do
            end do

! Base de Dubiner (hat_F)
            do F = 1, 3
                do j = 1, Quad_ref%N_G_Face
                    gamma_F_xi = gamma_F(3, F, Quad_ref%Coord_Ref_Face(j))

                    N_base = -1
                    do k = 0, p
                        do m = 0, p-k
                            N_base = N_base + 1
                            Base_ref%Phi_F_Tri(N_base, j, F) = Dubiner_2D(k, m, gamma_F_xi(1), gamma_F_xi(2))
                        end do
                    end do
                end do
            end do

        end subroutine make_Base_ref

        subroutine make_Base_phy(Mesh, Base_ref, i, Base_phy, xi, eta, j, nabla_Tri, nabla_Quad, F)
! ==================================================================================================== !
! Subroutine construisant les operateurs de passage associes au changement de variable
! (x, y) = hat_phi_i(xi, eta)

! Entrees de la subroutine :
!   - Mesh : Type derive contenant les informations de maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Sommets_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des sommets de
! la maille Omega_i
!   - Coord_Sommets(1:2, 1:N_noeuds) : Tableau contenant les coordonnes dans l'espace physique
! du sommet s de Omega_i
!   - Long_arete(1:N_aretes) : Tableau contenant la longueur de l'arete a
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - Grad_xi_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la derivee par rapport a xi de la base
! de Legendre evaluee aux points de Gauss-Legendre
!   - Grad_eta_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la derivee par rapport a eta de la base
! de Legendre evaluee aux points de Gauss-Legendre
!   - Grad_xi_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la derivee par rapport a xi de la base
! de Dubiner evaluee aux points de Gauss-Legendre
!   - Grad_eta_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la derivee par rapport a eta de la base
! de Dubiner evaluee aux points de Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - j : Entier referant au j-eme point de Gauss
!   - F : Entier referant a la face F de la maille Omega_i
!   - xi : Coordonnee x des elements de refrence
!   - eta : Coordonnee y des elements de reference

! Entrees/Sorties de la subroutine :
!   - Base_phy : Type derive contenant les operateurs de passage issu du changement de variable de
! la maille physique vers la maille elementaire
!   - Inv_J_i_T(2, 2) : Transposee de l'inverse de la matrice Jacobienne associee a hat_phi_i
!   - Det_J_i : Dertiminant de la matrice Jacobienne associee a hat_phi_i
!   - Det_J_F : Determinant de la matrice Jacobienne associee a hat_phi_F
!   - nabla_Tri(0:N_ddl_i_Tri-1, 1:N_G_Tri, 1:2) : Operateur gradient issu du changement de 
! variable lorsque la maille elementaire est T
!   - nabla_Quad(0:N_ddl_i_Quad-1, 1:N_G_Quad, 1:2) : Operateur gradient issu du changement de 
! variable lorsque la maille elementaire est Q
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i
            integer, intent(in), optional                               :: j, F

            real(kind = pr), intent(in), optional                       :: xi, eta

            type(Maillage), intent(in)                                  :: Mesh
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entrees/Sorties de la fonction
            real(pr), dimension(:, :, :), intent(inout), optional       :: nabla_Tri
            real(pr), dimension(:, :, :), intent(inout), optional       :: nabla_Quad

            type(Base_Omega_i), intent(inout)                           :: Base_phy

! Variables locales
            integer                                                     :: s_k, s_l, s_n, s_m, a

            real(kind = pr), dimension(2, 2)                            :: J_i_xi_eta

            if (present(xi) .AND. present(eta)) then
                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    s_k = Mesh%Sommets_Omega_i(1, i) ; s_l = Mesh%Sommets_Omega_i(2, i) ; s_n = Mesh%Sommets_Omega_i(3, i)

                    J_i_xi_eta(1, 1) = Mesh%Coord_Sommets(1, s_l) - Mesh%Coord_Sommets(1, s_k)
                    J_i_xi_eta(1, 2) = Mesh%Coord_Sommets(1, s_n) - Mesh%Coord_Sommets(1, s_k)

                    J_i_xi_eta(2, 1) = Mesh%Coord_Sommets(2, s_l) - Mesh%Coord_Sommets(2, s_k)
                    J_i_xi_eta(2, 2) = Mesh%Coord_Sommets(2, s_n) - Mesh%Coord_Sommets(2, s_k)
                
                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    s_k = Mesh%Sommets_Omega_i(1, i) ; s_l = Mesh%Sommets_Omega_i(2, i)
                    s_n = Mesh%Sommets_Omega_i(3, i) ; s_m = Mesh%Sommets_Omega_i(4, i)

                    J_i_xi_eta(1, 1) = 1._pr/4*((1._pr - eta)*(Mesh%Coord_Sommets(1, s_l) - Mesh%Coord_Sommets(1, s_k))         &
                    &                         + (1._pr + eta)*(Mesh%Coord_Sommets(1, s_n) - Mesh%Coord_Sommets(1, s_m)))
                    J_i_xi_eta(1, 2) = 1._pr/4*((1._pr - xi)*(Mesh%Coord_Sommets(1, s_m) - Mesh%Coord_Sommets(1, s_k))          &
                    &                         + (1._pr + xi)*(Mesh%Coord_Sommets(1, s_n) - Mesh%Coord_Sommets(1, s_l)))

                    J_i_xi_eta(2, 1) = 1._pr/4*((1._pr - eta)*(Mesh%Coord_Sommets(2, s_l) - Mesh%Coord_Sommets(2, s_k))         &
                    &                         + (1._pr + eta)*(Mesh%Coord_Sommets(2, s_n) - Mesh%Coord_Sommets(2, s_m)))
                    J_i_xi_eta(2, 2) = 1._pr/4*((1._pr - xi)*(Mesh%Coord_Sommets(2, s_m) - Mesh%Coord_Sommets(2, s_k))          &
                    &                         + (1._pr + xi)*(Mesh%Coord_Sommets(2, s_n) - Mesh%Coord_Sommets(2, s_l)))

                else
                    print *, "make_J_i : Forme non prise en charge"
                    return
                end if

! Determinant de la matrice Jacobienne locale
                Base_phy%Det_J_i = J_i_xi_eta(1, 1)*J_i_xi_eta(2, 2) - J_i_xi_eta(2, 1)*J_i_xi_eta(1, 2)

! Transposee de l'inverse de la matrice Jacobienne locale
                Base_phy%Inv_J_i_T(1, 1) = 1._pr/Base_phy%Det_J_i*J_i_xi_eta(2, 2)
                Base_phy%Inv_J_i_T(1, 2) = -1._pr/Base_phy%Det_J_i*J_i_xi_eta(2, 1)

                Base_phy%Inv_J_i_T(2, 1) = -1._pr/Base_phy%Det_J_i*J_i_xi_eta(1, 2)
                Base_phy%Inv_J_i_T(2, 2) = 1._pr/Base_phy%Det_J_i*J_i_xi_eta(1, 1)

! Gradient numerique associe au changement de variable
                if (Mesh%N_sommets_Omega_i(i) == 3 .AND. present(nabla_Tri)) then
! Composante x pour une maille triangulaire
                    nabla_Tri(:, j, 1) = Base_phy%Inv_J_i_T(1, 1)*Base_ref%Grad_xi_Tri(:, j)                                &
                &                                                      + Base_phy%Inv_J_i_T(1, 2)*Base_ref%Grad_eta_Tri(:, j)

! Composante y pour une maille triangulaire
                    nabla_Tri(:, j, 2) = Base_phy%Inv_J_i_T(2, 1)*Base_ref%Grad_xi_Tri(:, j)                                &
                &                                                      + Base_phy%Inv_J_i_T(2, 2)*Base_ref%Grad_eta_Tri(:, j)

                else if (Mesh%N_sommets_Omega_i(i) == 4 .AND. present(nabla_Quad)) then
! Composante x pour une maille quadrangulaire
                    nabla_Quad(:, j, 1) = Base_phy%Inv_J_i_T(1, 1)*Base_ref%Grad_xi_Quad(:, j)                              &
                &                                                     + Base_phy%Inv_J_i_T(1, 2)*Base_ref%Grad_eta_Quad(:, j)

! Composante y pour une maillage quadrangulaire
                    nabla_Quad(:, j, 2) = Base_phy%Inv_J_i_T(2, 1)*Base_ref%Grad_xi_Quad(:, j)                              &
                &                                                     + Base_phy%Inv_J_i_T(2, 2)*Base_ref%Grad_eta_Quad(:, j)
                end if
            end if

! Determinant de la Jacobienne liee a la transformation vers la face de reference
            if (present(F)) then
                a = Mesh%Aretes_Omega_i(F, i)
                Base_phy%Det_J_F = Mesh%Long_arete(a)/2._pr
            end if

        end subroutine make_Base_phy

end module mod_base