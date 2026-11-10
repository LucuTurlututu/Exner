module mod_cl

    use mod_parametres
    use mod_qte_physiques
    use mod_tri_maillage
    use mod_types_derives
    use mod_base

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_cl :
!   - Cholesky
!   - get_Coord_Sommets_ghost
!   - hat_phi_i_ghost
!   - Dirichlet
!   - Neumann
!   - Periodique
! ==================================================================================================== !

    contains

        subroutine Cholesky(Base_ref, i, A_i, b_i, X_i)
! ==================================================================================================== !
! Subroutine implementant l'algorithme de Cholesky

! Entrees de la subroutine :
!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i

!   - i : Entier rferant a la maille Omega_i
!   - A_i(0:N_ddl_i-1, 0:N_ddl_i-1) : Bloc i d'une matrice A
!   - b_i(0:N_ddl_i-1) : Composante i d'un vecteur b

! Sortie de la subroutine :
!   - X_i(0:N_ddl_i-1) : Vecteur solution du probleme A_i*X_i = b_i
! ==================================================================================================== !

! Entrees de la subroutine :
            type(Base_hat_Omega), intent(in)                            :: Base_ref

            integer, intent(in)                                         :: i

            real(pr), dimension(0:), intent(in)                         :: b_i
            real(pr), dimension(0:, 0:), intent(in)                     :: A_i

! Sortie de la subroutine
            real(pr), dimension(0:Base_ref%N_ddl_i(i)-1), intent(out)   :: X_i

! Variables locales
            integer                                                     :: k, m, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            BLOCK
                real(kind = pr)                                         :: z
                real(kind = pr), dimension(0:N_ddl_i-1)                 :: Y, X_temp
                real(kind = pr), dimension(0:N_ddl_i-1, 0:N_ddl_i-1)    :: L

                L = 0._pr

                do k = 0, N_ddl_i-1
                    do m = 0, k
                        if (k == m) then
                            if (m == 0) then
                                z = A_i(k, k)
                            else
                                z = A_i(k, k) - DOT_PRODUCT(L(k, 0:m-1), L(k, 0:m-1))
                            end if
                            
                            L(k, k) = SQRT(MAX(eps_var, z))

                        else
                            if (m == 0) then
                                L(k, m) = A_i(k, m)/L(m, m)
                            else
                                L(k, m) = (A_i(k, m) - DOT_PRODUCT(L(k, 0:m-1), L(m, 0:m-1)))/L(m, m)
                            end if

                        end if
                    end do
                end do

! Descente
                do k = 0, N_ddl_i-1
                    if (k == 0) then
                        Y(k) = b_i(k)/L(k, k)
                    else
                        Y(k) = (b_i(k) - DOT_PRODUCT(L(k, 0:k-1), Y(0:k-1)))/L(k, k)
                    end if
                end do

! Remontee
                do k = N_ddl_i-1, 0, -1
                    if (k == N_ddl_i-1) then
                        X_temp(k) = Y(k)/L(k, k)
                    else
                        X_temp(k) = (Y(k) - DOT_PRODUCT(L(k+1:N_ddl_i-1, k), X_temp(k+1:N_ddl_i-1)))/L(k, k)
                    end if
                end do
                X_i = X_temp
            END BLOCK

        end subroutine Cholesky

        function get_Coord_Sommets_ghost(Mesh, i, F) result(Coord_Sommets_ghost)
! ==================================================================================================== !
! Fonction renvoyant les coordonnees des sommets de la maille fantome Omega_r, adjacente
! a Omega_i par la face F

! Entrees de la fonction :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i

! Sortie de la fonction :
!   - Coord_Sommets_ghost(1:2, 1:N_Sommets_Omega_i(i)) : Tableau contenant les coordonnees
! des sommets de la maille fantome Omega_r
! ==================================================================================================== !

! Entrees de la fonction :
            integer, intent(in)                                         :: i, F

            type(Maillage), intent(in)                                  :: Mesh

! Sortie de la fonction :
            real(pr), dimension(1:2, 1:Mesh%N_Sommets_Omega_i(i))       :: Coord_Sommets_ghost

! Variables locales
            integer                                                     :: N_sommets, s, s_k

            real(kind = pr)                                             :: d_e
            real(kind = pr), dimension(2)                               :: Coord_s_k, Coord_s_1

            N_sommets = Mesh%N_sommets_Omega_i(i)

            s_k = Mesh%Sommets_Omega_i(F, i)
            Coord_s_k = Mesh%Coord_Sommets(:, s_k)

            do s = 1, N_sommets
                Coord_s_1 = Mesh%Coord_Sommets(:, Mesh%Sommets_Omega_i(s, i))

                d_e = SUM((Coord_s_1 - Coord_s_k)*Mesh%n_F(:, F, i))

                Coord_Sommets_ghost(:, s) = Coord_s_1 - 2._pr*d_e*Mesh%n_F(:, F, i)
            end do

        end function get_Coord_Sommets_ghost

        function hat_phi_i_ghost(N_sommets, Coord_Sommets_ghost, xi, eta) result(hat_phi_i_ghost_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant l'application envoyant l'element de reference sur la maille physique

! Entrees de la fonction :
!   - N_Sommets : Nombre de sommets de la maille Omega_i
!   - Coord_Sommets_ghost(1:2, 1:N_Sommets_Omega_i(i)) : Tableau contenant les coordonnees
! des sommets de la maille fantome Omega_r
!   - xi : Coordonnee x des elements de refrence
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction :
!   - hat_phi_i_ghost_xi_eta(1:2) : Vecteur (x, y) des coordonnes physiques de la maille fantome
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: N_sommets

            real(kind = pr), intent(in)                                 :: xi, eta
            real(kind = pr), dimension(1:2, 1:N_sommets)                :: Coord_Sommets_ghost

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: hat_phi_i_ghost_xi_eta

            if (N_Sommets == 3) then
                hat_phi_i_ghost_xi_eta = (1._pr - xi - eta)*Coord_Sommets_ghost(:, 1)                   &
            &                          + xi*Coord_Sommets_ghost(:, 2) + eta*Coord_Sommets_ghost(:, 3)

            else if (N_sommets == 4) then
                hat_phi_i_ghost_xi_eta = 1._pr/4*(1 - xi)*(1 - eta)*Coord_Sommets_ghost(:, 1)           &
            &                          + 1._pr/4*(1 + xi)*(1 - eta)*Coord_Sommets_ghost(:, 2)           &
            &                          + 1._pr/4*(1 + xi)*(1 + eta)*Coord_Sommets_ghost(:, 3)           &
            &                          + 1._pr/4*(1 - xi)*(1 + eta)*Coord_Sommets_ghost(:, 4)

            else
                print *, "hat_phi_i_ghost : Forme non prise en charge"
                return
            end if

        end function hat_phi_i_ghost

        function Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t) result(B_r)
! ==================================================================================================== !
! Fonction imposant des conditions aux limites de Dirichlet

! Entrees de la fonction :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i
!   - t : Temps courant

! Sortie de la fonction :
!   - B_r(0:N_ddl_i-1) : Vecteur de modes de la bathymetrie dans la maille Omega_r
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i, F

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Calcul), intent(in)                                    :: Calc_Exner

! Sortie de la fonction
            real(kind = pr), dimension(0:Base_ref%N_ddl_i(i)-1)         :: B_r

! Variables locales
            integer                                                     :: j, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            B_r = 0._pr

            BLOCK
                real(kind = pr)                                         :: xi_j, eta_j
                real(kind = pr), dimension(2)                           :: Coord_phy_ghost
                real(pr), dimension(1:2, 1:Mesh%N_sommets_Omega_i(i))   :: Coord_Sommets_ghost

                type(Base_Omega_i)                                      :: Base_phy

                Coord_Sommets_ghost = get_Coord_Sommets_ghost(Mesh, i, F)

                B_r = 0._pr

                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)

                        Coord_phy_ghost = hat_phi_i_ghost(Mesh%N_sommets_Omega_i(i), Coord_Sommets_ghost, xi_j, eta_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        B_r(0:N_ddl_i-1) = B_r(0:N_ddl_i-1) + Quad_ref%w_j_Ref_Tri(j)*b(Coord_phy_ghost, t)                 &
                    &                                                         *Base_ref%Phi_i_Tri(:, j)*ABS(Base_phy%Det_J_i)
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)

                        Coord_phy_ghost = hat_phi_i_ghost(Mesh%N_sommets_Omega_i(i), Coord_Sommets_ghost, xi_j, eta_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        B_r(0:N_ddl_i-1) = B_r(0:N_ddl_i-1) + Quad_ref%w_j_Ref_Quad(j)*b(Coord_phy_ghost, t)                &
                    &                                                        *Base_ref%Phi_i_Quad(:, j)*ABS(Base_phy%Det_J_i)
                    end do
                end if

                call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i_std, B_r(0:N_ddl_i-1), B_r(0:N_ddl_i-1))

            END BLOCK

        end function Dirichlet

        function Neumann(Base_ref, Calc_Exner, i) result(B_r)
! ==================================================================================================== !
! Fonction imposant des conditions aux limites de Neumann

! Entrees de la fonction :
!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner

!   - i : Entier referant a la maille Omega_i

! Sortie de la fonction :
!   - B_r(0:N_ddl_i-1) : Vecteur de modes de la bathymetrie dans la maille Omega_r
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i

            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Calcul), intent(in)                                    :: Calc_Exner

! Sortie de la fonction
            real(kind = pr), dimension(0:Base_ref%N_ddl_i(i)-1)         :: B_r

            B_r = Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i)

        end function Neumann

        function Periodique(Mesh, Base_ref, Calc_Exner, i, F) result(B_r)
! ==================================================================================================== !
! Fonction imposant des conditions aux limites periodiques
    
! Entrees de la fonction :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i

! Sortie de la fonction :
!   - B_r(0:N_ddl_i-1) : Vecteur de modes de la bathymetrie dans la maille Omega_r
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: i, F

            type(Maillage), intent(in)                                  :: Mesh
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Calcul), intent(in)                                    :: Calc_Exner

! Sortie de la fonction
            real(kind = pr), dimension(0:Base_ref%N_ddl_i(i)-1)         :: B_r

! Variables locales
            integer                                                     :: r

            r = Mesh%N_Omega_r_P(F, i)

            B_r = Calc_Exner%B_i(0:Base_ref%N_ddl_i(r)-1, r)

        end function Periodique

end module mod_cl
