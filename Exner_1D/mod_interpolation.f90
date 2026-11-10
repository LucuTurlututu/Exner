module mod_interpolation

    use mod_parametres
    use mod_fonctions
    use mod_legendre

    implicit none

    contains

        subroutine make_B_i(wi_legendre, Phi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine faisant la projection a l'instant initial de B_i sur la base
! de Legendre

! Entrees de la subroutine :
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids
!   - Phi_i(0:p, 1:N_G) : Matrice Phi_i
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss

! Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Sortie de la fonction
            real(kind = pr), dimension(0:p, 1:N_x), intent(out)         :: B_i

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr)                         :: interpolation_b_i

            do i = 1, N_x
                do k = 0, p
                    interpolation_b_i = 0._pr
                    do j = 1, N_G
                        interpolation_b_i = interpolation_b_i + b(x_gauss(j, i), 0._pr)*                &
                    &                                               Phi_i(k, j)*wi_legendre(j)
                    end do
                    B_i(k, i) = real(2*k+1, pr)/2*interpolation_b_i
                end do
            end do

        end subroutine make_B_i

        subroutine reconstruction_locale(i, B_i, b_L_i, b_R_i)
! ============================================================================ !
! Subroutine renvoyant la valeur de b_h(x_{i-1/2}^+, t) si xi = -1 et 
! b_h(x_{i+1/2}^-, t) si xi = 1

! Entrees de la subroutine :
!   - i : Entier concernant le numero i de la maille
!   - B_i(0:p, 1:N_x) : Matrice B_i

! Sorties de la subroutine :
!   - b_L_i(2) : Vecteur contentant les valeurs de la bathymetrie a gauche de
! l'interface des mailles i-1 et i
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie a droite de
! l'interface des mailles i-1 et i
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: B_i

! Sorties de la subroutine
            real(kind = pr), dimension(2), intent(out)                  :: b_L_i, b_R_i

! Variable locale
            integer                                 :: k

            b_L_i = 0._pr ; b_R_i = 0._pr
            
            if (i == 1) then
                select case (cas)

                case default
! Condition periodique
                    do k = 0, p
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x)*Li(k, 1._pr)
                        b_R_i(1) = b_R_i(1) + B_i(k, 1)*Li(k, -1._pr)

                        b_L_i(2) = b_L_i(2) + B_i(k, 1)*Li(k, 1._pr)
                        b_R_i(2) = b_R_i(2) + B_i(k, 2)*Li(k, -1._pr)
                    end do

                case (2, 8, 9, 10)
! Condition de Dirichlet
                    b_L_i(1) = b(0._pr, 0._pr)
                    do k = 0, p
                        b_R_i(1) = b_R_i(1) + B_i(k, 1)*Li(k, -1._pr)

                        b_L_i(2) = b_L_i(2) + B_i(k, 1)*Li(k, 1._pr)
                        b_R_i(2) = b_R_i(2) + B_i(k, 2)*Li(k, -1._pr)
                    end do

                end select

            else if (i == N_x) then
                select case (cas)

                case default
! Condition periodique
                    do k = 0, p
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x-1)*Li(k, 1._pr)
                        b_R_i(1) = b_R_i(1) + B_i(k, N_x)*Li(k, -1._pr)

                        b_L_i(2) = b_L_i(2) + B_i(k, N_x)*Li(k, 1._pr)
                        b_R_i(2) = b_R_i(2) + B_i(k, 1)*Li(k, -1._pr)
                    end do

                case (2, 8, 9, 10)
! Condition de Dirichlet
                    do k = 0, p
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x-1)*Li(k, 1._pr)
                        b_R_i(1) = b_R_i(1) + B_i(k, N_x)*Li(k, -1._pr)

                        b_L_i(2) = b_L_i(2) + B_i(k, N_x)*Li(k, 1._pr)
                    end do
                    b_R_i(2) = b(L, 0._pr)

                end select
            
            else
                do k = 0, p
                    b_L_i(1) = b_L_i(1) + B_i(k, i-1)*Li(k, 1._pr)
                    b_R_i(1) = b_R_i(1) + B_i(k, i)*Li(k, -1._pr)

                    b_L_i(2) = b_L_i(2) + B_i(k, i)*Li(k, 1._pr)
                    b_R_i(2) = b_R_i(2) + B_i(k, i+1)*Li(k, -1._pr)
                end do
            end if

        end subroutine reconstruction_locale

        subroutine make_Beta_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i_loc, b_L_i, b_R_i, Beta_i)
! ============================================================================ !
! Subroutine renvoyant le vecteur local de pente pour le terme b_x

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice Phi_xi_i contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss
!   - B_i_loc(0:p) : Colonne de B_i pour chaque point de discretisation
!   - b_L_i(2) : Vecteur contenant les valeurs de la bathymetrie b sur les
! interfaces gauche
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie b sur les
! interfaces droite
            
! Sortie de la subroutine :
!   - Beta_i(0:p) : Vecteur local de pente de b_x
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx
            real(kind = pr), dimension(2), intent(in)                   :: b_L_i, b_R_i
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p), intent(in)                 :: B_i_loc
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i

! Sortie de la subroutine
            real(kind = pr), dimension(0:p), intent(out)                :: Beta_i

! Variables locales
            integer                                 :: j, k, n
            real(kind = pr)                         :: b_hat_L, b_hat_R, b_h_loc
            real(kind = pr), dimension(0:p)         :: integrale_b
            
            integrale_b = 0._pr

            b_hat_L = (b_L_i(1) + b_R_i(1))/2
            b_hat_R = (b_L_i(2) + b_R_i(2))/2

            do j = 1, N_G
                b_h_loc = 0._pr
                do n = 0, p
                    b_h_loc = b_h_loc + B_i_loc(n)*Phi_i(n, j)
                end do
                
                do k = 0, p
                    integrale_b(k) = integrale_b(k) + b_h_loc*Phi_xi_i(k, j)*wi_legendre(j)
                end do
            end do

            do k = 0, p
                Beta_i(k) = real(2*k+1, pr)/dx*(-integrale_b(k)                                         &
        &                                   + b_hat_R*Li(k, 1._pr) - b_hat_L*Li(k, -1._pr))
            end do

        end subroutine make_Beta_i

        subroutine make_M_i(dx, wi_legendre, Phi_i, Beta_i, M_i)
! ============================================================================ !
! Subroutine fabriquant la matrice de masse locale M_i

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids
!   - Phi_i(0:p, 1:N_G) : Matrice Phi_i contenant la valeur des polynomes de 
! Legendre en chaque point de Gauss
!   - Beta_i(0:p) : Vecteur local de pente de b_x

! Sortie de la subroutine :
!   - M_i(0:p, 0:p) : Matrice de masse locale M_i a l'instant t
! ============================================================================ !

! Entrees de la subroutine :
            real(kind = pr), intent(in)             :: dx
            real(kind = pr), dimension(0:p), intent(in)                 :: Beta_i
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i

! Sortie de la subroutine
            real(kind = pr), dimension(0:p, 0:p), intent(out)           :: M_i

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr)                         :: b_x_loc, m_b_x, det_J
    
            det_J = dx/2._pr
            M_i = 0._pr

            do i = 1, N_G
                b_x_loc = 0._pr
                do k = 0, p
                    b_x_loc = b_x_loc + Beta_i(k)*Phi_i(k, i)
                end do

                m_b_x = m(b_x_loc)
                
                do k = 0, p
                    do j = 0, p
                        M_i(k, j) = M_i(k, j) + det_J*wi_legendre(i)*m_b_x*Phi_i(k, i)*Phi_i(j, i)
                    end do
                end do
            end do

        end subroutine make_M_i

        subroutine make_omega_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i, omega_i)
! ============================================================================ !
! Fonction renvoyant la masse omega_i associee a la maille i afin de verifier
! la bonne conservation de la quantite de masse

! Entrees de la subroutine
!   - dx : Pas d'espace
!   - wi_legendre(1:N_G) : Poids d'interpolation
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss
!   - B_i(0:p, 1:N_x) : Matrice de bathymetrie

! Sortie de la subroutine
!   - omega_i(1:N_x) : Vecteur des quantites de masse liees a la maille i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: B_i

! Sortie de la subroutine
            real(kind = pr), dimension(1:N_x), intent(out)              :: omega_i

! Variables locales
            integer                                 :: i, k
            real(kind = pr), dimension(2)           :: b_L_i, b_R_i
            real(kind = pr), dimension(0:p)         :: Beta_i
            real(kind = pr), dimension(0:p, 0:p)    :: M_i

            omega_i = 0._pr

            do i = 1, N_x
                call reconstruction_locale(i, B_i, b_L_i, b_R_i)
                
                call make_Beta_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i(:, i), b_L_i, b_R_i, Beta_i)

                call make_M_i(dx, wi_legendre, Phi_i, Beta_i, M_i)

                do k = 0, p
                    omega_i(i) = omega_i(i) + B_i(k, i)*M_i(0, k)
                end do
            end do

        end subroutine make_omega_i

end module mod_interpolation