module mod_legendre

    use mod_parametres

    implicit none

    contains

        pure function Li(k, x) result(Li_k_x)
! ============================================================================ !
! Fonction renvoyant le polynome de Legendre de degre k <= p

! Entrees de la fonction :
!   - k : degre du polynome de Legendre
!   - x : point de discretisation

! Sortie de la fonction :
!   - Li_k_x : Polynome de Legendre de degre k evalue au point x
! ============================================================================ !
            
! Entrees de la fonction
            integer, intent(in)                     :: k
            real(kind = pr), intent(in)             :: x

! Sortie de la fonction
            real(kind = pr)                         :: Li_k_x

! Variables locales
            integer                                 :: n
            real(kind = pr)                         :: Li_n, Li_nm1

            Li_k_x = 0._pr

            if (k == 0) then
                Li_k_x = 1._pr
            else if (k == 1) then
                Li_k_x = x
            else
! Formule de Bonnet pour la base de Legendre
                Li_nm1 = 1._pr ; Li_n = x
                do n = 1, k-1
                    Li_k_x = real(2*n+1, pr)/(n+1)*x*Li_n - real(n, pr)/(n+1)*Li_nm1
                    Li_nm1 = Li_n ; Li_n = Li_k_x
                end do
            end if

        end function Li

        pure function Li_x(k, x) result(Li_x_kx)
! ============================================================================ !
! Fonction renvoyant la derivee du polynome de Legendre de degre k <= p

! Entrees de la fonction :
!   - k : degre de la derivee du polynome de Legendre
!   - x : point de discretisation

! Sortie de la fonction :
!   - Li_x_kx : Derivee du polynome de Legendre de degre k evalue au point x
! ============================================================================ !

! Entrees de la fonction
            integer, intent(in)                     :: k
            real(kind = pr), intent(in)             :: x

! Sortie de la fonction
            real(kind = pr)                         :: Li_x_kx

! Variables locales
            integer                                 :: n
            real(kind = pr)                         :: Li_nm1, Li_x_nm2, Li_x_nm1

            Li_x_kx = 0._pr

            if (k == 0) then
                Li_x_kx = 0._pr
            else if (k == 1) then
                Li_x_kx = 1._pr
            else
                Li_x_nm2 = 0._pr ; Li_x_nm1 = 1._pr
                do n = 2, k
                    Li_nm1 = Li(n-1, x)
                    Li_x_kx = (2*n-1)*Li_nm1 + Li_x_nm2
                    Li_x_nm2 = Li_x_nm1 ; Li_x_nm1 = Li_x_kx
                end do
            end if

        end function Li_x

        subroutine racines_legendre(Li_racines)
! ============================================================================ !
! Subroutine renvoyant un tableau Li_racines(1:N_G) contenant les racines de
! chaque polynome de Legendre Li_n, n <= N_G

! Sortie de la subroutine :
!   Li_racines(1:N_G) : Tableau contenant les N_G racines de Li
! ============================================================================ !

! Sortie de la subroutine
            real(kind = pr), dimension(1:N_G), intent(out)              :: Li_racines

! Variables locales
            integer                                 :: m, j
            real(kind = pr)                         :: x_0, x_q, x_qp1, Li_q, Li_x_q

! On ne cherche que les racines positives par symetrie des polynomes de legendre
            m = INT((N_G + 1)/2._pr)

            do j = 1, m
! On initialise avec une approximation relativement bonne des racines
                x_0 = COS(pi*(j - 0.25_pr)/(N_G + 0.5_pr))

                x_q = 1._pr ; x_qp1 = x_0
                do while (ABS(x_qp1 - x_q) > eps_racine)
                    x_q = x_qp1
                    Li_q = Li(N_G, x_q)
                    Li_x_q = Li_x(N_G, x_q)

                    x_qp1 = x_q - Li_q/Li_x_q
                end do

                Li_racines(j) = -x_q ; Li_racines(N_G-j+1) = x_q
            end do

        end subroutine racines_legendre

        subroutine poids_legendre(Li_racines, wi_legendre)
! ============================================================================ !
! Subroutine renvoyant un tableau wi_legendre(1:N_G) contenant les poids de
! d'interpolation

! Entree de la subroutine
!   Li_racines(1:N_G) : Tableau contenant les N_G racines de Li

! Sortie de la subroutine :
!   wi_legendre(1:N_G) : Tableau contenant les N_G poids
! ============================================================================ !

! Entree de la subroutine
            real(kind = pr), dimension(1:N_G), intent(in)           :: Li_racines

! Sortie de la subroutine
            real(kind = pr), dimension(1:N_G), intent(out)          :: wi_legendre

! Variables locales
            integer                                                 :: j
            real(kind = pr)                                         :: x_j, Li_x_xj

            do j = 1, N_G
                x_j = Li_racines(j) ; Li_x_xj = Li_x(N_G, x_j)
                wi_legendre(j) = 2/((1 - x_j**2)*Li_x_xj**2)
            end do

        end subroutine poids_legendre

        subroutine points_gauss(dx, Li_racines, x_gauss)
! ============================================================================ !
! Subroutine mappant l'intervalle de reference vers les intervalles reels du
! maillage

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - Li_racines(1:N_G) : Tableau contenant les N_G racines de Li

! Sortie de la subroutine :
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx
            real(kind = pr), dimension(1:N_G), intent(in)               :: Li_racines

! Sortie de la subroutine
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(out)       :: x_gauss

! Variables locales
            integer                                 :: i, j
            real(kind = pr)                         :: x_mdemi

            do i = 1, N_x
                x_mdemi = (i-1)*dx
                do j = 1, N_G
                    x_gauss(j, i) = x_mdemi + (Li_racines(j) + 1)*dx/2
                end do
            end do

        end subroutine points_gauss

        subroutine make_Phi_i(Li_racines, Phi_i)
! ============================================================================ !
! Subroutine creant la matrice Phi_i contenant la valeur des polynomes de
! Legendre en chaque point de Gauss

! Entree de la subroutine :
!   - Li_racines(1:N_G) : Tableau contenant les N_G racines de Li

! Sortie de la subroutine :
!   - Phi_i(0:p, 1:N_G) : Matrice Phi_i
! ============================================================================ !

! Entree de la subroutine
            real(kind = pr), dimension(1: N_G), intent(in)              :: Li_racines

! Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_G), intent(out)         :: Phi_i

! Variables locales
            integer                                 :: j, k

            do j = 1, N_G
                do k = 0, p
                    Phi_i(k, j) = Li(k, Li_racines(j))
                end do
            end do

        end subroutine make_Phi_i

        subroutine make_Phi_xi_i(Li_racines, Phi_xi_i)
! ============================================================================ !
! Subroutine creant la matrice Phi_xi_i contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree de la subroutine :
!   - Li_racines(1:N_G) : Tableau contenant les N_G racines de Li

! Sortie de la subroutine :
!   - Phi_xi_i(0:p, 1:N_G) : Matrice Phi_xi_i
! ============================================================================ !

! Entree de la subroutine
            real(kind = pr), dimension(1: N_G), intent(in)              :: Li_racines

! Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_G), intent(out)         :: Phi_xi_i

! Variables locales
            integer                                 :: j, k

            do j = 1, N_G
                do k = 0, p
                    Phi_xi_i(k, j) = Li_x(k, Li_racines(j))
                end do
            end do

        end subroutine make_Phi_xi_i

end module mod_legendre