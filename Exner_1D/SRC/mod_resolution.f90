module mod_resolution

    use mod_parametres
    use mod_fonctions
    use mod_legendre
    use mod_interpolation
    use mod_limiteur

    implicit none

    contains

        subroutine flux_rusanov(i, dx, t, b_L_i, b_R_i, F_pdemi)
! ============================================================================ !
! Subroutine calculant le flux de Rusanov F_{i+1/2} pour une maille i

! Entrees de la subroutine :
!   - i : Entier concernant la maille i
!   - dx : Pas d'espace
!   - t : Temps courant
!   - b_L_i(2) : Vecteur contentant les valeurs de la bathymetrie a gauche de
! l'interface des mailles i-1 et i
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie a droite de
! l'interface des mailles i-1 et i

! Sortie de la subroutine :
!   - F_pdemi(2) : Flux de Rusanov aux interfaces d'une maille i
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: i
            real(kind = pr), intent(in)             :: dx, t
            real(kind = pr), dimension(2), intent(in)                   :: b_L_i, b_R_i

! Sortie de la subroutine
            real(kind = pr), dimension(2), intent(out)                  :: F_pdemi

! Variables locales
            integer                                 :: n
            real(kind = pr)                         :: x_int, h_L, h_R, u_L, u_R, tau_L, tau_R, alpha_pdemi
            real(kind = pr), dimension(2)           :: F_L, F_R

            F_pdemi = 0._pr
            
            do n = 1, 2
                x_int = (i-(2-n))*dx
                h_L = h(x_int, t) ; u_L = u(x_int, t) ; tau_L = tau(u_L)
                h_R = h(x_int, t) ; u_R = u(x_int, t) ; tau_R = tau(u_R)

                F_L(n) = F(b_L_i(n), h_L, u_L, tau_L)
                F_R(n) = F(b_R_i(n), h_R, u_R, tau_R)

                alpha_pdemi = MAX(ABS(F_b(b_L_i(n), h_L, u_L, tau_L)), ABS(F_b(b_R_i(n), h_R, u_R, tau_R)))

                F_pdemi(n) = 1._pr/2*((F_L(n) + F_R(n)) - alpha_pdemi*(b_R_i(n) - b_L_i(n)))
            end do

        end subroutine flux_rusanov

        subroutine make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i, R_i)
! ============================================================================ !
! Subroutine constituant le vecteur R_i tel que
!       M_i*B_t = R_i

! Entrees de la subroutine :
!   - i : Entier referant a la maille i
!   - dx : Pas d'espace
!   - t : Temps courant
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss
!   - B_i(0:p, 1:N_x) : Matrice B_i

! Sortie de la subroutine :
!   - R_i(0:p) : Vecteur second membre R_i
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: i
            real(kind = pr), intent(in)             :: dx, t
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: B_i

! Sorties de la subroutine
            real(kind = pr), dimension(0:p), intent(out)                :: R_i

! Variables locales
            integer                                 :: j, k
            real(kind = pr)                         :: h_h, u_h, tau_h
            real(kind = pr), dimension(2)           :: b_L_i, b_R_i, F_pdemi
            real(kind = pr), dimension(1:N_G)       :: b_h_gauss, F_gauss, S_gauss, R_tilde_gauss

            call reconstruction_locale(i, B_i, b_L_i, b_R_i)
            call flux_rusanov(i, dx, t, b_L_i, b_R_i, F_pdemi)

            b_h_gauss = MATMUL(TRANSPOSE(Phi_i), B_i(:, i))

            do j = 1, N_G
                h_h = h(x_gauss(j, i), t)
                u_h = u(x_gauss(j, i), t)
                tau_h = tau(u_h)

                F_gauss(j) = F(b_h_gauss(j), h_h, u_h, tau_h)*wi_legendre(j)
                S_gauss(j) = S(b_h_gauss(j), h_h, u_h, tau_h)*wi_legendre(j)
                R_tilde_gauss(j) = R_tilde(x_gauss(j, i), t)*wi_legendre(j)
            end do

            R_i = MATMUL(Phi_xi_i, F_gauss)

            do k = 0, p
                R_i(k) = R_i(k) - (F_pdemi(2)*Li(k, 1._pr) - F_pdemi(1)*Li(k, -1._pr))                  &
            &                                 + dx/2*(DOT_PRODUCT(S_gauss + R_tilde_gauss, Phi_i(k, :)))
            end do

        end subroutine make_R_i

        subroutine Cholesky(M_i, R_i_loc, X_i)
! ============================================================================ !
! Subroutine effectuant la decomposition de Cholesky de la matrice M_i pour
! obtenir le vecteur X_i solution du probleme M_i*X_i = R_i_loc

! Entrees de la subroutine
!   - M_i(0:p, 0:p) : Matrice de masse locale
!   - R_i_loc(0:p) : Vecteur second membre correspondant a la composante i
! du RHS

! Sortie de la subroutine
!   - X_i : Solution du probleme M_i*X_i = R_i_loc
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), dimension(0:p), intent(in)                 :: R_i_loc
            real(kind = pr), dimension(0:p, 0:p), intent(in)            :: M_i

! Sortie de la subroutine
            real(kind = pr), dimension(0:p), intent(out)                :: X_i

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr)                         :: z
            real(kind = pr), dimension(0:p)         :: Y
            real(kind = pr), dimension(0:p, 0:p)    :: L

            L = 0._pr
            do i = 0, p
                do j = 0, i
                    z = 0._pr
                    do k = 0, j - 1
                        z = z + L(i, k)*L(j, k)
                    end do

                    if (i == j) then
                        L(i, i) = SQRT(MAX(M_i(i, i) - z, eps_var))
                    else
                        L(i, j) = (M_i(i, j) - z)/L(j, j)
                    end if
                end do
            end do

! Descente
            do i = 0, p
                z = 0._pr
                do j = 0, i - 1
                    z = z + L(i, j)*Y(j)
                end do
                Y(i) = (R_i_loc(i) - z)/L(i, i)
            end do

! Remontee
            do i = p, 0, -1
                z = 0._pr
                do j = i + 1, p
                    z = z + L(j, i)*X_i(j)
                end do
                X_i(i) = (Y(i) - z)/L(i, i)
            end do

        end subroutine Cholesky

        subroutine make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i, b_R_i, R_i, B_t_i)
! ============================================================================ !
! Subroutine construisant la derivee de B par rapport au temps telle que
! B_t_i = M_i^{-1}*R_i

! Entrees de la subroutine :
!   - i : Entier referant a la maille i
!   - dx : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss
!   - B_i(0:p, 1:N_x) : Matrice B_i
!   - b_L_i(2) : Vecteur contentant les valeurs de la bathymetrie a gauche de
! l'interface des mailles i-1 et i
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie a droite de
! l'interface des mailles i-1 et i
!   - R_i(0:p) : Vecteur second membre R_i

! Sortie de la subroutine :
!   - B_t_i(0:p) : Vecteur B_t tel que B_t = M_i^{-1}*R_i
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: i
            real(kind = pr), intent(in)             :: dx
            real(kind = pr), dimension(2), intent(in)                   :: b_L_i, b_R_i
            real(kind = pr), dimension(0:p), intent(in)                 :: R_i
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: B_i

! Sorties de la subroutine
            real(kind = pr), dimension(0:p), intent(out)                :: B_t_i

! Variables locales
            real(kind = pr), dimension(0:p)         :: Beta_i
            real(kind = pr), dimension(0:p, 0:p)    :: M_i

            call make_Beta_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i(:, i), b_L_i, b_R_i, Beta_i)

            call make_M_i(dx, wi_legendre, Phi_i, Beta_i, M_i)

            call Cholesky(M_i, R_i, B_t_i)

        end subroutine make_B_t_i

        subroutine find_dt(dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i, dt)
! ============================================================================ !
! Subroutine renvoyant le pas de temps optimal module par le facteur de surete
! CFL

! Entrees de la subroutine
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice Phi_i
!   - Phi_xi_i(0:p, 1:N_G) : Matrice Phi_xi_i
!   - B_i(0:p, 1:N_x) : Bathymetrie au temps actuel

! Sortie de la subroutine
!   - dt : Pas de temps optimal
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: B_i

! Sortie de la subroutine
            real(kind = pr), intent(out)            :: dt

! Variables locales
            integer                                 :: i, j, k, n
            real(kind = pr)                         :: lambda_i, mu_i, h_ij, u_ij, tau_ij, inv_m_b_x_h, &
        &   S_h, S_h_eps, mu_ij, x_int, b_h_int, b_x_h_L, b_x_h_R, b_x_h_int, h_int, u_int, tau_int,    &
        &   lambda_tilde, lambda_int, lambda_max, mu_max, dt_adv, dt_sou
            real(kind = pr), dimension(2)           :: b_L_i, b_R_i
            real(kind = pr), dimension(0:p)         :: Beta_i
            real(kind = pr), dimension(1:N_G)       :: b_h_gauss, b_x_h_gauss

            lambda_i = 0._pr ; mu_i = 0._pr
            lambda_max = 0._pr ; mu_max = 0._pr

            do i = 1, N_x
                call reconstruction_locale(i, B_i, b_L_i, b_R_i)

                call make_Beta_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i(:, i), b_L_i, b_R_i, Beta_i)

                b_h_gauss = MATMUL(TRANSPOSE(Phi_i), B_i(:, i))
                b_x_h_gauss = MATMUL(TRANSPOSE(Phi_i), Beta_i)

! Pas de temps lie a l'advection
                b_x_h_L = 0._pr ; b_x_h_R = 0._pr
                do k = 0, p
                    b_x_h_L = b_x_h_L + Beta_i(k)*Li(k, -1._pr)
                    b_x_h_R = b_x_h_R + Beta_i(k)*Li(k, 1._pr)
                end do

                do n = 1, 2
                    if (n == 1) then
! On est sur l'interface de gauche de la maille
                        x_int = real(i-1._pr/2, pr)*dx
                        b_h_int = b_R_i(1)
                        b_x_h_int = b_x_h_L
                    else
! On est sur l'interface de droite
                        x_int = real(i+1._pr/2, pr)*dx
                        b_h_int = b_L_i(2)
                        b_x_h_int = b_x_h_R
                    end if

                    h_int = h(x_int, t)
                    u_int = u(x_int, t)
                    tau_int = tau(u_int)

                    inv_m_b_x_h = 1/m(b_x_h_int)

! Terme lie aux solutions manufacturees
                    lambda_tilde = 0._pr
                    if (1 <= cas .AND. cas <= 6) then
                        S_h = S(b_h_int, h_int, u_int, tau_int) + R_tilde(x_int, t)
                        if (eps_var < b_x_h_int) then
                            lambda_tilde = ABS(S_h)/ABS(b_x_h_int + eps_var)
                        end if
                    end if

                    select case (cas)

                    case (11)
                        lambda_int = MAX(ABS(u_int - SQRT(g_part*h_int)), ABS(u_int + SQRT(g_part*h_int)),                  &
                    &                                                   inv_m_b_x_h*ABS(F_b(b_h_int, h_int, u_int, tau_int)))

                    case default
                        lambda_int = inv_m_b_x_h*ABS(F_b(b_h_int, h_int, u_int, tau_int))

                    end select
                    lambda_i = MAX(lambda_i, lambda_int + lambda_tilde)
                end do

! Pas de temps lie au terme source
                do j = 1, N_G
                    h_ij = h(x_gauss(j, i), t)
                    u_ij = u(x_gauss(j, i), t)
                    tau_ij = tau(u_ij)

                    inv_m_b_x_h = 1/m(b_x_h_gauss(j))

                    S_h = S(b_h_gauss(j), h_ij, u_ij, tau_ij)
                    S_h_eps = S(b_h_gauss(j) + eps_dt, h_ij, u_ij, tau_ij)

                    mu_ij = inv_m_b_x_h*ABS((S_h_eps - S_h)/eps_dt)
                    mu_i = MAX(mu_i, mu_ij)
                end do
                lambda_max = MAX(lambda_max, lambda_i)
                mu_max = MAX(mu_max, mu_i)
            end do

            dt_adv = t_final ; dt_sou = t_final
            if (eps_var < lambda_max) then
                dt_adv = dx/((2*p+1)*lambda_max)
            end if

            if (eps_var < mu_max) then
                dt_sou = 1._pr/mu_max
            end if

            dt = CFL*MIN(alpha_RK*dt_adv, dt_sou)
            
        end subroutine find_dt

        subroutine RK1(dx, t, dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i)
! ============================================================================ !
! Subroutine pour le schema RK1 (ie Euler Explicite)

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps t^n
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
!   - b_L(0:N_x) : Vecteur contenant b_h(x_{i-1/2}^+, t) pour toutes les
! mailles
!   - b_R(0:N_x) : Vecteur contenant b_h(x_{i+1/2}^-, t) pour toutes les
! mailles

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i, b_R_i
            real(kind = pr), dimension(0:p)         :: R_i, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_t

            do i = 1, N_x
                call reconstruction_locale(i, B_i, b_L_i, b_R_i)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i, R_i)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i, b_R_i, R_i, B_t_i)

                B_t(:, i) = B_t_i
            end do

            B_i = B_i + dt*B_t

        end subroutine RK1

        subroutine RK3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema RK3

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_n, B_n_1, B_n_2, B_t_n, B_t_1, B_t_2

! 1e etape
            B_n = B_i
            do i = 1, N_x
                call reconstruction_locale(i, B_n, b_L_i_n, b_R_i_n)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t_n(:, i) = B_t_i
            end do
            B_n_1 = B_n + 1._pr/2*dt*B_t_n

            select case (cas)

            case (7)
                call krivodonova_2007(t, B_n_1)

            end select

! 2e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_1, b_L_i_1, b_R_i_1)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t_1(:, i) = B_t_i
            end do
            B_n_2 = B_n + dt*(-B_t_n + 2._pr*B_t_1)

            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_n_2)

            end select

! 3e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_2, b_L_i_2, b_R_i_2)

                call make_R_i(i, dx, t + dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t_2(:, i) = B_t_i
            end do
            B_i = B_n + 1._pr/6*dt*(B_t_n + 4._pr*B_t_1 + B_t_2)
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt, B_i)

            end select

        end subroutine RK3

        subroutine SSPRK_3_3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema SSPRK(3, 3)

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_t, B_n, B_n_1, B_n_2

! 1e etape
            B_n = B_i
            do i = 1, N_x
                call reconstruction_locale(i, B_n, b_L_i_n, b_R_i_n)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_1 = B_n + dt*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t, B_n_1)

            end select

! 2e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_1, b_L_i_1, b_R_i_1)

                call make_R_i(i, dx, t + dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_2 = 3._pr/4*B_n + 1._pr/4*(B_n_1 + dt*B_t)

            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_n_2)

            end select

! 3e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_2, b_L_i_2, b_R_i_2)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_i = 1._pr/3*B_n + 2._pr/3*(B_n_2 + dt*B_t)
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt, B_i)

            end select

        end subroutine SSPRK_3_3

        subroutine SSPRK_4_3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema SSPRK(4, 3)

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2, b_L_i_3, b_R_i_3
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, R_i_3, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_t, B_n, B_n_1, B_n_2, B_n_3

! 1e etape
            B_n = B_i
            do i = 1, N_x
                call reconstruction_locale(i, B_n, b_L_i_n, b_R_i_n)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_1 = B_n + dt/2*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t, B_n_1)

            end select

! 2e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_1, b_L_i_1, b_R_i_1)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_2 = B_n_1 + 1._pr/2*dt*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_n_2)

            end select

! 3e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_2, b_L_i_2, b_R_i_2)

                call make_R_i(i, dx, t + dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_3 = 2._pr/3*B_n + 1._pr/3*B_n_2 + 1._pr/6*dt*B_t
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_n_3)

            end select

! 4e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_3, b_L_i_3, b_R_i_3)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_3, R_i_3)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_3, b_L_i_3, b_R_i_3, R_i_3, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_i = B_n_3 + 1._pr/2*dt*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t + dt, B_i)

            end select

        end subroutine SSPRK_4_3

        subroutine RK4(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema RK4

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2, b_L_i_3, b_R_i_3
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, R_i_3, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_n, B_n_1, B_n_2, B_n_3, B_t_n, B_t_1, B_t_2, B_t_3

! 1e etape
            B_n = B_i
            do i = 1, N_x
                call reconstruction_locale(i, B_n, b_L_i_n, b_R_i_n)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t_n(:, i) = B_t_i
            end do
            B_n_1 = B_n + 1._pr/2*dt*B_t_n

            select case (cas)

            case (7)
                call krivodonova_2007(t, B_n_1)

            end select

! 2e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_1, b_L_i_1, b_R_i_1)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t_1(:, i) = B_t_i
            end do
            B_n_2 = B_n + 1._pr/2*dt*B_t_1

            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_n_2)

            end select

! 3e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_2, b_L_i_2, b_R_i_2)

                call make_R_i(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t_2(:, i) = B_t_i
            end do
            B_n_3 = B_n + dt*B_t_2
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt/2, B_i)

            end select
    
! 4e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_3, b_L_i_3, b_R_i_3)

                call make_R_i(i, dx, t + dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_3, R_i_3)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_3, b_L_i_3, b_R_i_3, R_i_3, B_t_i)

                B_t_3(:, i) = B_t_i
            end do
            B_i = B_n + 1._pr/6*dt*(B_t_n + 2._pr*B_t_1 + 2._pr*B_t_2 + B_t_3)
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt, B_i)

            end select

        end subroutine RK4

        subroutine SSPRK_5_4(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema SSPRK(5, 4)

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2,                &
        &                                                                               b_L_i_3, b_R_i_3, b_L_i_4, b_R_i_4
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, R_i_3, R_i_4, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_t, B_n, B_n_1, B_n_2, B_n_3, B_n_4, B_t_4

! 1e etape
            B_n = B_i
            do i = 1, N_x
                call reconstruction_locale(i, B_n, b_L_i_n, b_R_i_n)

                call make_R_i(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_1 = B_n + 0.39175222700392_pr*dt*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t + 0.39175222700392_pr*dt, B_n_1)

            end select

! 2e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_1, b_L_i_1, b_R_i_1)

                call make_R_i(i, dx, t + 0.39175222700392_pr*dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_2 = 0.44437049406734_pr*B_n + 0.55562950593266_pr*B_n_1 + 0.36841059262959_pr*dt*B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t + 0.58607968896779_pr*dt, B_n_2)

            end select

! 3e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_2, b_L_i_2, b_R_i_2)

                call make_R_i(i, dx, t + 0.58607968896779_pr*dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_3 = 0.62010185138540_pr*B_n + 0.37989814861460_pr*B_n_2 + 0.25189177424738_pr*dt*B_t
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + 0.47454236302687_pr*dt, B_n_3)

            end select

! 4e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_3, b_L_i_3, b_R_i_3)

                call make_R_i(i, dx, t + 0.47454236302687_pr*dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_3, R_i_3)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_3, b_L_i_3, b_R_i_3, R_i_3, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_4 = 0.17807995410773_pr*B_n + 0.82192004589227_pr*B_n_3 + 0.54497475021237_pr*dt*B_t
            B_t_4 = B_t

            select case (cas)

            case (7)
                call krivodonova_2007(t + 0.93501063100924_pr*dt, B_n_4)

            end select

! 5e etape
            do i = 1, N_x
                call reconstruction_locale(i, B_n_4, b_L_i_4, b_R_i_4)

                call make_R_i(i, dx, t + 0.93501063100924_pr*dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_n_4, R_i_4)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_4, b_L_i_4, b_R_i_4, R_i_4, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_i = 0.00683325884039_pr*B_n + 0.51723167208978_pr*B_n_2 + 0.12759831133288_pr*B_n_3 +                         &
        &                              0.34833675773694_pr*B_n_4 + 0.08460416338212_pr*dt*B_t_4 + 0.22600748319395_pr*dt*B_t
            
            select case (cas)

            case (7)
                call krivodonova_2007(t + dt, B_i)

            end select

        end subroutine SSPRK_5_4

end module mod_resolution