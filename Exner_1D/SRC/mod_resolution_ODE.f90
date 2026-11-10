module mod_resolution_ODE

! ============================================================================ !
! Module permettant la resolution de systemes d'EDO pour avoir une solution
! semi-analytique du probleme couple Shallow-Water-Exner

! Contient les routines :
!   - flux_rusanov_ODE
!   - make_R_i_ODE
!   - RK3_ODE
! Les equivalents des routines usuelles mais acceptant en plus un couple
! de grandeurs (h, u) calcule numeriquement
! ============================================================================ !

    use mod_parametres
    use mod_fonctions
    use mod_legendre
    use mod_interpolation
    use mod_limiteur
    use mod_resolution

    implicit none

    contains

        subroutine make_h_u_b_ini(wi_legendre, x_gauss, Phi_i, h_ini, u_ini, b_ini)
! ============================================================================ !
! Subroutine construisant les vecteurs de hauteur d'eau, vitesse et bathymetrie
! initiale pour la solution semi-analytique

! Entrees de la subroutine :
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss

! Sorties de la subroutine :
!   - h_ini(0:p, 1:N_x) : Vecteur (de modes) contenant la hauteur d'eau 
! a l'instant initial
!   - u_ini(0:p, 1:N_x) : Vecteur (de modes) contenant la vitesse a 
! l'instant initial
!   - b_ini(0:p, 1:N_x) : Vecteur contenant (de modes) la bathymetrie 
! a l'instant initial
! ============================================================================ !

! Entree de la subroutine
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Sorties de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(out)         :: h_ini, u_ini, b_ini

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr)                         :: int_h, int_u, int_b

            do i = 1, N_x
                do k = 0, p
                    int_h = 0._pr ; int_u = 0._pr ; int_b = 0._pr
                    do j = 1, N_G
                        int_h = int_h + h(x_gauss(j, i), 0._pr)*Phi_i(k, j)*wi_legendre(j)
                        int_u = int_u + u(x_gauss(j, i), 0._pr)*Phi_i(k, j)*wi_legendre(j)
                        int_b = int_b + b(x_gauss(j, i), 0._pr)*Phi_i(k, j)*wi_legendre(j)
                    end do
                    h_ini(k, i) = real(2*k+1, pr)/2*int_h
                    u_ini(k, i) = real(2*k+1, pr)/2*int_u
                    b_ini(k, i) = real(2*k+1, pr)/2*int_b
                end do
            end do

        end subroutine make_h_u_b_ini

        subroutine make_h_u_b_tilde(s, dt, wi_legendre, Phi_i, h_tilde_n, u_tilde_n, b_tilde_n,                             &
    &                                                h_tilde_sm1, u_tilde_sm1, b_tilde_sm1, h_tilde_s, u_tilde_s, b_tilde_s)
! ============================================================================ !
! Subroutine resolvant un systeme couple permettant d'obtenir le triplet
! (h_tilde, u_tilde, b_tilde) tel que
!   (h_semi, u_semi, b_semi) = (h_ini*h_tilde, u_ini*u_tilde, b_ini*b_tilde)

! Entrees de la subroutine :
!   - s : Numero de l'etape de RK3_ODE
!   - dt : Pas de temps
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - h_tilde_n(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau h_tilde 
! au temps t^n
!   - u_tilde_n(0:p, 1:N_x) : Vecteur contenant la vitesse u_tilde 
! au temps t^n
!   - b_tilde_n(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique
! b_tilde au temps t^n
!   - h_tilde_sm1(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau h_tilde 
! a l'etape s-1 de RK3_ODE
!   - u_tilde_sm1(0:p, 1:N_x) : Vecteur contenant la vitesse u_tilde 
! a l'etape s-1 de RK3_ODE
!   - b_tilde_sm1(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique 
! b_tilde a l'etape s-1 de RK3_ODE

! Sorties de la subroutine :
!   - h_tilde_s(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau h_tilde
! a l'etape s de RK3_ODE
!   - u_tilde_s(0:p, 1:N_x) : Vecteur contenant la vitesse u_tilde 
! a l'etape s de RK3_ODE
!   - b_tilde_s(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique 
! b_tilde a l'etape s de RK3_ODE
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: s
            real(kind = pr), intent(in)             :: dt
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_tilde_n, u_tilde_n, b_tilde_n
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_tilde_sm1, u_tilde_sm1, b_tilde_sm1

! Sorties de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(out)         :: h_tilde_s, u_tilde_s, b_tilde_s

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr), dimension(0:p)         :: h_tilde_t, u_tilde_t, b_tilde_t
            real(kind = pr), dimension(1:N_G)       :: h_gauss_sm1, u_gauss_sm1, b_gauss_sm1, h_gauss_t, u_gauss_t, b_gauss_t

            do i = 1, N_x
! ============================================================================ !
! A cause de la non-lineatite de l'EDO on doit passer par les points de Gauss
! pour resoudre le systeme
! ============================================================================ !
! Projection sur les points de Gauss
                h_gauss_sm1 = MATMUL(TRANSPOSE(Phi_i), h_tilde_sm1(:, i))
                u_gauss_sm1 = MATMUL(TRANSPOSE(Phi_i), u_tilde_sm1(:, i))
                b_gauss_sm1 = MATMUL(TRANSPOSE(Phi_i), b_tilde_sm1(:, i))

                do j = 1, N_G
                    h_gauss_t(j) = -a_FN*u_gauss_sm1(j)*h_gauss_sm1(j)
                    u_gauss_t(j) = -a_FN*u_gauss_sm1(j)**2 - 6._pr*g_part*(1._pr/(1._pr - p_d))*A_Grass*                   &
                &                                                                            a_FN**2/c_FN*b_gauss_sm1(j)
                    b_gauss_t(j) = -c_FN*u_gauss_sm1(j)**3
                end do

! Retour a l'espace DG
                h_tilde_t = 0._pr ; u_tilde_t = 0._pr ; b_tilde_t = 0._pr
                do k = 0, p
                    do j = 1, N_G
                        h_tilde_t(k) = h_tilde_t(k) + h_gauss_t(j)*Phi_i(k, j)*wi_legendre(j)
                        u_tilde_t(k) = u_tilde_t(k) + u_gauss_t(j)*Phi_i(k, j)*wi_legendre(j)
                        b_tilde_t(k) = b_tilde_t(k) + b_gauss_t(j)*Phi_i(k, j)*wi_legendre(j)
                    end do
                    h_tilde_t(k) = real(2*k+1, pr)/2*h_tilde_t(k)
                    u_tilde_t(k) = real(2*k+1, pr)/2*u_tilde_t(k)
                    b_tilde_t(k) = real(2*k+1, pr)/2*b_tilde_t(k)
                end do

! ============================================================================ !
! Une fois que l'on a resolu le systeme d'EDO, on peut simplement appliquer
! la transformation de l'etape s du schema SSP-RK3
! ============================================================================ !
                select case (s)

                case (1)
                    h_tilde_s(:, i) = h_tilde_n(:, i) + dt*h_tilde_t
                    u_tilde_s(:, i) = u_tilde_n(:, i) + dt*u_tilde_t
                    b_tilde_s(:, i) = b_tilde_n(:, i) + dt*b_tilde_t

                case (2)
                    h_tilde_s(:, i) = 3._pr/4*h_tilde_n(:, i) + 1._pr/4*(h_tilde_sm1(:, i) + dt*h_tilde_t)
                    u_tilde_s(:, i) = 3._pr/4*u_tilde_n(:, i) + 1._pr/4*(u_tilde_sm1(:, i) + dt*u_tilde_t)
                    b_tilde_s(:, i) = 3._pr/4*b_tilde_n(:, i) + 1._pr/4*(b_tilde_sm1(:, i) + dt*b_tilde_t)

                case (3)
                    h_tilde_s(:, i) = 1._pr/3*h_tilde_n(:, i) + 2._pr/3*(h_tilde_sm1(:, i) + dt*h_tilde_t)
                    u_tilde_s(:, i) = 1._pr/3*u_tilde_n(:, i) + 2._pr/3*(u_tilde_sm1(:, i) + dt*u_tilde_t)
                    b_tilde_s(:, i) = 1._pr/3*b_tilde_n(:, i) + 2._pr/3*(b_tilde_sm1(:, i) + dt*b_tilde_t)

                end select
            end do

        end subroutine make_h_u_b_tilde

        subroutine make_h_u_b_semi(wi_legendre, Phi_i, h_ini, u_ini, b_ini, h_tilde, u_tilde, b_tilde,  &
    &                                                                            h_semi, u_semi, b_semi)
! ============================================================================ !
! Subroutine construisant le triplet (h_semi, u_semi, b_semi) de solution
! semi-analytique du cas test 10

! Entrees de la subroutine :
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - h_ini(0:p, 1:N_x) : Vecteur (de modes) contenant la hauteur d'eau 
! a l'instant initial
!   - u_ini(0:p, 1:N_x) : Vecteur (de modes) contenant la vitesse a 
! l'instant initial
!   - b_ini(0:p, 1:N_x) : Vecteur contenant (de modes) la bathymetrie 
! a l'instant initial
!   - h_tilde(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau h_tilde 
!   - u_tilde(0:p, 1:N_x) : Vecteur contenant la vitesse u_tilde 
!   - b_tilde(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique

! Sorties de la subroutine :
!   - h_semi(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau semi-analytique
! en tous les points du domaine
!   - u_semi(0:p, 1:N_x) : Vecteur contenant la vitesse semi-analytique
! en tous les points du domaine
!   - b_semi(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique
! en tous les points du domaine
! ============================================================================ !

! Entrees de la subroutine :
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_ini, u_ini, b_ini
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_tilde, u_tilde, b_tilde

! Sorties de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(out)         :: h_semi, u_semi, b_semi

! Variables locales
            integer                                 :: i, k
            real(kind = pr), dimension(1:N_G)       :: h_ini_gauss, u_ini_gauss, b_ini_gauss
            real(kind = pr), dimension(1:N_G)       :: h_h_gauss, u_h_gauss, b_h_gauss

            h_semi = 0._pr ; u_semi = 0._pr ; b_semi = 0._pr

            do i = 1, N_x
                h_ini_gauss = MATMUL(TRANSPOSE(Phi_i), h_ini(:, i))
                u_ini_gauss = MATMUL(TRANSPOSE(Phi_i), u_ini(:, i))
                b_ini_gauss = MATMUL(TRANSPOSE(Phi_i), b_ini(:, i))

                h_h_gauss = MATMUL(TRANSPOSE(Phi_i), h_tilde(:, i))
                u_h_gauss = MATMUL(TRANSPOSE(Phi_i), u_tilde(:, i))
                b_h_gauss = MATMUL(TRANSPOSE(Phi_i), b_tilde(:, i))

                do k = 0, p
                    h_semi(k, i) = real(2*k+1, pr)/2*SUM(h_ini_gauss*h_h_gauss*Phi_i(k, :)*wi_legendre)
                    u_semi(k, i) = real(2*k+1, pr)/2*SUM(u_ini_gauss*u_h_gauss*Phi_i(k, :)*wi_legendre)
                    b_semi(k, i) = real(2*k+1, pr)/2*SUM(b_ini_gauss*b_h_gauss*Phi_i(k, :)*wi_legendre)
                end do
            end do

        end subroutine make_h_u_b_semi

        subroutine reconstruction_locale_ODE(i, h_semi, u_semi, b_semi, B_i, h_L_i, u_L_i, b_L_i, h_R_i, u_R_i, b_R_i)
! ============================================================================ !
! Subroutine renvoyant la valeur de X_h(x_{i-1/2}^+, t) si xi = -1 et 
! X_h(x_{i+1/2}^-, t) si xi = 1 où X = {h, u, b}

! Entrees de la subroutine :
!   - i : Entier referant a la maille i
!   - h_semi(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau semi-analytique
! en tous les points du domaine
!   - u_semi(0:p, 1:N_x) : Vecteur contenant la vitesse semi-analytique
! en tous les points du domaine
!   - b_semi(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique
! en tous les points du domaine
!   - B_i(0:p, 1:N_x) : Vecteur de bathymetrie (calcule numeriquement et non
! avec le systeme d'ODE !)

! Sorties de la subroutine :
!   - h_L_i(2) : Vecteur contentant les valeurs de la hauteur d'eau a gauche de
! l'interface des mailles i-1 et i
!   - u_L_i(2) : Vecteur contentant les valeurs de la vitesse a gauche de
! l'interface des mailles i-1 et i
!   - b_L_i(2) : Vecteur contentant les valeurs de la bathymetrie a gauche de
! l'interface des mailles i-1 et i
!   - h_R_i(2) : Vecteur contentant les valeurs de la hauteur d'eau a gauche de
! l'interface des mailles i-1 et i
!   - u_R_i(2) : Vecteur contentant les valeurs de la vitesse a gauche de
! l'interface des mailles i-1 et i
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie a droite de
! l'interface des mailles i-1 et i
! ============================================================================ !

! Entrees de la subroutine
            integer, intent(in)                     :: i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_semi, u_semi, b_semi, B_i

! Sorties de la subroutine
            real(kind = pr), dimension(2), intent(out)                  :: h_L_i, u_L_i, b_L_i
            real(kind = pr), dimension(2), intent(out)                  :: h_R_i, u_R_i, b_R_i

! Variables locales
            integer                                 :: k

            h_L_i = 0._pr ; h_R_i = 0._pr
            u_L_i = 0._pr ; u_R_i = 0._pr
            b_L_i = 0._pr ; b_R_i = 0._pr
            
            if (i == 1) then
                select case (cas)

                case default
! Condition periodique
                    do k = 0, p
                        h_L_i(1) = h_L_i(1) + h_semi(k, N_x)*Li(k, 1._pr) ; h_R_i(1) = h_R_i(1) + h_semi(k, 1)*Li(k, -1._pr)
                        u_L_i(1) = u_L_i(1) + u_semi(k, N_x)*Li(k, 1._pr) ; u_R_i(1) = u_R_i(1) + u_semi(k, 1)*Li(k, -1._pr)
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x)*Li(k, 1._pr) ; b_R_i(1) = b_R_i(1) + B_i(k, 1)*Li(k, -1._pr)

                        h_L_i(2) = h_L_i(2) + h_semi(k, 1)*Li(k, 1._pr) ; h_R_i(2) = h_R_i(2) + h_semi(k, 2)*Li(k, -1._pr)
                        u_L_i(2) = u_L_i(2) + u_semi(k, 1)*Li(k, 1._pr) ; u_R_i(2) = u_R_i(2) + u_semi(k, 2)*Li(k, -1._pr)
                        b_L_i(2) = b_L_i(2) + B_i(k, 1)*Li(k, 1._pr) ; b_R_i(2) = b_R_i(2) + B_i(k, 2)*Li(k, -1._pr)
                    end do

                case (11)
! Condition de Dirichlet
                    h_L_i(1) = 0._pr ; u_L_i(1) = 0._pr ; b_L_i(1) = 0._pr
                    do k = 0, p
                        h_L_i(1) = h_L_i(1) + h_semi(k, 1)*Li(k, -1._pr)
                        u_L_i(1) = u_L_i(1) + u_semi(k, 1)*Li(k, -1._pr)
                        b_L_i(1) = b_L_i(1) + b_semi(k, 1)*Li(k, -1._pr)
                    end do

                    do k = 0, p
                        h_R_i(1) = h_R_i(1) + h_semi(k, 1)*Li(k, -1._pr)
                        u_R_i(1) = u_R_i(1) + u_semi(k, 1)*Li(k, -1._pr)
                        b_R_i(1) = b_R_i(1) + B_i(k, 1)*Li(k, -1._pr)

                        h_L_i(2) = h_L_i(2) + h_semi(k, 1)*Li(k, 1._pr) ; h_R_i(2) = h_R_i(2) + h_semi(k, 2)*Li(k, -1._pr)
                        u_L_i(2) = u_L_i(2) + u_semi(k, 1)*Li(k, 1._pr) ; u_R_i(2) = u_R_i(2) + u_semi(k, 2)*Li(k, -1._pr)
                        b_L_i(2) = b_L_i(2) + B_i(k, 1)*Li(k, 1._pr) ; b_R_i(2) = b_R_i(2) + B_i(k, 2)*Li(k, -1._pr)
                    end do

                end select

            else if (i == N_x) then
                select case (cas)

                case default
! Condition periodique
                    do k = 0, p
                        h_L_i(1) = h_L_i(1) + h_semi(k, N_x-1)*Li(k, 1._pr) ; h_R_i(1) = h_R_i(1) + h_semi(k, N_x)*Li(k, -1._pr)
                        u_L_i(1) = u_L_i(1) + u_semi(k, N_x-1)*Li(k, 1._pr) ; u_R_i(1) = u_R_i(1) + u_semi(k, N_x)*Li(k, -1._pr)
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x-1)*Li(k, 1._pr) ; b_R_i(1) = b_R_i(1) + B_i(k, N_x)*Li(k, -1._pr)

                        h_L_i(2) = h_L_i(2) + h_semi(k, N_x)*Li(k, 1._pr) ; h_R_i(2) = h_R_i(2) + h_semi(k, 1)*Li(k, -1._pr)
                        u_L_i(2) = u_L_i(2) + u_semi(k, N_x)*Li(k, 1._pr) ; u_R_i(2) = u_R_i(2) + u_semi(k, 1)*Li(k, -1._pr)
                        b_L_i(2) = b_L_i(2) + B_i(k, N_x)*Li(k, 1._pr) ; b_R_i(2) = b_R_i(2) + B_i(k, 1)*Li(k, -1._pr)
                    end do

                case (11)
! Condition de Dirichlet ou sortie libre
                    do k = 0, p
                        h_L_i(1) = h_L_i(1) + h_semi(k, N_x-1)*Li(k, 1._pr) ; h_R_i(1) = h_R_i(1) + h_semi(k, N_x)*Li(k, -1._pr)
                        u_L_i(1) = u_L_i(1) + u_semi(k, N_x-1)*Li(k, 1._pr) ; u_R_i(1) = u_R_i(1) + u_semi(k, N_x)*Li(k, -1._pr)
                        b_L_i(1) = b_L_i(1) + B_i(k, N_x-1)*Li(k, 1._pr) ; b_R_i(1) = b_R_i(1) + B_i(k, N_x)*Li(k, -1._pr)

                        h_L_i(2) = h_L_i(2) + h_semi(k, N_x)*Li(k, 1._pr)
                        u_L_i(2) = u_L_i(2) + u_semi(k, N_x)*Li(k, 1._pr)
                        b_L_i(2) = b_L_i(2) + B_i(k, N_x)*Li(k, 1._pr)
                    end do

! Condition de sortie libre
                    h_R_i(2) = h_L_i(2) ; u_R_i(2) = u_L_i(2) ; b_R_i(2) = b_L_i(2)

! Condition de Dirichlet
                    ! do k = 0, p
                    !     h_R_i(2) = h_R_i(2) + h_semi(k, N_x)*Li(k, 1._pr)
                    !     u_R_i(2) = u_R_i(2) + u_semi(k, N_x)*Li(k, 1._pr)
                    !     b_R_i(2) = b_R_i(2) + b_semi(k, N_x)*Li(k, 1._pr)
                    ! end do

                end select
            
            else
                do k = 0, p
                    h_L_i(1) = h_L_i(1) + h_semi(k, i-1)*Li(k, 1._pr) ; h_R_i(1) = h_R_i(1) + h_semi(k, i)*Li(k, -1._pr)
                    u_L_i(1) = u_L_i(1) + u_semi(k, i-1)*Li(k, 1._pr) ; u_R_i(1) = u_R_i(1) + u_semi(k, i)*Li(k, -1._pr)
                    b_L_i(1) = b_L_i(1) + B_i(k, i-1)*Li(k, 1._pr) ; b_R_i(1) = b_R_i(1) + B_i(k, i)*Li(k, -1._pr)

                    h_L_i(2) = h_L_i(2) + h_semi(k, i)*Li(k, 1._pr) ; h_R_i(2) = h_R_i(2) + h_semi(k, i+1)*Li(k, -1._pr)
                    u_L_i(2) = u_L_i(2) + u_semi(k, i)*Li(k, 1._pr) ; u_R_i(2) = u_R_i(2) + u_semi(k, i+1)*Li(k, -1._pr)
                    b_L_i(2) = b_L_i(2) + B_i(k, i)*Li(k, 1._pr) ; b_R_i(2) = b_R_i(2) + B_i(k, i+1)*Li(k, -1._pr)
                end do
            end if

        end subroutine reconstruction_locale_ODE

        subroutine flux_rusanov_ODE(h_L_i, u_L_i, b_L_i, h_R_i, u_R_i, b_R_i, F_pdemi)
! ============================================================================ !
! Subroutine calculant le flux de Rusanov F_{i+1/2} pour une maille i

! Entrees de la subroutine :
!   - h_L_i(2) : Vecteur contentant les valeurs de la hauteur d'eau a gauche de
! l'interface des mailles i-1 et i
!   - u_L_i(2) : Vecteur contentant les valeurs de la vitesse a gauche de
! l'interface des mailles i-1 et i
!   - b_L_i(2) : Vecteur contentant les valeurs de la bathymetrie a gauche de
! l'interface des mailles i-1 et i
!   - h_R_i(2) : Vecteur contentant les valeurs de la hauteur d'eau a gauche de
! l'interface des mailles i-1 et i
!   - u_R_i(2) : Vecteur contentant les valeurs de la vitesse a gauche de
! l'interface des mailles i-1 et i
!   - b_R_i(2) : Vecteur contenant les valeurs de la bathymetrie a droite de
! l'interface des mailles i-1 et i

! Sortie de la subroutine :
!   - F_pdemi(2) : Flux de Rusanov aux interfaces d'une maille i
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), dimension(2), intent(in)                   :: h_L_i, u_L_i, b_L_i
            real(kind = pr), dimension(2), intent(in)                   :: h_R_i, u_R_i, b_R_i

! Sortie de la subroutine
            real(kind = pr), dimension(2), intent(out)                  :: F_pdemi

! Variables locales
            integer                                 :: n
            real(kind = pr)                         :: h_L, h_R, u_L, u_R, tau_L, tau_R, alpha_pdemi
            real(kind = pr), dimension(2)           :: F_L, F_R

            F_pdemi = 0._pr
            
            do n = 1, 2
                h_L = h_L_i(n) ; u_L = u_L_i(n) ; tau_L = tau(u_L)
                h_R = h_R_i(n) ; u_R = u_R_i(n) ; tau_R = tau(u_R)

                F_L(n) = F(b_L_i(n), h_L, u_L, tau_L)
                F_R(n) = F(b_R_i(n), h_R, u_R, tau_R)

                alpha_pdemi = MAX(ABS(F_b(b_L_i(n), h_L, u_L, tau_L)), ABS(F_b(b_R_i(n), h_R, u_R, tau_R)))

                F_pdemi(n) = 1._pr/2*((F_L(n) + F_R(n)) - alpha_pdemi*(b_R_i(n) - b_L_i(n)))
            end do

        end subroutine flux_rusanov_ODE

        subroutine make_R_i_ODE(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, h_semi, u_semi, b_semi, B_i, R_i)
! ============================================================================ !
! Subroutine constituant le vecteur R_i tel que
!       M_i*B_t = R_i

! Entrees de la subroutine :
!   - i : Entier referant a la maille i
!   - dx : Pas d'espace
!   - t : Temps courant
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss
!   - h_semi(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau semi-analytique
! en tous les points du domaine
!   - u_semi(0:p, 1:N_x) : Vecteur contenant la vitesse semi-analytique
! en tous les points du domaine
!   - b_semi(0:p, 1:N_x) : Vecteur contenant la bathymetrie semi-analytique
! en tous les points du domaine
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
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_semi, u_semi, b_semi, B_i

! Sorties de la subroutine
            real(kind = pr), dimension(0:p), intent(out)                :: R_i

! Variables locales
            integer                                 :: j, k
            real(kind = pr)                         :: tau_h
            real(kind = pr), dimension(2)           :: h_L_i, u_L_i, b_L_i, h_R_i, u_R_i, b_R_i, F_pdemi
            real(kind = pr), dimension(1:N_G)       :: h_h_gauss, u_h_gauss, b_h_gauss, F_gauss, S_gauss, R_tilde_gauss

            call reconstruction_locale_ODE(i, h_semi, u_semi, b_semi, B_i, h_L_i, u_L_i, b_L_i, h_R_i, u_R_i, b_R_i)
            call flux_rusanov_ODE(h_L_i, u_L_i, b_L_i, h_R_i, u_R_i, b_R_i, F_pdemi)

            h_h_gauss = MATMUL(TRANSPOSE(Phi_i), h_semi(:, i))
            u_h_gauss = MATMUL(TRANSPOSE(Phi_i), u_semi(:, i))
            b_h_gauss = MATMUL(TRANSPOSE(Phi_i), B_i(:, i))

            do j = 1, N_G
                tau_h = tau(u_h_gauss(j))

                F_gauss(j) = F(b_h_gauss(j), h_h_gauss(j), u_h_gauss(j), tau_h)*wi_legendre(j)
                S_gauss(j) = S(b_h_gauss(j), h_h_gauss(j), u_h_gauss(j), tau_h)*wi_legendre(j)
                R_tilde_gauss(j) = R_tilde(x_gauss(j, i), t)*wi_legendre(j)
            end do

            R_i = MATMUL(Phi_xi_i, F_gauss)

            do k = 0, p
                R_i(k) = R_i(k) - (F_pdemi(2)*Li(k, 1._pr) - F_pdemi(1)*Li(k, -1._pr))                  &
            &                                 + dx/2*(DOT_PRODUCT(S_gauss + R_tilde_gauss, Phi_i(k, :)))
            end do

        end subroutine make_R_i_ODE

        subroutine SSPRK_3_3_ODE(dx, t, dt, h_ini, u_ini, b_ini, h_tilde, u_tilde, b_tilde, wi_legendre,&
    &                                                                     Phi_i, Phi_xi_i, x_gauss, B_i)
! ============================================================================ !
! Subroutine pour le schema SSPRK(3,3)

! Entrees de la subroutine :
!   - dx : Pas d'espace
!   - t : Temps actuel
!   - dt : Pas de temps
!   - h_ini(1:N_x) : Vecteur contenant la hauteur d'eau a l'instant initial
!   - u_ini(1:N_x) : Vecteur contenant la vitesse a l'instant initial
!   - b_ini(1:N_x) : Vecteur contenant la bathymetrie a l'instant initial
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - Phi_xi_i(0:p, 1:N_G) : Matrice contenant la valeur de la derivee des
! polynomes de Legendre en chaque point de Gauss

! Entrees/Sorties de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
!   - h_tilde(1:N_x) : Vecteur contenant la hauteur d'eau semi-analytique
! en tous les points du domaine au temps t^n en entree, t^{n+1} en sortie
!   - u_tilde(1:N_x) : Vecteur contenant la vitesse semi-analytique
! en tous les points du domaine au temps t^n en entree, t^{n+1} en sortie
!   - b_tilde(1:N_x) : Vecteur contenant la bathymetrie semi-analytique
! en tous les points du domaine au temps t^n en entree, t^{n+1} en sortie
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), intent(in)             :: dx, t, dt
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_ini, u_ini, b_ini
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i, Phi_xi_i
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss

! Entrees/Sorties de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: h_tilde, u_tilde, b_tilde
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i
            real(kind = pr), dimension(2)           :: h_L_i_n, h_R_i_n, h_L_i_1, h_R_i_1, h_L_i_2, h_R_i_2
            real(kind = pr), dimension(2)           :: u_L_i_n, u_R_i_n, u_L_i_1, u_R_i_1, u_L_i_2, u_R_i_2
            real(kind = pr), dimension(2)           :: b_L_i_n, b_R_i_n, b_L_i_1, b_R_i_1, b_L_i_2, b_R_i_2
            real(kind = pr), dimension(0:p)         :: R_i_n, R_i_1, R_i_2, B_t_i
            real(kind = pr), dimension(0:p, 1:N_x)  :: h_tilde_1, u_tilde_1, b_tilde_1, h_tilde_2, u_tilde_2, b_tilde_2,    &
        & h_tilde_3, u_tilde_3, b_tilde_3, h_semi, u_semi, b_semi, h_semi_1, u_semi_1, b_semi_1, h_semi_2, u_semi_2, b_semi_2
            real(kind = pr), dimension(0:p, 1:N_x)  :: B_t, B_n, B_n_1, B_n_2

! 1e etape
            B_n = B_i

            call make_h_u_b_tilde(1, dt, wi_legendre, Phi_i, h_tilde, u_tilde, b_tilde,                                     &
        &                                                        h_tilde, u_tilde, b_tilde, h_tilde_1, u_tilde_1, b_tilde_1)
            
            call make_h_u_b_semi(wi_legendre, Phi_i, h_ini, u_ini, b_ini, h_tilde, u_tilde, b_tilde, h_semi, u_semi, b_semi)

            do i = 1, N_x
                call reconstruction_locale_ODE(i, h_semi, u_semi, b_semi, B_n, h_L_i_n, u_L_i_n, b_L_i_n,                   &
            &                                                                                     h_R_i_n, u_R_i_n, b_R_i_n)

                call make_R_i_ODE(i, dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, h_semi, u_semi, b_semi, B_n, R_i_n)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_i, b_L_i_n, b_R_i_n, R_i_n, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_1 = B_n + dt*B_t

            select case (cas)

            case (7, 11)
                ! call krivodonova_2007(t, B_n_1)

            end select

! 2e etape
            call make_h_u_b_tilde(2, dt, wi_legendre, Phi_i, h_tilde, u_tilde, b_tilde,                                     &
        &                                                  h_tilde_1, u_tilde_1, b_tilde_1, h_tilde_2, u_tilde_2, b_tilde_2)

            call make_h_u_b_semi(wi_legendre, Phi_i, h_ini, u_ini, b_ini, h_tilde_1, u_tilde_1, b_tilde_1,                  &
        &                                                                                      h_semi_1, u_semi_1, b_semi_1)

            do i = 1, N_x
                call reconstruction_locale_ODE(i, h_semi_1, u_semi_1, b_semi_1, B_n_1, h_L_i_1, u_L_i_1, b_L_i_1,           &
            &                                                                                     h_R_i_1, u_R_i_1, b_R_i_1)

                call make_R_i_ODE(i, dx, t + dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, h_semi_1, u_semi_1, b_semi_1,       &
            &                                                                                                  B_n_1, R_i_1)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_1, b_L_i_1, b_R_i_1, R_i_1, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_n_2 = 3._pr/4*B_n + 1._pr/4*(B_n_1 + dt*B_t)

            select case (cas)

            case (7, 11)
                ! call krivodonova_2007(t + dt, B_n_2)

            end select

! 3e etape
            call make_h_u_b_tilde(3, dt, wi_legendre, Phi_i, h_tilde, u_tilde, b_tilde,                                     &
        &                                                  h_tilde_2, u_tilde_2, b_tilde_2, h_tilde_3, u_tilde_3, b_tilde_3)

            call make_h_u_b_semi(wi_legendre, Phi_i, h_ini, u_ini, b_ini, h_tilde_2, u_tilde_2, b_tilde_2,                  &
        &                                                                                      h_semi_2, u_semi_2, b_semi_2)

            do i = 1, N_x
                call reconstruction_locale_ODE(i, h_semi_2, u_semi_2, b_semi_2, B_n_2, h_L_i_2, u_L_i_2, b_L_i_2,           &
            &                                                                                     h_R_i_2, u_R_i_2, b_R_i_2)

                call make_R_i_ODE(i, dx, t + dt/2, wi_legendre, x_gauss, Phi_i, Phi_xi_i, h_semi_2, u_semi_2, b_semi_2,     &
            &                                                                                                  B_n_2, R_i_2)

                call make_B_t_i(i, dx, wi_legendre, Phi_i, Phi_xi_i, B_n_2, b_L_i_2, b_R_i_2, R_i_2, B_t_i)

                B_t(:, i) = B_t_i
            end do
            B_i = 1._pr/3*B_n + 2._pr/3*(B_n_2 + dt*B_t)

            h_tilde = h_tilde_3 ; u_tilde = u_tilde_3 ; b_tilde = b_tilde_3
            
            select case (cas)

            case (7, 11)
                ! call krivodonova_2007(t + dt/2, B_i)

            end select

        end subroutine SSPRK_3_3_ODE

        elemental function phi_ABC(x) result(phi_ABC_x)
! ============================================================================ !
! Fonction phi_ABC pour les absorbing boundary conditions (ABCs)

! Entree de la fonction
!   - x : Point de Gauss du maillage
    
! Sortie de la fonction
!   - phi_ABC_x : Fonction phi_ABC evaluee en x
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)                 :: x

! Sortie de la fonction
            real(kind = pr)                             :: phi_ABC_x

            phi_ABC_x = 0._pr

            if (L_dom <= x .AND. x <= L) then
                phi_ABC_x = ((x - L_dom)/(L - L_dom))
            else if (L < x) then
                phi_ABC_x = 1._pr
            end if

        end function phi_ABC

        elemental function Gamma_ABC(x) result(Gamma_ABC_x)
! ============================================================================ !
! Fonction Gamma_ABC pour les ABC definie par Engsig-Karup et al.

! Entree de la fonction
!   - x : Point de Gauss du maillage
    
! Sortie de la fonction
!   - Gamma_ABC_x : Fonction Gamma_ABC evaluee en x
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: x

! Sortie de la fonction
            real(kind = pr)                         :: Gamma_ABC_x

            Gamma_ABC_x = -2._pr*(1._pr - phi_ABC(x))**3 + 3._pr*(1._pr - phi_ABC(x))**2

        end function Gamma_ABC

        subroutine apply_ABC(wi_legendre, x_gauss, Phi_i, h_ini, b_ini, h_tilde, u_tilde, B_i)
! ============================================================================ !
! Subroutine appliquant la relaxation sur les vecteurs h_tilde, u_tilde et
! B_i

! Entrees de la subroutine
!   - wi_legendre(1:N_G) : Tableau contenant les N_G poids d'interpolation
!   - x_gauss(1:N_G, 1:N_x) : Matrice contenant les coordonnees des
! points de Gauss
!   - Phi_i(0:p, 1:N_G) : Matrice contenant la valeur des polynomes de
! Legendre en chaque point de Gauss
!   - h_ini(0:p, 1:N_x) : Vecteur (de modes) contenant la hauteur d'eau 
! a l'instant initial
!   - b_ini(0:p, 1:N_x) : Vecteur contenant (de modes) la bathymetrie 
! a l'instant initial

! Entrees/Sorties de la subroutine :
!   - h_tilde(0:p, 1:N_x) : Vecteur contenant la hauteur d'eau semi-analytique
! en tous les points du domaine
!   - u_tilde(0:p, 1:N_x) : Vecteur contenant la vitesse semi-analytique
! en tous les points du domaine
!   - B_i(0:p, 1:N_x) : Vecteur contenant la bathymetrie en tous les points
! du domaine
! ============================================================================ !

! Entrees de la subroutine
            real(kind = pr), dimension(1:N_G), intent(in)               :: wi_legendre
            real(kind = pr), dimension(1:N_G, 1:N_x), intent(in)        :: x_gauss
            real(kind = pr), dimension(0:p, 1:N_G), intent(in)          :: Phi_i
            real(kind = pr), dimension(0:p, 1:N_x), intent(in)          :: h_ini, b_ini

! Entrees/Sorties de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: h_tilde, u_tilde, B_i

! Variables locales
            integer                                 :: i, j, k
            real(kind = pr)                         :: Gamma_xj
            real(kind = pr), dimension(1:N_G)       :: h_h_gauss, u_h_gauss, b_h_gauss, h_ini_gauss, b_ini_gauss
            real(kind = pr), dimension(1:N_G)       :: h_eq_gauss, u_eq_gauss, b_eq_gauss

            do i = 1, N_x
                if (L_dom < x_gauss(1, i)) then
                    h_h_gauss = MATMUL(TRANSPOSE(Phi_i), h_tilde(:, i))
                    u_h_gauss = MATMUL(TRANSPOSE(Phi_i), u_tilde(:, i))
                    b_h_gauss = MATMUL(TRANSPOSE(Phi_i), B_i(:, i))

                    h_ini_gauss = MATMUL(TRANSPOSE(Phi_i), h_ini(:, i))
                    b_ini_gauss = MATMUL(TRANSPOSE(Phi_i), b_ini(:, i))

                    do j = 1, N_G
                        Gamma_xj = Gamma_ABC(x_gauss(j, i))

                        h_eq_gauss(j) = h_ini_gauss(j)
                        u_eq_gauss(j) = 0._pr
                        b_eq_gauss(j) = 0._pr

                        h_h_gauss(j) = Gamma_xj*h_h_gauss(j) + (1._pr - Gamma_xj)*h_eq_gauss(j)
                        ! u_h_gauss(j) = Gamma_xj*u_h_gauss(j) + (1._pr - Gamma_xj)*u_eq_gauss(j)
                        b_h_gauss(j) = Gamma_xj*b_h_gauss(j) + (1._pr - Gamma_xj)*b_eq_gauss(j)
                    end do

                    do k = 0, p
                        h_tilde(k, i) = real(2*k+1, pr)/2*SUM(h_h_gauss*Phi_i(k, :)*wi_legendre)
                        ! u_tilde(k, i) = real(2*k+1, pr)/2*SUM(u_h_gauss*Phi_i(k, :)*wi_legendre)
                        B_i(k, i) = real(2*k+1, pr)/2*SUM(b_h_gauss*Phi_i(k, :)*wi_legendre)
                    end do
                end if
            end do

        end subroutine apply_ABC

end module mod_resolution_ODE