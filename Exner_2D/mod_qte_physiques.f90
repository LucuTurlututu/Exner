module mod_qte_physiques

    use mod_parametres

    implicit none

! ==================================================================================================== !
! Parametres physiques
! ==================================================================================================== !
! Coefficient de trainee laminaire
    real(kind = pr), parameter                      :: C_lam = 1.e-3_pr

! Coefficient de trainee turbulente
    real(kind = pr), parameter                      :: C_tur = 1.e-2_pr

! Masse volumique de l'eau
    real(kind = pr), parameter                      :: rho_f = 1.e3_pr

! Viscosite cinematique de l'eau
    real(kind = pr), parameter                      :: nu_f = 1.308e-6_pr

! Masse volumique du sol
    real(kind = pr), parameter                      :: rho_s = 2.650e3_pr

! Diametre median des grains
    real(kind = pr), parameter                      :: d_50 = 5.e-4_pr

! Gravite specifique des sediments
    real(kind = pr), parameter                      :: s_rel = rho_s/rho_f

! Constante de gravitation
    real(kind = pr), parameter                      :: g_const = 9.81_pr

! Porosite du sol
    real(kind = pr), parameter                      :: p_sol = 0.3_pr

! Coefficient de debit de Manning-Srickler
    real(kind = pr), parameter                      :: K_s = 35._pr

! ============================================================================ !
! Parametres des lois
! ============================================================================ !
! Exposant de la loi de Grass
    integer, parameter                              :: m_Grass = 3
    
! Constante de la loi de Grass
    real(kind = pr), parameter                      :: A_Grass = 1.e-3_pr

! Diametre adimentionnel
    real(kind = pr), parameter                      :: D_etoile = d_50*(((s_rel-1._pr)*g_const)/nu_f**2)**(1._pr/3)

! Vitesse de sedimentation
    real(kind = pr), parameter                      :: w_0 = SQRT(((13.95_pr*nu_f)/d_50)**2 +                               &
&                                                                1.09_pr*(s_rel - 1._pr)*g_const*d_50) - (13.95_pr*nu_f)/d_50

! Cas test de Berthon et al. (2012)
    real(kind = pr), parameter                      :: alpha_Berthon = 0.005_pr
    real(kind = pr), parameter                      :: beta_Berthon = 0.005_pr
    real(kind = pr), parameter                      :: C_Berthon = 1._pr

    real(kind = pr), parameter                      :: f_Berthon = 0.37415_pr
    real(kind = pr), parameter                      :: f_Berthon_leq = 73.65_pr
    real(kind = pr), parameter                      :: f_Berthon_geq = 236.91_pr

    real(kind = pr), parameter                      :: A_Berthon_G = 0.005_pr
    real(kind = pr), parameter                      :: A_Berthon_MPM = 8._pr*(f_Berthon/(8._pr*(s_rel-1._pr)*g_const*d_50))     &
&                                                                                 **(3._pr/2)*SQRT((s_rel-1._pr)*g_const*d_50**3)
    real(kind = pr), parameter                      :: A_Berthon_VR_leq = 0.053_pr*d_etoile**(-0.3_pr)*theta_c**(2.1_pr)*       &
&                                 (f_Berthon_leq/(8._pr*g_const*(s_rel-1._pr)*d_50))**(2.1_pr)*SQRT((s_rel-1._pr)*g_const*d_50**3)
    real(kind = pr), parameter                      :: A_Berthon_VR_geq = 0.1_pr*d_etoile**(-0.3_pr)*theta_c**(1.5_pr)*         &
&                                 (f_Berthon_geq/(8._pr*g_const*(s_rel-1._pr)*d_50))**(1.5_pr)*SQRT((s_rel-1._pr)*g_const*d_50**3)

    real(kind = pr), parameter                      :: u_cr_Berthon_G = 0._pr
    real(kind = pr), parameter                      :: u_cr_Berthon_MPM = theta_c*(f_Berthon/(8._pr*(s_rel-1._pr)*g_const*d_50))&
&                                                                                                                     **(-1._pr)
    real(kind = pr), parameter                      :: u_cr_Berthon_VR = theta_c*(f_Berthon/(8._pr*(s_rel-1._pr)*g_const*d_50)) &
&                                                                                                                     **(-1._pr)
    real(kind = pr), parameter                      :: Q_Berthon = 1._pr

    contains

        pure function norme_2(vect) result(norme_2_vect)
! ==================================================================================================== !
! Fonction retournant la norme 2 d'un vecteur

! Entree de la fonction :
!   - vect(:) : Vecteur quelqconque

! Sortie de la fonction :
!   - norme_2_vect : Norme euclidienne du vecteur vect
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), dimension(:), intent(in)                   :: vect

! Sortie de la fonction
            real(kind = pr)                                             :: norme_2_vect

            norme_2_vect = SQRT(SUM(vect**2))

        end function norme_2

        pure function b(Coord_phy, t) result(b_x_y_t)
! ==================================================================================================== !
! Fonction modelisant la hauteur de la bathymetrie b

! Entrees de la fonction :
!   - Coord_phy(2) : Vecteur (x, y) des coordonnees physiques
!   - t : Temps courant

! Sortie de la fonction :
!   - b_x_y_t : Hauteur de la bathymetrie evaluee en (x, y, t)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: t
            real(kind = pr), dimension(2), intent(in)                   :: Coord_phy

! Sortie de la fonction
            real(kind = pr)                                             :: b_x_y_t

! Variables locales
            real(kind = pr)                                             :: u_e_Berthon, u_Berthon, u_ini, u_1, b_0_Berthon

            select case (cas)

            case (1, 2, 3)
                b_x_y_t = SIN(2._pr*pi*(Coord_phy(1) - t))

            case (4, 5, 6)
                b_x_y_t = SIN(2._pr*pi*(Coord_phy(1) - t))*SIN(2._pr*pi*(Coord_phy(2) - t))

            case (7, 8)
! Cas test de Berthon et al. (2012)
                u_1 = DOT_PRODUCT(u(Coord_phy, 0._pr), [1._pr, 0._pr])

                b_0_Berthon = -(u_1**3 + 2._pr*g_const*Q_Berthon)/(2._pr*g_const*u_1) + C_Berthon

                b_x_y_t = -alpha_Berthon*t + b_0_Berthon

            case (9)
! Cas test de Berthon et al. (2012) avec VR
                u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_VR_leq)**(1._pr/2.1_pr)
                u_Berthon = SQRT(u_e_Berthon + u_cr_Berthon_VR)

                if (T_s(Q_Berthon/u_Berthon, [u_Berthon, 0._pr]) >= 3._pr) then
                    u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_VR_geq)**(2._pr/3._pr)
                    u_ini = SQRT(u_e_Berthon + u_cr_Berthon_VR)
                else
                    u_ini = u_Berthon
                end if

                b_0_Berthon = -(u_ini**3 + 2._pr*g_const*Q_Berthon)/(2._pr*g_const*u_ini) + C_Berthon

                b_x_y_t = -alpha_Berthon*t + b_0_Berthon

            case default
                b_x_y_t = 0._pr

            end select

        end function b

        pure function h(Coord_phy, t) result(h_x_y_t)
! ==================================================================================================== !
! Fonction modelisant la hauteur d'eau h

! Entrees de la fonction :
!   - Coord_phy(2) : Vecteur (x, y) des coordonnees physiques
!   - t : Temps courant

! Sortie de la fonction :
!   - h_x_y_t : Hauteur d'eau evaluee en (x, y, t)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: t
            real(kind = pr), dimension(2), intent(in)                   :: Coord_phy

! Sortie de la fonction
            real(kind = pr)                                             :: h_x_y_t

            select case (cas)
            
            case (1, 2, 3, 4, 5, 6)
                h_x_y_t = 1._pr/norme_2(u(Coord_phy, t))

            case (7, 8, 9)
                h_x_y_t = Q_Berthon/norme_2(u(Coord_phy, t))

            case default
                h_x_y_t = 1._pr

            end select

        end function h

        pure function u(Coord_phy, t) result(u_x_y_t)
! ==================================================================================================== !
! Fonction modelisant la vitesse de l'eau h

! Entrees de la fonction :
!   - Coord_phy(2) : Vecteur (x, y) des coordonnees physiques
!   - t : Temps courant

! Sortie de la fonction :
!   - u_x_y_t(1:2) : Vitesse de l'eau evaluee en (x, y, t)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: t
            real(kind = pr), dimension(2), intent(in)                   :: Coord_phy

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: u_x_y_t

! Variables locales
            real(kind = pr)                                             :: u_e_Berthon, u_Berthon

            select case (cas)

            case (1, 2, 3)
                u_x_y_t(1) = 1._pr + 1._pr/2._pr*b(Coord_phy, t)
                u_x_y_t(2) = 0._pr

            case (4, 5, 6)
                u_x_y_t(1) = 1._pr/2._pr*(1._pr + 1._pr/2._pr*b(Coord_phy, t))
                u_x_y_t(2) = 1._pr/2._pr*(1._pr + 1._pr/2._pr*b(Coord_phy, t))

            case (7)
! Cas test de Berthon et al. (2012) avec la loi de Grass
                u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_G)**(2._pr/3)

                u_x_y_t(1) = SQRT(u_e_Berthon + u_cr_Berthon_G)
                u_x_y_t(2) = 0._pr

            case (8)
! Cas test de Berthon et al. (2012) avec la loi de MPM
                u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_MPM)**(2._pr/3)

                u_x_y_t(1) = SQRT(u_e_Berthon + u_cr_Berthon_MPM)
                u_x_y_t(2) = 0._pr

            case (9)
! Cas test de Berthon et al. (2012) avec la loi de VR
                u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_VR_leq)**(1._pr/2.1_pr)

                u_Berthon = SQRT(u_e_Berthon + u_cr_Berthon_VR)

                if (T_s(Q_Berthon/u_Berthon, [u_Berthon, 0._pr]) >= 3._pr) then
                    u_e_Berthon = ((alpha_Berthon*Coord_phy(1) + beta_Berthon)/A_Berthon_VR_geq)**(2._pr/3._pr)

                    u_x_y_t(1) = SQRT(u_e_Berthon + u_cr_Berthon_VR)
                    u_x_y_t(2) = 0._pr
                else
                    u_x_y_t(1) = u_Berthon
                    u_x_y_t(2) = 0._pr
                end if

            case default
                u_x_y_t = 0._pr

            end select

        end function u

        pure function tau(u) result(tau_u)
! ==================================================================================================== !
! Fonction modelisant la contrainte de frottement a la paroi tau

! Entree de la fonction :
!   - u(1:2) : Vecteur vitesse d'eau

! Sortie de la fonction :
!   - tau_u(1:2) : Contrainte de frottement au fond
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: tau_u

            tau_u = rho_f*(C_lam + C_tur*norme_2(u))*u

        end function tau

        pure function C_Chezy(h) result(C_Chezy_h)
! ==================================================================================================== !
! Fonction renvoyant le coefficient de Chezy

! Entree de la fonction :
!   - h : Hauteur d'eau

! Sortie de la fonction :
!   - C_Chezy_h : Coefficient de Chezy
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), intent(in)                                     :: h

! Sortie de la fonction
            real(kind = pr)                                                 :: C_Chezy_h

            C_Chezy_h = K_s*h**(1._pr/6)

        end function C_Chezy

        pure function C_prime_Chezy(h) result(C_prime_Chezy_h)
! ==================================================================================================== !
! Fonction renvoyant le coefficient de peau de Chezy

! Entree de la fonction :
!   - h : Hauteur d'eau

! Sortie de la fonction :
!   - C_Chezy_h : Coefficient de peau de Chezy
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), intent(in)                                 :: h

! Sortie de la fonction
            real(kind = pr)                                             :: C_prime_Chezy_h

            C_prime_Chezy_h = SQRT(g_const)*(6.25_pr + 5.75_pr*LOG10(MAX(1._pr, h/(2.5_pr*d_50))))
        
        end function C_prime_Chezy

        pure function mu_corr(h) result(mu_corr_h)
! ==================================================================================================== !
! Fonction renvoyant le facteur de correction

! Entree de la fonction :
!   - h : Hauteur d'eau

! Sortie de la fonction :
!   - mu_corr_h : Facteur de correction
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), intent(in)                                 :: h

! Sortie de la fonction
            real(kind = pr)                                             :: mu_corr_h

            mu_corr_h = (C_Chezy(h)/C_prime_Chezy(h))**(3._pr/2)

        end function mu_corr

        pure function theta(u) result(theta_u)
! ==================================================================================================== !
! Fonction modelisant la contrainte de frottement a la paroi adimensionnee theta

! Entree de la fonction :
!   - u(1:2) : Vecteur vitesse d'eau

! Sortie de la fonction :
!   - theta_u(1:2) : Contrainte de frottement adimensionnee au fond
! ==================================================================================================== !
            
! Entree de la fonction
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: theta_u

            theta_u = tau(u)/(rho_f*g_const*(s_rel - 1._pr)*d_50)

        end function theta

        pure function theta_prime(h, u) result(theta_prime_h_u)
! ==================================================================================================== !
! Fonction renvoyant la contrainte de peau adimensionnee theta_prime

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - theta_prime_h_u : Contrainte de peau adimensionnee
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sorties de la fonction
            real(kind = pr)                                             :: theta_prime_h_u

            theta_prime_h_u = mu_corr(h)*norme_2(theta(u))

        end function theta_prime

        pure function T_s(h, u) result(T_s_h_u)
! ==================================================================================================== !
! Fonction renvoyant l'ecart relatif de Van Rijn

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - T_s_h_u : Ecart relatif de Van Rijn
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: T_s_h_u

            T_s_h_u = (theta_prime(h, u) - theta_c)/theta_c

        end function T_s

        pure function Grass(u) result(flux_Grass)
! ==================================================================================================== !
! Fonction renvoyant le flux sedimentaire de Grass

! Entree de la fonction :
!   u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   flux_Grass(1:2) : Flux sedimentaire suivant la loi de Grass
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: flux_Grass

            flux_Grass = A_Grass*(norme_2(u))**(m_Grass-1)*u

        end function Grass

        pure function MPM(h, u) result(flux_MPM)
! ==================================================================================================== !
! Fonction renvoyant le flux de charriage de MPM

! Entrees de la fonction :
!   h : Hauteur d'eau
!   u(1:2) : Vitesse d'eau

! Sorties de la fonction :
!   flux_MPM(1:2) : Flux de charriage suivant la loi de MPM
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: flux_MPM

            flux_MPM = SQRT(g_const*(s_rel-1._pr)*d_50**3)*8._pr*MAX(theta_prime(h, u) - theta_c, 0._pr)**(1.5_pr)          &
        &                                                                                                       *u/norme_2(u)

        end function MPM

        pure function Van_Rijn(h, u) result(flux_Van_Rijn)
! ==================================================================================================== !
! Fonction renvoyant le flux de charriage de Van Rijn

! Entrees de la fonction :
!   h : Hauteur d'eau
!   u(1:2) : Vitesse d'eau

! Sorties de la fonction :
!   flux_Van_Rijn(1:2) : Flux de charriage suivant la loi de Van Rijn
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: flux_Van_Rijn

            if (T_s(h, u) <= 0._pr) then
                flux_Van_Rijn = 0._pr
            else if (T_s(h, u) < 3._pr) then
                flux_Van_Rijn = 0.053_pr*D_etoile**(-0.3_pr)*T_s(h, u)**(2.1_pr)
            else if (T_s(h, u) >= 3._pr) then
                flux_Van_Rijn = 0.1_pr*D_etoile**(-0.3_pr)*T_s(h, u)**(1.5_pr)
            end if

            flux_Van_Rijn = SQRT(g_const*(s_rel-1._pr)*d_50**3)*flux_Van_Rijn*u/norme_2(u)

        end function Van_Rijn

end module mod_qte_physiques