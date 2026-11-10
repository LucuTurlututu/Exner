module mod_qte_physiques

    use mod_parametres

    implicit none

! ============================================================================ !
! Parametres physiques
! ============================================================================ !
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

! Diametre tel que 90% des grains ont un diametre inferieur
    real(kind = pr), parameter                      :: d_90 = 2.5_pr*d_50

! Gravite specifique des sediments
    real(kind = pr), parameter                      :: s_rel = rho_s/rho_f

! Constante de gravitation
    real(kind = pr), parameter                      :: g_part = 9.81_pr

! Concentration en sediments du sol
    real(kind = pr), parameter                      :: c_sol = 1._pr

! Porosite du sol
    real(kind = pr), parameter                      :: p_d = 0.3_pr

! Coefficient d'erosion
    real(kind = pr), parameter                      :: k_er = 1.e-3_pr

! Coefficient de debit de Manning-Srickler
    real(kind = pr), parameter                      :: K_s = 35._pr
! ============================================================================ !


! ============================================================================ !
! Parametres des lois
! ============================================================================ !
! Parametres de la loi de Grass
! Exposant de la loi de Grass
    integer, parameter                              :: m_Grass = 3
! Constante de la loi de Grass
    real(kind = pr), parameter                      :: A_Grass = 1.e-3_pr

! Seuil de Shields
    real(kind = pr), parameter                      :: theta_c = 0.047 ! Reel
    ! real(kind = pr), parameter                      :: theta_c = tau_c/(rho_f*(s_rel-1._pr)*g_part*d_50) ! Defini a partir de tau_c

! Diametre adimentionnel
    real(kind = pr), parameter                      :: d_etoile = d_50*(((s_rel-1._pr)*g_part)/nu_f**2)**(1._pr/3)

! Vitesse de sedimentation
    real(kind = pr), parameter                      :: w_0 = SQRT(((13.95_pr*nu_f)/d_50)**2 +                               &
            &                                                    1.09_pr*(s_rel - 1._pr)*g_part*d_50) - (13.95_pr*nu_f)/d_50
! ============================================================================ !


! ============================================================================ !
! Parametres des cas tests
! ============================================================================ !
! Debit utilise pour construire la hauteur d'eau
    real(kind = pr), parameter                      :: Q_SW = 1._pr

! Cas test 3
    real(kind = pr), parameter                      :: c_grass = 1._pr
    real(kind = pr), dimension(5), parameter        :: A_k = [0.10_pr, 0.05_pr, 0.12_pr, 0.08_pr, 0.15_pr]
    real(kind = pr), dimension(5), parameter        :: x_c_k = [0.20_pr, 0.35_pr, 0.50_pr, 0.65_pr, 0.80_pr]
    real(kind = pr), dimension(5), parameter        :: sigma_k = [0.03_pr, 0.01_pr, 0.05_pr, 0.02_pr, 0.04_pr]

! Cas test de Berthon et al. (2012) (cas 8 et 9)
    real(kind = pr), parameter                      :: alpha_Berthon = 0.005_pr
    real(kind = pr), parameter                      :: beta_Berthon = 0.005_pr
    real(kind = pr), parameter                      :: C_Berthon = 1._pr

    real(kind = pr), parameter                      :: f_Berthon = 0.37415_pr
    real(kind = pr), parameter                      :: f_Berthon_leq = 73.65_pr
    real(kind = pr), parameter                      :: f_Berthon_geq = 236.91_pr

    real(kind = pr), parameter                      :: A_Berthon_G = 0.005_pr
    real(kind = pr), parameter                      :: A_Berthon_MPM = 8._pr*(f_Berthon/(8._pr*(s_rel-1._pr)*g_part*d_50))      &
&                                                                                 **(3._pr/2)*SQRT((s_rel-1._pr)*g_part*d_50**3)
    real(kind = pr), parameter                      :: A_Berthon_VR_leq = 0.053_pr*d_etoile**(-0.3_pr)*theta_c**(2.1_pr)*       &
&                                 (f_Berthon_leq/(8._pr*g_part*(s_rel-1._pr)*d_50))**(2.1_pr)*SQRT((s_rel-1._pr)*g_part*d_50**3)
    real(kind = pr), parameter                      :: A_Berthon_VR_geq = 0.1_pr*d_etoile**(-0.3_pr)*theta_c**(1.5_pr)*         &
&                                 (f_Berthon_geq/(8._pr*g_part*(s_rel-1._pr)*d_50))**(1.5_pr)*SQRT((s_rel-1._pr)*g_part*d_50**3)

    real(kind = pr), parameter                      :: u_cr_Berthon_G = 0._pr
    real(kind = pr), parameter                      :: u_cr_Berthon_MPM = theta_c*(f_Berthon/(8._pr*(s_rel-1._pr)*g_part*d_50)) &
&                                                                                                                     **(-1._pr)
    real(kind = pr), parameter                      :: u_cr_Berthon_VR = theta_c*(f_Berthon/(8._pr*(s_rel-1._pr)*g_part*d_50)) &
&                                                                                                                     **(-1._pr)
    real(kind = pr), parameter                      :: Q_Berthon = 1._pr

! Cas test de Fernandez-Nieto et al. (2025) (cas 10)
    real(kind = pr), parameter                      :: a_FN = 0.01_pr
    real(kind = pr), parameter                      :: b_FN = 0.01_pr
    real(kind = pr), parameter                      :: c_FN = 1._pr/(1._pr - p_d)
    real(kind = pr), parameter                      :: h_0_FN = 10._pr
! ============================================================================ !

    contains

        elemental function b(x, t) result(b_xt)
! ============================================================================ !
! Fonction b modelisant la bathymetrie                                          

! Entrees de la fonction :
!   - x : point de discretisation
!   - t : temps de discretisation

! Sortie de la fonction :
!   - b_xt : fonction b evaluee en (x, t)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: x, t

! Sortie de la fonction
            real(kind = pr)                         :: b_xt

! Variables locales
            integer                                 :: k
            real(kind = pr)                         :: x_c_t, delta_x, h_0, Beta, b_0_Berthon, u_e_Berthon, u_Berthon, u_ini

            b_xt = 0._pr

            select case (cas)

            case (1, 4, 5, 6)
! Solution manufacturee
                b_xt = SIN(2*pi*(x-t))

            case (2)
! Solution manufacturee stationnaire
                b_xt = b_0 + 1._pr/2*SIN(2*pi*x)

            case (3)
! Cas test bathymetrie bossue
                do k = 1, 5
                    x_c_t = MODULO(x_c_k(k) + c_grass*t, 1._pr)
        
                    delta_x = x - x_c_t
                    if (delta_x >  0.5_pr) delta_x = delta_x - 1._pr
                    if (delta_x < -0.5_pr) delta_x = delta_x + 1._pr
                    
                    b_xt = b_xt + A_k(k) * EXP(-(delta_x**2)/(2._pr*sigma_k(k)**2))
                end do

            case (7)
! Cas test de Long et al. (2008)
                h_0 = 6._pr ; Beta = 0.01_pr
                x_c_t = 150._pr

                b_xt = -h_0 + 2._pr*EXP(-Beta*(x - x_c_t)**2)

            case (8, 9)
! Cas test de Berthon et al. (2012)
                b_0_Berthon = -(u(x, 0._pr)**3 + 2._pr*g_part*Q_Berthon)/(2._pr*g_part*u(x, 0._pr)) + C_Berthon

                b_xt = -alpha_Berthon*t + b_0_Berthon
            
            case (10)
! Cas test de Berthon et al. (2012) avec VR
                u_e_Berthon = ((alpha_Berthon*x + beta_Berthon)/A_Berthon_VR_leq)**(1._pr / 2.1_pr)
                u_Berthon = SQRT(u_e_Berthon + u_cr_Berthon_VR)

                if ((theta_prime(Q_Berthon/u_Berthon, u_Berthon) - theta_c) / theta_c >= 3._pr) then
                    u_e_Berthon = ((alpha_Berthon*x + beta_Berthon) / A_Berthon_VR_geq)**(2._pr / 3._pr)
                    u_ini = SQRT(u_e_Berthon + u_cr_Berthon_VR)
                else
                    u_ini = u_Berthon
                end if

                b_0_Berthon = -(u_ini**3 + 2._pr*g_part*Q_Berthon)/(2._pr*g_part*u_ini) + C_Berthon

                b_xt = -alpha_Berthon*t + b_0_Berthon

            case (11)
! Cas test de Fernandez-Nieto et al. (2025)
                b_xt = 3._pr*(1._pr/(1._pr - p_d))*A_Grass*a_FN/c_FN*u(x, 0._pr)**2
            
            end select

        end function b

        elemental function h(x, t) result(h_xt)
! ============================================================================ !
! Fonction h renvoyant la hauteur de SW

! Entrees de la fonction :
!   - x : Point du maillage
!   - t : Temps actuel
    
! Sortie de la fonction
!   - h_xt : Hauteur d'eau de SW
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: x, t

! Sortie de la fonction
            real(kind = pr)                         :: h_xt

            h_xt = 1._pr

            select case (cas)

            case (5, 6)
                h_xt = Q_SW/u(x, t)

            case (8, 9, 10)
! Cas de Berthon et al. (2012)
                h_xt = Q_Berthon/u(x, t)

            case (11)
! Cas de Fernandez-Nieto et al. (2025)
                h_xt = h_0_FN

            end select

        end function h

        elemental function u(x, t) result(u_xt)
! ============================================================================ !
! Fonction u modelisant la vitesse de l'eau

! Entrees de la fonction :
!   - x : Point du maillage
!   - t : Temps actuel

! Sortie de la fonction :
!   - u_xt : Fonction u evaluee (x, t)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: x, t

! Sortie de la fonction
            real(kind = pr)                         :: u_xt

! Variables locales
            real(kind = pr)                         :: Q, u_e_Berthon, u_Berthon

            u_xt = 0._pr

            select case(cas)

            case (1, 4, 5, 6)
! Solution manufacturee
                u_xt = 1._pr + 1._pr/2*b(x, t)

            case (2, 3)
! Solution manufacturee stationnaire et cas test bathymetrie bossue
                u_xt = 1._pr + 1._pr/2*SIN(2*pi*x)

            case (7)
! Cas test de Long et al. (2008)
                Q = 10._pr

                u_xt = Q/(-b(x, t))

            case (8)
! Cas test de Berthon et al. (2012) avec la loi de Grass
                u_e_Berthon = ((alpha_Berthon*x + beta_Berthon)/A_Berthon_G)**(2._pr/3)

                u_xt = SQRT(u_e_Berthon + u_cr_Berthon_G)

            case (9)
! Cas test de Berthon et al. (2012) avec la loi de MPM
                u_e_Berthon = ((alpha_Berthon*x + beta_Berthon)/A_Berthon_MPM)**(2._pr/3)

                u_xt = SQRT(u_e_Berthon + u_cr_Berthon_MPM)

            case (10)
! Cas test de Berthon et al. (2012) avec la loi de VR
                u_e_Berthon = ((alpha_Berthon*x + beta_Berthon)/A_Berthon_VR_leq)**(1._pr/2.1_pr)

                u_Berthon = SQRT(u_e_Berthon + u_cr_Berthon_VR)

                if ((theta_prime(Q_Berthon/u_Berthon, u_Berthon) - theta_c)/theta_c >= 3._pr) then
                    u_e_Berthon = ((alpha_Berthon*x + beta_Berthon)/A_Berthon_VR_geq)**(2._pr/3._pr)

                    u_xt = SQRT(u_e_Berthon + u_cr_Berthon_VR)
                else
                    u_xt = u_Berthon
                end if

            case (11)
! Cas test de Fernandez-Nieto et al. (2025)
                u_xt = a_FN*x + b_FN

            end select

        end function u

        elemental function tau(u) result(tau_u)
! ============================================================================ !
! Fonction tau modelisant la contrainte au fond

! Entree de la fonction :
!   - u : Vitesse de l'eau
        
! Sortie de la fonction
!   - tau_u : Contrainte au fond
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: u

! Sortie de la fonction
            real(kind = pr)                         :: tau_u

! Variables locales
            real(kind = pr)                         :: K_d

            select case (cas)

            case (1, 2, 3)
                tau_u = u*ABS(u)

            case default
                K_d = C_lam + C_tur*ABS(u)

                tau_u = rho_f*K_d*ABS(u)

            end select

        end function tau

        elemental function C_Chezy(h) result(C_Chezy_h)
! ============================================================================ !
! Fonction C_Chezy renvoyant le coefficient de Chezy suivant la loi de
! Manning-Strickler

! Entree de la fonction :
!   - h : Hauteur d'eau de SW
    
! Sortie de la fonction :
!   - C_Chezy_h : Coefficient de Chezy
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: h

! Sortie de la fonction
            real(kind = pr)                         :: C_Chezy_h

            C_Chezy_h = K_s*h**(1._pr/6)

        end function C_Chezy

        elemental function C_prime_Chezy(h) result(C_prime_Chezy_h)
! ============================================================================ !
! Fonction C_prime_Chezy renvoyant le coefficient de peau de Chezy utilisant
! la rugosite de Nikuradse

! Entree de la fonction :
!   - h : Hauteur d'eau de SW
    
! Sortie de la fonction :
!   - C_prime_Chezy_h : Coefficient de peau de Chezy
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: h

! Sortie de la fonction
            real(kind = pr)                         :: C_prime_Chezy_h

            C_prime_Chezy_h = SQRT(g_part)*(6.25_pr + 5.75_pr*LOG10(MAX(1._pr, h/(2.5_pr*d_50))))

        end function C_prime_Chezy

        elemental function mu(h) result (mu_h)
! ============================================================================ !
! Fonction mu renvoyant le facteur de forme permettant de passer de la
! contrainte totale adimensionee a la contrainte de peau adimensionnee
!       theta_prime = mu*theta

! Entree de la fonction :
!   - h : Hauteur d'eau de SW
    
! Sortie de la fonction
!   - mu_h : Facteur de forme pour la hauteur h
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: h

! Sortie de la fonction
            real(kind = pr)                         :: mu_h

            mu_h = (C_Chezy(h)/C_prime_Chezy(h))**(3._pr/2)

        end function mu

        elemental function theta_prime(h, u) result(theta_prime_hu)
! ============================================================================ !
! Fonction theta_prime renvoyant la contrainte de peau adimensionnee

! Entrees de la fonction :
!   - h : Hauteur de l'eau
!   - u : Vitesse de l'eau
    
! Sortie de la fonction
!   - theta_prime_hu : Contrainte de peau adimensionnee pour le couple
! de variables de SW (h, u)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: h, u

! Sortie de la fonction
            real(kind = pr)                         :: theta_prime_hu, theta_Berthon

            select case (cas)

            case (8, 9)
! Cas de Berthon et al. (2012)
                theta_prime_hu = (f_Berthon*u**2)/(8._pr*(s_rel-1._pr)*g_part*d_50)

            case (10)
! Cas de Berthon et al. avec la loi de Van Rijn
                theta_Berthon = (f_Berthon_leq*u**2)/(8._pr*(s_rel-1._pr)*g_part*d_50)
                if ((theta_Berthon - theta_c)/theta_c <= 3._pr) then
                    theta_prime_hu = theta_Berthon
                else
                    theta_prime_hu = (f_Berthon_geq*u**2)/(8._pr*(s_rel-1._pr)*g_part*d_50)
                end if

            case default
                theta_prime_hu = mu(h)*(tau(u)/(rho_f*(s_rel-1._pr)*g_part*d_50))

            end select

        end function theta_prime

        elemental function Grass(u) result(flux_Grass)
! ============================================================================ !
! Fonction implementant le flux de Grass

! Entree de la fonction :
!   - u : Vitesse issue de SW
        
! Sortie de la fonction :
!   - flux_Grass : Flux de charriage suivant la loi de Grass
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)                 :: u

! Sortie de la fonction
            real(kind = pr)                             :: flux_Grass

            flux_Grass = A_Grass*u*ABS(u)**(m_Grass-1)

        end function Grass

        elemental function MPM(h, u) result(flux_MPM)
! ============================================================================ !
! Fonction implementant le flux de Meyer-Peter et Müller

! Entrees de la fonction :
!   - h : Hauteur issue de SW
!   - u : Vitesse issue de SW
        
! Sortie de la fonction :
!   - flux_MPM : Flux de charriage suivant la loi de MPM
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)                 :: h, u

! Sortie de la fonction
            real(kind = pr)                             :: flux_MPM

            flux_MPM = SQRT(g_part*(s_rel-1._pr)*d_50**3)*8._pr*MAX(0._pr, theta_prime(h, u) - theta_c)**(3._pr/2)

        end function MPM

        elemental function Van_Rijn(h, u) result(flux_Van_Rijn)
! ============================================================================ !
! Fonction implementant le flux de Van-Rijn

! Entrees de la fonction :
!   - h : Hauteur issue de SW
!   - u : Vitesse issue de SW
        
! Sortie de la fonction :
!   - flux_Van_Rijn : Flux de charriage suivant la loi de Van_Rijn
! ============================================================================ !
        
! Entrees de la fonction
            real(kind = pr), intent(in)             :: h, u

! Sortie de la fonction
            real(kind = pr)                         :: flux_Van_Rijn

! Variable locale
            real(kind = pr)                         :: T

            T = (theta_prime(h, u) - theta_c)/theta_c
            flux_Van_Rijn = 0._pr
            
            if (0._pr < T .AND. T < 3._pr) then
                flux_Van_Rijn = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.053_pr*d_etoile**(-0.3_pr)*T**2.1_pr
            else if (3._pr <= T) then
                flux_Van_Rijn = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.1_pr*d_etoile**(-0.3_pr)*T**(3._pr/2)
            end if

        end function Van_Rijn

end module mod_qte_physiques
