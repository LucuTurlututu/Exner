module mod_fonctions

    use mod_parametres
    use mod_qte_physiques

    implicit none
    
    contains

        elemental function k_s_total(h, u) result(k_s_total_hu)
! ============================================================================ !
! Fonction k_s_total pour le coefficient de rugosite totale equivalente
! de Nikuradse

! Entrees de la fonction :
!   - h : Hateur d'eau issue de SW
!   - u : Vitesse issue de SW

! Sortie de la fonction :
!   - k_s_total : Rugosite equivalente de Nikuradse
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: h, u

! Sortie de la fonction
            real(kind = pr)                         :: k_s_total_hu

! Variables locales
            real(kind = pr)                         :: Delta_r, lambda_r, Delta_d, lambda_d, T_s

! Grandeurs liees aux rides
            Delta_r = 250._pr*d_50
            lambda_r = 1000._pr*d_50

! Grandeurs liees aux dunes
            T_s = (theta_prime(h, u) - theta_c)/theta_c
            Delta_d = 0.11_pr*h*(d_50/h)**(0.3_pr)*(1._pr - EXP(-0.5_pr*T_s))*(25._pr - T_s)
            lambda_d = 7.3_pr*h

            k_s_total_hu = 2.5_pr*d_50 + 28._pr*(Delta_r**2/lambda_r) + 1.1_pr*Delta_d*(1-EXP(-25._pr*Delta_d/lambda_d))

        end function k_s_total

        elemental function c_Delta(h, u) result(c_Delta_hu)
! ============================================================================ !
! Fonction c_Delta representant la concentration volumique de reference a
! l'interface entre la couche de charriage et la zone de suspension

! Entrees de la fonction :
!   - h : Hateur d'eau issue de SW
!   - u : Vitesse issue de SW

! Sortie de la fonction :
!   - c_Delta_hu : Concentration volumique de reference
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: h, u

! Sortie de la fonction
            real(kind = pr)                         :: c_Delta_hu

! Variables locales
            real(kind = pr)                         :: T_s, Delta_ref

            T_s = (theta_prime(h, u) - theta_c)/theta_c
            Delta_ref = MAX(k_s_total(h, u), 0.01_pr*h)

            c_Delta_hu = 0.015_pr*(d_50/Delta_ref)*(T_s**(1.5_pr)/d_etoile**(0.3_pr))

        end function c_Delta

        elemental function Qe(b, h, u) result(Qe_bhu)
! ============================================================================ !
! Fonction Q_e modelisant le flux de d'erosion

! Entrees de la fonction :
!   - b : bathymetrie au point du maillage ou debuterait l'erosion
!   - h : Hauteur d'eau issue de SW
!   - u : vitesse au point du maillage ou debuterait l'erosion

! Sortie de la fonction :
!   - Qe_bhu : fonction Q_e evaluee en (b, h, u)
! ============================================================================ !
        
! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, h, u

! Sortie de la fonction
            real(kind = pr)                         :: Qe_bhu

! Variables locales
            real(kind = pr)                         :: a, b_exp

            Qe_bhu = 0._pr

            select case (cas)

            case (1, 2, 3)
                Qe_bhu = b*u

            case (4)
! Loi de Grass
                Qe_bhu = Grass(u)

            case (5)
! Loi de MPM
                Qe_bhu = MPM(h, u)

            case (6)
! Loi de Van Rijn
                Qe_bhu = Van_Rijn(h, u)

            case (7)
                a = 0.001_pr ; b_exp = 3._pr

                Qe_bhu = a*u**b_exp

            case (11)
                Qe_bhu = 1._pr/(1._pr - p_d)*Grass(u)

            end select

        end function Qe

        elemental function Qd(b, u) result(Qd_bu)
! ============================================================================ !
! Fonction Q_d modelisant le flux de depot

! Entrees de la fonction :
!   - b : bathymetrie au point du maillage ou debuterait le depot
!   - u : vitesse au point du maillage ou debuterait le depot

! Sortie de la fonction :
!   - Qd_bu : fonction Q_d evaluee en (b, u)
! ============================================================================ !
        
! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, u

! Sortie de la fonction
            real(kind = pr)                         :: Qd_bu

! Variable locale
            real(kind = pr)                         :: a, b_exp

            Qd_bu = 0._pr

            select case (cas)

            case (1, 2, 3)
                Qd_bu = b*u*EXP(-ABS(u))

            case (4)
! Loi de Grass
                Qd_bu = Grass(u)

            case (5, 6)
! Loi MPM
                Qd_bu = 0._pr

            case (7)
                a = 0.001_pr ; b_exp = 3._pr

                Qd_bu = a*u**b_exp

            case (11)
                Qd_bu = 1._pr/(1._pr - p_d)*Grass(u)

            end select

        end function Qd

        elemental function D(b, h, u) result (D_bhu)
! ============================================================================ !
! Fonction D modelisant le terme source de depot

! Entrees de la fonction :
!   - b : Bathymetrie au point du maillage ou debuterait le depot
!   - h : Hauteur d'eau issue de SW
!   - u : vitesse au point du maillage ou debuterait le depot

! Sortie de la fonction :
!   - D_bhu : fonction D evluee en (b, h, u)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, h, u

! Sortie de la fonction
            real(kind = pr)                         :: D_bhu

            D_bhu = 0._pr

            select case (cas)

            case (1, 2, 3)
                D_bhu = EXP(-ABS(u))

            case (4, 5, 6)
                D_bhu = w_0*c_Delta(h, u)*MAX(1._pr - theta_prime(h, u)/theta_c, 0._pr)

            end select

        end function D

        elemental function E(h, u) result (E_hu)
! ============================================================================ !
! Fonction E modelisant le terme source d'erosion

! Entrees de la fonction :
!   - h : Hauteur d'eau de SW
!   - u : vitesse au point du maillage ou debuterait l'erosion

! Sortie de la fonction :
!   - E_hu : fonction D evluee en (h, u)
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: h, u

! Sortie de la fonction
            real(kind = pr)                         :: E_hu

            E_hu = 0._pr

            select case (cas)

            case (1, 2, 3)
                E_hu = ABS(u)

            case (4, 5, 6)
                E_hu = w_0*c_Delta(h, u)

            end select

        end function E

        elemental function I_lisse(tau) result(indi_lisse_tau)
! ============================================================================ !
! Fonction lissant les indicatrices

! Entree de la fonction :
!   - tau : Contrainte au point x_i

! Sortie de la fonction :
!   - indi_lisse_tau : Fonction sigmoide lissant l'indicatrice
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: tau

! Sortie de la fonction
            real(kind = pr)                         :: indi_lisse_tau

            indi_lisse_tau = 0.5_pr*(1._pr + TANH((tau - tau_c)/eps_tau))

        end function I_lisse

        elemental function I_leq(tau) result(indi_leq)
! ============================================================================ !
! Fonction indicatrice pour la contrainte tau <= tau_c

! Entree de la fonction :
!   - tau : contrainte au point x_i

! Sortie de la fonction :
!   - indi_leq : indicatrice de tau <= tau_c
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: tau

! Sortie de la fonction
            real(kind = pr)                         :: indi_leq

            indi_leq = 1._pr - I_lisse(tau)

        end function I_leq

        elemental function I_geq(tau) result(indi_geq)
! ============================================================================ !
! Fonction indicatrice pour la contrainte tau > tau_c

! Entree de la fonction :
!   - tau : contrainte au point x_i

! Sortie de la fonction :
!   - indi_geq : indicatrice de tau > tau_c
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: tau

! Sortie de la fonction
            real(kind = pr)                         :: indi_geq

            indi_geq = I_lisse(tau)

        end function I_geq

        elemental function F(b, h, u, tau) result(F_bhu)
! ============================================================================ !
! Fonction F pour le terme de flux

! Entrees de la fonction :
!   - b : bathymetrie au point x_i
!   - h : Hauteur d'eau issue de SW
!   - u : vitesse au point x_i
!   - tau : contrainte au point x_i

! Sortie de la fonction :
!   - F_bhu : terme source F evalue en (b, u)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, h, u, tau

! Sortie de la fonction
            real(kind = pr)                         :: F_bhu

! Variables locales
            real(kind = pr)                         :: u_Long, u_e_Berthon

            select case (cas)

            case (7)
! Cas de Long et al. (2008)
                u_Long = 10._pr/(0._pr - b)

                F_bhu = (0.001_pr*u_Long**3)

            case (8)
! Cas de Berthon et al. (2012) avec la loi de Grass
                u_e_Berthon = MAX(0._pr, u**2 - u_cr_Berthon_G)

                F_bhu = A_Berthon_G*u_e_Berthon**(3._pr/2)

            case (9)
! Cas de Berthon et al. (2012) avec la loi de MPM
                u_e_Berthon = MAX(0._pr, u**2 - u_cr_Berthon_MPM)

                F_bhu = A_Berthon_MPM*u_e_Berthon**(3._pr/2)

            case (10)
! Cas de Berthon et al. (2012) avec la loi de Van Rijn
                u_e_Berthon = MAX(0._pr, u**2 - u_cr_Berthon_VR)

                if ((theta_prime(h, u) - theta_c) / theta_c < 3._pr) then
                    F_bhu = A_Berthon_VR_leq*u_e_Berthon**(2.1_pr)
                else
                    F_bhu = A_Berthon_VR_geq*u_e_Berthon**(1.5_pr)
                end if

            case default
                F_bhu = Qe(b, h, u)*I_geq(tau) + Qd(b, u)*I_leq(tau)

            end select

        end function F

        elemental function F_b(b, h, u, tau) result(F_b_bhu)
! ============================================================================ !
! Derivee de la fonction F par rapport a b

! Entrees de la fonction :
!   - b : bathymetrie au point x_i
!   - h : Hauteur d'eau issue de SW
!   - u : vitesse au point x_i
!   - tau : contrainte au point x_i

! Sortie de la fonction :
!   - F_b_bhu : Derivee de F par rapport a b
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, h, u, tau

! Sortie de la fonction
            real(kind = pr)                         :: F_b_bhu

! Variables locales
            real(kind = pr)                         :: u_b, h_b, K_d, tau_b, C, C_prime, C_C_prime_h, mu_b,  &
        &  theta_prime_b, I_b, Q_Long, a_Long, T, Fr, F_q

            F_b_bhu = u

            select case (cas)

            case (1, 2, 3)
                F_b_bhu = u*EXP(-ABS(u))*I_leq(tau) + u*I_geq(tau)

            case (4)
                ! u_b = 1._pr/2

                ! F_b_bhu = A_Grass*u_b*ABS(u)**(m_Grass-2)*(ABS(u) + (m_Grass-1)*SIGN(1._pr, u)*u)

                if (eps_var < h) then
                    Fr = ABS(u)/SQRT(g_part*h)
                    
                    F_q = A_Grass*1._pr/h*m_Grass*ABS(u)**(m_Grass - 1)

                    F_b_bhu = F_q*u/(1._pr - Fr**2)
                else
                    F_b_bhu = 0._pr
                end if

            case (5)
                K_d = C_lam + C_tur*ABS(u)
                u_b = 1._pr/2
                h_b = -h*u_b/u
                tau_b = rho_f*SIGN(1._pr, u)*u_b*(C_tur*ABS(u) + K_d)

                C = C_Chezy(h)
                C_prime = C_prime_Chezy(h)
                C_C_prime_h = 1._pr/C_prime**2*1._pr/6*K_s*h**(-5._pr/6)*(C_prime - 6._pr*18._pr/LOG(10._pr))
                mu_b = 3._pr/2*h_b*C_C_prime_h*(C/C_prime)**(1._pr/2)
                theta_prime_b = 1._pr/(rho_f*g_part*(s_rel-1._pr)*d_50)*(tau*mu_b + mu(h)*tau_b)
                
                I_b = 1._pr/(2*eps_tau)*(1._pr - TANH((tau - tau_c)/eps_tau)**2)*tau_b

                F_b_bhu = 8._pr*SQRT(g_part*(s_rel-1._pr)*d_50**3)*(3._pr/2)*theta_prime_b*              &
            &           (MAX(0._pr, theta_prime(h, u) - theta_c))**(1._pr/2)*I_geq(tau) + Qe(b, h, u)*I_b

            case (6)
                K_d = C_lam + C_tur*ABS(u)
                u_b = 1._pr/2
                h_b = -h*u_b/u
                tau_b = rho_f*SIGN(1._pr, u)*u_b*(C_tur*ABS(u) + K_d)

                C = C_Chezy(h)
                C_prime = C_prime_Chezy(h)
                C_C_prime_h = 1._pr/C_prime**2*1._pr/6*K_s*h**(-5._pr/6)*(C_prime - 6._pr*18._pr/LOG(10._pr))
                mu_b = 3._pr/2*h_b*C_C_prime_h*(C/C_prime)**(1._pr/2)
                theta_prime_b = 1._pr/(rho_f*g_part*(s_rel-1._pr)*d_50)*(tau*mu_b + mu(h)*tau_b)
                
                I_b = 1._pr/(2*eps_tau)*(1._pr - TANH((tau - tau_c)/eps_tau)**2)*tau_b

                T = (theta_prime(h, u) - theta_c)/theta_c
                if (0._pr < T .AND. T < 3._pr) then
                    F_b_bhu = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.053_pr*d_etoile**(-0.3_pr)*           &
                &                    2.1_pr/theta_c*theta_prime_b*T**1.1_pr*I_geq(tau) + Qe(b, h, u)*I_b
                else if (3._pr <= T) then
                    F_b_bhu = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.1_pr*d_etoile**(-0.3_pr)*             &
                &          3._pr/(2._pr*theta_c)*theta_prime_b*T**(1._pr/2)*I_geq(tau) + Qe(b, h, u)*I_b
                end if

            case (7)
                if (eps_var < h) then
                    Fr = ABS(u)/SQRT(g_part*h)
                    
                    F_q = A_Grass*1._pr/h*m_Grass*ABS(u)**(m_Grass - 1)

                    F_b_bhu = F_q*u/(1._pr - Fr**2)
                else
                    F_b_bhu = 0._pr
                end if

            case (8, 9, 10)
                F_b_bhu = 0._pr

            case (11)
! Valeur propre de De Vries
                if (eps_var < h) then
                    Fr = ABS(u)/SQRT(g_part*h)
                    
                    F_q = A_Grass*1._pr/h*m_Grass*ABS(u)**(m_Grass - 1)

                    F_b_bhu = F_q*u/(1._pr - Fr**2)
                else
                    F_b_bhu = 0._pr
                end if

            end select

        end function F_b

        elemental function m(b_x) result(m_b_x)
! ============================================================================ !
! Fonction renvoyant la fonction m liee a la metrique

! Entree de la fonction :
!   - b_x : Derivee de la bathymetrie par rapport a x

! Sortie de la fonction :
!   - m_b_x : Fonction m = 1/SQRT(1+b_x**2)
! ============================================================================ !

! Entree de la fonction
            real(kind = pr), intent(in)             :: b_x

! Sortie de la fonction
            real(kind = pr)                         :: m_b_x

            select case (cas)

            case (7)
                m_b_x = 1._pr - 0.4_pr

            case (8, 9, 10, 11)
                m_b_x = 1._pr
                
            case default
                m_b_x = 1._pr/SQRT(1._pr + b_x**2)

            end select

        end function m

        elemental function S(b, h, u, tau) result(S_bhu)
! ============================================================================ !
! Fonction S pour le terme source

! Entrees de la fonction :
!   - b : bathymetrie au point x_i
!   - h : Hauteur d'eau issue de SW
!   - u : vitesse au point x_i
!   - tau : contrainte au point x_i

! Sortie de la fonction :
!   - S_bhu : terme source S evalue en (b, h, u)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: b, h, u, tau

! Sortie de la fonction
            real(kind = pr)                         :: S_bhu

            select case (cas)

            case (7, 8, 9, 10, 11)
                S_bhu = 0._pr

            case default
                S_bhu = D(b, h, u)*I_leq(tau) - E(h, u)*I_geq(tau)

            end select

        end function S

        elemental function R_tilde(x, t) result(R_tilde_xt)
! ============================================================================ !
! Fonction R_tilde pour la solution manufacturee associee a la solution exacte
!       b_tilde(x, t) = SIN(2*pi*(x-t))
! a la vitesse et la contrainte
!       u(x, t) = 1 + 1/2*SIN(2*pi*(x-t))
!       tau(u) = u*ABS(u)
! au flux
!       F(b_tilde, u, tau) = Q_e(b_tilde, u) I_leq(tau) + Q_d(b_tilde, u) I_geq(tau)
! ou :
!       Q_e = b*u
!       Q_d = b*u*EXP(-ABS(u))
! et au terme source
!       S(u, tau) = D(u) I_leq(tau) - E(u) I_geq(tau)
! avec :
!       D(u) = EXP(-ABS(u))
!       E(u) = ABS(u)
! Ainsi R_tilde est definie comme
!       R_tilde = m(b_tilde_x) b_tilde_t + F_x(b_tilde) - S(b_tilde)
! ou :
!       m(b_tilde_x) = 1/SQRT(1+b_tilde_x**2)

! Entrees de la fonction :
!   - x : Point d'interpolation de gauss
!   - t : Temps actuel

! Sortie de la fonction :
!   - R_tilde_xt : Fonction R_tilde evaluee aux entrees de la fonction
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: x, t

! Sortie de la fonction
            real(kind = pr)                         :: R_tilde_xt

! Variables locales
            integer                                 :: k
            real(kind = pr)                         :: m_b_x, b_t, b_tilde, Q_e_x, Q_d_x, I_x,          &
        &   Qe_Qd_I_x, F_b_x, S_b, x_c_t, b_x, b_k_xt, delta_x, u_x, u_tilde, tau_x, tau_tilde,         &
        &   C, C_prime, C_C_prime_h, h_x, mu_x, theta_prime_x, T_Van_Rijn, h_tilde

            select case (cas)

            case (1)
! Solution manufacturee
! Terme temporel
                m_b_x = 1._pr/SQRT(1._pr + 4*pi**2*COS(2*pi*(x-t))**2)
                b_t = -2*pi*COS(2*pi*(x-t))

! Terme de flux
                b_tilde = SIN(2*pi*(x-t))
                h_tilde = h(x, t)
                Q_e_x = 2*pi*COS(2*pi*(x-t))*u(x, t) + b_tilde*pi*COS(2*pi*(x-t))
                Q_d_x = (Q_e_x - Qe(b_tilde, h_tilde, u(x, t))*pi*COS(2*pi*(x-t)))*EXP(-u(x, t))

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau(u(x, t)) - tau_c)/eps_tau)**2)*(2*u(x, t)        &
            &                                                                       *pi*COS(2*pi*(x-t)))
                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u(x, t)) - Qd(b_tilde, u(x, t)))*I_x

                F_b_x = Q_e_x*I_geq(tau(u(x, t))) + Q_d_x*I_leq(tau(u(x, t))) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u(x, t), tau(u(x, t)))

                R_tilde_xt = m_b_x*b_t + F_b_x - S_b

            case (2)
! Solution manufacturee stationnaire
! Terme de flux
                b_tilde = b(x, 0._pr)
                h_tilde = h(x, 0._pr)
                Q_e_x = (u(x, 0._pr) + b_tilde)*pi*COS(2*pi*x)
                Q_d_x = (Q_e_x - Qe(b_tilde, h_tilde, u(x, 0._pr))*pi*COS(2*pi*x))*EXP(-u(x,0._pr))

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau(u(x, 0._pr)) - tau_c)/eps_tau)**2)*              &
            &                                                               (2*u(x, 0._pr)*pi*COS(2*pi*x))
                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u(x, 0._pr)) - Qd(b_tilde, u(x, 0._pr)))*I_x

                F_b_x = Q_e_x*I_geq(tau(u(x, 0._pr))) + Q_d_x*I_leq(tau(u(x, 0._pr))) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u(x, 0._pr), tau(u(x, 0._pr)))

                R_tilde_xt = F_b_x - S_b

            case (3)
! Cas test bossu
                b_x = 0._pr
                do k = 1, 5
                    x_c_t = MODULO(x_c_k(k) + c_grass*t, 1._pr)

                    delta_x = x - x_c_t
                    if (delta_x >  0.5_pr) delta_x = delta_x - 1._pr
                    if (delta_x < -0.5_pr) delta_x = delta_x + 1._pr

                    b_k_xt = A_k(k)*EXP(-(delta_x)**2/(2._pr*sigma_k(k)**2))

                    b_x = b_x - (delta_x)/sigma_k(k)**2*b_k_xt
                end do

! Terme temporel
                m_b_x = m(b_x)
                b_t = -c_grass*b_x

! Terme de flux
                b_tilde = b(x, t)
                h_tilde = h(x, t)
                Q_e_x = b_x*u(x, t) + b_tilde*pi*COS(2*pi*x)
                Q_d_x = (Q_e_x - Qe(b_tilde, h_tilde, u(x, t))*pi*COS(2*pi*x))*EXP(-u(x, t))

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau(u(x, t)) - tau_c)/eps_tau)**2)*                  &
            &                                                               (2*u(x, t)*pi*COS(2*pi*x))
                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u(x, t)) - Qd(b_tilde, u(x, t)))*I_x

                F_b_x = Q_e_x*I_geq(tau(u(x, t))) + Q_d_x*I_leq(tau(u(x, t))) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u(x, t), tau(u(x, t)))

                R_tilde_xt = m_b_x*b_t + F_b_x - S_b

            case (4)
! Cas test loi de Grass
! Terme temporel
                b_x = 2*pi*COS(2*pi*(x-t))
                m_b_x = m(b_x)
                b_t = -2*pi*COS(2*pi*(x-t))

! Terme de flux
                b_tilde = b(x, t)
                h_tilde = h(x, t)
                u_tilde = u(x, t)
                tau_tilde = tau(u_tilde)
                u_x = pi*COS(2*pi*(x-t))
                tau_x = rho_f*SIGN(1._pr, u_tilde)*u_x*(C_lam + 2._pr*C_tur*ABS(u_tilde))

                Q_e_x = A_Grass*u_x*ABS(u_tilde)**(m_Grass-2)*(ABS(u_tilde) + (m_Grass-1)*SIGN(1._pr, u_tilde)*u_tilde)
                Q_d_x = Q_e_x

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau_tilde - tau_c)/eps_tau)**2)*tau_x

                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u_tilde) - Qd(b_tilde, u_tilde))*I_x

                F_b_x = Q_e_x*I_geq(tau_tilde) + Q_d_x*I_leq(tau_tilde) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u_tilde, tau_tilde)

                R_tilde_xt = m_b_x*b_t + F_b_x - S_b

            case (5)
! Cas MPM
! Terme temporel
                b_x = 2*pi*COS(2*pi*(x-t))
                m_b_x = m(b_x)
                b_t = -2*pi*COS(2*pi*(x-t))

! Terme de flux
                b_tilde = b(x, t)
                h_tilde = h(x, t)
                u_tilde = u(x, t)
                tau_tilde = tau(u_tilde)
                u_x = pi*COS(2*pi*(x-t))
                h_x = -h_tilde*u_x/u_tilde
                tau_x = rho_f*SIGN(1._pr, u_tilde)*u_x*(C_lam + 2._pr*C_tur*ABS(u_tilde))

                C = C_Chezy(h_tilde)
                C_prime = C_prime_Chezy(h_tilde)
                C_C_prime_h = 1._pr/C_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)*(C_prime -                                   &
            &                                                                        6._pr*SQRT(g_part)*5.75_pr/LOG(10._pr))
                mu_x = 3._pr/2*h_x*C_C_prime_h*(C/C_prime)**(1._pr/2)
                theta_prime_x = 1._pr/(rho_f*g_part*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu(h_tilde)*tau_x)

                Q_e_x = 8._pr*SQRT((s_rel-1._pr)*g_part*d_50**3)*1.5_pr*theta_prime_x*                                      &
            &                                               (MAX(0._pr, theta_prime(h_tilde, u_tilde) - theta_c))**(1._pr/2)
                Q_d_x = 0._pr

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau_tilde - tau_c)/eps_tau)**2)*tau_x

                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u_tilde) - Qd(b_tilde, u_tilde))*I_x

                F_b_x = Q_e_x*I_geq(tau_tilde) + Q_d_x*I_leq(tau_tilde) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u_tilde, tau_tilde)

                R_tilde_xt = m_b_x*b_t + F_b_x - S_b

            case (6)
! Cas Van Rijn
! Terme temporel
                b_x = 2*pi*COS(2*pi*(x-t))
                m_b_x = m(b_x)
                b_t = -2*pi*COS(2*pi*(x-t))

! Terme de flux
                b_tilde = b(x, t)
                h_tilde = h(x, t)
                u_tilde = u(x, t)
                tau_tilde = tau(u_tilde)
                u_x = pi*COS(2*pi*(x-t))
                h_x = -h_tilde*u_x/u_tilde
                tau_x = rho_f*SIGN(1._pr, u_tilde)*u_x*(C_lam + 2._pr*C_tur*ABS(u_tilde))

                C = C_Chezy(h_tilde)
                C_prime = C_prime_Chezy(h_tilde)
                C_C_prime_h = 1._pr/C_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)*(C_prime -                                   &
            &                                                                        6._pr*SQRT(g_part)*5.75_pr/LOG(10._pr))
                mu_x = 3._pr/2*h_x*C_C_prime_h*(C/C_prime)**(1._pr/2)
                theta_prime_x = 1._pr/(rho_f*g_part*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu(h_tilde)*tau_x)

                T_Van_Rijn = (theta_prime(h_tilde, u_tilde) - theta_c)/theta_c

                Q_e_x = 0._pr
                if (0._pr < T_Van_Rijn .AND. T_Van_Rijn < 3._pr) then
                    Q_e_x = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.053_pr*d_etoile**(-0.3_pr)*            &
                &                                        2.1_pr/theta_c*theta_prime_x*T_Van_Rijn**1.1_pr
                else if (3._pr <= T_Van_Rijn) then
                    Q_e_x = SQRT(g_part*(s_rel-1._pr)*d_50**3)*0.1_pr*d_etoile**(-0.3_pr)*              &
                &                              3._pr/(2._pr*theta_c)*theta_prime_x*T_Van_Rijn**(1._pr/2)
                end if
                Q_d_x = 0._pr

                I_x = 1._pr/(2*eps_tau)*(1._pr - TANH((tau_tilde - tau_c)/eps_tau)**2)*tau_x

                Qe_Qd_I_x = (Qe(b_tilde, h_tilde, u_tilde) - Qd(b_tilde, u_tilde))*I_x

                F_b_x = Q_e_x*I_geq(tau_tilde) + Q_d_x*I_leq(tau_tilde) + Qe_Qd_I_x

! Terme source
                S_b = S(b_tilde, h_tilde, u_tilde, tau_tilde)

                R_tilde_xt = m_b_x*b_t + F_b_x - S_b

            case default
                R_tilde_xt = 0._pr

            end select

        end function R_tilde

end module mod_fonctions