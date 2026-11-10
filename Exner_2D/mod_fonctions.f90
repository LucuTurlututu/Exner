module mod_fonctions

    use mod_parametres
    use mod_qte_physiques
    
    implicit none

    contains

        pure function m(nabla_b) result(m_nabla_b)
! ==================================================================================================== !
! Fonction renvoyant la metrique geometrique

! Entree de la fonction :
!   - nabla_b(1:2) : Gradient de la hauteur de la batymetrie b

! Sortie de la fonction :
!   - m_nabla_b : Metique evaluee en nabla_b
! ==================================================================================================== !

! Entree de la fonction
            real(kind = pr), dimension(1:2), intent(in)                 :: nabla_b

! Sortie de la fonction
            real(kind = pr)                                             :: m_nabla_b

            select case (cas)

            case (1, 2, 3, 4, 5, 6)
                m_nabla_b = 1._pr/SQRT(1._pr + norme_2(nabla_b)**2)

            case default
                m_nabla_b = 1._pr

            end select

        end function m

        pure function Q_d(b, h, u) result(Q_d_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le flux de depot
        
! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
    
! Sortie de la fonction :
!   - Q_d_b_h_u(1:2) : Fonction Q_d evaluee en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: Q_d_b_h_u

            select case (cas)

            case (1, 4, 7)
! Loi de Grass
                Q_d_b_h_u = Grass(u)

            case (2, 5, 8)
! Loi de MPM
                Q_d_b_h_u = 0._pr

            case (3, 6, 9)
! Loi de Van Rijn
                Q_d_b_h_u = 0._pr

            case default
                Q_d_b_h_u = 0._pr

            end select
    
        end function Q_d

        pure function Q_e(b, h, u) result(Q_e_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le flux d'erosion
        
! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
    
! Sortie de la fonction :
!   - Q_e_b_h_u(1:2) : Fonction Q_e evaluee en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: Q_e_b_h_u

            select case (cas)

            case (1, 4, 7)
! Loi de Grass
                Q_e_b_h_u = Grass(u)

            case (2, 5, 8)
! Loi de MPM
                Q_e_b_h_u = MPM(h, u)

            case (3, 6, 9)
! Loi de Van Rijn
                Q_e_b_h_u = Van_Rijn(h, u)

            case default
                Q_e_b_h_u = 0._pr

            end select
    
        end function Q_e

        pure function Sigmoide(val, ref, larg) result(Sigmoide_val_ref_larg)
! ==================================================================================================== !
! Fonction renvoyant la Sigmoide de lissage
        
! Entrees de la fonction :
!   - val : Grandeur sur laquelle on souhaite effectuer le lissage
!   - ref : Valeur de reference de la grandeur val
!   - larg : Parametere controlant la largeur de la transition
    
! Sortie de la fonction :
!   - Sigmoide_val_ref_larg : Sigmoide evaluee en (val, ref, larg)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: val, ref, larg

! Sortie de la fonction
            real(kind = pr)                                             :: Sigmoide_val_ref_larg

            Sigmoide_val_ref_larg = 1._pr/2*(1._pr + TANH((val - ref)/larg))

        end function Sigmoide

        pure function I_leq(val, ref, larg) result(I_leq_val_ref_larg)
! ==================================================================================================== !
! Fonction renvoyant l'indicatrice inferieure ou egale lissee
        
! Entrees de la fonction :
!   - val : Grandeur sur laquelle on souhaite effectuer le lissage
!   - ref : Valeur de reference de la grandeur val
!   - larg : Parametere controlant la largeur de la transition
    
! Sortie de la fonction :
!   - I_leq_val_ref_larg : Indicatrice inferieure ou egale lissee evaluee en (val, ref, larg)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: val, ref, larg

! Sortie de la fonction
            real(kind = pr)                                             :: I_leq_val_ref_larg

            I_leq_val_ref_larg = 1._pr - Sigmoide(val, ref, larg)

        end function I_leq

        pure function I_geq(val, ref, larg) result(I_geq_val_ref_larg)
! ==================================================================================================== !
! Fonction renvoyant l'indicatrice superieure ou egale lissee
        
! Entrees de la fonction :
!   - val : Grandeur sur laquelle on souhaite effectuer le lissage
!   - ref : Valeur de reference de la grandeur val
!   - larg : Parametere controlant la largeur de la transition
    
! Sortie de la fonction :
!   - I_geq_val_ref_larg : Indicatrice superieure ou egale lissee evaluee en (val, ref, larg)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: val, ref, larg

! Sortie de la fonction
            real(kind = pr)                                             :: I_geq_val_ref_larg

            I_geq_val_ref_larg = Sigmoide(val, ref, larg)

        end function I_geq

        pure function G(b, h, u) result(G_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le flux sedimentaire G

! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - G_b_h_u(1:2) : Flux evalue en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr), dimension(2)                               :: G_b_h_u

! Variable locale
            real(kind = pr)                                             :: u_e_Berthon

            select case (cas)

            case (7)
! Cas de Berthon et al. (2012) avec la loi de Grass
                u_e_Berthon = MAX(0._pr, DOT_PRODUCT(u, u) - u_cr_Berthon_G)

                G_b_h_u = A_Berthon_G*u_e_Berthon**(3._pr/2)

            case (8)
! Cas de Berthon et al. (2012) avec la loi de MPM
                u_e_Berthon = MAX(0._pr, DOT_PRODUCT(u, u) - u_cr_Berthon_MPM)

                G_b_h_u = A_Berthon_MPM*u_e_Berthon**(3._pr/2)

            case (9)
! Cas de Berthon et al. (2012) avec la loi de Van Rijn
                u_e_Berthon = MAX(0._pr, DOT_PRODUCT(u, u) - u_cr_Berthon_VR)

                if ((theta_prime(h, u) - theta_c) / theta_c < 3._pr) then
                    G_b_h_u = A_Berthon_VR_leq*u_e_Berthon**(2.1_pr)
                else
                    G_b_h_u = A_Berthon_VR_geq*u_e_Berthon**(1.5_pr)
                end if

            case default
                G_b_h_u = Q_d(b, h, u)*I_leq(theta_prime(h, u), theta_c, eps_theta)                     &
            &                                 + Q_e(b, h, u)*I_geq(theta_prime(h, u), theta_c, eps_theta)

            end select

        end function G

        pure function partial_h_n_G(h, u, n_F) result(partial_h_G_h_u_n_F)
! ==================================================================================================== !
! Fonction renvoyant la derivee par rapport a h du flux sedimentaire suivant la normale n_F

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
!   - n_F(1:2) : Normale sortante unitaire a la face F

! Sortie de la fonction :
!   - partial_h_G_h_u_n_F : Derivee de G par rapport a h evaluee en (h, u) suivant n_F
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u, n_F

! Sortie de la fonction
            real(kind = pr)                                             :: partial_h_G_h_u_n_F

! Variables locales
            real(kind = pr)                                             :: u_n, theta_prime_h_u, partial_h_theta_prime
            real(kind = pr)                                             :: A, C, T_VR, dh, h_plus, h_moins
            real(kind = pr), dimension(2)                               :: u_plus, u_moins

            u_n = DOT_PRODUCT(u, n_F)

            select case (cas)

            case (1, 4, 7)
! Loi de Grass
                partial_h_G_h_u_n_F = -m_Grass/h*A_Grass*norme_2(u)**(m_Grass-1)*u_n

            case (2, 5, 8)
! Loi de MPM
                theta_prime_h_u = theta_prime(h, u)

                dh = MAX(eps_DF*h, eps_DF)
                h_plus = h + dh ; h_moins = MAX(h - dh, eps_DF)

                u_plus  = (h/h_plus)*u ; u_moins = (h/h_moins)*u

                partial_h_theta_prime = (theta_prime(h_plus, u_plus) - theta_prime(h_moins, u_moins))                       &
            &                                                                              /MAX(2._pr*eps_DF*h, 2._pr*eps_DF)

                A = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u-theta_c, 0._pr)**0.5_pr            &
            &                                                                           *partial_h_theta_prime*u_n/norme_2(u)

                C = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                       &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                            *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*partial_h_theta_prime

                partial_h_G_h_u_n_F = A*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case (3, 6, 9)
! Loi de Van Rijn
                theta_prime_h_u = theta_prime(h, u)

                dh = MAX(eps_DF*h, eps_DF)
                h_plus = h + dh ; h_moins = MAX(h - dh, eps_DF)

                u_plus  = (h/h_plus)*u ; u_moins = (h/h_moins)*u

                partial_h_theta_prime = (theta_prime(h_plus, u_plus) - theta_prime(h_moins, u_moins))                       &
            &                                                                              /MAX(2._pr*eps_DF*h, 2._pr*eps_DF)

                T_VR = T_s(h, u)

                if (0._pr < T_VR .AND. T_VR < 3._pr) then
                    A = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*2.1_pr*T_VR**1.1_pr                &
                &                                                               *partial_h_theta_prime/theta_c*u_n/norme_2(u)

                    C = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)            &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**2.1_pr*partial_h_theta_prime

                else if (T_VR > 3._pr) then
                    A = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*T_VR**0.5_pr                  &
                &                                                               *partial_h_theta_prime/theta_c*u_n/norme_2(u)

                    C = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)              &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**1.5_pr*partial_h_theta_prime

                else
                    A = 0._pr ; C = 0._pr
                
                end if

                partial_h_G_h_u_n_F = A*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case default
                partial_h_G_h_u_n_F = 0._pr

            end select

        end function partial_h_n_G

        pure function partial_q_n_G(h, u, n_F) result(partial_q_G_h_u_n_F)
! ==================================================================================================== !
! Fonction renvoyant la derivee par rapport a q du flux sedimentaire suivant la normale n_F

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
!   - n_F(1:2) : Normale sortante unitaire a la face F

! Sortie de la fonction :
!   - partial_q_G_h_u_n_F : Derivee de G par rapport a q evaluee en (h, u) suivant n_F
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u, n_F

! Sortie de la fonction
            real(kind = pr)                                             :: partial_q_G_h_u_n_F

! Variables locales
            real(kind = pr)                                             :: u_n, theta_prime_h_u, partial_u_theta_prime
            real(kind = pr)                                             :: A, B, C, T_VR

            u_n = DOT_PRODUCT(u, n_F)

            select case (cas)

            case (1, 4, 7)
! Loi de Grass
                partial_q_G_h_u_n_F = 1._pr/h*A_Grass*norme_2(u)**(m_Grass-1)*(1._pr + (m_Grass-1)*u_n**2/norme_2(u)**2)

            case (2, 5, 8)
! Loi de MPM
                theta_prime_h_u = theta_prime(h, u)
                partial_u_theta_prime = mu_corr(h)/(rho_f*g_const*(s_rel-1._pr)*d_50)                                       &
            &                                                        *(rho_f*(C_lam + 2._pr*C_tur*norme_2(u))*u_n/norme_2(u))

                A = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u-theta_c, 0._pr)**0.5_pr            &
            &                                                                           *partial_u_theta_prime*u_n/norme_2(u)

                B = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr                 &
            &                                                                *1._pr/norme_2(u)*(1._pr - u_n**2/norme_2(u)**2)

                C = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                       &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                            *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*partial_u_theta_prime

                partial_q_G_h_u_n_F = (A + B)*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case (3, 6, 9)
                theta_prime_h_u = theta_prime(h, u)
                partial_u_theta_prime = mu_corr(h)/(rho_f*g_const*(s_rel-1._pr)*d_50)                                       &
            &                                                        *(rho_f*(C_lam + 2._pr*C_tur*norme_2(u))*u_n/norme_2(u))
                T_VR = T_s(h, u)

                if (0._pr < T_VR .AND. T_VR < 3._pr) then
                    A = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*2.1_pr*T_VR**1.1_pr                &
                &                                                               *partial_u_theta_prime/theta_c*u_n/norme_2(u)

                    B = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*T_VR**2.1_pr                       &
                &                                                            *1._pr/norme_2(u)*(1._pr - u_n**2/norme_2(u)**2)

                    C = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)            &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**2.1_pr*partial_u_theta_prime

                else if (T_VR > 3._pr) then
                    A = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*T_VR**0.5_pr                  &
                &                                                               *partial_u_theta_prime/theta_c*u_n/norme_2(u)

                    B = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*T_VR**1.5_pr                         &
                &                                                            *1._pr/norme_2(u)*(1._pr - u_n**2/norme_2(u)**2)

                    C = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)              &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**1.5_pr*partial_u_theta_prime

                else
                    A = 0._pr ; B = 0._pr ; C = 0._pr
                
                end if

                partial_q_G_h_u_n_F = (A + B)*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case default
                partial_q_G_h_u_n_F = 0._pr

            end select

        end function partial_q_n_G

        pure function partial_q_t_G(h, u, n_F, t_F) result(partial_q_G_h_u_t_F)
! ==================================================================================================== !
! Fonction renvoyant la derivee par rapport a q du flux sedimentaire suivant la tangente t_F

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
!   - n_F(1:2) : Normale sortante unitaire a la face F
!   - t_F(1:2) : Tangente a la normale sortante unitaire a la face F

! Sortie de la fonction :
!   - partial_q_G_h_u_t_F : Derivee de G par rapport a q evaluee en (h, u) suivant t_F
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u, n_F, t_F

! Sortie de la fonction
            real(kind = pr)                                             :: partial_q_G_h_u_t_F

! Variables locales
            real(kind = pr)                                             :: u_n, u_t, theta_prime_h_u, partial_u_theta_prime
            real(kind = pr)                                             :: A, B, C, T_VR

            u_n = DOT_PRODUCT(u, n_F)
            u_t = DOT_PRODUCT(u, t_F)

            select case (cas)

            case (1, 4, 7)
! Loi de Grass
                partial_q_G_h_u_t_F = 1._pr/h*A_Grass*norme_2(u)**(m_Grass-1)*(m_Grass-1)*u_n*u_t/norme_2(u)**2

            case (2, 5, 8)
! Loi de MPM
                theta_prime_h_u = theta_prime(h, u)
                partial_u_theta_prime = mu_corr(h)/(rho_f*g_const*(s_rel-1._pr)*d_50)                                       &
            &                                                        *(rho_f*(C_lam + 2._pr*C_tur*norme_2(u))*u_t/norme_2(u))

                A = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u-theta_c, 0._pr)**0.5_pr            &
            &                                                                           *partial_u_theta_prime*u_n/norme_2(u)

                B = -8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr                &
            &                                                                       *1._pr/norme_2(u)*(u_n*u_t/norme_2(u)**2)

                C = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                       &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                            *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*partial_u_theta_prime

                partial_q_G_h_u_t_F = (A + B)*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case (3, 6, 9)
! Loi de Van Rijn
                theta_prime_h_u = theta_prime(h, u)
                partial_u_theta_prime = mu_corr(h)/(rho_f*g_const*(s_rel-1._pr)*d_50)                                       &
            &                                                        *(rho_f*(C_lam + 2._pr*C_tur*norme_2(u))*u_t/norme_2(u))
                T_VR = T_s(h, u)

                if (0._pr < T_VR .AND. T_VR < 3._pr) then
                    A = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*2.1_pr*T_VR**1.1_pr                &
                &                                                               *partial_u_theta_prime/theta_c*u_n/norme_2(u)

                    B = -0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*T_VR**2.1_pr                      &
                &                                                                   *1._pr/norme_2(u)*(u_n*u_t/norme_2(u)**2)

                    C = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)            &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**2.1_pr*partial_u_theta_prime

                else if (T_VR > 3._pr) then
                    A = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*T_VR**0.5_pr                  &
                &                                                               *partial_u_theta_prime/theta_c*u_n/norme_2(u)

                    B = -0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*T_VR**1.5_pr                        &
                &                                                                   *1._pr/norme_2(u)*(u_n*u_t/norme_2(u)**2)

                    C = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)              &
                &                *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**1.5_pr*partial_u_theta_prime

                else
                    A = 0._pr ; B = 0._pr ; C = 0._pr
                
                end if

                partial_q_G_h_u_t_F = (A + B)*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

            case default
                partial_q_G_h_u_t_F = 0._pr

            end select

        end function partial_q_t_G

        function lambda_DV(h, u, n_F) result(lambda_DV_h_u_n_F)
! ==================================================================================================== !
! Fonction renvoyant la valeur propre de De Vries

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau
!   - n_F(1:2) : Normale sortante unitaire a la face F

! Sortie de la fonction :
!   - lambda_DV_h_u_n_F : Valeur propre de De Vries
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u, n_F

! Sortie de la fonction
            real(kind = pr)                                             :: lambda_DV_h_u_n_F

! Variable locale
            real(kind = pr)                                             :: u_n, u_t, Fr, num, denom
            real(kind = pr)                                             :: alpha_n, delta_n, delta_t
            real(kind = pr), dimension(2)                               :: t_F

            t_F = [-n_F(2), n_F(1)]

            if (h > eps_var) then
                u_n = DOT_PRODUCT(u, n_F) ; u_t = DOT_PRODUCT(u, t_F)

                alpha_n = partial_h_n_G(h, u, n_F)
                delta_n = partial_q_n_G(h, u, n_F) ; delta_t = partial_q_t_G(h, u, n_F, t_F)

                Fr = ABS(u_n)/SQRT(g_const*h)

                Fr_max = MAX(Fr_max, Fr)

                num = (alpha_n + u_t*delta_t)*u_n
                denom = (-u_n*(1._pr - Fr**2) + alpha_n + u_t*delta_t - u_n*delta_n)

                if (ABS(denom) > eps_var) then
                    lambda_DV_h_u_n_F = num/denom
                
                else
                    lambda_DV_h_u_n_F = 0._pr

                end if

            else
                lambda_DV_h_u_n_F = 0._pr

            end if

        end function lambda_DV

        pure function k_s_total(h, u) result(k_s_total_h_u)
! ==================================================================================================== !
! Fonction renvoyant la rugosite equivalente totale de Nikuradse

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction
!   - k_s_total_h_u : Rugosite equivalente de Nikuradse
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: k_s_total_h_u

! Variables locales
            real(kind = pr)                                             :: Delta_r, lambda_r, Delta_d, lambda_d

! Grandeurs liees aux rides
            Delta_r = 250._pr*d_50
            lambda_r = 1000._pr*d_50

! Grandeurs liees aux dunes
            Delta_d = 0.11_pr*h*(d_50/h)**(0.3_pr)*(1._pr - EXP(-0.5_pr*T_s(h, u)))*(25._pr - T_s(h, u))
            lambda_d = 7.3_pr*h

            k_s_total_h_u = 2.5_pr*d_50 + 28._pr*(Delta_r**2/lambda_r) + 1.1_pr*Delta_d*(1._pr - EXP(-25._pr*Delta_d/lambda_d))

        end function k_s_total

        pure function c_Delta(h, u) result(c_Delta_h_u)
! ==================================================================================================== !
! Fonction representant la concentration volumique de reference a l'interface entre la couche de 
! charriage et la zone de suspension

! Entrees de la fonction :
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction
!   - c_Delta_h_u : Concentration volumique de reference
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: c_Delta_h_u

! Variable locale
            real(kind = pr)                                             :: Delta_ref

            Delta_ref = MAX(k_s_total(h, u), 0.01_pr*h)

            c_Delta_h_u = 0.015_pr*(d_50/Delta_ref)*(T_s(h, u)**(1.5_pr)/d_etoile**(0.3_pr))

        end function c_Delta

        pure function D(b, h, u) result(D_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le terme source de depot D

! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - D_b_h_u : Terme source D evalue en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: D_b_h_u

            D_b_h_u = w_0*c_Delta(h, u)*MAX(1._pr - theta_prime(h, u)/theta_c, 0._pr)

        end function D

        pure function E(b, h, u) result(E_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le terme source de d'erosion E

! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - E_b_h_u : Terme source E evalue en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: E_b_h_u

            E_b_h_u = w_0*c_Delta(h, u)

        end function E

        pure function S(b, h, u) result(S_b_h_u)
! ==================================================================================================== !
! Fonction renvoyant le terme source S

! Entrees de la fonction :
!   - b : Hauteur de la bathymetrie
!   - h : Hauteur d'eau
!   - u(1:2) : Vitesse d'eau

! Sortie de la fonction :
!   - S_b_h_u : Terme source S evalue en (b, h, u)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: b, h
            real(kind = pr), dimension(2), intent(in)                   :: u

! Sortie de la fonction
            real(kind = pr)                                             :: S_b_h_u

            select case (cas)

            case (7, 8, 9, 10)
                S_b_h_u = 0._pr

            case default
                S_b_h_u = D(b, h, u)*I_leq(theta_prime(h, u), theta_c, eps_theta)                       &
            &                                   - E(b, h, u)*I_geq(theta_prime(h, u), theta_c, eps_theta)

            end select

        end function S

        function R_tilde(Coord_phy, t) result(R_tilde_x_y_t)
! ==================================================================================================== !
! Fonction renvoyant le terme source supplementaire pour la solutin manufacturee

! Entrees de la fonction :
!   - Coord_phy(1:2) : Coordonnees (x, y)
!   - t : Temps courant

! Sortie de la fonction :
!   - R_tilde_x_y_t : Fonction R_tilde evaluee en (x, y, t)

! N.B. Il ne faut pas inclure la metrique directement ici ! Il est necessaire d'utiliser
! l'approximation numerique du gradient de b pour que la matrice de masse implicite
! ait un sens pour le solveur de Newton
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)                                 :: t
            real(kind = pr), dimension(2), intent(in)                   :: Coord_phy

! Sorties de la fonction
            real(kind = pr)                                             :: R_tilde_x_y_t

! Variables locales
            real(kind = pr)                                             :: X, Y, b_x, b_y, b_t, theta_prime_h_u, theta_prime_x
            real(kind = pr)                                             :: u_1, u_x, u_y, h_tilde, div_G, coeff_2D, b_tilde, A, C
            real(kind = pr)                                             :: C_h, C_h_prime, partial_h_C_h_C_h_prime, mu_x, T_VR
            real(kind = pr)                                             :: tau_tilde, h_x, h_y, tau_x, theta_prime_y, norme_u_y
            real(kind = pr)                                             :: A_x, A_y, C_x, C_y, norme_u, tau_y, mu_y, norme_u_x
            real(kind = pr), dimension(2)                               :: nabla_b_tilde, u_tilde

            select case (cas)

            case (1)
! Cas manufacture 1D : Loi de Grass
                X = 2._pr*pi*(Coord_phy(1) - t)
                b_tilde = b(Coord_phy, t)
                b_x = 2._pr*pi*COS(X)
                b_y = 0._pr
                b_t = -b_x

                u_1 = 1._pr + 1._pr/2*b_tilde
                u_x = 0.5_pr*b_x

                div_G = A_Grass*u_x*ABS(u_1)**(m_Grass-2)*(ABS(u_1) + (m_Grass-1)*SIGN(1._pr, u_1)*u_1)

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))

            case (2)
! Cas manufacture 1D : Loi de MPM
                X = 2._pr*pi*(Coord_phy(1) - t)
                b_tilde = b(Coord_phy, t)
                b_x = 2._pr*pi*COS(X)
                b_y = 0._pr
                b_t = -b_x

                u_1 = 1._pr + 1._pr/2*b_tilde
                u_x = 1._pr/2*b_x

                h_tilde = h(Coord_phy, t)
                h_x = -SIGN(1._pr, u_1)*u_x/(ABS(u_1))**2

                tau_tilde = rho_f*(C_lam + C_tur*u_1)*u_1
                tau_x = rho_f*(C_lam + 2._pr*C_tur*u_1)*u_x

                theta_prime_h_u = theta_prime(h(Coord_phy, t), u(Coord_phy, t))

                C_h = C_Chezy(h_tilde)
                C_h_prime = C_prime_Chezy(h_tilde)
                partial_h_C_h_C_h_prime = 1._pr/C_h_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)                                &
            &                                                          *(C_h_prime - 6._pr*SQRT(g_const)*5.75_pr/LOG(10._pr))
                mu_x = 3._pr/2*h_x*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)
                theta_prime_x = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu_corr(h_tilde)*tau_x)

                A = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u-theta_c, 0._pr)**0.5_pr*theta_prime_x

                C = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                       &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                                    *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*theta_prime_x

                div_G = A*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))

             case (3)
! Cas manufacture 1D : Loi de Van Rijn
                X = 2._pr*pi*(Coord_phy(1) - t)
                b_tilde = b(Coord_phy, t)
                b_x = 2._pr*pi*COS(X)
                b_y = 0._pr
                b_t = -b_x

                u_1 = 1._pr + 1._pr/2*b_tilde
                u_x = 1._pr/2*b_x

                h_tilde = h(Coord_phy, t)
                h_x = -SIGN(1._pr, u_1)*u_x/(ABS(u_1))**2

                tau_tilde = rho_f*(C_lam + C_tur*u_1)*u_1
                tau_x = rho_f*(C_lam + 2._pr*C_tur*u_1)*u_x

                theta_prime_h_u = theta_prime(h(Coord_phy, t), u(Coord_phy, t))

                C_h = C_Chezy(h_tilde)
                C_h_prime = C_prime_Chezy(h_tilde)
                partial_h_C_h_C_h_prime = 1._pr/C_h_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)                                &
            &                                                          *(C_h_prime - 6._pr*SQRT(g_const)*5.75_pr/LOG(10._pr))

                mu_x = 3._pr/2*h_x*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)
                theta_prime_x = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu_corr(h_tilde)*tau_x)

                T_VR = T_s(h_tilde, [u_1, 0._pr])


                if (0._pr < T_VR .AND. T_VR < 3._pr) then
                    A = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*2.1_pr*1._pr/theta_c               &
                &                                                                               *theta_prime_x*T_VR**(1.1_pr)

                    C = 0.053_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)            &
                &                      *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*theta_prime_x*T_VR**(2.1_pr)

                else if (T_VR > 3._pr) then
                    A = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*1._pr/theta_c                 &
                &                                                                               *theta_prime_x*T_VR**(0.5_pr)

                    C = 0.1_pr*D_etoile**(-0.3_pr)*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)              &
                &                      *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*theta_prime_x*T_VR**(1.5_pr)

                else
                    A = 0._pr ; C = 0._pr

                end if

                div_G = A*I_geq(theta_prime_h_u, theta_c, eps_theta) + C

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))

            case (4)
! Cas manufacture 2D : Loi de Grass
                X = 2._pr*pi*(Coord_phy(1) - t)
                Y = 2._pr*pi*(Coord_phy(2) - t)

                b_x = 2._pr*pi*COS(X)*SIN(Y)
                b_y = 2._pr*pi*SIN(X)*COS(Y)
                b_t = -(b_x + b_y)

                u_1 = 1._pr/2._pr*(1._pr + 1._pr/2._pr*b(Coord_phy, t))
                
                coeff_2D = 2._pr**(0.5_pr*(m_Grass-1._pr))

                div_G = A_Grass*m_Grass*coeff_2D*(u_1**(m_Grass-1._pr))*(0.25_pr*(b_x + b_y))

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))

            case (5)
! Cas manufacture 2D : Loi de MPM
                X = 2._pr*pi*(Coord_phy(1) - t)
                Y = 2._pr*pi*(Coord_phy(2) - t)

                b_tilde = b(Coord_phy, t)
                b_x = 2._pr*pi*COS(X)*SIN(Y)
                b_y = 2._pr*pi*SIN(X)*COS(Y)
                b_t = -(b_x + b_y)
                nabla_b_tilde = [b_x, b_y]

                u_tilde = u(Coord_phy, t)
                u_1 = 1._pr/2._pr*(1._pr + 0.5_pr*b(Coord_phy, t))
                u_x = 1._pr/4*b_x ; u_y = 1._pr/4*b_y
                norme_u = norme_2(u(Coord_phy, t))
                norme_u_x = SQRT(2._pr)*u_x ; norme_u_y = SQRT(2._pr)*u_y

                h_tilde = h(Coord_phy, t)
                h_x = -norme_u_x/(norme_u**2)
                h_y = -norme_u_y/(norme_u**2)

                C_h = C_Chezy(h_tilde)
                C_h_prime = C_prime_Chezy(h_tilde)
                partial_h_C_h_C_h_prime = 1._pr/C_h_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)                                &
            &                                                          *(C_h_prime - 6._pr*SQRT(g_const)*5.75_pr/LOG(10._pr))
                mu_x = 3._pr/2*h_x*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)
                mu_y = 3._pr/2*h_y*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)

                tau_tilde = rho_f*(C_lam + C_tur*norme_u)*norme_u
                tau_x = rho_f*(C_lam + 2._pr*C_tur*norme_u)*norme_u_x
                tau_y = rho_f*(C_lam + 2._pr*C_tur*norme_u)*norme_u_y

                theta_prime_h_u = theta_prime(h_tilde, u_tilde)

                theta_prime_x = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu_corr(h_tilde)*tau_x)
                theta_prime_y = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_y + mu_corr(h_tilde)*tau_y)

                A_x = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u - theta_c, 0._pr)**0.5_pr*theta_prime_x
                C_x = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                     &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                                    *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*theta_prime_x

                A_y = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1.5_pr*MAX(theta_prime_h_u-theta_c, 0._pr)**0.5_pr*theta_prime_y
                C_y = 8._pr*SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)                                     &
            &                               *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)                       &
            &                                                    *MAX(theta_prime_h_u - theta_c, 0._pr)**1.5_pr*theta_prime_y

                div_G = (A_x*I_geq(theta_prime_h_u, theta_c, eps_theta) + C_x)*(u_1/norme_u)                                &
            &                                          + (A_y*I_geq(theta_prime_h_u, theta_c, eps_theta) + C_y)*(u_1/norme_u)

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))

            case (6)
! Cas manufacture 2D : Loi de Van Rijn
                X = 2._pr*pi*(Coord_phy(1) - t)
                Y = 2._pr*pi*(Coord_phy(2) - t)

                b_tilde = b(Coord_phy, t)
                b_x = 2._pr*pi*COS(X)*SIN(Y)
                b_y = 2._pr*pi*SIN(X)*COS(Y)
                b_t = -(b_x + b_y)

                u_tilde = u(Coord_phy, t)
                u_1 = 1._pr/2._pr*(1._pr + 0.5_pr*b(Coord_phy, t))
                u_x = 1._pr/4*b_x ; u_y = 1._pr/4*b_y
                norme_u = norme_2(u(Coord_phy, t))
                norme_u_x = SQRT(2._pr)*u_x ; norme_u_y = SQRT(2._pr)*u_y

                h_tilde = h(Coord_phy, t)
                h_x = -norme_u_x/(norme_u**2)
                h_y = -norme_u_y/(norme_u**2)

                C_h = C_Chezy(h_tilde)
                C_h_prime = C_prime_Chezy(h_tilde)
                partial_h_C_h_C_h_prime = 1._pr/C_h_prime**2*1._pr/6*K_s*h_tilde**(-5._pr/6)                                &
            &                                                          *(C_h_prime - 6._pr*SQRT(g_const)*5.75_pr/LOG(10._pr))
                mu_x = 3._pr/2*h_x*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)
                mu_y = 3._pr/2*h_y*partial_h_C_h_C_h_prime*(C_h/C_h_prime)**(1._pr/2)

                tau_tilde = rho_f*(C_lam + C_tur*norme_u)*norme_u
                tau_x = rho_f*(C_lam + 2._pr*C_tur*norme_u)*norme_u_x
                tau_y = rho_f*(C_lam + 2._pr*C_tur*norme_u)*norme_u_y

                theta_prime_h_u = theta_prime(h_tilde, u_tilde)

                theta_prime_x = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_x + mu_corr(h_tilde)*tau_x)
                theta_prime_y = 1._pr/(rho_f*g_const*(s_rel-1._pr)*d_50)*(tau_tilde*mu_y + mu_corr(h_tilde)*tau_y)

                T_VR = T_s(h_tilde, u_tilde)

                if (0._pr < T_VR .AND. T_VR < 3._pr) then
                    A_x = SQRT(g_const*(s_rel-1._pr)*d_50**3)*0.053_pr*D_etoile**(-0.3_pr)*2.1_pr*theta_prime_x/theta_c     &
                &                                                                                               *T_VR**1.1_pr
                    C_x = SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)*0.053_pr*D_etoile**(-0.3_pr)          &
                &                        *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**2.1_pr*theta_prime_x

                    A_y = SQRT(g_const*(s_rel-1._pr)*d_50**3)*0.053_pr*D_etoile**(-0.3_pr)*2.1_pr*theta_prime_y/theta_c     &
                &                                                                                               *T_VR**1.1_pr
                    C_y = SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)*0.053_pr*D_etoile**(-0.3_pr)          &
                &                        *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**2.1_pr*theta_prime_y

                else if (T_VR > 3._pr) then
                    A_x = SQRT(g_const*(s_rel-1._pr)*d_50**3)*0.1_pr*D_etoile**(-0.3_pr)*1.5_pr*theta_prime_x/theta_c       &
                &                                                                              *T_VR**0.5_pr
                    C_x = SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)*0.1_pr*D_etoile**(-0.3_pr)            &
                &                        *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**1.5_pr*theta_prime_x

                    A_y = SQRT(g_const*(s_rel-1._pr)*d_50**3)*0.1_pr*D_etoile**(-0.3_pr)*1.5_pr*theta_prime_y/theta_c       &
                &                                                                                               *T_VR**0.5_pr
                    C_y = SQRT(g_const*(s_rel-1._pr)*d_50**3)*1._pr/(2._pr*eps_theta)*0.1_pr*D_etoile**(-0.3_pr)            &
                &                        *(1._pr - TANH((theta_prime_h_u - theta_c)/eps_theta)**2)*T_VR**1.5_pr*theta_prime_y

                else
                    A_x = 0._pr ; C_x = 0._pr ; A_y = 0._pr ; C_y = 0._pr
                
                end if

                div_G = (A_x*I_geq(theta_prime_h_u, theta_c, eps_theta) + C_x)*(u_1/norme_u)                                &
            &                                          + (A_y*I_geq(theta_prime_h_u, theta_c, eps_theta) + C_y)*(u_1/norme_u)

                nabla_b_tilde(1) = b_x ; nabla_b_tilde(2) = b_y
                R_tilde_x_y_t = m(nabla_b_tilde)*b_t + div_G - S(b(Coord_phy, t), h(Coord_phy, t), u(Coord_phy, t))
    
            case default
                R_tilde_x_y_t = 0._pr

            end select

        end function R_tilde

end module mod_fonctions