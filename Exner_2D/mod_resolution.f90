module mod_resolution

    use mod_parametres
    use mod_qte_physiques
    use mod_fonctions
    use mod_types_derives
    use mod_base
    use mod_interpolation
    use mod_limiteur

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_resolution :
!   - 
! ==================================================================================================== !

    contains

        subroutine Point_Fixe(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I, R_i, coeff_I,             &
        &                     B_i_q, B_i_qp1, P_i_q)
! ==================================================================================================== !
! Subroutine construisant le point fixe utilise pour l'algorithme de Picard et l'acceleration 
! d'Anderson

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
    
!   - i : Entier referant a la maille Omega_i
!   - t_s_I : Temps courant implicite
!   - R_i(0:N_ddl_i-1) : Second membre de la methode IMEX a l'etape s
!   - coeff_I : Coefficient associe au terme implicite dans le membre de gauche
!   - B_i_q(0:N_ddl_i-1) : Bathymetrie a l'etape q de l'algorithme iteratif

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner

! Sorties de la subroutine :
!   - B_i_qp1(0:N_ddl_i-1) : Bathymetrie a l'etape q+1 de l'algorithme iteratif
!   - P_i_q(0:N_ddl_i-1) : Vecteur residu de l'etape q
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t_s_I, coeff_I
            real(kind = pr), dimension(:), intent(in)                   :: R_i, B_i_q

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Sorties de la subroutine
            real(kind = pr), dimension(:), intent(out)                  :: B_i_qp1, P_i_q

! Variable locale
            integer                                                     :: N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%B_i(0:N_ddl_i-1, i) = B_i_q

            call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I)

! Point fixe classique
            B_i_qp1 = R_i + coeff_I*Calc_Exner%N_Omega_i(i)%L_i_I

! Residu
            P_i_q = B_i_qp1 - B_i_q

        end subroutine Point_Fixe

        subroutine Picard(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, i, t_s_I, R_i, coeff_I, e_q)
! ==================================================================================================== !
! Subroutine implementant l'algorithme de Picard pour la resolution implicite

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
    
!   - i : Entier referant a la maille Omega_i
!   - t_s_I : Temps courant implicite
!   - R_i(0:N_ddl_i-1) : Second membre de la methode IMEX a l'etape s
!   - coeff_I : Coefficient associe au terme implicite dans le membre de gauche

! Entrees/Sorties de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi

! Sortie de la subroutine :
!   - e_q : Critere d'arret de la methode de Newton, i.e. l'ecart relatif en norme L^2
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t_s_I, coeff_I
            real(kind = pr), dimension(:), intent(in)                   :: R_i

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
    
! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            type(Schema_Temporel), intent(inout)                        :: RK

! Sortie de la subroutine
            real(kind = pr), intent(out)                                :: e_q

! Variables locales
            integer                                                     :: q, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            BLOCK
                real(kind = pr), dimension(0:N_ddl_i-1)                 :: B_i_q, B_i_qp1, P_i_q

                B_i_q = Calc_Exner%B_i(0:N_ddl_i-1, i)

                do q = 1, N_max_Picard
                    call Point_Fixe(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I, R_i, coeff_I, B_i_q, B_i_qp1, P_i_q)

                    e_q = SQRT(DOT_PRODUCT(P_i_q, P_i_q))

                    B_i_q = B_i_qp1

                    if (e_q < eps_Newton) exit
                end do
                Calc_Exner%B_i(0:N_ddl_i-1, i) = B_i_q

                RK%q_max_Picard = MAX(RK%q_max_Picard, q)
            END BLOCK

        end subroutine Picard

        subroutine Newton(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, s, i, t_s_I, R_i, coeff_I)
! ==================================================================================================== !
! Subroutine implementant l'algorithme de Newton pour la resolution implicite

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - s : Entier referant a l'etape s du schema IMEX
!   - i : Entier referant a la maille Omega_i
!   - t_s_I : Temps courant implicite
!   - R_i(0:N_ddl_i-1) : Second membre de la methode IMEX a l'etape s
!   - coeff_I : Coefficient associe au terme implicite dans le membre de gauche

! Entrees/Sorties de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: s, i

            real(kind = pr), intent(in)                                 :: t_s_I, coeff_I
            real(kind = pr), dimension(:), intent(in)                   :: R_i

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            
! Entrees/Sorties de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            type(Schema_Temporel), intent(inout)                        :: RK

! Variables locales
            integer                                                     :: q, k, info, nrhs, seuil, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            BLOCK
                logical                                                 :: switch_Picard

                real(kind = pr)                                         :: e_q, q_moyen_Newton
                real(kind = pr), dimension(0:N_ddl_i-1)                 :: B_i_q, B_peps, B_meps, N_i_q, N_peps, N_meps, dB_i_q

                nrhs = 1
                switch_Picard = .FALSE.

                if (RK%N_appels_Newton > 0) then
                    q_moyen_Newton = real(RK%q_Newton, pr)/real(RK%N_appels_Newton, pr)
                    seuil = MAX(3, INT(q_moyen_Newton) + 2)
                else
                    seuil = N_max_Newton
                end if

                RK%N_appels_Newton = RK%N_appels_Newton + 1

                do q = 1, N_max_Newton

! Calcul du residu a l'etat courant B_i_q
                    call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I)
                    N_i_q = Calc_Exner%B_i(0:N_ddl_i-1, i) - coeff_I*Calc_Exner%N_Omega_i(i)%L_i_I - R_i

                    e_q = SQRT(DOT_PRODUCT(N_i_q, N_i_q))
                    if (e_q < eps_Newton) exit

                    if (q > seuil) then
                        switch_Picard = .TRUE.
                        exit
                    end if

! Construction de la matrice Jacobienne
! N.B. Comme on utilise un schema SDIRK, la Jacobienne reste identique pour toutes les etapes
! s du schema
                    if (s == 1 .AND. q == 1) then 
                        B_i_q = Calc_Exner%B_i(0:N_ddl_i-1, i)

                        do k = 0, N_ddl_i-1
! Perturbation positive
                            B_peps = B_i_q
                            B_peps(k) = B_peps(k) + eps_DF
                            Calc_Exner%B_i(0:N_ddl_i-1, i) = B_peps

                            call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I)
                            N_peps = B_peps - coeff_I*Calc_Exner%N_Omega_i(i)%L_i_I - R_i

! Perutbation negative
                            B_meps = B_i_q
                            B_meps(k) = B_meps(k) - eps_DF
                            Calc_Exner%B_i(0:N_ddl_i-1, i) = B_meps

                            call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I)
                            N_meps = B_meps - coeff_I*Calc_Exner%N_Omega_i(i)%L_i_I - R_i

! Matrice Jacobienne
                            Calc_Exner%N_Omega_i(i)%J_i_q(:, k) = (N_peps - N_meps)/(2._pr*eps_DF)
                        end do
                        Calc_Exner%B_i(0:N_ddl_i-1, i) = B_i_q

                        call DGETRF(N_ddl_i, N_ddl_i, Calc_Exner%N_Omega_i(i)%J_i_q, N_ddl_i, Calc_Exner%N_Omega_i(i)%ipiv, info)

                        if (info /= 0) then
                            switch_Picard = .TRUE.
                            exit
                        end if
                    end if

                    dB_i_q = -N_i_q

                    call DGETRS("N", N_ddl_i, nrhs, Calc_Exner%N_Omega_i(i)%J_i_q, N_ddl_i,                                     &
                &                                                            Calc_Exner%N_Omega_i(i)%ipiv, dB_i_q, N_ddl_i, info)

                    if (info /= 0) then
                        switch_Picard = .TRUE.
                        exit
                    end if

! Mise a jour a l'etape q
                    Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%B_i(0:N_ddl_i-1, i) + dB_i_q
                end do

                if (.NOT. switch_Picard) RK%q_max_Newton = MAX(RK%q_max_Newton, q) ; RK%q_Newton = RK%q_Newton + q

! Algorithme de Picard
                if (switch_Picard) then
                    call Picard(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, i, t_s_I, R_i, coeff_I, e_q)
                end if

                if (e_q >= eps_Newton) then
                    print *, "Attention : Le solveur Hybride n'a pas converge a la maille ", i, " Erreur = ", e_q
                end if

            END BLOCK

        end subroutine Newton

        subroutine make_L_i_E(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant l'operateur de la methode des lignes explicite L_i_E

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - L_i_E(0:N_ddl_i-1) : Composante i de l'operateur explicite de la methode des maille
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i
            
            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

            Calc_Exner%N_Omega_i(i)%L_i_E = 0._pr

! Construction de M_i
            call make_M_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de I_i_EF
            call make_I_i_EF(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de I_i_VF
            call make_I_i_VF(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de I_i_R_tilde
            call make_I_i_R_tilde(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de L_i_E
            call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i,                                                         &
        &     Calc_Exner%N_Omega_i(i)%I_i_EF - Calc_Exner%N_Omega_i(i)%I_i_VF + Calc_Exner%N_Omega_i(i)%I_i_R_tilde,        &
        &                                                                                      Calc_Exner%N_Omega_i(i)%L_i_E)

        end subroutine make_L_i_E

        subroutine make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant l'operateur de la methode des lignes implicite L_i_I

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - L_i_I(0:N_ddl_i-1) : Composante i de l'operateur implicite de la methode des maille
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i
            
            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

            Calc_Exner%N_Omega_i(i)%L_i_I = 0._pr

! Construction de M_i
            call make_M_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de I_i_S
            call make_I_i_S(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

! Construction de L_i_I
            call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i,                                                         &
        &                                                       Calc_Exner%N_Omega_i(i)%I_i_S, Calc_Exner%N_Omega_i(i)%L_i_I)

        end subroutine make_L_i_I

        subroutine find_dt(Mesh, Quad_ref, Base_ref, Calc_Exner, t)
! ==================================================================================================== !
! Subroutine calculant le pas de temps optimal

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
! ==================================================================================================== !

! Entrees de la subroutine
            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: i, j, F, N_ddl_i

            BLOCK
                real(kind = pr)                                         :: lambda_i_n, lambda_tilde, mu_i_n
                real(kind = pr)                                         :: b_j, b_F, h_j, h_F, xi_j, eta_j, S_j, S_j_eps, dt_A, dt_S
                real(kind = pr), dimension(2)                           :: Coord_phy, Coord_phy_F, u_j, u_F, nabla_b_j

                dt_A = T_final ; dt_S = T_final

                !$omp parallel do private(i, j, F, N_ddl_i, lambda_i_n, lambda_tilde, mu_i_n) &
                !$omp private(b_j, b_F, h_j, h_F, xi_j, eta_j, S_j, S_j_eps, Coord_phy, Coord_phy_F, u_j, u_F, nabla_b_j) &
                !$omp shared(Mesh, Quad_ref, Base_ref, Calc_Exner, t) &
                !$omp reduction(min: dt_A, dt_S)

                do i = 1, Mesh%N_elements
                    N_ddl_i = Base_ref%N_ddl_i(i)

                    lambda_i_n = 0._pr ; lambda_tilde = 0._pr ; mu_i_n = 0._pr

                    BLOCK
                        real(kind = pr), dimension(1:Quad_ref%N_G_Face)         :: B_i_F, B_r_F

! Pas de temps lie a l'advection
                        call make_Beta_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)

                        do F = 1, Mesh%N_Sommets_Omega_i(i)
                            call reconstruction_locale(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, B_i_F, B_r_F)

                            do j = 1, Quad_ref%N_G_Face
                                xi_j = Quad_ref%Coord_Ref_Face(j)

                                Coord_phy_F = hat_phi_F(Mesh, i, F, xi_j)

                                b_F = B_i_F(j)
                                h_F = h(Coord_phy_F, t)
                                u_F = u(Coord_phy_F, t)

                                if (Mesh%N_Sommets_Omega_i(i) == 3) then
                                    nabla_b_j(1) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j, F),                                     &
                                &                                                 Calc_Exner%N_Omega_i(i)%Beta_i(0:N_ddl_i-1, 1))
                                    nabla_b_j(2) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j, F),                                     &
                                &                                                 Calc_Exner%N_Omega_i(i)%Beta_i(0:N_ddl_i-1, 2))

                                else if (Mesh%N_Sommets_Omega_i(i) == 4) then
                                    nabla_b_j(1) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j, F),                                    &
                                &                                                 Calc_Exner%N_Omega_i(i)%Beta_i(0:N_ddl_i-1, 1))
                                    nabla_b_j(2) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j, F),                                    &
                                &                                                 Calc_Exner%N_Omega_i(i)%Beta_i(0:N_ddl_i-1, 2))
                                end if

! Compensation liee au terme manufacture
                                if (norme_2(nabla_b_j) > eps_var) then
                                    lambda_tilde = (S(b(Coord_phy_F, t), h_F, u_F) + R_tilde(Coord_phy_F, t))/norme_2(nabla_b_j)
                                else
                                    lambda_tilde = 0._pr
                                end if

                                lambda_i_n = MAX(lambda_i_n, 1._pr/m(nabla_b_j)*(ABS(lambda_DV(h_F, u_F, Mesh%n_F(:, F, i))))   &
                            &                                                                                + ABS(lambda_tilde))
                            end do
                        end do

                        if (lambda_i_n > eps_var) then
                            dt_A = MIN(dt_A, 1._pr/(2*p+1)*Mesh%Diam_Omega_i(i)/lambda_i_n)
                        end if

! Pas de temps lie au terme source
                        if (schema == 1) then
                            if (Mesh%N_sommets_Omega_i(i) == 3) then
                                do j = 1, Quad_ref%N_G_Tri
                                    xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j =  Quad_ref%Coord_Ref_Tri(2, j)

                                    Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                                    b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                                    h_j = h(Coord_phy, t)
                                    u_j = u(Coord_phy, t)

                                    S_j = S(b_j, h_j, u_j)
                                    S_j_eps = S(b_j + eps_DF, h_j, u_j)

                                    mu_i_n = MAX(mu_i_n, (S_j_eps - S_j)/eps_DF)
                                end do

                            else if (Mesh%N_sommets_Omega_i(i) == 4) then
                                do j = 1, Quad_ref%N_G_Quad
                                    xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j =  Quad_ref%Coord_Ref_Quad(2, j)

                                    Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                                    b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                                    h_j = h(Coord_phy, t)
                                    u_j = u(Coord_phy, t)

                                    S_j = S(b_j, h_j, u_j)
                                    S_j_eps = S(b_j + eps_DF, h_j, u_j)

                                    mu_i_n = MAX(mu_i_n, (S_j_eps - S_j)/eps_DF)
                                end do
                            end if

                            if (mu_i_n > eps_var) then
                                dt_S = MIN(dt_S, 1._pr/mu_i_n)
                            end if
                        end if
                    END BLOCK
                end do

                select case (schema)

                case (3)
                    Calc_Exner%dt = alpha_CFL*dt_A
                
                case default
                    Calc_Exner%dt = alpha_CFL*MIN(dt_A, dt_S)
                
                end select

            END BLOCK

        end subroutine find_dt

        subroutine make_Butcher_SSP_2_2(RK)
! ==================================================================================================== !
! Subroutine initialisant les coefficients du tableau de Butcher pour le schema SSP(2,2)
            
! Entree/Sortie de la subroutine :
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
!   - A_Butcher(1:m_etapes, 1:m_etapes) : Matrice A de Butcher explicite triangulaire 
! strictement inferieure
!   - A_tilde_Butcher(1:m_etapes, 1:m_etapes) : Matrice A_tilde de Butcher implicite triangulaire
! inferieure
!   - b_Butcher(1:m_etapes) : Vecteur b de Butcher explicite
!   - b_tilde_Butcher(1:m_etapes) : Vecteur b_tilde de Butcher implicite
!   - c_Butcher(1:m_etapes) : Vecteur c de Butcher explicite
!   - c_tilde_Butcher(1:m_etapes) : Vecteur c_tilde de Butcher implicite
! ==================================================================================================== !

! Entree/Sortie de la subroutine
            type(Schema_Temporel), intent(inout)                        :: RK

            RK%m_etapes = 2

            allocate(RK%A_Butcher(1:2, 1:2), RK%b_Butcher(1:2), RK%c_Butcher(1:2))
            allocate(RK%A_tilde_Butcher(1:2, 1:2), RK%b_tilde_Butcher(1:2), RK%c_tilde_Butcher(1:2))

! Allocation des tableaux explicites
            RK%A_Butcher = 0._pr
            RK%A_Butcher(2, 1) = 1._pr

            RK%b_Butcher = [1._pr/2._pr, 1._pr/2._pr]

            RK%c_Butcher = [0._pr, 1._pr]

! Allocation des tableaux implicites
            RK%A_tilde_Butcher = 0._pr
            RK%b_tilde_Butcher = 0._pr
            RK%c_tilde_Butcher = 0._pr

        end subroutine make_Butcher_SSP_2_2

        subroutine make_Butcher_SSP_3_3(RK)
! ==================================================================================================== !
! Subroutine initialisant les coefficients du tableau de Butcher pour le schema SSP(3,3)
            
! Entree/Sortie de la subroutine :
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
!   - A_Butcher(1:m_etapes, 1:m_etapes) : Matrice A de Butcher explicite triangulaire 
! strictement inferieure
!   - A_tilde_Butcher(1:m_etapes, 1:m_etapes) : Matrice A_tilde de Butcher implicite triangulaire
! inferieure
!   - b_Butcher(1:m_etapes) : Vecteur b de Butcher explicite
!   - b_tilde_Butcher(1:m_etapes) : Vecteur b_tilde de Butcher implicite
!   - c_Butcher(1:m_etapes) : Vecteur c de Butcher explicite
!   - c_tilde_Butcher(1:m_etapes) : Vecteur c_tilde de Butcher implicite
! ==================================================================================================== !

! Entree/Sortie de la subroutine
            type(Schema_Temporel), intent(inout)                        :: RK

            RK%m_etapes = 3

            allocate(RK%A_Butcher(1:3, 1:3), RK%b_Butcher(1:3), RK%c_Butcher(1:3))
            allocate(RK%A_tilde_Butcher(1:3, 1:3), RK%b_tilde_Butcher(1:3), RK%c_tilde_Butcher(1:3))

! Allocation des tableaux explicites
            RK%A_Butcher = 0._pr
            RK%A_Butcher(2, 1) = 1._pr
            RK%A_Butcher(3, 1) = 1._pr/4._pr ; RK%A_Butcher(3, 2) = 1._pr/4._pr

            RK%b_Butcher = [1._pr/6._pr, 1._pr/6._pr, 2._pr/3._pr]

            RK%c_Butcher = [0._pr, 1._pr, 1._pr/2._pr]

! Allocation des tableaux implicites
            RK%A_tilde_Butcher = 0._pr
            RK%b_tilde_Butcher = 0._pr
            RK%c_tilde_Butcher = 0._pr

        end subroutine make_Butcher_SSP_3_3

        subroutine make_Butcher_IMEX_4_3_3(RK)
! ==================================================================================================== !
! Subroutine initialisant les coefficients du tableau de Butcher pour le schema IMEX-SSP3(4,3,3)
            
! Entree/Sortie de la subroutine :
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
!   - A_Butcher(1:m_etapes, 1:m_etapes) : Matrice A de Butcher explicite triangulaire 
! strictement inferieure
!   - A_tilde_Butcher(1:m_etapes, 1:m_etapes) : Matrice A_tilde de Butcher implicite triangulaire
! inferieure
!   - b_Butcher(1:m_etapes) : Vecteur b de Butcher explicite
!   - b_tilde_Butcher(1:m_etapes) : Vecteur b_tilde de Butcher implicite
!   - c_Butcher(1:m_etapes) : Vecteur c de Butcher explicite
!   - c_tilde_Butcher(1:m_etapes) : Vecteur c_tilde de Butcher implicite
! ==================================================================================================== !

! Entree/Sortie de la subroutine
            type(Schema_Temporel), intent(inout)                        :: RK

! Variables locales
            real(kind = pr)                                             :: alpha_I, beta_I, eta_I

            RK%m_etapes = 4

            allocate(RK%A_Butcher(1:4, 1:4), RK%b_Butcher(1:4), RK%c_Butcher(1:4))
            allocate(RK%A_tilde_Butcher(1:4, 1:4), RK%b_tilde_Butcher(1:4), RK%c_tilde_Butcher(1:4))

! Allocation des tableaux explicites
            RK%A_Butcher = 0._pr
            RK%A_Butcher(3, 2) = 1._pr
            RK%A_Butcher(4, 2) = 1._pr/4._pr ; RK%A_Butcher(4, 3) = 1._pr/4._pr

            RK%b_Butcher = [0._pr, 1._pr/6._pr, 1._pr/6._pr, 2._pr/3._pr]

            RK%c_Butcher = [0._pr, 0._pr, 1._pr, 1._pr/2._pr]

! Allocation des tableaux implicites
            alpha_I = 0.24169426078821_pr ; beta_I = 0.06042356519705_pr ; eta_I = 0.12915286960590_pr

            RK%A_tilde_Butcher = 0._pr
            RK%A_tilde_Butcher(1, 1) = alpha_I
            RK%A_tilde_Butcher(2, 1) = -alpha_I ; RK%A_tilde_Butcher(2, 2) = alpha_I
            RK%A_tilde_Butcher(3, 2) = 1._pr - alpha_I ; RK%A_tilde_Butcher(3, 3) = alpha_I
            RK%A_tilde_Butcher(4, :) = [beta_i, eta_I, 1._pr/2._pr - beta_I - eta_I - alpha_I, alpha_I]

            RK%b_tilde_Butcher = [0._pr, 1._pr/6._pr, 1._pr/6._pr, 2._pr/3._pr]

            RK%c_tilde_Butcher = [alpha_I, 0._pr, 1._pr, 1._pr/2._pr]

        end subroutine make_Butcher_IMEX_4_3_3

        subroutine make_Butcher_SSP_5_4(RK)
! ==================================================================================================== !
! Subroutine initialisant les coefficients du tableau de Butcher pour le schema SSP(3,3)
            
! Entree/Sortie de la subroutine :
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
!   - A_Butcher(1:m_etapes, 1:m_etapes) : Matrice A de Butcher explicite triangulaire 
! strictement inferieure
!   - A_tilde_Butcher(1:m_etapes, 1:m_etapes) : Matrice A_tilde de Butcher implicite triangulaire
! inferieure
!   - b_Butcher(1:m_etapes) : Vecteur b de Butcher explicite
!   - b_tilde_Butcher(1:m_etapes) : Vecteur b_tilde de Butcher implicite
!   - c_Butcher(1:m_etapes) : Vecteur c de Butcher explicite
!   - c_tilde_Butcher(1:m_etapes) : Vecteur c_tilde de Butcher implicite
! ==================================================================================================== !

            type(Schema_Temporel), intent(inout)                        :: RK

            RK%m_etapes = 5

            allocate(RK%A_Butcher(1:5, 1:5), RK%b_Butcher(1:5), RK%c_Butcher(1:5))
            allocate(RK%A_tilde_Butcher(1:5, 1:5), RK%b_tilde_Butcher(1:5), RK%c_tilde_Butcher(1:5))

! Allocation des tableaux explicites
            RK%A_Butcher = 0._pr
            RK%A_Butcher(2, 1) = 0.39175222700392_pr
            RK%A_Butcher(3, 1) = 0.21766909633821_pr ; RK%A_Butcher(3, 2) = 0.36841059262959_pr
            RK%A_Butcher(4, 1) = 0.08269208670950_pr ; RK%A_Butcher(4, 2) = 0.13995850206999_pr
            RK%A_Butcher(4, 3) = 0.25189177424738_pr
            RK%A_Butcher(5, 1) = 0.06796628370320_pr ; RK%A_Butcher(5, 2) = 0.11503469844438_pr
            RK%A_Butcher(5, 3) = 0.20703489864929_pr ; RK%A_Butcher(5, 4) = 0.54497475021237_pr

            RK%b_Butcher = [0.14681187618661_pr, 0.24848290924556_pr, 0.10425883036650_pr, 0.27443890091960_pr, 0.22600748319395_pr]

            RK%c_Butcher = [0._pr, 0.39175222700392_pr, 0.58607968896779_pr, 0.47454236302687_pr, 0.93501063100924_pr]

! Allocation des tableaux implicites
            RK%A_tilde_Butcher = 0._pr
            RK%b_tilde_Butcher = 0._pr
            RK%c_tilde_Butcher = 0._pr

        end subroutine make_Butcher_SSP_5_4

        subroutine SSP(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t_n)
! ==================================================================================================== !
! Subroutine effectuant un pas de temps des schemas SSP

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi

!   - t_n : Temps apres n pas de temps

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
! ==================================================================================================== !

! Entrees de la subroutine
            real(kind = pr), intent(in)                                 :: t_n

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Schema_Temporel), intent(in)                           :: RK

! Entrees/Sorties de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: i, s, N_ddl_i
            real(kind = pr)                                             :: t_s

            !$omp parallel private(i, s, t_s, N_ddl_i) shared(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t_n)

! Stockage de la bathymetrie a l'instant t_n
            !$omp do
            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%N_Omega_i(i)%B_i_n = Calc_Exner%B_i(0:N_ddl_i-1, i)
            end do
            !$omp end do

            do s = 1, RK%m_etapes
                t_s = t_n + RK%c_Butcher(s)*Calc_Exner%dt

! Calcul et stockage des operateurs de la methode des lignes pour l'etape s
                !$omp do
                do i = 1, Mesh%N_elements
                    call make_L_i_E(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s)
                    call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s)

                    Calc_Exner%N_Omega_i(i)%L_i_s_E(:, s) = Calc_Exner%N_Omega_i(i)%L_i_E + Calc_Exner%N_Omega_i(i)%L_i_I
                end do
                !$omp end do

! Calcul de l'etat intermediaire a l'etape (s+1)
                if (s < RK%m_etapes) then
                    !$omp do
                    do i = 1, Mesh%N_elements
                        N_ddl_i = Base_ref%N_ddl_i(i)

                        Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%N_Omega_i(i)%B_i_n(0:N_ddl_i-1)                             &
                    &           + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_E(0:N_ddl_i-1, 1:s), RK%A_Butcher(s+1, 1:s))
                    end do
                    !$omp end do

                    ! call Krivodonova_2007(Mesh, Quad_ref, Base_ref, Calc_Exner, t_s)
                end if
            end do

! Calcul de la solution a l'instant t_np1
            !$omp do
            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%N_Omega_i(i)%B_i_n(0:N_ddl_i-1)                                     &
            &    + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_E(0:N_ddl_i-1, 1:RK%m_etapes), RK%b_Butcher(1:RK%m_etapes))
            end do
            !$omp end do

            !$omp end parallel

        end subroutine SSP

        subroutine IMEX(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t_n)
! ==================================================================================================== !
! Subroutine effectuant un pas de temps des schemas IMEX-SSP

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - t_n : Temps apres n pas de temps

! Entrees/Sorties de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - RK : Type derive contenant les tableaux de Butcher du schema RK choisi
! ==================================================================================================== !

! Entrees de la subroutine
            real(kind = pr), intent(in)                                 :: t_n

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entrees/Sorties de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            type(Schema_Temporel), intent(inout)                        :: RK

! Variables locales
            integer                                                     :: i, s, N_ddl_i
            real(kind = pr)                                             :: t_s_E, t_s_I

            !$omp parallel private(i, s, t_s_E, t_s_I, N_ddl_i) shared(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t_n)

! Stockage de la bathymetrie a l'instant t_n
            !$omp do
            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%N_Omega_i(i)%B_i_n = Calc_Exner%B_i(0:N_ddl_i-1, i)
            end do
            !$omp end do

            do s = 1, RK%m_etapes
                t_s_E = t_n + RK%c_Butcher(s)*Calc_Exner%dt
                t_s_I = t_n + RK%c_tilde_Butcher(s)*Calc_Exner%dt

                !$omp do
                do i = 1, Mesh%N_elements
                    N_ddl_i = Base_ref%N_ddl_i(i)
                    
                    BLOCK
                        real(pr), dimension(0:N_ddl_i-1)                :: R_i

! Construction du second membre de l'algorithme IMEX
                        if (s == 1) then
                            R_i = Calc_Exner%N_Omega_i(i)%B_i_n(0:N_ddl_i-1)
                        else
                            R_i = Calc_Exner%N_Omega_i(i)%B_i_n(0:N_ddl_i-1)                                                    &
                        &       + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_E(:, 1:s-1), RK%A_Butcher(s, 1:s-1))       &
                        &       + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_I(:, 1:s-1), RK%A_tilde_Butcher(s, 1:s-1))
                        end if

                        call Newton(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, s, i, t_s_I, R_i,                                 &
                    &                                                                     Calc_Exner%dt*RK%A_tilde_Butcher(s, s))
                    END BLOCK
                end do
                !$omp end do

! Calcul et stockage des operateurs de la methode des lignes pour l'etape s
                !$omp do
                do i = 1, Mesh%N_elements
                    call make_L_i_E(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_E)
                    call make_L_i_I(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t_s_I)

                    Calc_Exner%N_Omega_i(i)%L_i_s_E(:, s) = Calc_Exner%N_Omega_i(i)%L_i_E
                    Calc_Exner%N_Omega_i(i)%L_i_s_I(:, s) = Calc_Exner%N_Omega_i(i)%L_i_I
                end do
                !$omp end do
            end do

! Calcul de la solution a l'instant t_np1
            !$omp do
            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%N_Omega_i(i)%B_i_n(0:N_ddl_i-1)                                 &
            & + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_E(0:N_ddl_i-1, 1:RK%m_etapes),                           &
            &                                                                                  RK%b_Butcher(1:RK%m_etapes)) &
            & + Calc_Exner%dt*MATMUL(Calc_Exner%N_Omega_i(i)%L_i_s_I(0:N_ddl_i-1, 1:RK%m_etapes),                           &
            &                                                                            RK%b_tilde_Butcher(1:RK%m_etapes))
            end do
            !$omp end do

            !$omp end parallel

        end subroutine IMEX

end module mod_resolution