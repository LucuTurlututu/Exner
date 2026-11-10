program exner_2d

    use mod_parametres
    use mod_types_derives
    use mod_maillage
    use mod_legendre
    use mod_base
    use mod_interpolation
    use mod_resolution
    use omp_lib
    
    implicit none

    character(len = 255)                            :: t_ch

    integer                                         :: i, j, it

    real(kind = pr)                                 :: t, t_sortie, dt_sortie, b_min, b_max, b_j, b_ex, t_debut, t_fin
    real(kind = pr)                                 :: xi_j, eta_j, erreur_L2, erreur_Linf
    real(kind = pr), dimension(2)                   :: Coord_phy
    real(kind = pr), dimension(:), allocatable      :: racines_Legendre, w_j_Legendre

    type(Maillage)                                  :: Mesh
    type(Quadrature_hat_Omega)                      :: Quad_ref
    type(Base_hat_Omega)                            :: Base_ref
    type(Calcul)                                    :: Calc_Exner
    type(Schema_Temporel)                           :: RK

! Nombre de threads d'OpenMP
    call omp_set_num_threads(4)

    print *, "================================================================"
    print *, "Veuillez choisir le cas test :"
    print *, "  1) Solution manufacturée 1D : Loi de Grass"
    print *, "  2) Solution manufacturée 1D : Loi de Meyer-Peter et Müller"
    print *, "  3) Solution manufacturée 1D : Loi de Van Rijn"
    print *, "  4) Solution manufacturée 2D : Loi de Grass"
    print *, "  5) Solution manufacturée 2D : Loi de Meyer-Peter et Müller"
    print *, "  6) Solution manufacturée 2D : Loi de Van Rijn"
    print *, "  7) Berthon et al. (2012) 1D : Loi de Grass"
    print *, "  8) Berthon et al. (2012) 1D : Loi de Meyer-Peter et Müller"
    print *, "  9) Berthon et al. (2012) 1D : Loi de Van Rijn"
    print *, "================================================================"
    read *, cas

    t = 0._pr
    it = 0
    select case (cas)

    case (1, 2, 3)
        t_sortie = 0.1_pr ; dt_sortie = 0.1_pr

    case (7, 8, 9)
        t_sortie = 0.7_pr ; dt_sortie = 0.7_pr

    case default
        t_sortie = 0.05_pr ; dt_sortie = 0.05_pr

    end select

    print *, "================================================================"
    print *, "Veuillez choisir le schéma temporel :"
    print *, "  1) SSP(2,2)"
    print *, "  2) SSP(3,3)"
    print *, "  3) IMEX-SSP3(4,3,3)"
    print *, "  4) SSP(5,4)"
    print *, "================================================================"
    read *, schema
    print *, "================================================================"

    t_debut = omp_get_wtime()

    select case (schema)

    case (1)
        call make_Butcher_SSP_2_2(RK)

    case (2)
        call make_Butcher_SSP_3_3(RK)

    case (3)
        call make_Butcher_IMEX_4_3_3(RK)

    case (4)
        call make_Butcher_SSP_5_4(RK)

    case default
        print *, "Schéma non défini"
        stop

    end select

    call make_Mesh(fichier, Mesh)

    call racines_Legendre_1D(racines_Legendre)
    call poids_Legendre_1D(racines_Legendre, w_j_Legendre)

! Creation des tableaux de quadrature sur l'element de reference
    call make_Quad_ref(racines_Legendre, w_j_Legendre, Quad_ref)
    call make_Base_ref(Mesh, Quad_ref, Base_ref)

    deallocate(racines_Legendre, w_j_Legendre)

! Allocation de tous les tableaux
    allocate(Calc_Exner%N_Omega_i(Mesh%N_elements))

    !$omp parallel do private(i) shared(Calc_Exner, Mesh, Base_ref, Quad_ref, RK)

    do i = 1, Mesh%N_elements
        if (Mesh%N_sommets_Omega_i(i) == 3) then
            allocate(Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri(0:Base_ref%N_ddl_i_Tri-1, Quad_ref%N_G_Tri, 2))
            Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri(:, :, :) = 0._pr

        else if (Mesh%N_sommets_Omega_i(i) == 4) then
            allocate(Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad(0:Base_ref%N_ddl_i_Quad-1, Quad_ref%N_G_Quad, 2))
            Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad(:, :, :) = 0._pr

        end if

        allocate(Calc_Exner%N_Omega_i(i)%M_i_std(0:Base_ref%N_ddl_i(i)-1, 0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%M_i_std = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%K_i_x_std(0:Base_ref%N_ddl_i(i)-1, 0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%K_i_x_std = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%K_i_y_std(0:Base_ref%N_ddl_i(i)-1, 0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%K_i_y_std = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%Beta_i(0:Base_ref%N_ddl_i(i)-1, 1:2))
        Calc_Exner%N_Omega_i(i)%Beta_i = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%M_i(0:Base_ref%N_ddl_i(i)-1, 0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%M_i = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%I_i_EF(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%I_i_EF = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%I_i_VF(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%I_i_VF = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%I_i_S(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%I_i_S = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%I_i_R_tilde(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%I_i_R_tilde = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%L_i_E(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%L_i_E = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%L_i_I(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%L_i_I = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%B_i_n(0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%B_i_n = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%L_i_s_E(0:Base_ref%N_ddl_i(i)-1, 1:RK%m_etapes))
        Calc_Exner%N_Omega_i(i)%L_i_s_E = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%L_i_s_I(0:Base_ref%N_ddl_i(i)-1, 1:RK%m_etapes))
        Calc_Exner%N_Omega_i(i)%L_i_s_I = 0._pr

        allocate(Calc_Exner%N_Omega_i(i)%ipiv(1:Base_ref%N_ddl_i(i)))
        Calc_Exner%N_Omega_i(i)%ipiv = 0

        allocate(Calc_Exner%N_Omega_i(i)%J_i_q(0:Base_ref%N_ddl_i(i)-1, 0:Base_ref%N_ddl_i(i)-1))
        Calc_Exner%N_Omega_i(i)%J_i_q = 0._pr
    end do

    !$omp end parallel do

! Construction de la matrice de masse standard
    call make_M_i_std(Mesh, Quad_ref, Base_ref, Calc_Exner)

! Construction de la matrice de rigidite standard
    call make_K_i_std(Mesh, Quad_ref, Base_ref, Calc_Exner)

! Construction du vecteur de bathymetrie a l'instant initial
    call make_B_i(Mesh, Quad_ref, Base_ref, Calc_Exner)

! Sortie dans un fichier
    write(t_ch, "(F8.3)") t
    open(10, file="./SORTIE/b_"//trim(adjustl(t_ch))//".dat")

    do i = 1, Mesh%N_elements
        if (Mesh%N_sommets_Omega_i(i) == 3) then
            do j = 1, Quad_ref%N_G_Tri
                xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                write(10, *) Coord_phy(1), Coord_phy(2), b_j
            end do

        else if (Mesh%N_sommets_Omega_i(i) == 4) then
            do j = 1, Quad_ref%N_G_Quad
                xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                write(10, *) Coord_phy(1), Coord_phy(2), b_j
            end do

        end if
    end do

    close(10)

! Calcul des extremas a l'instant initial
    b_min = 1.e10_pr ; b_max = -1.e10_pr

    RK%q_max_Newton = 0 ; RK%q_Newton = 0 ; RK%N_appels_Newton = 0
    RK%q_max_Picard = 0

    !$omp parallel do private(i, j, b_j) &
    !$omp shared(Mesh, Quad_ref, Base_ref, Calc_Exner) &
    !$omp reduction(min: b_min) reduction(max: b_max)

    do i = 1, Mesh%N_elements
        if (Mesh%N_sommets_Omega_i(i) == 3) then
            do j = 1, Quad_ref%N_G_Tri
                b_j = 0._pr

                b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                if (b_j < b_min) b_min = b_j
                if (b_max < b_j) b_max = b_j
            end do

        else if (Mesh%N_sommets_Omega_i(i) == 4) then
            do j = 1, Quad_ref%N_G_Quad
                b_j = 0._pr

                b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                if (b_j < b_min) b_min = b_j
                if (b_max < b_j) b_max = b_j
            end do
        end if
    end do

    !$omp end parallel do

    print *, "================================================================"
    print *, "b_min(0) = ", b_min
    print *, "b_max(0) = ", b_max
    print *, "================================================================"

! Boucle principale
    do while (t < T_final)
        call find_dt(Mesh, Quad_ref, Base_ref, Calc_Exner, t)

        if (t_sortie < t + Calc_Exner%dt) then
            Calc_Exner%dt = t_sortie - t
        end if

        if (T_final < t + Calc_Exner%dt) then
            Calc_Exner%dt = t_final - t
        end if

        select case (schema)

        case (1, 2, 4)
            call SSP(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t)

        case (3)
            call IMEX(Mesh, Quad_ref, Base_ref, Calc_Exner, RK, t)

        end select

        t = t + Calc_Exner%dt
        it = it + 1

        if (ABS(t - t_sortie) < eps_var) then
            b_min = 1.e10_pr ; b_max = -1.e10_pr
            print *, "t = ", t

! Sortie dans un fichier
            write(t_ch, "(F8.3)") t
            open(10, file="./SORTIE/b_"//trim(adjustl(t_ch))//".dat")

            do i = 1, Mesh%N_elements
                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        write(10, *) Coord_phy(1), Coord_phy(2), b_j
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        write(10, *) Coord_phy(1), Coord_phy(2), b_j
                    end do

                end if
            end do

            close(10)

! Calcul des extremas a t_sortie
            !$omp parallel do private(i, j, b_j) &
            !$omp shared(Mesh, Quad_ref, Base_ref, Calc_Exner) &
            !$omp reduction(min: b_min) reduction(max: b_max)

            do i = 1, Mesh%N_elements
                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        b_j = 0._pr

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        if (b_j < b_min) b_min = b_j
                        if (b_max < b_j) b_max = b_j
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        b_j = 0._pr

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        if (b_j < b_min) b_min = b_j
                        if (b_max < b_j) b_max = b_j
                    end do
                end if
            end do

            !$omp end parallel do

! Calcul des erreurs en norme L^2 et L^inf a t_sortie
            erreur_L2 = 0._pr ; erreur_Linf = 0._pr

            !$omp parallel do private(i, j, b_j, xi_j, eta_j, Coord_phy, b_ex) &
            !$omp shared(Mesh, Quad_ref, Base_ref, Calc_Exner, t) &
            !$omp reduction(+: erreur_L2) reduction(max: erreur_Linf)

            do i = 1, Mesh%N_elements
                BLOCK
                    type(Base_Omega_i)                              :: Base_phy

                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        b_j = 0._pr

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_ex = b(Coord_phy, t)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        erreur_L2 = erreur_L2 + Quad_ref%w_j_Ref_Tri(j)*(b_j - b_ex)**2*Base_phy%Det_J_i

                        erreur_Linf = MAX(erreur_Linf, ABS(b_j - b_ex))
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        b_j = 0._pr

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:Base_ref%N_ddl_i(i)-1, i))

                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_ex = b(Coord_phy, t)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        erreur_L2 = erreur_L2 + Quad_ref%w_j_Ref_Quad(j)*(b_j - b_ex)**2*Base_phy%Det_J_i

                        erreur_Linf = MAX(erreur_Linf, ABS(b_j - b_ex))
                    end do
                end if

                END BLOCK
            end do

            !$omp end parallel do

            print *, "Erreur L2 = ", SQRT(erreur_L2)
            print *, "Erreur Linf = ", erreur_Linf
            print *, "b_min(t_sortie) = ", b_min
            print *, "b_max(t_sortie) = ", b_max
            print *, "================================================================"

            t_sortie = t_sortie + dt_sortie

        end if
    end do

    t_fin = omp_get_wtime()

    print *, "Temps de simulation (s) : ", t_fin - t_debut
    print *, "Nombre de pas de temps : ", it
    print *, "Pas de temps moyen (s) : ", T_final/real(it, pr)
    print *, "Nombre maximal d'étapes de Newton : ", RK%q_max_Newton
    print *, "Nombre moyen d'étapes de Newton : ", RK%q_Newton/real(RK%m_etapes*it*Mesh%N_elements, pr)
    print *, "Nombre maximal d'étapes de Picard : ", RK%q_max_Picard
    print *, "Nombre de Froude maximal : ", Fr_max
    print *, "================================================================"
        
end program exner_2d
