program exner

    use mod_parametres
    use mod_fonctions
    use mod_legendre
    use mod_interpolation
    use mod_resolution
    use mod_resolution_ODE

    implicit none

    character(len = 50)                             :: t_ch
    integer                                         :: i, j, k, n_dt
    real(kind = pr)                                 :: dx, t, dt, t_sortie, dt_sortie, z, xi, b_xi,         &
&            erreur_L2, erreur_Linf, xi_gauss, b_h, b_ex, t_debut, t_fin, b_min, b_max
    real(kind = pr), dimension(2)                   :: b_L_i, b_R_i
    real(kind = pr), dimension(0:p)                 :: B_i_loc, Beta_i
    real(kind = pr), dimension(1:N_G)               :: test
    real(kind = pr), dimension(:), allocatable      :: Li_racines, wi_legendre, omega_i
    real(kind = pr), dimension(:, :), allocatable   :: Phi_i, Phi_xi_i, x_gauss, B_i, h_ini, u_ini, b_ini,  &
&            h_tilde, u_tilde, b_tilde

! Selection du cas test dans le terminal
    print *, "================================================================"
    print *, "Résolution explicite RK3"
    print *, "================================================================"
    print *, "Veuillez choisir le cas test :"
    print *, "  1) Solution manufacturée"
    print *, "  2) Solution manufacturée stationnaire"
    print *, "  3) Cas test bathymétrie bossue"
    print *, "  4) Solution manufacturée - Loi de Grass"
    print *, "  5) Solution manufacturée - Loi de Meyer-Peter et Müller"
    print *, "  6) Solution manufacturée - Loi de Van Rijn"
    print *, "  7) Solution analytique de Long et al. (2008)"
    print *, "  8) Solution analytique de Berthon et al. (2012) utilisant la loi de Grass"
    print *, "  9) Solution analytique de Berthon et al. (2012) utilisant la loi de MPM"
    print *, "  10) Solution analytique de Berthon et al. (2012) utilisant la loi de VR"
    print *, "  11) Solution semi-analytique de Fernandez-Nieto et al. (2025)"
    read(*, *) cas
    print *, "Veuillez choisir le schéma :"
    print *, "  1) RK3"
    print *, "  2) SSPRK(3, 3)"
    print *, "  3) SSPRK(4, 3)"
    print *, "  4) RK4"
    print *, "  5) SSPRK(5, 4)"
    read(*, *) schema

    select case (schema)

    case (3)
        alpha_RK = 1.3_pr
    
    case (5)
        alpha_RK = 1.35_pr

    case default
        alpha_RK = 1._pr

    end select

    allocate(Li_racines(1:N_G), wi_legendre(1:N_G))
    allocate(omega_i(1:N_x))
    allocate(Phi_i(0:p, 1:N_G), Phi_xi_i(0:p, 1:N_G))
    allocate(x_gauss(1:N_G, 1:N_x))
    allocate(B_i(0:p, 1:N_x))
    allocate(h_ini(0:p, 1:N_x), u_ini(0:p, 1:N_x), b_ini(0:p, 1:N_x))
    allocate(h_tilde(0:p, 1:N_x), u_tilde(0:p, 1:N_x), b_tilde(0:p, 1:N_x))

    t = 0._pr

    select case (cas)

    case (7, 11)
        t_sortie = 200._pr ; dt_sortie = 200._pr

    case (8, 9, 10)
        t_sortie = 0.7_pr ; dt_sortie = 0.7_pr

    case default
        t_sortie = 0.1_pr ; dt_sortie = 0.1_pr

    end select
    
    dx = L/N_x
    n_dt = 0
    b_min = 1.e10_pr ; b_max = -1.e10_pr

    call cpu_time(t_debut)

    call racines_legendre(Li_racines)
    call poids_legendre(Li_racines, wi_legendre)
    call points_gauss(dx, Li_racines, x_gauss)

    call make_Phi_i(Li_racines, Phi_i)
    call make_Phi_xi_i(Li_racines, Phi_xi_i)

    call make_B_i(wi_legendre, Phi_i, x_gauss, B_i)

    call reconstruction_locale(1, B_i, b_L_i, b_R_i)

    B_i_loc = B_i(:, 1)

    call make_Beta_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i_loc, b_L_i, b_R_i, Beta_i)

    test = 0.0_pr

! Boucle sur les points de Gauss (on veut la valeur à chaque point j)
    do j = 1, N_G
        ! Somme sur toutes les fonctions de base (de 0 à p)
        do k = 0, p
            test(j) = test(j) + Beta_i(k) * Phi_i(k, j)
        end do
        
        ! Ici, test(j) contient la vraie valeur physique de nabla_b au point j
        print *, "Au point de Gauss ", j, " : nabla_b = ", test(j)
    end do

    call make_omega_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i, omega_i)
    print *, "================================================================"
    print *, "Début : Masse totale = ", SUM(omega_i)
    print *, "================================================================"

    write(t_ch, "(F8.3)") t
    open(10, file="./SORTIE/b_"//trim(adjustl(t_ch))//".dat")
    do i = 1, N_x
        do j = 0, p
            if (0 < p) then
                z = -1._pr + j*2._pr/p
            else
                z = 0._pr
            end if
            xi = (i - 0.5_pr)*dx + z*dx/2

            b_xi = 0._pr
            do k = 0, p
                b_xi = b_xi + B_i(k, i)*Li(k, z)
            end do

            if (b_xi < b_min) b_min = b_xi
            if (b_max < b_xi) b_max = b_xi

            write(10, *) xi, b_xi
        end do
    end do
    close(10)

    print *, "b_min(0) = ", b_min
    print *, "b_max(0) = ", b_max
    print *, "================================================================"

    select case (cas)

    case (11)
        call make_h_u_b_ini(wi_legendre, x_gauss, Phi_i, h_ini, u_ini, b_ini)

        h_tilde = 0._pr ; u_tilde = 0._pr ; b_tilde = 0._pr
        h_tilde(0, :) = h_0_FN ; u_tilde(0, :) = 1._pr ; b_tilde(0, :) = 1._pr

    end select

    do while (t < t_final)
        call find_dt(dx, t, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i, dt)

        if (t_sortie < t + dt) then
            dt = t_sortie - t
        end if

        if (t_final < t + dt) then
            dt = t_final - t
        end if

        select case (cas)

        case (11)
            call SSPRK_3_3_ODE(dx, t, dt, h_ini, u_ini, b_ini, h_tilde, u_tilde, b_tilde, wi_legendre,  &
    &                                                                     Phi_i, Phi_xi_i, x_gauss, B_i)

            call apply_ABC(wi_legendre, x_gauss, Phi_i, h_ini, b_ini, h_tilde, u_tilde, B_i)

        case default
            ! call RK1(dx, t, dt, wi_legendre, x_gauss, Phi_i, Phi_xi_i, B_i)

            select case (schema)
            
            case (1)
                call RK3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)

            case (2)
                call SSPRK_3_3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)

            case (3)
                call SSPRK_4_3(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)

            case (4)
                call SSPRK_5_4(dx, t, dt, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i)

            end select

        end select

        t = t + dt
        n_dt = n_dt + 1
        
        if (ABS(t - t_sortie) < eps_var) then
            b_min = 1.e10_pr ; b_max = -1.e10_pr
            print *, "t = ", t
            write(t_ch, "(F8.3)") t
            open(10, file="./SORTIE/b_"//trim(adjustl(t_ch))//".dat")
            do i = 1, N_x
                do j = 0, p
                    if (0 < p) then
                        z = -1._pr + j*2._pr/p
                    else
                        z = 0._pr
                    end if
                    xi = (i - 0.5_pr)*dx + z*dx/2

                    b_xi = 0._pr ; b_ex = 0._pr
                    do k = 0, p
                        b_xi = b_xi + B_i(k, i)*Li(k, z)

                        if (cas == 11) then
                            b_ex = b_ex + b_tilde(k, i)*Li(k, z)
                        end if
                    end do

                    if (cas == 11) then
                        b_ex = b(xi, 0._pr)*b_ex
                    else
                        b_ex = b(xi, t)
                    end if

                    if (b_xi < b_min) b_min = b_xi
                    if (b_max < b_xi) b_max = b_xi

                    if (cas == 11) then
                        write(10, *) xi, b_xi, b_ex
                    else
                        write(10, *) xi, b_xi
                    end if
                end do
            end do
            close(10)

            erreur_L2 = 0._pr
            erreur_Linf = 0._pr
            do i = 1, N_x
                do j = 1, N_G
                    xi_gauss = x_gauss(j, i)
                    b_h = 0._pr
                    do k = 0, p
                        b_h = b_h + B_i(k, i)*Phi_i(k, j)
                    end do

                    select case (cas)

                    case (11)
                        b_ex = 0._pr
                        do k = 0, p
                            b_ex = b_ex + b_tilde(k, i)*Phi_i(k, j)
                        end do
                        b_ex = b(xi_gauss, 0._pr)*b_ex
                    
                    case default
                        b_ex = b(xi_gauss, t)

                    end select

                    erreur_L2 = erreur_L2 + wi_legendre(j)*(b_h - b_ex)**2*dx/2
                    
                    if(erreur_Linf < ABS(b_h - b_ex)) then
                        erreur_Linf = ABS(b_h - b_ex)
                    end if
                end do
            end do

            print *, "Erreur L2 = ", SQRT(erreur_L2)
            print *, "Erreur Linf = ", erreur_Linf
            print *, "b_min(t_sortie) = ", b_min
            print *, "b_max(t_sortie) = ", b_max
            print *, "================================================================"

            t_sortie = t_sortie + dt_sortie
            
        end if
    end do

    call make_omega_i(dx, wi_legendre, Phi_i, Phi_xi_i, B_i, omega_i)
    print *, "Fin : Masse totale = ", SUM(omega_i)
    print *, "================================================================"

    call cpu_time(t_fin)

    print *, "Temps de simulation (s) : ", t_fin - t_debut
    print *, "Nombre de pas de temps : ", n_dt
    print *, "Pas de temps moyen (s) : ", t_final/n_dt
    print *, "================================================================"

    deallocate(h_ini, u_ini, b_ini, h_tilde, u_tilde, b_tilde)
    deallocate(Li_racines, wi_legendre, Phi_i, Phi_xi_i, x_gauss, B_i, omega_i)
    
end program exner
