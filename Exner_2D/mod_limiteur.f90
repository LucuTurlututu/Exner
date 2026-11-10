module mod_limiteur

    use mod_parametres
    use mod_types_derives
    use mod_cl

    implicit none

    contains

        pure function minmod(a, b, c, d, e) result(minmod_a_b_c_d_e)
! ==================================================================================================== !
! Fonction implementant la fonction minmod(a, b, c, d, e)

! Entrees de la fonction :
!   - a, b, c, d, e : Arguments de la fonction minmod

! Sortie de la fonction :
!   - minmod_a_b_c_d_e : Fonction minmod(a, b, c, d, e)
! ==================================================================================================== !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: a, b, c, d, e

! Sortie de la fonction
            real(kind = pr)                         :: minmod_a_b_c_d_e

            if (SIGN(1._pr, a) == SIGN(1._pr, b) .AND. SIGN(1._pr, a) == SIGN(1._pr, c) .AND.                           &
        &                                   SIGN(1._pr, a) == SIGN(1._pr, d) .AND. SIGN(1._pr, a) == SIGN(1._pr, e)) then
                minmod_a_b_c_d_e = SIGN(1._pr, a)*MIN(ABS(a), ABS(b), ABS(c), ABS(d), ABS(e))
            
            else
                minmod_a_b_c_d_e = 0._pr

            end if

        end function minmod

        subroutine Krivodonova_2007(Mesh, Quad_ref, Base_ref, Calc_Exner, t)
! ==================================================================================================== !
! Subroutine implementant le limiteur de Krivodonova (2007) pour les mailles rectangulaires

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
            logical                                                     :: limiteur, modif_1, modif_2

            integer                                                     :: i, k, N_ddl_i
            integer                                                     :: index_i, index_j, index_k_m, M, L

            real(kind = pr)                                             :: c_k_m, c_k_m_tilde
            real(kind = pr)                                             :: D_i_moins, D_i_plus, D_j_moins, D_j_plus
            real(kind = pr), dimension(0:p)                             :: alpha

            alpha(0) = 1.0_pr
            do k = 1, p
                alpha(k) = SQRT(real(2*k-1, pr)/real(2*k+1, pr))
            end do

            do i = 1, Mesh%N_elements

                N_ddl_i = Base_ref%N_ddl_i(i)

                BLOCK
                    real(kind = pr), dimension(0:N_ddl_i-1)             :: B_Sud, B_Est, B_Nord, B_Ouest

! Face Sud
                    if (Mesh%N_Omega_r(1, i) /= 0) then
                        B_Sud = Calc_Exner%B_i(0:N_ddl_i-1, Mesh%N_Omega_r(1, i))
                    
                    else
                        select case (Mesh%CL_Arete_Bord(Mesh%Aretes_Omega_i(1, i)))
                        case ("D"); B_Sud = Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, 1, t)
                        case ("N"); B_Sud = Neumann(Base_ref, Calc_Exner, i)
                        case ("P"); B_Sud = Periodique(Mesh, Base_ref, Calc_Exner, i, 1)
                        end select

                    end if

! Face Est
                    if (Mesh%N_Omega_r(2, i) /= 0) then
                        B_Est = Calc_Exner%B_i(0:N_ddl_i-1, Mesh%N_Omega_r(2, i))
                    
                    else
                        select case (Mesh%CL_Arete_Bord(Mesh%Aretes_Omega_i(2, i)))
                        case ("D"); B_Est = Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, 2, t)
                        case ("N"); B_Est = Neumann(Base_ref, Calc_Exner, i)
                        case ("P"); B_Est = Periodique(Mesh, Base_ref, Calc_Exner, i, 2)
                        end select

                    end if

! Face Nord
                    if (Mesh%N_Omega_r(3, i) /= 0) then
                        B_Nord = Calc_Exner%B_i(0:N_ddl_i-1, Mesh%N_Omega_r(3, i))
                    
                    else
                        select case (Mesh%CL_Arete_Bord(Mesh%Aretes_Omega_i(3, i)))
                        case ("D"); B_Nord = Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, 3, t)
                        case ("N"); B_Nord = Neumann(Base_ref, Calc_Exner, i)
                        case ("P"); B_Nord = Periodique(Mesh, Base_ref, Calc_Exner, i, 3)
                        end select

                    end if

! Face Ouest
                    if (Mesh%N_Omega_r(4, i) /= 0) then
                        B_Ouest = Calc_Exner%B_i(0:N_ddl_i-1, Mesh%N_Omega_r(4, i))
                    
                    else
                        select case (Mesh%CL_Arete_Bord(Mesh%Aretes_Omega_i(4, i)))
                        case ("D"); B_Ouest = Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, 4, t)
                        case ("N"); B_Ouest = Neumann(Base_ref, Calc_Exner, i)
                        case ("P"); B_Ouest = Periodique(Mesh, Base_ref, Calc_Exner, i, 4)
                        end select

                    end if

! ==================================================================================================== !
! Cascade de Krivodonova
! ==================================================================================================== !
                    limiteur = .TRUE.
                    M = p  ! M représente le sub-script principal (p, p-1, ..., 1)

                    do while (M >= 1 .and. limiteur)

! ==================================================================================================== !
! Terme diagonal
! ==================================================================================================== !
                        index_k_m = M*(p+1) + M
                        c_k_m = Calc_Exner%B_i(index_k_m, i)

                        index_i = M*(p+1) + (M-1)
                        D_i_plus  = alpha(M)*(B_Est(index_i) - Calc_Exner%B_i(index_i, i))
                        D_i_moins = alpha(M)*(Calc_Exner%B_i(index_i, i) - B_Ouest(index_i))

                        index_j = (M-1)*(p+1) + M
                        D_j_plus  = alpha(M)*(B_Nord(index_j) - Calc_Exner%B_i(index_j, i))
                        D_j_moins = alpha(M)*(Calc_Exner%B_i(index_j, i) - B_Sud(index_j))

                        c_k_m_tilde = minmod(c_k_m, D_j_plus, D_j_moins, D_i_plus, D_i_moins)

                        if (ABS(c_k_m_tilde - c_k_m) > eps_var) then
                            Calc_Exner%B_i(index_k_m, i) = c_k_m_tilde

                        else
! Si le terme diagonal n'est pas limite on sort de la routine
                            limiteur = .FALSE.
                            exit

                        end if

! ==================================================================================================== !
! Paires extra-diagonales
! ==================================================================================================== !
                        do L = M - 1, 0, -1
                            index_k_m = M*(p+1) + L
                            c_k_m = Calc_Exner%B_i(index_k_m, i)

                            if (L >= 1) then
                                index_i = M*(p+1) + (L-1)
                                D_i_plus = alpha(L)*(B_Est(index_i) - Calc_Exner%B_i(index_i, i))
                                D_i_moins = alpha(L)*(Calc_Exner%B_i(index_i, i) - B_Ouest(index_i))

                            else
                                D_i_plus = 0.0_pr ; D_i_moins = 0.0_pr

                            end if

                            index_j = (M-1)*(p+1) + L
                            D_j_plus  = alpha(M)*(B_Nord(index_j) - Calc_Exner%B_i(index_j, i))
                            D_j_moins = alpha(M)*(Calc_Exner%B_i(index_j, i) - B_Sud(index_j))

                            c_k_m_tilde = minmod(c_k_m, D_j_plus, D_j_moins, D_i_plus, D_i_moins)

                            modif_1 = (ABS(c_k_m_tilde - c_k_m) > eps_var)
                            if (modif_1) Calc_Exner%B_i(index_k_m, i) = c_k_m_tilde

                            index_k_m = L*(p+1) + M
                            c_k_m = Calc_Exner%B_i(index_k_m, i)

                            index_i = L*(p+1) + (M-1)
                            D_i_plus  = alpha(M)*(B_Est(index_i) - Calc_Exner%B_i(index_i, i))
                            D_i_moins = alpha(M)*(Calc_Exner%B_i(index_i, i) - B_Ouest(index_i))
                            
                            if (L >= 1) then
                                index_j = (L-1)*(p+1) + M
                                D_j_plus  = alpha(L)*(B_Nord(index_j) - Calc_Exner%B_i(index_j, i))
                                D_j_moins = alpha(L)*(Calc_Exner%B_i(index_j, i) - B_Sud(index_j))

                            else
                                D_j_plus = 0.0_pr ; D_j_moins = 0.0_pr

                            end if

                            c_k_m_tilde = minmod(c_k_m, D_j_plus, D_j_moins, D_i_plus, D_i_moins)

                            modif_2 = (ABS(c_k_m_tilde - c_k_m) > eps_var)
                            if (modif_2) Calc_Exner%B_i(index_k_m, i) = c_k_m_tilde

                            if (.NOT. (modif_1 .OR. modif_2)) then
! Si les paires ne sont pas limitees on sort de la boucle
                                limiteur = .FALSE.
                                exit
                            end if

                        end do
                        M = M - 1
                    end do
                END BLOCK
            end do

        end subroutine Krivodonova_2007

end module mod_limiteur