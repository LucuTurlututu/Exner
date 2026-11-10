module mod_limiteur

    use mod_parametres
    use mod_qte_physiques

    implicit none

    contains

        pure function minmod(a, b, c) result(minmod_abc)
! ============================================================================ !
! Fonction implementant la fonction minmod(a, b, c)

! Entrees de la fonction :
!   - a, b, c : Arguments de la fonction minmod

! Sortie de la fonction :
!   - minmod_abc : Fonction minmod(a, b, c)
! ============================================================================ !

! Entrees de la fonction
            real(kind = pr), intent(in)             :: a, b, c

! Sortie de la fonction
            real(kind = pr)                         :: minmod_abc

            if ((a > 0._pr .AND. b > 0._pr .AND. c > 0._pr) .OR.                                        &
        &                                              (a < 0._pr .AND. b < 0._pr .AND. c < 0._pr)) then
                minmod_abc = SIGN(1._pr, a) * MIN(ABS(a), ABS(b), ABS(c))
            else
                minmod_abc = 0._pr
            end if

        end function minmod

        subroutine krivodonova_2007(t, B_i)
! ============================================================================ !
! Subroutine implementant le limiteur de Krivodonova (2007)

! Entree de la subroutine :
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - B_i(0:p, 1:N_x) : Matrice B_i
! ============================================================================ !

! Entree de la subroutine
            real(kind = pr), intent(in)             :: t

! Entree/Sortie de la subroutine
            real(kind = pr), dimension(0:p, 1:N_x), intent(inout)       :: B_i

! Variables locales
            integer                                 :: i, k
            logical                                 :: limiteur
            real(kind = pr)                         :: c_k, c_k_tilde, D_moins, D_plus, B_i_L, B_i_R, alpha_i

            do i = 1, N_x
                limiteur = .TRUE.
                k = p
                alpha_i = 1._pr
                do while (1 <= k .AND. limiteur)
                    if (i == 1) then
! Gestion des conditions aux limites sur la 1e maille
                        select case (cas)

                        case (2, 8, 9, 10)
! On a des conditions aux limites de Dirichlet
                            if (k-1 == 0) then
! L'ordre 0 de B_i contient la moyenne, i.e. la condition de Dirichlet
                                B_i_L = b(0._pr, t)
                            else
! L'ordre 1 <= k de B_i est nul pour ne pas creer de pente artificielle
                                B_i_L = B_i(k-1, 1)
                            end if

                        case default
! On a des conditions aux limites periodiques
                            B_i_L = B_i(k-1, N_x)
                        
                        end select
                    else
                        B_i_L = B_i(k-1, i-1)
                    end if

                    if (i == N_x) then
! Gestion des conditions aux limites sur la derniere maille
                        select case (cas)

                        case (2, 8, 9, 10)
! On a des conditions aux limites de Dirichlet
                            if (k-1 == 0) then
! L'ordre 0 de B_i contient la moyenne, i.e. la condition de Dirichlet
                                B_i_R = b(L, t)
                            else
! L'ordre 1 <= k de B_i est nul pour ne pas creer de pente artificielle
                                B_i_R = B_i(k-1, N_x)
                            end if

                        case default
! On a des conditions aux limites periodiques
                            B_i_R = B_i(k-1, 1)
                        
                        end select
                    else
                        B_i_R = B_i(k-1, i+1)
                    end if

! Application du limiteur de Krivodonova
                    c_k = B_i(k, i)
                    D_moins = B_i(k-1, i) - B_i_L
                    D_plus = B_i_R - B_i(k-1, i)

                    c_k_tilde = minmod(c_k, alpha_i*D_moins, alpha_i*D_plus)

                    if (ABS(c_k_tilde - c_k) < eps_var) then
                        limiteur = .FALSE.
                    else
                        B_i(k, i) = c_k_tilde
                        k = k-1
                    end if
                end do
            end do

        end subroutine krivodonova_2007

end module mod_limiteur