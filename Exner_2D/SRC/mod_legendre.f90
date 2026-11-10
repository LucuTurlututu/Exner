module mod_legendre

    use mod_parametres

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_legendre :
!   - Legendre_1D
!   - Legendre_xi_1D
!   - Jacobi_1D
!   - Jacobi_xi_1D
!   - racines_Legendre_1D
!   - poids_Legendre_1D
!   - Legendre_2D
!   - Dubiner_2D
!   - Dubiner_xi_2D
!   - Dubiner_eta_2D
! ==================================================================================================== !

    contains

        pure function Legendre_1D(k, xi) result(Legendre_k_xi)
! ==================================================================================================== !
! Fonction renvoyant les polynomes de Legendre 1D par la formule de Bonnet

! Entrees de la fonction :
!   - k : Degre du polynome de Legendre 1D
!   - xi : Coordonnee x des elements de reference

! Sortie de la fonction :
!   - Legendre_k_xi : Polynome de Legendre de degre k evalue en xi
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k

            real(kind = pr), intent(in)                                 :: xi

! Sorties de la fonction
            real(kind = pr)                                             :: Legendre_k_xi

! Variables locales
            integer                                                     :: n

            real(kind = pr)                                             :: Legendre_nm1_xi, Legendre_n_xi

            Legendre_k_xi = 0._pr

            if (k == 0) then
                Legendre_k_xi = 1._pr
            else if (k == 1) then
                Legendre_k_xi = xi
            else
                Legendre_nm1_xi = 1._pr ; Legendre_n_xi = xi
                do n = 1, k-1
                    Legendre_k_xi = real(2*n+1, pr)/(n+1)*xi*Legendre_n_xi - real(n, pr)/(n+1)*Legendre_nm1_xi
                    Legendre_nm1_xi = Legendre_n_xi ; Legendre_n_xi = Legendre_k_xi
                end do
            end if

        end function Legendre_1D

        pure function Legendre_xi_1D(k, xi) result(Legendre_xi_k_xi)
! ==================================================================================================== !
! Fonction renvoyant la derivee par rapport a xi du polynome de Legendre 1D de degre k

! Entrees de la fonction :
!   - k : Degre du polynome de Legendre 1D
!   - xi : Coordonnee x des elements de reference

! Sortie de la fonction :
!   - Legendre_xi_k_xi : Derivee du polynome de Legendre de degre k
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k

            real(kind = pr), intent(in)                                 :: xi

! Sortie de la fonction
            real(kind = pr)                                             :: Legendre_xi_k_xi

! Variables locales
            integer                                                     :: n

            real(kind = pr)                                             :: Legendre_nm1_xi, Legendre_xi_nm2_xi, Legendre_xi_nm1_xi

            Legendre_xi_k_xi = 0._pr

            if (k == 0) then
                Legendre_xi_k_xi = 0._pr
            else if (k == 1) then
                Legendre_xi_k_xi = 1._pr
            else
                Legendre_xi_nm2_xi = 0._pr ; Legendre_xi_nm1_xi = 1._pr
                do n = 2, k
                    Legendre_nm1_xi = Legendre_1D(n-1, xi)
                    Legendre_xi_k_xi = (2*n - 1)*Legendre_nm1_xi + Legendre_xi_nm2_xi
                    Legendre_xi_nm2_xi = Legendre_xi_nm1_xi
                    Legendre_xi_nm1_xi = Legendre_xi_k_xi
                end do
            end if

        end function Legendre_xi_1D

        pure function Jacobi_1D(k, alpha_J, beta_J, xi) result(Jacobi_k_xi)
! ==================================================================================================== !
! Fonction renvoyant le polynome de Jacobi 1D de degre k et de parametres (alpha_J, beta_J)

! Entrees de la fonction :
!   - k : Degre du polynome de Jacobi 1D
!   - alpha_J : Premier parametre de Jacobi
!   - beta_J : Second parametre de Jacobi
!   - xi : Coordonne x des elements de reference

! Sortie de la fonction :
!   - Jacobi_k_xi : Polynome de Jacobi de degre k et de parametres (alpha_J, beta_J) evalue en xi
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, alpha_J, beta_J

            real(kind = pr), intent(in)                                 :: xi

! Sortie de la fonction
            real(kind = pr)                                             :: Jacobi_k_xi

! Variables locales
            integer                                                     :: n

            real(kind = pr)                                             :: Jacobi_nm1_xi, Jacobi_n_xi
            real(kind = pr)                                             :: c_J_1, c_J_2, c_J_3

            Jacobi_k_xi = 0._pr

            if (k == 0) then
                Jacobi_k_xi = 1._pr
            else if (k == 1) then
                Jacobi_k_xi = 1._pr/2*(alpha_J - beta_J + (alpha_J + beta_J + 2)*xi)
            else
                Jacobi_nm1_xi = 1._pr ; Jacobi_n_xi = 1._pr/2*(alpha_J - beta_J + (alpha_J + beta_J + 2)*xi)
                do n = 1, k-1
                    c_J_1 = real(2*(n+1)*(n + alpha_J + beta_J + 1)*(2*n + alpha_J + beta_J), pr)
                    c_J_2 = real(2*n + alpha_J + beta_J + 1, pr)*(real((2*n + alpha_J + beta_J)*                            &
                    &                               (2*n + alpha_J + beta_J + 2), pr)*xi + real(alpha_J**2 - beta_J**2, pr))
                    c_J_3 = real(-2*(n + alpha_J)*(n + beta_J)*(2*n + alpha_J + beta_J + 2), pr)

                    Jacobi_k_xi = c_J_2/c_J_1*Jacobi_n_xi + c_J_3/c_J_1*Jacobi_nm1_xi
                    Jacobi_nm1_xi = Jacobi_n_xi ; Jacobi_n_xi = Jacobi_k_xi
                end do
            end if

        end function Jacobi_1D

        pure function Jacobi_xi_1D(k, alpha_J, beta_J, xi) result(Jacobi_xi_k_xi)
! ==================================================================================================== !
! Fonction renvoyant la derivee du polynome de Jacobi 1D de degre k et de parametres (alpha_J, beta_J)

! Entrees de la fonction :
!   - k : Degre du polynome de Jacobi 1D
!   - alpha_J : Premier parametre de Jacobi
!   - beta_J : Second parametre de Jacobi
!   - xi : Coordonne x des elements de reference

! Sortie de la fonction :
!   - Jacobi_xi_k_xi : Derivee du polynome de Jacobi de degre k et de parametres (alpha_J, beta_J)
! evaluee en xi
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, alpha_J, beta_J

            real(kind = pr), intent(in)                                 :: xi

! Sortie de la fonction
            real(kind = pr)                                             :: Jacobi_xi_k_xi

            if (k == 0) then
                Jacobi_xi_k_xi = 0._pr
            else
                Jacobi_xi_k_xi = 1._pr/2*(alpha_J + beta_J + k + 1)*Jacobi_1D(k-1, alpha_J+1, beta_J+1, xi)
            end if

        end function Jacobi_xi_1D

        subroutine racines_Legendre_1D(racines_Legendre)
! ==================================================================================================== !
! Subroutine renvoyant un tableau contenant les racines des polynomes de Legendre 1D

! Sortie de la subroutine :
!   - racines_Legendre(1:N_G) : Tableau contenant les N_G racines du polynome de Legendre 1D de
! de degre N_G
! ==================================================================================================== !

! Sortie de la subroutine
            real(kind = pr), dimension(:), allocatable, intent(out)     :: racines_Legendre

! Variables locales
            integer                                                     :: j, m

            real(kind = pr)                                             :: xi_0, xi_q, xi_qp1, Legendre_q_xi, Legendre_xi_q_xi 

            allocate(racines_Legendre(1:N_G))

! On ne cherche que les racines positives par symetrie des polynomes de legendre
            m = INT((N_G + 1)/2._pr)

            do j = 1, m
! On initialise avec une approximation relativement bonne des racines
                xi_0 = COS(pi*(j - 0.25_pr)/(N_G + 0.5_pr))

                xi_q = 1._pr ; xi_qp1 = xi_0
                do while (ABS(xi_qp1 - xi_q) > eps_racine)
                    xi_q = xi_qp1
                    Legendre_q_xi = Legendre_1D(N_G, xi_q)
                    Legendre_xi_q_xi = Legendre_xi_1D(N_G, xi_q)

                    xi_qp1 = xi_q - Legendre_q_xi/Legendre_xi_q_xi
                end do

                racines_Legendre(j) = -xi_q
                racines_Legendre(N_G-j+1) = xi_q
            end do

        end subroutine racines_Legendre_1D

        subroutine poids_Legendre_1D(racines_Legendre, w_j_Legendre)
! ==================================================================================================== !
! Subroutine renvoyant les poids de la quadrature de Gauss-Legendre

! Entree de la subroutine :
!   - racines_Legendre(1:N_G) : Tableau contenant les N_G racines du polynome de Legendre 1D de
! de degre N_G

! Sorties de la subroutine :
!   - w_j_Legendre(1:N_G) : Tableau contenant les poids de quadrature associes aux points de
! Gauss-Legendre
! ==================================================================================================== !

! Entree de la subroutine
            real(kind = pr), dimension(:), intent(in)                   :: racines_Legendre

! Sortie de la subroutine
            real(kind = pr), dimension(:), allocatable, intent(out)     :: w_j_Legendre

! Variables locales
            integer                                                     :: j

            real(kind = pr)                                             :: xi_j, Legendre_xi_N_G_xi_j

            allocate(w_j_Legendre(1:N_G))

            do j = 1, N_G
                xi_j = racines_Legendre(j)
                Legendre_xi_N_G_xi_j = Legendre_xi_1D(N_G, xi_j)

                w_j_Legendre(j) = 2._pr/((1._pr - xi_j**2)*Legendre_xi_N_G_xi_j**2)
            end do

        end subroutine poids_Legendre_1D

        subroutine quadrature_Dunavant(N_Dunavant, racines_Dunavant, w_j_Dunavant)
! ==================================================================================================== !
! Subroutine renvoyant le nombre de points, les coordonnes des points de Dunavant et les poids
! associes

! Sorties de la subroutine :
!   - N_Dunavant : Nombre de points de quadrature
!   - racines_Dunavant(1:2, 1:N_Dunavant) : Coordonnes de la quadrature de Dunavant
!   - w_j_Dunavant(1:N_Dunavant) : Poids de quadrature de Dunavant
! ==================================================================================================== !

! Sorties de la subroutine
            integer, intent(out)                                        :: N_Dunavant
            real(kind = pr), dimension(:), allocatable, intent(out)     :: w_j_Dunavant
            real(kind = pr), dimension(:, :), allocatable, intent(out)  :: racines_Dunavant

            select case (N_G)

            case (3)
                N_Dunavant = 7
                allocate(w_j_Dunavant(1:N_Dunavant))
                allocate(racines_Dunavant(1:2, 1:N_Dunavant))

! Point central
                racines_Dunavant(:, 1) = [1._pr/3._pr, 1._pr/3._pr]
                w_j_Dunavant(1) = 0.225000000000000_pr*1._pr/2._pr

! Premier groupe symetrique
                racines_Dunavant(:, 2) = [0.059715871789770_pr, 0.470142064105115_pr]
                racines_Dunavant(:, 3) = [0.470142064105115_pr, 0.059715871789770_pr]
                racines_Dunavant(:, 4) = [0.470142064105115_pr, 0.470142064105115_pr]
                w_j_Dunavant(2:4)   = 0.132394152788506_pr*1._pr/2._pr

! Second groupe symetrique
                racines_Dunavant(:, 5) = [0.797426985353087_pr, 0.101286507323456_pr]
                racines_Dunavant(:, 6) = [0.101286507323456_pr, 0.797426985353087_pr]
                racines_Dunavant(:, 7) = [0.101286507323456_pr, 0.101286507323456_pr]
                w_j_Dunavant(5:7) = 0.125939180544827_pr*1._pr/2._pr

            case (4)
                N_Dunavant = 13
                allocate(w_j_Dunavant(1:N_Dunavant))
                allocate(racines_Dunavant(1:2, 1:N_Dunavant))

! Point central
                racines_Dunavant(:, 1) = [1._pr/3._pr, 1._pr/3._pr]
                w_j_Dunavant(1) = -0.149570044467682_pr*1._pr/2._pr

! Premier groupe
                racines_Dunavant(:, 2) = [0.479308067841920_pr, 0.260345966079040_pr]
                racines_Dunavant(:, 3) = [0.260345966079040_pr, 0.479308067841920_pr]
                racines_Dunavant(:, 4) = [0.260345966079040_pr, 0.260345966079040_pr]
                w_j_Dunavant(2:4) = 0.175615257433208_pr*1._pr/2._pr

! Deuxieme groupe
                racines_Dunavant(:, 5) = [0.869739794195568_pr, 0.065130102902216_pr]
                racines_Dunavant(:, 6) = [0.065130102902216_pr, 0.869739794195568_pr]
                racines_Dunavant(:, 7) = [0.065130102902216_pr, 0.065130102902216_pr]
                w_j_Dunavant(5:7) = 0.053347235608838_pr*1._pr/2._pr

! Troisieme groupe
                racines_Dunavant(:, 8) = [0.048690315425316_pr, 0.312865496004874_pr]
                racines_Dunavant(:, 9) = [0.048690315425316_pr, 0.638444188569810_pr]
                racines_Dunavant(:, 10) = [0.312865496004874_pr, 0.048690315425316_pr]
                racines_Dunavant(:, 11) = [0.312865496004874_pr, 0.638444188569810_pr]
                racines_Dunavant(:, 12) = [0.638444188569810_pr, 0.048690315425316_pr]
                racines_Dunavant(:, 13) = [0.638444188569810_pr, 0.312865496004874_pr]
                w_j_Dunavant(8:13) = 0.077113760890257_pr*1._pr/2._pr

            case (5)
                N_Dunavant = 19
                allocate(w_j_Dunavant(1:N_Dunavant))
                allocate(racines_Dunavant(1:2, 1:N_Dunavant))

! Point central
                racines_Dunavant(:, 1) = [1._pr/3._pr, 1._pr/3._pr]
                w_j_Dunavant(1) = 0.097135796282799_pr*1._pr/2._pr

! Premier groupe
                racines_Dunavant(:, 2) = [0.020634961602525_pr, 0.489682519198738_pr]
                racines_Dunavant(:, 3) = [0.489682519198738_pr, 0.020634961602525_pr]
                racines_Dunavant(:, 4) = [0.489682519198738_pr, 0.489682519198738_pr]
                w_j_Dunavant(2:4) = 0.031334700227139_pr*1._pr/2._pr

! Deuxieme groupe
                racines_Dunavant(:, 5) = [0.125820817014127_pr, 0.437089591492937_pr]
                racines_Dunavant(:, 6) = [0.437089591492937_pr, 0.125820817014127_pr]
                racines_Dunavant(:, 7) = [0.437089591492937_pr, 0.437089591492937_pr]
                w_j_Dunavant(5:7) = 0.077827541004774_pr*1._pr/2._pr

! Troisieme groupe
                racines_Dunavant(:, 8) = [0.623592928761935_pr, 0.188203535619033_pr]
                racines_Dunavant(:, 9) = [0.188203535619033_pr, 0.623592928761935_pr]
                racines_Dunavant(:, 10)= [0.188203535619033_pr, 0.188203535619033_pr]
                w_j_Dunavant(8:10) = 0.079647738927210_pr*1._pr/2._pr

! Quatrieme groupe
                racines_Dunavant(:, 11) = [0.910540973211095_pr, 0.044729513394453_pr]
                racines_Dunavant(:, 12) = [0.044729513394453_pr, 0.910540973211095_pr]
                racines_Dunavant(:, 13) = [0.044729513394453_pr, 0.044729513394453_pr]
                w_j_Dunavant(11:13) = 0.025577675658698_pr*1._pr/2._pr

! Cinquieme groupe
                racines_Dunavant(:, 14) = [0.036838412054736_pr, 0.221962989160766_pr] 
                racines_Dunavant(:, 15) = [0.036838412054736_pr, 0.741198598784498_pr] 
                racines_Dunavant(:, 16) = [0.221962989160766_pr, 0.036838412054736_pr] 
                racines_Dunavant(:, 17) = [0.221962989160766_pr, 0.741198598784498_pr] 
                racines_Dunavant(:, 18) = [0.741198598784498_pr, 0.036838412054736_pr] 
                racines_Dunavant(:, 19) = [0.741198598784498_pr, 0.221962989160766_pr] 
                w_j_Dunavant(14:19) = 0.043283539377289_pr*1._pr/2._pr

            case (6)
                N_Dunavant = 27
                allocate(w_j_Dunavant(1:N_Dunavant))
                allocate(racines_Dunavant(1:2, 1:N_Dunavant))

! Type n1
                racines_Dunavant(:, 1) = [-0.069222096541517_pr, 0.534611048270758_pr]
                racines_Dunavant(:, 2) = [0.534611048270758_pr, -0.069222096541517_pr]
                racines_Dunavant(:, 3) = [0.534611048270758_pr, 0.534611048270758_pr]
                w_j_Dunavant(1:3) = 0.000927006328961_pr*1._pr/2._pr

! Premier groupe
                racines_Dunavant(:, 4) = [0.202061394068290_pr, 0.398969302965855_pr]
                racines_Dunavant(:, 5) = [0.398969302965855_pr, 0.202061394068290_pr]
                racines_Dunavant(:, 6) = [0.398969302965855_pr, 0.398969302965855_pr]
                w_j_Dunavant(4:6) = 0.077149534914813_pr*1._pr/2._pr

! Deuxieme groupe
                racines_Dunavant(:, 7) = [0.593380199137435_pr, 0.203309900431282_pr]
                racines_Dunavant(:, 8) = [0.203309900431282_pr, 0.593380199137435_pr]
                racines_Dunavant(:, 9) = [0.203309900431282_pr, 0.203309900431282_pr]
                w_j_Dunavant(7:9) = 0.059322977380774_pr*1._pr/2._pr

! Troisieme groupe
                racines_Dunavant(:, 10) = [0.761298175434837_pr, 0.119350912282581_pr]
                racines_Dunavant(:, 11) = [0.119350912282581_pr, 0.761298175434837_pr]
                racines_Dunavant(:, 12)= [0.119350912282581_pr, 0.119350912282581_pr]
                w_j_Dunavant(10:12) = 0.036184540503418_pr*1._pr/2._pr

! Quatrieme groupe
                racines_Dunavant(:, 13) = [0.935270103777448_pr, 0.032364948111276_pr]
                racines_Dunavant(:, 14) = [0.032364948111276_pr, 0.935270103777448_pr]
                racines_Dunavant(:, 15) = [0.032364948111276_pr, 0.032364948111276_pr]
                w_j_Dunavant(13:15) = 0.013659731002678_pr*1._pr/2._pr

! Cinquieme groupe
                racines_Dunavant(:, 16) = [0.050178138310495_pr, 0.356620648261293_pr] 
                racines_Dunavant(:, 17) = [0.050178138310495_pr, 0.593201213428213_pr] 
                racines_Dunavant(:, 18) = [0.356620648261293_pr, 0.050178138310495_pr] 
                racines_Dunavant(:, 19) = [0.356620648261293_pr, 0.593201213428213_pr] 
                racines_Dunavant(:, 20) = [0.593201213428213_pr, 0.050178138310495_pr] 
                racines_Dunavant(:, 21) = [0.593201213428213_pr, 0.356620648261293_pr]
                w_j_Dunavant(16:21) = 0.052337111962204_pr*1._pr/2._pr

! Sixieme groupe
                racines_Dunavant(:, 22) = [0.021022016536166_pr, 0.171488980304042_pr] 
                racines_Dunavant(:, 23) = [0.021022016536166_pr, 0.807489003159792_pr] 
                racines_Dunavant(:, 24) = [0.171488980304042_pr, 0.021022016536166_pr] 
                racines_Dunavant(:, 25) = [0.171488980304042_pr, 0.807489003159792_pr] 
                racines_Dunavant(:, 26) = [0.807489003159792_pr, 0.021022016536166_pr] 
                racines_Dunavant(:, 27) = [0.807489003159792_pr, 0.171488980304042_pr]
                w_j_Dunavant(22:27) = 0.020707659639141_pr*1._pr/2._pr

            case default
                print*, "Erreur : Règle de Dunavant non implémentée pour N_G =", N_G
                stop

            end select

        end subroutine quadrature_Dunavant

        pure function Legendre_2D(k, m, xi, eta) result(Legendre_k_m_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant les polynomes de Legendre 2D normalises

! Entrees de la fonction :
!   - k, m : Degres des polynomes de Legendre 1D
!   - xi : Coordonnee x des elements de reference
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction :
!   - Legendre_k_m_xi_eta : Polynome de Legendre 2D, de degre k suivant xi, m suivant eta, evalue
! en (xi, eta)
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, m

            real(kind = pr), intent(in)                                 :: xi, eta

! Sortie de la fonction
            real(kind = pr)                                             :: Legendre_k_m_xi_eta

            Legendre_k_m_xi_eta = SQRT(real(2*k+1, pr)/2)*Legendre_1D(k, xi)*SQRT(real(2*m+1, pr)/2)*Legendre_1D(m, eta)

        end function Legendre_2D

        pure function Dubiner_2D(k, m, xi, eta) result(Dubiner_k_m_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant la base de Dubiner normalisee
    
! Entrees de la fonction :
!   - k : Degre du polynome de Legendre 1D
!   - m : Degre du polynome de Jacobi 1D
!   - xi : Coordonnee x des elements de reference
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction
!   - Dubiner_k_m_xi_eta : Base de Dubiner 2D evaluee en (xi, eta)
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, m

            real(kind = pr), intent(in)                                 :: xi, eta

! Sortie de la fonction
            real(kind = pr)                                             :: Dubiner_k_m_xi_eta

! Variable locale
            real(kind = pr)                                             :: C_k_m, arg_Legendre

            C_k_m = SQRT(real(2*(2*k+1)*(k + m + 1), pr))/real(2**k, pr)

            if (ABS(1._pr - xi) < eps_var) then
                if (k == 0) then
                    arg_Legendre = 0._pr
                else
                    Dubiner_k_m_xi_eta = 0._pr
                    return
                end if
            else
                arg_Legendre = (2._pr*eta)/(1._pr - xi) - 1._pr
            end if

            Dubiner_k_m_xi_eta = C_k_m*2._pr**k*Legendre_1D(k, arg_Legendre)*(1._pr - xi)**k                                &
        &                                                                             *Jacobi_1D(m, 2*k+1, 0, 2._pr*xi-1._pr)

        end function Dubiner_2D

        pure function Dubiner_xi_2D(k, m, xi, eta) result(Dubiner_xi_k_m_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant la derivee de la base de Dubiner normalisee par rapport a xi
    
! Entrees de la fonction :
!   - k : Degre du polynome de Legendre 1D
!   - m : Degre du polynome de Jacobi 1D
!   - xi : Coordonnee x des elements de reference
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction
!   - Dubiner_xi_k_m_xi_eta : Derivee de la base de Dubiner 2D par rapport a xi evaluee 
! en (xi, eta)
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, m

            real(kind = pr), intent(in)                                 :: xi, eta

! Sortie de la fonction
            real(kind = pr)                                             :: Dubiner_xi_k_m_xi_eta

! Variables locales
            real(kind = pr)                                             :: C_k_m, Li, Li_xi, arg_Li, Ji, Ji_xi, arg_Ji

            C_k_m = SQRT(real(2*(2*k+1)*(k + m + 1), pr))/real(2**k, pr)

            if (ABS(1._pr - xi) < eps_var) then
                Dubiner_xi_k_m_xi_eta = 0._pr
                return
            end if

            arg_Li = (2._pr*eta)/(1._pr - xi) - 1._pr
            arg_Ji = 2._pr*xi - 1._pr

            Li = Legendre_1D(k, arg_Li)
            Li_xi = Legendre_xi_1D(k, arg_Li)
            Ji = Jacobi_1D(m, 2*k+1, 0, arg_Ji)
            Ji_xi = Jacobi_xi_1D(m, 2*k+1, 0, arg_Ji)

            Dubiner_xi_k_m_xi_eta = C_k_m*2._pr**k*(2._pr*eta*Li_xi*(1._pr - xi)**(k-2)*Ji                                  &
        &                                            + Li*(-k)*(1._pr - xi)**MAX(0, k-1)*Ji + Li*(1._pr - xi)**k*2._pr*Ji_xi)

        end function Dubiner_xi_2D

        pure function Dubiner_eta_2D(k, m, xi, eta) result(Dubiner_eta_k_m_xi_eta)
! ==================================================================================================== !
! Fonction renvoyant la derivee de la base de Dubiner normalisee par rapport a eta
    
! Entrees de la fonction :
!   - k : Degre du polynome de Legendre 1D
!   - m : Degre du polynome de Jacobi 1D
!   - xi : Coordonnee x des elements de reference
!   - eta : Coordonnee y des elements de reference

! Sortie de la fonction
!   - Dubiner_eta_k_m_xi_eta : Derivee de la base de Dubiner 2D par rapport a eta evaluee 
! en (xi, eta)
! ==================================================================================================== !

! Entrees de la fonction
            integer, intent(in)                                         :: k, m

            real(kind = pr), intent(in)                                 :: xi, eta

! Sortie de la fonction
            real(kind = pr)                                             :: Dubiner_eta_k_m_xi_eta

! Variables locales
            real(kind = pr)                                             :: C_k_m, Li_xi, arg_Li, Ji, arg_Ji

            C_k_m = SQRT(real(2*(2*k+1)*(k + m + 1), pr))/real(2**k, pr)

            if (ABS(1._pr - xi) < eps_var) then
                Dubiner_eta_k_m_xi_eta = 0._pr
                return
            end if

            arg_Li = (2._pr*eta)/(1._pr - xi) - 1._pr
            arg_Ji = 2._pr*xi - 1._pr

            Li_xi = Legendre_xi_1D(k, arg_Li)
            Ji = Jacobi_1D(m, 2*k+1, 0, arg_Ji)

            Dubiner_eta_k_m_xi_eta = C_k_m*2._pr**k*(2._pr/(1._pr - xi)*Li_xi*(1._pr - xi)**k*Ji)

        end function Dubiner_eta_2D

end module mod_legendre