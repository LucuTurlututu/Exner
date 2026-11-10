module mod_interpolation

    use mod_parametres
    use mod_qte_physiques
    use mod_fonctions
    use mod_types_derives
    use mod_base
    use mod_cl

    implicit none

! ==================================================================================================== !
! Module contenant les routines pour mod_interpolation :
!   - make_B_i
!   - reconstruction_locale
!   - flux_Rusanov
!   - make_Beta_i
!   - make_M_i
!   - make_I_i_EF
!   - make_I_i_VF
!   - make_I_i_S
! ==================================================================================================== !

    contains

        subroutine make_M_i_std(Mesh, Quad_ref, Base_ref, Calc_Exner)
! ==================================================================================================== !
! Subroutine calculant la matrice de masse standard

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
! ==================================================================================================== !

! Entrees de la subroutine
            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: i, j, k, N_ddl_i

            !$omp parallel do private(i, j, k, N_ddl_i) shared(Mesh, Quad_ref, Base_ref, Calc_Exner)

            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%N_Omega_i(i)%M_i_std = 0._pr

                BLOCK
                    real(kind = pr)                                             :: xi_j, eta_j
                    real(kind = pr), dimension(2)                               :: Coord_phy

                    type(Base_Omega_i)                                          :: Base_phy

                    if (Mesh%N_Sommets_Omega_i(i) == 3) then
                        do j = 1, Quad_ref%N_G_Tri
                            xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            !$omp simd
                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%M_i_std(k, :) = Calc_Exner%N_Omega_i(i)%M_i_std(k, :)                   &
                            & + Quad_ref%w_j_Ref_Tri(j)*Base_ref%Phi_i_Tri(k, j)                                                &
                            &                                           *Base_ref%Phi_i_Tri(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                            end do
                        end do
                        
                    else if (Mesh%N_sommets_Omega_i(i) == 4) then
                        do j = 1, Quad_ref%N_G_Quad
                            xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            !$omp simd
                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%M_i_std(k, :) = Calc_Exner%N_Omega_i(i)%M_i_std(k, :)                   &
                            & + Quad_ref%w_j_Ref_Quad(j)*Base_ref%Phi_i_Quad(k, j)                                              &
                            &                                          *Base_ref%Phi_i_Quad(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                            end do
                        end do

                    end if

                END BLOCK
            end do

            !$omp end parallel do

        end subroutine make_M_i_std

        subroutine make_K_i_std(Mesh, Quad_ref, Base_ref, Calc_Exner)
! ==================================================================================================== !
! Subroutine calculant la matrice de rigidite standard

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
! ==================================================================================================== !

! Entrees de la subroutine
            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: i, j, k, N_ddl_i

            !$omp parallel do private(i, j, k, N_ddl_i) shared(Mesh, Quad_ref, Base_ref, Calc_Exner)

            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                Calc_Exner%N_Omega_i(i)%K_i_x_std = 0._pr
                Calc_Exner%N_Omega_i(i)%K_i_y_std = 0._pr

                BLOCK
                    real(kind = pr)                                             :: xi_j, eta_j
                    real(kind = pr), dimension(2)                               :: Coord_phy

                    type(Base_Omega_i)                                          :: Base_phy

                    if (Mesh%N_Sommets_Omega_i(i) == 3) then
                        do j = 1, Quad_ref%N_G_Tri
                            xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j,                             &
                        &                                           j = j, nabla_Tri = Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri)

                            !$omp simd
                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%K_i_x_std(k, :) = Calc_Exner%N_Omega_i(i)%K_i_x_std(k, :)               &
                            & + Quad_ref%w_j_Ref_Tri(j)*ABS(Base_phy%Det_J_i)                                                   &
                            &                        *Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri(k, j, 1)*Base_ref%Phi_i_Tri(:, j)

                                Calc_Exner%N_Omega_i(i)%K_i_y_std(k, :) = Calc_Exner%N_Omega_i(i)%K_i_y_std(k, :)               &
                            & + Quad_ref%w_j_Ref_Tri(j)*ABS(Base_phy%Det_J_i)                                                   &
                            &                        *Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri(k, j, 2)*Base_ref%Phi_i_Tri(:, j)
                            end do
                        end do
                    
                    else if (Mesh%N_sommets_Omega_i(i) == 4) then
                        do j = 1, Quad_ref%N_G_Quad
                            xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j,                             &
                        &                                         j = j, nabla_Quad = Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad)

                            !$omp simd
                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%K_i_x_std(k, :) = Calc_Exner%N_Omega_i(i)%K_i_x_std(k, :)               &
                            & + Quad_ref%w_j_Ref_Quad(j)*ABS(Base_phy%Det_J_i)                                                  &
                            &                      *Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad(k, j, 1)*Base_ref%Phi_i_Quad(:, j)

                                Calc_Exner%N_Omega_i(i)%K_i_y_std(k, :) = Calc_Exner%N_Omega_i(i)%K_i_y_std(k, :)               &
                            & + Quad_ref%w_j_Ref_Quad(j)*ABS(Base_phy%Det_J_i)                                                  &
                            &                      *Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad(k, j, 2)*Base_ref%Phi_i_Quad(:, j)
                            end do
                        end do

                    end if

                END BLOCK
            end do

            !$omp end parallel do

        end subroutine make_K_i_std

        subroutine make_B_i(Mesh, Quad_ref, Base_ref, Calc_Exner)
! ==================================================================================================== !
! Subroutine projetant la bathymetrie initiale sur l'espace DG

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_elements : Nombre de mailles
!   - N_max_sommets : Nombre maximal de sommets d'une maille Omega_i
    
!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Tri : Nombre de points de Gauss une maille triangulaire
!   - N_G_Quad : Nombre de points de Gauss pour une maille quadrangulaire
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnes des points de quadrature de Gauss-Legendre pour
! une maille triangulaire
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnes des points de quadrature de Gauss-Legendre
! pour une maille quadrangulaire
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids de Gauss-Legendre associes pour une maille triangulaire
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids de Gauss-Legendre associes pour une maille
! quadrangulaire

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Tri : Nombre de degre de liberte de la maille triangulaire Omega_i
!   - N_ddl_i_Quad : Nombre de degre de liberte de la maille quadrangulaire Omega_i
!   - Phi_i_Tri(0:N_ddl_i_Tri-1, 1:N_G_Tri) : Tableau contenant les vecteurs de la base de
! Dubiner evalues aux points de Gauss-Legendre (xi_j, eta_j)
!   - Phi_i_Quad(0:N_ddl_i_Quad-1, 1:N_G_Quad) : Tableau contenant les vecteurs de base de 
! Legendre evalues aux points de Gauss-Legendre (xi_j, eta_j)

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - B_i(0:N_ddl_i-1, 1:N_elements) : Tableau contenant les modes de la bathymetrie sur
! l'ensemble des mailles
! ==================================================================================================== !

! Entrees de la subroutine
            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: i, j, N_ddl_i

! B_i(0, :) represente la moyenne sur la maille Omega_i

            if (Mesh%N_max_sommets == 3) then
                if (.NOT. allocated(Calc_Exner%B_i)) allocate(Calc_Exner%B_i(0:Base_ref%N_ddl_i_Tri-1, Mesh%N_elements))
            else if (Mesh%N_max_sommets == 4) then
                if (.not. allocated(Calc_Exner%B_i)) allocate(Calc_Exner%B_i(0:Base_ref%N_ddl_i_Quad-1, Mesh%N_elements))
            else
                print *, "make_B_i : Nombre maximal de sommets non pris en charge"
                return
            end if

            Calc_Exner%B_i = 0._pr

            !$omp parallel do private(i, j, N_ddl_i) shared(Mesh, Quad_ref, Base_ref, Calc_Exner)

            do i = 1, Mesh%N_elements
                N_ddl_i = Base_ref%N_ddl_i(i)

                BLOCK
                    real(kind = pr)                                     :: xi_j, eta_j
                    real(kind = pr), dimension(2)                       :: Coord_phy

                    type(Base_Omega_i)                                  :: Base_phy

                    if (Mesh%N_sommets_Omega_i(i) == 3) then
                        do j = 1, Quad_ref%N_G_Tri
                            xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%B_i(0:N_ddl_i-1, i)                                     &
                        &            + Quad_ref%w_j_Ref_Tri(j)*b(Coord_phy, 0._pr)*Base_ref%Phi_i_Tri(:, j)*ABS(Base_phy%Det_J_i)
                        end do

                    else if (Mesh%N_sommets_Omega_i(i) == 4) then
                        do j = 1, Quad_ref%N_G_Quad
                            xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)

                            Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)
                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            Calc_Exner%B_i(0:N_ddl_i-1, i) = Calc_Exner%B_i(0:N_ddl_i-1, i)                                     &
                        &          + Quad_ref%w_j_Ref_Quad(j)*b(Coord_phy, 0._pr)*Base_ref%Phi_i_Quad(:, j)*ABS(Base_phy%Det_J_i)
                        end do

                    end if

                    call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i_std, Calc_Exner%B_i(0:N_ddl_i-1, i),                 &
                &                                                                                 Calc_Exner%B_i(0:N_ddl_i-1, i))

                END BLOCK
            end do

            !$omp end parallel do

        end subroutine make_B_i

        subroutine reconstruction_locale(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, B_i_F, B_r_F)
! ==================================================================================================== !
! Subroutine renvoyant la projection de la hauteur de la bathymetrie b sur l'espace DG evaluee
! a l'interface F de la maille Omega_i

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - Aretes_Omega_i(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero des aretes de la
! maille Omega_i
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille adjacente
! a la maille Omega_i par la face F
!   - Faces_Omega_r_P(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la face F pour
! la maille Omega_r, i.e. F_r dans le cas de conditions periodiques
! De plus, F_r est negatif si les faces F et F_r ont le meme sens de parcours et positif sinon
!   - CL_Arete_bord(1:N_aretes) : Tableau contenant des lettres referant a la condition aux limites
! appliquee sur le bord :
!       "D" : Condition aux limites de Dirichlet
!       "N" : Conditions aux limites de Neumann
!       "P" : Coniditions aux limites periodiques
!       "O" : L'arete a n'est pas une arete de bord

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Face : Nombre de points de Gauss d'une face F

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Tri : Nombre de degre de liberte de la maille triangulaire Omega_i
!   - N_ddl_i_Quad : Nombre de degre de liberte de la maille quadrangulaire Omega_i
!   - Phi_F_Tri(0:N_ddl_i_Tri-1, 1:N_G_Face, 1:3) : Tableau contenant les vecteurs de la base de
! Dubiner evalues aux points de Gauss-Legendre de la face F
!   - Phi_F_Quad(0:N_ddl_i_Quad-1, 1:N_G_Face, 1:4) : Tableau contenant les vecteurs de la base de 
! Legendre evalues aux points de Gauss-Legendre de la face F

!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - B_i(0:N_ddl_i-1, 1:N_elements) : Tableau contenant les modes de la bathymetrie sur
! l'ensemble des mailles

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i
!   - t : Temps courant

! Sorties de la subroutine :
!   - B_i_F(1:N_G_Face) : Hauteur de la bathymetrie a la face F de la maille Omega_i
!   - B_r_F(1:N_G_Face) : Hauteur de la bathymetrie a la face F de la maille Omega_r
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i, F

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Calcul), intent(in)                                    :: Calc_Exner

! Sorties de la subroutine
            real(kind = pr), dimension(Quad_ref%N_G_Face), intent(out)  :: B_i_F, B_r_F

! Variables locales
            integer                                                     :: a, r, F_r, j, j_r

            real(kind = pr), dimension(0:Base_ref%N_ddl_i(i)-1)         :: B_r

! Interface de la maille Omega_i
            if (Mesh%N_sommets_Omega_i(i) == 3) then
                B_i_F = MATMUL(TRANSPOSE(Base_ref%Phi_F_Tri(:, :, F)), Calc_Exner%B_i(0:Base_ref%N_ddl_i_Tri-1, i))

            else if (Mesh%N_sommets_Omega_i(i) == 4) then
                B_i_F = MATMUL(TRANSPOSE(Base_ref%Phi_F_Quad(:, :, F)), Calc_Exner%B_i(0:Base_ref%N_ddl_i_Quad-1, i))

            end if

! Interface de la maille Omega_r
            B_r_F = 0._pr
            r = Mesh%N_Omega_r(F, i) ; F_r = Mesh%Faces_Omega_r(F, i)

            if (r /= 0) then
! On est sur une arete interieure
                do j = 1, Quad_ref%N_G_Face
! Le sens de parcours de F pour Omega_r est inverse par rapport au sens de parcours de la maille Omega_i
! Ainsi, on inverse les points de Gauss dans la maille adjacente
                    j_r = Quad_ref%N_G_Face - j + 1
                    if(Mesh%N_sommets_Omega_i(r) == 3) then
                        B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j_r, F_r), Calc_Exner%B_i(0:Base_ref%N_ddl_i_Tri-1, r))

                    else if (Mesh%N_sommets_Omega_i(r) == 4) then
                        B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j_r, F_r), Calc_Exner%B_i(0:Base_ref%N_ddl_i_Quad-1, r))
                    end if
                end do

            else
! On est sur une arete de bord
                a = Mesh%Aretes_Omega_i(F, i)
                select case (Mesh%CL_Arete_Bord(a))

                case ("D")
                    B_r = Dirichlet(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t)

                    do j = 1, Quad_ref%N_G_Face
                        j_r = Quad_ref%N_G_Face - j + 1
                        if (Mesh%N_sommets_Omega_i(i) == 3) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j_r, F), B_r)

                        else if (Mesh%N_sommets_Omega_i(i) == 4) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j_r, F), B_r)

                        end if
                    end do

                case ("N")
                    B_r = Neumann(Base_ref, Calc_Exner, i)

                    do j = 1, Quad_ref%N_G_Face
                        j_r = Quad_ref%N_G_Face - j + 1
                        if (Mesh%N_sommets_Omega_i(i) == 3) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j_r, F), B_r)

                        else if (Mesh%N_sommets_Omega_i(i) == 4) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j_r, F), B_r)

                        end if
                    end do

                case ("P")
                    B_r = Periodique(Mesh, Base_ref, Calc_Exner, i, F)

                    r = Mesh%N_Omega_r_P(F, i) ; F_r = ABS(Mesh%Faces_Omega_r_P(F, i))

                    do j = 1, Quad_ref%N_G_Face
                        ! Gestion du sens de parcours selon le signe de Faces_Omega_r_P
                        if (Mesh%Faces_Omega_r_P(F, i) < 0) then
                            j_r = j
                        else
                            j_r = Quad_ref%N_G_Face - j + 1
                        end if

                        if (Mesh%N_sommets_Omega_i(r) == 3) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Tri(:, j_r, F_r), B_r)

                        else if (Mesh%N_sommets_Omega_i(r) == 4) then
                            B_r_F(j) = DOT_PRODUCT(Base_ref%Phi_F_Quad(:, j_r, F_r), B_r)

                        end if
                    end do
                    
                end select
            end if

        end subroutine reconstruction_locale

        subroutine flux_Rusanov(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, hat_G_F)
! ==================================================================================================== !
! Subroutine construisant le flux de Rusanov a la face F

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_Omega_r(1:N_max_sommets, 1:N_elements) : Tableau contenant le numero de la maille
! adjacente a la maille Omega_i par la face F

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Face : Nombre de points de Gauss pour la face de reference
!   - Coord_Ref_Face(1:2, 1:N_G_Face) : Coordonnees des points de Gauss pour la face de
! reference

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega

!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner

!   - i : Entier referant a la maille Omega_i
!   - F : Entier referant a la face F de la maille Omega_i
!   - t : Temps courant

! Sortie de la subroutine :
!   - hat_G_F(1:N_G_Face) : Flux de Rusanov pour la face F de la maille Omega_i en chaque point
! de Gauss
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i, F

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref
            type(Calcul), intent(in)                                    :: Calc_Exner

! Sortie de la subroutine
            real(pr), dimension(1:Quad_ref%N_G_Face), intent(out)       :: hat_G_F

! Variables locales
            integer                                                     :: j

            BLOCK
                real(kind = pr)                                         :: xi_j, alpha_F, h_i_F, h_r_F
                real(kind = pr), dimension(1:2)                         :: Coord_phy_F, u_i_F, u_r_F, G_L, G_R
                real(kind = pr), dimension(1:Quad_ref%N_G_Face)         :: B_i_F, B_r_F

                call reconstruction_locale(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, B_i_F, B_r_F)

                do j = 1, Quad_ref%N_G_Face
                    xi_j = Quad_ref%Coord_Ref_Face(j)
                    Coord_phy_F = hat_phi_F(Mesh, i, F, xi_j)

                    h_i_F = h(Coord_phy_F, t) ; u_i_F = u(Coord_phy_F, t)
                    h_r_F = h(Coord_phy_F, t) ; u_r_F = u(Coord_phy_F, t)

                    G_L = G(B_i_F(j), h_i_F, u_i_F) ; G_R = G(B_r_F(j), h_r_F, u_r_F)

                    alpha_F = MAX(ABS(lambda_DV(h_i_F, u_i_F, Mesh%n_F(:, F, i))),                                          &
                &                                                            ABS(lambda_DV(h_r_F, u_r_F, Mesh%n_F(:, F, i))))

                    hat_G_F(j) = DOT_PRODUCT(1._pr/2*(G_L + G_R), Mesh%n_F(:, F, i)) - alpha_F/2*(B_r_F(j) - B_i_F(j))
                end do
            END BLOCK

        end subroutine flux_Rusanov

        subroutine make_Beta_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant le vecteur de modes Beta_i = (Beta_i_x, Beta_i_y) associe a la
! projection de grad(b) sur l'espace DG

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_faces_Omega_i(1:N_elements) : Tableau contenant le nombre de faces de la maille Omega_i
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i
!   - n_F(1:2, 1:N_max_sommets, 1:N_elements) : Tableau contenant les composantes (n_x, n_y) de la 
! normale unitaire sortante de la face F de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T
!   - w_j_Ref_Face(1:N_G_Face) : Poids d'interpolation associes aux points de Gauss pour la
! face de reference

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - Beta_i(0:N_ddl_i-1) : Vecteur contenant les modes de la projection de grad(b) sur l'espace DG
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: F, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            BLOCK
                real(kind = pr), dimension(1:Quad_ref%N_G_Face)         :: B_i_F, B_r_F, hat_b_F, hat_b_x_F, hat_b_y_F
                real(kind = pr), dimension(0:N_ddl_i-1)                 :: I_i_x_VF, I_i_y_VF
                real(kind = pr), dimension(0:N_ddl_i-1)                 :: Beta_i_x, Beta_i_y

                type(Base_Omega_i)                                      :: Base_phy

! Flux numeriques
                I_i_x_VF = 0._pr ; I_i_y_VF = 0._pr

                do F = 1, Mesh%N_faces_Omega_i(i)
                    call reconstruction_locale(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, B_i_F, B_r_F)
                    
                    hat_b_F = 1._pr/2*(B_i_F + B_r_F)

                    call make_Base_phy(Mesh, Base_ref, i, Base_phy, F = F)

                    hat_b_x_F = Quad_ref%w_j_Ref_Face*hat_b_F*Mesh%n_F(1, F, i)*Base_phy%Det_J_F
                    hat_b_y_F = Quad_ref%w_j_Ref_Face*hat_b_F*Mesh%n_F(2, F, i)*Base_phy%Det_J_F

                    if (Mesh%N_sommets_Omega_i(i) == 3) then
                        I_i_x_VF = I_i_x_VF + MATMUL(Base_ref%Phi_F_Tri(:, :, F), hat_b_x_F)

                        I_i_y_VF = I_i_y_VF + MATMUL(Base_ref%Phi_F_Tri(:, :, F), hat_b_y_F)

                    else if (Mesh%N_sommets_Omega_i(i) == 4) then
                        I_i_x_VF = I_i_x_VF + MATMUL(Base_ref%Phi_F_Quad(:, :, F), hat_b_x_F)

                        I_i_y_VF = I_i_y_VF + MATMUL(Base_ref%Phi_F_Quad(:, :, F), hat_b_y_F)

                    end if
                end do

! Construction du vecteur Beta_i
                Beta_i_x = 0._pr ; Beta_i_y = 0._pr ; Calc_Exner%N_Omega_i(i)%Beta_i = 0._pr

                Beta_i_x = -MATMUL(Calc_Exner%N_Omega_i(i)%K_i_x_std, Calc_Exner%B_i(0:N_ddl_i-1, i)) + I_i_x_VF
                Beta_i_y = -MATMUL(Calc_Exner%N_Omega_i(i)%K_i_y_std, Calc_Exner%B_i(0:N_ddl_i-1, i)) + I_i_y_VF

                call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i_std, Beta_i_x, Calc_Exner%N_Omega_i(i)%Beta_i(:, 1))
                call Cholesky(Base_ref, i, Calc_Exner%N_Omega_i(i)%M_i_std, Beta_i_y, Calc_Exner%N_Omega_i(i)%Beta_i(:, 2))
            END BLOCK

        end subroutine make_Beta_i

        subroutine make_M_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant la matrice de masse locale M_i

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - M_i(0:N_ddl_i-1, 0:N_ddl_i-1) : Matrice de masse locale M_i
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: j, k, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%N_Omega_i(i)%M_i = 0._pr

            BLOCK
                real(kind = pr)                                         :: xi_j, eta_j
                real(kind = pr), dimension(2)                           :: nabla_b_j

                type(Base_Omega_i)                                      :: Base_phy

                call make_Beta_i(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
                
                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    BLOCK
                        real(kind = pr), dimension(Quad_ref%N_G_Tri, 2) :: nabla_b

! Calcul du produit matriciel D_i^T*beta_i
                        nabla_b(:, 1) = MATMUL(TRANSPOSE(Base_ref%Phi_i_Tri), Calc_Exner%N_Omega_i(i)%Beta_i(:, 1))
                        nabla_b(:, 2) = MATMUL(TRANSPOSE(Base_ref%Phi_i_Tri), Calc_Exner%N_Omega_i(i)%Beta_i(:, 2))

! Construction de la matrice de masse
                        do j = 1, Quad_ref%N_G_Tri
                            xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                            nabla_b_j = nabla_b(j, :)

                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%M_i(k, :) = Calc_Exner%N_Omega_i(i)%M_i(k, :)                       &
                            & + Quad_ref%w_j_Ref_Tri(j)*m(nabla_b_j)                                                        &
                            &                        *Base_ref%Phi_i_Tri(k, j)*Base_ref%Phi_i_Tri(:, j)*ABS(Base_phy%Det_J_i)
                            end do
                        end do
                    END BLOCK

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    BLOCK
                        real(kind = pr), dimension(Quad_ref%N_G_Quad, 2) :: nabla_b

! Calcul du produit matriciel D_i^T*beta_i
                        nabla_b(:, 1) = MATMUL(TRANSPOSE(Base_ref%Phi_i_Quad), Calc_Exner%N_Omega_i(i)%Beta_i(:, 1))
                        nabla_b(:, 2) = MATMUL(TRANSPOSE(Base_ref%Phi_i_Quad), Calc_Exner%N_Omega_i(i)%Beta_i(:, 2))

! Construction de la matrice de masse
                        do j = 1, Quad_ref%N_G_Quad
                            xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                            nabla_b_j = nabla_b(j, :)

                            call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                            do k = 0, N_ddl_i-1
                                Calc_Exner%N_Omega_i(i)%M_i(k, :) = Calc_Exner%N_Omega_i(i)%M_i(k, :)                       &
                            & + Quad_ref%w_j_Ref_Quad(j)*m(nabla_b_j)                                                       &
                            &                      *Base_ref%Phi_i_Quad(k, j)*Base_ref%Phi_i_Quad(:, j)*ABS(Base_phy%Det_J_i)
                            end do
                        end do
                    END BLOCK

                end if
            END BLOCK

        end subroutine make_M_i

        subroutine make_I_i_EF(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant le vecteur de contribution volumique du flux I_i_EF

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Quad : Nombre de DDL pour une maille quadrangulaire
!   - N_ddl_i_Tri : Nombre de DDL pour une maille triangulaire
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - I_i_EF(0:N_ddl_i-1) : Vecteur contenant la contribution volumique du flux
! ==================================================================================================== !
            
! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner

! Variables locales
            integer                                                     :: j, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%N_Omega_i(i)%I_i_EF = 0._pr

            BLOCK
                real(kind = pr)                                         :: xi_j, eta_j, b_j, h_j
                real(kind = pr), dimension(2)                           :: Coord_phy, u_j, G_j

                type(Base_Omega_i)                                      :: Base_phy

                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                        h_j = h(Coord_phy, t)
                        u_j = u(Coord_phy, t)

                        G_j = G(b_j, h_j, u_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j,                             &
                    &                                           j = j, nabla_Tri = Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri)

                        Calc_Exner%N_Omega_i(i)%I_i_EF(:) = Calc_Exner%N_Omega_i(i)%I_i_EF(:)                               &
                    & + Quad_ref%w_j_Ref_Tri(j)*ABS(Base_phy%Det_J_i)                                                       &
                    &                              *MATMUL(Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Tri(0:N_ddl_i-1, j, :), G_j)
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                        h_j = h(Coord_phy, t)
                        u_j = u(Coord_phy, t)

                        G_j = G(b_j, h_j, u_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j,                             &
                    &                                         j = j, nabla_Quad = Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad)

                        Calc_Exner%N_Omega_i(i)%I_i_EF(:) = Calc_Exner%N_Omega_i(i)%I_i_EF(:)                               &
                    & + Quad_ref%w_j_Ref_Quad(j)*ABS(Base_phy%Det_J_i)                                                      &
                    &                             *MATMUL(Calc_Exner%N_Omega_i(i)%nabla_i_Phi_i_Quad(0:N_ddl_i-1, j, :), G_j)
                    end do
                end if
            END BLOCK

        end subroutine make_I_i_EF

        subroutine make_I_i_VF(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant le vecteur de contribution surfacique du flux I_i_VF

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_faces_Omega_i(1:N_elements) : Tableau contenant le nombre de faces de la maille Omega_i
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Face : Nombre de points de Gauss pour la face de reference
!   - Coord_Ref_Face(1:2, 1:N_G_Face) : Coordonnees des points de Gauss pour la face de
! reference
!   - w_j_Ref_Face(1:N_G_Face) : Poids d'interpolation associes aux points de Gauss pour la
! face de reference

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Quad : Nombre de DDL pour une maille quadrangulaire
!   - N_ddl_i_Tri : Nombre de DDL pour une maille triangulaire
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_F_Quad(0:N_ddl_i-1, 1:N_G_Face, 1:4) : Vecteur contenant la base de Legendre evaluee aux
! points de Gauss-Legendre des 4 faces F
!   - Phi_F_Tri(0:N_ddl_i-1, 1:N_G_Face, 1:3) : Vecteur contenant la base de Legendre evaluee aux
! points de Gauss-Legendre des 3 faces F

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - I_i_VF : Vecteur contenant la contribution surfacique du flux
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            
! Variables locales
            integer                                                     :: j, F, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%N_Omega_i(i)%I_i_VF = 0._pr

            BLOCK
                real(kind = pr), dimension(1:Quad_ref%N_G_Face)         :: hat_G_F

                type(Base_Omega_i)                                      :: Base_phy

                do F = 1, Mesh%N_faces_Omega_i(i)
                    call make_Base_phy(Mesh, Base_ref, i, Base_phy, F = F)

                    call flux_Rusanov(Mesh, Quad_ref, Base_ref, Calc_Exner, i, F, t, hat_G_F)

                    if (Mesh%N_sommets_Omega_i(i) == 3) then
                        do j = 1, Quad_ref%N_G_Face
                            Calc_Exner%N_Omega_i(i)%I_i_VF(:) = Calc_Exner%N_Omega_i(i)%I_i_VF(:)                           &
                        &                   + Quad_ref%w_j_Ref_Face(j)*hat_G_F(j)                                           &
                        &                                        *Base_ref%Phi_F_Tri(0:N_ddl_i-1, j, F)*ABS(Base_phy%Det_J_F)
                        end do

                    else if (Mesh%N_sommets_Omega_i(i) == 4) then
                        do j = 1, Quad_ref%N_G_Face
                            Calc_Exner%N_Omega_i(i)%I_i_VF(:) = Calc_Exner%N_Omega_i(i)%I_i_VF(:)                           &
                        &                   + Quad_ref%w_j_Ref_Face(j)*hat_G_F(j)                                           &
                        &                                       *Base_ref%Phi_F_Quad(0:N_ddl_i-1, j, F)*ABS(Base_phy%Det_J_F)
                        end do

                    end if

                end do
            END BLOCK

        end subroutine make_I_i_VF

        subroutine make_I_i_S(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant le vecteur de contribution volumique du terme source I_i_S

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Quad : Nombre de DDL pour une maille quadrangulaire
!   - N_ddl_i_Tri : Nombre de DDL pour une maille triangulaire
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - I_i_S(0:N_ddl_i-1) : Vecteur contenant la contribution volumique du terme source
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            
! Variables locales
            integer                                                     :: j, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%N_Omega_i(i)%I_i_S = 0._pr

            BLOCK
                real(kind = pr)                                         :: xi_j, eta_j, b_j, h_j, S_j
                real(kind = pr), dimension(2)                           :: Coord_phy, u_j
                
                type(Base_Omega_i)                                      :: Base_phy

                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Tri(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                        h_j = h(Coord_phy, t)
                        u_j = u(Coord_phy, t)

                        S_j = S(b_j, h_j, u_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        Calc_Exner%N_Omega_i(i)%I_i_S(:) = Calc_Exner%N_Omega_i(i)%I_i_S(:)                                 &
                    &                  + Quad_ref%w_j_Ref_Tri(j)*S_j*Base_ref%Phi_i_Tri(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        b_j = DOT_PRODUCT(Base_ref%Phi_i_Quad(:, j), Calc_Exner%B_i(0:N_ddl_i-1, i))
                        h_j = h(Coord_phy, t)
                        u_j = u(Coord_phy, t)

                        S_j = S(b_j, h_j, u_j)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        Calc_Exner%N_Omega_i(i)%I_i_S(:) = Calc_Exner%N_Omega_i(i)%I_i_S(:)                                 &
                    &                + Quad_ref%w_j_Ref_Quad(j)*S_j*Base_ref%Phi_i_Quad(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                    end do

                end if
            END BLOCK

        end subroutine make_I_i_S

        subroutine make_I_i_R_tilde(Mesh, Quad_ref, Base_ref, Calc_Exner, i, t)
! ==================================================================================================== !
! Subroutine construisant le vecteur de contribution volumique du terme source pour la 
! solution manufacturee I_i_R_tilde

! Entrees de la subroutine :
!   - Mesh : Type derive contenant l'ensemble des informations liees au maillage
!   - N_sommets_Omega_i(1:N_elements) : Tableau contenant le nombre de sommets de la maille Omega_i

!   - Quad_ref : Type derive contenant les tableaux lies a la quadrature de Gauss sur l'element
! de reference hat_Omega
!   - N_G_Quad : Nombre de points de Gauss pour la maille de reference Q
!   - N_G_Tri : Nombre de points de Gauss pour la maille de reference T
!   - Coord_Ref_Quad(1:2, 1:N_G_Quad) : Coordonnees des points de Gauss pour la maille de
! reference Q
!   - Coord_Ref_Tri(1:2, 1:N_G_Tri) : Coordonnees des points de Gauss pour la maille de
! reference T
!   - w_j_Ref_Quad(1:N_G_Quad) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference Q
!   - w_j_Ref_Tri(1:N_G_Tri) : Poids d'interpolation associes aux points de Gauss pour la
! maille de reference T

!   - Base_ref : Type derive contenant les vecteurs de base sur l'element de reference hat_Omega
!   - N_ddl_i_Quad : Nombre de DDL pour une maille quadrangulaire
!   - N_ddl_i_Tri : Nombre de DDL pour une maille triangulaire
!   - N_ddl_i(1:N_elements) : Tableau contenant le nombre de DDL de la maille Omega_i
!   - Phi_i_Quad(0:N_ddl_i-1, 1:N_G_Quad) : Vecteur contenant la base de Legendre evaluee aux points de
! Gauss-Legendre
!   - Phi_i_Tri(0:N_ddl_i-1, 1:N_G_Tri) : Vecteur contenant la base de Dubiner evaluee aux points de
! Gauss-Legendre

!   - i : Entier referant a la maille Omega_i
!   - t : Temps courant

! Entree/Sortie de la subroutine :
!   - Calc_Exner : Type derive contenant les tableaux lies a la resolution numerique de l'equation
! d'Exner
!   - I_i_R_tilde(0:N_ddl_i-1) : Vecteur contenant la contribution volumique du terme source
! pour la solution manufacturee
! ==================================================================================================== !

! Entrees de la subroutine
            integer, intent(in)                                         :: i

            real(kind = pr), intent(in)                                 :: t

            type(Maillage), intent(in)                                  :: Mesh
            type(Quadrature_hat_Omega), intent(in)                      :: Quad_ref
            type(Base_hat_Omega), intent(in)                            :: Base_ref

! Entree/Sortie de la subroutine
            type(Calcul), intent(inout)                                 :: Calc_Exner
            
! Variables locales
            integer                                                     :: j, N_ddl_i

            N_ddl_i = Base_ref%N_ddl_i(i)

            Calc_Exner%N_Omega_i(i)%I_i_R_tilde = 0._pr

            BLOCK
                real(kind = pr)                                         :: xi_j, eta_j, R_tilde_j
                real(kind = pr), dimension(2)                           :: Coord_phy
                
                type(Base_Omega_i)                                      :: Base_phy

                if (Mesh%N_sommets_Omega_i(i) == 3) then
                    do j = 1, Quad_ref%N_G_Tri
                        xi_j = Quad_ref%Coord_Ref_Tri(1, j) ; eta_j = Quad_ref%Coord_Ref_Tri(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        R_tilde_j = R_tilde(Coord_phy, t)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        Calc_Exner%N_Omega_i(i)%I_i_R_tilde(:) = Calc_Exner%N_Omega_i(i)%I_i_R_tilde(:)                         &
                    &                + Quad_ref%w_j_Ref_Tri(j)*R_tilde_j*Base_ref%Phi_i_Tri(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                    end do

                else if (Mesh%N_sommets_Omega_i(i) == 4) then
                    do j = 1, Quad_ref%N_G_Quad
                        xi_j = Quad_ref%Coord_Ref_Quad(1, j) ; eta_j = Quad_ref%Coord_Ref_Quad(2, j)
                        Coord_phy = hat_phi_i(Mesh, i, xi_j, eta_j)

                        R_tilde_j = R_tilde(Coord_phy, t)

                        call make_Base_phy(Mesh, Base_ref, i, Base_phy, xi = xi_j, eta = eta_j)

                        Calc_Exner%N_Omega_i(i)%I_i_R_tilde(:) = Calc_Exner%N_Omega_i(i)%I_i_R_tilde(:)                         &
                    &              + Quad_ref%w_j_Ref_Quad(j)*R_tilde_j*Base_ref%Phi_i_Quad(0:N_ddl_i-1, j)*ABS(Base_phy%Det_J_i)
                    end do

                end if
            END BLOCK

        end subroutine make_I_i_R_tilde

end module mod_interpolation