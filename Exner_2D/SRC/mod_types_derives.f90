module mod_types_derives

    use mod_parametres

    implicit none

! ==================================================================================================== !
! Module contenant les types derives :
!   - Quadrature_hat_Omega : Type derive contenant les points et les poids de quadrature sur
! l'element de reference
!   - Base_hat_Omega : Type derive contenant les vecteurs de la base de reference
!   - Base_Omega_i : Type derive contenant le Jacobien et les transformations pour la maille
! Omega_i
!   - Calcul : Type derive de calcul pour l'equation d'Exner
! ==================================================================================================== !

    type :: Maillage
        character(len = 1), dimension(:), allocatable                   :: CL_Arete_bord

        integer                                                         :: N_elements
        integer                                                         :: N_aretes
        integer                                                         :: N_max_sommets
        
        integer, dimension(:), allocatable                              :: N_faces_Omega_i
        integer, dimension(:), allocatable                              :: N_sommets_Omega_i
        integer, dimension(:, :), allocatable                           :: Sommets_Omega_i
        integer, dimension(:, :), allocatable                           :: Aretes_Omega_i

        integer, dimension(:, :), allocatable                           :: N_Omega_r
        integer, dimension(:, :), allocatable                           :: Faces_Omega_r
        integer, dimension(:, :), allocatable                           :: N_Omega_r_P
        integer, dimension(:, :), allocatable                           :: Faces_Omega_r_P

        real(kind = pr), dimension(:), allocatable                      :: Aire_Omega_i
        real(kind = pr), dimension(:), allocatable                      :: Long_arete
        real(kind = pr), dimension(:), allocatable                      :: Diam_Omega_i
        real(kind = pr), dimension(:, :), allocatable                   :: Coord_Sommets
        real(kind = pr), dimension(:, :, :), allocatable                :: n_F
    end type Maillage

    type :: Quadrature_hat_Omega
        integer                                                         :: N_G_Quad
        integer                                                         :: N_G_Tri
        integer                                                         :: N_G_Face

        real(kind = pr), dimension(:), allocatable                      :: w_j_Ref_Quad
        real(kind = pr), dimension(:), allocatable                      :: w_j_Ref_Tri

        real(kind = pr), dimension(:), allocatable                      :: w_j_Ref_Face
        real(kind = pr), dimension(:), allocatable                      :: Coord_Ref_Face

        real(kind = pr), dimension(:, :), allocatable                   :: Coord_Ref_Quad
        real(kind = pr), dimension(:, :), allocatable                   :: Coord_Ref_Tri
    end type Quadrature_hat_Omega

    type :: Base_hat_Omega
        integer                                                         :: N_ddl_i_Quad
        integer                                                         :: N_ddl_i_Tri
        integer, dimension(:), allocatable                              :: N_ddl_i

        real(kind = pr), dimension(:, :), allocatable                   :: Phi_i_Quad
        real(kind = pr), dimension(:, :), allocatable                   :: Phi_i_Tri

        real(kind = pr), dimension(:, :), allocatable                   :: Grad_xi_Quad
        real(kind = pr), dimension(:, :), allocatable                   :: Grad_eta_Quad
        real(kind = pr), dimension(:, :), allocatable                   :: Grad_xi_Tri
        real(kind = pr), dimension(:, :), allocatable                   :: Grad_eta_Tri

        real(kind = pr), dimension(:, :, :), allocatable                :: Phi_F_Quad
        real(kind = pr), dimension(:, :, :), allocatable                :: Phi_F_Tri
    end type Base_hat_Omega

    type :: Base_Omega_i
        real(kind = pr)                                                 :: Det_J_i
        real(kind = pr)                                                 :: Det_J_F

        real(kind = pr), dimension(2, 2)                                :: Inv_J_i_T
    end type Base_Omega_i

    type :: Cellule
        integer, dimension(:), allocatable                              :: ipiv

        real(kind = pr), dimension(:), allocatable                      :: I_i_EF
        real(kind = pr), dimension(:), allocatable                      :: I_i_VF
        real(kind = pr), dimension(:), allocatable                      :: I_i_S
        real(kind = pr), dimension(:), allocatable                      :: I_i_R_tilde

        real(kind = pr), dimension(:), allocatable                      :: L_i_E
        real(kind = pr), dimension(:), allocatable                      :: L_i_I

        real(kind = pr), dimension(:), allocatable                      :: B_i_n

        real(kind = pr), dimension(:, :), allocatable                   :: L_i_s_E
        real(kind = pr), dimension(:, :), allocatable                   :: L_i_s_I

        real(kind = pr), dimension(:, :), allocatable                   :: Beta_i
        real(kind = pr), dimension(:, :), allocatable                   :: M_i

        real(kind = pr), dimension(:, :), allocatable                   :: M_i_std
        real(kind = pr), dimension(:, :), allocatable                   :: K_i_x_std
        real(kind = pr), dimension(:, :), allocatable                   :: K_i_y_std

        real(kind = pr), dimension(:, :), allocatable                   :: J_i_q

        real(kind = pr), dimension(:, :, :), allocatable                :: nabla_i_Phi_i_Quad
        real(kind = pr), dimension(:, :, :), allocatable                :: nabla_i_Phi_i_Tri
    end type Cellule

    type :: Schema_Temporel
        integer                                                         :: m_etapes

        integer                                                         :: q_Newton
        integer                                                         :: q_max_Newton
        integer                                                         :: N_appels_Newton

        integer                                                         :: q_max_Picard

        real(kind = pr), dimension(:), allocatable                      :: b_Butcher
        real(kind = pr), dimension(:), allocatable                      :: b_tilde_Butcher

        real(kind = pr), dimension(:), allocatable                      :: c_Butcher
        real(kind = pr), dimension(:), allocatable                      :: c_tilde_Butcher

        real(kind = pr), dimension(:, :), allocatable                   :: A_Butcher
        real(kind = pr), dimension(:, :), allocatable                   :: A_tilde_Butcher
    end type Schema_Temporel
    
    type :: Calcul
        real(kind = pr)                                                 :: dt

        real(kind = pr), dimension(:, :), allocatable                   :: B_i

        type(Cellule), dimension(:), allocatable                        :: N_Omega_i
    end type Calcul

end module mod_types_derives