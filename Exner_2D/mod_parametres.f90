module mod_parametres

    implicit none

! Parametre de precision pour les reels
    integer, parameter                              :: pr = 8

! Degre polynomial de l'espace DG
    integer, parameter                              :: p = 2

! Nombre de points de Gauss dans une direction de l'espace, i.e. sur les faces
    integer, parameter                              :: N_G = p + 3

! Nom du fichier de maillage
!   - "mesh" : Maillages du carre [0, 1]^2
!   - "berthon" : Maillages du carre [0, 7]^2 pour les cas test de Berthon
!   - "[...]1" : Maillages de Delaunay
!   - "[...]2" : Maillages cartésiens
!   - "[...]_1:5" : Niveau du raffinement du maillage ; 1 = + Grossier, 5 = + Fin
    character(len = *), parameter                   :: fichier = TRIM(ADJUSTL("./MESH/mesh2_3.msh"))

! Temps final
!   - Cas test 1:3 : T_final = 0.5_pr
!   - Cas test 4:6 : T_final = 1._pr
!   - Cas test 7:9 : T_final = 7._pr
    real(kind = pr), parameter                      :: T_final = 0.5_pr

! Facteur de surete CFL
    real(kind = pr), parameter                      :: alpha_CFL = 0.9_pr

! Seuil de Shields
    real(kind = pr), parameter                      :: theta_c = 0.047_pr

! Parametre de lissage des indicatrices
    real(kind = pr), parameter                      :: eps_theta = 0.1_pr*theta_c

! Seuil de tolerance en dessous duquel on considere que nous sommes au zero marchine
    real(kind = pr), parameter                      :: eps_var = 1.e-14_pr

! Crietere d'arret de la methode de Newton trouvant racines des polynomes de Legendre
! en 1D
    real(kind = pr), parameter                      :: eps_racine = 1.e-12_pr

! Parametre de variance pour effectuer une DF numerique
    real(kind = pr), parameter                      :: eps_DF = 1.e-8_pr

! Critere d'arret de la methode de Newton pour la resolution IMEX
    real(kind = pr), parameter                      :: eps_Newton = 1.e-12_pr

! Seuil en deca duquel on considre qu'on a trouve la face periodique
    real(kind = pr), parameter                      :: eps_P = 1.e-8_pr

! Nombre maximal d'iterations pour une boucle de l'algorithme de Newton
    integer, parameter                              :: N_max_Newton = 20

! Nombre maximal d'iterations pour une boucle de l'algorithme de Picard
    integer, parameter                              :: N_max_Picard = 30

! Nombre de Froude maximal
    real(kind = pr)                                 :: Fr_max = 0._pr

! Numero du cas test
    integer                                         :: cas

! Numero du schema
    integer                                         :: schema

! Longueur du domaine
    real(kind = pr)                                 :: Long

! Largeur du domaine
    real(kind = pr)                                 :: Larg

! Constante pi
    real(kind = pr), parameter                      :: pi = 4._pr*ATAN(1._pr)

end module mod_parametres