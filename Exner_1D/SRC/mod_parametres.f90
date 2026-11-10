module mod_parametres

    implicit none

    integer, parameter                              :: pr = 8

! Nombre pi
    real(kind = pr), parameter                      :: pi = 4._pr*ATAN(1._pr)

! Temps final de simulation
! t_fnal = 1 pour les cas tests 1, 2, 3, 4, 5, 6
! t_final = 6000 pour le cas test 7
! t_final = 7 pour le cas test 8, 9
! t_final = 2000 pour le test 10
    real(kind = pr), parameter                      :: t_final = 6000._pr

! Longueur du domaine physique
! L_dom = 1 pour les cas tests 1, 2, 3, 4, 5, 6
! L_dom = 300 pour le cas test 7
! L_dom = 7 pour le cas test 8, 9
! L_dom = 400 pour le cas test 10
    real(kind = pr), parameter                      :: L_dom = 400._pr

! Longueur du domaine d'absorption
! L_spo = 50 pour le cas test 10
! L_spo = 0 sinon
    real(kind = pr), parameter                      :: L_spo = 0._pr

! Longueur totale du domaine
    real(kind = pr), parameter                      :: L = L_dom + L_spo

! Nombre de points de discretisation du domaine physique
    integer, parameter                              :: N_dom = 300

! Nombre de points de discretation du domaine d'absorption
! N_spo = 50 pour le cas test 10
! N_spo = 0 sinon
    integer, parameter                              :: N_spo = 0

! Nombre de points total du domaine
    integer, parameter                              :: N_x = N_dom + N_spo

! Contrainte critique
! tau_c = 1 pour les cas tests 1, 2, 3
! tau_c = 10 pour le cas test 4, 5, 6
! tau_c n'intervient pas dans le cas test 7, 8, 9
    real(kind = pr), parameter                      :: tau_c = 10._pr

! Degre polynomial
    integer, parameter                              :: p = 2

! Nombre de points pour la quadrature de Gauss
    integer, parameter                              :: N_G = p + 3

! Seuil de tolerence pour les racines des polynomes de Legendre
    real(kind = pr), parameter                      :: eps_racine = 1e-12_pr

! Ecart pour la difference finie calculant la derivee de S dans find_dt
    real(kind = pr), parameter                      :: eps_dt = 1e-8_pr

! Facteur de surete CFL
    real(kind = pr), parameter                      :: CFL = 0.9_pr

! Seuil en dessous duquel on considere que u_max ou S_b sont a l'erreur
! machine
    real(kind = pr), parameter                      :: eps_var = 1e-14_pr

! Parametre de lissage
    real(kind = pr), parameter                      :: eps_tau = 0.1_pr*tau_c

! Hauteur moyenne pour la bathymetrie
    real(kind = pr), parameter                      :: b_0 = 0.5_pr

! Entier pour le cas test
    integer                                         :: cas

! Entier pour le schema
    integer                                         :: schema

! Coefficient d'optimisation du schema SSPRK
    real(kind = pr)                                 :: alpha_RK

end module mod_parametres