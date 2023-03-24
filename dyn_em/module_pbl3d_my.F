  module module_pbl3d_my

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                         !
    ! Purpose: 1) Calculates turbulent fluxes of momentum, heat and moisture  !
    !                                                                         !
    !          2) Calculates tendencies due to 3D turbulent mixing            !
    !                                                                         !
    ! Author: Pedro A. Jimenez & Timothy W. Juliano                           !
    !                                                                         !
    ! Methods: Uses the Mellor and Yamada turbulence closure model            !
    !                                                                         !
    ! References: M73: Mellor, G. L., 1973: Analytic prediction of the        !
    !                  properties of stratified planetary surface layers. J.  !
    !                  Atmos. Sci., 30, 1061-1069.                            !
    !                                                                         !
    !             MY74: Mellor, G. L., and T. Yamada, 1974: A hierarchy of    !
    !                   turbulence closure models for planetary boundary      !
    !                   layers. J. Atmos. Sci., 31, 1791-1806.                !
    !                                                                         !
    !             YM75: Yamada, T., and G. Mellor, 1975: A simulation of the  !
    !                   Wangara atmospheric boundary layer data. J. Atmos.    !
    !                   Sci., 32, 2309-2329.                                  !
    !                                                                         !
    !              Y75: Yamada, T., The critical richardson number and the    !
    !                   ratio of the eddy transport coefficients obtained     !
    !                   from a turbulence closure model. JAS, 32, 926 - 933   !
    !                                                                         !
    !             MY82: Mellor, G. L., and T. Yamada, 1982: Development of a  !
    !                   turbulent closure model for geophysical fluid         !
    !                   problems. Rev. Geophys.                               !
    !                   Space Phys., 20, 851-875.                             !
    !                                                                         !
    !              Y83: Yamada, T., 1983: Simulations of nocturnal drainage   !
    !                   flows by a q2/l turbulence closure model              !
    !                                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use, intrinsic :: iso_fortran_env

    use module_model_constants, only: KARMAN, G, P608, T0, CP, RCP
    use module_state_description, only: p_qc, p_qr, p_qi, p_qs, p_qg, p_qnc, p_qnr, p_qni, p_qnwfa, p_qnifa, param_first_scalar
    use module_big_step_utilities_em, only: grid_config_rec_type

    implicit none

    private

    public :: Set_init_turb_state_my, Calc_turb_fluxes_my

    real, parameter :: G_OVER_TREF = G / T0
    real, parameter :: DWIND_DZ_MIN = 1.0E-10
    real, public    :: Q_SQ_MIN, TURB_FLUX_MIN
    logical, parameter :: USE_HIGH_ORDER = .true.

    integer, parameter :: SP = REAL32
    integer, parameter :: DP = REAL64
    integer, parameter :: SOLVER_PREC = DP
    integer, parameter :: WORKING_PREC = SP

      ! Model constants
    real         :: a_1, b_2, c_1, c_2, c_3, alpha, lambda
    real, public :: a_2, b_1

      ! For the level 2 model
    real :: rif_c, rf1, rf2, cm, ch2
    real :: fac1, fac2, fac3, fac4

    integer, parameter :: I_TO_TEST = 4, J_TO_TEST = 4, K_TO_TEST = 3

  contains

    subroutine Calc_turb_fluxes_my (config_flags, dz, z_at_mass, z_at_w, u, v, w, th, th2, qv, thetav, tsk, t2, q2, psfc, rho, ust, hfx, qfx, rmol, &
              du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dqv_dx, dqv_dy, dqv_dz, dthetav_dx, dthetav_dy, dthetav_dz, &
              dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz, itimestep, ss, n_tracer, q_sq, q_sq_prog, &
              l_master, l_master_at_mass, l_boulac, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
              turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv, &
              turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqc, turb_flux_vqc, &
              turb_flux_wqc, turb_flux_uqr, turb_flux_vqr, turb_flux_wqr, turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, &
              turb_flux_uqs, turb_flux_vqs, turb_flux_wqs, turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, turb_flux_uqnc, turb_flux_vqnc, &
              turb_flux_wqnc, turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_wqni, &
              turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa, &
              mat_cond_heat, mat_cond_moist, rpgf_heat, rpgf_moist, msfux, msfuy, msfvx, msfvy, msftx, msfty, xland, pblh, &
              dx, dy, rdx, rdy, dn, dnw, rdz, dt, id, rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, mut, c1h, c2h, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Calculates turbulent fluxes                                  !
      !                                                                       !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                         !
      !                                                                       !
      ! Method: Mellor and Yamada model                                       !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags
      integer, intent(in) :: itimestep, id, ss, n_tracer

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, z_at_mass, z_at_w, u, v, &
          w, th, th2, qv, rho, zx, zy, rdz, rdzw
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq, q_sq_prog
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, &
          dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, dqc_dz, dqr_dz,    &
          dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master, l_master_at_mass, l_boulac
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, jms:jme),          intent (in)    :: tsk, t2, thetav, q2, psfc, ust, hfx, qfx, rmol, xland
      real, dimension (ims:ime, jms:jme),          intent (inout) :: pblh
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1h, c2h
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: turb_flux_u2, turb_flux_v2, turb_flux_w2,  &
          turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta,  &
          turb_flux_wtheta, turb_flux_uqc, turb_flux_vqc, turb_flux_wqc, turb_flux_uqr, turb_flux_vqr,          &
          turb_flux_wqr, turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, turb_flux_uqs, turb_flux_vqs,             &
          turb_flux_wqs, turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, turb_flux_uqnc, turb_flux_vqnc,           &
          turb_flux_wqnc, turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, turb_flux_uqni, turb_flux_vqni,       &
          turb_flux_wqni, turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, turb_flux_uqnifa,               &
          turb_flux_vqnifa, turb_flux_wqnifa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: mat_cond_heat, mat_cond_moist, rpgf_heat, rpgf_moist
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp, dn, dnw
      real,                                        intent (in)    :: dx, dy, rdx, rdy, cf1, cf2, cf3, dt
      integer,                                     intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! Local vars
      logical, parameter :: DEBUG = .false.
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf

      if (config_flags%pbl3d_opt < 2) then  ! PBL approx analytical solution only
        call Calc_fluxes_pbl_approx (config_flags, ss, l_master, l_boulac, du_dz, dv_dz, dthetav_dz, dqv_dz, &
            dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz, &
            th2, qv, dn, dnw, fnm, fnp, dx, dy, dz, z_at_w, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
            jte, kts, kte, q_sq, q_sq_prog, turb_flux_u2, turb_flux_v2, turb_flux_w2,     &
            turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                 &
            turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv,    &
            turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
            turb_flux_uqc, turb_flux_vqc, turb_flux_wqc, turb_flux_uqr, turb_flux_vqr, turb_flux_wqr, &
            turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, turb_flux_uqs, turb_flux_vqs, turb_flux_wqs, &
            turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, turb_flux_uqnc, turb_flux_vqnc, turb_flux_wqnc, &
            turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_wqni, &
            turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa, &
            u, v, w, thetav, rho, ust, hfx, qfx, rmol)
      else if (config_flags%pbl3d_opt == 2) then  ! Full 3D numerical solution
        call Calc_fluxes (config_flags, id, config_flags%specified, q_sq, q_sq_prog, l_master, l_boulac, &
            du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, &
            dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, &
            turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
            turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, mat_cond_heat, mat_cond_moist, rpgf_heat, rpgf_moist, &
            th2, thetav, qv, rho, hfx, qfx, rmol, u, v, ust, dz, rdzw, z_at_w, fnm, fnp, cf1, cf2, cf3, dn, dnw, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Diagnose PBL height
      call Diagnose_pbl_height (thetav, q_sq, z_at_w, dz, xland, pblh, ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      if (DEBUG) print *, 'Done Calc_turb_fluxes_my...'

    end subroutine Calc_turb_fluxes_my


    subroutine Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_prog
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real,                                        intent (in)  :: cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq

        ! Local vars
      integer :: i, j, k
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real    :: cft1, cft2

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq(i, k, j) = Max (Q_SQ_MIN, fnm(k) * q_sq_prog(i, k, j) + fnp(k) * q_sq_prog(i, k - 1, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          q_sq(i, kts, j) = Max (Q_SQ_MIN, q_sq(i, kts, j))
            ! Top face
          q_sq(i, ktf + 1, j) = Q_SQ_MIN
        end do
      end do

    end subroutine Fill_q_sq_with_q_sq_prog


    subroutine Fill_q_sq_prog_with_q_sq (q_sq_prog, q_sq, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_prog

        ! Local vars
      integer :: i, j, k
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real    :: cft1, cft2

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_prog(i, k, j) = 0.5 * ( q_sq(i, k + 1, j) + q_sq(i, k, j) )
          end do
        end do
      end do

    end subroutine Fill_q_sq_prog_with_q_sq


    subroutine Calc_l_master_algebra (pbl3d_l_opt, q_sq, dz, rdzw, rmol, wthv, dthetav_dz, l_master, &
              l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the master length scale following MY74  !
        !          Also included MYNN length scale formulation        !
        !                                                             !
        ! Author: Pedro A. Jimenez & Timothy W. Juliano               !
        !                                                             !
        ! Comments: Use an algebraic equation to calculate L          !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer,                                     intent (in) :: pbl3d_l_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: dz, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in) :: rmol
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: wthv
      real, dimension (its:ite, kts:kte, jts:jte), intent (in) :: dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: l_master, l_boulac
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local vars
      real, dimension (kts:kte) :: z_at_walls_1d
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, parameter :: SMALL_VAL = 0.00001
!      real, parameter :: alpha_1 = 0.23
        ! Original MYNN vals
!      real, parameter :: alpha_2 = 1.0
!      real, parameter :: alpha_3 = 5.0
!      real, parameter :: alpha_4 = 100.0
!      real, parameter :: cns = 2.7
        ! Operational MYNN vals
!      real, parameter :: alpha_2 = 0.65
!      real, parameter :: alpha_3 = 3.0
!      real, parameter :: alpha_4 = 20.0
!      real, parameter :: cns = 2.3
        ! Experimental MYNN vals
!      real, parameter :: alpha_2 = 0.3
!      real, parameter :: alpha_3 = 2.0
!      real, parameter :: alpha_4 = 10.0
!      real, parameter :: cns = 3.5
        ! 3DTKE MYNN vals
      real, parameter :: alpha_1 = 0.8
      real, parameter :: alpha_2 = 1.0
      real, parameter :: alpha_3 = 1.0
      real, parameter :: alpha_4 = 100.0
      real, parameter :: cns = 2.7
      real :: l0, l0_num, l0_den, qdz
      real :: q, q_dz, dz_wall
      real :: l_s, l_b, l_f, zet, kz, N, q_c

      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

         ! Calculates l0 (Eq. 72 MY74)
      do j = j_start, j_end
        do i = i_start, i_end
          l0_num = SMALL_VAL
          l0_den = SMALL_VAL
          z_at_walls_1d(kts) = 0.0
          do k = kts + 1, ktf
            z_at_walls_1d(k) = z_at_walls_1d(k - 1) + 1.0 / rdzw(i, k - 1, j)
            dz_wall = 0.5 * (dz(i, k, j) + dz(i, k - 1, j))
            q = Sqrt (q_sq(i, k, j))
            q_dz = q * dz_wall
            l0_num = l0_num + q_dz * z_at_walls_1d(k)
            l0_den = l0_den + q_dz
          end do
          z_at_walls_1d(ktf + 1) = z_at_walls_1d(ktf) + 1.0 / rdzw(i, ktf, j)

          if (pbl3d_l_opt == 1) then  ! MY alpha
            l0 = alpha * l0_num / l0_den
          else if (pbl3d_l_opt == 2) then  ! MYNN alpha
            l0 = alpha_1 * l0_num / l0_den
            q_c = (G_OVER_TREF * Max (wthv(i, kts, j), 0.0) * l0) ** ( 1.0 / 3.0 )
          end if

          do k = kts, ktf
            if (pbl3d_l_opt == 1) then
               ! Calculates master length scale (Eq. 71 MY74)
              l_master(i, k, j) = l0 * KARMAN * z_at_walls_1d(k) / (KARMAN * z_at_walls_1d(k) + l0)
            else if (pbl3d_l_opt == 2) then
               ! Calculates master length scale (Nakanishi 2001)
               ! l surface
              kz  = KARMAN * z_at_walls_1d(k)
              zet = rmol(i, j) * z_at_walls_1d(k)
              if (zet > 0.0) then
                l_s = kz / (1.0 + cns * MIN(1.0, zet))
              else
                l_s = kz * (1.0 - alpha_4 * zet) ** 0.2
              end if
               ! l buoyancy
              if (dthetav_dz(i, k, j) > 0.0) then
                N = Sqrt (G_OVER_TREF * dthetav_dz(i, k, j))
                 ! note q_c is a velocity scale defined on the ground hence kts index
                l_b = alpha_2 * Sqrt (q_sq(i, k, j)) / N * ( 1.0 + alpha_3 / alpha_2 * Sqrt (q_c / (l0 * N)) )
                l_f = alpha_2 * Sqrt (q_sq(i, k, j)) / N
              else
                l_b = 1.0E10
                l_f = l_b
              end if

              l_master(i, k, j) = MIN( (l_s * l0 * l_b) / (l_s * l0 + l_s * l_b + l0 * l_b), l_f )
              l_boulac(i, k, j) = l_master(i, k, j)  ! here we are just filling in dissipation boulac length scale

            end if
          end do
          l_master(i, kte, j) = l_master(i, ktf, j)
          l_boulac(i, kte, j) = l_boulac(i, ktf, j)

          if (DEBUG) call Debug_l_my_algebra

        end do
      end do

      contains

        subroutine Debug_l_master_algebra

          if (i == I_TO_TEST .and. j == J_TO_TEST) then
              print *, 'l_master = '
              print *, l_master(I_TO_TEST, :, J_TO_TEST)
              print *, 'l0 =', l0
              print *, 'l0_num =', l0_num
              print *, 'l0_den =', l0_den
              print *, 'z_at_walls = '
              print *, z_at_walls_1d
              print *, 'dz = '
              print *, dz(i, :, j)
              print *, 'rdzw = '
              print *, rdzw(i, :, j)
              print *, 'dzw = '
              print *, 1.0 / rdzw(i, :, j)
            end if

        end subroutine Debug_l_master_algebra

    end subroutine Calc_l_master_algebra


    subroutine Calc_l_messinger_master_algebra (q_sq, dz, rdzw, rmol, wthv, dthetav_dz, l_master, &
              l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the master length scale following MY74  !
        !          Also included MYNN length scale formulation        !
        !                                                             !
        ! Author: Pedro A. Jimenez & Timothy W. Juliano               !
        !                                                             !
        ! Comments: Use an algebraic equation to calculate L          !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: dz, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in) :: rmol
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: wthv
      real, dimension (its:ite, kts:kte, jts:jte), intent (in) :: dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: l_master, l_boulac
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local vars
      real, dimension (kts:kte) :: z_at_walls_1d
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, parameter :: SMALL_VAL = 0.00001
      real, parameter :: alpha_1 = 0.23
      real :: l0, l0_num, l0_den, qdz
      real :: q, q_dz, dz_wall
      real :: l_b, l_d, l_del, N, zi
      real, parameter :: c_r = 0.25
      real, parameter :: c_g = 10.
      real, parameter :: c_s = 0.23

      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

         ! Calculates l0 (Eq. 72 MY74)
      do j = j_start, j_end
        do i = i_start, i_end
          l0_num = SMALL_VAL
          l0_den = SMALL_VAL
          z_at_walls_1d(kts) = 0.0
          do k = kts + 1, ktf
            z_at_walls_1d(k) = z_at_walls_1d(k - 1) + 1.0 / rdzw(i, k - 1, j)
            dz_wall = 0.5 * (dz(i, k, j) + dz(i, k - 1, j))
            q = Sqrt (q_sq(i, k, j))
            q_dz = q * dz_wall
            l0_num = l0_num + q_dz * z_at_walls_1d(k)
            l0_den = l0_den + q_dz
          end do
          z_at_walls_1d(ktf + 1) = z_at_walls_1d(ktf) + 1.0 / rdzw(i, ktf, j)

          l0 = alpha * l0_num / l0_den

          zi = c_g * l0

          do k = kts + 1, ktf
            N = Sqrt (G_OVER_TREF * dthetav_dz(i, k, j))
            l_d = c_r * Sqrt (q_sq(i, k, j)) / N

            l_b = l0 * KARMAN * z_at_walls_1d(k) / (KARMAN * z_at_walls_1d(k) + l0)

            if (z_at_walls_1d(k) <= zi) then
              l_master(i, k, j) = Min (l_b, l_d)
            else
              l_del = c_s * (1.0 / rdzw(i, k, j))
              l_master(i, k, j) = Min (l_del, l_d, KARMAN * z_at_walls_1d(k))
            end if
            l_boulac(i, k, j) = l_master(i, k, j)  ! here we are just filling in dissipation boulac length scale
          end do

          l_master(i, kts, j) = 0.0
          l_boulac(i, kts, j) = 0.0
          l_master(i, kte, j) = l_master(i, ktf, j)
          l_boulac(i, kte, j) = l_boulac(i, ktf, j)

        end do
      end do

    end subroutine Calc_l_messinger_master_algebra


    subroutine Calc_l_boulac_master_algebra (q_sq, dz, th, dlu, dld, dls, dlk, l_master, l_boulac, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the master length scale following       !
        !          Bougeault and Lacarrere (1989) and                 ! 
        !          Martilli et al. (2002)                             !
        !                                                             !
        ! Author: Timothy W. Juliano                                  !
        !                                                             !
        ! Comments: Adapted from Boulac PBL scheme code               !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: th
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dlu, dld, dls, dlk
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: l_master, l_boulac
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k, kk, i_start, i_end, j_start, j_end, ktf
      real :: beta, bbb, zup, zdo, zdo_sup, zzz, zup_inf, dzt, te, tl
      real, dimension (its:ite, kts:kte, jts:jte) :: z

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate height agl
      do j = j_start, j_end
        do i = i_start, i_end
          z(i, kts, j) = 0.0
          do k = kts + 1, kte
            z(i, k, j) = z(i, k - 1, j) + dz(i, k - 1, j)
          end do
        end do
      end do

        ! Calculate up and down length scales according to BouLac
      do j = j_start, j_end
        do i = i_start, i_end
          do k = kts, kte
            zup = 0.
            dlu(i, k, j) = z(i, kte, j) - z(i, k, j) - dz(i, k, j)/2.
            zzz = 0.
            zup_inf = 0.
            beta = G_OVER_TREF      !Buoyancy coefficient
            te = Max (q_sq(i, k, j) / 2., Q_SQ_MIN)
            do kk = k, kte - 1
              dzt = (dz(i, kk + 1, j) + dz(i, kk, j)) / 2.
              zup = zup - beta * th(i, k, j) * dzt
              zup = zup + beta * (th(i, kk + 1, j) + th(i, kk, j)) * dzt / 2.
              zzz = zzz + dzt
              if (te .lt. zup .and. te .ge. zup_inf) then
                bbb = (th(i, kk + 1, j) - th(i, kk, j)) / dzt
                if (bbb .ne. 0) then
                  tl = (-beta * (th(i, kk, j) - th(i, k, j)) + sqrt( max( 0.0, (beta * (th(i, kk, j) - th(i, k, j)))**2. + &
                            2. * bbb * beta * (te - zup_inf)))) / bbb / beta
                else
                  if (th(i, kk, j) .ne. th(i, k, j)) then
                    tl = (te - zup_inf) / (beta * (th(i, kk, j) - th(i, k, j)))
                  else
                    tl = 0.
                  end if
                end if
                dlu(i, k, j) = zzz - dzt + tl
              end if
              zup_inf = zup
            end do

            zdo = 0.
            zdo_sup = 0.
            dld(i, k, j) = z(i, k, j) + dz(i, k, j) / 2.
            zzz = 0.
            do kk = k, kts + 1, -1
              dzt = (dz(i, kk - 1, j) + dz(i, kk, j)) / 2.
              zdo = zdo + beta * th(i, k, j) * dzt
              zdo = zdo - beta * (th(i, kk - 1, j) + th(i, kk, j)) * dzt / 2.
              zzz = zzz + dzt
              if (te .lt. zdo .and. te .ge. zdo_sup) then
                bbb = (th(i, kk, j) - th(i, kk - 1, j)) / dzt
                if (bbb .ne. 0.0) then
                  tl = (beta * (th(i, kk, j) - th(i, k, j)) + sqrt( max( 0.0, (beta * (th(i, kk, j) - th(i, k, j)))**2. + &
                           2. * bbb * beta * (te - zdo_sup)))) / bbb / beta
                else
                  if (th(i, kk, j) .ne. th(i, k, j)) then
                    tl = (te - zdo_sup) / (beta * (th(i, kk, j) - th(i, k, j)))
                  else
                    tl = 0.
                  end if
                end if
                dld(i, k, j) = zzz - dzt + tl
              end if
              zdo_sup = zdo
            end do
          end do

            ! Compute turbulent and dissipation length scales
          do k = kts, ktf
            dlu(i, k, j) = max( dlu(i, k, j), 0.0 )
            dld(i, k, j) = max( min( dld(i, k, j), z(i, k, j) ), 0.0 )
            dls(i, k, j) = max( sqrt( dlu(i, k, j) * dld(i, k, j) ), 0.1 )
            dlk(i, k, j) = min( dlu(i, k, j), dld(i, k, j) )
            l_master(i, k, j) = dlk(i, k, j)  ! turbulent length scale
            l_boulac(i, k, j) = dls(i, k, j)  ! dissipation length scale
          end do
          l_master(i, kte, j) = l_master(i, ktf, j)
          l_boulac(i, kte, j) = l_boulac(i, ktf, j)
        end do
      end do

    end subroutine Calc_l_boulac_master_algebra


    subroutine Calc_fluxes_pbl_approx (config_flags, ss, l_master, l_boulac, du_dz, dv_dz, dthetav_dz, dqv_dz, &
              dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz, &
              t, qv, dn, dnw, fnm, fnp, dx, dy, dz, z_at_w, zx, zy, msftx, msfty, rdz, rdzw, cf1, cf2, cf3, dt, id, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts,    &
              jte, kts, kte, q_sq, q_sq_prog, turb_flux_u2, turb_flux_v2, turb_flux_w2,     &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                 &
              turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv,    &
              turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, &
              turb_flux_uqc, turb_flux_vqc, turb_flux_wqc, turb_flux_uqr, turb_flux_vqr, turb_flux_wqr, &
              turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, turb_flux_uqs, turb_flux_vqs, turb_flux_wqs, &
              turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, turb_flux_uqnc, turb_flux_vqnc, turb_flux_wqnc, &
              turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_wqni, &
              turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa, &
              u, v, w, thetav, rho, ust, hfx, qfx, rmol)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Diagnose q_sq and turb fluxes                            !
      !                                                                   !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                     !
      !                                                                   !
      ! Comments: q_sq and fluxes are calculated at the faces             !
      !                                                                   !
      ! References: Y75                                                   !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags
      integer,                                     intent (in)    :: ss

      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, zx, zy, rdz, rdzw, z_at_w
      real,                                        intent (in)    :: cf1, cf2, cf3, dt, dx, dy
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz, dqv_dz, &
          dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master, l_boulac
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, &
          turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw,        &
          turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v,                  &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv,             &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,                        &
          turb_flux_uqc, turb_flux_vqc, turb_flux_wqc,                                 &
          turb_flux_uqr, turb_flux_vqr, turb_flux_wqr,                                 &
          turb_flux_uqi, turb_flux_vqi, turb_flux_wqi,                                 &
          turb_flux_uqs, turb_flux_vqs, turb_flux_wqs,                                 &
          turb_flux_uqg, turb_flux_vqg, turb_flux_wqg,                                 &
          turb_flux_uqnc, turb_flux_vqnc, turb_flux_wqnc,                              &
          turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr,                              &
          turb_flux_uqni, turb_flux_vqni, turb_flux_wqni,                              &
          turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa,                        &
          turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, t, thetav, rho, qv
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, ust, hfx, qfx, rmol
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: rif, sm_l2, sh_l2, sm, sh, q_ratio, q_sq_hl88
      real, dimension (its:ite, kts:kte, jts:jte) :: dlu, dld, dls, dlk
      integer :: i, j

        ! Calculate q_sq (level 2) or Fill q_sq with prognostic q_sq (level 2.5)
      if ( config_flags%pbl3d_prog .eq. 0 ) then   ! level 2 model
        call Prep_for_fluxes_l2_pbl_approx (config_flags, q_sq, l_master, l_boulac, rif, sm, sh, du_dz, dv_dz, &
            dthetav_dz, id, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else if ( config_flags%pbl3d_prog .gt. 0 ) then   ! level 2.5 model
        call Prep_for_fluxes_l2p5_pbl_approx (config_flags, q_sq, q_sq_prog, l_master, l_boulac, rif, sm_l2, sh_l2, sm, sh, &
            q_ratio, q_sq_hl88, du_dz, dv_dz, dthetav_dz, id, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Diagnose master length scale at faces
      if (config_flags%pbl3d_l_opt < 3) then
        call Calc_l_master_algebra (config_flags%pbl3d_l_opt, q_sq, dz, rdzw, rmol, turb_flux_wtheta_v, dthetav_dz, l_master, &
            l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else if (config_flags%pbl3d_l_opt == 3) then
        call Calc_l_boulac_master_algebra (q_sq, dz, t, dlu, dld, dls, dlk, l_master, l_boulac, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else
        call Calc_l_messinger_master_algebra (q_sq, dz, rdzw, rmol, turb_flux_wtheta_v, dthetav_dz, l_master, &
            l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Diagnose momentum, heat, moisture fluxes for pbl approx
      if ( config_flags%pbl3d_prog .eq. 0 ) then  ! For level 2 model, no need for realizability
        call Diagnose_fluxes_pbl_approx (do_scalar_mix=config_flags%pbl3d_scalar_mix, ss=ss, nsteps=config_flags%pbl3d_nsteps, &
            fnm=fnm, fnp=fnp, t=t, qv=qv, ust=ust, hfx=hfx, qfx=qfx, sm_tilde=sm, sh_tilde=sh, rif=rif, &
            q_sq=q_sq, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, &
            dqc_dz=dqc_dz, dqr_dz=dqr_dz, dqi_dz=dqi_dz, dqs_dz=dqs_dz, dqg_dz=dqg_dz, dqnc_dz=dqnc_dz, &
            dqnr_dz=dqnr_dz, dqni_dz=dqni_dz, dqnwfa_dz=dqnwfa_dz, dqnifa_dz=dqnifa_dz, u2=turb_flux_u2, &
            v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, &
            wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, &
            uqc=turb_flux_uqc, vqc=turb_flux_vqc, wqc=turb_flux_wqc, uqr=turb_flux_uqr, vqr=turb_flux_vqr, wqr=turb_flux_wqr, &
            uqi=turb_flux_uqi, vqi=turb_flux_vqi, wqi=turb_flux_wqi, uqs=turb_flux_uqs, vqs=turb_flux_vqs, wqs=turb_flux_wqs, &
            uqg=turb_flux_uqg, vqg=turb_flux_vqg, wqg=turb_flux_wqg, uqnc=turb_flux_uqnc, vqnc=turb_flux_vqnc, &
            wqnc=turb_flux_wqnc, uqnr=turb_flux_uqnr, vqnr=turb_flux_vqnr, wqnr=turb_flux_wqnr, uqni=turb_flux_uqni, &
            vqni=turb_flux_vqni, wqni=turb_flux_wqni, uqnwfa=turb_flux_uqnwfa, vqnwfa=turb_flux_vqnwfa, wqnwfa=turb_flux_wqnwfa, &
            uqnifa=turb_flux_uqnifa, vqnifa=turb_flux_vqnifa, wqnifa=turb_flux_wqnifa, uth=turb_flux_utheta, &
            vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      else if ( config_flags%pbl3d_prog .gt. 0 ) then  ! For level 2.5 model, include information for HL88 mod
        call Diagnose_fluxes_pbl_approx (do_scalar_mix=config_flags%pbl3d_scalar_mix, ss=ss, nsteps=config_flags%pbl3d_nsteps, &
            fnm=fnm, fnp=fnp, t=t, qv=qv, ust=ust, hfx=hfx, qfx=qfx, sm_tilde=sm, sh_tilde=sh, sm_l2=sm_l2, sh_l2=sh_l2, rif=rif, &
            q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, &
            dqc_dz=dqc_dz, dqr_dz=dqr_dz, dqi_dz=dqi_dz, dqs_dz=dqs_dz, dqg_dz=dqg_dz, dqnc_dz=dqnc_dz, &
            dqnr_dz=dqnr_dz, dqni_dz=dqni_dz, dqnwfa_dz=dqnwfa_dz, dqnifa_dz=dqnifa_dz, u2=turb_flux_u2, &
            v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, &
            wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, &
            uqc=turb_flux_uqc, vqc=turb_flux_vqc, wqc=turb_flux_wqc, uqr=turb_flux_uqr, vqr=turb_flux_vqr, wqr=turb_flux_wqr, &
            uqi=turb_flux_uqi, vqi=turb_flux_vqi, wqi=turb_flux_wqi, uqs=turb_flux_uqs, vqs=turb_flux_vqs, wqs=turb_flux_wqs, &
            uqg=turb_flux_uqg, vqg=turb_flux_vqg, wqg=turb_flux_wqg, uqnc=turb_flux_uqnc, vqnc=turb_flux_vqnc, &
            wqnc=turb_flux_wqnc, uqnr=turb_flux_uqnr, vqnr=turb_flux_vqnr, wqnr=turb_flux_wqnr, uqni=turb_flux_uqni, &
            vqni=turb_flux_vqni, wqni=turb_flux_wqni, uqnwfa=turb_flux_uqnwfa, vqnwfa=turb_flux_vqnwfa, wqnwfa=turb_flux_wqnwfa, &
            uqnifa=turb_flux_uqnifa, vqnifa=turb_flux_vqnifa, wqnifa=turb_flux_wqnifa, uth=turb_flux_utheta, &
            vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
      end if

      call Diagnose_fluxes_surface (config_flags%pbl3d_sfc_opt, u, v, t, rho, dz, qv, ust, hfx, qfx, &
          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

    end subroutine Calc_fluxes_pbl_approx


    subroutine Prep_for_fluxes_l2_pbl_approx (config_flags, q_sq, l_master, l_boulac, rif, sm, sh, du_dz, dv_dz, &
              dthetav_dz, id, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Calculate stability functions for q_sq diagnosis          !
      !                                                                    !
      ! Author: Timothy W. Juliano                                         !
      !                                                                    !
      ! Comments: Stability functions and q_sq are calculated at the faces !
      !           The q_sq equation is for the PBL approximation           !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: rif, sm, sh
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: l_dissip
      character (4) :: my_level

      my_level = 'l2'
        ! Calc stability functions for level 2 model
      call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
         ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
         its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)
      call Fill_dissip_length_scale (config_flags%pbl3d_l_opt, l_master, l_boulac, l_dissip, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call Calc_q_sq_l2_pbl_approx (l_master=l_master, l_dissip=l_dissip, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
          jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
          jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, sm=sm, rif=rif)
        ! Relax q_sq at the boundaries (only for parent domain and specified B.C. runs)
      if ( id .eq. 1 .and. config_flags%specified) then
        call Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

    end subroutine Prep_for_fluxes_l2_pbl_approx


    subroutine Prep_for_fluxes_l2p5_pbl_approx (config_flags, q_sq, q_sq_prog, l_master, l_boulac, rif, sm_l2, sh_l2, sm, sh, &
              q_ratio, q_sq_hl88, du_dz, dv_dz, dthetav_dz, id, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                           !
      ! Purpose: Calculate stability functions for q_sq diagnosis (level 2 model) !
      !          so that we can apply HL88 realizability criterion                !
      !          Also calculate stability functions for level 2.5 model           !
      !                                                                           !
      ! Author: Timothy W. Juliano                                                !
      !                                                                           !
      ! Comments: Stability functions and q_sq are calculated at the faces        !
      !           The q_sq equation is for the PBL approximation                  !
      !                                                                           !
      ! References: Helfand and Labraga (1988, JAS)                               !
      !                                                                           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real,                                        intent (in)    :: cf1, cf2, cf3
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: rif, sm_l2, sh_l2, sm, sh, q_ratio, q_sq_hl88
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: l_dissip
      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level

      if (use_hl88) then   ! calc sm and sh then q_sq for level 2 model to prepare for HL88 fix
        my_level = 'l2'
        call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm_l2, sh=sh_l2, rif=rif)
        call Fill_dissip_length_scale (config_flags%pbl3d_l_opt, l_master, l_boulac, l_dissip, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        call Calc_q_sq_l2_pbl_approx (l_master=l_master, l_dissip=l_dissip, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
            jts=jts, jte=jte, kts=kts, kte=kte, q_sq_hl88=q_sq_hl88, sm=sm_l2, rif=rif)
      end if
      call Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      my_level = 'l2p5'
        ! Calc stability functions for level 2.5 model
      call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
        jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, &
        kts=kts, kte=kte, sm_l2=sm_l2, sh_l2=sh_l2, sm=sm, sh=sh, q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master_at_face=l_master)

    end subroutine Prep_for_fluxes_l2p5_pbl_approx


    subroutine Calc_fluxes (config_flags, id, specified, q_sq, q_sq_prog, l_master, l_boulac, du_dx, du_dy, du_dz, dv_dx, dv_dy, &
        dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
        turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, &
        turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, mat_cond_heat, mat_cond_moist, &
        rpgf_heat, rpgf_moist, t, thetav, qv, rho, hfx, qfx, rmol, u, v, ust, dz, rdzw, z_at_w, fnm, fnp, cf1, cf2, cf3, dn, dnw, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                !
      ! Purpose: Diagnoses turb fluxes                                 !
      !                                                                !
      ! Author: Pedro A. Jimenez                                       !
      !                                                                !
      ! Comments: The fluxes are calculated at the faces.              !
      !                                                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      logical,                                     intent (in)    :: specified
      integer,                                     intent (in)    :: id

      integer,                                     intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                                     dw_dx, du_dz, dw_dy, dv_dz, dthetav_dx, &
                                                                     dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, &
                                                                     dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_utheta_v, turb_flux_vtheta_v,   &
                                                                     turb_flux_wtheta_v, turb_flux_theta2_v,   &
                                                                     turb_flux_utheta, turb_flux_vtheta,       &
                                                                     turb_flux_wtheta, turb_flux_uqv,          &
                                                                     turb_flux_vqv, turb_flux_wqv,             &
                                                                     q_sq, l_master, l_boulac
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: mat_cond_heat, mat_cond_moist, rpgf_heat, rpgf_moist
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: u, v, rho, t, thetav, qv
      real, dimension (ims:ime, jms:jme),          intent (in)    :: ust, hfx, qfx, rmol
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: dz, rdzw, z_at_w
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp, dn, dnw
      real,                                        intent (in)    :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf, iter
      real, dimension (its:ite, kts:kte, jts:jte) :: sm, sh, sm_l2, sh_l2, rif, q_sq_hl88, q_ratio
      real, dimension (its:ite, kts:kte, jts:jte) :: dlu, dld, dls, dlk

      real (kind = SOLVER_PREC) :: tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv

      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate q_sq (level 2) or Fill q_sq with prognostic q_sq (level 2.5)
        ! Full 3D level 2 model is unstable right now
!      if ( config_flags%pbl3d_prog .eq. 0 ) then   ! level 2 model
!        call Prep_for_fluxes_l2 (config_flags, q_sq, l_master, l_boulac, rif, sm, sh, du_dx, dv_dy, dw_dz, &
!            du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
!            turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, id, ids, ide, jds, jde, kds, kde, &
!            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!      else if ( config_flags%pbl3d_prog .gt. 0 ) then   ! level 2.5 model
!        call Prep_for_fluxes_l2p5 (config_flags, q_sq, q_sq_prog, l_master, l_boulac, rif, sm_l2, sh_l2, sm, sh, &
!            q_ratio, q_sq_hl88, du_dx, dv_dy, dw_dz, du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz, &
!            turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, &
!            id, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
!            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!      end if

        ! Calling level 2 model PBL approx for now
      if ( config_flags%pbl3d_prog .eq. 0 ) then   ! level 2 model
        call Prep_for_fluxes_l2_pbl_approx (config_flags, q_sq, l_master, l_boulac, rif, sm, sh, du_dz, dv_dz, &
            dthetav_dz, id, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else if ( config_flags%pbl3d_prog .gt. 0 ) then   ! level 2.5 model
        call Prep_for_fluxes_l2p5_pbl_approx (config_flags, q_sq, q_sq_prog, l_master, l_boulac, rif, sm_l2, sh_l2, sm, sh, &
            q_ratio, q_sq_hl88, du_dz, dv_dz, dthetav_dz, id, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Diagnose master length scale at faces
      if (config_flags%pbl3d_l_opt < 3) then
        call Calc_l_master_algebra (config_flags%pbl3d_l_opt, q_sq, dz, rdzw, rmol, turb_flux_wtheta_v, dthetav_dz, l_master, &
            l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else if (config_flags%pbl3d_l_opt == 3) then
        call Calc_l_boulac_master_algebra (q_sq, dz, t, dlu, dld, dls, dlk, l_master, l_boulac, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else
        call Calc_l_messinger_master_algebra (q_sq, dz, rdzw, rmol, turb_flux_wtheta_v, dthetav_dz, l_master, &
            l_boulac, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      mat_cond_heat(its:ite,kts:kte,jts:jte) = -9999.
      mat_cond_moist(its:ite,kts:kte,jts:jte) = -9999.

      loop_j: do j = j_start, j_end
        loop_k: do k = kts + 1, ktf
          loop_i: do i = i_start, i_end

            if_no_turb: if (q_sq(i, k, j) == Q_SQ_MIN) then

              call Set_fluxes_to_min (q_sq(i, k, j), turb_flux_u2(i, k, j),   &
                 turb_flux_v2(i, k, j), turb_flux_w2(i, k, j), turb_flux_uv(i, k, j), &
                 turb_flux_uw(i, k, j), turb_flux_vw(i, k, j), turb_flux_utheta_v(i, k, j), &
                 turb_flux_vtheta_v(i, k, j), turb_flux_wtheta_v(i, k, j), turb_flux_theta2_v(i, k, j), &
                 turb_flux_utheta(i, k, j), turb_flux_vtheta(i, k, j), turb_flux_wtheta(i, k, j),  &
                 turb_flux_uqv(i, k, j), turb_flux_vqv(i, k, j), turb_flux_wqv(i, k, j) )

            else

              call Diagnose_fluxes (q_sq(i, k, j), l_master(i, k, j), du_dx(i, k, j), du_dy(i, k, j),                                &
                  du_dz(i, k, j), dv_dx(i, k, j), dv_dy(i, k, j), dv_dz(i, k, j), dw_dx(i, k, j), dw_dy(i, k, j), dw_dz(i, k, j),    &
                  dthetav_dx(i, k, j), dthetav_dy(i, k, j), dthetav_dz(i, k, j), dqv_dx(i, k, j), dqv_dy(i, k, j), dqv_dz(i, k, j),  &
                  turb_flux_u2(i, k, j), turb_flux_v2(i, k, j), turb_flux_w2(i, k, j), turb_flux_uv(i, k, j), turb_flux_uw(i, k, j), &
                  turb_flux_vw(i, k, j), turb_flux_utheta_v(i, k, j), turb_flux_vtheta_v(i, k, j), turb_flux_wtheta_v(i, k, j),      &
                  turb_flux_theta2_v(i, k, j), turb_flux_uqv(i, k, j), turb_flux_vqv(i, k, j), turb_flux_wqv(i, k, j), tf_utv,       &
                  tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, mat_cond_heat(i, k, j), mat_cond_moist(i, k, j), rpgf_heat(i, k, j),       &
                  rpgf_moist(i, k, j), config_flags%pbl3d_solver, config_flags%pbl3d_precond, config_flags%pbl3d_out_iter,           &
                  config_flags%pbl3d_in_iter, config_flags%pbl3d_tol_abs, config_flags%pbl3d_tol_rel                                 )

                ! Calc heat fluxes
              call Calc_heat_flux (t, qv, fnm, fnp, i, j, k, tf_utv, tf_vtv, tf_wtv,  tf_uqv, tf_vqv, tf_wqv, &
                  turb_flux_utheta(i, k, j), turb_flux_vtheta(i, k, j), turb_flux_wtheta(i, k, j), ims, ime, jms, jme, kms, kme)

            end if if_no_turb

          end do loop_i
        end do loop_k
      end do loop_j

      call Diagnose_fluxes_surface (config_flags%pbl3d_sfc_opt, u, v, t, rho, dz, qv, ust, hfx, qfx, &
          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

    end subroutine Calc_fluxes


    subroutine Prep_for_fluxes_l2 (config_flags, q_sq, l_master, l_boulac, rif, sm, sh, du_dx, dv_dy, dw_dz, &
              du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_wtv, &
              id, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Calculate stability functions for q_sq diagnosis          !
      !                                                                    !
      ! Author: Timothy W. Juliano                                         !
      !                                                                    !
      ! Comments: Stability functions and q_sq are calculated at the faces !
      !           The q_sq equation is for the full 3D model               !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx, dv_dy, dw_dz, &
                                    du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: tf_u2, tf_v2, tf_w2, &
                                                            tf_uv, tf_uw, tf_vw, tf_wtv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: rif, sm, sh
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: l_dissip
      character (4) :: my_level

      my_level = 'l2'
        ! Calc stability functions for level 2 model
      call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
         ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
         its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)
      call Fill_dissip_length_scale (config_flags%pbl3d_l_opt, l_master, l_boulac, l_dissip, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      call Calc_q_sq_l2 (l_master=l_master, l_dissip=l_dissip, du_dx=du_dx, dv_dy=dv_dy, dw_dz=dw_dz, du_dy=du_dy, &
          dv_dx=dv_dx, dw_dx=dw_dx, du_dz=du_dz, dw_dy=dw_dy, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
          tf_u2=tf_u2, tf_v2=tf_v2, tf_w2=tf_w2, tf_uv=tf_uv, tf_uw=tf_uw, tf_vw=tf_vw, tf_wtv=tf_wtv, &
          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, rif=rif)
        ! Relax q_sq at the boundaries (only for parent domain and specified B.C. runs)
      if ( id .eq. 1 .and. config_flags%specified) then
        call Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

    end subroutine Prep_for_fluxes_l2


    subroutine Prep_for_fluxes_l2p5 (config_flags, q_sq, q_sq_prog, l_master, l_boulac, rif, sm_l2, sh_l2, sm, sh, &
              q_ratio, q_sq_hl88, du_dx, dv_dy, dw_dz, du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz, &
              tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_wtv, id, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                           !
      ! Purpose: Calculate stability functions for q_sq diagnosis (level 2 model) !
      !          so that we can apply HL88 realizability criterion                !
      !          Also calculate stability functions for level 2.5 model           !
      !                                                                           !
      ! Author: Timothy W. Juliano                                                !
      !                                                                           !
      ! Comments: Stability functions and q_sq are calculated at the faces        !
      !           The q_sq equation is for the full 3D model                      !
      !                                                                           !
      ! References: Helfand and Labraga (1988, JAS)                               !
      !                                                                           !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)    :: config_flags

      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real,                                        intent (in)    :: cf1, cf2, cf3
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx, dv_dy, dw_dz, &
                                    du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: tf_u2, tf_v2, tf_w2, &
                                                            tf_uv, tf_uw, tf_vw, tf_wtv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: rif, sm_l2, sh_l2, sm, sh, q_ratio, q_sq_hl88
      integer, intent (in)    :: id
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: l_dissip
      logical, parameter :: use_hl88 = .true.
      character (4) :: my_level

      if (use_hl88) then   ! calc sm and sh then q_sq for level 2 model to prepare for HL88 fix
        my_level = 'l2'
        call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm_l2, sh=sh_l2, rif=rif)
        call Fill_dissip_length_scale (config_flags%pbl3d_l_opt, l_master, l_boulac, l_dissip, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        call Calc_q_sq_l2 (l_master=l_master, l_dissip=l_dissip, du_dx=du_dx, dv_dy=dv_dy, dw_dz=dw_dz, du_dy=du_dy, &
            dv_dx=dv_dx, dw_dx=dw_dx, du_dz=du_dz, dw_dy=dw_dy, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
            tf_u2=tf_u2, tf_v2=tf_v2, tf_w2=tf_w2, tf_uv=tf_uv, tf_uw=tf_uw, tf_vw=tf_vw, tf_wtv=tf_wtv, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, q_sq_hl88=q_sq_hl88, rif=rif)
      end if
      call Fill_q_sq_with_q_sq_prog (q_sq_prog, q_sq, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      my_level = 'l2p5'
        ! Calc stability functions for level 2.5 model
      call Calc_sm_sh_l2_or_l2p5 (my_level=my_level, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
        jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, &
        kts=kts, kte=kte, sm_l2=sm_l2, sh_l2=sh_l2, sm=sm, sh=sh, q_sq=q_sq, q_sq_hl88=q_sq_hl88, q_ratio=q_ratio, l_master_at_face=l_master)

    end subroutine Prep_for_fluxes_l2p5


    subroutine Fill_dissip_length_scale (l_opt, l_master, l_boulac, l_dissip, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      integer,                                     intent (in)  :: l_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: l_dissip

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end
      real :: l
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            if (l_opt < 3) then
              l_dissip(i, k, j) = l_master(i, k, j)
            else
              l_dissip(i, k, j) = l_boulac(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Fill_dissip_length_scale


    subroutine Set_fluxes_to_min (q_sq, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, &
        tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, tf_ut, tf_vt, tf_wt, tf_uqv, &
        tf_vqv, tf_wqv)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Set the turbulent fluxes to the minimum val  !
        !                                                       !
        ! Author: Pedro A. Jimenez                              !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, intent (inout) :: q_sq, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw,  &
          tf_utv, tf_vtv, tf_wtv, tf_t2v, tf_ut, tf_vt, tf_wt, tf_uqv, tf_vqv, &
          tf_wqv

      q_sq = Q_SQ_MIN
      tf_u2 = TURB_FLUX_MIN
      tf_v2 = TURB_FLUX_MIN
      tf_w2 = TURB_FLUX_MIN
      tf_uv = TURB_FLUX_MIN
      tf_uw = TURB_FLUX_MIN
      tf_vw = TURB_FLUX_MIN
      tf_utv = TURB_FLUX_MIN
      tf_vtv = TURB_FLUX_MIN
      tf_wtv = TURB_FLUX_MIN
      tf_t2v = TURB_FLUX_MIN
      tf_ut = TURB_FLUX_MIN
      tf_vt = TURB_FLUX_MIN
      tf_wt = TURB_FLUX_MIN
      tf_uqv = TURB_FLUX_MIN
      tf_vqv = TURB_FLUX_MIN
      tf_wqv = TURB_FLUX_MIN

    end subroutine Set_fluxes_to_min


    subroutine Diagnose_fluxes (q_sq_at_walls_1dpbl, l_master, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz,      &
        dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,        &
        turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv, &
        turb_flux_vqv, turb_flux_wqv, tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, mat_cond_heat, mat_cond_moist, rpgf_heat,    &
        rpgf_moist, pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter, pbl3d_tol_abs, pbl3d_tol_rel)

      implicit none

      real, intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                              dw_dx, du_dz, dw_dy, dv_dz,             &
                              dthetav_dx, dthetav_dy, dthetav_dz,     &
                              dqv_dx, dqv_dy, dqv_dz
      real, intent (in)    :: q_sq_at_walls_1dpbl, l_master
      real, intent (inout) :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                              turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                              turb_flux_utheta_v, turb_flux_vtheta_v,   &
                              turb_flux_wtheta_v, turb_flux_theta2_v,   &
                              turb_flux_uqv, turb_flux_vqv, turb_flux_wqv
      real, intent (out)   :: mat_cond_heat, mat_cond_moist, rpgf_heat, rpgf_moist
      real (kind = SOLVER_PREC), intent (out) :: tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv

!!! FOR INTERATIVE SOLVER
      integer, intent(in)                       :: pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter
      real, intent(in)                          :: pbl3d_tol_abs, pbl3d_tol_rel

        ! Local vars
      real (kind = SOLVER_PREC) :: q, q_sq, l, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_t2v
      real (kind = SOLVER_PREC) :: du_dx_tmp, dv_dy_tmp, dw_dz_tmp, du_dy_tmp, dv_dx_tmp,      &
                                   dw_dx_tmp, du_dz_tmp, dw_dy_tmp, dv_dz_tmp, dthetav_dx_tmp, &
                                   dthetav_dy_tmp, dthetav_dz_tmp, dqv_dx_tmp, dqv_dy_tmp, dqv_dz_tmp

      logical, parameter :: DEBUG = .false.

        ! Convert vars from single prec to double prec
      q_sq = Working_to_solver_prec(q_sq_at_walls_1dpbl)
      q = Sqrt ( q_sq )
      l = Working_to_solver_prec(l_master)
      du_dx_tmp = Working_to_solver_prec(du_dx)
      du_dy_tmp = Working_to_solver_prec(du_dy)
      du_dz_tmp = Working_to_solver_prec(du_dz)
      dv_dx_tmp = Working_to_solver_prec(dv_dx)
      dv_dy_tmp = Working_to_solver_prec(dv_dy)
      dv_dz_tmp = Working_to_solver_prec(dv_dz)
      dw_dx_tmp = Working_to_solver_prec(dw_dx)
      dw_dy_tmp = Working_to_solver_prec(dw_dy)
      dw_dz_tmp = Working_to_solver_prec(dw_dz)
      dthetav_dx_tmp = Working_to_solver_prec(dthetav_dx)
      dthetav_dy_tmp = Working_to_solver_prec(dthetav_dy)
      dthetav_dz_tmp = Working_to_solver_prec(dthetav_dz)
      dqv_dx_tmp = Working_to_solver_prec(dqv_dx)
      dqv_dy_tmp = Working_to_solver_prec(dqv_dy)
      dqv_dz_tmp = Working_to_solver_prec(dqv_dz)

        ! Solve system of equations for momentum fluxes
      call Solve_turb_system (q, q_sq, l, du_dx_tmp, du_dy_tmp, du_dz_tmp, dv_dx_tmp, dv_dy_tmp, dv_dz_tmp, dw_dx_tmp, dw_dy_tmp, dw_dz_tmp, &
          dthetav_dx_tmp, dthetav_dy_tmp, dthetav_dz_tmp, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, &
          turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, &
          tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, mat_cond_heat, rpgf_heat, &
          pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter, pbl3d_tol_abs, pbl3d_tol_rel)

        ! Solve system of equations for moisture fluxes
      call Solve_turb_system_moist (q, l, du_dx_tmp, du_dy_tmp, du_dz_tmp, dv_dx_tmp, dv_dy_tmp, dv_dz_tmp, dw_dx_tmp, dw_dy_tmp, dw_dz_tmp, dthetav_dx_tmp, &
          dthetav_dy_tmp, dthetav_dz_tmp, dqv_dx_tmp, dqv_dy_tmp, dqv_dz_tmp, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, &
          tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, mat_cond_moist, rpgf_moist, pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter, pbl3d_tol_abs, pbl3d_tol_rel)

        ! Convert vars from double prec to single prec
      turb_flux_u2 = Solver_to_working_prec(tf_u2)
      turb_flux_v2 = Solver_to_working_prec(tf_v2)
      turb_flux_w2 = Solver_to_working_prec(tf_w2)
      turb_flux_uv = Solver_to_working_prec(tf_uv)
      turb_flux_uw = Solver_to_working_prec(tf_uw)
      turb_flux_vw = Solver_to_working_prec(tf_vw)
      turb_flux_utheta_v = Solver_to_working_prec(tf_utv)
      turb_flux_vtheta_v = Solver_to_working_prec(tf_vtv)
      turb_flux_wtheta_v = Solver_to_working_prec(tf_wtv)
      turb_flux_theta2_v = Solver_to_working_prec(tf_t2v)
      turb_flux_uqv = Solver_to_working_prec(tf_uqv)
      turb_flux_vqv = Solver_to_working_prec(tf_vqv)
      turb_flux_wqv = Solver_to_working_prec(tf_wqv)

    end subroutine Diagnose_fluxes


    pure function Working_to_solver_prec (var)

      implicit none

      real (kind = WORKING_PREC), intent(in) :: var
      real (kind = SOLVER_PREC)              :: Working_to_solver_prec

      working_to_solver_prec = Real (var, kind = SOLVER_PREC)

    end function Working_to_solver_prec


    pure function Solver_to_working_prec (var)

      implicit none

      real (kind = SOLVER_PREC), intent(in) :: var
      real (kind = WORKING_PREC)            :: Solver_to_working_prec

      solver_to_working_prec = Real (var, kind = WORKING_PREC)

    end function Solver_to_working_prec


    subroutine Solve_turb_system (q, q_sq, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
              dthetav_dx, dthetav_dy, dthetav_dz, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, &
              turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, &
              tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v, mat_cond_heat, rpgf_heat, &
              pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter, pbl3d_tol_abs, pbl3d_tol_rel)

      implicit none

      real (kind = SOLVER_PREC), intent (in)    :: q, q_sq, l
      real (kind = SOLVER_PREC), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                   dw_dx, du_dz, dw_dy, dv_dz,             &
                                                   dthetav_dx, dthetav_dy, dthetav_dz
      real,                      intent (in)    :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                   turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                   turb_flux_utheta_v, turb_flux_vtheta_v,   &
                                                   turb_flux_wtheta_v, turb_flux_theta2_v
      real,                      intent (out)   :: mat_cond_heat, rpgf_heat
      real (kind = SOLVER_PREC), intent (out)   :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v

      integer, intent(in)                       :: pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter
      real, intent(in)                          :: pbl3d_tol_abs, pbl3d_tol_rel

        ! Local vars
      character (len = 1) :: equed
      integer, parameter   :: N_VARS = 10
      integer :: info, ipiv(N_VARS)
      real (kind = SOLVER_PREC) :: rcond, cond, anorm, ferr, berr
      real (kind = SOLVER_PREC), dimension(4*N_VARS) :: work
      real (kind = SOLVER_PREC) :: swork(N_VARS*(N_VARS+N_VARS))
      real (kind = SOLVER_PREC), dimension(N_VARS) :: iwork, rsf, csf
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: b(:), x(:)
      real (kind = SOLVER_PREC), allocatable, dimension(:,:) :: a(:, :), af(:, :)

!!! FOR INTERATIVE SOLVER
      integer, parameter   :: N_NZ = 61
      integer, parameter   :: ia(N_NZ) = (/1,1,1,1,1,1,1,&
                                           2,2,2,2,2,2,2,&
                                           3,3,3,3,3,3,3,&
                                           4,4,4,4,4,&
                                           5,5,5,5,5,5,&
                                           6,6,6,6,6,6,&
                                           7,7,7,7,7,7,&
                                           8,8,8,8,8,8,&
                                           9,9,9,9,9,9,9,&
                                           10,10,10,10/)
      integer, parameter   :: ia_comp(N_VARS+1) = (/1,8,15,22,27,33,39,45,51,58,62/)
      integer, parameter   :: ja(N_NZ) = (/1,2,3,4,5,6,9,&
                                           1,2,3,4,5,6,9,&
                                           1,2,3,4,5,6,9,&
                                           1,2,4,5,6,&
                                           1,3,4,5,6,7,&
                                           2,3,4,5,6,8,&
                                           1,4,5,7,8,9,&
                                           2,4,6,7,8,9,&
                                           3,5,6,7,8,9,10,&
                                           7,8,9,10/)
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: a_nz(:), b_fg(:)
      real (kind = SOLVER_PREC) :: tf_u2_tmp, tf_v2_tmp, tf_w2_tmp, tf_uv_tmp, tf_uw_tmp, tf_vw_tmp, &
                                   tf_utv_tmp, tf_vtv_tmp, tf_wtv_tmp, tf_t2v_tmp
      real (kind = SOLVER_PREC) :: tol_abs, tol_rel
!!! FOR INTERATIVE SOLVER

      external dgesv
      external dgesvx

      allocate (a(N_VARS, N_VARS), af(N_VARS, N_VARS), b(N_VARS), x(N_VARS))
      allocate (a_nz(N_NZ), b_fg(N_VARS))

        ! Fills in the B matrix
      call Fill_in_b_vector (b, q, q_sq, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)

      if (pbl3d_solver == 1) then  ! Direct method
          ! Fills in the A matrix
        call Fill_in_a_matrix (a, q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            dthetav_dx, dthetav_dy, dthetav_dz)

          ! Compute solution to system of linear equations
        call dgesvx ('N', 'N', N_VARS, 1, a, N_VARS, af, N_VARS, ipiv, equed, rsf, csf, b, N_VARS, x, N_VARS, rcond, ferr, berr, work, iwork, info)
        mat_cond_heat = 1. / rcond
        rpgf_heat = work(1)

        tf_u2 = x(1)
        tf_v2 = x(2)
        tf_w2 = x(3)
        tf_uv = x(4)
        tf_uw = x(5)
        tf_vw = x(6)
        tf_utv = x(7)
        tf_vtv = x(8)
        tf_wtv = x(9)
        tf_t2v = x(10)

      else  ! Iterative method

          ! Convert vars from single prec to double prec
        tf_u2_tmp = Working_to_solver_prec(turb_flux_u2)
        tf_v2_tmp = Working_to_solver_prec(turb_flux_v2)
        tf_w2_tmp = Working_to_solver_prec(turb_flux_w2)
        tf_uv_tmp = Working_to_solver_prec(turb_flux_uv)
        tf_uw_tmp = Working_to_solver_prec(turb_flux_uw)
        tf_vw_tmp = Working_to_solver_prec(turb_flux_vw)
        tf_utv_tmp = Working_to_solver_prec(turb_flux_utheta_v)
        tf_vtv_tmp = Working_to_solver_prec(turb_flux_vtheta_v)
        tf_wtv_tmp = Working_to_solver_prec(turb_flux_wtheta_v)
        tf_t2v_tmp = Working_to_solver_prec(turb_flux_theta2_v)

          ! Fills in the first guess RHS
        call Fill_in_b_first_guess (b_fg, tf_u2_tmp, tf_v2_tmp, tf_w2_tmp, tf_uv_tmp, tf_uw_tmp, tf_vw_tmp, tf_utv_tmp, tf_vtv_tmp, tf_wtv_tmp, tf_t2v_tmp)

          ! Fills in the A matrix, non-zero
        call Fill_in_a_matrix_non_zero (a_nz, q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            dthetav_dx, dthetav_dy, dthetav_dz)

          ! Call MGMRES
        tol_abs = Working_to_solver_prec(pbl3d_tol_abs)
        tol_rel = Working_to_solver_prec(pbl3d_tol_rel)

        if (pbl3d_precond == 0) then  ! No preconditioner
          call mgmres_st ( N_VARS, N_NZ, ia, ja, a_nz, b_fg, b, pbl3d_out_iter, pbl3d_in_iter, tol_abs, tol_rel )
        else
          call pmgmres_ilu_cr ( N_VARS, N_NZ, ia_comp, ja, a_nz, b_fg, b, pbl3d_out_iter, pbl3d_in_iter, tol_abs, tol_rel )
        end if

        tf_u2 = b_fg(1)
        tf_v2 = b_fg(2)
        tf_w2 = b_fg(3)
        tf_uv = b_fg(4)
        tf_uw = b_fg(5)
        tf_vw = b_fg(6)
        tf_utv = b_fg(7)
        tf_vtv = b_fg(8)
        tf_wtv = b_fg(9)
        tf_t2v = b_fg(10)
      end if

      deallocate (a, af, b, x)
      deallocate (a_nz, b_fg)

    end subroutine Solve_turb_system


    subroutine Solve_turb_system_moist (q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dx, &
              dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, &
              tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, &
              mat_cond_moist, rpgf_moist, pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter, pbl3d_tol_abs, pbl3d_tol_rel)

      implicit none

      real (kind = SOLVER_PREC), intent (in)    :: q, l
      real (kind = SOLVER_PREC), intent (in)    :: du_dx, dv_dy, dw_dz, du_dy, dv_dx,      &
                                                   dw_dx, du_dz, dw_dy, dv_dz,             &
                                                   dthetav_dx, dthetav_dy, dthetav_dz,     &
                                                   dqv_dx, dqv_dy, dqv_dz
      real,                      intent (in)    :: turb_flux_uqv, turb_flux_vqv, turb_flux_wqv
      real (kind = SOLVER_PREC), intent (in)    :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv
      real (kind = SOLVER_PREC), intent (out)   :: tf_uqv, tf_vqv, tf_wqv
      real,                      intent (out)   :: mat_cond_moist, rpgf_moist

      integer, intent(in)                       :: pbl3d_solver, pbl3d_precond, pbl3d_out_iter, pbl3d_in_iter
      real, intent(in)                          :: pbl3d_tol_abs, pbl3d_tol_rel

        ! Local vars
      character (len = 1) :: equed
      integer, parameter   :: N_VARS = 4
      integer :: info, ipiv(N_VARS)
      real (kind = SOLVER_PREC) :: rcond, cond, anorm, ferr, berr
      real (kind = SOLVER_PREC), dimension(4*N_VARS) :: work
      real (kind = SOLVER_PREC) :: swork(N_VARS*(N_VARS+N_VARS))
      real (kind = SOLVER_PREC), dimension(N_VARS) :: iwork, rsf, csf
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: b(:), x(:)
      real (kind = SOLVER_PREC), allocatable, dimension(:,:) :: a(:, :), af(:, :)

!!! FOR INTERATIVE SOLVER
      integer, parameter   :: N_NZ = 14
      integer, parameter   :: ia(N_NZ) = (/1,1,1,&
                                           2,2,2,&
                                           3,3,3,3,&
                                           4,4,4,4/)
      integer, parameter   :: ia_comp(N_VARS+1) = (/1,4,7,12,15/)
      integer, parameter   :: ja(N_NZ) = (/1,2,3,&
                                           1,2,3,&
                                           1,2,3,4,&
                                           1,2,3,4/)
      real (kind = SOLVER_PREC), allocatable, dimension(:) :: a_nz(:), b_fg(:)
      real (kind = SOLVER_PREC) :: tf_uqv_tmp, tf_vqv_tmp, tf_wqv_tmp
      real (kind = SOLVER_PREC) :: tol_abs, tol_rel
!!! FOR INTERATIVE SOLVER

      external dgesv
      external dgesvx

      allocate (a(N_VARS, N_VARS), af(N_VARS, N_VARS), b(N_VARS), x(N_VARS))
      allocate (a_nz(N_NZ), b_fg(N_VARS))

        ! Fills in the B matrix
      call Fill_in_b_vector_moist (b, dqv_dx, dqv_dy, dqv_dz, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv)

      if (pbl3d_solver == 1) then  ! Direct method
          ! Fills in the A matrix
        call Fill_in_a_matrix_moist (a, q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            dthetav_dx, dthetav_dy, dthetav_dz)

          ! Compute solution to system of linear equations
        call dgesvx ('N', 'N', N_VARS, 1, a, N_VARS, af, N_VARS, ipiv, equed, rsf, csf, b, N_VARS, x, N_VARS, rcond, ferr, berr, work, iwork, info)
        mat_cond_moist = 1. / rcond

        tf_uqv = x(1)
        tf_vqv = x(2)
        tf_wqv = x(3)

      else  ! Iterative method

          ! Convert vars from single prec to double prec
        tf_uqv_tmp = Working_to_solver_prec(turb_flux_uqv)
        tf_vqv_tmp = Working_to_solver_prec(turb_flux_vqv)
        tf_wqv_tmp = Working_to_solver_prec(turb_flux_wqv)

          ! Fills in the first guess RHS
        call Fill_in_b_first_guess_moist (b_fg, tf_uqv_tmp, tf_vqv_tmp, tf_wqv_tmp)

          ! Fills in the A matrix, non-zero
        call Fill_in_a_matrix_non_zero_moist (a_nz, q, l, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            dthetav_dx, dthetav_dy, dthetav_dz)

          ! Call MGMRES
        tol_abs = Working_to_solver_prec(pbl3d_tol_abs)
        tol_rel = Working_to_solver_prec(pbl3d_tol_rel)

        if (pbl3d_precond == 0) then  ! No preconditioner
          call mgmres_st ( N_VARS, N_NZ, ia, ja, a_nz, b_fg, b, pbl3d_out_iter, pbl3d_in_iter, tol_abs, tol_rel )
        else
          call pmgmres_ilu_cr ( N_VARS, N_NZ, ia_comp, ja, a_nz, b_fg, b, pbl3d_out_iter, pbl3d_in_iter, tol_abs, tol_rel )
        end if

        tf_uqv = b_fg(1)
        tf_vqv = b_fg(2)
        tf_wqv = b_fg(3)
      end if

      deallocate (a, af, b, x)
      deallocate (a_nz, b_fg)

    end subroutine Solve_turb_system_moist


    subroutine Fill_in_b_vector (b, q, q_sq, l, u_x, u_y, u_z, v_x, &
        v_y, v_z, w_x, w_y, w_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the vector b with the indenpendent terms         !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      real (kind = SOLVER_PREC), intent(in) :: q, q_sq, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b(:)

      b = 0.0
      b(1) = Real (q ** 3.0 / (6.0 * a_1 * l) + 3.0 * c_1 * q_sq * u_x, kind = SOLVER_PREC)
      b(2) = Real (q ** 3.0 / (6.0 * a_1 * l) + 3.0 * c_1 * q_sq * v_y, kind = SOLVER_PREC)
      b(3) = Real (q ** 3.0 / (6.0 * a_1 * l) + 3.0 * c_1 * q_sq * w_z, kind = SOLVER_PREC)
      b(4) = Real (c_1 * q_sq * (v_x + u_y), kind = SOLVER_PREC)
      b(5) = Real (c_1 * q_sq * (w_x + u_z), kind = SOLVER_PREC)
      b(6) = Real (c_1 * q_sq * (w_y + v_z), kind = SOLVER_PREC)

    end subroutine Fill_in_b_vector


    subroutine Fill_in_b_first_guess (b_fg, tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the vector b with the indenpendent terms         !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      real (kind = SOLVER_PREC), intent(in) :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv, tf_t2v
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b_fg(:)

      b_fg = 0.0
      b_fg(1) = Real (tf_u2, kind = SOLVER_PREC)
      b_fg(2) = Real (tf_v2, kind = SOLVER_PREC)
      b_fg(3) = Real (tf_w2, kind = SOLVER_PREC)
      b_fg(4) = Real (tf_uv, kind = SOLVER_PREC)
      b_fg(5) = Real (tf_uw, kind = SOLVER_PREC)
      b_fg(6) = Real (tf_vw, kind = SOLVER_PREC)
      b_fg(7) = Real (tf_utv, kind = SOLVER_PREC)
      b_fg(8) = Real (tf_vtv, kind = SOLVER_PREC)
      b_fg(9) = Real (tf_wtv, kind = SOLVER_PREC)
      b_fg(10) = Real (tf_t2v, kind = SOLVER_PREC)

    end subroutine Fill_in_b_first_guess


    subroutine Fill_in_b_vector_moist (b, qv_x, qv_y, qv_z, tf_u2, tf_v2, tf_w2, tf_uv, &
        tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv)

        implicit none

      real (kind = SOLVER_PREC), intent(in) :: qv_x, qv_y, qv_z
      real (kind = SOLVER_PREC), intent(in) :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_utv, tf_vtv, tf_wtv
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b(:)

      b = 0.0
      b(1) = Real (- (qv_x * tf_u2 + qv_y * tf_uv + qv_z * tf_uw), kind = SOLVER_PREC)
      b(2) = Real (- (qv_y * tf_v2 + qv_x * tf_uv + qv_z * tf_vw), kind = SOLVER_PREC)
      b(3) = Real (- (qv_z * tf_w2 + qv_x * tf_uw + qv_y * tf_vw), kind = SOLVER_PREC)
      b(4) = Real (- (qv_x * tf_utv + qv_y * tf_vtv + qv_z * tf_wtv), kind = SOLVER_PREC)

    end subroutine Fill_in_b_vector_moist


    subroutine Fill_in_b_first_guess_moist (b_fg, tf_uqv, tf_vqv, tf_wqv)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the vector b with the indenpendent terms         !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      implicit none

      real (kind = SOLVER_PREC), intent(in) :: tf_uqv, tf_vqv, tf_wqv
      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: b_fg(:)

      b_fg = 0.0
      b_fg(1) = Real (tf_uqv, kind = SOLVER_PREC)
      b_fg(2) = Real (tf_vqv, kind = SOLVER_PREC)
      b_fg(3) = Real (tf_wqv, kind = SOLVER_PREC)
      b_fg(4) = Real (0.0, kind = SOLVER_PREC)

    end subroutine Fill_in_b_first_guess_moist


    subroutine Fill_in_a_matrix (a, q, l, u_x, u_y, u_z, v_x, v_y, &
        v_z, w_x, w_y, w_z, tv_x, tv_y, tv_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the coefficients of the A Matrix                 !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~!!!!!!

      implicit none

      real (kind = SOLVER_PREC), dimension(:, :), intent(inout) :: a(:, :)
      real (kind = SOLVER_PREC), intent(in) :: q, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, &
          tv_x, tv_y, tv_z

      a = 0.0

         ! u2 equation
      a(1, 1) = Real (q / (2.0 * a_1 * l) + 2.0 * u_x, kind = SOLVER_PREC)
      a(1, 2) = Real (- v_y, kind = SOLVER_PREC)
      a(1, 3) = Real (- w_z, kind = SOLVER_PREC)
      a(1, 4) = Real (2.0 * u_y - v_x, kind = SOLVER_PREC)
      a(1, 5) = Real (2.0 * u_z - w_x, kind = SOLVER_PREC)
      a(1, 6) = Real (- v_z - w_y, kind = SOLVER_PREC)
      a(1, 9) = Real (G_OVER_TREF, kind = SOLVER_PREC)

         ! v2 equation
      a(2, 1) = Real (- u_x, kind = SOLVER_PREC)
      a(2, 2) = Real (q / (2.0 * a_1 * l) + 2.0 * v_y, kind = SOLVER_PREC)
      a(2, 3) = Real (- w_z, kind = SOLVER_PREC)
      a(2, 4) = Real (2.0 * v_x - u_y, kind = SOLVER_PREC)
      a(2, 5) = Real (- u_z - w_x, kind = SOLVER_PREC)
      a(2, 6) = Real (2.0 * v_z - w_y, kind = SOLVER_PREC)
      a(2, 9) = Real (G_OVER_TREF, kind = SOLVER_PREC)

         ! w2 equation
      a(3, 1) = Real (- u_x, kind = SOLVER_PREC)
      a(3, 2) = Real (- v_y, kind = SOLVER_PREC)
      a(3, 3) = Real (q / (2.0 * a_1 * l) + 2.0 * w_z, kind = SOLVER_PREC)
      a(3, 4) = Real (- u_y - v_x, kind = SOLVER_PREC)
      a(3, 5) = Real (2.0 * w_x - u_z, kind = SOLVER_PREC)
      a(3, 6) = Real (2.0 * w_y - v_z, kind = SOLVER_PREC)
      a(3, 9) = Real (- 2.0 * G_OVER_TREF, kind = SOLVER_PREC)

         ! uv equation
      a(4, 1) = Real (v_x, kind = SOLVER_PREC)
      a(4, 2) = Real (u_y, kind = SOLVER_PREC)
      a(4, 4) = Real (q / (3.0 * a_1 * l) + u_x + v_y, kind = SOLVER_PREC)
      a(4, 5) = Real (v_z, kind = SOLVER_PREC)
      a(4, 6) = Real (u_z, kind = SOLVER_PREC)

       ! uw equation
      a(5, 1) = Real (w_x, kind = SOLVER_PREC)
      a(5, 3) = Real (u_z, kind = SOLVER_PREC)
      a(5, 4) = Real (w_y, kind = SOLVER_PREC)
      a(5, 5) = Real (q / (3.0 * a_1 * l) + w_z + u_x, kind = SOLVER_PREC)
      a(5, 6) = Real (u_y, kind = SOLVER_PREC)
      a(5, 7) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! vw equation
      a(6, 2) = Real (w_y, kind = SOLVER_PREC)
      a(6, 3) = Real (v_z, kind = SOLVER_PREC)
      a(6, 4) = Real (w_x, kind = SOLVER_PREC)
      a(6, 5) = Real (v_x, kind = SOLVER_PREC)
      a(6, 6) = Real (q / (3.0 * a_1 * l) + w_z + v_y, kind = SOLVER_PREC)
      a(6, 8) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! ut_v equation
      a(7, 1) = Real (tv_x, kind = SOLVER_PREC)
      a(7, 4) = Real (tv_y, kind = SOLVER_PREC)
      a(7, 5) = Real (tv_z, kind = SOLVER_PREC)
      a(7, 7) = Real (q / (3.0 * a_2 * l) + u_x, kind = SOLVER_PREC)
      a(7, 8) = Real (u_y, kind = SOLVER_PREC)
      a(7, 9) = Real (u_z, kind = SOLVER_PREC)

         ! vt_v equation
      a(8, 2) = Real (tv_y, kind = SOLVER_PREC)
      a(8, 4) = Real (tv_x, kind = SOLVER_PREC)
      a(8, 6) = Real (tv_z, kind = SOLVER_PREC)
      a(8, 7) = Real (v_x, kind = SOLVER_PREC)
      a(8, 8) = Real (q / (3.0 * a_2 * l) + v_y, kind = SOLVER_PREC)
      a(8, 9) = Real (v_z, kind = SOLVER_PREC)

         ! wt_v equation
      a(9, 3) = Real (tv_z, kind = SOLVER_PREC)
      a(9, 5) = Real (tv_x, kind = SOLVER_PREC)
      a(9, 6) = Real (tv_y, kind = SOLVER_PREC)
      a(9, 7) = Real (w_x, kind = SOLVER_PREC)
      a(9, 8) = Real (w_y, kind = SOLVER_PREC)
      a(9, 9) = Real (q / (3.0 * a_2 * l) + w_z, kind = SOLVER_PREC)
      a(9, 10) = Real (- G_OVER_TREF, kind = SOLVER_PREC)

         ! tt_v equation
      a(10, 7) = Real (tv_x, kind = SOLVER_PREC)
      a(10, 8) = Real (tv_y, kind = SOLVER_PREC)
      a(10, 9) = Real (tv_z, kind = SOLVER_PREC)
      a(10, 10) = Real (q / (b_2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix


    subroutine Fill_in_a_matrix_non_zero (a_nz, q, l, u_x, u_y, u_z, v_x, v_y, &
        v_z, w_x, w_y, w_z, tv_x, tv_y, tv_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the coefficients of the A Matrix                 !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~!!!!!!

      implicit none

      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: a_nz(:)
      real (kind = SOLVER_PREC), intent(in) :: q, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, &
          tv_x, tv_y, tv_z

      a_nz = 0.0

      a_nz(1) = Real (q / (2.0 * a_1 * l) + 2.0 * u_x, kind = SOLVER_PREC)
      a_nz(2) = Real (- v_y, kind = SOLVER_PREC)
      a_nz(3) = Real (- w_z, kind = SOLVER_PREC)
      a_nz(4) = Real (2.0 * u_y - v_x, kind = SOLVER_PREC)
      a_nz(5) = Real (2.0 * u_z - w_x, kind = SOLVER_PREC)
      a_nz(6) = Real (- v_z - w_y, kind = SOLVER_PREC)
      a_nz(7) = Real (G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(8) = Real (- u_x, kind = SOLVER_PREC)
      a_nz(9) = Real (q / (2.0 * a_1 * l) + 2.0 * v_y, kind = SOLVER_PREC)
      a_nz(10) = Real (- w_z, kind = SOLVER_PREC)
      a_nz(11) = Real (2.0 * v_x - u_y, kind = SOLVER_PREC)
      a_nz(12) = Real (- u_z - w_x, kind = SOLVER_PREC)
      a_nz(13) = Real (2.0 * v_z - w_y, kind = SOLVER_PREC)
      a_nz(14) = Real (G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(15) = Real (- u_x, kind = SOLVER_PREC)
      a_nz(16) = Real (- v_y, kind = SOLVER_PREC)
      a_nz(17) = Real (q / (2.0 * a_1 * l) + 2.0 * w_z, kind = SOLVER_PREC)
      a_nz(18) = Real (- u_y - v_x, kind = SOLVER_PREC)
      a_nz(19) = Real (2.0 * w_x - u_z, kind = SOLVER_PREC)
      a_nz(20) = Real (2.0 * w_y - v_z, kind = SOLVER_PREC)
      a_nz(21) = Real (- 2.0 * G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(22) = Real (v_x, kind = SOLVER_PREC)
      a_nz(23) = Real (u_y, kind = SOLVER_PREC)
      a_nz(24) = Real (q / (3.0 * a_1 * l) + u_x + v_y, kind = SOLVER_PREC)
      a_nz(25) = Real (v_z, kind = SOLVER_PREC)
      a_nz(26) = Real (u_z, kind = SOLVER_PREC)
      a_nz(27) = Real (w_x, kind = SOLVER_PREC)
      a_nz(28) = Real (u_z, kind = SOLVER_PREC)
      a_nz(29) = Real (w_y, kind = SOLVER_PREC)
      a_nz(30) = Real (q / (3.0 * a_1 * l) + w_z + u_x, kind = SOLVER_PREC)
      a_nz(31) = Real (u_y, kind = SOLVER_PREC)
      a_nz(32) = Real (- G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(33) = Real (w_y, kind = SOLVER_PREC)
      a_nz(34) = Real (v_z, kind = SOLVER_PREC)
      a_nz(35) = Real (w_x, kind = SOLVER_PREC)
      a_nz(36) = Real (v_x, kind = SOLVER_PREC)
      a_nz(37) = Real (q / (3.0 * a_1 * l) + w_z + v_y, kind = SOLVER_PREC)
      a_nz(38) = Real (- G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(39) = Real (tv_x, kind = SOLVER_PREC)
      a_nz(40) = Real (tv_y, kind = SOLVER_PREC)
      a_nz(41) = Real (tv_z, kind = SOLVER_PREC)
      a_nz(42) = Real (q / (3.0 * a_2 * l) + u_x, kind = SOLVER_PREC)
      a_nz(43) = Real (u_y, kind = SOLVER_PREC)
      a_nz(44) = Real (u_z, kind = SOLVER_PREC)
      a_nz(45) = Real (tv_y, kind = SOLVER_PREC)
      a_nz(46) = Real (tv_x, kind = SOLVER_PREC)
      a_nz(47) = Real (tv_z, kind = SOLVER_PREC)
      a_nz(48) = Real (v_x, kind = SOLVER_PREC)
      a_nz(49) = Real (q / (3.0 * a_2 * l) + v_y, kind = SOLVER_PREC)
      a_nz(50) = Real (v_z, kind = SOLVER_PREC)
      a_nz(51) = Real (tv_z, kind = SOLVER_PREC)
      a_nz(52) = Real (tv_x, kind = SOLVER_PREC)
      a_nz(53) = Real (tv_y, kind = SOLVER_PREC)
      a_nz(54) = Real (w_x, kind = SOLVER_PREC)
      a_nz(55) = Real (w_y, kind = SOLVER_PREC)
      a_nz(56) = Real (q / (3.0 * a_2 * l) + w_z, kind = SOLVER_PREC)
      a_nz(57) = Real (- G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(58) = Real (tv_x, kind = SOLVER_PREC)
      a_nz(59) = Real (tv_y, kind = SOLVER_PREC)
      a_nz(60) = Real (tv_z, kind = SOLVER_PREC)
      a_nz(61) = Real (q / (b_2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix_non_zero


    subroutine Fill_in_a_matrix_moist (a, q, l, u_x, u_y, u_z, v_x, v_y, v_z, &
        w_x, w_y, w_z, tv_x, tv_y, tv_z)

      implicit none

      real (kind = SOLVER_PREC), dimension(:, :), intent(inout) :: a(:, :)
      real (kind = SOLVER_PREC), intent (in) :: q, l
      real (kind = SOLVER_PREC), intent (in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, &
          w_y, w_z, tv_x, tv_y, tv_z

      a = 0.0

        ! ur equation
      a(1, 1) = Real (q / (3.0 * a_2 * l) + u_x, kind = SOLVER_PREC)
      a(1, 2) = Real (u_y, kind = SOLVER_PREC)
      a(1, 3) = Real (u_z, kind = SOLVER_PREC)
      a(1, 4) = 0.0

        ! vr equation
      a(2, 1) = Real (v_x, kind = SOLVER_PREC)
      a(2, 2) = Real (q / (3.0 * a_2 * l) + v_y, kind = SOLVER_PREC)
      a(2, 3) = Real (v_z, kind = SOLVER_PREC)
      a(2, 4) = 0.0

        ! wr equation
      a(3, 1) = Real (w_x, kind = SOLVER_PREC)
      a(3, 2) = Real (w_y, kind = SOLVER_PREC)
      a(3, 3) = Real (q / (3.0 * a_2 * l) + w_z, kind = SOLVER_PREC)
      a(3, 4) = Real (-G_OVER_TREF, kind = SOLVER_PREC)

        ! rthetav equation
      a(4, 1) = Real (tv_x, kind = SOLVER_PREC)
      a(4, 2) = Real (tv_y, kind = SOLVER_PREC)
      a(4, 3) = Real (tv_z, kind = SOLVER_PREC)
      a(4, 4) = Real (2.0 * q / (b_2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix_moist


    subroutine Fill_in_a_matrix_non_zero_moist (a_nz, q, l, u_x, u_y, u_z, v_x, v_y, &
        v_z, w_x, w_y, w_z, tv_x, tv_y, tv_z)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Fills in the coefficients of the A Matrix                 !
      !                                                                    !
      ! Author: Pedro A. Jimenez                                           !
      !                                                                    !
      ! Comments: A * X = B wherein X are the tubulent fluxes of momentum, !
      !           theta_v, and theta2_v                                    !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~!!!!!!

      implicit none

      real (kind = SOLVER_PREC), dimension(:), intent(inout) :: a_nz(:)
      real (kind = SOLVER_PREC), intent(in) :: q, l
      real (kind = SOLVER_PREC), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, &
          tv_x, tv_y, tv_z

      a_nz = 0.0

      a_nz(1) = Real (q / (3.0 * a_2 * l) + u_x, kind = SOLVER_PREC)
      a_nz(2) = Real (u_y, kind = SOLVER_PREC)
      a_nz(3) = Real (u_z, kind = SOLVER_PREC)
      a_nz(4) = Real (v_x, kind = SOLVER_PREC)
      a_nz(5) = Real (q / (3.0 * a_2 * l) + v_y, kind = SOLVER_PREC)
      a_nz(6) = Real (v_z, kind = SOLVER_PREC)
      a_nz(7) = Real (w_x, kind = SOLVER_PREC)
      a_nz(8) = Real (w_y, kind = SOLVER_PREC)
      a_nz(9) = Real (q / (3.0 * a_2 * l) + w_z, kind = SOLVER_PREC)
      a_nz(10) = Real (-G_OVER_TREF, kind = SOLVER_PREC)
      a_nz(11) = Real (tv_x, kind = SOLVER_PREC)
      a_nz(12) = Real (tv_y, kind = SOLVER_PREC)
      a_nz(13) = Real (tv_z, kind = SOLVER_PREC)
      a_nz(14) = Real (2.0 * q / (b_2 * l), kind = SOLVER_PREC)

    end subroutine Fill_in_a_matrix_non_zero_moist


    subroutine Calc_heat_flux (th2, qv, fnm, fnp, i, j, k, tf_utv, &
        tf_vtv, tf_wtv, tf_uqv, tf_vqv, tf_wqv, tf_ut, tf_vt, tf_wt, &
        ims, ime, jms, jme, kms, kme)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                             !
      ! Purpose: Calc heat fluxes based on the moist fluses and     !
      !          virtual potential temperature fluxes               !
      !                                                             !
      ! Author: Pedro A. Jimenez                                    !
      !                                                             !
      ! Reference: YM75 (Eq. A13)                                   !
      !                                                             !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: i, k, j, ims, ime, jms, jme, kms, kme
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: th2, qv
      real, dimension (kms:kme), intent (in) :: fnm, fnp
      real (kind = SOLVER_PREC), intent(in) :: tf_utv, tf_vtv, tf_wtv, &
          tf_uqv, tf_vqv, tf_wqv
      real, intent(out) :: tf_ut, tf_vt, tf_wt

        ! Local vars
      real :: th_wall, qv_wall

      th_wall = T0 + fnm(k) * th2(i, k, j) + fnp(k) * th2(i, k - 1, j)
      qv_wall = fnm(k) * qv(i, k, j) + fnp(k) * qv(i, k - 1, j)

      tf_ut = Real ((tf_utv - P608 * th_wall * tf_uqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

      tf_vt = Real ((tf_vtv - P608 * th_wall * tf_vqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

      tf_wt = Real ((tf_wtv - P608 * th_wall * tf_wqv) / (1.0 + P608 * qv_wall), &
          kind = WORKING_PREC)

    end subroutine Calc_heat_flux


    subroutine Diagnose_fluxes_pbl_approx (do_scalar_mix, ss, nsteps, fnm, fnp, t, qv, ust, hfx, qfx, sm_tilde, sh_tilde, &
              sm_l2, sh_l2, rif, q_sq, q_sq_hl88, q_ratio, l_master, du_dz, dv_dz, dthetav_dz, dqv_dz, dqc_dz, &
              dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz, &
              u2, v2, w2, uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth, &
              uqc, vqc, wqc, uqr, vqr, wqr, uqi, vqi, wqi, uqs, vqs, wqs, uqg, vqg, wqg, uqnc, vqnc, wqnc, &
              uqnr, vqnr, wqnr, uqni, vqni, wqni, uqnwfa, vqnwfa, wqnwfa, uqnifa, vqnifa, wqnifa, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent (in)                                        :: do_scalar_mix, ss, nsteps
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz, dqv_dz, &
          dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: t, qv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: u2, v2, w2, &
          uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth, &
          uqc, vqc, wqc, uqr, vqr, wqr, uqi, vqi, wqi, uqs, vqs, wqs, uqg, vqg, wqg, &
          uqnc, vqnc, wqnc, uqnr, vqnr, wqnr, uqni, vqni, wqni, uqnwfa, vqnwfa, wqnwfa, uqnifa, vqnifa, wqnifa
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: sm_tilde, sh_tilde, rif
      real, dimension (its:ite, kts:kte, jts:jte), optional, intent (in)  :: q_ratio, q_sq_hl88, sm_l2, sh_l2
      real, dimension (ims:ime, jms:jme),          intent (in)  :: ust, hfx, qfx
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: qsq, l, sm, sh, ri, dudz, dvdz, dthvdz, dqvdz, dqcdz, dqrdz, dqidz, dqsdz, dqgdz, &
              dqncdz, dqnrdz, dqnidz, dqnwfadz, dqnifadz
      real :: t_wall, qv_wall
      real :: dxm, dym, tmpzx, tmpzy, sf_alpha
      real :: c_p, u1, v1, wind_speed, ust_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts + 1, ktf
          i_loop: do i = i_start, i_end

            if_low_q_sq: if (q_sq(i, k, j) > Q_SQ_MIN) then

                ! Save temporary scalars
              qsq = q_sq(i, k, j)
              sm = sm_tilde(i, k, j)
              sh = sh_tilde(i, k, j)
              l = l_master(i, k, j)
              ri = rif(i, k, j)
              dudz = du_dz(i, k, j)
              dvdz = dv_dz(i, k, j)
              dthvdz = dthetav_dz(i, k, j)
              dqvdz = dqv_dz(i, k, j)

                ! Eq. 27 MY82
              uw(i, k, j) = Diag_uw_pbl_approx (qsq, l, sm, dudz)
              vw(i, k, j) = Diag_vw_pbl_approx (qsq, l, sm, dvdz)
              uv(i, k, j) = Diag_uv_pbl_approx (qsq, l, dudz, dvdz, uw(i, k, j), vw(i, k, j))

                ! Eq. 28 MY82
              uthv(i, k, j) = Diag_uthv_pbl_approx (l, sh, sm, dthvdz, dudz)
              vthv(i, k, j) = Diag_vthv_pbl_approx (l, sh, sm, dthvdz, dvdz)
              wthv(i, k, j) = Diag_wthv_pbl_approx (qsq, l, sh, dthvdz)
              wqv(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqvdz)

                ! Eq. 29 MY82
              th2v(i, k, j) = Diag_th2v_pbl_approx (l, sh, dthvdz)

                ! Eq. 26 MY82
              u2(i, k, j) = Diag_u2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
              v2(i, k, j) = Diag_v2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))
              w2(i, k, j) = Diag_w2_pbl_approx (qsq, l, ri, dudz, dvdz, uw(i, k, j), vw(i, k, j), wthv(i, k, j))

              if (DEBUG) call Debug_negative_flux

              if (DEBUG) call Debug_big_qsq

                ! Diagnose the rest of qv fluxes
              uqv(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqvdz, uw(i, k, j), wqv(i, k, j))
              vqv(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqvdz, vw(i, k, j), wqv(i, k, j))

                ! Diagnose heat fluxes
              t_wall = T0 + fnm(k) * t(i, k, j) + fnp(k) * t(i, k - 1, j)
              qv_wall = fnm(k) * qv(i, k, j) + fnp(k) * qv(i, k - 1, j)

              uth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, uthv(i, k, j), uqv(i, k, j))
              vth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, vthv(i, k, j), vqv(i, k, j))
              wth(i, k, j) = Diag_xth_pbl_approx (t_wall, qv_wall, wthv(i, k, j), wqv(i, k, j))

              if ( do_scalar_mix .eq. 1 .and. ss .eq. nsteps) then
                  ! Cloud mass
                if ( p_qc .ge. param_first_scalar ) then
                  dqcdz = dqc_dz(i, k, j)
                  wqc(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqcdz)
                  uqc(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqcdz, uw(i, k, j), wqc(i, k, j))
                  vqc(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqcdz, vw(i, k, j), wqc(i, k, j))
                else
                  wqc(i, k, j) = TURB_FLUX_MIN
                  uqc(i, k, j) = TURB_FLUX_MIN
                  vqc(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Rain mass
                if ( p_qr .ge. param_first_scalar ) then
                  dqrdz = dqr_dz(i, k, j)
                  wqr(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqrdz)
                  uqr(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqrdz, uw(i, k, j), wqr(i, k, j))
                  vqr(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqrdz, vw(i, k, j), wqr(i, k, j))
                else
                  wqr(i, k, j) = TURB_FLUX_MIN
                  uqr(i, k, j) = TURB_FLUX_MIN
                  vqr(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Ice mass
                if ( p_qi .ge. param_first_scalar ) then
                  dqidz = dqi_dz(i, k, j)
                  wqi(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqidz)
                  uqi(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqidz, uw(i, k, j), wqi(i, k, j))
                  vqi(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqidz, vw(i, k, j), wqi(i, k, j))
                else
                  wqi(i, k, j) = TURB_FLUX_MIN
                  uqi(i, k, j) = TURB_FLUX_MIN
                  vqi(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Snow mass
                if ( p_qs .ge. param_first_scalar ) then
                  dqsdz = dqs_dz(i, k, j)
                  wqs(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqsdz)
                  uqs(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqsdz, uw(i, k, j), wqs(i, k, j))
                  vqs(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqsdz, vw(i, k, j), wqs(i, k, j))
                else
                  wqs(i, k, j) = TURB_FLUX_MIN
                  uqs(i, k, j) = TURB_FLUX_MIN
                  vqs(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Graupel mass
                if ( p_qg .ge. param_first_scalar ) then
                  dqgdz = dqg_dz(i, k, j)
                  wqg(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqgdz)
                  uqg(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqgdz, uw(i, k, j), wqg(i, k, j))
                  vqg(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqgdz, vw(i, k, j), wqg(i, k, j))
                else
                  wqg(i, k, j) = TURB_FLUX_MIN
                  uqg(i, k, j) = TURB_FLUX_MIN
                  vqg(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Cloud number
                if ( p_qnc .ge. param_first_scalar ) then
                  dqncdz = dqnc_dz(i, k, j)
                  wqnc(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqncdz)
                  uqnc(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqncdz, uw(i, k, j), wqnc(i, k, j))
                  vqnc(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqncdz, vw(i, k, j), wqnc(i, k, j))
                else
                  wqnc(i, k, j) = TURB_FLUX_MIN
                  uqnc(i, k, j) = TURB_FLUX_MIN
                  vqnc(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Rain number
                if ( p_qnr .ge. param_first_scalar ) then
                  dqnrdz = dqnr_dz(i, k, j)
                  wqnr(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqnrdz)
                  uqnr(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqnrdz, uw(i, k, j), wqnr(i, k, j))
                  vqnr(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqnrdz, vw(i, k, j), wqnr(i, k, j))
                else
                  wqnr(i, k, j) = TURB_FLUX_MIN
                  uqnr(i, k, j) = TURB_FLUX_MIN
                  vqnr(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Ice number
                if ( p_qni .ge. param_first_scalar ) then
                  dqnidz = dqni_dz(i, k, j)
                  wqni(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqnidz)
                  uqni(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqnidz, uw(i, k, j), wqni(i, k, j))
                  vqni(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqnidz, vw(i, k, j), wqni(i, k, j))
                else
                  wqni(i, k, j) = TURB_FLUX_MIN
                  uqni(i, k, j) = TURB_FLUX_MIN
                  vqni(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Water friendly aerosol number
                if ( p_qnwfa .ge. param_first_scalar ) then
                  dqnwfadz = dqnwfa_dz(i, k, j)
                  wqnwfa(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqnwfadz)
                  uqnwfa(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqnwfadz, uw(i, k, j), wqnwfa(i, k, j))
                  vqnwfa(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqnwfadz, vw(i, k, j), wqnwfa(i, k, j))
                else
                  wqnwfa(i, k, j) = TURB_FLUX_MIN
                  uqnwfa(i, k, j) = TURB_FLUX_MIN
                  vqnwfa(i, k, j) = TURB_FLUX_MIN
                end if
                  ! Ice friendly aerosol number
                if ( p_qnifa .ge. param_first_scalar ) then
                  dqnifadz = dqnifa_dz(i, k, j)
                  wqnifa(i, k, j) = Diag_wqv_pbl_approx (qsq, l, sh, dqnifadz)
                  uqnifa(i, k, j) = Diag_uqv_pbl_approx (qsq, l, dudz, dqnifadz, uw(i, k, j), wqnifa(i, k, j))
                  vqnifa(i, k, j) = Diag_vqv_pbl_approx (qsq, l, dvdz, dqnifadz, vw(i, k, j), wqnifa(i, k, j))
                else
                  wqnifa(i, k, j) = TURB_FLUX_MIN
                  uqnifa(i, k, j) = TURB_FLUX_MIN
                  vqnifa(i, k, j) = TURB_FLUX_MIN
                end if
              else
                  ! Cloud mas
                wqc(i, k, j) = TURB_FLUX_MIN
                uqc(i, k, j) = TURB_FLUX_MIN
                vqc(i, k, j) = TURB_FLUX_MIN
                  ! Rain mass
                wqr(i, k, j) = TURB_FLUX_MIN
                uqr(i, k, j) = TURB_FLUX_MIN
                vqr(i, k, j) = TURB_FLUX_MIN
                  ! Ice mass
                wqi(i, k, j) = TURB_FLUX_MIN
                uqi(i, k, j) = TURB_FLUX_MIN
                vqi(i, k, j) = TURB_FLUX_MIN
                  ! Snow mass
                wqs(i, k, j) = TURB_FLUX_MIN
                uqs(i, k, j) = TURB_FLUX_MIN
                vqs(i, k, j) = TURB_FLUX_MIN
                  ! Graupel mass
                wqg(i, k, j) = TURB_FLUX_MIN
                uqg(i, k, j) = TURB_FLUX_MIN
                vqg(i, k, j) = TURB_FLUX_MIN
                  ! Cloud number
                wqnc(i, k, j) = TURB_FLUX_MIN
                uqnc(i, k, j) = TURB_FLUX_MIN
                vqnc(i, k, j) = TURB_FLUX_MIN
                  ! Rain number
                wqnr(i, k, j) = TURB_FLUX_MIN
                uqnr(i, k, j) = TURB_FLUX_MIN
                vqnr(i, k, j) = TURB_FLUX_MIN
                  ! Ice number
                wqni(i, k, j) = TURB_FLUX_MIN
                uqni(i, k, j) = TURB_FLUX_MIN
                vqni(i, k, j) = TURB_FLUX_MIN
                  ! Water friendly aerosol number
                wqnwfa(i, k, j) = TURB_FLUX_MIN
                uqnwfa(i, k, j) = TURB_FLUX_MIN
                vqnwfa(i, k, j) = TURB_FLUX_MIN
                  ! Ice friendly aerosol number
                wqnifa(i, k, j) = TURB_FLUX_MIN
                uqnifa(i, k, j) = TURB_FLUX_MIN
                vqnifa(i, k, j) = TURB_FLUX_MIN
              end if

            else

              call Set_turb_flux_min

            end if if_low_q_sq

          end do i_loop
        end do k_loop
      end do j_loop

        ! Zero fluxes at the suface
      call Set_turb_flux_min_surface

        ! Zero fluxes at the top
      call Set_turb_flux_min_top

      contains

        subroutine Debug_negative_flux

          if (u2(i, k, j) < 0 .or. v2(i, k, j) < 0 .or. w2(i, k, j) < 0) then
            print *, '------------------------------------------------------'
            print *, 'NEGATIVE FLUX'
            print *, 'i, k, j = ', i, k, j
            if (u2(i, k, j) < 0) then
              print *, 'u2, qsq, u2/qsq = ', u2(i, k, j), qsq, ABS( u2(i, k, j) / qsq )
            else if (v2(i, k, j) < 0) then
              print *, 'v2, qsq, v2/qsq = ', v2(i, k, j), qsq, ABS( v2(i, k, j) / qsq )
            else if (w2(i, k, j) < 0) then
              print *, 'w2, qsq, w2/qsq = ', w2(i, k, j), qsq, ABS( w2(i, k, j) / qsq )
            end if
          end if

        end subroutine Debug_negative_flux


        subroutine Debug_big_qsq

          if (u2(i, k, j) + v2(i, k, j) + w2(i, k, j) > (1.0 + epsilon(qsq)) * qsq) then
            print *, '------------------------------------------------------'
            print *, 'BIG QSQ'
            print *, 'i, k, j = ', i, k, j
            print *, 'u2, v2, w2, sum, = ', u2(i, k, j), v2(i, k, j), w2(i, k, j), (u2(i, k, j) + v2(i, k, j) + w2(i, k, j))
            print *, 'qsq, sum/qsq = ', qsq, (u2(i, k, j) + v2(i, k, j) + w2(i, k, j)) / qsq
            u2(i, k, j) = qsq / 3.0
            v2(i, k, j) = qsq / 3.0
            w2(i, k, j) = qsq / 3.0
          end if

        end subroutine Debug_big_qsq


        subroutine Set_turb_flux_min

          u2(i, k, j) = TURB_FLUX_MIN
          v2(i, k, j) = TURB_FLUX_MIN
          w2(i, k, j) = TURB_FLUX_MIN
          uv(i, k, j) = TURB_FLUX_MIN
          uw(i, k, j) = TURB_FLUX_MIN
          vw(i, k, j) = TURB_FLUX_MIN
          uqv(i, k, j) = TURB_FLUX_MIN
          vqv(i, k, j) = TURB_FLUX_MIN
          wqv(i, k, j) = TURB_FLUX_MIN
          uqc(i, k, j) = TURB_FLUX_MIN
          vqc(i, k, j) = TURB_FLUX_MIN
          wqc(i, k, j) = TURB_FLUX_MIN
          uqr(i, k, j) = TURB_FLUX_MIN
          vqr(i, k, j) = TURB_FLUX_MIN
          wqr(i, k, j) = TURB_FLUX_MIN
          uqi(i, k, j) = TURB_FLUX_MIN
          vqi(i, k, j) = TURB_FLUX_MIN
          wqi(i, k, j) = TURB_FLUX_MIN
          uqs(i, k, j) = TURB_FLUX_MIN
          vqs(i, k, j) = TURB_FLUX_MIN
          wqs(i, k, j) = TURB_FLUX_MIN
          uqg(i, k, j) = TURB_FLUX_MIN
          vqg(i, k, j) = TURB_FLUX_MIN
          wqg(i, k, j) = TURB_FLUX_MIN
          uqnc(i, k, j) = TURB_FLUX_MIN
          vqnc(i, k, j) = TURB_FLUX_MIN
          wqnc(i, k, j) = TURB_FLUX_MIN
          uqnr(i, k, j) = TURB_FLUX_MIN
          vqnr(i, k, j) = TURB_FLUX_MIN
          wqnr(i, k, j) = TURB_FLUX_MIN
          uqni(i, k, j) = TURB_FLUX_MIN
          vqni(i, k, j) = TURB_FLUX_MIN
          wqni(i, k, j) = TURB_FLUX_MIN
          uqnwfa(i, k, j) = TURB_FLUX_MIN
          vqnwfa(i, k, j) = TURB_FLUX_MIN
          wqnwfa(i, k, j) = TURB_FLUX_MIN
          uqnifa(i, k, j) = TURB_FLUX_MIN
          vqnifa(i, k, j) = TURB_FLUX_MIN
          wqnifa(i, k, j) = TURB_FLUX_MIN
          uth(i, k, j) = TURB_FLUX_MIN
          uthv(i, k, j) = TURB_FLUX_MIN
          vth(i, k, j) = TURB_FLUX_MIN
          vthv(i, k, j) = TURB_FLUX_MIN
          wth(i, k, j) = TURB_FLUX_MIN
          wthv(i, k, j) = TURB_FLUX_MIN
          th2v(i, k, j) = TURB_FLUX_MIN

        end subroutine Set_turb_flux_min


        subroutine Set_turb_flux_min_surface

          do j = j_start, j_end
            do i = i_start, i_end
              u2(i, kts, j) = TURB_FLUX_MIN
              v2(i, kts, j) = TURB_FLUX_MIN
              w2(i, kts, j) = TURB_FLUX_MIN
              uv(i, kts, j) = TURB_FLUX_MIN
              uw(i, kts, j) = TURB_FLUX_MIN
              vw(i, kts, j) = TURB_FLUX_MIN
              uqv(i, kts, j) = TURB_FLUX_MIN
              vqv(i, kts, j) = TURB_FLUX_MIN
              wqv(i, kts, j) = TURB_FLUX_MIN
              uqc(i, kts, j) = TURB_FLUX_MIN
              vqc(i, kts, j) = TURB_FLUX_MIN
              wqc(i, kts, j) = TURB_FLUX_MIN
              uqr(i, kts, j) = TURB_FLUX_MIN
              vqr(i, kts, j) = TURB_FLUX_MIN
              wqr(i, kts, j) = TURB_FLUX_MIN
              uqi(i, kts, j) = TURB_FLUX_MIN
              vqi(i, kts, j) = TURB_FLUX_MIN
              wqi(i, kts, j) = TURB_FLUX_MIN
              uqs(i, kts, j) = TURB_FLUX_MIN
              vqs(i, kts, j) = TURB_FLUX_MIN
              wqs(i, kts, j) = TURB_FLUX_MIN
              uqg(i, kts, j) = TURB_FLUX_MIN
              vqg(i, kts, j) = TURB_FLUX_MIN
              wqg(i, kts, j) = TURB_FLUX_MIN
              uqnc(i, kts, j) = TURB_FLUX_MIN
              vqnc(i, kts, j) = TURB_FLUX_MIN
              wqnc(i, kts, j) = TURB_FLUX_MIN
              uqnr(i, kts, j) = TURB_FLUX_MIN
              vqnr(i, kts, j) = TURB_FLUX_MIN
              wqnr(i, kts, j) = TURB_FLUX_MIN
              uqni(i, kts, j) = TURB_FLUX_MIN
              vqni(i, kts, j) = TURB_FLUX_MIN
              wqni(i, kts, j) = TURB_FLUX_MIN
              uqnwfa(i, kts, j) = TURB_FLUX_MIN
              vqnwfa(i, kts, j) = TURB_FLUX_MIN
              wqnwfa(i, kts, j) = TURB_FLUX_MIN
              uqnifa(i, kts, j) = TURB_FLUX_MIN
              vqnifa(i, kts, j) = TURB_FLUX_MIN
              wqnifa(i, kts, j) = TURB_FLUX_MIN
              uth(i, kts, j) = TURB_FLUX_MIN
              uthv(i, kts, j) = TURB_FLUX_MIN
              vth(i, kts, j) = TURB_FLUX_MIN
              vthv(i, kts, j) = TURB_FLUX_MIN
              wth(i, kts, j) = TURB_FLUX_MIN
              wthv(i, kts, j) = TURB_FLUX_MIN
              th2v(i, kts, j) = TURB_FLUX_MIN
            end do
          end do

        end subroutine Set_turb_flux_min_surface

        subroutine Set_turb_flux_min_top

          do j = j_start, j_end
            do i = i_start, i_end
              u2(i, ktf + 1, j) = TURB_FLUX_MIN
              v2(i, ktf + 1, j) = TURB_FLUX_MIN
              w2(i, ktf + 1, j) = TURB_FLUX_MIN
              uv(i, ktf + 1, j) = TURB_FLUX_MIN
              uw(i, ktf + 1, j) = TURB_FLUX_MIN
              vw(i, ktf + 1, j) = TURB_FLUX_MIN
              uqv(i, ktf + 1, j) = TURB_FLUX_MIN
              vqv(i, ktf + 1, j) = TURB_FLUX_MIN
              wqv(i, ktf + 1, j) = TURB_FLUX_MIN
              uqc(i, ktf + 1, j) = TURB_FLUX_MIN
              vqc(i, ktf + 1, j) = TURB_FLUX_MIN
              wqc(i, ktf + 1, j) = TURB_FLUX_MIN
              uqr(i, ktf + 1, j) = TURB_FLUX_MIN
              vqr(i, ktf + 1, j) = TURB_FLUX_MIN
              wqr(i, ktf + 1, j) = TURB_FLUX_MIN
              uqi(i, ktf + 1, j) = TURB_FLUX_MIN
              vqi(i, ktf + 1, j) = TURB_FLUX_MIN
              wqi(i, ktf + 1, j) = TURB_FLUX_MIN
              uqs(i, ktf + 1, j) = TURB_FLUX_MIN
              vqs(i, ktf + 1, j) = TURB_FLUX_MIN
              wqs(i, ktf + 1, j) = TURB_FLUX_MIN
              uqg(i, ktf + 1, j) = TURB_FLUX_MIN
              vqg(i, ktf + 1, j) = TURB_FLUX_MIN
              wqg(i, ktf + 1, j) = TURB_FLUX_MIN
              uqnc(i, ktf + 1, j) = TURB_FLUX_MIN
              vqnc(i, ktf + 1, j) = TURB_FLUX_MIN
              wqnc(i, ktf + 1, j) = TURB_FLUX_MIN
              uqnr(i, ktf + 1, j) = TURB_FLUX_MIN
              vqnr(i, ktf + 1, j) = TURB_FLUX_MIN
              wqnr(i, ktf + 1, j) = TURB_FLUX_MIN
              uqni(i, ktf + 1, j) = TURB_FLUX_MIN
              vqni(i, ktf + 1, j) = TURB_FLUX_MIN
              wqni(i, ktf + 1, j) = TURB_FLUX_MIN
              uqnwfa(i, ktf + 1, j) = TURB_FLUX_MIN
              vqnwfa(i, ktf + 1, j) = TURB_FLUX_MIN
              wqnwfa(i, ktf + 1, j) = TURB_FLUX_MIN
              uqnifa(i, ktf + 1, j) = TURB_FLUX_MIN
              vqnifa(i, ktf + 1, j) = TURB_FLUX_MIN
              wqnifa(i, ktf + 1, j) = TURB_FLUX_MIN
              uth(i, ktf + 1, j) = TURB_FLUX_MIN
              uthv(i, ktf + 1, j) = TURB_FLUX_MIN
              vth(i, ktf + 1, j) = TURB_FLUX_MIN
              vthv(i, ktf + 1, j) = TURB_FLUX_MIN
              wth(i, ktf + 1, j) = TURB_FLUX_MIN
              wthv(i, ktf + 1, j) = TURB_FLUX_MIN
              th2v(i, ktf + 1, j) = TURB_FLUX_MIN
            end do
          end do

        end subroutine Set_turb_flux_min_top

    end subroutine Diagnose_fluxes_pbl_approx


    subroutine Diagnose_fluxes_surface (pbl3d_sfc_opt, u, v, t, rho, dz, qv, ust, hfx, qfx, rmol, &
              q_sq, l_master, u2, v2, w2, uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                             !
        ! Purpose: Calculates the surface fluxes                      !
        !                                                             !
        ! Author: M. Eghdami & Pedro A. Jimenez & Timothy W. Juliano  !
        !                                                             !
        ! Comments: Nakanishi 2001 for 3D                             !
        !                                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent (in)                                        :: pbl3d_sfc_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: t, qv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: u2, v2, w2, &
          uv, uw, vw, uthv, vthv, wthv, th2v, uqv, vqv, wqv, uth, vth, wth
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, rho, dz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: ust, hfx, qfx, rmol
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real :: zet, pmz1, pmz2, phz
      real :: qss, prt, gamma1, kz, elv
      real :: c_p, u1, v1, wind_speed
      real, dimension (its:ite, jts:jte) :: cosa, sina, ust_sq
      ! The parameters below depend on stability functions of module_sf_mynn
      REAL, PARAMETER :: cphm_st=5.0, cphm_unst=16.0, &
                         cphh_st=5.0, cphh_unst=16.0

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          c_p = CP * (1.0 + 0.84 * qv(i, kts, j))
          wth(i, kts, j) = hfx(i, j) / (rho(i, kts, j) * c_p)
          wthv(i, kts, j) = hfx(i, j) / (rho(i, kts, j) * CP)
          wqv(i, kts, j) = qfx(i, j) / rho(i, kts, j)
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          u1 = 0.5 * (u(i, 1, j) + u(i + 1, 1, j))
          v1 = 0.5 * (v(i, 1, j) + v(i, 1, j + 1))
          wind_speed = Max (Sqrt ((u1 ** 2.0) + (v1 ** 2.0)), 0.01)
          ust_sq(i, j) = ust(i, j) * ust(i, j)
          cosa(i, j) = u1 / wind_speed
          sina(i, j) = v1 / wind_speed
          uw(i, kts, j) = - ust_sq(i, j) * cosa(i, j)
          vw(i, kts, j) = - ust_sq(i, j) * sina(i, j)
        end do
      end do

        ! Calculate additional surface turbulent fluxes
      if (pbl3d_sfc_opt > 0) then
         ! Prandtl number
        gamma1 = 1.0 / 3.0 - 2.0 * (a_1 / b_1)
        prt = a_1 * (gamma1 - c_1) /  (a_2 * gamma1)
        do j = j_start, j_end
          do i = i_start, i_end
            if (pbl3d_sfc_opt == 1) then
              zet  = 0.0
              pmz1 = 1.0
              pmz2 = 1.0
              phz  = prt
              elv  = 1.0 ! L/kz
            else if (pbl3d_sfc_opt == 2) then
              kz  = 0.5 * dz(i, kts, j) * KARMAN
              zet = 0.5 * dz(i, kts, j) * rmol(i, j)
              elv = 0.5 * l_master(i, kts + 1, j) / kz

               ! MYNN similarity functions (from Businger)
              if ( zet >= 0.0 ) then
                pmz1 = 1.0 + (cphm_st - 1.0) * zet ! This will always appear in the denomenator
                pmz2 = 1.0 + (cphm_st) * zet
                phz = prt + cphh_st * zet
              else
                pmz1 = 1.0 / (1.0 - cphm_unst * zet) ** 0.25 - zet
                pmz2 = 1.0 / (1.0 - cphm_unst * zet) ** 0.25
                phz = prt / sqrt(1.0 - cphh_unst * zet) ! cphh_unst does not match businger
              end if

!              pmz1 = psim(i, j) - zet ! This will always appear in the denomenator
!              pmz2 = psim(i, j)
!              phz = psih(i, j)

!              if (pmz1 .le. 0.0) then
!                write(90, *) '-----------------------------'
!                write(90, *) 'i =', i, 'j = ', j, 'k = ', k
!                write(90, *) 'kz = ', kz
!                write(90, *) 'zet = ', zet
!                write(90, *) 'elv = ', elv
!                write(90, *) 'pmz1 = ', pmz1
!                write(90, *) 'pmz1 - zet = ', pmz1 - zet
!                write(90, *) 'pmz2 = ', pmz2
!                write(90, *) 'phz = ', phz
!              end if

            end if

             ! Following equations derived from MY74
             ! Eq. 46
            qss = Max(Q_SQ_MIN, (b_1 * elv * pmz1) ** (2.0 / 3.0) * ust_sq(i, j))
            q_sq(i, kts, j) = qss

             ! Other momentum Eq. 47a, b, c, d
            u2(i, kts, j) = qss * (1.0 / 3.0 + (a_1 / b_1) * &
                (4.0 * cosa(i, j) * cosa(i, j) * pmz2 - 2.0 * sina(i, j) * sina(i, j) * pmz2 + 2.0 * zet) / pmz1)
            v2(i, kts, j) = qss * (1.0 / 3.0 + (a_1 / b_1) * &
                (-2.0 * cosa(i, j) * cosa(i, j) * pmz2 + 4.0 * sina(i, j) * sina(i, j) * pmz2 + 2.0 * zet) / pmz1)
            w2(i, kts, j) = qss * (1.0 / 3.0 + (a_1 / b_1) * &
                (-2.0 * pmz2 - 4.0 * zet) / pmz1)
            uv(i, kts, j) = qss * (a_1 / b_1) * (6.0 * sina(i, j) * cosa(i, j) * pmz2) / pmz1

             ! Heat flux Eq. 48a, b
            uth(i, kts, j) = - wth(i, kts, j) * (3.0 * a_2 / b_1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * cosa(i, j)
            vth(i, kts, j) = - wth(i, kts, j) * (3.0 * a_2 / b_1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * sina(i, j)

             ! Moisture Eq. 48a, b
            uqv(i, kts, j) = - wqv(i, kts, j) * (3.0 * a_2 / b_1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * cosa(i, j)
            vqv(i, kts, j) = - wqv(i, kts, j) * (3.0 * a_2 / b_1 ** (1.0 / 3.0)) * &
                (elv ** (2.0 / 3.0)) * (phz + pmz2) / pmz1 ** (1.0 / 3.0) * sina(i, j)

             ! Heat variance Eq. 49
            th2v(i, kts, j) = wth(i, kts, j) ** 2.0 / MAX(ust_sq(i, j), 0.001) * (elv ** (2.0 / 3.0)) * &
                (b_2 / b_1 ** (1.0 / 3.0)) * phz / pmz1 ** (1.0 / 3.0)
          end do
        end do
      end if

    end subroutine Diagnose_fluxes_surface


    pure function Diag_uw_pbl_approx (q_sq, l_master, sm, du_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sm, du_dz
      real             :: diag_uw_pbl_approx

      diag_uw_pbl_approx = - Sqrt (q_sq) * l_master * sm * du_dz

    end function Diag_uw_pbl_approx


    pure function Diag_vw_pbl_approx (q_sq, l_master, sm, dv_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sm, dv_dz
      real             :: diag_vw_pbl_approx

      diag_vw_pbl_approx = - Sqrt (q_sq) * l_master * sm * dv_dz

    end function Diag_vw_pbl_approx


    pure function Diag_uv_pbl_approx (q_sq, l_master, du_dz, dv_dz, tf_uw, tf_vw)

      implicit none

      real, intent(in) :: q_sq, l_master, du_dz, dv_dz, tf_uw, tf_vw
      real             :: diag_uv_pbl_approx

      diag_uv_pbl_approx = 3.0 * a_1 * l_master / Sqrt (q_sq) * (- tf_uw * dv_dz - tf_vw * du_dz)

    end function Diag_uv_pbl_approx


    pure function Diag_uthv_pbl_approx (l_master, sh, sm, dthetav_dz, du_dz)

      implicit none

      real, intent(in) :: l_master, sh, sm, dthetav_dz, du_dz
      real             :: diag_uthv_pbl_approx

      diag_uthv_pbl_approx = 3.0 * a_2 * l_master ** 2.0 * (sh + sm) * dthetav_dz * du_dz

    end function Diag_uthv_pbl_approx


    pure function Diag_vthv_pbl_approx (l_master, sh, sm, dthetav_dz, dv_dz)

      implicit none

      real, intent(in) :: l_master, sh, sm, dthetav_dz, dv_dz
      real             :: diag_vthv_pbl_approx

      diag_vthv_pbl_approx = 3.0 * a_2 * l_master ** 2.0 * (sh + sm) * dthetav_dz * dv_dz

    end function Diag_vthv_pbl_approx


    pure function Diag_wthv_pbl_approx (q_sq, l_master, sh, dthetav_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sh, dthetav_dz
      real             :: diag_wthv_pbl_approx

      diag_wthv_pbl_approx = - Sqrt (q_sq) * l_master * sh * dthetav_dz

    end function Diag_wthv_pbl_approx


    pure function Diag_wqv_pbl_approx (q_sq, l_master, sh, dqv_dz)

      implicit none

      real, intent(in) :: q_sq, l_master, sh, dqv_dz
      real             :: diag_wqv_pbl_approx

      diag_wqv_pbl_approx = - Sqrt (q_sq) * l_master * sh * dqv_dz

    end function Diag_wqv_pbl_approx


    pure function Diag_th2v_pbl_approx (l_master, sh, dthetav_dz)

      implicit none

      real, intent(in) :: l_master, sh, dthetav_dz
      real             :: diag_th2v_pbl_approx

      diag_th2v_pbl_approx = b_2 * l_master ** 2.0 * sh * dthetav_dz ** 2.0

    end function Diag_th2v_pbl_approx


    pure function Diag_u2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_u2_pbl_approx

      diag_u2_pbl_approx = q_sq / 3.0 + 2.0 * a_1 * l_master / Sqrt (q_sq) *  &
          (- 2.0 * tf_uw * du_dz + tf_vw * dv_dz - G_OVER_TREF * tf_wtheta_v)

    end function Diag_u2_pbl_approx


    pure function Diag_v2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_v2_pbl_approx

      diag_v2_pbl_approx = q_sq / 3.0 + 2.0 * a_1 * l_master / Sqrt (q_sq) *  &
          (tf_uw * du_dz - 2.0 * tf_vw * dv_dz - G_OVER_TREF * tf_wtheta_v)

    end function Diag_v2_pbl_approx


    pure function Diag_w2_pbl_approx (q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v)

      implicit none

      real, intent(in) :: q_sq, l_master, rif, du_dz, dv_dz, tf_uw, tf_vw, tf_wtheta_v
      real             :: diag_w2_pbl_approx

      diag_w2_pbl_approx = q_sq / 3.0 + 2.0 * a_1 * l_master / Sqrt (q_sq) *  &
          (tf_uw * du_dz + tf_vw * dv_dz + 2.0 * G_OVER_TREF * tf_wtheta_v)

    end function Diag_w2_pbl_approx


    pure function Diag_uqv_pbl_approx (q_sq, l_master, du_dz, dqv_dz, tf_uw, tf_wqv)

      implicit none

      real, intent(in) :: q_sq, l_master, du_dz, dqv_dz, tf_uw, tf_wqv
      real             :: diag_uqv_pbl_approx

      diag_uqv_pbl_approx = - (3.0 * a_2 * l_master / Sqrt (q_sq)) * (tf_uw * dqv_dz + tf_wqv * du_dz)

    end function Diag_uqv_pbl_approx


    pure function Diag_vqv_pbl_approx (q_sq, l_master, dv_dz, dqv_dz, tf_vw, tf_wqv)

      implicit none

      real, intent(in) :: q_sq, l_master, dv_dz, dqv_dz, tf_vw, tf_wqv
      real             :: diag_vqv_pbl_approx

      diag_vqv_pbl_approx = - (3.0 * a_2 * l_master / Sqrt (q_sq)) * (tf_vw * dqv_dz + tf_wqv * dv_dz)

    end function Diag_vqv_pbl_approx


    pure function Diag_xth_pbl_approx (t_wall, qv_wall, tf_xtheta_v, tf_xqv)

      implicit none

      real, intent(in) :: t_wall, qv_wall, tf_xtheta_v, tf_xqv
      real             :: diag_xth_pbl_approx

      diag_xth_pbl_approx = (tf_xtheta_v - P608 * t_wall * tf_xqv) / (1.0 + P608 * qv_wall)

    end function Diag_xth_pbl_approx


    subroutine Calc_sm_sh_l2_or_l2p5 (my_level, du_dz, dv_dz, dthetav_dz, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, &
              kms, kme, its, ite, jts, jte, kts, kte, sm_l2, sh_l2, sm, sh, q_sq, q_sq_hl88, q_ratio, l_master_at_mass, l_master_at_face, rif)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates Sm and Sh for level 2 model following Yamada (1975)                !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in)                             :: my_level
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)   :: sm, sh
      real, dimension (its:ite, kts:kte, jts:jte), intent (out), optional :: rif, q_ratio
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: q_sq
      real, dimension (its:ite, kts:kte, jts:jte), intent (in),  optional :: sm_l2, sh_l2, q_sq_hl88
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in),  optional :: l_master_at_face
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf, k_start
      real :: ri, turb_prtl_num, gm, gh, phi1, phi2, phi3, phi4, phi5, N, D
      real, dimension (its:ite, kts:kte, jts:jte) :: l_master

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)
      k_start = kts + 1

      level: if (my_level == 'l2') then   ! level 2 model

        j_loop_l2: do j = j_start, j_end
          k_loop_l2: do k = k_start, ktf
            i_loop_l2: do i = i_start, i_end

              if (abs(du_dz(i, k, j)) < DWIND_DZ_MIN .and. abs(dv_dz(i, k, j)) < DWIND_DZ_MIN) then
                ri = G_OVER_TREF * dthetav_dz(i, k, j) / (2.0 * DWIND_DZ_MIN * DWIND_DZ_MIN)
              else
                ri = G_OVER_TREF * dthetav_dz(i, k, j) / (du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0)
              end if

                ! Eq. on page 929 Y75
              rif(i, k, j) = fac1 * (ri + fac2 - sqrt(ri * ri + fac3 * ri + fac4))
              if (rif(i, k, j) /= rif(i, k, j)) then
                write(90, *) 'i =', i, 'j = ', j, 'k = ', k
                write(90, *) '    dt_dz = ', dthetav_dz(i, k, j)
                write(90, *) '    du_dz = ', du_dz(i, k, j)
                write(90, *) '    dv_dz = ', dv_dz(i, k, j)
              end if

              if_turb: if (rif(i, k, j) < rif_c) then

                  ! Eq. 9a Y75
                sm(i, k, j) = cm * (rif_c - rif(i, k, j)) * (rf1 - rif(i, k, j)) / ((1.0 - rif(i, k, j)) * (rf2 - rif(i, k, j)))
                  ! Eq. 11 Y75
                turb_prtl_num = (ch2 / cm) * (rf2 - rif(i, k, j)) / (rf1 - rif(i, k, j))
                sh(i, k, j) = turb_prtl_num * sm(i, k, j)

              else

                sm(i, k, j) = TURB_FLUX_MIN
                sh(i, k, j) = TURB_FLUX_MIN

              end if if_turb

            end do i_loop_l2
          end do k_loop_l2
        end do j_loop_l2

        do j = j_start, j_end
          do i = i_start, i_end
            sm(i, kts, j) = TURB_FLUX_MIN
            sh(i, kts, j) = TURB_FLUX_MIN
            sm(i, ktf + 1, j) = TURB_FLUX_MIN
            sh(i, ktf + 1, j) = TURB_FLUX_MIN
          end do
        end do

      else if (my_level == 'l2p5') then   ! level 2.5 model

        if (PRESENT(l_master_at_mass)) then

          k_start = kts

          do j = j_start, j_end
            do k = k_start, ktf
              do i = i_start, i_end
                l_master(i, k, j) = l_master_at_mass(i, k, j)
              end do
            end do
          end do

        else if (PRESENT(l_master_at_face)) then

          k_start = kts + 1

          do j = j_start, j_end
            do k = k_start, ktf
              do i = i_start, i_end
                l_master(i, k, j) = l_master_at_face(i, k, j)
              end do
            end do
          end do

        end if

        j_loop_l2p5: do j = j_start, j_end
          k_loop_l2p5: do k = k_start, ktf
            i_loop_l2p5: do i = i_start, i_end

                  ! Realizability criteria from HL88
                if ( q_sq(i, k, j) < q_sq_hl88(i, k, j) ) then
                  q_ratio(i, k, j) = sqrt( q_sq(i, k, j) / q_sq_hl88(i, k, j) )
                  sm(i, k, j) = sm_l2(i, k, j) * q_ratio(i, k, j)
                  sh(i, k, j) = sh_l2(i, k, j) * q_ratio(i, k, j)
                else
                  q_ratio(i, k, j) = sqrt( q_sq(i, k, j) / q_sq_hl88(i, k, j) )
                    ! Eqs. 33a and 33b on page 855 MY82
                  gh = - ( (l_master(i, k, j) ** 2.0) / Max (q_sq(i, k, j), Q_SQ_MIN) ) * G_OVER_TREF * dthetav_dz(i, k, j)
                  gm = ( (l_master(i, k, j) ** 2.0) / Max (q_sq(i, k, j), Q_SQ_MIN) ) * ( Max( du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0, 1.0E-10 ) )

                    ! Eqs. 2.31, 2.33, 2.35 - 2.41 from HL88
                  phi1 = 1.0 - 3.0 * a_2 * b_2 * gh
                  phi2 = 1.0 - 9.0  * a_1 * a_2 * gh
                  phi3 = phi1 + 9.0 * a_2 * a_2 * gh
                  phi4 = phi1 - 12.0 * a_1 * a_2 * gh
                  phi5 = phi1 - 18.0 * a_1 * a_2 * gh
                  N = phi2 + 18.0 * c_1 * a_1 * a_1 * gm
                  D = phi2 * phi4 + 6.0 * a_1 * a_1 * phi3 * gm

                  sm(i, k, j) = a_1 * (phi3 - 3.0 * c_1 * phi4) / D
                  sh(i, k, j) = a_2 * N / D

                    ! Eqs. 3.17, 3.19, 3.20 from HL88
                  sm(i, k, j) = max ( min ( sm(i, k, j), sqrt ( 0.44 ** 2.0 / gm ) ), 0.0 )
                  sh(i, k, j) = max ( min ( sh(i, k, j), 0.76 * b_2 ), 0.0 )
                end if

            end do i_loop_l2p5
          end do k_loop_l2p5
        end do j_loop_l2p5

        do j = j_start, j_end
          do i = i_start, i_end
            if (PRESENT(l_master_at_face)) then
              sm(i, kts, j) = 0.0 !TURB_FLUX_MIN
              sh(i, kts, j) = 0.0 !TURB_FLUX_MIN
            end if
            sm(i, ktf + 1, j) = 0.0 !TURB_FLUX_MIN
            sh(i, ktf + 1, j) = 0.0 !TURB_FLUX_MIN
          end do
        end do

      end if level

    end subroutine Calc_sm_sh_l2_or_l2p5


    subroutine Calc_q_sq_l2_pbl_approx (l_master, l_dissip, du_dz, dv_dz, dthetav_dz, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, q_sq_face, q_sq_hl88, sm, rif)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates diagnostic TKE equation (level 2) for PBL approximation            !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out), optional :: q_sq_face
      real, dimension (its:ite, kts:kte, jts:jte), intent (out), optional :: q_sq_hl88
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: sm, rif
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: l_dissip
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, kts1, ktf
      real :: ri, turb_prtl_num
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      kts1 = kts + 1
      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts1, ktf
          i_loop: do i = i_start, i_end

            if_turb: if (rif(i, k, j) < rif_c) then
 
              q_sq(i, k, j) =  b_1 * l_master(i, k, j) ** 2.0 *        &
                  (du_dz(i, k, j) ** 2.0 + dv_dz(i, k, j) ** 2.0) *  &
                  (1.0 - rif(i, k, j)) * sm(i, k, j)

              q_sq(i, k, j) = Max (q_sq(i, k, j), Q_SQ_MIN)

              if (DEBUG) call Debug_big_qsq

            else

              q_sq(i, k, j) = Q_SQ_MIN

            end if if_turb

          end do i_loop
        end do k_loop
      end do j_loop

        ! Zero q_sq at the suface and top
      do j = j_start, j_end
        do i = i_start, i_end
          q_sq(i, kts, j) = Q_SQ_MIN
          q_sq(i, ktf + 1, j) = Q_SQ_MIN
        end do
      end do

      if (PRESENT(q_sq_face)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_face(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      if (PRESENT(q_sq_hl88)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_hl88(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      contains

        subroutine Debug_big_qsq

          if (q_sq(i, k, j) > 20.0) then
            print *, '-----------------------------'
            print *, 'i, k, j = ', i, k, j
            print *, 'q_sq = ', q_sq(i, k, j)
            print *, 'du_dz = ', du_dz(i, k, j)
            print *, 'dv_dz = ', dv_dz(i, k, j)
            print *, 'dthetav_dz = ', dthetav_dz(i, k, j)
            print *, 'l dissip = ', l_dissip(i, k, j)
            print *, 'Ri = ', ri
            print *, 'Ri flux = ', rif(i, k, j)
            print *, 'sm = ', sm(i, k, j)
          end if

        end subroutine Debug_big_qsq

    end subroutine Calc_q_sq_l2_pbl_approx


    subroutine Calc_q_sq_l2 (l_master, l_dissip, du_dx, dv_dy, dw_dz, du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz, &
              tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_wtv, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, q_sq_face, q_sq_hl88, rif)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates diagnostic TKE equation (level 2) for full 3D model                !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx, dv_dy, dw_dz, &
                                    du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: tf_u2, tf_v2, tf_w2, tf_uv, tf_uw, tf_vw, tf_wtv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out), optional :: q_sq_face
      real, dimension (its:ite, kts:kte, jts:jte), intent (out), optional :: q_sq_hl88
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: rif
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: l_dissip
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, kts1, ktf
      real :: ri, turb_prtl_num, q
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq
      logical, parameter :: DEBUG = .false.

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      kts1 = kts + 1
      ktf = Min (kte, kde - 1)

      j_loop: do j = j_start, j_end
        k_loop: do k = kts1, ktf
          i_loop: do i = i_start, i_end

            if_turb: if (rif(i, k, j) < rif_c) then

              q = b_1 * l_master(i, k, j) * (- tf_u2(i, k, j) * du_dx(i, k, j) - tf_v2(i, k, j) * dv_dy(i, k, j) &
                      - tf_w2(i, k, j) * dw_dz(i, k, j) - tf_uv(i, k, j) * (du_dy(i, k, j) + dv_dx(i, k, j)) &
                      - tf_uw(i, k, j) * (du_dz(i, k, j) + dw_dx(i, k, j)) &
                      - tf_vw(i, k, j) * (dv_dz(i, k, j) + dw_dy(i, k, j)) + tf_wtv(i, k, j) * G_OVER_TREF)
              q = q ** 1./3.
              q_sq(i, k, j) = Max (q * q, Q_SQ_MIN)

              if (DEBUG) call Debug_big_qsq

            else

              q_sq(i, k, j) = Q_SQ_MIN

            end if if_turb

          end do i_loop
        end do k_loop
      end do j_loop

        ! Zero q_sq at the suface and top
      do j = j_start, j_end
        do i = i_start, i_end
          q_sq(i, kts, j) = Q_SQ_MIN
          q_sq(i, ktf + 1, j) = Q_SQ_MIN
        end do
      end do

      if (PRESENT(q_sq_face)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_face(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      if (PRESENT(q_sq_hl88)) then

        do j = j_start, j_end
          do k = kts, ktf + 1
            do i = i_start, i_end
              q_sq_hl88(i, k, j) = q_sq(i, k, j)
            end do
          end do
        end do

      end if

      contains

        subroutine Debug_big_qsq

          if (q_sq(i, k, j) > 20.0) then
            print *, '-----------------------------'
            print *, 'i, k, j = ', i, k, j
            print *, 'q_sq = ', q_sq(i, k, j)
            print *, 'du_dz = ', du_dz(i, k, j)
            print *, 'dv_dz = ', dv_dz(i, k, j)
            print *, 'dthetav_dz = ', dthetav_dz(i, k, j)
            print *, 'l dissip = ', l_dissip(i, k, j)
            print *, 'Ri = ', ri
            print *, 'Ri flux = ', rif(i, k, j)
          end if

        end subroutine Debug_big_qsq

    end subroutine Calc_q_sq_l2


    subroutine Taper_q_sq (q_sq, ids, jds, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq
      integer,                                     intent (in)    :: ids, jds, ide, jde, kde,      &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      integer :: bdy_dist, taper_dist
      real    :: taper

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            bdy_dist = min( min((i-ids),(ide-i)) , min((j-jds),(jde-j)) )
            bdy_dist = max( 1, min( bdy_dist, 5 ) )
            taper_dist = 5
            if (bdy_dist .LE. taper_dist) then
              taper = max(0., min( 1., real(bdy_dist) / real(taper_dist) ) )
              q_sq(i, k, j) = q_sq(i, k, j) * taper
            end if
          end do
        end do
      end do

    end subroutine Taper_q_sq


    subroutine Diagnose_pbl_height (thetav, q_sq, z_at_w, dz, xland, zi, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                      !
      !                  NOTES ON THE PBLH FORMULATION                       !
      !                  ADOPTED FROM MYNN PBL SCHEME                        !
      !                                                                      !
      ! The 1.5-theta-increase method defines PBL heights as the level at    !
      ! which the potential temperature first exceeds the minimum potential  !
      ! temperature within the boundary layer by 1.5 K. When applied to      !
      ! observed temperatures, this method has been shown to produce PBL-    !
      ! height estimates that are unbiased relative to profiler-based        !
      ! estimates (Nielsen-Gammon et al. 2008). However, their study did not !
      ! include LLJs. Banta and Pichugina (2008) show that a TKE-based       !
      ! threshold is a good estimate of the PBL height in LLJs. Therefore,   !
      ! a hybrid definition is implemented that uses both methods, weighting !
      ! the TKE-method more during stable conditions (PBLH < 400 m).         !
      ! A variable tke threshold (TKEeps) is used since no hard-wired        !
      ! value could be found to work best in all conditions.                 !
      !                                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: thetav, q_sq, z_at_w, dz
      real, dimension (ims:ime, jms:jme),          intent (in)    :: xland
      real, dimension (ims:ime, jms:jme),          intent (out)   :: zi
      integer,                                     intent (in)    :: ids, jds, ide, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end
      integer :: kthv, ktke
      real    :: minthv, delt_thv, tke_eps, maxqke, qtke, qtkem1, pblh_tke, wt
      real, parameter :: sbl_lim  = 200.0 ! Upper limit of stable BL height (m)
      real, parameter :: sbl_damp = 400.0 ! Transition length for blending (m)

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Find minimum thetav in the lowest 200 m AGL
      do j = j_start, j_end
        do i = i_start, i_end
          k = kts + 1
          kthv = 1
          minthv = 9.E9

          do while (z_at_w(i, k, j) .le. 200.)
            if (minthv > thetav(i, k, j)) then
              minthv = thetav(i, k, j)
              kthv = k
            end if
            k = k + 1
          end do

            ! Find thetav-based PBLH (best for daytime)
          k = kthv + 1
          if ((xland(i, j) - 1.5) .ge. 0) then
              ! Water
            delt_thv = 1.0
          else
              ! Land
            delt_thv = 1.25
          ENDIF

          zi(i, j) = 0.
          k = kthv + 1
          do k = kts + 1, kte - 1
            if (thetav(i, k, j) .ge. (minthv + delt_thv)) then
              zi(i, j) = z_at_w(i, k, j) - dz(i, k - 1, j) * &
                     Min ((thetav(i, k, j) - (minthv + delt_thv)) / Max (thetav(i, k, j) - thetav(i, k - 1, j), 1E-6), 1.0)
            end if
            if (k .eq. kte - 1) zi(i, j) = z_at_w(i, kts + 1, j)
            if (zi(i, j) .ne. 0.0) exit
          end do

            ! For stable boundary layers, use TKE method to complement the
            ! thetav-based definition (when the thetav-based PBLH is below ~0.5 km).
            ! The tanh weighting function will make the TKE-based definition negligible
            ! when the thetav-based definition is above ~1 km.
          ktke = 1
          maxqke = Max (q_sq(i, kts, j), Q_SQ_MIN)
            ! Use 5% of tke max (Kosovic and Curry, 2000; JAS)
          tke_eps = maxqke / 40.
          tke_eps = Max (tke_eps, 0.02)
          pblh_tke = 0.

          k = ktke + 1
          do k = kts + 1, kte - 1
            qtke = Max (q_sq(i, k, j) / 2., Q_SQ_MIN)
            qtkem1 = Max (q_sq(i, k - 1, j) / 2., Q_SQ_MIN)
            if (qtke .le. tke_eps) then
              pblh_tke = z_at_w(i, k, j) - dz(i, k - 1, j) * &
                           Min ((tke_eps - qtke) / Max (qtkem1 - qtke, 1E-6), 1.0)
                ! In case of near zero TKE, set PBLH = lowest level
              pblh_tke = Max (pblh_tke, z_at_w(i, kts + 1, j))
            end if
            if (k .eq. kte - 1) pblh_tke = z_at_w(i, kts + 1, j)
            if (pblh_tke .ne. 0.) exit
          end do

            ! With TKE advection turned on, the TKE-based PBLH can be very large
            ! in grid points with convective precipitation (> 8 km!),
            ! so an artificial limit is imposed to not let pblh_tke exceed the
            ! theta_v-based PBL height +/- 350 m.
            ! This has no impact on 98-99% of the domain, but is the simplest patch
            ! that adequately addresses these extremely large PBLHs.
          pblh_tke = Min (pblh_tke, zi(i, j) + 350.)
          pblh_tke = Max (pblh_tke, Max (zi(i, j) - 350., 10.))

          wt = 0.5 * Tanh ((zi(i, j) - sbl_lim) / sbl_damp) + .5
          if (maxqke <= 0.05) then
            ! Cold pool situation - default to thetav-based def
          else
            ! Blend the two PBLH types here
            zi(i, j) = pblh_tke * (1. - wt) + zi(i, j) * wt
          end if
        end do
      end do

    end subroutine Diagnose_pbl_height


    subroutine Set_model_constants (pbl3d_constants)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Assigns values to the model constants                        !
      !                                                                       !
      ! Author: Pedro A. Jimenez                                              !
      !                                                                       !
      ! Method: 1) M73 constants                                              !
      !                                                                       !
      !         2) MY74 constants                                             !
      !                                                                       !
      !         3) MY82 constants                                             !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      character (len = *), intent(in) :: pbl3d_constants

      logical, parameter :: DEBUG = .false.

      real, parameter :: A1_M73 = 0.78, &
                         A2_M73 = 0.79, &
                         B1_M73 = 15.0, &
                         B2_M73 = 8.0,  &
                         C1_M73 = 0.056,&
                         C2_M73 = 0.0  ,&
                         C3_M73 = 0.0  ,&
                         ALPHA_M73 = 0.1, &
                         LAMBDA_M73 = 0.23

      real, parameter :: A1_MY74 = 0.78, &
                         A2_MY74 = 0.78, &
                         B1_MY74 = 15.0, &
                         B2_MY74 = 8.0,  &
                         C1_MY74 = 0.056,&
                         C2_MY74 = 0.0  ,&
                         C3_MY74 = 0.0  ,&
                         ALPHA_MY74 = 0.1, &
                         LAMBDA_MY74 = 0.23

      real, parameter :: A1_MY82 = 0.92, &
                         A2_MY82 = 0.74, &
                         B1_MY82 = 16.6, &
                         B2_MY82 = 10.1, &
                         C1_MY82 = 0.08 ,&
                         C2_MY82 = 0.0  ,&
                         C3_MY82 = 0.0  ,&
                         ALPHA_MY82 = 0.1, &
                         LAMBDA_MY82 = 0.23

      real, parameter :: A1_MYJ = 0.659888514560862645,    &
                         A2_MYJ = 0.6574209922667784586,   &
                         B1_MYJ = 11.87799326209552761,    &
                         B2_MYJ = 7.226971804046074028,    &
                         C1_MYJ = 0.000830955950095854396, &
                         C2_MYJ = 0.0,                     &
                         C3_MYJ = 0.0,                     &
                         ALPHA_MYJ = 0.3,                  &
                         LAMBDA_MYJ = 0.23

      real, parameter :: A1_Boulac = 0.3,     &
                         A2_Boulac = 0.33,    &
                         B1_Boulac = 8.4,     &
                         B2_Boulac = 6.4,     &
                         C1_Boulac = 0.08,    &
                         C2_Boulac = 0.0,     &
                         C3_Boulac = 0.0,     &
                         ALPHA_Boulac = 0.1,  &
                         LAMBDA_Boulac = 0.23

        ! Set model contants
      select case (trim(pbl3d_constants))

        case('M73')

          a_1 = A1_M73
          a_2 = A2_M73
          b_1 = B1_M73
          b_2 = B2_M73
          c_1 = C1_M73
          c_2 = C2_M73
          c_3 = C3_M73
          alpha = ALPHA_M73
          lambda = LAMBDA_M73

        case('MY74')

          a_1 = A1_MY74
          a_2 = A2_MY74
          b_1 = B1_MY74
          b_2 = B2_MY74
          c_1 = C1_MY74
          c_2 = C2_MY74
          c_3 = C3_MY74
          alpha = ALPHA_MY74
          lambda = LAMBDA_MY74

        case('MY82')

          a_1 = A1_MY82
          a_2 = A2_MY82
          b_1 = B1_MY82
          b_2 = B2_MY82
          c_1 = C1_MY82
          c_2 = C2_MY82
          c_3 = C3_MY82
          alpha = ALPHA_MY82
          lambda = LAMBDA_MY82

        case('MYJ')

          a_1 = A1_MYJ
          a_2 = A2_MYJ
          b_1 = B1_MYJ
          b_2 = B2_MYJ
          c_1 = C1_MYJ
          c_2 = C2_MYJ
          c_3 = C3_MYJ
          alpha = ALPHA_MYJ
          lambda = LAMBDA_MYJ

        case('Boulac')

          a_1 = A1_Boulac
          a_2 = A2_Boulac
          b_1 = B1_Boulac
          b_2 = B2_Boulac
          c_1 = C1_Boulac
          c_2 = C2_Boulac
          c_3 = C3_Boulac
          alpha = ALPHA_Boulac
          lambda = LAMBDA_Boulac

        case default

          write (OUTPUT_UNIT, *) 'ERROR: Unknown definition of model constants'
          write (OUTPUT_UNIT, *) 'Valid values of turbmix3d_constants are Y73, MY74 and MY82'
          stop

      end select

      if (DEBUG) then
        write (OUTPUT_UNIT, *) 'Inside subroutine Set_model_constants:'
        write (OUTPUT_UNIT, *) 'Model constants:'
        write (OUTPUT_UNIT, *) 'a1 = ', a_1
        write (OUTPUT_UNIT, *) 'a2 = ', a_2
        write (OUTPUT_UNIT, *) 'b1 = ', b_1
        write (OUTPUT_UNIT, *) 'b2 = ', b_2
        write (OUTPUT_UNIT, *) 'c1 = ', c_1
        write (OUTPUT_UNIT, *) 'c2 = ', c_2
        write (OUTPUT_UNIT, *) 'c3 = ', c_3
        write (OUTPUT_UNIT, *) 'Alpha = ', alpha
        write (OUTPUT_UNIT, *) 'Lambda = ', lambda
      end if

      if (DEBUG) print *, 'Done Set_model_constants...'

    end subroutine Set_model_constants


    subroutine Set_my_l2_global_vars

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Set global variables for Level 2                             !
      !                                                                       !
      ! Author: Pedro A. Jimenez (adapted from MYNN)                          !
      !                                                                       !
      ! Method: 1) M75 global vars                                            !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, parameter :: DEBUG = .false.

      real              :: e1, e2, e3, e4, e5

        ! For the level 2 model (Y75)
      e1 = b_1 - 6.0 * a_1
      e2 = b_1 + 12.0 * a_1 * (1.0 - c_2) + 3.0 * b_2 * (1.0 - c_3)
      e3 = b_1 * (1.0 - 3.0 * c_1) - 6.0 * a_1
      e4 = b_1 * (1.0 - 3.0 * c_1) + 12.0 * a_1 * (1.0 - c_2) + 9.0 * a_2 * (1.0 - c_2)
      e5 = b_1 + 3.0 * a_1 * (1.0 - c_2) + 3.0 * b_2 * (1.0 - c_3)

      rif_c = e1 / e2
      rf1   = e3 / e4
      rf2   = e1 / e5
      cm    = a_1 * e2 * e4 / (b_1 * e5)
      ch2   = a_2 * e2 / b_1

      fac1 = a_2 * e5 / (2.0 * a_1 * e4)
      fac2 = a_1 * e3 / (a_2 * e5)
      fac3 = 2.0 * (a_1 / a_2) * (e3 * e5 - 2.0 * e1 * e4 ) / (e5 * e5)
      fac4 = (a_1 * e3 / (a_2 * e5) ) ** 2.0

      if (DEBUG) print *, 'Done Set_my_l2_global_vars...'

    end subroutine Set_my_l2_global_vars


    subroutine Set_my_qsq_min (pbl3d_prog)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Set minimum values for Q_SQ and TURB_FLUX                    !
      !                                                                       !
      ! Author: Timothy W. Juliano                                            !
      !                                                                       !
      ! Method: Depending on level 2 or 2.5 model                             !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: pbl3d_prog

      logical, parameter :: DEBUG = .false.

      if (pbl3d_prog == 0) then
        Q_SQ_MIN = 1.0E-5
        TURB_FLUX_MIN = 1.0E-12
      else
        Q_SQ_MIN = 1.0E-5
        TURB_FLUX_MIN = 0.0
      end if

      if (DEBUG) print *, 'Done Set_my_qsq_min...'

    end subroutine Set_my_qsq_min


    subroutine Set_init_turb_state_my (restart, pbl3d_constants,                &
              pbl3d_nsteps, pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt,  &
              ids, ide, jds, jde, kds, kde,                                     &
              ims, ime, jms, jme, kms, kme,                                     &
              its, ite, jts, jte, kts, kte,                                     &
              itf, jtf, ktf, i_start, i_end, j_start, j_end,                    &
              fnm, fnp, dnw, ph, phb,                                           &
              rdzw, z_at_mass, dz,                                              &
              t_2, qv, alt, u, v, ust, hfx, qfx, rmol, rho,                     &
              du_dz, dv_dz, dthetav_dz, dqv_dz,                                 &
              q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear,              &
              q_sq_buoyancy, q_sq_dissip,                                       &
              l_master, l_master_at_mass,                                       &
              turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,           &
              turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                   &
              turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,       &
              turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,             &
              turb_flux_uqv, turb_flux_vqv, turb_flux_wqv                       )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Initialize turbulence state                                  !
      !                                                                       !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano (adapted from MYNN)     !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, parameter :: DEBUG = .true.

      logical, intent(in) :: restart
      integer, intent(in) :: pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt, pbl3d_nsteps
      character (len = *), intent(in) :: pbl3d_constants
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                  intent(in)  :: fnm, fnp, dnw
      real, dimension(ims:ime, jms:jme),          intent(in)  :: ust, hfx, qfx, rmol
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: t_2, qv, &
                                                                 u, v, alt, ph, phb
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: rho
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: q_sq_prog, q_sq,     &
          q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: l_master, l_master_at_mass, &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv

      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: rdzw, z_at_mass, dz
      real, dimension(its:ite, kts:kte, jts:jte), intent(in)  :: du_dz, dv_dz, &
                                                                 dthetav_dz, dqv_dz
      integer, intent(in)  :: itf, jtf, ktf, i_start, i_end, j_start, j_end

        ! Local vars
      integer :: i, j, k
      integer :: init_iter
      real, dimension (its:ite, kts:kte, jts:jte) :: sm, sh, rif

      if (DEBUG) then
        write (OUTPUT_UNIT, *) 'Inside subroutine Init_turb_mixing_3d:'
        write (OUTPUT_UNIT, *) 'PBL3D option = ', pbl3d_opt
        write (OUTPUT_UNIT, *) 'Model constants = ', trim(pbl3d_constants)
        write (OUTPUT_UNIT, *) 'PBL3D number of substeps = ', pbl3d_nsteps
      end if

      call Set_model_constants (pbl3d_constants)

      call Set_my_l2_global_vars

      call Set_my_qsq_min (pbl3d_prog)

      if (restart) return

        ! Initializes thetav, TKE, L, and the turbulent fluxes
      j_loop_init: do j = jts, jtf
        i_loop_init: do i = its, itf
          q_sq(i, kts:kte, j) = Q_SQ_MIN
          l_master(i, kts:kte, j) = Q_SQ_MIN
          if (pbl3d_prog == 1) then
            q_sq_prog(i, kts:ktf, j) = Q_SQ_MIN
            q_sq_vdiff(i, kts:ktf, j) = 0.0
            q_sq_hdiff(i, kts:ktf, j) = 0.0
            q_sq_shear(i, kts:ktf, j) = 0.0
            q_sq_buoyancy(i, kts:ktf, j) = 0.0
            q_sq_dissip(i, kts:ktf, j) = 0.0
            l_master_at_mass(i, kts:ktf, j) = Q_SQ_MIN
          end if
          turb_flux_u2(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_v2(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_w2(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_uv(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_uw(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_vw(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_utheta_v(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_vtheta_v(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_wtheta_v(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_theta2_v(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_utheta(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_vtheta(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_wtheta(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_uqv(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_vqv(i, kts:kte, j) = TURB_FLUX_MIN
          turb_flux_wqv(i, kts:kte, j) = TURB_FLUX_MIN
          rho(i, kts:ktf, j) = 1. / alt(i, kts:ktf, j)
        end do i_loop_init
      end do j_loop_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FOR NOW, INITIALIZE TO MIN VALS !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        ! Calc stability functions
!      call Calc_sm_sh_l2_or_l2p5 (my_level='l2', du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, &
!          ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
!          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, sm=sm, sh=sh, rif=rif)

!      k = 1
!      do j = j_start, j_end
!        do i = i_start, i_end
!          q_sq(i, k, j) = Max(Q_SQ_MIN, b_1 ** (2.0 / 3.0) * ust(i, j) * ust(i, j))
!        end do
!      end do

!      do init_iter = 1, 2
!          ! Calc q_sq
!        call Calc_q_sq_l2_pbl_approx (l_master=l_master, du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, ids=ids, ide=ide, &
!            jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, &
!            jts=jts, jte=jte, kts=kts, kte=kte, q_sq_face=q_sq, sm=sm, rif=rif)

!          ! Diagnose master length scale
!         call Calc_l_master_algebra (config_flags%pbl3d_l_opt, q_sq, dz, rdzw, rmol, turb_flux_wtheta_v, dthetav_dz, l_master, &
!             ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!      end do

!          ! Diagnose momentum, heat, moisture fluxes for pbl approx
!      call Diagnose_fluxes_pbl_approx (fnm=fnm, fnp=fnp, t=t_2, qv=qv, ust=ust, hfx=hfx, qfx=qfx, &
!          sm_tilde=sm, sh_tilde=sh, rif=rif, q_sq=q_sq, l_master=l_master, &
!          du_dz=du_dz, dv_dz=dv_dz, dthetav_dz=dthetav_dz, dqv_dz=dqv_dz, u2=turb_flux_u2, v2=turb_flux_v2, w2=turb_flux_w2, uv=turb_flux_uv, &
!          uw=turb_flux_uw, vw=turb_flux_vw, uthv=turb_flux_utheta_v, vthv=turb_flux_vtheta_v, wthv=turb_flux_wtheta_v, th2v=turb_flux_theta2_v, &
!          uqv=turb_flux_uqv, vqv=turb_flux_vqv, wqv=turb_flux_wqv, uth=turb_flux_utheta, vth=turb_flux_vtheta, wth=turb_flux_wtheta, ids=ids, ide=ide, &
!          jds=jds, jde=jde, kds=kds, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

!      call Diagnose_fluxes_surface (pbl3d_sfc_opt, u, v, t_2, rho, dz, qv, ust, hfx, qfx, &
!          rmol, q_sq, l_master, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
!          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
!          turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, &
!          turb_flux_wtheta, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

!      if (pbl3d_prog == 1) then
!        call Fill_q_sq_prog_with_q_sq (q_sq_prog, q_sq, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!          ! Diagnose master length scale at mass points for prognostic TKE
!        call Fill_l_mass_with_l_face (l_master_at_face, l_master_at_mass, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
!      end if

      if (DEBUG) print *, 'Done Set_init_turb_state_my...'

    end subroutine Set_init_turb_state_my

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! atx_cr() computes A'*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) W(N), the value of A'*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(n+1)
  integer ja(nz_num)
  integer k1
  integer k2
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(ja(k1:k2)) = w(ja(k1:k2)) + a(k1:k2) * x(i)
  end do

  return
end
subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! ATX_ST computes A'*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A'.
!
!    Output, real ( kind = rk ) W(N), the value of A'*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(nz_num)
  integer j
  integer ja(nz_num)
  integer k
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(j) = w(j) + a(k) * x(i)
  end do

  return
end
subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_CR computes A*x for a matrix stored in sparse compressed row form.
!
!  Discussion:
!
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) W(N), the value of A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(n+1)
  integer ja(nz_num)
  integer k1
  integer k2
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end
subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

!*****************************************************************************80
!
!! AX_ST computes A*x for a matrix stored in sparset triplet form.
!
!  Discussion:
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input, real ( kind = rk ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = rk ) W(N), the value of A*X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(nz_num)
  integer j
  integer ja(nz_num)
  integer k
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  w(1:n) = 0.0D+00

  do k = 1, nz_num
    i = ia(k)
    j = ja(k)
    w(i) = w(i) + a(k) * x(j)
  end do

  return
end
subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

!*****************************************************************************80
!
!! DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!    The array UA can be used to locate the diagonal elements of the matrix.
!
!    It is assumed that every row of the matrix includes a diagonal element,
!    and that the elements of each row have been ascending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!    On output, the order of the entries of JA may have changed because of
!    the sorting.
!
!    Output, integer UA(N), the index of the diagonal element
!    of each row.
!
  implicit none

  integer n
  integer nz_num

  integer i
  integer ia(n+1)
  integer k
  integer ja(nz_num)
  integer ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

!*****************************************************************************80
!
!! ILU_CR computes the incomplete LU factorization of a matrix.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input, integer UA(N), the index of the diagonal element
!    of each row.
!
!    Output, real ( kind = rk ) L(NZ_NUM), the ILU factorization of A.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(n+1)
  integer iw(n)
  integer j
  integer ja(nz_num)
  integer jj
  integer jrow
  integer jw
  integer k
  real ( kind = rk ) l(nz_num)
  real ( kind = rk ) tl
  integer ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end
subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

!*****************************************************************************80
!
!! LUS_CR applies the incomplete LU preconditioner.
!
!  Discussion:
!
!    The linear system M * Z = R is solved for Z.  M is the incomplete
!    LU preconditioner matrix, and R is a vector supplied by the user.
!    So essentially, we're solving L * U * Z = R.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = rk ) L(NZ_NUM), the matrix values.
!
!    Input, integer UA(N), the index of the diagonal element
!    of each row.
!
!    Input, real ( kind = rk ) R(N), the right hand side.
!
!    Output, real ( kind = rk ) Z(N), the solution of the system M * Z = R.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  integer i
  integer ia(n+1)
  integer j
  integer ja(nz_num)
  real ( kind = rk ) l(nz_num)
  real ( kind = rk ) r(n)
  integer ua(n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end
subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, &
  tol_rel )

!*****************************************************************************80
!
!! MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
!
!  Discussion:
!
!    The linear system A*X=B is solved iteratively.
!
!    The matrix A is assumed to be stored in sparse triplet form.  Only
!    the nonzero entries of A are stored.  For instance, the K-th nonzero
!    entry in the matrix is stored by:
!
!      A(K) = value of entry,
!      IA(K) = row of entry,
!      JA(K) = column of entry.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, integer NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the matrix values.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = rk ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer MR, the maximum number of (inner) iterations
!    to take.  0 < MR <= N.
!
!    Input, real ( kind = rk ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = rk ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mr
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) av
  real ( kind = rk ) c(1:mr)
  real ( kind = rk ), parameter :: delta = 1.0D-03
  real ( kind = rk ) g(1:mr+1)
  real ( kind = rk ) h(1:mr+1,1:mr)
  real ( kind = rk ) htmp
  integer i
  integer ia(nz_num)
  integer itr
  integer itr_max
  integer itr_used
  integer j
  integer ja(nz_num)
  integer k
  integer k_copy
  real ( kind = rk ) mu
  real ( kind = rk ) r(1:n)
  real ( kind = rk ) rho
  real ( kind = rk ) rho_tol
  real ( kind = rk ) rhs(1:n)
  real ( kind = rk ) s(1:mr)
  real ( kind = rk ) tol_abs
  real ( kind = rk ) tol_rel
  real ( kind = rk ) v(1:n,1:mr+1)
  logical, parameter :: verbose = .true.
  real ( kind = rk ) x(1:n)
  real ( kind = rk ) y(1:mr+1)

  itr_used = 0

  if ( n < mr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
    write ( *, '(a)' ) '  N < MR.'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  MR = ', mr
    stop
  end if

  do itr = 1, itr_max

    call ax_st ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    rho = sqrt ( dot_product ( r(1:n), r(1:n) ) )

    if ( verbose ) then
      write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_st ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - h(j,k) * v(1:n,j)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( av + delta * h(k+1,k) == av ) then

        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do

        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then

        y(1:k+1) = h(1:k+1,k)

        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y(1:k+1) )
        end do

        h(1:k+1,k) = y(1:k+1)

      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g(1:k+1) )
      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i8,a,g14.6)' ) '  K =   ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)'       ) ' '
    write ( *, '(a)'       ) 'MGMRES_ST:'
    write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine mult_givens ( c, s, k, g )

!*****************************************************************************80
!
!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
!
!  Discussion:
!
!    In order to make it easier to compare this code with the Original C,
!    the vector indexing is 0-based.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, real ( kind = rk ) C, S, the cosine and sine of a Givens
!    rotation.
!
!    Input, integer K, indicates the location of the first
!    vector entry.
!
!    Input/output, real ( kind = rk ) G(1:K+1), the vector to be modified.
!    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k

  real ( kind = rk ) c
  real ( kind = rk ) g(1:k+1)
  real ( kind = rk ) g1
  real ( kind = rk ) g2
  real ( kind = rk ) s

  g1 = c * g(k) - s * g(k+1)
  g2 = s * g(k) + c * g(k+1)

  g(k)   = g1
  g(k+1) = g2

  return
end
subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel )

!*****************************************************************************80
!
!! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!
!  Discussion:
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!
!    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
!    corrections to the code on 31 May 2007.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 2012
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the linear system.
!
!    Input, integer NZ_NUM, the number of nonzero matrix values.
!
!    Input, integer IA(N+1), JA(NZ_NUM), the row and column indices
!    of the matrix values.  The row vector has been compressed.
!
!    Input, real ( kind = rk ) A(NZ_NUM), the matrix values.
!
!    Input/output, real ( kind = rk ) X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real ( kind = rk ) RHS(N), the right hand side of the linear system.
!
!    Input, integer ITR_MAX, the maximum number of (outer) 
!    iterations to take.
!
!    Input, integer MR, the maximum number of (inner) iterations 
!    to take.  MR must be less than N.
!
!    Input, real ( kind = rk ) TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real ( kind = rk ) TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer mr
  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  real ( kind = rk ) av
  real ( kind = rk ) c(mr+1)
  real ( kind = rk ), parameter :: delta = 1.0D-03
  real ( kind = rk ) g(mr+1)
  real ( kind = rk ) h(mr+1,mr)
  real ( kind = rk ) htmp
  integer i
  integer ia(n+1)
  integer itr
  integer itr_max
  integer itr_used
  integer j
  integer ja(nz_num)
  integer k
  integer k_copy
  real ( kind = rk ) l(ia(n+1)+1)
  real ( kind = rk ) mu
  real ( kind = rk ) r(n)
  real ( kind = rk ) rho
  real ( kind = rk ) rho_tol
  real ( kind = rk ) rhs(n)
  real ( kind = rk ) s(mr+1)
  real ( kind = rk ) tol_abs
  real ( kind = rk ) tol_rel
  integer ua(n)
  real ( kind = rk ) v(n,mr+1);
  logical, parameter :: verbose = .true.
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(mr+1)

  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end
subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Input, integer NZ_NUM, the number of nonzeros.
!
!    Input, integer IA(N+1), the compressed row indices.
!
!    Input/output, integer JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = rk ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nz_num

  real ( kind = rk ) a(nz_num)
  integer i
  integer ia(n+1)
  integer i4temp
  integer ja(nz_num)
  integer k
  integer l
  real ( kind = rk ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end module module_pbl3d_my
