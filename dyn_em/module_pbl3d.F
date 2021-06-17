  module module_pbl3d

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
    use module_pbl3d_my, ONLY : Set_init_turb_state_my, Calc_turb_fluxes_my
    use module_pbl3d_my, ONLY : Q_SQ_MIN, TURB_FLUX_MIN, a_2, b_1

    implicit none

    private

    public :: Init_pbl3d, &
              Init_temp_arrays_for_substep, Update_wrf_tends_temp_state_and_zero_tends, Allocate_temp_tends, Deallocate_temp_tends, &
              Fill_l_mass_with_l_face, Calc_turb_fluxes_driver, Vertical_turb_mix, Horizontal_turb_mix

    real, parameter :: G_OVER_TREF = G / T0
    logical, parameter :: USE_HIGH_ORDER = .true.

    integer, parameter :: I_TO_TEST = 4, J_TO_TEST = 4, K_TO_TEST = 3

  contains

    subroutine Init_temp_arrays_for_substep (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, &
              u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Prepare for sub-stepping                                     !
      !              1. Save copy of state variables                          !
      !              2. Allocate temporary tendency arrays                    !
      !              3. Zero out temporary tendency arrays                    !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u_2, v_2, w_2, t_2, moist, rho, q_sq_2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real, dimension (:,:,:), allocatable, intent (inout)      :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                   t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                                     intent (in)  :: ide, jde, kde,                &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

      call Save_copy_state (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, &
              rho_tmp, q_sq_tmp, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Allocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      call Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, its, ite, jts, jte, kts, kte)

    end subroutine Init_temp_arrays_for_substep


    subroutine Save_copy_state (u_2, v_2, w_2, t_2, moist, rho, q_sq_2, u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, &
              rho_tmp, q_sq_tmp, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Save copy of state arrays for sub-stepping                   !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u_2, v_2, w_2, t_2, moist, rho, q_sq_2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      integer,                                     intent (in)  :: ide, jde, kde,                &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes
      i_start = its - 1
      i_end = ite + 1
      j_start = jts - 1
      j_end = jte + 1
      k_start = kts

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            u_tmp(i,k,j) = u_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            v_tmp(i,k,j) = v_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = kte

      DO j = j_start - 1, j_end + 1
        DO k = k_start, k_end
          DO i = i_start - 1, i_end + 1
            w_tmp(i,k,j) = w_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            th_tmp(i,k,j) = t_2(i,k,j)
            qv_tmp(i,k,j) = moist(i,k,j)
            rho_tmp(i,k,j) = rho(i,k,j)
            q_sq_tmp(i,k,j) = q_sq_2(i,k,j)
          ENDDO
        ENDDO
      ENDDO

    end subroutine Save_copy_state


    subroutine Allocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Allocate temporary tendencies for sub-stepping               !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (:,:,:), allocatable, intent (inout)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                  t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)        :: its, ite, jts, jte, kts, kte

      IF (.NOT. ALLOCATED(ru_tendf_tmp))   ALLOCATE ( ru_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(rv_tendf_tmp))   ALLOCATE ( rv_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(rw_tendf_tmp))   ALLOCATE ( rw_tendf_tmp(its:ite,kts:kte,jts:jte)   )
      IF (.NOT. ALLOCATED(t_tendf_tmp))    ALLOCATE ( t_tendf_tmp(its:ite,kts:kte,jts:jte)    )
      IF (.NOT. ALLOCATED(moist_tend_tmp)) ALLOCATE ( moist_tend_tmp(its:ite,kts:kte,jts:jte) )
      IF (.NOT. ALLOCATED(q_sq_tend_tmp))  ALLOCATE ( q_sq_tend_tmp(its:ite,kts:kte,jts:jte) )

    end subroutine Allocate_temp_tends


    subroutine Deallocate_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Deallocate temporary tendencies for sub-stepping             !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (:,:,:), allocatable, intent (inout)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                  t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)        :: its, ite, jts, jte, kts, kte

      IF (ALLOCATED(ru_tendf_tmp))    DEALLOCATE ( ru_tendf_tmp   )
      IF (ALLOCATED(rv_tendf_tmp))    DEALLOCATE ( rv_tendf_tmp   )
      IF (ALLOCATED(rw_tendf_tmp))    DEALLOCATE ( rw_tendf_tmp   )
      IF (ALLOCATED(t_tendf_tmp))     DEALLOCATE ( t_tendf_tmp    )
      IF (ALLOCATED(moist_tend_tmp))  DEALLOCATE ( moist_tend_tmp )
      IF (ALLOCATED(q_sq_tend_tmp))   DEALLOCATE ( q_sq_tend_tmp  )

    end subroutine Deallocate_temp_tends


    subroutine Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
              t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, ide, jde, kde, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Zero out temporary tendencies for sub-stepping               !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (out)  :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                    t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                              intent (in)          :: ide, jde, kde, its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            ru_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rv_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rw_tendf_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            t_tendf_tmp(i,k,j) = 0.0
            moist_tend_tmp(i,k,j) = 0.0
            q_sq_tend_tmp(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

    end subroutine Zero_temp_tends


    subroutine Update_wrf_tends_temp_state_and_zero_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, &
              moist_tend, q_sq_tend, ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, dt, c1h, c2h, c1f, c2f, muu, muv, mut, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                          !
      ! Purpose: Update arrays for sub-stepping                                  !
      !              1. Update WRF tendencies with temporary tendencies          !
      !              2. Update temporary state arrays with temporary tendencies  !
      !              3. Zero out temporary tendency arrays                       !
      !                                                                          !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                            !
      !                                                                          !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout)  :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real,                                        intent (in)     :: dt
      real, dimension (kms:kme),                   intent (in)     :: c1h, c2h, c1f, c2f
      real, dimension (ims:ime, jms:jme),          intent (in)     :: muu, muv, mut
      integer,                                     intent (in)     :: ide, jde, kde,                &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

      call Update_wrf_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend, &
          ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Update_temp_state (u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
          ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
          dt, c1h, c2h, c1f, c2f, muu, muv, mut, ide, jde, kde, ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte)

      call Zero_temp_tends (ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, &
          q_sq_tend_tmp, ide, jde, kde, its, ite, jts, jte, kts, kte)

    end subroutine Update_wrf_tends_temp_state_and_zero_tends


    subroutine Update_wrf_tends (ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Update WRF tendency arrays with temporary tendencies         !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: ru_tendf, rv_tendf, rw_tendf, t_tendf, moist_tend, q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      integer,                                     intent (in)     :: ide, jde, kde,                &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            ru_tendf(i,k,j) = ru_tendf(i,k,j) + ru_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rv_tendf(i,k,j) = rv_tendf(i,k,j) + rv_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            rw_tendf(i,k,j) = rw_tendf(i,k,j) + rw_tendf_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            t_tendf(i,k,j) = t_tendf(i,k,j) + t_tendf_tmp(i,k,j)
            moist_tend(i,k,j) = moist_tend(i,k,j) + moist_tend_tmp(i,k,j)
            q_sq_tend(i,k,j) = q_sq_tend(i,k,j) + q_sq_tend_tmp(i,k,j)
          ENDDO
        ENDDO
      ENDDO

    end subroutine Update_wrf_tends


    subroutine Update_temp_state (u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp, &
              ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp, &
              dt, c1h, c2h, c1f, c2f, muu, muv, mut, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                       !
      ! Purpose: Updates temporary state arrays for sub-stepping              !
      !                                                                       !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                         !
      !                                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout)  :: u_tmp, v_tmp, w_tmp, th_tmp, qv_tmp, rho_tmp, q_sq_tmp
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)     :: ru_tendf_tmp, rv_tendf_tmp, rw_tendf_tmp, &
                                                                      t_tendf_tmp, moist_tend_tmp, q_sq_tend_tmp
      real,                                        intent (in)     :: dt
      real, dimension (kms:kme),                   intent (in)     :: c1h, c2h, c1f, c2f
      real, dimension (ims:ime, jms:jme),          intent (in)     :: muu, muv, mut 
      integer,                                     intent (in)     :: ide, jde, kde, &
                                                                      ims, ime, jms, jme, kms, kme, &
                                                                      its, ite, jts, jte, kts, kte

        ! Local
      integer :: i, j, k, i_start, i_end, j_start, j_end, k_start, k_end

        ! Calculate indexes for u
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            u_tmp(i,k,j) = u_tmp(i,k,j) + (ru_tendf_tmp(i,k,j) * (dt / (c1h(k) * muu(i,j) + c2h(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for v
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = jte
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            v_tmp(i,k,j) = v_tmp(i,k,j) + (rv_tendf_tmp(i,k,j) * (dt / (c1h(k) * muv(i,j) + c2h(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for w
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = kte

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            w_tmp(i,k,j) = w_tmp(i,k,j) + (rw_tendf_tmp(i,k,j) * (dt / (c1f(k) * mut(i,j) + c2f(k))))
          ENDDO
        ENDDO
      ENDDO

        ! Calculate indexes for theta/qv
      i_start = its
      i_end = Min ( ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      k_start = kts
      k_end = Min (kte, kde - 1)

      DO j = j_start, j_end
        DO k = k_start, k_end
          DO i = i_start, i_end
            th_tmp(i,k,j) = th_tmp(i,k,j) + (t_tendf_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k))))
            qv_tmp(i,k,j) = qv_tmp(i,k,j) + (moist_tend_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k))))
            rho_tmp(i,k,j) = rho_tmp(i,k,j)
            q_sq_tmp(i,k,j) = max ( Q_SQ_MIN , q_sq_tmp(i,k,j) + (q_sq_tend_tmp(i,k,j) * (dt / (c1h(k) * mut(i,j) + c2h(k)))) )
          ENDDO
        ENDDO
      ENDDO

    end subroutine Update_temp_state


    subroutine Calc_turb_fluxes_driver (config_flags, dz, z_at_mass, z_at_w, u, v, w, th, th2, qv, qc, qr, qi, qs, qg, &
              qnc, qnr, qni, qnwfa, qnifa, &
              thetav, tsk, t2, q2, psfc, rho, ust, hfx, qfx, rmol, &
              itimestep, ss, n_tracer, q_sq, q_sq_prog, q_sq_tend, &
              q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, &
              khx, khy, khz, l_master, l_master_at_mass, l_boulac, turb_flux_u2, turb_flux_v2, turb_flux_w2,   &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
              turb_flux_theta2_v, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv, turb_flux_uqc, turb_flux_vqc, turb_flux_wqc, &
              turb_flux_uqr, turb_flux_vqr, turb_flux_wqr, turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, &
              turb_flux_uqs, turb_flux_vqs, turb_flux_wqs, turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, &
              turb_flux_uqnc, turb_flux_vqnc, turb_flux_wqnc, turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, &
              turb_flux_uqni, turb_flux_vqni, turb_flux_wqni, turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, &
              turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa, &
              turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, turb_flux_u2_mass, turb_flux_v2_mass, &
              turb_flux_w2_mass, turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, turb_flux_wtheta_v_mass, &
              msfux, msfuy, msfvx, msfvy, msftx, msfty, dx, dy, rdx, rdy, dn, dnw, rdz, dt, id, rdzw, fnm, fnp, cf1, cf2, cf3, &
              zx, zy, mut, c1h, c2h, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

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
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: qc, qr, qi, qs, qg, qnc, qnr, qni, qnwfa, qnifa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: q_sq, q_sq_prog
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_hdiff_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_shear_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_buoyancy_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: khx, khy, khz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: thetav
      real, dimension (ims:ime, kms:kme, jms:jme), intent (inout) :: l_master, l_master_at_mass, l_boulac
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, jms:jme),          intent (in)    :: tsk, t2, q2, psfc, ust, hfx, qfx, rmol
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
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2_mass, turb_flux_v2_mass, turb_flux_w2_mass, &
                                                                     turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, &
                                                                     turb_flux_wtheta_v_mass
      real, dimension (kms:kme),                   intent (in)    :: fnm, fnp, dn, dnw
      real,                                        intent (in)    :: dx, dy, rdx, rdy, cf1, cf2, cf3, dt
      integer,                                     intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                                                     ims, ime, jms, jme, kms, kme, &
                                                                     its, ite, jts, jte, kts, kte

        ! Local vars
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, &
          dw_dx, dw_dy, dw_dz, dthetav_dx, dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz,     &
          dqc_dx, dqc_dy, dqc_dz, dqr_dx, dqr_dy, dqr_dz, dqi_dx, dqi_dy, dqi_dz, dqs_dx, dqs_dy, dqs_dz, &
          dqg_dx, dqg_dy, dqg_dz, dqnc_dx, dqnc_dy, dqnc_dz, dqnr_dx, dqnr_dy, dqnr_dz, &
          dqni_dx, dqni_dy, dqni_dz, dqnwfa_dx, dqnwfa_dy, dqnwfa_dz, dqnifa_dx, dqnifa_dy, dqnifa_dz
      logical, parameter :: DEBUG = .false.
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf

        ! Calculate derivatives at the faces to prepare for fluxes
      call Calc_wind_derivatives (config_flags, u, v, w, du_dx, dv_dy, dw_dz,  &
          du_dy, dv_dx, dw_dx, du_dz, dw_dy, dv_dz,                            &
          msftx, msfty, msfux, msfuy, msfvx, msfvy,                            &
          rdx, rdy, dn, dnw, dz, rdz,                                          &
          rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy,                               &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its,     &
          ite, jts, jte, kts, kte)

      call Calc_thetav (th2, qv, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, thetav)

      call Calc_mass_var_derivatives (thetav, qv,                      &
          msftx, msfty, msfux, msfuy, msfvx, msfvy, rdx, rdy,          &
          fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn,         &
          dthetav_dx, dthetav_dy, dthetav_dz,                          &
          dqv_dx, dqv_dy, dqv_dz,                                      &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,  &
          its, ite, jts, jte, kts, kte)

      if ( config_flags%pbl3d_scalar_mix .eq. 1 .and. ss .eq. config_flags%pbl3d_nsteps ) then
        call Calc_scalar_var_derivatives (qc, qr, qi, qs, qg, qnc, qnr, qni, qnwfa, qnifa, msftx, msfty, msfux, msfuy, msfvx, msfvy, &
            rdx, rdy, fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn, dqc_dx, dqc_dy, dqc_dz, dqr_dx, dqr_dy, dqr_dz, &
            dqi_dx, dqi_dy, dqi_dz, dqs_dx, dqs_dy, dqs_dz, dqg_dx, dqg_dy, dqg_dz, dqnc_dx, dqnc_dy, dqnc_dz, &
            dqnr_dx, dqnr_dy, dqnr_dz, dqni_dx, dqni_dy, dqni_dz, dqnwfa_dx, dqnwfa_dy, dqnwfa_dz, dqnifa_dx, dqnifa_dy, dqnifa_dz, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Calculate indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (DEBUG .and. itimestep == 1) then
        print *, 'its = ', its, 'ite = ', ite
        print *, 'ims = ', ims, 'ime = ', ime
        print *, 'jts = ', jts, 'jte = ', jte
        print *, 'jms = ', jms, 'jme = ', jme
      end if

      call Calc_turb_fluxes_my (config_flags, dz, z_at_mass, z_at_w, u, v, w, th, th2, qv, thetav, tsk, t2, q2, psfc, rho, ust, hfx, qfx, rmol, &
          du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dqv_dx, dqv_dy, dqv_dz, dthetav_dx, dthetav_dy, dthetav_dz, &
          dqc_dz, dqr_dz, dqi_dz, dqs_dz, dqg_dz, dqnc_dz, dqnr_dz, dqni_dz, dqnwfa_dz, dqnifa_dz, itimestep, ss, n_tracer, q_sq, q_sq_prog, &
          l_master, l_master_at_mass, l_boulac, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_uqv, &
          turb_flux_vqv, turb_flux_wqv, turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqc, turb_flux_vqc, &
          turb_flux_wqc, turb_flux_uqr, turb_flux_vqr, turb_flux_wqr, turb_flux_uqi, turb_flux_vqi, turb_flux_wqi, &
          turb_flux_uqs, turb_flux_vqs, turb_flux_wqs, turb_flux_uqg, turb_flux_vqg, turb_flux_wqg, turb_flux_uqnc, turb_flux_vqnc, &
          turb_flux_wqnc, turb_flux_uqnr, turb_flux_vqnr, turb_flux_wqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_wqni, &
          turb_flux_uqnwfa, turb_flux_vqnwfa, turb_flux_wqnwfa, turb_flux_uqnifa, turb_flux_vqnifa, turb_flux_wqnifa, &
          msfux, msfuy, msfvx, msfvy, msftx, msfty, &
          dx, dy, rdx, rdy, dn, dnw, rdz, dt, id, rdzw, fnm, fnp, cf1, cf2, cf3, zx, zy, mut, c1h, c2h, ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calculate eddy diffusivities for tracer mixing
      if ( ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar ) then
        call Calc_kh (du_dz, dv_dz, dthetav_dz, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, &
            l_master, q_sq, khx, khy, khz, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

        ! Predict q_sq if level 2.5 model
      if ( config_flags%pbl3d_prog .gt. 0 ) then   ! Prognostic TKE
        call Calc_q_sq_l2p5 (config_flags, u, v, w, zx, zy, th2, qv, rho, thetav, tsk, t2, q2, psfc, ust, hfx, &
            du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, dthetav_dz, q_sq, q_sq_prog, q_sq_tend, &
            q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, l_master, l_master_at_mass, l_boulac, &
            turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, &
            turb_flux_u2_mass, turb_flux_v2_mass, turb_flux_w2_mass, turb_flux_uv_mass, turb_flux_uw_mass, turb_flux_vw_mass, &
            turb_flux_wtheta_v_mass, dn, dnw, fnm, fnp, msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1h, c2h, z_at_mass, rdx, rdy, &
            rdz, rdzw, cf1, cf2, cf3, dt, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      if (DEBUG) print *, 'Done Calc_turb_fluxes_driver...'

    end subroutine Calc_turb_fluxes_driver


    subroutine Calc_wind_derivatives (config_flags, u, v, w, du_dx, dv_dy, dw_dz, du_dy, dv_dx, dw_dx, &
        du_dz, dw_dy, dv_dz, msftx, msfty, msfux, msfuy, msfvx, msfvy, rdx, rdy, dn, dnw, dz, rdz, rdzw, fnm, fnp, &
        cf1, cf2, cf3, zx, zy, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                                         !
        ! Purpose: Calculate spatial derivatives of the wind variables at the center of the faces !
        !                                                                                         !
        ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                          !
        !                                                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, zx, zy, dz, rdz, rdzw
      real, dimension (kms:kme), intent (in)                    :: fnm, fnp, dn, dnw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx, dv_dy, dw_dz, du_dy, dv_dx, &
                                                                   dw_dx, du_dz, dw_dy, dv_dz

      logical, parameter :: DEBUG = .false.

      call Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dw_dz (w, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      call Calc_du_dx (u, zx, rdz, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3,     &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx  )

      call Calc_dv_dx (v, zx, rdz, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, ide, jde,  &
          kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dx     )

      call Calc_dw_dx (w, zx, rdz, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, config_flags, dw_dx    )

      call Calc_du_dy (u, zy, rdz, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3,     &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy  )

      call Calc_dv_dy (v, zy, rdz, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds,    &
          jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy  )

      call Calc_dw_dy (w, zy, rdz, msfty, rdy, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy                  )

      if (DEBUG) print *, 'Done Calc_wind_derivatives...'

    end subroutine Calc_wind_derivatives


    subroutine Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, rdz
      real, dimension (kms:kme),                   intent (in)  :: dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: u_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Unstagger U
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            u_at_mass(i, k, j) = 0.5 * (u(i, k, j) + u(i + 1, k, j)) 
          end do
        end do
      end do

        ! Calculate derivative at the center of the faces
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dz(i, k, j) = (u_at_mass(i, k, j) - u_at_mass(i, k - 1, j)) * rdz(i, k, j) 
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dz(i, kts, j) = 0.0
            ! Top
          du_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dz


    subroutine Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, rdz
      real, dimension (kms:kme),                   intent (in)  :: dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: v_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Unstagger V
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            v_at_mass(i, k, j) = 0.5 * (v(i, k, j) + v(i, k, j + 1))
          end do
        end do
      end do

        ! Calculate derivative at the center of the faces
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dz(i, k, j) = (v_at_mass(i, k, j) - v_at_mass(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dv_dz(i, kts, j) = 0.0
            ! Top face
          dv_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dz


    subroutine Calc_dw_dz (w, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dz                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, rdz
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dz

        ! Local vars
      integer :: i, j, k
      real, dimension (its:ite, kts:kte, jts:jte) :: w_at_mass
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate w at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            w_at_mass(i, k, j) = 0.5 * (w(i, k, j) + w(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate partial w / partial z
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dw_dz(i, k, j) = (w_at_mass(i, k, j) - w_at_mass(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dw_dz(i, kts, j) =  0.0
            ! Top face
          dw_dz(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dw_dz


    subroutine Calc_du_dx (u, zx, rdz, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3,         &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx      )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zx, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfux
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_ux
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dx, u_mux, u_mux_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at u points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_ux(i, j) = msftx(i, j) * msfux(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            u_mux(i, k, j) = u(i, k, j) / msfux(i, j)
          end do
        end do
      end do

        ! Interpolate u_mux from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_mux_at_full(i, k, j) = fnm(k) * u_mux(i, k, j) + fnp(k) * u_mux(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_mux_at_full(i, 1, j) =  cf1 * u_mux(i, 1, j) + cf2 * u_mux(i, 2, j) + cf3 * u_mux(i, 3, j) 
            ! Top face
          u_mux_at_full(i, kte, j) =  cft1 *  u_mux(i, ktes1, j) +  cft2 * u_mux(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! du_dz_dz_dx = (partial u_mux / partial z)|sigma (du_dz) * 
              ! (partial z / partial x)|sigma (dz_dx)
            du_dz_dz_dx(i, k, j) = (0.5 * (u_mux(i, k, j) + u_mux(i + 1, k, j)) - &
                0.5 * (u_mux(i, k - 1, j) + u_mux(i + 1, k - 1, j))) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial u / partial x)|z = msf_tx_ux * (partial u_mux_at_full / partial x)|sigma - du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dx(i, k, j) = msf_tx_ux(i, j) * ((u_mux_at_full(i + 1, k, j) - u_mux_at_full(i, k, j)) * rdx -  &
                du_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dx(i, kts, j) = 0.0
            ! Top face
          du_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dx


    subroutine Calc_dv_dx (v, zx, rdz, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, ide, jde,         &
              kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dx            )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zx, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfvx
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_vx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dx, v_mvx, v_mvx_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_vx(i, j) = msftx(i, j) * msfvx(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection (at u points)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            v_mvx(i, k, j) = 0.25 * (v(i, k, j) + v(i, k, j + 1) + v(i - 1, k, j) + v(i - 1, k, j + 1)) / msfvx(i, j)
          end do
        end do
      end do

        ! Interpolate v_mvx from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            v_mvx_at_full(i, k, j) = fnm(k) * v_mvx(i, k, j) + fnp(k) * v_mvx(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          v_mvx_at_full(i, 1, j) = cf1 * v_mvx(i, 1, j) + cf2 * v_mvx(i, 2, j) + cf3 * v_mvx(i, 3, j)
            ! Top face
          v_mvx_at_full(i, kte, j) = cft1 * v_mvx(i, ktes1, j) + cft2 * v_mvx(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dv_dz_dz_dx = (partial v_mvx / partial z)|sigma (dv_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dv_dz_dz_dx(i, k, j) = (0.5 * (v_mvx(i, k, j) + v_mvx(i + 1, k, j)) - &
                0.5 * (v_mvx(i, k - 1, j) + v_mvx(i + 1, k - 1, j))) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial v / partial x)|z = msf_tx_vx * (partial v_mvx_at_full / partial x)|sigma - dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dx(i, k, j) = msf_tx_vx(i, j) * ((v_mvx_at_full(i + 1, k, j) - v_mvx_at_full(i, k, j)) * rdx -  &
                dv_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dv_dx(i, 1, j) =  0.0
          dv_dx(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dv_dx


    subroutine Calc_dw_dx (w, zx, rdz, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, config_flags, dw_dx           )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dx                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zx, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real,                                        intent (in)  :: rdx
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dx

        ! Local vars
      integer :: i, j, k, i_offset
      real    :: dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_tx
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dx, w_mtx, w_mtx_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_tx_tx(i, j) = msftx(i, j) * msftx(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        i_offset = 2
      else
        i_offset = 1
      end if
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start - i_offset, i_end + i_offset
            w_mtx(i, k, j) = w(i, k, j) / msftx(i, j)
          end do
        end do
      end do

        ! Average w_mtx to the mass points (half levels)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start - 1, i_end + 1
            w_mtx_at_half(i, k, j) = 0.5 * (w_mtx(i, k, j) + w_mtx(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dw_dz_dz_dx = (partial w_mtx / partial z)|sigma (dw_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dw_dz_dz_dx(i, k, j) = (w_mtx_at_half(i, k, j) - w_mtx_at_half(i, k - 1, j)) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial w / partial x)|z = msf_tx_tx * (partial w_mtx / partial x)|sigma - dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * (8.0 * (w_mtx(i + 1, k, j) - w_mtx(i - 1, k, j)) - (w_mtx(i + 2, k, j) - w_mtx(i - 2, k, j)) ) / 12.0 - dw_dz_dz_dx(i, k, j))
            else
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * 0.5 * (w_mtx(i + 1, k, j) - w_mtx(i - 1, k, j)) - dw_dz_dz_dx(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dw_dx(i, 1, j) = 0.0
          dw_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dx


    subroutine Calc_du_dy (u, zy, rdz, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3,       &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy    )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zy, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfuy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_uy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dy, u_muy, u_muy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_uy(i, j) = msfty(i, j) * msfuy(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection (at v points)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            u_muy(i, k, j) = 0.25 * (u(i, k, j) + u(i + 1, k, j) + u(i, k, j - 1) + u(i + 1, k, j - 1)) / msfuy(i, j) 
          end do
        end do
      end do

        ! Interpolate u_muy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            u_muy_at_full(i, k, j) = fnm(k) * u_muy(i, k, j) + fnp(k) * u_muy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          u_muy_at_full(i, 1, j) =  cf1 * u_muy(i, 1, j) + cf2 * u_muy(i, 2, j) + cf3 * u_muy(i, 3, j)
            ! Top face
          u_muy_at_full(i, kte, j) = cft1 * u_muy(i, ktes1, j) +  cft2 * u_muy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dy)
      do j = j_start, j_end 
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces 
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! du_dz_dz_dy = (partial u_muy / partial z)|sigma (du_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            du_dz_dz_dy(i, k, j) = (0.5 * (u_muy(i, k, j) + u_muy(i, k, j + 1)) - &
                0.5 * (u_muy(i, k - 1, j) + u_muy(i, k - 1, j + 1))) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial u / partial y)|z = msf_ty_uy * (partial u_muy_at_full / partial y)|sigma - du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            du_dy(i, k, j) = msf_ty_uy(i, j) * ((u_muy_at_full(i, k, j + 1) - u_muy_at_full(i, k, j)) * rdy -  &
                du_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          du_dy(i, kts, j) = 0.0
            ! Top face
          du_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dy


    subroutine Calc_dv_dy (v, zy, rdz, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds,     &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy   )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zy, rdz
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfvy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jds, jde, kde,            &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_vy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dy, v_mvy, v_mvy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_vy(i, j) = msfty(i, j) * msfvy(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! Because msfvy at the poles will be undefined (1./0.),
              ! we will have trouble. But we are OK since v at the
              ! poles is 0., and that takes precedence in this case
            if ((config_flags%polar) .and. ((j == jds) .or. (j == jde))) then
              v_mvy(i, k, j) = 0.0
            else
              v_mvy(i, k, j) = v(i, k, j) / msfvy(i, j)
            end if
          end do
        end do
      end do

        ! Interpolate v_mvy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_mvy_at_full(i, k, j) = fnm(k) * v_mvy(i, k, j) + fnp(k) * v_mvy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_mvy_at_full(i, 1, j) =  cf1 * v_mvy(i, 1, j) + cf2 * v_mvy(i, 2, j) + cf3 * v_mvy(i, 3, j)
            ! Top face
          v_mvy_at_full(i, kte, j) = cft1 * v_mvy(i, ktes1, j) + cft2 * v_mvy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dv_dz_dz_dy = (partial v_mvy / partial z)|sigma (dv_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dv_dz_dz_dy(i, k, j) = (0.5 * (v_mvy(i, k, j) + v_mvy(i, k, j + 1)) - &
                0.5 * (v_mvy(i, k - 1, j) + v_mvy(i, k - 1, j + 1))) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial v / partial y)|z = msf_ty_vy * (partial v_mvy_at_full / partial y)|sigma - dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dv_dy(i, k, j) = msf_ty_vy(i, j) * ((v_mvy_at_full(i, k, j + 1) - v_mvy_at_full(i, k, j)) * rdy -  &
                dv_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          dv_dy(i, 1, j) = 0.0
            ! Top face
          dv_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dy


    subroutine Calc_dw_dy (w, zy, rdz, msfty, rdy, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy                  )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dy                             !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zy, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real,                                        intent (in)  :: rdy
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dy

        ! Local vars
      integer :: i, j, k, j_offset
      real    :: dz_dy
      real, dimension (its:ite, jts:jte) :: msf_ty_ty
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dy, w_mty, w_mty_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1) 

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_ty_ty(i, j) = msfty(i, j) * msfty(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        j_offset = 2
      else 
        j_offset = 1
      end if
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, kte
          do i = i_start, i_end
            w_mty(i, k, j) = w(i, k, j) / msfty(i, j)
          end do
        end do
      end do

        ! Average w_mty to the mass points (half levels)
      do j = j_start - 1, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            w_mty_at_half(i,k,j) = 0.5 * (w_mty(i, k, j) + w_mty(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dy)
      do j = j_start, j_end 
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dw_dz_dz_dy = (partial w_mty / partial z)|sigma (dw_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dw_dz_dz_dy(i, k, j) = (w_mty_at_half(i, k, j) - w_mty_at_half(i, k - 1, j)) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial w / partial y)|z = msf_ty_ty * (partial w_mty / partial y)|sigma - dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * (8.0 * (w_mty(i, k, j + 1) - w_mty(i, k, j - 1)) - (w_mty(i, k, j + 2) - w_mty(i, k, j - 2)) ) / 12.0 - dw_dz_dz_dy(i, k, j))
            else
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * 0.5 * (w_mty(i, k, j + 1) - w_mty(i, k, j - 1)) - dw_dz_dz_dy(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dw_dy(i, kts, j) = 0.0
          dw_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dy


    subroutine Calc_thetav (th, qv, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, thetav)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Calculate virtual potential temperature at mass points   !
      !                                                                   !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                    !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: th, qv
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: thetav

        ! Local vars
      integer :: i, j, k
      real :: sqv
      integer :: ktf, i_start, i_end, j_start, j_end

      ktf = Min (kte, kde - 1)

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Calculate virtual potential temperature
      do j = j_start - 1, j_end + 1
        do k = kts, ktf
          do i = i_start - 1, i_end + 1
            sqv = qv(i, k, j) / (1.0 + qv(i, k, j))
            thetav(i, k, j) = th(i, k, j) * (1.0 + P608 * sqv)
          end do
        end do
      end do

    end subroutine Calc_thetav


    subroutine Calc_mass_var_derivatives (thetav, qv, msftx, msfty, msfux, msfuy, msfvx, msfvy,  &
              rdx, rdy, fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn, dthetav_dx,         &
              dthetav_dy, dthetav_dz, dqv_dx, dqv_dy, dqv_dz, ids, ide, jds, jde, kds, kde,      &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte                         )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                     !
      ! Purpose: Calculate spatial derivatives of mass variables at the center of the faces !
      !                                                                                     !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                      !
      !                                                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: thetav, qv, zx, zy
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdz, rdzw
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dthetav_dx, dthetav_dy, dthetav_dz, &
                                                                   dqv_dx,     dqv_dy,     dqv_dz

      logical, parameter :: DEBUG = .false.

        ! Calculate thetav derivatives
      call Calc_dmassvar_dx (thetav, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dx)

      call Calc_dmassvar_dy (thetav, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dy)

      call Calc_dmassvar_dz (thetav, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dz)

        ! Calculate qv derivatives
      call Calc_dmassvar_dx (qv, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dx)

      call Calc_dmassvar_dy (qv, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dy)

      call Calc_dmassvar_dz (qv, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dz)

      if (DEBUG) print *, 'Done Calc_mass_var_derivatives...'

    end subroutine Calc_mass_var_derivatives


    subroutine Calc_scalar_var_derivatives (qc, qr, qi, qs, qg, qnc, qnr, qni, qnwfa, qnifa, msftx, msfty, msfux, msfuy, msfvx, msfvy, &
              rdx, rdy, fnm, fnp, cf1, cf2, cf3, zx, zy, rdz, rdzw, dnw, dn, dqc_dx, dqc_dy, dqc_dz, dqr_dx, dqr_dy, dqr_dz, &
              dqi_dx, dqi_dy, dqi_dz, dqs_dx, dqs_dy, dqs_dz, dqg_dx, dqg_dy, dqg_dz, dqnc_dx, dqnc_dy, dqnc_dz, &
              dqnr_dx, dqnr_dy, dqnr_dz, dqni_dx, dqni_dy, dqni_dz, dqnwfa_dx, dqnwfa_dy, dqnwfa_dz, dqnifa_dx, dqnifa_dy, dqnifa_dz, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                                     !
      ! Purpose: Calculate spatial derivatives of hydrometeor & scalar variables at the center of the faces !
      !                                                                                                     !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez                                                      !
      !                                                                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfux, msfuy, msfvx, msfvy, msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: zx, zy
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: qc, qr, qi, qs, qg, qnc, qnr, qni, qnwfa, qnifa
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdz, rdzw
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dqc_dx, dqc_dy, dqc_dz, &
         dqr_dx, dqr_dy, dqr_dz, dqi_dx, dqi_dy, dqi_dz, dqs_dx, dqs_dy, dqs_dz, dqg_dx, dqg_dy, dqg_dz, &
         dqnc_dx, dqnc_dy, dqnc_dz, dqnr_dx, dqnr_dy, dqnr_dz, dqni_dx, dqni_dy, dqni_dz, &
         dqnwfa_dx, dqnwfa_dy, dqnwfa_dz, dqnifa_dx, dqnifa_dy, dqnifa_dz

      logical, parameter :: DEBUG = .false.

        ! Calculate cloud mass derivatives
      if ( p_qc .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qc, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqc_dx)

        call Calc_dmassvar_dy (qc, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqc_dy)

        call Calc_dmassvar_dz (qc, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqc_dz)
      else
        dqc_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqc_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqc_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate rain mass derivatives
      if ( p_qr .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qr, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqr_dx)

        call Calc_dmassvar_dy (qr, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqr_dy)

        call Calc_dmassvar_dz (qr, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqr_dz)
      else
        dqr_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqr_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqr_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate ice mass derivatives
      if ( p_qi .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qi, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqi_dx)

        call Calc_dmassvar_dy (qi, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqi_dy)

        call Calc_dmassvar_dz (qi, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqi_dz)
      else
        dqi_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqi_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqi_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate snow mass derivatives
      if ( p_qs .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qs, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqs_dx)

        call Calc_dmassvar_dy (qs, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqs_dy)

        call Calc_dmassvar_dz (qs, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqs_dz)
      else
        dqs_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqs_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqs_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate graupel mass derivatives
      if ( p_qg .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qg, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqg_dx)

        call Calc_dmassvar_dy (qg, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqg_dy)

        call Calc_dmassvar_dz (qg, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqg_dz)
      else
        dqg_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqg_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqg_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate cloud number derivatives
      if ( p_qnc .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qnc, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnc_dx)

        call Calc_dmassvar_dy (qnc, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnc_dy)

        call Calc_dmassvar_dz (qnc, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnc_dz)
      else
        dqnc_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqnc_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqnc_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate rain number derivatives
      if ( p_qnr .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qnr, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnr_dx)

        call Calc_dmassvar_dy (qnr, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnr_dy)

        call Calc_dmassvar_dz (qnr, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnr_dz)
      else
        dqnr_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqnr_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqnr_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate ice number derivatives
      if ( p_qni .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qni, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqni_dx)

        call Calc_dmassvar_dy (qni, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqni_dy)

        call Calc_dmassvar_dz (qni, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqni_dz)
      else
        dqni_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqni_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqni_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate water friendly aerosol number derivatives
      if ( p_qnwfa .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qnwfa, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnwfa_dx)

        call Calc_dmassvar_dy (qnwfa, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnwfa_dy)

        call Calc_dmassvar_dz (qnwfa, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnwfa_dz)
      else
        dqnwfa_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqnwfa_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqnwfa_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

        ! Calculate ice friendly aerosol number derivatives
      if ( p_qnifa .ge. param_first_scalar ) then
        call Calc_dmassvar_dx (qnifa, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnifa_dx)

        call Calc_dmassvar_dy (qnifa, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
            ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnifa_dy)

        call Calc_dmassvar_dz (qnifa, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqnifa_dz)
      else
        dqnifa_dx(its:ite, kts:kte, jts:jte) = 0.0
        dqnifa_dy(its:ite, kts:kte, jts:jte) = 0.0
        dqnifa_dz(its:ite, kts:kte, jts:jte) = 0.0
      end if

      if (DEBUG) print *, 'Done Calc_scalar_var_derivatives...'

    end subroutine Calc_scalar_var_derivatives


    subroutine Calc_dmassvar_dx (var, zx, rdz, msftx, fnm, fnp, dn, dnw, rdx, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial x        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zx, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dx

        ! Local vars
      integer :: i, j, k
      real :: dz_dx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_full, dvar_dz_dz_dx
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i - 1, k, j) + var(i, k, j)) + &
                fnp(k) * (var(i - 1, k - 1, j) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i - 1, 1, j) +  &
              cf2 * var(i - 1, 2, j) + cf3 * var(i - 1, 3, j))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i - 1, ktes1, j) + (var(i - 1, ktes1, j) - &
              var(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at center of faces
            dz_dx = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
              ! dvar_dz_dz_dx = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dvar_dz_dz_dx(i, k, j) =  (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial var / partial x)|z = msftx * (partial var_at_full / partial x)|sigma - dvar_dz_dz_dx)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dx(i, k, j) = msftx(i, j) * (var_at_full(i + 1, k, j) -  var_at_full(i, k, j)) * rdx - dvar_dz_dz_dx(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dx(i, 1, j) = 0.0
          dvar_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dx


    subroutine Calc_dmassvar_dy (var, zy, rdz, msfty, fnm, fnp, dn, dnw, rdy, cf1, cf2, cf3,  &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial y        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, zy, rdz
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real, dimension (kms:kme),                   intent (in)  :: fnm, fnp, dn, dnw
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dy

        ! Local vars
      integer :: i, j, k
      real :: dz_dy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: var_at_full, dvar_dz_dz_dy
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      ktes1 = kte - 1
      ktes2 = kte - 2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1,   ktf
          do i = i_start, i_end
            var_at_full(i, k, j) = 0.5 * (fnm(k) * (var(i, k, j - 1) + var(i, k, j)) + &
                fnp(k) * (var(i, k - 1, j - 1) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          var_at_full(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
              cf3 * var(i, 3, j) + cf1 * var(i, 1, j - 1) +  &
              cf2 * var(i, 2, j - 1) + cf3 * var(i, 3, j - 1))
          var_at_full(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
              var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              var(i, ktes1, j - 1) + (var(i, ktes1, j - 1) - &
              var(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dvar_dz_dz_dy)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at center of faces
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
              ! dvar_dz_dz_dy = (partial var / partial z)|sigma (dvar_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dvar_dz_dz_dy(i, k, j) = (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial var / partial y)|z = msfty * (partial var_at_full / partial y)|sigma - dvar_dz_dz_dy)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dy(i, k, j) = msfty(i, j) * (var_at_full(i, k, j + 1) - var_at_full(i, k, j)) * rdy - dvar_dz_dz_dy(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dy(i, 1, j) = 0.0
          dvar_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dy


    subroutine Calc_dmassvar_dz (var, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial z        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: var, rdz
      integer,                                     intent (in)  :: ide, jde, kde,                  &
                                                                   ims, ime, jms, jme, kms, kme,   &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dvar_dz

        ! Local vars
      integer :: i, j, k
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktf = Min (kte, kde - 1)

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Interpolate mass var from half levels to full levels (center of faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dvar_dz(i, k, j) = (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dvar_dz(i, kts, j) = 0.0
          dvar_dz(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dz


    subroutine Calc_kh (du_dz, dv_dz, dthetav_dz, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, &
              turb_flux_vtheta_v, turb_flux_wtheta_v, l_master, q_sq, khx, khy, khz, ids, ide, &
              jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                   !
      ! Purpose: Diagnose 3D eddy diffusivities for tracer mixing         !
      !                                                                   !
      ! Author: Timothy W. Juliano                                        !
      !                                                                   !
      ! Comments: kh* are calculated at the faces                         !
      !                                                                   !
      ! References: MY82                                                  !
      !                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: du_dz, dv_dz, dthetav_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_uw, turb_flux_vw, &
                                                                   turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: khx, khy, khz
      integer, intent (in)    :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

        ! local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf

        ! Define dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          khx(i, kts, j) = 0.0
          khy(i, kts, j) = 0.0
          khz(i, kts, j) = 0.0
        end do
      end do

        ! MY82, Eqs. 28 and 31
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
              ! Khx
            if (dthetav_dz(i, k, j) == 0.0 .or. du_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khx(i, k, j) = 0.0
            else
              khx(i, k, j) = Max( ( (turb_flux_utheta_v(i, k, j) / (3.0 * a_2 * l_master(i, k, j) / sqrt(q_sq(i, k, j)))) &
                            + (turb_flux_uw(i, k, j) * dthetav_dz(i, k, j)) ) / (dthetav_dz(i, k, j) * du_dz(i, k, j)), 0.0 )
            end if
              ! Khy
            if (dthetav_dz(i, k, j) == 0.0 .or. dv_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khy(i, k, j) = 0.0
            else
              khy(i, k, j) = Max( ( (turb_flux_vtheta_v(i, k, j) / (3.0 * a_2 * l_master(i, k, j) / sqrt(q_sq(i, k, j)))) &
                            + (turb_flux_vw(i, k, j) * dthetav_dz(i, k, j)) ) / (dthetav_dz(i, k, j) * dv_dz(i, k, j)), 0.0 )
            end if
              ! Khz
            if (dthetav_dz(i, k, j) == 0.0 .or. q_sq(i, k, j) == Q_SQ_MIN) then
              khz(i, k, j) = 0.0
            else
              khz(i, k, j) = Max( -turb_flux_wtheta_v(i, k, j) / dthetav_dz(i, k, j), 0.0 )
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          khx(i, ktf + 1, j) = 0.0
          khy(i, ktf + 1, j) = 0.0
          khz(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Calc_kh


    subroutine Vertical_turb_mix (dth, du, dv, dw, dqv, dqc, dqr, dqi, dqs, dqg, dqnc, dqnr, dqni, dqnwfa, dqnifa, &
              khz, tracer, tracer_tend, pbl3d_opt, ss, config_flags, turb_flux_wtheta, turb_flux_uw, turb_flux_vw, &
              turb_flux_w2, turb_flux_wqv, turb_flux_wqc, turb_flux_wqr, turb_flux_wqi, turb_flux_wqs, turb_flux_wqg, &
              turb_flux_wqnc, turb_flux_wqnr, turb_flux_wqni, turb_flux_wqnwfa, turb_flux_wqnifa, rdz, dnw, dn, rho, &
              fnm, fnp, n_tracer, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                               !
      ! Purpose: Performs vertical turbulent mixing   !
      !                                               !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano !
      !                                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: pbl3d_opt
      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ss, n_tracer
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, dn, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, khz, rho, &
          turb_flux_wtheta, turb_flux_uw, turb_flux_vw, turb_flux_w2, turb_flux_wqv, &
          turb_flux_wqc, turb_flux_wqr, turb_flux_wqi, turb_flux_wqs, turb_flux_wqg, &
          turb_flux_wqnc, turb_flux_wqnr, turb_flux_wqni, turb_flux_wqnwfa, turb_flux_wqnifa
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: dth, &
          du, dv, dw, dqv
      real, dimension (ims:ime, kms:kme, jms:jme), intent(inout) :: dqc, dqr, dqi, dqs, dqg, &
          dqnc, dqnr, dqni, dqnwfa, dqnifa
      real, dimension (ims:ime, kms:kme, jms:jme, n_tracer), intent(in) ::  tracer
      real, dimension (ims:ime, kms:kme, jms:jme, n_tracer), intent(inout) ::  tracer_tend
      
        ! Local vars
      integer :: im
      logical, parameter :: ADD_TEND_WIND=.true., ADD_TEND_T=.true., &
          ADD_TEND_Q=.true.
      REAL, DIMENSION (its:ite, kts:kte, jts:jte) :: hydromet_tend
      
      if (ADD_TEND_WIND) then
        call Vertical_mixing_u (du, config_flags, turb_flux_uw, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)

        call Vertical_mixing_v (dv, config_flags, turb_flux_vw, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)

        if (pbl3d_opt > 0) then
          call Vertical_mixing_w (dw, config_flags, turb_flux_w2,   &
              dn, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
        end if

      end if

      if (ADD_TEND_T) then
        call Vertical_mixing_s_with_flux (dth, config_flags, turb_flux_wtheta, &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
            jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      if (ADD_TEND_Q) then
        call Vertical_mixing_s_with_flux (dqv, config_flags, turb_flux_wqv,  &
            dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
            jme, kms, kme, its, ite, jts, jte, kts, kte)
      end if

      if ( config_flags%pbl3d_scalar_mix .eq. 1 .and. ss .eq. config_flags%pbl3d_nsteps ) then
          ! Cloud mass
        if ( p_qc .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqc,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqc(its:ite, kts:kte, jts:jte) = dqc(its:ite, kts:kte, jts:jte) + &
                                           hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Rain mass
        if ( p_qr .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqr,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqr(its:ite, kts:kte, jts:jte) = dqr(its:ite, kts:kte, jts:jte) + &
                                           hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice mass
        if ( p_qi .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqi,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqi(its:ite, kts:kte, jts:jte) = dqi(its:ite, kts:kte, jts:jte) + &
                                           hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Snow mass
        if ( p_qs .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqs,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqs(its:ite, kts:kte, jts:jte) = dqs(its:ite, kts:kte, jts:jte) + &
                                           hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
           ! Graupel mass
        if ( p_qg .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqg,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqg(its:ite, kts:kte, jts:jte) = dqg(its:ite, kts:kte, jts:jte) + &
                                           hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Cloud number
        if ( p_qnc .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqnc,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqnc(its:ite, kts:kte, jts:jte) = dqnc(its:ite, kts:kte, jts:jte) + &
                                            hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Rain number
        if ( p_qnr .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqnr,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqnr(its:ite, kts:kte, jts:jte) = dqnr(its:ite, kts:kte, jts:jte) + &
                                            hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice number
        if ( p_qni .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqni,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqni(its:ite, kts:kte, jts:jte) = dqni(its:ite, kts:kte, jts:jte) + &
                                            hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Water friendly aerosol number
        if ( p_qnwfa .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqnwfa,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqnwfa(its:ite, kts:kte, jts:jte) = dqnwfa(its:ite, kts:kte, jts:jte) + &
                                              hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice friendly aerosol number
        if ( p_qnifa .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          call Vertical_mixing_s_with_flux (hydromet_tend, config_flags, turb_flux_wqnifa,  &
              dnw, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, &
              jme, kms, kme, its, ite, jts, jte, kts, kte)
          dqnifa(its:ite, kts:kte, jts:jte) = dqnifa(its:ite, kts:kte, jts:jte) + &
                                              hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
      end if

      if (ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar) then
        do im = param_first_scalar, n_tracer
          call Vertical_mixing_s (tracer_tend(ims, kms, jms, im), config_flags, tracer(ims, kms, jms, im), &
              khz, dnw, rdz, rho, fnm, fnp, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end do
      end if

    end subroutine Vertical_turb_mix


    subroutine Vertical_mixing_u (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  & 
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                       !
      ! Purpose: Calculates tendencies of resolved variables  !
      !          at u points due to vertical turbulent        !
      !          mixing using turbulent fluxes.               !
      !                                                       !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer, intent(in ) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency
      real, dimension (kms:kme), intent(in) ::    dnw, fnm, fnp

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_avg
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_u
      real (kind = kind(turb_flux_u)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min (jde - 2, jte)

        ! Calc flux at u point in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i - 1, kts, j) + rho(i, kts, j)) 
          turb_flux_u(i, kts, j) = rho_avg * 0.5 * (turb_flux(i - 1, kts, j) + turb_flux(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf 
          do i = i_start, i_end
            rho_avg = 0.5 * (fnm(k) * (rho(i, k, j) + rho(i - 1, k, j)) + &
                fnp(k) * (rho(i, k - 1, j) + rho(i - 1, k - 1, j)))
            turb_flux_u(i, k, j) = rho_avg * 0.5 * (turb_flux(i - 1, k, j) + turb_flux(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i - 1, ktf, j) + rho(i, ktf, j))
          turb_flux_u(i, ktf + 1, j) =  rho_avg * 0.5 * (turb_flux(i - 1, ktf + 1, j) + turb_flux(i, ktf + 1, j)) 
        end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G / dnw(k) * (turb_flux_u(i, k + 1, j) - turb_flux_u(i, k, j)) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_u


    subroutine Vertical_mixing_v (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                       !
      ! Purpose: Calculates tendencies of resolved variables  !
      !          at v points due to vertical turbulent        !
      !          mixing using turbulent fluxes.               !
      !                                                       !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags
      integer, intent (in) :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte
      real, dimension (ims:ime, kms:kme, jms:jme ), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) ::tendency
      real, dimension (kms:kme), intent(in) :: dnw, fnm, fnp

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_avg
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_v
      real (kind = kind(turb_flux_v)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end   = Min (jde - 1, jte)

        ! Calc flux at v points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i, kts, j - 1) + rho(i, kts, j))
          turb_flux_v(i, kts, j) = rho_avg * 0.5 * (turb_flux(i, kts, j - 1) + turb_flux(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            rho_avg = 0.5 * (fnm(k) * (rho(i, k, j - 1) + rho(i, k, j)) + &
                fnp(k) * (rho(i, k - 1, j - 1) + rho(i, k - 1, j)))
            turb_flux_v(i, k, j) = rho_avg * 0.5 * (turb_flux(i, k, j - 1) + turb_flux(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          rho_avg = 0.5 * (rho(i, ktf, j - 1) + rho(i, ktf, j))
          turb_flux_v(i, ktf + 1, j) =  rho_avg * 0.5 * (turb_flux(i, ktf + 1, j - 1) + turb_flux(i, ktf + 1, j))
        end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G / dnw(k) * (turb_flux_v(i, k + 1, j) - turb_flux_v(i, k, j)) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_v


    subroutine Vertical_mixing_w (tendency, config_flags, turb_flux, dn, rho,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates vertical velocity tendencies due  !
        !          to vertical turbulent mixing                 !
        !                                                       !
        ! Author: Pedro A. Jimenez                              !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, &
          ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent(in) :: dn
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real, dimension(its:ite, kts:kte, jts:jte) :: turb_flux_mass
      real (kind = kind(turb_flux_mass)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some dimensions
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end = Min (jde - 2, jte)

        ! Calc flux at mass points
      do j = j_start, j_end
        do k = kts + 1, ktf - 1
          do i = i_start, i_end
            turb_flux_mass(i, k, j) = 0.5 * (turb_flux(i, k + 1, j) + turb_flux(i, k, j)) * rho(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
          do i = i_start, i_end
            turb_flux_mass(i, kts, j) = 0.5 * (turb_flux(i, kts, j) + turb_flux(i, kts + 1, j)) * rho(i, kts, j)
            turb_flux_mass(i, ktf, j) = 0.5 * (turb_flux(i, ktf, j) + turb_flux(i, ktf + 1, j)) * rho(i, ktf, j)
          end do
      end do

        ! Calc tendency
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tendency(i, k, j) = ( G * (turb_flux_mass(i, k, j) - turb_flux_mass(i, k - 1, j)) / dn(k) ) / N_PBL3D
          end do
        end do
      end do

    end subroutine Vertical_mixing_w


    subroutine Vertical_mixing_s (tendency, config_flags, var, khz, dnw, rdz, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !          Code adapted from module_diffusion_em.F      !
        !                                                       !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: var, khz, rdz, rho
      real, dimension (ims:ime, kms:kme, jms:jme), intent(inout) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real, dimension(its:ite, kts:kte, jts:jte) :: H3, khz_avg

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            khz_avg(i, k, j) = khz(i, k, j) * (fnm(k) * rho(i, k, j) + fnp(k) * rho(i, k - 1, j))
            H3(i, k, j) = -khz_avg(i, k, j) * (var(i, k, j) - var(i, k - 1, j)) * rdz(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          H3(i, kts, j) = 0.
          H3(i, ktf + 1, j) = 0.
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tendency(i, k, j) = tendency(i, k, j) + ( G * (H3(i, k + 1, j) - H3(i, k, j)) / dnw(k) )
          end do
        end do
      end do

    end subroutine Vertical_mixing_s


    subroutine Vertical_mixing_s_with_flux (tendency, config_flags, turb_flux, dnw, rho, fnm, fnp,  &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !                                                       !
        ! Authors: Pedro A. Jimenez & Timothy W. Juliano        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                   intent(in) :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: turb_flux, rho
      real, dimension (its:ite, kts:kte, jts:jte), intent(out) :: tendency

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      real :: rho_kp1, rho_k
      real (kind = kind(turb_flux)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

        ! First layer
      do j = j_start, j_end
        do i = i_start, i_end
          rho_k = rho(i, kts, j)
          rho_kp1 = fnm(kts + 1) * rho(i, kts + 1, j) + fnp(kts + 1) * rho(i, kts, j)
          tendency(i, kts, j) = ( G * (rho_kp1 * turb_flux(i, kts + 1, j) - &
              rho_k * turb_flux(i, kts, j)) / dnw(kts) ) / N_PBL3D
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf - 1
          do i = i_start, i_end
            rho_k = fnm(k) * rho(i, k, j) + fnp(k) * rho(i, k - 1, j)
            rho_kp1 = fnm(k + 1) * rho(i, k + 1, j) + fnp(k + 1) * rho(i, k, j)
            tendency(i, k, j) = ( G * (rho_kp1 * turb_flux(i, k + 1, j) - &
                rho_k * turb_flux(i, k, j)) / dnw(k) ) / N_PBL3D
          end do
        end do
      end do

        ! Last layer
      do j = j_start, j_end
        do i = i_start, i_end
          rho_k = fnm(ktf) * rho(i, ktf, j) + fnp(ktf) * rho(i, ktf - 1, j)
          rho_kp1 = rho(i, ktf, j)
          tendency(i, ktf, j) = ( G * (rho_kp1 * turb_flux(i, ktf + 1, j) - &
              rho_k * turb_flux(i, ktf, j)) / dnw(ktf) ) / N_PBL3D
        end do
      end do

    end subroutine Vertical_mixing_s_with_flux


    subroutine Horizontal_turb_mix (rt_tendf, ru_tendf, rv_tendf, rw_tendf, moist_tendf, qc_tendf, qr_tendf, &
              qi_tendf, qs_tendf, qg_tendf, qnc_tendf, qnr_tendf, qni_tendf, qnwfa_tendf, qnifa_tendf, tracer_tend, &
              ss, config_flags, turb_flux_u2, turb_flux_v2, turb_flux_uv, &
              turb_flux_uw, turb_flux_vw, turb_flux_w2, turb_flux_utheta, turb_flux_vtheta, turb_flux_uqv, &
              turb_flux_vqv, turb_flux_uqc, turb_flux_vqc, turb_flux_uqr, turb_flux_vqr, turb_flux_uqi, turb_flux_vqi, &
              turb_flux_uqs, turb_flux_vqs, turb_flux_uqg, turb_flux_vqg, turb_flux_uqnc, turb_flux_vqnc, &
              turb_flux_uqnr, turb_flux_vqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_uqnwfa, turb_flux_vqnwfa, &
              turb_flux_uqnifa, turb_flux_vqnifa, khx, khy, n_tracer, tracer, &
              msfux, msfuy, msfvx, msfvy, msftx, msfty, dx, dy, rdx, rdy, rdz, rdzw, fnm, fnp, cf1, cf2, cf3, &
              zx, zy, dn, dnw, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                               !
      ! Purpose: Performs horizontal turbulent mixing !
      !                                               !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano !
      !                                               !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags
 
      INTEGER ,        INTENT(IN) :: ss 
      INTEGER ,        INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

      REAL ,           INTENT(IN) :: cf1, cf2, cf3, dx, dy

      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: fnm
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: fnp
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: dnw
      REAL , DIMENSION( kms:kme ) ,                 INTENT(IN) :: dn

      REAL , DIMENSION( ims:ime, jms:jme) ,         INTENT(IN) ::  msfux, &
                                                                   msfuy, &
                                                                   msfvx, &
                                                                   msfvy, &
                                                                   msftx, &
                                                                   msfty

      REAL , DIMENSION( its:ite, kts:kte, jts:jte), INTENT(INOUT) :: rt_tendf, &
                                                                     ru_tendf, &
                                                                     rv_tendf, &
                                                                     rw_tendf, &
                                                                     moist_tendf
      REAL , DIMENSION ( ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: qc_tendf, &
           qr_tendf, qi_tendf, qs_tendf, qg_tendf, qnc_tendf, qnr_tendf, qni_tendf, qnwfa_tendf, qnifa_tendf

      REAL , DIMENSION( ims:ime, kms:kme, jms:jme, n_tracer), INTENT(INOUT) :: tracer_tend
      REAL , DIMENSION( ims:ime, kms:kme, jms:jme, n_tracer), INTENT(IN) :: tracer
      INTEGER                                               , INTENT(IN) :: n_tracer

      REAL , DIMENSION( ims:ime, kms:kme, jms:jme), INTENT(IN) :: zx, zy, rdz, rdzw, &
                                                                  rho, khx, khy
                                                                  

      REAL, DIMENSION (ims:ime, kms:kme, jms:jme), intent(IN) :: turb_flux_u2,  &
          turb_flux_v2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_w2, &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_uqv, turb_flux_vqv, &
          turb_flux_uqc, turb_flux_vqc, turb_flux_uqr, turb_flux_vqr, turb_flux_uqi, turb_flux_vqi, &
          turb_flux_uqs, turb_flux_vqs, turb_flux_uqg, turb_flux_vqg, turb_flux_uqnc, turb_flux_vqnc, &
          turb_flux_uqnr, turb_flux_vqnr, turb_flux_uqni, turb_flux_vqni, turb_flux_uqnwfa, turb_flux_vqnwfa, &
          turb_flux_uqnifa, turb_flux_vqnifa

      REAL, INTENT(IN) :: rdx, rdy

        ! Local
      INTEGER :: im
      REAL, DIMENSION (its-1:ite, kts:kte, jts-1:jte) :: sf_alpha
      REAL, DIMENSION (its:ite, kts:kte, jts:jte) :: hydromet_tend

      call Calc_slope_factor (config_flags, msftx, msfty, dx, dy, zx, zy, rdzw, sf_alpha, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_u (ru_tendf, config_flags, turb_flux_u2, turb_flux_uv,  &
          msfux, msfuy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, rho, sf_alpha, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_v (rv_tendf, config_flags, turb_flux_v2, turb_flux_uv,  &
          msfvx, msfvy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, rho, sf_alpha, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_w (rw_tendf, config_flags, turb_flux_uw, turb_flux_vw, &
          msftx, msfty, rdx, rdy, fnm, fnp, dn, zx, zy, rdz, rho, sf_alpha, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_s_with_flux (rt_tendf, config_flags, msftx, msfty, rdx, rdy, &
          fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_utheta, turb_flux_vtheta, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Horizontal_mixing_s_with_flux (moist_tendf, config_flags, msftx, msfty, rdx, rdy, &
          fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqv, turb_flux_vqv, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      if ( config_flags%pbl3d_scalar_mix .eq. 1 .and. ss .eq. config_flags%pbl3d_nsteps ) then
          ! Cloud mass
        if ( p_qc .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qc_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqc, turb_flux_vqc, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qc_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Rain mass
        if ( p_qr .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qr_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqr, turb_flux_vqr, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qr_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice mass
        if ( p_qi .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qi_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqi, turb_flux_vqi, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qi_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Snow mass
        if ( p_qs .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qs_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqs, turb_flux_vqs, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qs_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Graupel mass
        if ( p_qg .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qg_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqg, turb_flux_vqg, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qg_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Cloud number
        if ( p_qnc .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qnc_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqnc, turb_flux_vqnc, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qnc_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
           ! Rain number
        if ( p_qnr .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qnr_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqnr, turb_flux_vqnr, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qnr_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice number
        if ( p_qni .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qni_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqni, turb_flux_vqni, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qni_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Water friendly aerosol number
        if ( p_qnwfa .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qnwfa_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqnwfa, turb_flux_vqnwfa, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qnwfa_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
          ! Ice friendly aerosol number
        if ( p_qnifa .ge. param_first_scalar ) then
          hydromet_tend(its:ite, kts:kte, jts:jte) = 0.0
          hydromet_tend(its:ite, kts:kte, jts:jte) = qnifa_tendf(its:ite, kts:kte, jts:jte)
          call Horizontal_mixing_s_with_flux (hydromet_tend, config_flags, msftx, msfty, rdx, rdy, &
              fnm, fnp, zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_uqnifa, turb_flux_vqnifa, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
          qnifa_tendf(its:ite, kts:kte, jts:jte) = hydromet_tend(its:ite, kts:kte, jts:jte)
        end if
      end if

      if (ss .eq. config_flags%pbl3d_nsteps .and. n_tracer .ge. param_first_scalar) then
        do im = param_first_scalar, n_tracer
          call Horizontal_mixing_s (tracer_tend(ims, kms, jms, im), config_flags, tracer(ims, kms, jms, im), &
              khx, khy, dn, dnw, rdzw, rho, fnm, fnp, cf1, cf2, cf3, rdx, rdy, zx, zy, msfux, msfvy, msftx, msfty, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        end do
      end if

    end subroutine Horizontal_turb_mix


    subroutine Calc_slope_factor (config_flags, msftx, msfty, dx, dy, zx, zy, rdzw, sf_alpha, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Calculate the slope factor to account for steep terrain   !
      !                                                                    !
      ! Author: Timothy W. Juliano                                         !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           Calculation is done at mass points                       !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent(in)  :: config_flags

      REAL, DIMENSION (ims:ime, jms:jme),          INTENT(IN)  :: msftx, msfty
      REAL,                                        INTENT(IN)  :: dx, dy
      REAL, DIMENSION (ims:ime, kms:kme, jms:jme), INTENT(IN)  :: zx, zy, rdzw
      REAL, DIMENSION (its-1:ite, kts:kte, jts-1:jte), INTENT(OUT) :: sf_alpha
      INTEGER,                                     INTENT(IN)  :: ide, jde, kde,                &
                                                                  ims, ime, jms, jme, kms, kme, &
                                                                  its, ite, jts, jte, kts, kte

        ! Local
      INTEGER :: i, j, k, i_start, i_end, j_start, j_end, ktf
      REAL    :: dxm, dym, tmpzx, tmpzy
      logical, parameter :: SLOPE_FACTOR = .true.

        ! Define dimensions
      i_start = its - 1
      i_end = ite
      j_start = jts - 1
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (SLOPE_FACTOR .and. config_flags%pbl3d_opt > 0) then
        do j = j_start, j_end
          do k = kts, ktf
            do i = i_start, i_end
              dxm = dx / msftx(i, j)
              dym = dy / msfty(i, j)
              tmpzx = (0.25 * ( abs(zx(i, k, j))+ abs(zx(i + 1, k, j))  &
                              + abs(zx(i, k + 1, j)) + abs(zx(i + 1, k + 1, j))) * rdzw(i, k, j) * dxm)
              tmpzy = (0.25 * ( abs(zy(i, k, j))+ abs(zy(i, k, j + 1))  &
                              + abs(zy(i, k + 1, j)) + abs(zy(i, k + 1, j + 1))) * rdzw(i, k, j) * dym)
              sf_alpha(i, k, j) = max(sqrt(tmpzx * tmpzx + tmpzy * tmpzy), 1.0)
            end do
          end do
        end do
      else
        do j = j_start, j_end
          do k = kts, ktf
            do i = i_start, i_end
              sf_alpha(i, k, j) = 1.0
            end do
          end do
        end do
      end if

    end subroutine Calc_slope_factor


    subroutine Horizontal_mixing_u (tendency, config_flags, turb_flux_u2,  turb_flux_uv, &
        msfux, msfuy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, rho, sf_alpha, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the u tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             du / dt = - [d (u'u') / dx + d (u'v') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfux, msfuy 
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_u2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local data
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      integer :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_u2_at_mass, tf_uv_at_corners

        ! Calc some indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max(ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min(ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max(jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min(jde - 2, jte)

        ! tf_u2_at_mass = rho * turb_flux_u2 (move flux from center of faces to mass points)
      call Prep_turb_flux_u2 (config_flags, tf_u2_at_mass, turb_flux_u2, rho, ids, ide, &
          jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! tf_uv_at_corners = rho * turb_flux_uv (move flux from center of faces to half levels at corners)
      ie_ext = 0
      je_ext = 1
      call Prep_turb_flux_uv (config_flags, tf_uv_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calculate the u tendency due to horizontal turbulent mixing
      call Calc_htend_du (config_flags, tendency, turb_flux_u2, turb_flux_uv, tf_u2_at_mass, tf_uv_at_corners, &
          msfux, msfuy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_u


    subroutine Horizontal_mixing_v (tendency, config_flags, turb_flux_v2, turb_flux_uv, &
        msfvx, msfvy, rdx, rdy, fnm, fnp, dnw, zx, zy, rdzw, rho, sf_alpha, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the v tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             dv / dt = - [d (u'v') / dx + d (v'v') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) :: msfvx, msfvy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_v2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local vars
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end
      integer :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) ::  &
          tf_uv_at_corners, tf_v2_at_mass

        ! Indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

        ! tf_uv_at_corners = rho * turb_flux_uv (move flux from center of faces to half levels at corners)
      ie_ext = 1
      je_ext = 0
      call Prep_turb_flux_uv (config_flags, tf_uv_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! tf_v2_at_mass = rho * turb_flux_v2 (move flux from center of faces to mass points)
      call Prep_turb_flux_v2 (config_flags, tf_v2_at_mass, turb_flux_v2, rho, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calculate the v tendency due to horizontal turbulent mixing
      call Calc_htend_dv (config_flags, tendency, turb_flux_v2, turb_flux_uv, tf_v2_at_mass, tf_uv_at_corners, &
          msfvx, msfvy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_v


    subroutine Horizontal_mixing_w (tendency, config_flags, turb_flux_uw, turb_flux_vw, &
        msftx, msfty, rdx, rdy, fnm, fnp, dn, zx, zy, rdz, rho, sf_alpha, &
        ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the w tendency due to horizontal turb mix             !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             dw / dt = - [d (u'w') / dx + d (v'w') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent(in) :: fnm, fnp, dn
      real, dimension (ims:ime, jms:jme), intent(in) ::  msftx, msfty

      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, rho, &
          zx, zy, turb_flux_uw, turb_flux_vw
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha

        ! Local data
      integer :: i, j, k, ktf
      integer :: i_start, i_end, j_start, j_end, n_offset

      real, dimension(its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: &
          tf_uw_at_full, tf_vw_at_full

        ! Calc indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
          config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
          config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
          config_flags%nested) j_start = Max (jds + 1,jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
          config_flags%nested) j_end = Min (jde - 2,jte)

      if (.not. USE_HIGH_ORDER) then
          ! tf_uw_at_full = rho * turb_flux_uw (move flux from center of faces to u points in faces)
        call Prep_turb_flux_uw (config_flags, tf_uw_at_full, turb_flux_uw, fnm, fnp, rho, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

          ! tf_vw_at_full = rho * turb_flux_vw (move flux from center of faces to v points in faces)
        call Prep_turb_flux_vw (config_flags, tf_vw_at_full, turb_flux_vw, fnm, fnp, rho, &
            ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
      else
        n_offset = 2
        do j = j_start - n_offset, j_end + n_offset
          do k = kts + 1, ktf
            do i = i_start - n_offset, i_end + n_offset
                ! tf_uw_at_full = rho * turb_flux_uw at the center of the faces
              tf_uw_at_full(i, k, j) = turb_flux_uw(i, k, j) * (fnm(k) * rho(i, k, j) + &
                  fnp(k) * rho(i, k - 1, j))
                ! tf_vw_at_full = rho * turb_flux_vw at the center of the faces
              tf_vw_at_full(i, k, j) = turb_flux_vw(i, k, j) * (fnm(k) * rho(i, k, j) + &
                  fnp(k) * rho(i, k - 1, j))
            end do
          end do
        end do
      end if

        ! Calculate the w tendency due to horizontal turbulent mixing
      call Calc_htend_dw (config_flags, tendency, turb_flux_uw, turb_flux_vw, tf_uw_at_full, tf_vw_at_full, &
          msftx, msfty, rdx, rdy, dn, zx, zy, rdz, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_w


    subroutine Horizontal_mixing_s_with_flux (tendency, config_flags, msftx, msfty, rdx, rdy, fnm, fnp,   &
              zx, zy, rdz, rdzw, dnw, rho, sf_alpha, turb_flux_x, turb_flux_y, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                    !
      ! Purpose: Add the scalar tendency due to horizontal turb mix        !
      !                                                                    !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano                      !
      !                                                                    !
      ! Comments: Adapted subroutine from the diffusion module             !
      !           The tendency is:                                         !
      !             ds / dt = - [d (u's') / dx + d (v's') / dy]            !
      !                                                                    !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, intent (in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) ::  &
          msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: rdz, &
          rdzw, zx, zy, turb_flux_x, turb_flux_y, rho
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: tendency

        ! Local vars
      integer :: ktf, i_start, i_end, j_start, j_end
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) ::    &
          tf_x_at_u, tf_y_at_v

        ! Define indexes
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde-1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
           config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
           config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
           config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
           config_flags%nested) j_end = Min (jde - 2, jte)

        ! Scalar turbulent flux in x direction (at u points)
      call Prep_turb_flux_sx (turb_flux_x, tf_x_at_u, rho, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        ! Scalar turbulent flux in y direction (at v points)
      call Prep_turb_flux_sy (turb_flux_y, tf_y_at_v, rho, ims, ime, jms, jme, kms, kme, &
          its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        ! Calculate scalar tendency due to horizontal turbulent mixing (at mass points)
      call Calc_htend_ds (config_flags, tendency, turb_flux_x, turb_flux_y, tf_x_at_u, tf_y_at_v, rho, sf_alpha, &
          fnm, fnp, zx, zy, dnw, rdzw, msfty, msftx, rdx, rdy, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, &
          kts, kte, i_start, i_end, j_start, j_end, ktf)

    end subroutine Horizontal_mixing_s_with_flux


    subroutine Horizontal_mixing_s (tendency, config_flags, var, khx, khy, dn, dnw, rdzw, rho, fnm, fnp,  &
              cf1, cf2, cf3, rdx, rdy, zx, zy, msfux, msfvy, msftx, msfty, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                       !
        ! Purpose: Calculates tendencies of resolved variables  !
        !          at mass points due to vertical turbulent     !
        !          mixing using turbulent fluxes.               !
        !          Code adapted from module_diffusion_em.F      !
        !                                                       !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez        !
        !                                                       !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent(in) :: config_flags
      integer,                     intent(in) :: ids, ide, jds, jde, kds, kde, &
                                                 ims, ime, jms, jme, kms, kme, &
                                                 its, ite, jts, jte, kts, kte

      real,                                        intent(in) :: cf1, cf2, cf3, rdx, rdy
      real, dimension (kms:kme),                   intent(in) :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent(in) :: msfux, msfvy, msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: var, khx, khy, zx, zy, rdzw, rho
      real, dimension (ims:ime, kms:kme, jms:jme), intent(inout) :: tendency

        ! Local vars
      integer :: i, j, k, ktf, ktes1, ktes2
      integer :: i_start, i_end, j_start, j_end
      real :: mrdx, mrdy, tmpzx, tmpzy, rdzu, rdzv
      real, dimension(its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: H1, H2, H1avg, H2avg, &
                                                                    khx_avg, khy_avg, zx_at_m, zy_at_m

        ! Define some dimensions
      i_start = its
      i_end   = Min (ite, ide - 1)
      j_start = jts
      j_end   = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)
      ktes1 = kte - 1
      ktes2 = kte - 2

      if (config_flags%open_xs .or. config_flags%specified .or. &
        config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
        config_flags%nested) i_end   = Min (ide - 2, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
        config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
        config_flags%nested) j_end   = Min (jde - 2, jte)

! H1 = partial var over partial x

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            khx_avg(i, k, j) = 0.25 * (khx(i - 1, k, j) + khx(i, k, j) + khx(i - 1, k + 1, j) + khx(i, k + 1, j)) * &
                                            0.5 * (rho(i - 1, k, j) + rho(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            H1avg(i, k, j) = 0.5 * (fnm(k) * (var(i - 1, k, j) + var(i, k, j)) + fnp(k) * (var(i - 1, k - 1, j) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          H1avg(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
                                    cf3 * var(i, 3, j) + cf1 * var(i - 1, 1, j) + &
                                    cf2 * var(i - 1, 2, j) + cf3 * var(i - 1, 3, j))
          H1avg(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
                                        var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
                                        var(i - 1, ktes1, j) + (var(i - 1, ktes1, j) - &
                                        var(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            tmpzx = 0.5 * (zx(i, k, j) + zx(i, k + 1, j))
            rdzu = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i - 1, k, j))
            H1(i, k, j) = -msfux(i, j) * khx_avg(i, k, j) * ( rdx * (var(i, k, j) - var(i - 1, k, j)) - tmpzx * &
                           (H1avg(i, k + 1, j) - H1avg(i, k, j)) * rdzu )
          end do
        end do
      end do

! H2 = partial var over partial y

      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            khy_avg(i, k, j) = 0.25 * (khy(i, k, j - 1) + khy(i, k, j) + khy(i, k + 1, j - 1) + khy(i, k + 1, j)) * &
                                            0.5 * (rho(i, k, j - 1) + rho(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            H2avg(i, k, j) = 0.5 * (fnm(k) * (var(i, k, j - 1) + var(i, k, j)) + fnp(k) * (var(i, k - 1, j - 1) + var(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          H2avg(i, kts, j) = 0.5 * (cf1 * var(i, 1, j) + cf2 * var(i, 2, j) + &
                                    cf3 * var(i, 3, j) + cf1 * var(i, 1, j - 1) + &
                                    cf2 * var(i, 2, j - 1) + cf3 * var(i, 3, j - 1))
          H2avg(i, ktf + 1, j) = 0.5 * (var(i, ktes1, j) + (var(i, ktes1, j) - &
                                        var(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
                                        var(i, ktes1, j - 1) + (var(i, ktes1, j - 1) - &
                                        var(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            tmpzy = 0.5 * (zy(i, k, j) + zy(i, k + 1, j))
            rdzv = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i, k, j - 1))
            H2(i, k, j) = -msfvy(i, j) * khy_avg(i, k, j) * ( rdy * (var(i, k, j) - var(i, k, j - 1)) - tmpzy * &
                           (H2avg(i, k + 1, j) - H2avg(i, k, j)) * rdzv )
          end do
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            H1avg(i, k, j) = 0.5 * (fnm(k) * (H1(i + 1, k, j) + H1(i, k, j)) + &
                                    fnp(k) * (H1(i + 1, k - 1, j) + H1(i, k - 1, j)))
            H2avg(i, k, j) = 0.5 * (fnm(k) * (H2(i, k, j + 1) + H2(i, k, j)) + &
                                    fnp(k) * (H2(i, k - 1, j + 1) + H2(i, k - 1, j)))
            zx_at_m(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_m(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          H1avg(i, kts, j) = 0.
          H1avg(i, ktf + 1, j) = 0.
          H2avg(i, kts, j) = 0.
          H2avg(i, ktf + 1, j) = 0.
          zx_at_m(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts, j) + zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_m(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy

            tendency(i, k, j) = tendency(i, k, j) + g / (dnw(k) * rdzw(i, k, j)) * &
                 ( mrdx * (H1(i + 1, k, j) - H1(i, k, j)) + mrdy * (H2(i, k, j + 1) - H2(i, k, j)) - &
                   msftx(i, j) * zx_at_m(i, k, j) * (H1avg(i, k + 1, j) - H1avg(i, k, j)) * rdzw(i, k, j) - &
                   msfty(i, j) * zy_at_m(i, k, j) * (H2avg(i, k + 1, j) - H2avg(i, k, j)) * rdzw(i, k, j) )
          end do
        end do
      end do

    end subroutine Horizontal_mixing_s


    subroutine Prep_turb_flux_u2 (config_flags, tf_at_mass, turb_flux_u2, rho, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, &
         ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: &
          turb_flux_u2, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

      i_start = i_start - 1

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_at_mass(i, k, j) = rho(i, k, j) * 0.5 * (turb_flux_u2(i, k + 1, j) + turb_flux_u2(i, k, j))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_u2


    subroutine Prep_turb_flux_v2 (config_flags, tf_at_mass, turb_flux_v2, rho, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, &
         ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: &
          turb_flux_v2, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1, jte)

      j_start = j_start - 1

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_at_mass(i, k, j) = rho(i, k, j) * 0.5 * (turb_flux_v2(i, k + 1, j) + turb_flux_v2(i, k, j))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_v2


    subroutine Prep_turb_flux_uv (config_flags, tf_at_corners, turb_flux_uv, rho, ie_ext, je_ext, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

      integer, intent (in) :: ie_ext, je_ext

      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), &
          intent (inout) :: tf_at_corners
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in)  &
          ::  turb_flux_uv, rho

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real  :: turb_flux_avg

        ! Calc indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .or. config_flags%specified .or. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .or. config_flags%specified .or. &
         config_flags%nested) i_end = Min (ide - 1 - ie_ext, ite)
      if (config_flags%open_ys .or. config_flags%specified .or. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .or. config_flags%specified .or. &
         config_flags%nested) j_end = Min (jde - 1 - je_ext, jte)

      i_end = i_end + ie_ext
      j_end = j_end + je_ext

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.125 * (turb_flux_uv(i - 1, k, j) + turb_flux_uv(i, k, j) +  &
                turb_flux_uv(i - 1, k, j - 1) + turb_flux_uv(i, k, j - 1) + turb_flux_uv(i - 1, k + 1, j) &
                + turb_flux_uv(i, k + 1, j) + turb_flux_uv(i - 1, k + 1, j - 1) + turb_flux_uv(i, k + 1, j - 1))
            tf_at_corners(i, k, j) = turb_flux_avg * 0.25 * (rho(i - 1, k, j) + rho(i, k, j) + &
                rho(i - 1, k, j - 1) + rho(i, k, j - 1))
          end do
        end do
      end do

    end subroutine Prep_turb_flux_uv


    subroutine Prep_turb_flux_uw (config_flags, tf_at_full, turb_flux_uw,  &
        fnm, fnp, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms,   &
        jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent (in) :: ids, ide, jds, jde, kds, kde, ims, ime, &
          jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent (in) :: fnm, fnp
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), & 
          intent (inout) :: tf_at_full
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in):: rho, &
          turb_flux_uw

        ! Local vars
      integer  :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real :: turb_flux_avg

        ! Define indexes
      i_start = its
      i_end = ite
      j_start = jts
      j_end = Min (jte, jde - 1 )
      ktf = Min (kte, kde - 1)

      if (config_flags%open_xs .OR. config_flags%specified .OR. &
         config_flags%nested) i_start = Max (ids + 1, its)
      if (config_flags%open_xe .OR. config_flags%specified .OR. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      if (config_flags%open_ys .OR. config_flags%specified .OR. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      if (config_flags%open_ye .OR. config_flags%specified .OR. &
         config_flags%nested) j_end = Min (jde - 2, jte)

      i_end = i_end + 1

       ! Bottom face
      do j = j_start, j_end
        do i = i_start, i_end
          turb_flux_avg = 0.5 * (turb_flux_uw(i, kts, j) +  &
              turb_flux_uw(i - 1, kts, j))
          tf_at_full(i, kts, j) = turb_flux_avg * 0.5 * &
              (rho(i - 1, kts, j) + rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.5 * (turb_flux_uw(i, k, j) +  &
                turb_flux_uw(i - 1, k, j)) 
            tf_at_full(i, k, j) = turb_flux_avg * 0.5 * (fnm(k) *  &
                (rho(i - 1, k, j) + rho(i, k, j) ) + fnp(k) * &
                (rho(i - 1, k - 1, j) + rho(i, k - 1, j)))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end
        do i = i_start, i_end
          tf_at_full(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Prep_turb_flux_uw


    subroutine Prep_turb_flux_vw (config_flags, tf_at_full, turb_flux_vw, &
        fnm, fnp, rho, ids, ide, jds, jde, kds, kde, ims, ime, jms,  &
        jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      type (grid_config_rec_type), intent(in)  :: config_flags

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, ims, ime, &
         jms, jme, kms, kme, its, ite, jts, jte, kts, kte

      real, dimension (kms:kme), intent(in) :: fnm, fnp
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), &
          intent (inout) :: tf_at_full
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rho, &
          turb_flux_vw

        ! Local vars
      integer :: i, j, k, ktf, i_start, i_end, j_start, j_end
      real :: turb_flux_avg

        ! Define indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = jte
      ktf = Min (kte, kde - 1)

      IF (config_flags%open_xs .OR. config_flags%specified .OR. &
         config_flags%nested) i_start = Max (ids + 1, its)
      IF (config_flags%open_xe .OR. config_flags%specified .OR. &
         config_flags%nested) i_end = Min (ide - 2, ite)
      IF (config_flags%open_ys .OR. config_flags%specified .OR. &
         config_flags%nested) j_start = Max (jds + 1, jts)
      IF (config_flags%open_ye .OR. config_flags%specified .OR. &
         config_flags%nested) j_end = Min (jde - 2, jte)

      j_end = j_end + 1

       ! Bottom face
      do j = j_start, j_end
        do i = i_start, i_end
          turb_flux_avg = 0.5 * (turb_flux_vw(i, kts, j) +  &
              turb_flux_vw(i, kts, j - 1))
          tf_at_full(i, kts, j) = turb_flux_avg * 0.5 * &
              (rho(i, kts, j) + rho(i, kts, j - 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            turb_flux_avg = 0.5 * (turb_flux_vw(i, k, j) + turb_flux_vw(i, k, j - 1)) 
            tf_at_full(i, k, j) = turb_flux_avg * 0.5 * (fnm(k) * (rho(i, k, j) + &
                rho(i, k, j - 1)) + fnp(k) * (rho(i, k - 1, j) + rho(i, k - 1, j - 1)))
          end do
        end do
      end do

        ! Top face
      do j = j_start, j_end
        do i = i_start, i_end
          tf_at_full(i, ktf + 1, j) = 0.0
        end do
      end do

    end subroutine Prep_turb_flux_vw


    subroutine Prep_turb_flux_sx (turb_flux_x, tf_x_at_u, rho, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_x, rho
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_x_at_u

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: turb_flux_average

        ! tf_x_at_u = turb flux in x direction (move flux from center of faces to u points)
      do j = j_start, j_end
        do k = kts, ktf - 1
          do i = i_start, i_end + 1
            turb_flux_average = 0.25 * (turb_flux_x(i - 1, k, j) + turb_flux_x(i, k, j) +  &
               turb_flux_x(i - 1, k + 1, j) + turb_flux_x(i, k + 1, j))
            tf_x_at_u(i, k, j) = turb_flux_average  * 0.5 * (rho(i - 1, k, j) + rho(i, k, j))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end
          do i = i_start, i_end + 1
            tf_x_at_u(i, ktf, j) = 0.25 * (turb_flux_x(i - 1, ktf, j) + turb_flux_x(i, ktf, j)) * 0.5 * (rho(i - 1, ktf, j) + rho(i, ktf, j))
          end do
      end do

    end subroutine Prep_turb_flux_sx


    subroutine Prep_turb_flux_sy (turb_flux_y, tf_y_at_v, rho, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_y, rho
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (inout) :: tf_y_at_v

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: turb_flux_average

        ! tf_y_at_v = turb flux in y direction (move flux from center of faces to v points)
      do j = j_start, j_end + 1
        do k = kts, ktf - 1
          do i = i_start, i_end
            turb_flux_average = 0.25 * (turb_flux_y(i, k, j - 1) + turb_flux_y(i, k, j) + &
                turb_flux_y(i, k + 1, j - 1) + turb_flux_y(i, k + 1, j))
            tf_y_at_v(i, k, j) = turb_flux_average * 0.5 * (rho(i, k, j - 1) + rho(i, k, j))
          end do
        end do
      end do

       ! Top face
      do j = j_start, j_end + 1
        do i = i_start, i_end
          tf_y_at_v(i, ktf, j) = 0.25 * (turb_flux_y(i, ktf, j - 1) + turb_flux_y(i, ktf, j)) * 0.5 * (rho(i, ktf, j - 1) + rho(i, ktf, j))
        end do
      end do

    end subroutine Prep_turb_flux_sy


    subroutine Calc_htend_du (config_flags, tendency, turb_flux_u2, turb_flux_uv, tf_u2_at_mass, tf_uv_at_corners, &
              msfux, msfuy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfux, msfuy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_u2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), intent(in) :: &
          tf_u2_at_mass, tf_uv_at_corners

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpdz, tmpalpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_u2_at_u, tf_uv_at_u, zy_at_u, zx_at_u
      real (kind = kind(tf_u2_at_u)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_u2_at_u = rho * turb_flux_u2 at u points in the faces
        ! tf_uv_at_u = rho * turb_flux_uv at u points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_u2_at_u(i, kts, j) = 0.5 * (turb_flux_u2(i - 1, kts, j) + &
              turb_flux_u2(i, kts, j)) * 0.5 * (rho(i - 1, kts, j) +   &
              rho(i, kts, j))
          tf_uv_at_u(i, kts, j) = 0.5 * (turb_flux_uv(i - 1, kts, j) + &
              turb_flux_uv(i, kts, j)) * 0.5 * (rho(i - 1, kts, j) +   &
              rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_u2_at_u(i, k, j) = 0.5 * (turb_flux_u2(i - 1, k, j) + &
                turb_flux_u2(i, k, j)) * 0.25 * (rho(i - 1, k, j) +   &
                rho(i, k, j) + rho(i - 1, k - 1, j) + rho(i, k - 1, j))
            tf_uv_at_u(i, k, j) = 0.5 * (turb_flux_uv(i - 1, k, j) + &
                turb_flux_uv(i, k, j)) * 0.25 * (rho(i - 1, k, j) +   &
                rho(i, k, j) + rho(i - 1, k - 1, j) + rho(i, k - 1, j))
              ! Horizontal derivatives of z at u points
            zy_at_u(i, k, j) = 0.125 * (zy(i - 1, k, j) + zy(i, k, j) + &
                zy(i - 1, k, j + 1) + zy(i, k, j + 1) + zy(i - 1, k + 1, j) + &
                zy(i, k + 1, j) + zy(i - 1, k + 1, j + 1) + zy(i, k + 1, j + 1))
            zx_at_u(i, k, j) = 0.5 * (zx(i, k, j) + zx(i, k + 1 , j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_u2_at_u(i, ktf + 1, j) = 0.0
          tf_uv_at_u(i, ktf + 1, j) = 0.0
          zx_at_u(i, kts, j) = 0.5 * (zx(i, kts, j) + zx(i, kts + 1 , j))
          zy_at_u(i, kts, j) = 0.125 * (zy(i - 1, kts, j) + zy(i, kts, j) + &
              zy(i - 1, kts, j + 1) + zy(i, kts, j + 1) + zy(i - 1, kts + 1, j) + &
              zy(i, kts + 1, j) + zy(i - 1, kts + 1, j + 1) + zy(i, kts + 1, j + 1))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msfux(i, j) * rdx
            mrdy = msfuy(i, j) * rdy
            tmpdz = (1.0 / rdzw(i, k, j) + 1.0 / rdzw(i - 1, k, j)) / 2.0
            tmpalpha = 0.5 * (sf_alpha(i - 1, k, j) + sf_alpha(i, k, j))
            tendency(i, k, j) = tendency(i, k, j) + ( G * tmpdz / dnw(k) *             &
               (mrdx * (tf_u2_at_mass(i, k, j) - tf_u2_at_mass(i - 1, k, j)) +                      &
                mrdy * (tf_uv_at_corners(i, k, j + 1) - tf_uv_at_corners(i, k, j)) -                      &
                msfux(i, j) * zx_at_u(i, k, j) * (tf_u2_at_u(i, k + 1, j) - tf_u2_at_u(i, k, j)) / tmpdz - &
                msfuy(i, j) * zy_at_u(i, k, j) * (tf_uv_at_u(i, k + 1, j) - tf_uv_at_u(i, k, j)) / tmpdz) ) / (N_PBL3D * tmpalpha)
          end do
        end do
      end do

    end subroutine Calc_htend_du


    subroutine Calc_htend_dv (config_flags, tendency, turb_flux_v2, turb_flux_uv, tf_v2_at_mass, tf_uv_at_corners, &
              msfvx, msfvy, rdx, rdy, dnw, zx, zy, rdzw, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dnw
      real, dimension (ims:ime, jms:jme), intent(in) ::  msfvx, msfvy
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdzw, rho, &
          zx, zy, turb_flux_v2, turb_flux_uv
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1), intent(in) :: &
          tf_v2_at_mass, tf_uv_at_corners

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpdz, tmpalpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: &
          tf_v2_at_v, tf_uv_at_v, zy_at_v, zx_at_v
      real (kind = kind(tf_v2_at_v)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_uv_at_v = rho * turb_flux_uv at v points in the faces
        ! tf_v2_at_v = rho * turb_flux_v2 at v points in the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_uv_at_v(i, kts, j) = 0.5 * (turb_flux_uv(i, kts, j - 1) + &
              turb_flux_uv(i, kts, j)) * 0.5 * (rho(i, kts, j - 1) +   &
              rho(i, kts, j))
          tf_v2_at_v(i, kts, j) = 0.5 * (turb_flux_v2(i, kts, j - 1) + &
              turb_flux_v2(i, kts, j)) * 0.5 * (rho(i, kts, j - 1) +   &
              rho(i, kts, j))
        end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_uv_at_v(i, k, j) = 0.5 * (turb_flux_uv(i, k, j - 1) + &
                turb_flux_uv(i, k, j)) * 0.25 * (rho(i, k, j - 1) +   &
                rho(i, k, j) + rho(i, k - 1, j - 1) + rho(i, k - 1, j))
            tf_v2_at_v(i, k, j) = 0.5 * (turb_flux_v2(i, k, j - 1) + &
                turb_flux_v2(i, k, j)) * 0.25 * (rho(i, k, j - 1) +   &
                rho(i, k, j) + rho(i, k - 1, j - 1) + rho(i, k - 1, j))
              ! Horiontal derivatives of z at v points
            zx_at_v(i, k, j) = 0.125 * (zx(i, k, j) + zx(i + 1, k, j) +     &
                zx(i, k, j - 1) + zx(i + 1, k, j - 1) + zx(i, k + 1, j) +   &
                zx(i + 1, k + 1, j) + zx(i, k + 1, j - 1) + zx(i + 1, k + 1, j - 1))
            zy_at_v(i, k, j) = 0.5 * (zy(i, k, j) + zy(i, k + 1 , j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_uv_at_v(i, ktf + 1,j)= 0.0
          tf_v2_at_v(i, ktf + 1, j) = 0.0
          zx_at_v(i, kts, j) = 0.125 * (zx(i, kts, j) + zx(i + 1, kts, j) +     &
              zx(i, kts, j - 1) + zx(i + 1, kts, j - 1) + zx(i, kts + 1, j) + &
              zx(i + 1, kts + 1, j) + zx(i, kts + 1, j - 1) + zx(i + 1, kts + 1, j - 1))
          zy_at_v(i, kts, j) = 0.5 * (zy(i, kts, j) + zy(i, kts + 1 , j))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msfvx(i, j) * rdx
            mrdy = msfvy(i, j) * rdy
            tmpdz = (1.0 / rdzw(i, k, j) + 1.0 / rdzw(i, k, j - 1)) / 2.0
            tmpalpha = 0.5 * (sf_alpha(i, k, j - 1) + sf_alpha(i, k, j))
            tendency(i, k, j) = tendency(i, k, j) + ( G * tmpdz / dnw(k) * &
                (mrdy * (tf_v2_at_mass(i, k, j) - tf_v2_at_mass(i, k, j - 1)) +       &
                 mrdx * (tf_uv_at_corners(i + 1, k, j) - tf_uv_at_corners(i, k, j)) -       &
                 msfvx(i, j) * zx_at_v(i, k, j) * (tf_uv_at_v(i, k + 1, j) - tf_uv_at_v(i, k, j)) / tmpdz - &
                 msfvy(i, j) * zy_at_v(i, k, j) * (tf_v2_at_v(i, k + 1, j) - tf_v2_at_v(i, k, j)) / tmpdz) ) / (N_PBL3D * tmpalpha)
          end do
        end do
      end do

    end subroutine Calc_htend_dv


    subroutine Calc_htend_dw (config_flags, tendency, turb_flux_uw, turb_flux_vw, tf_uw_at_full, tf_vw_at_full, &
              msftx, msfty, rdx, rdy, dn, zx, zy, rdz, rho, sf_alpha, ims, ime, jms, jme, kms, kme, &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf
      real,    intent(in) :: rdx, rdy

      real, dimension (kms:kme), intent (in) :: dn
      real, dimension (ims:ime, jms:jme), intent(in) ::  msftx, msfty
      real, dimension (its:ite, kts:kte, jts:jte), intent(inout) :: tendency
      real, dimension (ims:ime, kms:kme, jms:jme), intent(in) :: rdz, rho, &
          zx, zy, turb_flux_uw, turb_flux_vw
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2), intent(in) :: &
          tf_uw_at_full, tf_vw_at_full

        ! Local data
      integer :: i, j, k
      real    :: mrdx, mrdy, tmpalpha
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: &
          tf_uw_at_mass, tf_vw_at_mass, zy_at_w, zx_at_w
      real (kind = kind(tf_uw_at_mass)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_uw_at_mass = turb_flux_uw * rho at mass points
        ! tf_vw_at_mass = turb_flux_vw * rho at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tf_uw_at_mass(i, k, j) = 0.5 * (turb_flux_uw(i, k, j) + turb_flux_uw(i, k + 1, j)) *  rho(i, k, j)
            tf_vw_at_mass(i, k, j) = 0.5 * (turb_flux_vw(i, k, j) + turb_flux_vw(i, k + 1, j)) *  rho(i, k, j)
              ! Horiontal derivatives of z at center of faces
            zx_at_w(i, k, j) = 0.5 * (zx(i, k, j) + zx(i + 1, k, j))
            zy_at_w(i, k, j) = 0.5 * (zy(i, k, j) + zy(i, k, j + 1))
          end do
        end do
      end do

        ! Lid
      do j = j_start, j_end
        do i = i_start, i_end
          tf_uw_at_mass(i, ktf + 1, j) = 0.0
          tf_vw_at_mass(i, ktf + 1, j) = 0.0
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy
            tmpalpha = 0.5 * (sf_alpha(i, k - 1, j) + sf_alpha(i, k, j))
            if (USE_HIGH_ORDER) then
              tendency(i, k, j) = tendency(i, k, j) + ( G / (dn(k) * rdz(i, k, j)) * (   &
                  mrdx * ((8.0 * (tf_uw_at_full(i + 1, k, j) - tf_uw_at_full(i - 1, k, j))) -  &
                  (tf_uw_at_full(i + 2, k, j) - tf_uw_at_full(i - 2, k, j))) / 12.0 +          &
                  mrdy * (8.0 * (tf_vw_at_full(i, k, j + 1) - tf_vw_at_full(i, k, j - 1)) -    &
                  (tf_vw_at_full(i, k, j + 2) - tf_vw_at_full(i, k, j - 2))) / 12.0 -          &
                  rdz(i, k, j) * (msftx(i, j) * zx_at_w(i, k, j) * (tf_uw_at_mass(i, k, j) - tf_uw_at_mass(i, k - 1, j)) +  &
                                  msfty(i, j) * zy_at_w(i, k, j) * (tf_vw_at_mass(i, k, j) - tf_vw_at_mass(i, k - 1, j))) ) ) / (N_PBL3D * tmpalpha)
            else
              tendency(i, k, j) = tendency(i, k, j) + ( G / (dn(k) * rdz(i, k, j)) * (   &
                  mrdx * (tf_uw_at_full(i + 1, k, j) - tf_uw_at_full(i, k, j)) +       &
                  mrdy * (tf_vw_at_full(i, k, j + 1) - tf_vw_at_full(i, k, j)) -       &
                  rdz(i, k, j) * (msftx(i, j) * zx_at_w(i, k, j) * (tf_uw_at_mass(i, k, j) - tf_uw_at_mass(i, k - 1, j)) +  &
                                  msfty(i, j) * zy_at_w(i, k, j) * (tf_vw_at_mass(i, k, j) - tf_vw_at_mass(i, k - 1, j))) ) ) / (N_PBL3D * tmpalpha)
            end if
          end do
        end do
      end do

    end subroutine Calc_htend_dw


    subroutine Calc_htend_ds (config_flags, tendency, turb_flux_x, turb_flux_y, tf_x_at_u, tf_y_at_v, rho, sf_alpha, &
              fnm, fnp, zx, zy, dnw, rdzw, msfty, msftx, rdx, rdy, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, &
              kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      integer, intent(in) :: ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, intent (in) :: rdx, rdy
      real, dimension (kms:kme), intent (in) :: fnm, fnp, dnw
      real, dimension (ims:ime, jms:jme), intent (in) :: msftx, msfty
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in) :: turb_flux_x, &
          turb_flux_y, rho, rdzw, zx, zy
      real, dimension (its-1:ite, kts:kte, jts-1:jte), intent(in) :: sf_alpha
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ),  &
          intent (in) :: tf_x_at_u, tf_y_at_v
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: tendency

      integer, intent(in) :: i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: mrdx, mrdy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1 ) ::  &
          tf_x_at_full, tf_y_at_full, zx_at_m, zy_at_m
      real (kind = kind(tf_x_at_full)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

        ! Calculate component related to slope in sigma surfaces
        ! tf_x_at_full = turb_flux_x * rho at the center of the faces
        ! tf_y_at_full = turb_flux_y * rho at the center of the faces
      do j = j_start, j_end
        do i = i_start, i_end
          tf_x_at_full(i, kts, j) = turb_flux_x(i, kts, j) * rho(i, kts, j)
          tf_y_at_full(i, kts, j) = turb_flux_y(i, kts, j) * rho(i, kts, j)
          end do
      end do

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            tf_x_at_full(i, k, j) = turb_flux_x(i, k, j) * (fnm(k) * rho(i, k, j) + &
                fnp(k) * rho(i, k - 1, j))
            tf_y_at_full(i, k, j) = turb_flux_y(i, k, j) * (fnm(k) * rho(i, k, j) + &
                fnp(k) * rho(i, k - 1, j))
            zx_at_m(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + &
                zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_m(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + &
                zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          tf_x_at_full(i, ktf + 1, j) = 0.0
          tf_y_at_full(i, ktf + 1, j) = 0.0
          zx_at_m(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts, j) + &
              zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_m(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + &
              zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

        ! Calc tendency using coordinate transform (sigma to cartesian)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy
            tendency(i, k, j) = tendency(i, k, j) + ( G / (dnw(k) * rdzw(i, k, j)) * (  &
                mrdx * (tf_x_at_u(i + 1, k, j) - tf_x_at_u(i, k, j)) +                &
                mrdy * (tf_y_at_v(i, k, j + 1) - tf_y_at_v(i, k, j)) -                &
                msftx(i, j) * zx_at_m(i, k, j) * (tf_x_at_full(i, k + 1, j) - tf_x_at_full(i, k, j)) * rdzw(i, k, j) - &
                msfty(i, j) * zy_at_m(i, k, j) * (tf_y_at_full(i, k + 1, j) - tf_y_at_full(i, k, j)) * rdzw(i, k, j)) ) / (N_PBL3D * sf_alpha(i, k, j))
          end do
        end do
      end do

    end subroutine Calc_htend_ds


    subroutine Calc_q_sq_l2p5 (config_flags, u, v, w, zx, zy, t, qv, rho, thetav, tsk, t2, q2, psfc, ust, hfx, &
              du_dx_faces, du_dy_faces, du_dz_faces, dv_dx_faces, dv_dy_faces, dv_dz_faces, dw_dx_faces, dw_dy_faces, dw_dz_faces, &
              dthetav_dz_faces, q_sq, q_sq_prog, q_sq_tend, q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, &
              q_sq_dissip_tend, l_master_at_face, l_master_at_mass, l_boulac, turb_flux_u2_face, turb_flux_v2_face, turb_flux_w2_face, &
              turb_flux_uv_face, turb_flux_uw_face, turb_flux_vw_face, turb_flux_wtheta_v_face, turb_flux_u2, turb_flux_v2, &
              turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, dn, dnw, fnm, fnp, &
              msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates derivatives and fluxes at mass points                              !
      !          to prepare for prognostic TKE equation (level 2.5)                            !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      ! Method: 1) Follows module_diffusion_em to calc shear, buoyancy, and dissipation terms  !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: u, v, w, zx, zy, t, qv, rho, thetav
      real, dimension (ims:ime, jms:jme),          intent (in)    :: tsk, t2, q2, psfc, ust, hfx
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dx_faces, du_dy_faces, du_dz_faces, &
                                                                     dv_dx_faces, dv_dy_faces, dv_dz_faces, &
                                                                     dw_dx_faces, dw_dy_faces, dw_dz_faces, &
                                                                     dthetav_dz_faces
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend,    &
                                                                     q_sq_hdiff_tend,    &
                                                                     q_sq_shear_tend,    &
                                                                     q_sq_buoyancy_tend, &
                                                                     q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_face, l_boulac
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_wtheta_v
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: turb_flux_u2_face, turb_flux_v2_face, turb_flux_w2_face, &
                                                                     turb_flux_uv_face, turb_flux_uw_face, turb_flux_vw_face, &
                                                                     turb_flux_wtheta_v_face
      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1, c2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: z_at_mass, rdz, rdzw
      real,                                        intent (in)    :: rdx, rdy, cf1, cf2, cf3, dt
      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dz_mass, dv_dz_mass, dw_dz_mass, du_dx_mass, dv_dx_mass, dw_dx_mass, &
           du_dy_mass, dv_dy_mass, dw_dy_mass, dthetav_dz_mass, dqv_dz_mass
      real, dimension (its:ite, jts:jte)          :: thetav_sk

      ktf     = MIN( kte, kde-1 )
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      call Calc_wind_derivatives_at_mass (config_flags, u, v, w, zx, zy, du_dz_mass, dv_dz_mass, dw_dz_mass, &
          du_dx_mass, dv_dx_mass, dw_dx_mass, du_dy_mass, dv_dy_mass, dw_dy_mass, rdzw, dn, dnw, fnm, fnp, msftx, msfty, &
          msfux, msfuy, msfvx, msfvy, rdx, rdy, cf1, cf2, cf3, ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Calc_surface_thetav (tsk, t2, q2, psfc, thetav_sk, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Calc_mass_var_derivatives_at_mass (thetav, thetav_sk, q2, qv, dthetav_dz_mass, dqv_dz_mass, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

!      call Fill_derivatives_mass_with_derivatives_face (du_dx_faces, du_dy_faces, du_dz_faces, dv_dx_faces, dv_dy_faces, dv_dz_faces, &
!          dw_dx_faces, dw_dy_faces, dw_dz_faces, dthetav_dz_faces, du_dx_mass, du_dy_mass, du_dz_mass, dv_dx_mass, dv_dy_mass, dv_dz_mass, &
!          dw_dx_mass, dw_dy_mass, dw_dz_mass, dthetav_dz_mass, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Fill_fluxes_mass_with_fluxes_face (turb_flux_u2_face, turb_flux_v2_face, turb_flux_w2_face, turb_flux_uv_face, turb_flux_uw_face, &
          turb_flux_vw_face, turb_flux_wtheta_v_face, turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, &
          turb_flux_wtheta_v, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      call Fill_l_mass_with_l_face (l_master_at_face, l_master_at_mass, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        ! Calc TKE tendency
      call Calc_q_sq_rhs (config_flags, zx, zy, du_dz_mass, dv_dz_mass, dw_dz_mass, du_dx_mass, dv_dx_mass, dw_dx_mass, &
          du_dy_mass, dv_dy_mass, dw_dy_mass, dthetav_dz_mass, q_sq, q_sq_prog, q_sq_tend, &
          q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, q_sq_buoyancy_tend, q_sq_dissip_tend, &
          l_master_at_face, l_master_at_mass, l_boulac, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
          turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, dn, dnw, fnm, fnp, &
          msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
          ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

    end subroutine Calc_q_sq_l2p5


    subroutine Calc_wind_derivatives_at_mass (config_flags, u, v, w, zx, zy, du_dz, dv_dz, dw_dz, &
              du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, rdzw, dn, dnw, fnm, fnp, msftx, msfty, &
              msfux, msfuy, msfvx, msfvy, rdx, rdy, cf1, cf2, cf3, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                                         !
        ! Purpose: Calculate spatial derivatives of the wind variables at mass points             !
        !                                                                                         !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez                                          !
        !                                                                                         !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, v, w, zx, zy, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz, dv_dz, dw_dz, &
                                                                   du_dx, dv_dx, dw_dx, &
                                                                   du_dy, dv_dy, dw_dy
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3

      logical, parameter :: DEBUG = .false.

      call Calc_du_dz_at_mass (u, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz_at_mass (v, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dw_dz_at_mass (w, rdzw, dn, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      call Calc_du_dx_at_mass (u, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx)

      call Calc_dv_dx_at_mass (v, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, &
          ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dx)

      call Calc_dw_dx_at_mass (w, zx, rdzw, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
          ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dx)

      call Calc_du_dy_at_mass (u, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy)

      call Calc_dv_dy_at_mass (v, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds, &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy)

      call Calc_dw_dy_at_mass (w, zy, rdzw, msfty, rdy, ims, ime, jms, jme, kms, kme, &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy)

      if (DEBUG) print *, 'Done Calc_wind_derivatives_at_mass...'

    end subroutine Calc_wind_derivatives_at_mass


    subroutine Calc_du_dz_at_mass (u, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dz
      real,                                        intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: u_at_full
      real, dimension (its:ite, kts:kte, jts:jte)   :: u_at_face
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_at_full(i, k, j) = fnm(k) * u(i, k, j) + fnp(k) * u(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_at_full(i, 1, j) =  0.0
            ! Top face
          u_at_full(i, kte, j) =  cft1 * u(i, ktes1, j) + cft2 * u(i, ktes2, j)
        end do
      end do

        ! Calculate u at the center of the faces
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            u_at_face(i, k, j) = 0.5 * (u_at_full(i, k, j) + u_at_full(i + 1, k, j))
          end do
        end do
      end do

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dz(i, k, j) = (u_at_face(i, k + 1, j) - u_at_face(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          du_dz(i, kts, j) = 0.0
            ! Top
          du_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dz_at_mass


    subroutine Calc_dv_dz_at_mass (v, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dz
      real,                                        intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: v_at_full
      real, dimension (its:ite, kts:kte, jts:jte)   :: v_at_face
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_at_full(i, k, j) = fnm(k) * v(i, k, j) + fnp(k) * v(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_at_full(i, 1, j) =  0.0
            ! Top face
          v_at_full(i, kte, j) =  cft1 * v(i, ktes1, j) + cft2 * v(i, ktes2, j)
        end do
      end do

        ! Calculate v at the center of the faces
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start, i_end
            v_at_face(i, k, j) = 0.5 * (v_at_full(i, k, j) + v_at_full(i, k, j + 1))
          end do
        end do
      end do

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dz(i, k, j) = (v_at_face(i, k + 1, j) - v_at_face(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dv_dz(i, kts, j) = 0.0
            ! Top
          dv_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dz_at_mass


    subroutine Calc_dw_dz_at_mass (w, rdzw, dn, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dw_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dz                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dz

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Calculate derivative at the mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dw_dz(i, k, j) = (w(i, k + 1, j) - w(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dw_dz(i, kts, j) = 0.0
            ! Top
          dw_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dz_at_mass


    subroutine Calc_du_dx_at_mass (u, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfux, rdx, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zx, rdzw
      real, dimension (kms:kme), intent (in)                    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfux
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dx

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2, dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_ux
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: du_dz_dz_dx, u_mux, u_mux_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at u points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_ux(i, j) = msftx(i, j) * msfux(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            u_mux(i, k, j) = u(i, k, j) / msfux(i, j)
          end do
        end do
      end do

        ! Interpolate u_mux from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            u_mux_at_full(i, k, j) = fnm(k) * u_mux(i, k, j) + fnp(k) * u_mux(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          u_mux_at_full(i, 1, j) =  cf1 * u_mux(i, 1, j) + cf2 * u_mux(i, 2, j) + cf3 * u_mux(i, 3, j)
            ! Top face
          u_mux_at_full(i, kte, j) =  cft1 *  u_mux(i, ktes1, j) +  cft2 * u_mux(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! du_dz_dz_dx = (partial u_mux / partial z)|sigma (du_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            du_dz_dz_dx(i, k, j) = (0.5 * (u_mux_at_full(i, k + 1, j) + u_mux_at_full(i + 1, k + 1, j)) - &
                0.5 * (u_mux_at_full(i, k, j) + u_mux_at_full(i + 1, k, j))) * rdzw(i, k, j) * dz_dx

            if (DEBUG) then
              if ((i == 49 .or. i == 50 .or. i == 51) .and. (k == 3 .or. k == 4) .and. j == 54) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'dz_dx: '
                print *, dz_dx
                print *, 'du_dz_dz_dx: '
                print *, du_dz_dz_dx(i, k, j)
              end if
            end if
          end do
        end do
      end do

        ! (partial u / partial x)|z = msf_tx_ux * (partial u_mux_at_full / partial x)|sigma - du_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dx(i, k, j) = msf_tx_ux(i, j) * ((u_mux(i + 1, k, j) - u_mux(i, k, j)) * rdx - du_dz_dz_dx(i, k, j))

            if (DEBUG) then
              if ((i == 49 .or. i == 50 .or. i == 51) .and. (k == 3 .or. k == 4) .and. j == 54) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'u_mux_at_full: '
                print *, u_mux_at_full(i + 1, k, j), u_mux_at_full(i, k, j)
                print *, 'rdx: '
                print *, rdx
                print *, 'du_dx: '
                print *, du_dx(i, k, j)
              end if
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          du_dx(i, kts, j) = 0.0
            ! Top
          du_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dx_at_mass


    subroutine Calc_dv_dx_at_mass (v, zx, rdzw, dn, dnw, fnm, fnp, msftx, msfvx, rdx, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zx, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfvx
      real,                                        intent (in)  :: rdx, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dx

        ! Local vars
      integer :: i, j, k
      real    :: dz_dx, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_tx_vx
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dx, v_mvx, v_mvx_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_tx_vx(i, j) = msftx(i, j) * msfvx(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection (at u points)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            v_mvx(i, k, j) = 0.25 * (v(i, k, j) + v(i, k, j + 1) + v(i - 1, k, j) + v(i - 1, k, j + 1)) / msfvx(i, j)
          end do
        end do
      end do

        ! Interpolate v_mvx from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            v_mvx_at_full(i, k, j) = fnm(k) * v_mvx(i, k, j) + fnp(k) * v_mvx(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
            ! Surface
          v_mvx_at_full(i, 1, j) = cf1 * v_mvx(i, 1, j) + cf2 * v_mvx(i, 2, j) + cf3 * v_mvx(i, 3, j)
            ! Top face
          v_mvx_at_full(i, kte, j) = cft1 * v_mvx(i, ktes1, j) + cft2 * v_mvx(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! dv_dz_dz_dx = (partial v_mvx / partial z)|sigma (dv_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dv_dz_dz_dx(i, k, j) = (0.5 * (v_mvx_at_full(i, k + 1, j) + v_mvx_at_full(i + 1, k + 1, j)) - &
                0.5 * (v_mvx_at_full(i, k, j) + v_mvx_at_full(i + 1, k, j))) * rdzw(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial v / partial x)|z = msf_tx_vx * (partial v_mvx_at_full / partial x)|sigma - dv_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dx(i, k, j) = msf_tx_vx(i, j) * ((v_mvx(i + 1, k, j) - v_mvx(i, k, j)) * rdx - dv_dz_dz_dx(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dv_dx(i, 1, j) =  0.0
          dv_dx(i, kte, j) =  0.0
        end do
      end do

    end subroutine Calc_dv_dx_at_mass


    subroutine Calc_dw_dx_at_mass (w, zx, rdzw, msftx, rdx, ims, ime, jms, jme, kms, kme,  &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dx)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dx                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zx, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx
      real,                                        intent (in)  :: rdx
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dx

        ! Local vars
      integer :: i, j, k, i_offset
      real    :: dz_dx
      real, dimension (its:ite, jts:jte) :: msf_tx_tx
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dx, w_mtx, w_mtx_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_tx_tx(i, j) = msftx(i, j) * msftx(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        i_offset = 2
      else
        i_offset = 1
      end if
      do j = j_start, j_end
        do k = kts, kte
          do i = i_start - i_offset, i_end + i_offset
            w_mtx(i, k, j) = w(i, k, j) / msftx(i, j)
          end do
        end do
      end do

        ! Average w_mtx to the mass points (half levels)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start - i_offset, i_end + i_offset
            w_mtx_at_half(i, k, j) = 0.5 * (w_mtx(i, k, j) + w_mtx(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial x)|sigma at mass points
            dz_dx = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
              ! dw_dz_dz_dx = (partial w_mtx / partial z)|sigma (dw_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dw_dz_dz_dx(i, k, j) = (w_mtx(i, k + 1, j) - w_mtx(i, k, j)) * rdzw(i, k, j) * dz_dx
          end do
        end do
      end do

        ! (partial w / partial x)|z = msf_tx_tx * (partial w_mtx / partial x)|sigma - dw_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end

            if (DEBUG) then
              if (i == 3 .and. k == 2 .and. j == 3) then
                print *, '-----------------------------'
                print *, 'ims = ', ims, 'ime = ', ime
                print *, 'i_start = ', i_start, 'i_end = ', i_end
                print *, 'jms = ', jms, 'jme = ', jme
                print *, 'j_start = ', j_start, 'j_end = ', j_end
                print *, 'Row: '
                print *, w(:, k, j)
                print *, 'Column: '
                print *, w(i, k, :)
              end if
            end if

            if (USE_HIGH_ORDER) then
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * (8.0 * (w_mtx_at_half(i + 1, k, j) - w_mtx_at_half(i - 1, k, j)) - (w_mtx_at_half(i + 2, k, j) - w_mtx_at_half(i - 2, k, j)) ) / 12.0 - dw_dz_dz_dx(i, k, j))
            else
              dw_dx(i, k, j) = msf_tx_tx(i, j) * (rdx * 0.5 * (w_mtx_at_half(i + 1, k, j) - w_mtx_at_half(i - 1, k, j)) - dw_dz_dz_dx(i, k, j))
            end if

          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dw_dx(i, 1, j) = 0.0
          dw_dx(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dx_at_mass


    subroutine Calc_du_dy_at_mass (u, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfuy, rdy, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate du/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: u, zy, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfuy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: du_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_uy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: du_dz_dz_dy, u_muy, u_muy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

      logical, parameter :: DEBUG = .false.

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_uy(i, j) = msfty(i, j) * msfuy(i, j)
        end do
      end do

        ! Scale u component of vector to account for map projection (at v points)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            u_muy(i, k, j) = 0.25 * (u(i, k, j) + u(i + 1, k, j) + u(i, k, j - 1) + u(i + 1, k, j - 1)) / msfuy(i, j)
          end do
        end do
      end do

        ! Interpolate u_muy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            u_muy_at_full(i, k, j) = fnm(k) * u_muy(i, k, j) + fnp(k) * u_muy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          u_muy_at_full(i, 1, j) =  cf1 * u_muy(i, 1, j) + cf2 * u_muy(i, 2, j) + cf3 * u_muy(i, 3, j)
            ! Top face
          u_muy_at_full(i, kte, j) = cft1 * u_muy(i, ktes1, j) +  cft2 * u_muy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! du_dz_dz_dy = (partial u_muy / partial z)|sigma (du_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            du_dz_dz_dy(i, k, j) = (0.5 * (u_muy_at_full(i, k + 1, j) + u_muy_at_full(i, k + 1, j + 1)) - &
                0.5 * (u_muy_at_full(i, k, j) + u_muy_at_full(i, k, j + 1))) * rdzw(i, k, j) * dz_dy

            if (DEBUG) then
              if (i == 54 .and. (k == 3 .or. k == 4) .and. (j == 59 .or. j == 60 .or. j == 61 .or. j == 62)) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'dz_dy: '
                print *, dz_dy
                print *, 'du_dz_dz_dy: '
                print *, du_dz_dz_dy(i, k, j)
              end if
            end if

          end do
        end do
      end do

        ! (partial u / partial y)|z = msf_ty_uy * (partial u_muy_at_full / partial y)|sigma - du_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dy(i, k, j) = msf_ty_uy(i, j) * ((u_muy(i, k, j + 1) - u_muy(i, k, j)) * rdy - du_dz_dz_dy(i, k, j))

            if (DEBUG) then
              if (i == 54 .and. (k == 3 .or. k == 4) .and. (j == 59 .or. j == 60 .or. j == 61 .or. j == 62)) then
                print *, '-----------------------------'
                print *, 'i = ', i, 'j = ', j, 'k = ', k
                print *, 'u_muy_at_full: '
                print *, u_muy_at_full(i, k, j + 1), u_muy_at_full(i, k, j)
                print *, 'rdy: '
                print *, rdy
                print *, 'du_dy: '
                print *, du_dy(i, k, j)
              end if
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
!          du_dy(i, kts, j) = 0.0
            ! Top face
          du_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_du_dy_at_mass


    subroutine Calc_dv_dy_at_mass (v, zy, rdzw, dn, dnw, fnm, fnp, msfty, msfvy, rdy, cf1, cf2, cf3, ide, jds, &
              jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, config_flags, dv_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dv/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: v, zy, rdzw
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty, msfvy
      real,                                        intent (in)  :: rdy, cf1, cf2, cf3
      integer,                                     intent (in)  :: ide, jds, jde, kde,            &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte
      type (grid_config_rec_type),                 intent (in)  :: config_flags

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dv_dy

        ! Local vars
      integer :: i, j, k
      real    :: dz_dy, cft1, cft2
      real, dimension (its:ite, jts:jte) :: msf_ty_vy
      real, dimension (its - 1:ite + 1, kts:kte, jts - 1:jte + 1) :: dv_dz_dz_dy, v_mvy, v_mvy_at_full
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points multiplied by map factor at v points
      do j = j_start, j_end
        do i = i_start, i_end
          msf_ty_vy(i, j) = msfty(i, j) * msfvy(i, j)
        end do
      end do

        ! Scale v component of vector to account for map projection
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! Because msfvy at the poles will be undefined (1./0.),
              ! we will have trouble. But we are OK since v at the
              ! poles is 0., and that takes precedence in this case
            if ((config_flags%polar) .and. ((j == jds) .or. (j == jde))) then
              v_mvy(i, k, j) = 0.0
            else
              v_mvy(i, k, j) = v(i, k, j) / msfvy(i, j)
            end if
          end do
        end do
      end do

        ! Interpolate v_mvy from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            v_mvy_at_full(i, k, j) = fnm(k) * v_mvy(i, k, j) + fnp(k) * v_mvy(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
            ! Surface
          v_mvy_at_full(i, 1, j) =  cf1 * v_mvy(i, 1, j) + cf2 * v_mvy(i, 2, j) + cf3 * v_mvy(i, 3, j)
            ! Top face
          v_mvy_at_full(i, kte, j) = cft1 * v_mvy(i, ktes1, j) + cft2 * v_mvy(i, ktes2, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! dv_dz_dz_dy = (partial v_mvy / partial z)|sigma (dv_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dv_dz_dz_dy(i, k, j) = (0.5 * (v_mvy_at_full(i, k + 1, j) + v_mvy_at_full(i, k + 1, j + 1)) - &
                0.5 * (v_mvy_at_full(i, k, j) + v_mvy_at_full(i, k, j + 1))) * rdzw(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial v / partial y)|z = msf_ty_vy * (partial v_mvy_at_full / partial y)|sigma - dv_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dv_dy(i, k, j) = msf_ty_vy(i, j) * ((v_mvy(i, k, j + 1) - v_mvy(i, k, j)) * rdy - dv_dz_dz_dy(i, k, j))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
!          dv_dy(i, 1, j) = 0.0
            ! Top face
          dv_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dv_dy_at_mass


    subroutine Calc_dw_dy_at_mass (w, zy, rdzw, msfty, rdy, ims, ime, jms, jme, kms, kme, &
              ide, jde, kde, its, ite, jts, jte, kts, kte, dw_dy)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate dw/dy                             !
      !                                                      !
      ! Authors: Timothy W. Juliano & Pedro A. Jimenez       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: w, zy, rdzw
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msfty
      real,                                        intent (in)  :: rdy
      integer,                                     intent (in)  :: ide, jde, kde,                 &
                                                                   ims, ime, jms, jme, kms, kme,  &
                                                                   its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dw_dy

        ! Local vars
      integer :: i, j, k, j_offset
      real    :: dz_dy
      real, dimension (its:ite, jts:jte) :: msf_ty_ty
      real, dimension (its - 2:ite + 2, kts:kte, jts - 2:jte + 2) :: dw_dz_dz_dy, w_mty, w_mty_at_half
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Map factor at mass points squared
      do j = jts, jte
        do i = its, ite
          msf_ty_ty(i, j) = msfty(i, j) * msfty(i, j)
        end do
      end do

        ! Scale w component of vector to account for map projection
      if (USE_HIGH_ORDER) then
        j_offset = 2
      else
        j_offset = 1
      end if
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, kte
          do i = i_start, i_end
            w_mty(i, k, j) = w(i, k, j) / msfty(i, j)
          end do
        end do
      end do

        ! Average w_mty to the mass points (half levels)
      do j = j_start - j_offset, j_end + j_offset
        do k = kts, ktf
          do i = i_start, i_end
            w_mty_at_half(i,k,j) = 0.5 * (w_mty(i, k, j) + w_mty(i, k + 1, j))
          end do
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dy = (partial z / partial y)|sigma at mass points
            dz_dy = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
              ! dw_dz_dz_dy = (partial w_mty / partial z)|sigma (dw_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dw_dz_dz_dy(i, k, j) = (w_mty(i, k + 1, j) - w_mty(i, k, j)) * rdzw(i, k, j) * dz_dy
          end do
        end do
      end do

        ! (partial w / partial y)|z = msf_ty_ty * (partial w_mty / partial y)|sigma - dw_dz_dz_dy)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            if (USE_HIGH_ORDER) then
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * (8.0 * (w_mty_at_half(i, k, j + 1) - w_mty_at_half(i, k, j - 1)) - (w_mty_at_half(i, k, j + 2) - w_mty_at_half(i, k, j - 2)) ) / 12.0 - dw_dz_dz_dy(i, k, j))
            else
              dw_dy(i, k, j) = msf_ty_ty(i, j) * (rdy * 0.5 * (w_mty_at_half(i, k, j + 1) - w_mty_at_half(i, k, j - 1)) - dw_dz_dz_dy(i, k, j))
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!          dw_dy(i, kts, j) = 0.0
          dw_dy(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dw_dy_at_mass


    subroutine Calc_surface_thetav (tsk, t2, q2, psfc, thetav_sk, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, jms:jme),          intent (in)  :: tsk, t2, q2, psfc
      real, dimension (its:ite, jts:jte),          intent (out) :: thetav_sk
      integer,                                     intent (in)  :: ide, jde, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
        ! Local vars
      integer :: i, j
      integer :: i_start, i_end, j_start, j_end
      real    :: sqv, theta_sk

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

      do j = j_start, j_end
        do i = i_start, i_end
          sqv = q2(i, j) / (1.0 + q2(i, j))
          theta_sk = t2(i, j) * (100000. / psfc(i, j)) ** rcp
          thetav_sk(i, j) = theta_sk * (1.0 + P608 * sqv)
        end do
      end do

    end subroutine Calc_surface_thetav


    subroutine Calc_mass_var_derivatives_at_mass (thetav, thetav_sk, q2, qv, dthetav_dz, dqv_dz, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                     !
      ! Purpose: Calculate spatial derivatives of mass variables at mass points             !
      !                                                                                     !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano                                      !
      !                                                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: thetav, qv
      real, dimension (its:ite, jts:jte),          intent (in)  :: thetav_sk
      real, dimension (ims:ime, jms:jme),          intent (in)  :: q2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdzw
      integer,                                     intent (in)  :: ids, ide, jds, jde, kds, kde, &
                                                                   ims, ime, jms, jme, kms, kme, &
                                                                   its, ite, jts, jte, kts, kte
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: dthetav_dz, dqv_dz
      real,                                        intent (in)  :: cf1, cf2, cf3

      logical, parameter :: DEBUG = .false.

        ! Calculate thetav derivatives
      call Calc_dmassvar_dz_at_mass (var=thetav, thetav_sk=thetav_sk, rdzw=rdzw, dn=dn, dnw=dnw, fnm=fnm, fnp=fnp, &
          cf1=cf1, cf2=cf2, cf3=cf3, ide=ide, jde=jde, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, dvar_dz=dthetav_dz)

        ! Calculate thetav derivatives
      call Calc_dmassvar_dz_at_mass (var=qv, q2=q2, rdzw=rdzw, dn=dn, dnw=dnw, fnm=fnm, fnp=fnp, &
          cf1=cf1, cf2=cf2, cf3=cf3, ide=ide, jde=jde, kde=kde, ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
          its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte, dvar_dz=dqv_dz)

      if (DEBUG) print *, 'Done Calc_mass_var_derivatives_at_mass...'

    end subroutine Calc_mass_var_derivatives_at_mass


    subroutine Calc_dmassvar_dz_at_mass (var, thetav_sk, q2, rdzw, dn, dnw, fnm, fnp, cf1, cf2, cf3, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dvar_dz)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Calculate partial var over partial z        !
      !                                                      !
      ! Authors: Pedro A. Jimenez & Timothy W. Juliano       !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme),  intent (in)  :: var, rdzw
      real, dimension (its:ite, jts:jte), optional, intent (in)  :: thetav_sk
      real, dimension (ims:ime, jms:jme), optional, intent (in)  :: q2
      real, dimension (kms:kme),                    intent (in)  :: dn, dnw, fnm, fnp
      integer,                                      intent (in)  :: ide, jde, kde,                  &
                                                                    ims, ime, jms, jme, kms, kme,   &
                                                                    its, ite, jts, jte, kts, kte

      real, dimension (its:ite, kts:kte, jts:jte),  intent (out) :: dvar_dz
      real,                                         intent (in)  :: cf1, cf2, cf3

        ! Local vars
      integer :: i, j, k
      real    :: cft1, cft2
      integer :: ktf, ktes1, ktes2, i_start, i_end, j_start, j_end
      real, dimension (its:ite, kts:kte, jts:jte) :: var_at_full

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

        ! Interpolate from half levels to full levels
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            var_at_full(i, k, j) = fnm(k) * var(i, k, j) + fnp(k) * var(i, k - 1, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
            ! Surface
          if (PRESENT(thetav_sk)) then
            var_at_full(i, 1, j) =  thetav_sk(i, j) - 300.
          end if
          if (PRESENT(q2)) then
            var_at_full(i, 1, j) =  q2(i, j)
          end if
            ! Top face
          var_at_full(i, kte, j) =  cft1 * var(i, ktes1, j) + cft2 * var(i, ktes2, j)
        end do
      end do

        ! Calculate derivative at mass points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            dvar_dz(i, k, j) = (var_at_full(i, k + 1, j) - var_at_full(i, k, j)) * rdzw(i, k, j)
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
!            ! Surface
!          dvar_dz(i, kts, j) = 0.0
            ! Top
          dvar_dz(i, kte, j) = 0.0
        end do
      end do

    end subroutine Calc_dmassvar_dz_at_mass


    subroutine Fill_derivatives_mass_with_derivatives_face (du_dx_faces, du_dy_faces, du_dz_faces, dv_dx_faces, dv_dy_faces, dv_dz_faces, &
              dw_dx_faces, dw_dy_faces, dw_dz_faces, dthetav_dz_faces, du_dx_mass, du_dy_mass, du_dz_mass, dv_dx_mass, dv_dy_mass, dv_dz_mass, &
              dw_dx_mass, dw_dy_mass, dw_dz_mass, dthetav_dz_mass, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (its:ite, kts:kte, jts:jte), intent (in)   :: du_dx_faces, du_dy_faces, du_dz_faces, &
                                                                    dv_dx_faces, dv_dy_faces, dv_dz_faces, &
                                                                    dw_dx_faces, dw_dy_faces, dw_dz_faces, &
                                                                    dthetav_dz_faces
      real, dimension (its:ite, kts:kte, jts:jte), intent (out)  :: du_dx_mass, du_dy_mass, du_dz_mass, &
                                                                    dv_dx_mass, dv_dy_mass, dv_dz_mass, &
                                                                    dw_dx_mass, dw_dy_mass, dw_dz_mass, &
                                                                    dthetav_dz_mass
      integer,                                     intent (in)   :: ide, jde, kde,                  &
                                                                    ims, ime, jms, jme, kms, kme,   &
                                                                    its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

!      k = kts
!      do j = j_start, j_end
!          do i = i_start, i_end
!            du_dx_mass(i, k, j) = du_dx_faces(i, k + 1, j)
!            du_dy_mass(i, k, j) = du_dy_faces(i, k + 1, j)
!            du_dz_mass(i, k, j) = du_dz_faces(i, k + 1, j)
!            dv_dx_mass(i, k, j) = dv_dx_faces(i, k + 1, j)
!            dv_dy_mass(i, k, j) = dv_dy_faces(i, k + 1, j)
!            dv_dz_mass(i, k, j) = dv_dz_faces(i, k + 1, j)
!            dw_dx_mass(i, k, j) = dw_dx_faces(i, k + 1, j)
!            dw_dy_mass(i, k, j) = dw_dy_faces(i, k + 1, j)
!            dw_dz_mass(i, k, j) = dw_dz_faces(i, k + 1, j)
!            dthetav_dz_mass(i, k, j) = dthetav_dz_faces(i, k + 1, j)
!          end do
!      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            du_dx_mass(i, k, j) = 0.5 * ( du_dx_faces(i, k + 1, j) + du_dx_faces(i, k, j) )
            du_dy_mass(i, k, j) = 0.5 * ( du_dy_faces(i, k + 1, j) + du_dy_faces(i, k, j) )
            du_dz_mass(i, k, j) = 0.5 * ( du_dz_faces(i, k + 1, j) + du_dz_faces(i, k, j) )
            dv_dx_mass(i, k, j) = 0.5 * ( dv_dx_faces(i, k + 1, j) + dv_dx_faces(i, k, j) )
            dv_dy_mass(i, k, j) = 0.5 * ( dv_dy_faces(i, k + 1, j) + dv_dy_faces(i, k, j) )
            dv_dz_mass(i, k, j) = 0.5 * ( dv_dz_faces(i, k + 1, j) + dv_dz_faces(i, k, j) )
            dw_dx_mass(i, k, j) = 0.5 * ( dw_dx_faces(i, k + 1, j) + dw_dx_faces(i, k, j) )
            dw_dy_mass(i, k, j) = 0.5 * ( dw_dy_faces(i, k + 1, j) + dw_dy_faces(i, k, j) )
            dw_dz_mass(i, k, j) = 0.5 * ( dw_dz_faces(i, k + 1, j) + dw_dz_faces(i, k, j) )
            dthetav_dz_mass(i, k, j) = 0.5 * ( dthetav_dz_faces(i, k + 1, j) + dthetav_dz_faces(i, k, j) )
          end do
        end do
      end do

    end subroutine Fill_derivatives_mass_with_derivatives_face


    subroutine Fill_fluxes_mass_with_fluxes_face (tff_u2, tff_v2, tff_w2, tff_uv, tff_uw, tff_vw, tff_wthv, &
              tfm_u2, tfm_v2, tfm_w2, tfm_uv, tfm_uw, tfm_vw, tfm_wthv, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: tff_u2, tff_v2, tff_w2, tff_uv, tff_uw, tff_vw, tff_wthv
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)  :: tfm_u2, tfm_v2, tfm_w2, tfm_uv, tfm_uw, tfm_vw, tfm_wthv
      integer,                                     intent (in)   :: ide, jde, kde,                  &
                                                                    ims, ime, jms, jme, kms, kme,   &
                                                                    its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            tfm_u2(i, k, j) = 0.5 * ( tff_u2(i, k + 1, j) + tff_u2(i, k, j) )
            tfm_v2(i, k, j) = 0.5 * ( tff_v2(i, k + 1, j) + tff_v2(i, k, j) )
            tfm_w2(i, k, j) = 0.5 * ( tff_w2(i, k + 1, j) + tff_w2(i, k, j) )
            tfm_uv(i, k, j) = 0.5 * ( tff_uv(i, k + 1, j) + tff_uv(i, k, j) )
            tfm_uw(i, k, j) = 0.5 * ( tff_uw(i, k + 1, j) + tff_uw(i, k, j) )
            tfm_vw(i, k, j) = 0.5 * ( tff_vw(i, k + 1, j) + tff_vw(i, k, j) )
            tfm_wthv(i, k, j) = 0.5 * ( tff_wthv(i, k + 1, j) + tff_wthv(i, k, j) )
          end do
        end do
      end do

    end subroutine Fill_fluxes_mass_with_fluxes_face


    subroutine Fill_l_mass_with_l_face (l_master_at_face, l_master_at_mass, &
              ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      implicit none

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)   :: l_master_at_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)  :: l_master_at_mass
      integer,                                     intent (in)   :: ide, jde, kde,                  &
                                                                    ims, ime, jms, jme, kms, kme,   &
                                                                    its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k
      integer :: ktf, i_start, i_end, j_start, j_end

        ! Define some indexes
      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)
      ktf = Min (kte, kde - 1)

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            l_master_at_mass(i, k, j) = 0.5 * ( l_master_at_face(i, k + 1, j) + l_master_at_face(i, k, j) )
          end do
        end do
      end do

    end subroutine Fill_l_mass_with_l_face


    subroutine Calc_q_sq_rhs (config_flags, zx, zy, du_dz, dv_dz, dw_dz, du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, &
              dthetav_dz, q_sq, q_sq_prog, q_sq_tend, q_sq_vdiff_tend, q_sq_hdiff_tend, q_sq_shear_tend, &
              q_sq_buoyancy_tend, q_sq_dissip_tend, l_master_at_face, l_master_at_mass, l_boulac, turb_flux_u2, turb_flux_v2, turb_flux_w2, &
              turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_wtheta_v, dn, dnw, fnm, fnp, &
              msftx, msfty, msfux, msfuy, msfvx, msfvy, mut, c1, c2, z_at_mass, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, dt, &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                                        !
      ! Purpose: Calculates RHS of prognostic TKE equation (level 2.5)                         !
      !                                                                                        !
      ! Author: Timothy W. Juliano & Pedro A. Jimenez                                          !
      !                                                                                        !
      ! Method: 1) Follows module_diffusion_em to calc shear, buoyancy, and dissipation terms  !
      !                                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (grid_config_rec_type), intent (in) :: config_flags

      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: zx, zy
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)    :: du_dz, dv_dz, dw_dz, &
                            du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, dthetav_dz
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: q_sq_vdiff_tend,    &
                                                                     q_sq_hdiff_tend,    &
                                                                     q_sq_shear_tend,    &
                                                                     q_sq_buoyancy_tend, &
                                                                     q_sq_dissip_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: q_sq_prog
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_face, l_boulac
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: l_master_at_mass
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out)   :: turb_flux_u2, turb_flux_v2, turb_flux_w2, &
                                                                     turb_flux_uv, turb_flux_uw, turb_flux_vw, &
                                                                     turb_flux_wtheta_v
      real, dimension (kms:kme),                   intent (in)    :: dn, dnw, fnm, fnp
      real, dimension (ims:ime, jms:jme),          intent (in)    :: msftx, msfty, msfux, msfuy, msfvx, msfvy
      real, dimension (ims:ime, jms:jme),          intent (in)    :: mut
      real, dimension (kms:kme),                   intent (in)    :: c1, c2
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)    :: z_at_mass, rdz, rdzw
      real,                                        intent (in)    :: rdx, rdy, cf1, cf2, cf3, dt
      integer, intent (in) :: ids, ide, jds, jde, kds, kde,  &
                              ims, ime, jms, jme, kms, kme,  &
                              its, ite, jts, jte, kts, kte

        ! Local vars
      integer :: i, j, k, i_start, i_end, j_start, j_end, ktf
      real, dimension (its:ite, kts:kte, jts:jte) :: l_dissip
      logical, parameter :: DEBUG = .false.
      real (kind = kind(q_sq_tend)) :: N_PBL3D

      N_PBL3D = config_flags%pbl3d_nsteps

      ktf     = MIN( kte, kde-1 )
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. config_flags%specified .or. &
           config_flags%nested) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. config_flags%specified .or. &
           config_flags%nested) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. config_flags%specified .or. &
           config_flags%nested) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. config_flags%specified .or. &
           config_flags%nested) j_end   = MIN(jde-2,jte)

        !!! Diffusion !!!
      call Calc_q_sq_vertical_diffusion (N_PBL3D, q_sq_vdiff_tend, q_sq_tend, l_master_at_face, q_sq_prog, q_sq, &
          fnm, fnp, dnw, rdz, rdzw, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      if ( config_flags%pbl3d_opt .gt. 1 ) then
        call Calc_q_sq_horizontal_diffusion (N_PBL3D, q_sq_hdiff_tend, q_sq_tend, l_master_at_face, q_sq_prog, q_sq, &
            msftx, msfty, msfux, msfvy, fnm, fnp, dn, dnw, zx, zy, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, &
            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      end if

        !!! Shear !!!
      if ( config_flags%pbl3d_opt .lt. 2 ) then
        call Calc_q_sq_shear_pbl_approx (N_PBL3D, l_master_at_mass, q_sq_prog, q_sq_shear_tend, q_sq_tend, du_dz, dv_dz, turb_flux_uw, turb_flux_vw, &
            mut, c1, c2, fnm, fnp, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      else
        call Calc_q_sq_shear (N_PBL3D, q_sq_shear_tend, q_sq_tend, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
            turb_flux_u2, turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv, mut, c1, c2, fnm, fnp, &
            ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)
      end if

        !!! Buoyancy !!!
      call Calc_q_sq_buoyancy (N_PBL3D, q_sq_buoyancy_tend, q_sq_tend, turb_flux_wtheta_v, dthetav_dz, q_sq_prog, l_master_at_mass, mut, c1, c2, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        !!! Dissipation !!!
      call Fill_dissip_length_scale (config_flags%pbl3d_l_opt, l_master_at_mass, l_boulac, l_dissip, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      call Calc_q_sq_dissip (N_PBL3D, q_sq_dissip_tend, q_sq_tend, l_master_at_mass, l_dissip, q_sq_prog, mut, c1, c2, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

        !!! Bound tendency !!!
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_tend(i, k, j) = max ( q_sq_tend(i, k, j), -(c1(k) * mut(i,j) + c2(k)) * max ( 0.0 , N_PBL3D * q_sq_prog(i, k, j) ) / dt )

            if (DEBUG) call Debug_q_sq_tend
          end do
        end do
      end do

      contains

      subroutine Debug_q_sq_tend

        if (i == I_TO_TEST .and. j == J_TO_TEST ) then
          print *, '------------------------------------------'
          print *, '  q_sq_tend final = ', q_sq_tend(i, k, j)
        end if

      end subroutine Debug_q_sq_tend

    end subroutine Calc_q_sq_rhs


    subroutine Calc_q_sq_vertical_diffusion (N_PBL3D, q_sq_vdiff_tend, q_sq_tend, l_master, &
              q_sq_mass, q_sq_face, fnm, fnp, dnw, rdz, rdzw, ims, ime, jms, jme, kms, kme,     &
              its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf             )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                !
        ! Purpose: Calculate vertical diffusion term for prognostic TKE  !
        !          Follows MY82 Eq. 24                                   !
        !                                                                !
        ! Authors: Timothy W. Juliano & Pedro A. Jimenez                 !
        !                                                                !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_vdiff_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_mass, q_sq_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master
      real, dimension (kms:kme),                   intent (in)  :: dnw, fnm, fnp
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: rdz, rdzw
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real, parameter :: Sq = 0.20
      real :: q_sq_vdiff_tend_tmp
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq_vdiff_hold
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq_vdiff_hold(i, k, j) = - ( Sq * l_master(i, k, j) * Sqrt(q_sq_face(i, k, j)) * &
                                         ( q_sq_mass(i, k, j) - q_sq_mass(i, k - 1, j) ) * rdz(i, k, j) )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  l_master = ', l_master(i, k, j)
              print *, '  q_sq_face = ', q_sq_face(i, k, j)
              print *, '  q_sq_mass = ', q_sq_mass(i, k, j), q_sq_mass(i, k - 1, j)
              print *, '  q_sq_vdiff_hold = ', q_sq_vdiff_hold(i, k, j)
            end if
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          q_sq_vdiff_hold(i, kts, j) = 0.0
          q_sq_vdiff_hold(i, ktf + 1, j) = 0.0
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_vdiff_tend(i, k, j) = - rdzw(i, k, j) * ( q_sq_vdiff_hold(i, k + 1, j) - q_sq_vdiff_hold(i, k, j) )
            q_sq_vdiff_tend_tmp = ( G / dnw(k) ) *  &
                                  ( q_sq_vdiff_hold(i, k + 1, j) - q_sq_vdiff_hold(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( q_sq_vdiff_tend_tmp / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '  q_sq_vdiff_tend = ', q_sq_vdiff_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if 
          end do
        end do
      end do

    end subroutine Calc_q_sq_vertical_diffusion


    subroutine Calc_q_sq_horizontal_diffusion (N_PBL3D, q_sq_hdiff_tend, q_sq_tend, l_master, q_sq_mass, q_sq_face, &
              msftx, msfty, msfux, msfvy, fnm, fnp, dn, dnw, zx, zy, rdx, rdy, rdz, rdzw, cf1, cf2, cf3, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_hdiff_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq_mass, q_sq_face
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master
      real, dimension (ims:ime, jms:jme),          intent (in)  :: msftx, msfty, msfux, msfvy
      real, dimension (kms:kme),                   intent (in)  :: dn, dnw, fnm, fnp
      real,                                        intent (in)  :: rdx, rdy, cf1, cf2, cf3
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: zx, zy, rdz, rdzw
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k, ktes1, ktes2
      real :: cft1, cft2, tmpzx, rdzu, rdzv, dz_dx, dz_dy, &
              dq_sq_dz_dz_dx, dq_sq_dz_dz_dy, mrdx, mrdy
      real :: q_sq_hdiff_tend_tmp
      real, parameter :: Sq = 0.20
      real, dimension (its:ite, kts:kte, jts:jte) :: q_sq_hdiff_hold
      real, dimension (its-1:ite+1, kts:kte, jts-1:jte+1) :: xkxavg, q_sq_at_full_u, q_sq_at_full_v, &
           dq_sq_dx_at_u, dq_sq_dy_at_v, dq_sq_dx_at_w, dq_sq_dy_at_w, zx_at_mass, zy_at_mass

        ! Define some indexes
      ktes1 = kte - 1
      ktes2 = kte - 2
      cft2 = - 0.5 * dnw(ktes1) / dn(ktes1)
      cft1 = 1.0 - cft2

        !!! Do dq_sq/dx first !!!

         ! Interpolate q_sq from half levels to full levels (u points at faces)
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end + 1
            q_sq_at_full_u(i, k, j) = 0.5 * (fnm(k) * (q_sq_mass(i - 1, k, j) + q_sq_mass(i, k, j)) + &
                fnp(k) * (q_sq_mass(i - 1, k - 1, j) + q_sq_mass(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          q_sq_at_full_u(i, kts, j) = 0.5 * (cf1 * q_sq_mass(i, 1, j) + cf2 * q_sq_mass(i, 2, j) + &
              cf3 * q_sq_mass(i, 3, j) + cf1 * q_sq_mass(i - 1, 1, j) +  &
              cf2 * q_sq_mass(i - 1, 2, j) + cf3 * q_sq_mass(i - 1, 3, j))
          q_sq_at_full_u(i, ktf + 1, j) = 0.5 * (q_sq_mass(i, ktes1, j) + (q_sq_mass(i, ktes1, j) - &
              q_sq_mass(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              q_sq_mass(i - 1, ktes1, j) + (q_sq_mass(i - 1, ktes1, j) - &
              q_sq_mass(i - 1, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate average eddy diffusivity at u points
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
            xkxavg(i, k, j) = Sq * 0.25 * ( l_master(i - 1, k, j) * Sqrt(q_sq_face(i - 1, k, j)) +  &
                                            l_master(i, k, j) * Sqrt(q_sq_face(i, k, j))         +  &
                                            l_master(i - 1, k + 1, j) * Sqrt(q_sq_face(i - 1, k + 1, j)) +  &
                                            l_master(i, k + 1, j) * Sqrt(q_sq_face(i, k + 1, j))         )
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end + 1
          xkxavg(i, ktf + 1, j) = xkxavg(i, ktf, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dq_sq_dz_dz_dx)
      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end + 1
              ! dz_dx = (partial z / partial x)|sigma at u points
            dz_dx = 0.5 * (zx(i, k, j) + zx(i, k + 1, j))
            rdzu = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i - 1, k, j))
              ! dq_sq_dz_dz_dx = (partial q_sq / partial z)|sigma (dq_sq_dz) *
              ! (partial z / partial x)|sigma (dz_dx)
            dq_sq_dz_dz_dx = (q_sq_at_full_u(i, k + 1, j) - q_sq_at_full_u(i, k, j)) * rdzu * dz_dx
              ! (partial q_sq / partial x)|z = msfux * (partial q_sq_at_mass / partial x)|sigma - dq_sq_dz_dz_dx)
              ! at u points
            dq_sq_dx_at_u(i, k, j) = - msfux(i, j) * xkxavg(i, k, j) *  &
                                       0.5 * ( q_sq_at_full_u(i, k + 1, j) + q_sq_at_full_u(i, k, j) ) * &
                                       (q_sq_mass(i, k, j) - q_sq_mass(i - 1, k, j)) * rdx - dq_sq_dz_dz_dx
          end do
        end do
      end do

        !!! Do dq_sq/dy second !!!

         ! Interpolate q_sq from half levels to full levels (v points at faces)
      do j = j_start, j_end + 1
        do k = kts + 1, ktf
          do i = i_start, i_end
            q_sq_at_full_v(i, k, j) = 0.5 * (fnm(k) * (q_sq_mass(i, k, j - 1) + q_sq_mass(i, k, j)) + &
                fnp(k) * (q_sq_mass(i, k - 1, j - 1) + q_sq_mass(i, k - 1, j)))
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          q_sq_at_full_v(i, kts, j) = 0.5 * (cf1 * q_sq_mass(i, 1, j) + cf2 * q_sq_mass(i, 2, j) + &
              cf3 * q_sq_mass(i, 3, j) + cf1 * q_sq_mass(i, 1, j - 1) +  &
              cf2 * q_sq_mass(i, 2, j - 1) + cf3 * q_sq_mass(i, 3, j - 1))
          q_sq_at_full_v(i, ktf + 1, j) = 0.5 * (q_sq_mass(i, ktes1, j) + (q_sq_mass(i, ktes1, j) - &
              q_sq_mass(i, ktes2, j)) * 0.5 * dnw(ktes1) / dn(ktes1) + &
              q_sq_mass(i, ktes1, j - 1) + (q_sq_mass(i, ktes1, j - 1) - &
              q_sq_mass(i, ktes2, j - 1)) * 0.5 * dnw(ktes1) / dn(ktes1))
        end do
      end do

        ! Calculate average eddy diffusivity at v points
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
            xkxavg(i, k, j) = Sq * 0.25 * ( l_master(i, k, j - 1) * Sqrt(q_sq_face(i, k, j - 1)) +  &
                                            l_master(i, k, j) * Sqrt(q_sq_face(i, k, j))         +  &
                                            l_master(i, k + 1, j - 1) * Sqrt(q_sq_face(i, k + 1, j - 1)) +  &
                                            l_master(i, k + 1, j) * Sqrt(q_sq_face(i, k + 1, j))         )
          end do
        end do
      end do

      do j = j_start, j_end + 1
        do i = i_start, i_end
          xkxavg(i, ktf + 1, j) = xkxavg(i, ktf, j)
        end do
      end do

        ! Calculate component associated with slope in sigma surfaces (dq_sq_dz_dz_dy)
      do j = j_start, j_end + 1
        do k = kts, ktf
          do i = i_start, i_end
              ! dz_dx = (partial z / partial y)|sigma at y points
            dz_dy = 0.5 * (zy(i, k, j) + zy(i, k + 1, j))
            rdzv = 2. / (1. / rdzw(i, k, j) + 1. / rdzw(i, k, j - 1))
              ! dq_sq_dz_dz_dy = (partial q_sq / partial z)|sigma (dq_sq_dz) *
              ! (partial z / partial y)|sigma (dz_dy)
            dq_sq_dz_dz_dy = (q_sq_at_full_v(i, k + 1, j) - q_sq_at_full_v(i, k, j)) * rdzv * dz_dy
              ! (partial q_sq / partial y)|z = msfvy * (partial q_sq_at_mass / partial y)|sigma - dq_sq_dz_dz_dy)
              ! at v points
            dq_sq_dy_at_v(i, k, j) = - msfvy(i, j) * xkxavg(i, k, j) * &
                                       0.5 * ( q_sq_at_full_v(i, k + 1, j) + q_sq_at_full_v(i, k, j) ) * &
                                       (q_sq_mass(i, k, j) - q_sq_mass(i, k, j - 1)) * rdy - dq_sq_dz_dz_dy
          end do
        end do
      end do

        !!! Do full horizontal diffusion !!!

        ! Interpolate dq_sq_dx and dq_sq_dy to w points
      do j = j_start, j_end
        do k = kts + 1, ktf
          do i = i_start, i_end
            dq_sq_dx_at_w(i, k, j) = 0.5 * (fnm(k) * (dq_sq_dx_at_u(i + 1, k, j) + dq_sq_dx_at_u(i, k, j)) +  &
                                            fnp(k) * (dq_sq_dx_at_u(i + 1, k - 1, j) + dq_sq_dx_at_u(i, k - 1, j)))
            dq_sq_dy_at_w(i, k, j) = 0.5 * (fnm(k) * (dq_sq_dy_at_v(i, k, j + 1) + dq_sq_dy_at_v(i, k, j)) +  &
                                            fnp(k) * (dq_sq_dy_at_v(i, k - 1, j + 1) + dq_sq_dy_at_v(i, k - 1, j)))

            zx_at_mass(i, k, j) = 0.25 * (zx(i, k, j) + zx(i + 1, k, j) + zx(i, k + 1, j) + zx(i + 1, k + 1, j))
            zy_at_mass(i, k, j) = 0.25 * (zy(i, k, j) + zy(i, k, j + 1) + zy(i, k + 1, j) + zy(i, k + 1, j + 1))
          end do
        end do
      end do

      do j = j_start, j_end
        do i = i_start, i_end
          dq_sq_dx_at_w(i, kts, j) = 0.
          dq_sq_dx_at_w(i, ktf + 1, j) = 0.
          dq_sq_dy_at_w(i, kts, j) = 0.
          dq_sq_dy_at_w(i, ktf + 1,j) = 0.

          zx_at_mass(i, kts, j) = 0.25 * (zx(i, kts, j) + zx(i + 1, kts ,j) + zx(i, kts + 1, j) + zx(i + 1, kts + 1, j))
          zy_at_mass(i, kts, j) = 0.25 * (zy(i, kts, j) + zy(i, kts, j + 1) + zy(i, kts + 1, j) + zy(i, kts + 1, j + 1))
        end do
      end do

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            mrdx = msftx(i, j) * rdx
            mrdy = msfty(i, j) * rdy

            q_sq_hdiff_tend(i, k, j) = 2. * ( mrdx * (dq_sq_dx_at_u(i + 1, k, j) - dq_sq_dx_at_u(i, k, j)) + &
                                              mrdy * (dq_sq_dy_at_v(i, k, j + 1) - dq_sq_dy_at_v(i, k, j)) - &
                                              msftx(i, j) * zx_at_mass(i, k, j) * (dq_sq_dx_at_w(i, k + 1, j) - dq_sq_dx_at_w(i, k, j)) * rdzw(i, k, j) - &
                                              msfty(i, j) * zy_at_mass(i, k, j) * (dq_sq_dy_at_w(i, k + 1, j) - dq_sq_dy_at_w(i, k, j)) * rdzw(i, k, j) )
            q_sq_hdiff_tend_tmp = G / (dnw(k) * rdzw(i, k, j)) * &
                                  2. * ( mrdx * (dq_sq_dx_at_u(i + 1, k, j) - dq_sq_dx_at_u(i, k, j)) + &
                                         mrdy * (dq_sq_dy_at_v(i, k, j + 1) - dq_sq_dy_at_v(i, k, j)) - &
                                         msftx(i, j) * zx_at_mass(i, k, j) * (dq_sq_dx_at_w(i, k + 1, j) - dq_sq_dx_at_w(i, k, j)) * rdzw(i, k, j) - &
                                         msfty(i, j) * zy_at_mass(i, k, j) * (dq_sq_dy_at_w(i, k + 1, j) - dq_sq_dy_at_w(i, k, j)) * rdzw(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( q_sq_hdiff_tend_tmp / N_PBL3D )
          end do
        end do
      end do

    end subroutine Calc_q_sq_horizontal_diffusion


    subroutine Calc_q_sq_shear_pbl_approx (N_PBL3D, l_master, q_sq, q_sq_shear_tend, q_sq_tend, du_dz, dv_dz, turb_flux_uw, turb_flux_vw, &
              mut, c1, c2, fnm, fnp, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_shear_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: du_dz, dv_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_uw, &
                                                                   turb_flux_vw

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2, fnm, fnp
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_shear_tend(i, k, j) = -2. * ( turb_flux_uw(i, k, j) * du_dz(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dv_dz(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_shear_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_shear_tend = ', q_sq_shear_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_shear_pbl_approx


    subroutine Calc_q_sq_shear (N_PBL3D, q_sq_shear_tend, q_sq_tend, du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz, &
              turb_flux_u2, turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv, mut, c1, c2, fnm, fnp, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_shear_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: du_dx, du_dy, du_dz,  &
                                                                   dv_dx, dv_dy, dv_dz,  &
                                                                   dw_dx, dw_dy, dw_dz
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: turb_flux_u2,         &
          turb_flux_uw, turb_flux_v2, turb_flux_vw, turb_flux_w2, turb_flux_uv

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2, fnm, fnp
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_shear_tend(i, k, j) = -2. * ( turb_flux_u2(i, k, j) * du_dx(i, k, j) +  &
                                               turb_flux_uv(i, k, j) * du_dy(i, k, j) +  &
                                               turb_flux_uw(i, k, j) * du_dz(i, k, j) +  &
                                               turb_flux_uv(i, k, j) * dv_dx(i, k, j) +  &
                                               turb_flux_v2(i, k, j) * dv_dy(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dv_dz(i, k, j) +  &
                                               turb_flux_uw(i, k, j) * dw_dx(i, k, j) +  &
                                               turb_flux_vw(i, k, j) * dw_dy(i, k, j) +  &
                                               turb_flux_w2(i, k, j) * dw_dz(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_shear_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_shear_tend = ', q_sq_shear_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_shear


    subroutine Calc_q_sq_buoyancy (N_PBL3D, q_sq_buoyancy_tend, q_sq_tend, turb_flux_wtheta_v, dthv_dz, q_sq, l_master, mut, c1, c2, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_buoyancy_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq, l_master, turb_flux_wtheta_v
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: dthv_dz

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            q_sq_buoyancy_tend(i, k, j) = 2. * G_OVER_TREF * turb_flux_wtheta_v(i, k, j)
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) + ( (c1(k) * mut(i, j) + c2(k)) * q_sq_buoyancy_tend(i, k, j) / N_PBL3D )

            if ( DEBUG .and. i == I_TO_TEST .and. j == J_TO_TEST ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_buoyancy_tend = ', q_sq_buoyancy_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_buoyancy


    subroutine Fill_dissip_length_scale (l_opt, l_master, l_boulac, l_dissip, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      integer,                                     intent (in)  :: l_opt
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master, l_boulac
      real, dimension (its:ite, kts:kte, jts:jte), intent (out) :: l_dissip

      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      real :: l
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            if (l_opt < 3) then
              l_dissip(i, k, j) = l_master(i, k, j)
            else
                ! Boulac length scale is at faces
              l_dissip(i, k, j) = (l_boulac(i, k, j) + l_boulac(i, k + 1, j)) / 2.
            end if
          end do
        end do
      end do

    end subroutine Fill_dissip_length_scale


    subroutine Calc_q_sq_dissip (N_PBL3D, q_sq_dissip_tend, q_sq_tend, l_master, l_dissip, q_sq, mut, c1, c2, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, i_start, i_end, j_start, j_end, ktf)

      implicit none

      real,                                        intent (in)  :: N_PBL3D
      real, dimension (ims:ime, kms:kme, jms:jme), intent (out) :: q_sq_dissip_tend
      real, dimension (its:ite, kts:kte, jts:jte), intent (inout) :: q_sq_tend
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: q_sq
      real, dimension (ims:ime, kms:kme, jms:jme), intent (in)  :: l_master
      real, dimension (its:ite, kts:kte, jts:jte), intent (in)  :: l_dissip

      real, dimension (ims:ime, jms:jme), intent (in) :: mut
      real, dimension (kms:kme),          intent (in) :: c1, c2
      integer, intent (in) :: ims, ime, jms, jme, kms, kme,        &
                              its, ite, jts, jte, kts, kte,        &
                              i_start, i_end, j_start, j_end, ktf

        ! Local vars
      integer :: i, j, k
      logical, parameter :: DEBUG = .false.

      do j = j_start, j_end
        do k = kts, ktf
          do i = i_start, i_end
            if ( DEBUG ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq = ', q_sq(i, k, j)
              print *, '  l_master = ', l_master(i, k, j)
              print *, '  b1 = ', b_1
            end if
            q_sq_dissip_tend(i, k, j) = 2. * ( q_sq(i, k, j) ** 1.5 ) /  &
                                             ( b_1 * l_dissip(i, k, j) )
            q_sq_tend(i, k, j) = q_sq_tend(i, k, j) - ( (c1(k) * mut(i, j) + c2(k)) * q_sq_dissip_tend(i, k, j) / N_PBL3D )

            if ( DEBUG ) then
              print *, '------------------------------------------'
              print *, '  i =', i, 'j = ', j, 'k = ', k
              print *, '  q_sq_dissip_tend = ', q_sq_dissip_tend(i, k, j)
              print *, '  q_sq_tend = ', q_sq_tend(i, k, j)
            end if
          end do
        end do
      end do

    end subroutine Calc_q_sq_dissip


    subroutine Calc_metric (ph, phb, rdz, rdzw, dz, z_at_mass, z_at_w, ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, itf, jtf, ktf, i_start, i_end, j_start, j_end)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                                        !
      ! Purpose: Calculate metric terms for initialization of turbulence state !
      !                                                                        !
      ! Author: Timothy W. Juliano                                             !
      !                                                                        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      integer, intent(in) :: itf, jtf, ktf, i_start, i_end, j_start, j_end

      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: ph, phb
      real, dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rdz, rdzw, dz, &
                                                                 z_at_mass, z_at_w

        ! Local vars
      integer :: i, j, k

      do j = j_start, j_end

        do k = kts, kte
          do i = i_start, i_end
            z_at_w(i, k, j) = ( ph(i, k, j) + phb(i, k, j) ) / g
          end do
        end do

        do k = kts, ktf
          do i = i_start, i_end
            z_at_mass(i, k, j) = ( z_at_w(i, k + 1, j) + z_at_w(i, k, j) ) / 2.
          end do
        end do

        do k = kts, ktf
          do i = i_start, i_end
            rdzw(i, k, j) = 1.0 / ( z_at_w(i, k + 1, j) - z_at_w(i, k, j) )
            dz(i, k, j) = z_at_w(i, k + 1, j) - z_at_w(i, k, j)
          end do
        end do

        do k = kts + 1, ktf
          do i = i_start, i_end
            rdz(i, k, j) = 2.0 / ( z_at_w(i, k + 1, j) - z_at_w(i, k - 1, j) )
          end do
        end do

        do i = i_start, i_end
          rdz(i, kts, j) = 2. / (z_at_w(i, kts + 1, j) - z_at_w(i, kts, j))
        end do

      end do

    end subroutine Calc_metric


    subroutine Init_pbl3d (restart, pbl3d_opt, pbl3d_sfc_opt, pbl3d_l_opt, pbl3d_prog, pbl3d_constants, &
              pbl3d_nsteps, fnm, fnp, dnw, ph, phb, t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol, &
              q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip,                       &
              ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,               &
              its, ite, jts, jte, kts, kte, l_master, l_master_at_mass, turb_flux_u2,                     &
              turb_flux_v2, turb_flux_w2, turb_flux_uv, turb_flux_uw,                   &
              turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v,                     &
              turb_flux_wtheta_v, turb_flux_theta2_v, turb_flux_utheta,                 &
              turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqv,                        &
              turb_flux_vqv, turb_flux_wqv)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      ! Purpose: Initializes the PBL3D                       !
      !                                                      !
      ! Author: Pedro A. Jimenez & Timothy W. Juliano        !
      !                                                      !
      ! Method: 1) User selects set of closure constants     !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      logical, intent(in) :: restart
      integer, intent(in) :: pbl3d_opt, pbl3d_sfc_opt, pbl3d_l_opt, pbl3d_prog, pbl3d_nsteps
      character (len = *), intent(in) :: pbl3d_constants

        ! Grid vars
      real, dimension(kms:kme),                   intent(in)  :: fnm, fnp, dnw
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: ph, phb

        ! Thermo/kinematic state vars
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: t_2, qv, u, v, alt
      real, dimension(ims:ime, jms:jme),          intent(in)  :: ust, hfx, qfx, rmol
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: thetav

        ! Turb state vars
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: q_sq_prog, q_sq,     &
          q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: l_master, l_master_at_mass, &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv

      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      call Set_init_turb_state_driver (restart, pbl3d_constants,            &
          pbl3d_nsteps, pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt,  &
          ids, ide, jds, jde, kds, kde,                                     &
          ims, ime, jms, jme, kms, kme,                                     &
          its, ite, jts, jte, kts, kte,                                     &
          fnm, fnp, dnw, ph, phb,                                           &
          t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol,                  &
          q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear,              &
          q_sq_buoyancy, q_sq_dissip,                                       &
          l_master, l_master_at_mass,                                       &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,           &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                   &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,       &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,             &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv                       )

    end subroutine Init_pbl3d


    subroutine Set_init_turb_state_driver (restart, pbl3d_constants,            &
              pbl3d_nsteps, pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt,  &
              ids, ide, jds, jde, kds, kde,                                     &
              ims, ime, jms, jme, kms, kme,                                     &
              its, ite, jts, jte, kts, kte,                                     &
              fnm, fnp, dnw, ph, phb,                                           &
              t_2, qv, thetav, alt, u, v, ust, hfx, qfx, rmol,                  &
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
      integer :: i, j, k, itf, jtf, ktf
      integer, intent(in) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      real, dimension (kms:kme),                  intent(in)  :: fnm, fnp, dnw
      real, dimension(ims:ime, jms:jme),          intent(in)  :: ust, hfx, qfx, rmol
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in)  :: t_2, qv, u, v, &
                                                                 alt, ph, phb
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: thetav
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: q_sq_prog, q_sq,     &
          q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip
      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: l_master, l_master_at_mass, &
          turb_flux_u2, turb_flux_v2, turb_flux_w2, turb_flux_uv,          &
          turb_flux_uw, turb_flux_vw, turb_flux_utheta_v,                  &
          turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v,      &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta,            &
          turb_flux_uqv, turb_flux_vqv, turb_flux_wqv

        ! Local vars
      real, dimension (kts:kte) :: sqv
      real, dimension (ims:ime, kms:kme, jms:jme) :: rdz, rdzw, z_at_mass, z_at_w, dz
      real, dimension (its:ite, kts:kte, jts:jte) :: sm, sh, rif
      real, dimension (its:ite, kts:kte, jts:jte) :: du_dz, dv_dz, &
                                                     dthetav_dz, dqv_dz
      real, dimension (ims:ime, kms:kme, jms:jme) :: rho
      integer :: i_start, i_end, j_start, j_end, init_iter

        ! Defines mass point at faces indexes
      jtf = Min (jte, jde - 1)
      itf = Min (ite, ide - 1)
      ktf = Min (kte, kde - 1)

      if (DEBUG) write (OUTPUT_UNIT, *) "ktf = ", ktf, "itf = ", itf, "jtf = ", jtf , "kte = ", kte, "kde = ", kde, "kme = ", kme

      i_start = its
      i_end = Min (ite, ide - 1)
      j_start = jts
      j_end = Min (jte, jde - 1)

        ! Calc metric terms
      call Calc_metric (ph, phb, rdz, rdzw, dz, z_at_mass, z_at_w, ids, ide, jds, jde, kds, kde, &
          ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, itf, jtf, ktf, i_start, i_end, j_start, j_end)

      call Calc_thetav (t_2, qv, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, thetav)

        ! Calc vertical derivatives
      call Calc_du_dz (u, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, du_dz)

      call Calc_dv_dz (v, rdz, dnw, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dv_dz)

      call Calc_dmassvar_dz (thetav, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dthetav_dz)

      call Calc_dmassvar_dz (qv, rdz, ide, jde, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, dqv_dz)

        ! Call MY model to set initial state
      call Set_init_turb_state_my (restart, pbl3d_constants, pbl3d_nsteps, pbl3d_opt, pbl3d_prog, pbl3d_sfc_opt, pbl3d_l_opt, ids, ide, &
          jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte, itf, jtf, ktf, i_start, i_end, j_start, j_end, &
          fnm, fnp, dnw, ph, phb, rdzw, z_at_mass, dz, t_2, qv, alt, u, v, ust, hfx, qfx, rmol, rho, du_dz, dv_dz, dthetav_dz, dqv_dz, &
          q_sq, q_sq_prog, q_sq_vdiff, q_sq_hdiff, q_sq_shear, q_sq_buoyancy, q_sq_dissip, l_master, l_master_at_mass, turb_flux_u2, turb_flux_v2, &
          turb_flux_w2, turb_flux_uv, turb_flux_uw, turb_flux_vw, turb_flux_utheta_v, turb_flux_vtheta_v, turb_flux_wtheta_v, turb_flux_theta2_v, &
          turb_flux_utheta, turb_flux_vtheta, turb_flux_wtheta, turb_flux_uqv, turb_flux_vqv, turb_flux_wqv)

      if (DEBUG) print *, 'Done Set_init_turb_state_driver...'

    end subroutine Set_init_turb_state_driver

  end module module_pbl3d
