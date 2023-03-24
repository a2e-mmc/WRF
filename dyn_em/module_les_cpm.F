! WRF:MODEL_LAYER:PHYSICS

MODULE module_les_cpm

  USE module_model_constants    
  
#ifdef DM_PARALLEL
  USE module_dm
#endif
  
  
CONTAINS

  !!====================================================================
  
  SUBROUTINE force_down_meso_pblh( meso_pblh, pblh,               &
                                   ids, ide, jds, jde, kds, kde,  &
                                   ims, ime, jms, jme, kms, kme,  &
                                   its, ite, jts, jte, kts, kte   )
    
    !-------------------------------------------------------------------
    !
    ! Assign PBL height from mesoscale PBL scheme (pblh) to variable
    ! meso_pblh which is forced down to nested LES domains using the
    ! CPM option to base the perturbation height on the PBL height.
    !
    !-------------------------------------------------------------------
   
    
    IMPLICIT NONE

    INTEGER, INTENT( IN ) ::            &                                    ! Starting and ending index values for
         ids, ide, jds, jde, kds, kde,  &                                    ! domain (d),
         ims, ime, jms, jme, kms, kme,  &                                    ! memory (m) and
         its, ite, jts, jte, kts, kte                                        ! tile (t) dimensions.
    
    REAL, DIMENSION( ims:ime, jms:jme ), INTENT( IN   ) :: pblh              ! PBL height.                            [m]
    REAL, DIMENSION( ims:ime, jms:jme ), INTENT( OUT  ) :: meso_pblh         ! PBL height for CPM scheme.             [m] 
    
    INTEGER :: i,j                                                           ! Loop index variables.
    INTEGER :: i_start, i_end, j_start, j_end                                ! Loop start and end variables.                

    
    ! End declarations
    !-------------------------------------------------------------------
    
    
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)
    
    DO j = j_start, j_end
       DO i = i_start,i_end
          
          meso_pblh(i,j) = pblh(i,j)
                    
       END DO
    END DO
    
    
  END SUBROUTINE force_down_meso_pblh
  
  
  !!====================================================================
  
  SUBROUTINE force_down_meso_hfx( meso_hfx, hfx, rho,            &
                                  ids, ide, jds, jde, kds, kde,  &
                                  ims, ime, jms, jme, kms, kme,  &
                                  its, ite, jts, jte, kts, kte   )
    
    !-------------------------------------------------------------------
    !
    ! Assign surface sensible heat flux from mesoscale surface layer
    ! scheme (hfx) to variable meso_hfx, which is forced down to nested
    ! LES domains using the CPM, with cpm_opt=4 (variance scaling
    ! for unstable conditions). Units converted from W/m^2 to K m/s.
    !
    !-------------------------------------------------------------------

    
    IMPLICIT NONE
    
    INTEGER, INTENT( IN ) ::            &                                    ! Starting and ending index values for
         ids, ide, jds, jde, kds, kde,  &                                    ! domain (d),
         ims, ime, jms, jme, kms, kme,  &                                    ! memory (m) and
         its, ite, jts, jte, kts, kte                                        ! tile (t) dimensions.

    REAL, DIMENSION( ims:ime, jms:jme ),         INTENT( OUT  ) :: meso_hfx  ! Surface sensible heat flux for CPM.               [K*m/s]
    REAL, DIMENSION( ims:ime, jms:jme ),         INTENT( IN   ) :: hfx       ! Surface sensible heat flux.                       [W/m^2]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme), INTENT( IN   ) :: rho       ! Air density.                                      [kg/m^3]
    
    INTEGER :: i,j                                                           ! Loop index variables
    INTEGER :: i_start, i_end, j_start, j_end                                ! Loop start and end variables                

    
    ! End declarations
    !-------------------------------------------------------------------
        
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)

    DO j = j_start, j_end
       DO i = i_start,i_end
          
          meso_hfx(i,j) = hfx(i,j)/(rho(i,kts,j)*cp)
          
       END DO
    END DO
    
    
  END SUBROUTINE force_down_meso_hfx

  
  !!====================================================================

  SUBROUTINE calc_cpm( cpm_opt, cpm_lim_z, cpm_amp,               &
                       cpm_pec, cpm_rim, cpm_dt, cpm_lambda,      &
                       cpm_nb, cpm_sb, cpm_eb, cpm_wb,            &
                       cpm_ngc_h, cpm_ncells_h, cpm_noff_h,       &
                       cpm_ngc_v, cpm_noff_v, spec_bdy_width,     &
                       prttms, prtdt, prtseed, cpm_3d,            &
                       cpm_meso_pblh, meso_pblh, cpm_pblh_check,  &
                       cpm_meso_hfx, meso_hfx, spec_hfx,          &
                       t, u, v, w, rho,                           &
                       z, dx, dt, fnm, fnp,                       &
                       ids, ide, jds, jde, kds, kde,              &
                       ims, ime, jms, jme, kms, kme,              &
                       its, ite, jts, jte, kts, kte               )

    !-------------------------------------------------------------------
    !
    ! Cell Perturbation Method (CPM) module. Applies instantaneous 
    ! stochastic perturbations, using a uniform random number 
    ! distribution scaled by a target amplitude window, to the
    ! potential temperature (cpm_opt = 1-4) or vertical velocity
    ! (cpm_opt = 5) field. Many perturbation characteristics can be 
    ! specified in the namelist.input file. For a complete description,
    ! see Registry/registry.les_cpm.
    
    ! Once the number of grid points per cell in the horizontal 
    ! (cpm_ngc_h) and vertical (cpm_ngc_v) directions, number of cells 
    ! per boundary (cpm,_ncells_h) and any offset from the outflow edge
    ! (cpm_noff_h) are determined, the following values of cpm_opt
    ! can be set to specify the target amplitude window, refresh time,
    ! and which boundaries to perturb.
    !
    ! cpm_opt = 1   Specify the thermal perturbation amplitude
    !               (cpm_amp > 0.0), refresh time (cpm_dt > 0.0) and
    !               which boundaries to perturb (cpm_*b = 1 to perturb,
    !               0 to ignore) ) in the namelist.
    !
    ! cpm_opt = 2   Use Eckert number scaling for the thermal 
    !               perturbation amplitude. Developed for neutral
    !               stability conditions, but can be used for all
    !               stabilities. Refresh time and which boundary(ies)
    !               to perturb can either be specified in the 
    !               namelist (cpm_dt > 0.0 and one or more of
    !               cpm_*b > 0), or determined automatically using the
    !               advective timescale at the height of the slab
    !               (cpm_dt = 0.0), and the wind direction at the
    !               estimated height of the geostrophic wind
    !               (1.2*zi, with zi the PBL height), specified by
    !               cpm_*b all = 0 in the namelist. Perturbations
    !               will be applied at any boundary with a
    !               component of the geostrophic wind (averaged across
    !               the domain perimeter) isoriented into the domain.
    !
    ! cpm_opt = 3   Use modified Richardson # scaling for the thermal
    !               perturbation amplitude. Developed for stable
    !               conditions. Refresh time and boundaries to perturb 
    !               set as for cpm_opt = 2.
    !     
    ! cpm_opt = 4   Use convective variance scaling for the thermal
    !               perturbation amplitude. Developed for unstable
    !               conditions. Refresh time and boundaries to perturb
    !               set as for cpm_opt = 2.
    !
    ! cpm_opt = 5   Perturb vertical velocity with value specified in
    !               namelist. Refresh time and boundaries to perturb
    !               set as for cpm_opt = 2.
    !            
    ! Potential temperature options are based on the procedures
    ! outlined in the following three manuscripts:
    ! 
    ! Munoz-Esparza, D. et al. 2014: Bridging the Transition from    
    ! Mesoscale to Microscale Turbulence in Numerical Weather Prediction
    ! Models. Boundary-Layer Meteorology, 153:409-440.  
    ! doi:10.1007/s10546-014-9956-9                                
    !                                                              
    ! Munoz-Esparza, D. et al. 2015: A stochastic perturbation method  
    ! to generate inflow turbulence in large-eddy simulation models:
    ! Application to neutrally stratified atmospheric boundary layers.    
    ! Physics of Fluids, 27:035102. doi:10.1063/1.4913572                                        
    !                                                              
    ! Munoz-Esparza, D. and Kosovic, B. 2018: Generation of inflow      
    ! turbulence in large-eddy simulations of nonneutral atmospheric   
    ! boundary layers with the cell perturbation method. Monthly     
    ! Weather Review, 146:1889-1909. doi:10.1175/MWR-D-18-0077.1              
    !
    ! Vertical velocity technique outlined in the following manuscript,
    ! however herein the w values are perturbed directly, rather than
    ! via the w tendency, as in the manuscript.
    !
    ! Mazzaro, L.J. et al. 2019: Random Force Perturbations: A New   
    ! Extension of the Cell Perturbation Method for Turbulence Generation
    ! in Multiscale Atmospheric Boundary Layer Simulations. Journal of          
    ! Advances in Modeling Earth Systems, 11, 2311–2329.
    ! doi:10.1029/2019MS001608        
    !
    !-------------------------------------------------------------------

    
    IMPLICIT NONE
    
#ifdef DM_PARALLEL
    INCLUDE 'mpif.h'
#endif
    
    INTEGER, INTENT( IN ) ::            &                                      ! Starting and ending index values for
         ids, ide, jds, jde, kds, kde,  &                                      ! domain (d),
         ims, ime, jms, jme, kms, kme,  &                                      ! memory (m) and
         its, ite, jts, jte, kts, kte                                          ! tile (t) dimensions.
    
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT ) :: t         ! Potential temperature.                       [K]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT ) :: cpm_3d    ! Instantaneous 3D perturbation value.  [K or m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: u         ! Zonal wind component.                      [m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: v         ! Meridional wind component.                 [m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT ) :: w         ! Vertical wind component.                   [m/s]
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: rho
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN    ) :: z         ! Height above surface @ u,v,t points.         [m]
    REAL, DIMENSION( ims:ime, jms:jme ),          INTENT( IN    ) :: meso_pblh ! Mesoscale PBL height.                        [m]
    REAL, DIMENSION( ims:ime, jms:jme ),          INTENT( IN    ) :: meso_hfx  ! Mesoscale surface sensible heat flux.    [K*m/s] 
    REAL, DIMENSION( kms:kme ),                   INTENT( INOUT ) :: prttms    ! Time since last perturbations applied.       [s]
    REAL, DIMENSION( kms:kme ),                   INTENT( INOUT ) :: prtdt     ! Time interval to refresh perturbations.      [s]
    REAL, DIMENSION( kms:kme ),                   INTENT( IN    ) :: fnm       ! Vertical interpolation coefficient.
    REAL, DIMENSION( kms:kme ),                   INTENT( IN    ) :: fnp       ! Vertical interpolation coefficient.

    INTEGER, DIMENSION( kms:kme ),                INTENT( INOUT ) :: prtseed   ! 1D variable to store variable-length seed array. Dimensioned as
                                                                               ! kms:kme to avoid the need to define a new dimension variable.

    INTEGER, INTENT( IN    ) :: cpm_opt                         ! Which perturbation option to use.
    INTEGER, INTENT( IN    ) :: cpm_nb, cpm_sb, cpm_eb, cpm_wb  ! Specified boundaries to apply perturbations to (1 to perturb, 0 to ignore).
    INTEGER, INTENT( IN    ) :: cpm_ngc_h, cpm_ngc_v            ! Number of grid points per CPM cell in horizontal and vertical directions.
    INTEGER, INTENT( IN    ) :: cpm_ncells_h                    ! Number of CPM cells in both horizontal directions.
    INTEGER, INTENT( IN    ) :: cpm_noff_h                      ! Number of grid points from any outflow boundary(ies) to end perturbations.
    INTEGER, INTENT( IN    ) :: cpm_noff_v                      ! Number of grid cells from bottom boundary to start perturbations.
    INTEGER, INTENT( IN    ) :: spec_bdy_width                  ! Width which nested domain solution is set or blended with bounding domain solution.
    INTEGER, INTENT( IN    ) :: cpm_meso_pblh                   ! Whether to use mesoscale PBL height (1) or not (0).
    INTEGER, INTENT( IN    ) :: cpm_meso_hfx                    ! Whether to use mesoscale sensible heat flux (1) or specified value (0) for cpm_opt=4.
    REAL, INTENT( IN    ) :: cpm_pblh_check                     ! Time interval to check if mesoscale pblh has changed.                                         [s]
    REAL, INTENT( IN    ) :: cpm_lim_z                          ! Max height to apply cpm when cpm_meso_pblh=0; lower bound when cpm_meso_pblh=1.               [m]
    REAL, INTENT( IN    ) :: cpm_amp                            ! Target amplitude window for perturbations.                                             [K or m/s]
    REAL, INTENT( IN    ) :: cpm_pec                            ! Perturbation Eckert number value.
    REAL, INTENT( IN    ) :: cpm_rim                            ! Modified Richardson number value.
    REAL, INTENT( IN    ) :: cpm_dt                             ! Time to reapply perturbations at all heights if > 0; = 0.0 to automate.                       [s]
    REAL, INTENT( IN    ) :: cpm_lambda                         ! Account for width which LES solution is set to or blended with the bounding solution.
    REAL, INTENT( IN    ) :: spec_hfx                           ! User specified surface sensible heat flux.                                                [W/m^2]
    REAL, INTENT( IN    ) :: dt                                 ! Model time step.                                                                              [s]
    REAL, INTENT( IN    ) :: dx                                 ! Horizontal grid spacing.                                                                      [m]

!   LOCAL VARIABLES
    
    INTEGER :: i, j, k, slab_k, m                               ! Loop index variables.
    INTEGER :: i_start, i_end, j_start, j_end                   ! Loop endpoints.
    INTEGER :: north, south, east, west                         ! Designate which boundaries to apply perturbations to.
    INTEGER :: i_seed, seedsum                                  ! Size of seed and variable for its summation.
    INTEGER :: ni, nj                               !           ! Number of perturbation cells in horizontal directions (for random number query).
    INTEGER :: n_slabs_k                                        ! Number of perturbations cells in vertical direction.
    INTEGER :: k_slab_start, k_slab_end                         ! Vertical indices defining endpoints of perturbation slabs.
    INTEGER :: k_geo                                            ! k-index at which geostrophic winds are computed.
    INTEGER :: idum, jdum                                       ! Dummy index values for intrinsic WRF functions.

    REAL    :: pertmag                                          ! Magnitude of perturbation to scale amplitude window.                            [K]
    REAL    :: pblhdavg                                         ! Domain perimeter average PBL height .                                           [m]
    REAL    :: hfxdavg                                          ! Domain perimeter average sensible heat flux height.                         [K*m/s]
    REAL    :: h_geo                                            ! Height at which to estimate geostrophic winds.                                  [m]
    REAL    :: ugeodavg, vgeodavg, wsgeodavg                    ! Domain perimeter average geostrophic wind components and speed.               [m/s]
    REAL    :: wsslabdavg                                       ! Domain perimeter average wind speed within slab.                              [m/s]
    REAL    :: angleslabdavg                                    ! Domain perimeter average wind direction within slab.                          [deg]
    REAL    :: tslabdavg                                        ! Domain perimeter average temperature within slab.                               [K]
    REAL    :: zslabdavg                                        ! Domain perimeter average height of slab.                                        [m] 
    REAL    :: uslabdavg, vslabdavg                             ! Domain perimeter average horizontal wind components within slab.              [m/s]
    REAL    :: duslabdavg, dvslabdavg                           ! Domain perimeter average vertical differences of wind component across slab.  [m/s]
    REAL    :: dtslabdavg                                       ! Domain perimeter average vertical differences of temperature across slab.       [K]
    REAL    :: dzslabdavg                                       ! Domain perimeter average depth of slab.                                         [m]
    
    REAL    :: num_pts_in_sum_1, num_pts_in_sum_2               ! Number of elements on each patch to compute lateral boundary averages.
    REAL    :: counter_sum_1, counter_sum_2                     ! Number of elements in domain to compute lateral boundary averages.
    REAL    :: u_slab_top, u_slab_bot                           ! Zonal velocity on top and bottom of slabs, to compute difference.             [m/s]
    REAL    :: v_slab_top, v_slab_bot                           ! Meridional velocity on top and bottom of slabs, to compute difference.        [m/s]
    REAL    :: t_slab_top, t_slab_bot                           ! Temperature on top and bottom of slabs, to compute difference.                  [K]
    REAL    :: z_slab_top, z_slab_bot                           ! Heigt of top and bottom of slabs, to compute difference.                        [m]
    
    ! The summation variables below used to compute mean values around lateral edges on patches (p) and domain (d) must be double precision
    ! to ensure obtaining the same value using different processor splits

    DOUBLE PRECISION :: pblhpsum, pblhdsum                      ! For PBL height                                                                [m]
    DOUBLE PRECISION :: hfxpsum, hfxdsum                        ! For sensible heat flux                                                    [K*m/s]
    DOUBLE PRECISION :: ugeopsum, vgeopsum                      ! For geostrophic wind components                                             [m/s]
    DOUBLE PRECISION :: ugeodsum, vgeodsum                      ! "                                                                           [m/s]
    DOUBLE PRECISION :: uslabpsum, vslabpsum                    ! For slab-average wind components                                            [m/s]
    DOUBLE PRECISION :: uslabdsum, vslabdsum                    ! "                                                                           [m/s]
    DOUBLE PRECISION :: tslabpsum, tslabdsum                    ! For temperature                                                               [K]
    DOUBLE PRECISION :: zslabpsum, zslabdsum                    ! For height                                                                    [m]
    DOUBLE PRECISION :: duslabpsum, duslabdsum                  ! For vertical differences of wind components                                 [m/s]
    DOUBLE PRECISION :: dvslabpsum, dvslabdsum                  ! "                                                                           [m/s]
    DOUBLE PRECISION :: dtslabpsum, dtslabdsum                  ! For vertical differences of temperature                                       [K]
    DOUBLE PRECISION :: dzslabpsum, dzslabdsum                  ! For vertical differences of height                                            [m]
    DOUBLE PRECISION :: rhoslabpsum, rhoslabdsum                ! For air density                                                          [kg/m^3]
    
    REAL    :: pio2 = piconst/2.0                               !
    REAL    :: pio4 = piconst/4.0                               !       
    REAL    :: pblh_fac = 0.75                                  ! Percentage of PBL height to assign perturbations to.                           
    REAL    :: h_geo_fac = 1.2                                  ! Percentage of maximum perturbation height to retreive geostrophic winds.
    REAL    :: max_pert_z                                       ! Maximum physical height of slab average to which to apply perturbations.      [m]    
    REAL    :: zi                                               ! PBL height used by convective variance formulation.                           [m]
    REAL    :: thetastar                                        ! Convective temperature scale for convective variance formulation.             [K]
    REAL    :: f1, f2, f3, f4                                   ! Dummy variables for convective variance scaling.

    INTEGER, dimension( : ), allocatable :: seed                                    ! random number seed
    REAL, dimension( :, : ), allocatable :: perts_wb, perts_eb, perts_sb, perts_nb  ! random number arrays for lateral boundaries
  
#ifdef DM_PARALLEL
    INTEGER :: ierr
    INTEGER :: tag
    INTEGER :: master
    INTEGER :: status(MPI_STATUS_SIZE)
#endif

#ifdef DM_PARALLEL
    master = 0
    tag    = 0
#endif

    ! End declarations
    !-------------------------------------------------------------------
    
    IF ( (cpm_opt .GT. 1 ) .AND. ( cpm_amp .GT. 0.0 ) ) THEN
       
       IF (cpm_opt .LT. 5 ) THEN
          
          PRINT*,'---------- WARNING ----------'
          PRINT*,'Specifying cpm_amp > 0.0 in namelist while also selecting automated amplitude, cpm_opt = 2, 3 or 4. Stopping.'
          
          STOP
          
       ENDIF
       
    ENDIF
    
    IF (cpm_opt .EQ. 5 ) THEN
       
       IF ( cpm_amp .GT. 0.0 ) THEN
          
          pertmag = cpm_amp
          
       ELSE
          
          PRINT*,'---------- WARNING ----------'
          PRINT*,'Specifying w-perturbations (cpm_opt = 5) with zero amplitude (cpm_amp = 0.0). Stopping.'
          
          STOP

       ENDIF
       
    ENDIF
    
    idum = 1
    jdum = 1
    
    north = 0
    south = 0
    east = 0
    west = 0
    
    i_start = its
    i_end   = MIN(ite,ide-1)
    j_start = jts
    j_end   = MIN(jte,jde-1)
 
    n_slabs_k = kde/cpm_ngc_v
    
! Loop over all vertical slabs to see if slab should be perturbed (PBL grew) or have perturbation value refreshed (refresh or advective timescale exceeded)

    DO slab_k = 1, n_slabs_k     ! Outer k-loop over all slabs spanning the domain depth that might receive perturbations
       
       prttms(slab_k) = prttms(slab_k) + dt
       
       IF ( prttms(slab_k) .GE. prtdt(slab_k) ) THEN !prdt should be zero at startup
          
          prttms(slab_k) = dt

! Check each cpm_pblh_check seconds to see if slab is below the maximum perturbation height, max_pert_z (i.e. the PBL grew)

          prtdt(slab_k) = cpm_pblh_check 

!          PRINT*,'Computing new perturbations, slab_k = ',slab_k

! Determine the starting and ending vertical k-index values within the slab

          k_slab_start = MIN( (slab_k - 1)*cpm_ngc_v + 1, kde-1 ) + cpm_noff_v        !lowest k-value in slab
          k_slab_end = MIN( k_slab_start + cpm_ngc_v - 1, kde-1 )                     !highest k-value in slab   

!          PRINT*,'k_slab_start,k_slab_end ',k_slab_start,k_slab_end
          
! Only compute and apply perturbations for slabs with average z-coordinate value < the maximum perturbation, height max_pert_z
                    
! Find max_pert_z

          IF (cpm_meso_pblh .EQ. 0 ) THEN !Set max perturbation height to value specified in the namelist

             max_pert_z = cpm_lim_z
             
          ENDIF
          
          IF (cpm_meso_pblh .EQ. 1 ) THEN !Use mesoscale PBLH to set max_pert_z. Use cpm_lim_z as a minimum value in case PBL shrinks too low.
             
             ! Compute mean value of meso_pblh around the domain perimeter. Ghost point values not needed, as meso_pblh is not influenced by perturbations.
             ! Therefore we use the first grid point value within the domain boundary as a good approximation of the ghost point value.

             pblhpsum = 0.0  
             num_pts_in_sum_1 = 0.0
             
             IF ( its .EQ. ids ) THEN
                
                DO j = j_start,j_end
                   
                   pblhpsum = pblhpsum + meso_pblh(ids,j) 

                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0

                END DO
                
             ENDIF
             
             IF ( ite .EQ. ide ) THEN
                
                DO j = j_start,j_end
                   
                   pblhpsum = pblhpsum + meso_pblh(ide-1,j)
                   
                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0

                END DO
                
             ENDIF
             
             IF ( jts .EQ. jds ) THEN
                
                DO i = i_start,i_end
                   
                   pblhpsum = pblhpsum + meso_pblh(i,jds)
                   
                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0

                END DO
                
             ENDIF
             
             IF ( jte .EQ. jde ) THEN
                
                DO i = i_start,i_end
                   
                   pblhpsum = pblhpsum + meso_pblh(i,jde-1)
                   
                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0

                END DO
                
             ENDIF
             
             pblhdsum = pblhpsum
             
             counter_sum_1 = DBLE(num_pts_in_sum_1)
                          
#ifdef DM_PARALLEL
             
             pblhdsum  = wrf_dm_sum_real8(pblhpsum)
             
             counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
             
#endif
             
             pblhdavg = pblhdsum/counter_sum_1       
             
             max_pert_z = MAX(cpm_lim_z,pblh_fac*pblhdavg) 
             
          ENDIF
          
! Find mean z-value of slab to determine if slab is below max_pert_z, and to use if convective scaling is selected (cpm_opt = 4).

          zslabpsum = 0.0d0
          num_pts_in_sum_1 = 0.0
          
          IF ( its .EQ. ids ) THEN
             
             DO j = j_start,j_end              

                DO k = k_slab_start,k_slab_end
                   
                   zslabpsum = zslabpsum + z(ids-1,k,j)
                   
                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                   
                END DO
                
             END DO
             
          END IF
          
          IF ( ite .EQ. ide ) THEN
             
             DO j = j_start,j_end
 
                DO k = k_slab_start,k_slab_end

                   zslabpsum = zslabpsum + z(ide,k,j)
                   
                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0

                END DO
                
             END DO
             
          END IF
                    
          IF ( jts .EQ. jds ) THEN
             
             DO i = i_start,i_end
                
                DO k = k_slab_start,k_slab_end

                   zslabpsum = zslabpsum + z(i,k,jds-1)

                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                   
                END DO

             END DO

          END IF
          
          IF ( jte .EQ. jde ) THEN
             
             DO i = i_start,i_end
                
                DO k = k_slab_start,k_slab_end

                   zslabpsum = zslabpsum + z(i,k,jde)

                   num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                   
                END DO

             END DO
           
          END IF

          zslabdsum = zslabpsum

          counter_sum_1 = num_pts_in_sum_1
                    
#ifdef DM_PARALLEL
          
          zslabdsum  = wrf_dm_sum_real8(zslabpsum)

          counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
          
#endif

          zslabdavg = zslabdsum/counter_sum_1
          
          IF ( zslabdavg .LE. max_pert_z ) THEN  ! If the average z value for the slab is < max_pert_z, perturb slab  
             
             ! Set boundary(ies) to perturb from the namelist if one or more of cpm_*b are specified ( = 1 ).
             
             IF (cpm_opt .EQ. 1 ) THEN ! Perturbation magnitude specified in the namelist
                
                pertmag = cpm_amp
                
             ENDIF
             
             IF ( ( cpm_nb + cpm_sb + cpm_eb + cpm_wb ) .GT. 0 ) THEN !Boundary(ies) to perturb specified in namelist.
                
                north = cpm_nb
                south = cpm_sb
                east = cpm_eb
                west = cpm_wb
                
             ENDIF
             
             IF ( ( ( cpm_nb + cpm_sb + cpm_eb + cpm_wb ) .EQ. 0 ) .OR. ( cpm_opt .EQ. 2 ) ) THEN
      
                ! If boundary(ies) to perturb are not specified (all cpm_*b = 0), and/or the Eckert # formulation is selected,
                ! compute a "geostrophic height", h_geo, and obtain the u and v wind components at that height, to determine
                ! which boundaries to perturb (those with wind components into the domain), and/or the perturbation amplitude.
                
                h_geo = h_geo_fac*max_pert_z    
                
                ugeopsum = 0.0
                vgeopsum = 0.0
                num_pts_in_sum_1 = 0.0
                
                IF ( its .EQ. ids ) THEN
                   
                   DO j = j_start,j_end
                      
                      k_geo = kde-1     ! Max possible value if loop below does not find it (PBL deep relative to domain)
                      
                      DO k=kts,kde-2                   
                         
                         IF ( ( z(ids-1,k,j) .LE. h_geo ) .AND. ( z(ids-1,k+1,j) .GT. h_geo ) )  k_geo = k
                         
                      END DO
                      
                      ugeopsum = ugeopsum + u(ids-1,k_geo,j)
                      vgeopsum = vgeopsum + v(ids-1,k_geo,j)
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( ite .EQ. ide ) THEN
                   
                   DO j = j_start,j_end
                      
                      k_geo = kde-1 
                      
                      DO k=kts,kde-2                   
                         
                         IF ( ( z(ide,k,j) .LE. h_geo ) .AND. ( z(ide,k+1,j) .GT. h_geo ) )  k_geo = k
                         
                      END DO
                      
                      ugeopsum = ugeopsum + u(ide+1,k_geo,j)
                      vgeopsum = vgeopsum + v(ide,k_geo,j)
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( jts .EQ. jds ) THEN
                   
                   DO i = i_start,i_end
                      
                      k_geo = kde-1
                      
                      DO k=kts,kde-2                   
                         
                         IF ( ( z(i,k,jds-1) .LE. h_geo ) .AND. ( z(i,k+1,jds-1) .GT. h_geo ) )  k_geo = k
                         
                      END DO
                      
                      ugeopsum = ugeopsum + u(i,k_geo,jds-1)
                      vgeopsum = vgeopsum + v(i,k_geo,jds-1)
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( jte .EQ. jde ) THEN
                   
                   DO i = i_start,i_end
                      
                      k_geo = kde-1
                      
                      DO k=kts,kde-2                   
                         
                         IF ( ( z(i,k,jde) .LE. h_geo ) .AND. ( z(i,k+1,jde) .GT. h_geo ) )  k_geo = k
                         
                      END DO
                      
                      ugeopsum = ugeopsum + u(i,k_geo,jde)
                      vgeopsum = vgeopsum + v(i,k_geo,jde+1)
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                ugeodsum = ugeopsum
                
                vgeodsum = vgeopsum
                
                counter_sum_1 = num_pts_in_sum_1
                
#ifdef DM_PARALLEL
                
                ugeodsum  = wrf_dm_sum_real8(ugeopsum)
                
                vgeodsum  = wrf_dm_sum_real8(vgeopsum)
                
                counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
                
#endif
                
                ugeodavg = ugeodsum/counter_sum_1
                
                vgeodavg = vgeodsum/counter_sum_1

                IF ( ( cpm_nb + cpm_sb + cpm_eb + cpm_wb ) .EQ. 0 ) THEN
                   
                   ! Determination of which boundary(ies) to perturb based on geostrophic wind components
                   
                   IF ( vgeodsum .LT. 0.0 ) north = 1
                   IF ( vgeodsum .GT. 0.0 ) south = 1
                   IF ( ugeodsum .LT. 0.0 ) east =  1
                   IF ( ugeodsum .GT. 0.0 ) west = 1
                   
                ENDIF
                
                wsgeodavg = sqrt( ugeodavg*ugeodavg + vgeodavg*vgeodavg ) !Needed for Eckert # scaling
                
             END IF ! Automatic determination of inflow boundaries and/or geostrophic wind speed for Eckert # scaling
             
             IF (cpm_opt .EQ. 2) THEN !Eckert # formulation
                
                pertmag = (wsgeodavg*wsgeodavg)/(cpm_pec*cp)
                
             ENDIF 
             
             IF (cpm_dt .GT. 0.0 ) THEN  !If cpm_dt > 0 is specified in the namelist, assign it here and do not compute it.
                
                prtdt(slab_k) = cpm_dt
                
             ELSE ! Find lateral average wind speed in the slab to determine advective timescale for refresh interval
                
                uslabpsum = 0.0
                vslabpsum = 0.0
                num_pts_in_sum_1 = 0.0
                
                IF ( its .EQ. ids ) THEN
                   
                   DO j = j_start,j_end
                      
                      DO k = k_slab_start, k_slab_end
                         
                         uslabpsum = uslabpsum + u(ids-1,k,j)
                         vslabpsum = vslabpsum + v(ids-1,k,j)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( ite .EQ. ide ) THEN
                      
                   DO j = j_start,j_end
                      
                      DO k = k_slab_start, k_slab_end                     
                         
                         uslabpsum = uslabpsum + u(ide+1,k,j)
                         vslabpsum = vslabpsum + v(ide,k,j)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jts .EQ. jds ) THEN
                   
                   DO i = i_start,i_end
                      
                      DO k = k_slab_start, k_slab_end                     
                         
                         uslabpsum = uslabpsum + u(i,k,jds-1)
                         vslabpsum = vslabpsum + v(i,k,jds-1)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jte .EQ. jde ) THEN
                   
                   DO i = i_start,i_end
                      
                      DO k = k_slab_start, k_slab_end                     
                         
                         uslabpsum = uslabpsum + u(i,k,jde)
                         vslabpsum = vslabpsum + v(i,k,jde+1)
                            
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                uslabdsum  = uslabpsum
                
                vslabdsum  = vslabpsum
                
                counter_sum_1 = num_pts_in_sum_1
                                
#ifdef DM_PARALLEL
                
                uslabdsum  = wrf_dm_sum_real8(uslabpsum)
                   
                vslabdsum  = wrf_dm_sum_real8(vslabpsum)
                
                counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
                
#endif
                
                uslabdavg = uslabdsum/counter_sum_1
                
                vslabdavg = vslabdsum/counter_sum_1
                
                wsslabdavg = sqrt( uslabdavg*uslabdavg + vslabdavg*vslabdavg )
                
                angleslabdavg = atan( abs( uslabdavg )/abs( vslabdavg ) )
                
                IF (angleslabdavg .GT. pio4 ) angleslabdavg = pio2 - angleslabdavg 
                
                IF (wsslabdavg .EQ. 0.0) THEN !Check for division by zero
                   
                   PRINT*,'something wrong calculating wind speed scaling in calc_cpm. Stopping.'
                   
                   STOP
                   
                ENDIF

                prtdt(slab_k) = (1.0/cos(angleslabdavg))*(cpm_ngc_h*cpm_ncells_h - cpm_lambda*spec_bdy_width)*dx/wsslabdavg
                
             ENDIF ! cpm_dt .EQ. 0; determine refresh time automatically using advective timescale
             
             IF ( cpm_opt .EQ. 3 ) THEN !Ri # formulation (stable)
                
                ! First differences across the slab (slab depth may be greater than one grid point)
                
                tslabpsum = 0.0
                duslabpsum = 0.0
                dvslabpsum = 0.0
                dtslabpsum = 0.0
                dzslabpsum = 0.0
                
                num_pts_in_sum_1 = 0.0
                num_pts_in_sum_2 = 0.0
                
                IF ( its .EQ. ids ) THEN
                   
                   DO j = j_start,j_end
                      
                      IF ( k_slab_start .EQ. 1 ) THEN !Linearly project to surface using two lowest model values (k=1,2)
                         
                         u_slab_bot = u(ids-1,1,j) - ( ( u(ids-1,2,j) - u(ids-1,1,j)  )/( z(ids,2,j) - z(ids,1,j)) )*z(ids,1,j)
                         
                         v_slab_bot = v(ids-1,1,j) - ( ( v(ids-1,2,j) - v(ids-1,1,j)  )/( z(ids,2,j) - z(ids,1,j)) )*z(ids,1,j)
                         
                         t_slab_bot = t(ids-1,1,j) - ( ( t(ids-1,2,j) - t(ids-1,1,j)  )/( z(ids,2,j) - z(ids,1,j)) )*z(ids,1,j)
                         
                         z_slab_bot = 0.0
                         
                      ELSE
                         
                         u_slab_bot = fnm(k_slab_start)*u(ids-1,k_slab_start,j) + fnp(k_slab_start)*u(ids-1,k_slab_start-1,j)
                         
                         v_slab_bot = fnm(k_slab_start)*v(ids-1,k_slab_start,j) + fnp(k_slab_start)*v(ids-1,k_slab_start-1,j)
                         
                         t_slab_bot = fnm(k_slab_start)*t(ids-1,k_slab_start,j) + fnp(k_slab_start)*t(ids-1,k_slab_start-1,j)
                         
                         z_slab_bot = fnm(k_slab_start)*z(ids-1,k_slab_start,j) + fnp(k_slab_start)*z(ids-1,k_slab_start-1,j)
                         
                      ENDIF
                      
                      u_slab_top = fnm(k_slab_end+1)*u(ids-1,k_slab_end+1,j) + fnp(k_slab_end+1)*u(ids-1,k_slab_end,j) 
                      
                      v_slab_top = fnm(k_slab_end+1)*v(ids-1,k_slab_end+1,j) + fnp(k_slab_end+1)*v(ids-1,k_slab_end,j) 
                      
                      t_slab_top = fnm(k_slab_end+1)*t(ids-1,k_slab_end+1,j) + fnp(k_slab_end+1)*t(ids-1,k_slab_end,j) 
                      
                      z_slab_top = fnm(k_slab_end+1)*z(ids-1,k_slab_end+1,j) + fnp(k_slab_end+1)*z(ids-1,k_slab_end,j) 
                      
                      duslabpsum = duslabpsum + u_slab_top - u_slab_bot
                      
                      dvslabpsum = dvslabpsum + v_slab_top - v_slab_bot
                                            
                      dtslabpsum = dtslabpsum + t_slab_top - t_slab_bot
                      
                      dzslabpsum = dzslabpsum + z_slab_top - z_slab_bot
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                      DO k = k_slab_start, k_slab_end   ! Also need T in slab. If slab > 1 grid cell deep, the sum will have more points than the differences.  
                                                        ! Therefore use a separate counter variable to compute the average.
                         
                         tslabpsum = tslabpsum + t(ids-1,k,j)
                         
                         num_pts_in_sum_2 = num_pts_in_sum_2 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( ite .EQ. ide ) THEN
                   
                   DO j = j_start,j_end
                      
                      IF ( k_slab_start .EQ. 1 ) THEN 
                         
                         u_slab_bot = u(ide+1,1,j) - ( ( u(ide+1,2,j) - u(ide+1,1,j)  )/( z(ide+1,2,j) - z(ide+1,1,j)) )*z(ide+1,1,j)
                         
                         v_slab_bot = v(ide,1,j) - ( ( v(ide,2,j) - v(ide,1,j)  )/( z(ide,2,j) - z(ide,1,j)) )*z(ide,1,j)
                         
                         t_slab_bot = t(ide,1,j) - ( ( t(ide,2,j) - t(ide,1,j)  )/( z(ide,2,j) - z(ide,1,j)) )*z(ide,1,j)
                         
                         z_slab_bot = 0.0
                         
                      ELSE
                         
                         u_slab_bot = fnm(k_slab_start)*u(ide+1,k_slab_start,j) + fnp(k_slab_start)*u(ide+1,k_slab_start-1,j)
                         
                         v_slab_bot = fnm(k_slab_start)*v(ide,k_slab_start,j) + fnp(k_slab_start)*v(ide,k_slab_start-1,j)
                         
                         t_slab_bot = fnm(k_slab_start)*t(ide,k_slab_start,j) + fnp(k_slab_start)*t(ide,k_slab_start-1,j)
                         
                         z_slab_bot = fnm(k_slab_start)*z(ide,k_slab_start,j) + fnp(k_slab_start)*z(ide,k_slab_start-1,j) 
                         
                      ENDIF
                      
                      u_slab_top = fnm(k_slab_end+1)*u(ide+1,k_slab_end+1,j) + fnp(k_slab_end+1)*u(ide+1,k_slab_end,j) 
                      
                      v_slab_top = fnm(k_slab_end+1)*v(ide,k_slab_end+1,j) + fnp(k_slab_end+1)*v(ide,k_slab_end,j) 
                      
                      t_slab_top = fnm(k_slab_end+1)*t(ide,k_slab_end+1,j) + fnp(k_slab_end+1)*t(ide,k_slab_end,j) 
                      
                      z_slab_top = fnm(k_slab_end+1)*z(ide,k_slab_end+1,j) + fnp(k_slab_end+1)*z(ide,k_slab_end,j)
                      
                      duslabpsum = duslabpsum + u_slab_top - u_slab_bot
                      
                      dvslabpsum = dvslabpsum + v_slab_top - v_slab_bot
                      
                      dtslabpsum = dtslabpsum + t_slab_top - t_slab_bot
                      
                      dzslabpsum = dzslabpsum + z_slab_top - z_slab_bot
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                      DO k = k_slab_start, k_slab_end   
                         
                         tslabpsum = tslabpsum + t(ide,k,j)
                         
                         num_pts_in_sum_2 = num_pts_in_sum_2 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jts .EQ. jds ) THEN
                   
                   DO i = i_start,i_end
                      
                      IF ( k_slab_start .EQ. 1 ) THEN 
                         
                         u_slab_bot = u(i,1,jds-1) - ( ( u(i,2,jds-1) - u(i,1,jds-1)  )/( z(i,2,jds-1) - z(i,1,jds-1)) )*z(i,1,jds-1)
                         
                         v_slab_bot = v(i,1,jds-1) - ( ( v(i,2,jds-1) - v(i,1,jds-1)  )/( z(i,2,jds-1) - z(i,1,jds-1)) )*z(i,1,jds-1)
                         
                         t_slab_bot = t(i,1,jds-1) - ( ( t(i,2,jds-1) - t(i,1,jds-1)  )/( z(i,2,jds-1) - z(i,1,jds-1)) )*z(i,1,jds-1)
                         
                         z_slab_bot = 0.0
                         
                      ELSE
                         
                         u_slab_bot = fnm(k_slab_start)*u(i,k_slab_start,jds-1) + fnp(k_slab_start)*u(i,k_slab_start-1,jds-1)
                         
                         v_slab_bot = fnm(k_slab_start)*v(i,k_slab_start,jds-1) + fnp(k_slab_start)*v(i,k_slab_start-1,jds-1)
                         
                         t_slab_bot = fnm(k_slab_start)*t(i,k_slab_start,jds-1) + fnp(k_slab_start)*t(i,k_slab_start-1,jds-1)
                         
                         z_slab_bot = fnm(k_slab_start)*z(i,k_slab_start,jds-1) + fnp(k_slab_start)*z(i,k_slab_start-1,jds-1)                           
                         
                      ENDIF
                      
                      u_slab_top = fnm(k_slab_end+1)*u(i,k_slab_end+1,jds-1) + fnp(k_slab_end+1)*u(i,k_slab_end,jds-1) 
                      
                      v_slab_top = fnm(k_slab_end+1)*v(i,k_slab_end+1,jds-1) + fnp(k_slab_end+1)*v(i,k_slab_end,jds-1) 
                      
                      t_slab_top = fnm(k_slab_end+1)*t(i,k_slab_end+1,jds-1) + fnp(k_slab_end+1)*t(i,k_slab_end,jds-1) 
                      
                      z_slab_top = fnm(k_slab_end+1)*z(i,k_slab_end+1,jds-1) + fnp(k_slab_end+1)*z(i,k_slab_end,jds-1)
                      
                      duslabpsum = duslabpsum + u_slab_top - u_slab_bot
                      
                      dvslabpsum = dvslabpsum + v_slab_top - v_slab_bot
                                            
                      dtslabpsum = dtslabpsum + t_slab_top - t_slab_bot
                      
                      dzslabpsum = dzslabpsum + z_slab_top - z_slab_bot
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                      DO k = k_slab_start, k_slab_end 
                         
                         tslabpsum = tslabpsum + t(i,k,jds-1)
                         
                         num_pts_in_sum_2 = num_pts_in_sum_2 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jte .EQ. jde ) THEN
                   
                   DO i = i_start,i_end
                      
                      IF ( k_slab_start .EQ. 1 ) THEN 
                         
                         u_slab_bot = u(i,1,jde) - ( ( u(i,2,jde) - u(i,1,jde)  )/( z(i,2,jde) - z(i,1,jde)) )*z(i,1,jde)
                         
                         v_slab_bot = v(i,1,jde+1) - ( ( v(i,2,jde+1) - v(i,1,jde+1)  )/( z(i,2,jde+1) - z(i,1,jde+1)) )*z(i,1,jde+1)
                         
                         t_slab_bot = t(i,1,jde) - ( ( t(i,2,jde) - t(i,1,jde)  )/( z(i,2,jde) - z(i,1,jde)) )*z(i,1,jde)
                         
                         z_slab_bot = 0.0
                         
                      ELSE
                         
                         u_slab_bot = fnm(k_slab_start)*u(i,k_slab_start,jde) + fnp(k_slab_start)*u(i,k_slab_start-1,jde)
                         
                         v_slab_bot = fnm(k_slab_start)*v(i,k_slab_start,jde+1) + fnp(k_slab_start)*v(i,k_slab_start-1,jde+1)
                         
                         t_slab_bot = fnm(k_slab_start)*t(i,k_slab_start,jde) + fnp(k_slab_start)*t(i,k_slab_start-1,jde)
                         
                         z_slab_bot = fnm(k_slab_start)*z(i,k_slab_start,jde) + fnp(k_slab_start)*z(i,k_slab_start-1,jde) 
                         
                      ENDIF
                      
                      u_slab_top = fnm(k_slab_end+1)*u(i,k_slab_end+1,jde) + fnp(k_slab_end+1)*u(i,k_slab_end,jde) 
                      
                      v_slab_top = fnm(k_slab_end+1)*v(i,k_slab_end+1,jde+1) + fnp(k_slab_end+1)*v(i,k_slab_end,jde+1) 
                      
                      t_slab_top = fnm(k_slab_end+1)*t(i,k_slab_end+1,jde) + fnp(k_slab_end+1)*t(i,k_slab_end,jde) 
                      
                      z_slab_top = fnm(k_slab_end+1)*z(i,k_slab_end+1,jde) + fnp(k_slab_end+1)*z(i,k_slab_end,jde)
                      
                      duslabpsum = duslabpsum + u_slab_top - u_slab_bot
                      
                      dvslabpsum = dvslabpsum + v_slab_top - v_slab_bot
                                            
                      dtslabpsum = dtslabpsum + t_slab_top - t_slab_bot
                      
                      dzslabpsum = dzslabpsum + z_slab_top - z_slab_bot
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                      DO k = k_slab_start, k_slab_end 
                         
                         tslabpsum = tslabpsum + t(i,k,jde)
                         
                         num_pts_in_sum_2 = num_pts_in_sum_2 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                tslabdsum = tslabpsum
                   
                duslabdsum  = duslabpsum
                
                dvslabdsum  = dvslabpsum
                 
                dtslabdsum  = dtslabpsum
                
                dzslabdsum  = dzslabpsum
                
                counter_sum_1 = num_pts_in_sum_1
                
                counter_sum_2 = num_pts_in_sum_2
                
#ifdef DM_PARALLEL
                
                tslabdsum = wrf_dm_sum_real8(tslabpsum)
                
                duslabdsum  = wrf_dm_sum_real8(duslabpsum)
                
                dvslabdsum  = wrf_dm_sum_real8(dvslabpsum)
                
                dtslabdsum  = wrf_dm_sum_real8(dtslabpsum)
                
                dzslabdsum  = wrf_dm_sum_real8(dzslabpsum)
                
                counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
                
                counter_sum_2 = wrf_dm_sum_real(num_pts_in_sum_2)
                
#endif
                
                duslabdavg = duslabdsum/counter_sum_1
                
                dvslabdavg = dvslabdsum/counter_sum_1
                
                dtslabdavg = dtslabdsum/counter_sum_1
                
                dzslabdavg = dzslabdsum/counter_sum_1
                
                tslabdavg = tslabdsum/counter_sum_2 + 300.0 !convert to real K
                
                pertmag = dtslabdavg - cpm_rim*tslabdavg*( duslabdavg*duslabdavg + dvslabdavg*dvslabdavg )/(g*dzslabdavg)
                
                
             ENDIF !cpm_opt = 3: Richardson number (stable) formulation
             
             
             IF ( cpm_opt .EQ. 4 ) THEN ! Unstable formulation
                
                ! Need average heat flux for scaling. Ghost point values do not exist, so use first interior
                ! grid point values since those will be very close.
                
                hfxpsum = 0.0  
                num_pts_in_sum_1 = 0.0   
                
                IF ( its .EQ. ids ) THEN
                   
                   DO j = j_start,j_end

                      IF ( cpm_meso_hfx .EQ. 1 ) THEN 
                         
                         hfxpsum = hfxpsum + meso_hfx(ids,j)
                         
                      ELSE
                         
                         hfxpsum = hfxpsum + spec_hfx/(rho(ids,kts,j)*cp)
                         
                      END IF
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( ite .EQ. ide ) THEN
                   
                   DO j = j_start,j_end

                      IF ( cpm_meso_hfx .EQ. 1 ) THEN 
                         
                         hfxpsum = hfxpsum + meso_hfx(ide-1,j)
                         
                      ELSE

                         hfxpsum = hfxpsum + spec_hfx/(rho(ide-1,kts,j)*cp)

                      END IF
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( jts .EQ. jds ) THEN
                   
                   DO i = i_start,i_end

                      IF ( cpm_meso_hfx .EQ. 1 ) THEN 
                         
                         hfxpsum = hfxpsum + meso_hfx(i,jds) 
                         
                      ELSE

                         hfxpsum = hfxpsum + spec_hfx/(rho(i,kts,jds)*cp)
                         
                      END IF
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                IF ( jte .EQ. jde ) THEN
                   
                   DO i = i_start,i_end

                      IF ( cpm_meso_hfx .EQ. 1 ) THEN 
                      
                         hfxpsum = hfxpsum + meso_hfx(i,jde-1) 
                      
                      ELSE

                         hfxpsum = hfxpsum + spec_hfx/(rho(i,kts,jde-1)*cp)

                      END IF
                      
                      num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                      
                   END DO
                   
                ENDIF
                
                hfxdsum = hfxpsum
                
                counter_sum_1 = num_pts_in_sum_1
                
#ifdef DM_PARALLEL
                
                hfxdsum  = wrf_dm_sum_real8(hfxpsum)
                
                counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
                
#endif
                
                hfxdavg = hfxdsum/counter_sum_1
                
                
                tslabpsum = 0.0
                num_pts_in_sum_1 = 0.0
                
                IF ( its .EQ. ids ) THEN
                   
                   DO j = j_start,j_end
                      
                      DO k = k_slab_start, k_slab_end
                         
                         tslabpsum = tslabpsum + t(ids-1,k,j)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( ite .EQ. ide ) THEN
                   
                   DO j = j_start,j_end
                         
                      DO k = k_slab_start, k_slab_end   
                         
                         tslabpsum = tslabpsum + t(ide,k,j)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jts .EQ. jds ) THEN
                   
                   DO i = i_start,i_end
                      
                      DO k = k_slab_start, k_slab_end 
                         
                         tslabpsum = tslabpsum + t(i,k,jds-1)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                IF ( jte .EQ. jde ) THEN
                   
                   DO i = i_start,i_end
                      
                      DO k = k_slab_start, k_slab_end 
                         
                         tslabpsum = tslabpsum + t(i,k,jde)
                         
                         num_pts_in_sum_1 = num_pts_in_sum_1 + 1.0
                         
                      END DO
                      
                   END DO
                   
                ENDIF
                
                tslabdsum = tslabpsum
                
                counter_sum_1 = num_pts_in_sum_1
                
#ifdef DM_PARALLEL
                
                tslabdsum = wrf_dm_sum_real8(tslabpsum)
                
                counter_sum_1 = wrf_dm_sum_real(num_pts_in_sum_1)
                
#endif
                
                tslabdavg = tslabdsum/counter_sum_1 + 300.0 !convert to real K
                
                zi = MAX(cpm_lim_z,pblhdavg)
                
                thetastar = (hfxdavg**0.6667)/( ( (g/tslabdavg)*zi )**0.3333 )
                
                f1 = (zslabdavg/zi)**(-0.6667)
                f2 = (1.0 - zslabdavg/zi)**(1.3333)
                f3 = (zslabdavg/zi)**(1.3333)
                f4 = (1.0 - zslabdavg/zi)**(-0.6667)
                
                pertmag = thetastar*( 2.0*f1*f2 + 0.94*f3*f4 )
                   
             ENDIF ! Unstable
             
             ! Calculate (new) perturbations  ==========================
                          
             ni = (ide-1)/cpm_ngc_h + 1 
             nj = (jde-1)/cpm_ngc_h + 1
             
             ! Allocate perturbation cell array for each lateral boundary 

             IF (west .EQ. 1) ALLOCATE( perts_wb(1:nj,1:cpm_ncells_h) )
             IF (east .EQ. 1) ALLOCATE( perts_eb(1:nj,1:cpm_ncells_h) )
             IF (south .EQ. 1) ALLOCATE( perts_sb(1:ni,1:cpm_ncells_h) )
             IF (north .EQ. 1) ALLOCATE( perts_nb(1:ni,1:cpm_ncells_h) )
             
             CALL RANDOM_SEED(size=i_seed)       ! Get the size of the seed, each processor.
             
             ALLOCATE( seed(1:i_seed) )          ! Allocate an array to hold the seed, each processor.
             
#ifdef DM_PARALLEL
             
             IF ( mytask .EQ. master ) THEN
                
                seedsum = 0                      ! Sum the stored seed (up to i_seed) to see if it has been assigned.
                DO k = 1,i_seed                  ! If so, now master has the current prtseed value in seed.
                   seed(k) = prtseed(k)
                   seedsum = seedsum + seed(k)
                END DO
                
                IF (seedsum .EQ. 0 ) THEN        ! If the seed is empty, get the first seed, which is now only on master.
                    
                   CALL RANDOM_SEED(get=seed)
                   
                ENDIF
                
             ENDIF
             
             CALL MPI_BCAST(seed,i_seed,MPI_REAL,master,MPI_COMM_WORLD,ierr) ! Broadcast seed from master to all
             
             CALL RANDOM_SEED(put=seed) !Use same seed value so all processors receive idential random numbers
             
             IF (west .EQ. 1) CALL RANDOM_NUMBER(perts_wb)
             IF (east .EQ. 1) CALL RANDOM_NUMBER(perts_eb)
             IF (south .EQ. 1) CALL RANDOM_NUMBER(perts_sb)
             IF (north .EQ. 1) CALL RANDOM_NUMBER(perts_nb)
             
             IF ( mytask .EQ. master ) THEN ! Must get a new seed and store into prtseed so the sequence at the time of the
                                            ! next perturbation update begins with new seed. 
                
                CALL RANDOM_SEED(get=seed)
                
                DO k = 1,i_seed
                   prtseed(k) = seed(k)
                ENDDO!
                
             ENDIF
             
#endif

             ! Apply perturbations along chosen lateral boundaries =====
             
             IF ( west .EQ. 1 ) THEN
                
                IF (its .LE. ids + cpm_ngc_h*cpm_ncells_h) THEN 
                   
                   DO j = MAX(jts,jds + (1-south)*cpm_noff_h), MIN(jte,jde - 1 - (1-north)*cpm_noff_h)
                      
                      DO i = its, MIN(ite, ide-1)
                         
                         m = (i-1)/cpm_ngc_h+1
                         
                         IF ( m .LE. cpm_ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end
                               
                               cpm_3d(i,k,j) = (perts_wb(((j-1)/cpm_ngc_h+1),m)-0.5)*2.0*pertmag !Output variable only, for debugging.
                               
                               IF (cpm_opt .LE. 4) THEN 
                                  
                                  t(i,k,j) = t(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                               IF (cpm_opt .EQ. 5) THEN  
                                  
                                  w(i,k,j) = w(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO! j
                   
                ENDIF
                
             ENDIF
             
             IF ( east .EQ. 1 ) THEN
                
                IF (ite .GE. ide - 1 - cpm_ngc_h*cpm_ncells_h) THEN 
                   
                   DO j = MAX(jts,jds + (1-south)*cpm_noff_h), MIN(jte,jde - 1 - (1-north)*cpm_noff_h)
                      
                      DO i = MIN(ite, ide-1), its, -1
                         
                         m = (ide-i-1)/cpm_ngc_h+1
                         
                         IF ( m .LE. cpm_ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end
                               
                               cpm_3d(i,k,j) = (perts_eb(((j-1)/cpm_ngc_h+1),m)-0.5)*2.0*pertmag
                               
                               IF (cpm_opt .LE. 4) THEN
                                  
                                  t(i,k,j) = t(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                               IF (cpm_opt .EQ. 5) THEN 
                                  
                                  w(i,k,j) = w(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             IF ( south .EQ. 1 ) THEN
                
                IF (jts .LE. jds + cpm_ngc_h*cpm_ncells_h) THEN 
                   
                   DO i = MAX(its,ids + cpm_noff_h*(1 - west)), MIN(ite,ide - 1 - cpm_noff_h*(1 - east))
                      
                      DO j = jts, MIN(jte, jde-1)
                         
                         m = (j-1)/cpm_ngc_h+1
                         
                         IF ( m .LE. cpm_ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end
                               
                               cpm_3d(i,k,j) = (perts_sb(((i-1)/cpm_ngc_h+1),m)-0.5)*2.0*pertmag
                               
                               IF (cpm_opt .LE. 4) THEN
                                  
                                  t(i,k,j) = t(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                               IF (cpm_opt .EQ. 5) THEN 
                                  
                                  w(i,k,j) = w(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             IF ( north .EQ. 1 ) THEN
                
                IF (jte .GE. jde - 1 - cpm_ngc_h*cpm_ncells_h) THEN 
                   
                   DO i = MAX(its,ids + cpm_noff_h*(1 - west)), MIN(ite,ide - 1 - cpm_noff_h*(1 - east))
                      
                      DO j = MIN(jte, jde-1),jts,-1
                         
                         m = (jde-j-1)/cpm_ngc_h+1
                         
                         IF ( m .LE. cpm_ncells_h ) THEN
                            
                            DO k = k_slab_start, k_slab_end
                               
                               cpm_3d(i,k,j) = (perts_nb(((i-1)/cpm_ngc_h+1),m)-0.5)*2.0*pertmag
                               
                               IF (cpm_opt .LE. 4) THEN
                                  
                                  t(i,k,j) = t(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                               IF (cpm_opt .EQ. 5) THEN 
                                  
                                  w(i,k,j) = w(i,k,j) + cpm_3d(i,k,j)
                                  
                               ENDIF
                               
                            END DO ! k
                            
                         ENDIF
                         
                      END DO ! i
                      
                   END DO ! j
                   
                ENDIF
                
             ENDIF
             
             DEALLOCATE(seed)
             
             IF (west .EQ. 1) DEALLOCATE( perts_wb)
             IF (east .EQ. 1) DEALLOCATE( perts_eb )
             IF (south .EQ. 1) DEALLOCATE( perts_sb )
             IF (north .EQ. 1) DEALLOCATE( perts_nb )
             
200          CONTINUE
             
             
          END IF ! IF ( min_slab_z .GE. max_pblh ) 
          
       END IF !(  prttms(slab_k) .GE. prtdt(slab_k) ) THEN
   
    END DO !Big outer k-loop over number of vertical perturbation slabs

    
    !===================================================================
    
    
  END SUBROUTINE calc_cpm
  
  !=====================================================================
  
  
END MODULE module_les_cpm

!=======================================================================
!=======================================================================
