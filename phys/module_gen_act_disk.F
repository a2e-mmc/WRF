
MODULE module_gen_act_disk

!RSA for converting lat/lon to x/y
USE module_llxy

!RSA for turbine-specific functions
USE module_gad_psu15, ONLY: calc_rpm_psu15, calc_pitch_psu15, calc_fnft_psu15
USE module_gad_nrel5, ONLY: calc_rpm_nrel5, calc_pitch_nrel5, calc_fnft_nrel5

#ifdef DM_PARALLEL
  USE MODULE_DM
#endif
  
  CONTAINS
    
  SUBROUTINE gen_act_disk ( ru_tendf, rv_tendf, rw_tendf, &
                            wp_opt,                       & 
                            turbine_opt,                  & !RSA
                            windturb_spec,                & !RSA
                            turb_loc_in_latlon,           & !RSA
                            truelat1, truelat2,           & !RSA
                            stand_lon,                    & !RSA
                            n_turbines,                   &
                            n_timeseries,                 &
                            x_turbine,                    &
                            y_turbine,                    &
                            x_inflow,                     & !RSA
                            y_inflow,                     & !RSA
                            hub_height,                   &
                            rotor_diameter,               &
                            blade_length,                 &
                            theta_turbine,                &
                            cutinspd, cutoutspd,          & !RSA
                            turb_power,                   & !RSA
                            acc_yaw_err,                  & 
                            wp_acc_u, wp_acc_v, wp_acc_w, & ! these variables likely not needed. 
                            wp_ts,                        &
                            thrust, torque,               &
                            rotrate, power,               &
                            u, v, w,                      &
                            muu, muv, mut,                &
                            rdzw,  dx, dy, dt,            &
                            v0t, d0t, itimestep,          &
                            xlat_v, xlong_u, hgt,         & !RSA
                            fnm, fnp,                     & !RSA
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            ips, ipe, jps, jpe, kps, kpe, &
                            its, ite, jts, jte, kts, kte  )

!-----------------------------------------------------------------------------
!
! Purpose: To calculate forces using the generalized actuator disk model. 
!
! See Mirocha et al. (2014), J. Renew. Sust. Energy:
! "Implementation of a generalized actuator disk wind turbine model into the
! weather research and forecasting model for large-eddy simulation applications"
!
! Lookup tables use generic Penn State 1.5 MW GE SLE turbine (approximate).
!
! Note that wind direction and turbine yaw are measured counterclockwise from
! 0 deg. for a westerly wind / west-facing turbine. Thus, 90 is a southerly
! wind / south-facing turbine, 180 is an easterly wind / east-facing turbine,
! 270 is a northerly wind / north-facing turbine, up to 360 deg.
!
! Last modified 10/2018 by RSA and JDM. 
! - Code will not compute GAD forces if flow is from behind the turbine; turbine
! must first yaw into mean flow direction (when yaw error exceeds threshold) 
! before forces will be calculated.
! - Added ability specify upwind distance for turbine inflow calculation relative 
! to turbine hub.
! - Added ability to read in turbine info from a file specifying turbine
! parameters, including position in either x,y or lat,lon.
!
!=============================================================================

  IMPLICIT NONE

#if defined(DM_PARALLEL) && !defined(STUBMPI)
  INCLUDE "mpif.h" ! MSG TESTING
#endif

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(INOUT) &
  :: ru_tendf    & ! u tendency * mass                   (kg s-4)
   , rv_tendf    & ! v
   , rw_tendf      ! w
  
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(INOUT) &
  :: wp_acc_u    & ! acceleration for u                  (m s-2)
   , wp_acc_v    & !                  v                  (m s-2)
   , wp_acc_w      !                  w                  (m s-2)
 
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT(INOUT) & 
  :: wp_ts         ! diagnostic for debuging

  INTEGER,                                       INTENT( IN  ) &  
  :: n_turbines    ! number of turbines

  INTEGER,                                       INTENT( IN  ) &  
  :: n_timeseries  ! length of timeseries for averages

  INTEGER,                                       INTENT( IN  ) &  
  :: wp_opt, turbine_opt

  CHARACTER*256,                                 INTENT( IN  ) &
  :: windturb_spec !RSA file for reading in turbine info

  LOGICAL,                                       INTENT( IN  ) &
  :: turb_loc_in_latlon !RSA flag for turbine location in lat/lon, else in meters

  REAL,                                          INTENT( IN  ) &
  :: truelat1, truelat2, stand_lon !RSA for converting lat/lon to x/y

  REAL, DIMENSION( ims:ime, jms:jme ),           INTENT( IN  ) &
  :: xlat_v, xlong_u !RSA for converting lat/lon to x/y

  REAL, DIMENSION( n_turbines ),                 INTENT(INOUT) &  !should be in if reading from aux file
  :: x_turbine      & ! x location of a turbine              (m)
   , y_turbine      & ! y location of a turbine              (m)
   , x_inflow       & ! x location of turbine inflow calcs   (m)
   , y_inflow       & ! y location of turbine inflow calcs   (m)
   , hub_height     & ! hub height                           (m)
   , rotor_diameter & ! hub height                           (m)
   , blade_length   & ! blade_length                         (m)
   , theta_turbine  & ! angle of rotor plane w.r.t. x axis (deg)
   , cutinspd       & ! turbine cutin speed                (m/s)
   , cutoutspd      & ! turbine cutout speed               (m/s)
   , turb_power       ! turbine power output          (kg*m2/s3)

  REAL, DIMENSION( n_turbines ),                 INTENT(INOUT) & 
  ::  acc_yaw_err

  REAL, DIMENSION( ims:ime, jms:jme ),            INTENT( OUT ) &
  :: thrust, torque, rotrate, power

   REAL, DIMENSION( n_timeseries, n_turbines ),   INTENT(INOUT) &
  :: v0t          & ! Wind speed timeseries               (m s-1)
   , d0t            ! Wind direction timeseries.            (deg)
                    ! 0 is north, positive rotating clockwise  

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT( IN  ) &
  :: u           & ! u-velocity                  (m s-1)
   , v           & ! v
   , w             ! w

  REAL, DIMENSION( ims:ime, jms:jme ),           INTENT( IN  ) &
  :: muu         & ! column pressure at u nodes  (kg m-1 s-2)
   , muv         & !                    v
   , mut         & !                    t
   , hgt           ! RSA terrain height          (m)

  REAL, DIMENSION( kms:kme ),                    INTENT( IN  ) &
  :: fnm, fnp      ! RSA for v0_w calculation

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT( IN  ) &
  :: rdzw          !inverse vertical grid spacing at w-nodes (m-1)

  REAL,                                          INTENT( IN  ) &
  :: dx          & !grid spacing in x-direction  (m) 
   , dy            !                y            (m)

  REAL,                                          INTENT( IN  ) & !JDM May 2, 2013
  :: dt

  INTEGER,                                       INTENT( IN  ) & !JDM May 2, 2013
  :: itimestep

  INTEGER,                                       INTENT( IN  ) &
  :: ids, ide, jds, jde, kds, kde   &
   , ims, ime, jms, jme, kms, kme   &
   , ips, ipe, jps, jpe, kps, kpe   &
   , its, ite, jts, jte, kts, kte

! LOCAL VARIABLES

  REAL, DIMENSION( its:ite+1 ) &
  :: u_x         & !x coordinate of u-nodes (m)
   , v_x         & !                v       (m)
   , w_x           !                w       (m)
  
  REAL, DIMENSION( jts:jte+1 ) &
  :: u_y         & !y coordinate of u-nodes (m)
   , v_y         & !                v       (m)
   , w_y           !                w       (m)

  REAL, DIMENSION( its:ite+1, kts:kte, jts:jte+1 ) &
  :: u_z         & !z coordinate of u-nodes (m)
   , v_z           !                v       (m)
   
  REAL, DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1 ) &
  :: w_z         & !z coordinate of w-nodes  (m)
   , c_z           !z coordinate of c-nodes  (m) 

  REAL, DIMENSION( its:ite, jts:jte ) &
  :: dru       & !normal distance from u-nodes to the rotor plane (m)
   , drv       & !                     v                          (m)
   , drw         !                     w                          (m)

  REAL, DIMENSION( its:ite, jts:jte ) &
  :: dau       & !normal distance from u-nodes to the rotor axis line (m)
   , dav       & !                     v                              (m)
   , daw         !                     w                              (m)

  REAL, DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1 ) &
  :: dzu       & !normal distance from u-nodes to the rotor plane (m)
   , dzv       & !                    v                           (m)
   , dzw       & !                    w                           (m)
   , dzc         !                    c                           (m)

  REAL, DIMENSION( its-1:ite+1, kts:kte, jts-1:jte+1 ) &
  :: ru       & !distance from u-nodes to the nacelle            (m)
   , rv       & !                    v                           (m)
   , rw       & !                    w                           (m)
   , rc         !                    c                           (m)

  REAL, DIMENSION( its:ite, kts:kte, jts:jte ) &
  :: dz !RSA thickness of each grid cell                         (m)

  REAL, DIMENSION( kts:kte ) &
  :: turb_z   & ! z coordinate at a turbine                      (m)  
   , inflow_z   !RSA z coordinate at inflow location             (m)

  REAL :: xmin_tile, xmax_tile, ymin_tile, ymax_tile
  REAL :: xmin_turb, xmax_turb, ymin_turb, ymax_turb
  REAL :: zmax_turb
  REAL :: start_fac

  INTEGER :: i, j, k, n, i_start, j_start, i_end, j_end, k_end, k_top, i_inflow, j_inflow, k_inflow, &
             i_turbine, j_turbine

  REAL :: dxo2, dyo2, zmin_k, zmin_kp1, zmax_k, zmax_kp1, eps

  REAL :: theta

  INTEGER :: naninf

#ifdef DM_PARALLEL
  INTEGER         :: ierr
  INTEGER         :: tag
  INTEGER         :: master
  INTEGER         :: status(MPI_STATUS_SIZE)
#endif

  REAL, DIMENSION(n_turbines) :: v0n
  REAL, DIMENSION(n_turbines) :: d0n
  INTEGER, DIMENSION(n_turbines) :: yawing


  REAL, DIMENSION(2) :: v0d0inst

!--------------------------------------------------------

  REAL, DIMENSION(2)    :: fnft
  REAL                  :: fn, ft
  REAL                  :: v0, v0_u, v0_v, v0_w, v0_c, d0, v0t_avg, d0t_avg !RSA
  REAL                  :: u_c, v_c, w_c !RSA
  REAL                  :: rho
  REAL                  :: rpm, blades, pitchdeg
  REAL                  :: rotor_radius
  REAL                  :: dscale, cscale
  REAL                  :: zeta
  REAL                  :: fx, fy, fz !RSA
  REAL                  :: cnst
  REAL                  :: v0inst
  INTEGER               :: nv0t

  INTEGER, DIMENSION(n_turbines) :: turbine_on_patch
  INTEGER, DIMENSION(n_turbines) :: inflow_on_patch
  REAL, DIMENSION(n_turbines) :: turb_power_local !turbine power calculation on local processor
  REAL, DIMENSION(n_turbines) :: omega !turbine rotation rate (rad/s)

  REAL                  :: time_fac
  REAL                  :: xtime
  REAL                  :: u_geo
  REAL                  :: v_geo
  REAL                  :: th_geo
  REAL                  :: ws_geo
  REAL                  :: nt

  REAL                  :: acc_u, acc_v, acc_w, acc_c
  REAL                  :: radmin, radmax

! variables for yawing

  REAL                  :: d0inst
  REAL                  :: yaw_thresh
  REAL                  :: yaw_err
  REAL                  :: yaw_inc
  REAL                  :: yaw_new
  REAL                  :: yaw_rate
  REAL                  :: delta
  REAL                  :: sgn
  REAL                  :: dum
  REAL                  :: operate
! variables for hub drag

  REAL                  :: cdh

! variables for tower drag

  REAL                  :: cdt, d_tower, Rtwrt, Rtwrb, Rtwr

! variables for tilt

  REAL                  :: tilt_deg 
  REAL                  :: tilt   
  REAL                  :: sigma
  REAL                  :: sigmar
  REAL                  :: adenom, bdenom
  REAL                  :: adenom2, bdenom2
  REAL                  :: hoff
  REAL                  :: drut, drvt, drwt, drct
  REAL                  :: radius

! RSA for reading in turbine info from file
  CHARACTER*512        :: inline
  INTEGER              :: istat
  REAL                 :: lat,lon,hubheight,diameter,bldlngth,thetaturb,cutin,cutout

! RSA for converting lat/lon to x/y
  TYPE(proj_info)      :: proj
  REAL                 :: lat1,lon1,knowni,knownj,inew,jnew!,truelat1,truelat2,stand_lon

! test
  REAL :: zd, zetazd
! temporary variable from MSG debugging
REAL :: temp
!---------------------------------------------------------

  CHARACTER(LEN=12), PARAMETER :: FMT1 = "(I5,1F18.12)"
  CHARACTER(LEN=12), PARAMETER :: FMT2 = "(I5,2F18.12)"
  CHARACTER(LEN=12), PARAMETER :: FMT3 = "(I5,3F18.12)"
  CHARACTER(LEN=12), PARAMETER :: FMT4 = "(3I5,2F18.12)"

  REAL, DIMENSION(10) :: turb_xy_locs

#ifdef DM_PARALLEL
  master = 0
  tag    = 0
#endif

!____________________________________________________________________________________________________________________________
!
! 0. Set turbine parameters (to be assigned via auxiliary input file in future version)
!____________________________________________________________________________________________________________________________

! PRINT *,'Starting GAD routine'

  hoff = 5.0 ! hub offset
  yaw_rate = 2.0 !degrees per second
  pitchdeg=  0.0
  dscale  =  1.0
  cscale  =  1.0
  rpm     = 16.0
  blades  =  3
  rho     =  1.225
  start_fac = 1.0
  yaw_thresh = 10000.0 !10.0*60.0*10.0**2
  ! yaw_thresh = 100.0 !10.0*60.0*10.0**2 !RSAtest
  yaw_new = 0.0
  tilt_deg = 4.0 ! in degrees
  tilt = tilt_deg*pi/180.0
  cdh = 0.0! 1.0
  cdt = 0.0!  1.2
  !JDM-04-2015  d_tower = 2.0 !m
  Rtwrt = 1.5 !NM radius is say 1.5 meters
  Rtwrb = 2.0 !NM radius of tower at base

  !RSA if this is the first time the GAD is turned on, then read in turbine info 
  !from file or hard coded values
  dum = 0.0
  DO k = 1,n_timeseries
      dum = dum + v0t(k,1)
  ENDDO

if (dum .LT. 0.0001) then
   print*,'use dum variable to read files'
else
   ! print*,'dum variable is',dum
endif

! print*,'dum variable is equal to',dum,'after loop'
! print*,'check if i am reading the turbine files every time',' maybe use x_turbine variable (if equal to zero, then continue)',x_turbine(1)
 
  IF ( x_turbine(1) .EQ. 0.0 ) THEN !MS changed this from dum .EQ. 0.0 
     !RSA for turbine info read from file (adapted from NM code)
     IF (TRIM(windturb_spec) .NE. 'none') THEN
        OPEN(file=TRIM(windturb_spec),unit=19,FORM='FORMATTED',STATUS='OLD',IOSTAT=istat)
        IF ( istat .EQ. 0 ) THEN !no error in file open
           n = 0
           !count number of lines in file
           DO WHILE (.true.)
              READ(19,'(A256)',END=30) inline !RSA read file and goto 30 if end of file
              IF ( index(inline,'!') .EQ. 0 ) n = n + 1
              print*,'file reading',n,'turbines'
           ENDDO

30         CONTINUE !RSA reached end of file
           IF (n_turbines .NE. n) THEN
              print*,'RSA n_turbines in namelist does not match number of turbines in file windturb_spec'
              STOP
              ! CALL wrf_error_fatal('RSA n_turbines in namelist does not match number of turbines in file windturb_spec')
           ENDIF

           REWIND(19) !RSA go back to beginning of file
           print*,'RSA reading input file windturb_spec for ',n,' turbines'
           i = 1
           DO WHILE (.true.)
              READ(19,'(A256)',END=120) inline !RSA read file and goto 120 if end of file
              IF ( i .LE. n_turbines .AND. index(inline,'!') .EQ. 0 ) THEN
                 !RSA read in each part of the line seperated by spaces and go to line 130 if there is an error
                 READ(inline,*,ERR=130) lat,lon,hubheight,diameter,bldlngth,thetaturb,cutin,cutout 
                 x_turbine(i) = lon   
                 y_turbine(i) = lat
                 hub_height(i) = hubheight
                 rotor_diameter(i) = diameter
                 blade_length(i) = bldlngth
                 theta_turbine(i) = thetaturb
                 cutinspd(i) = cutin
                 cutoutspd(i) = cutout
                 i = i + 1 !RSA counting lines
                ! PRINT *,'Reading file turbines(MSG)',x_turbine(i-1),y_turbine(i-1),'locs',hub_height(i-1),'hub height',rotor_diameter(i-1), 'diameter',blade_length(i-1),'blades',theta_turbine(i-1),'angle',cutinspd(i-1),'cutin',cutoutspd(i-1),'cutout'
              ENDIF
           ENDDO

120	     CONTINUE !RSA reached end of file
           CLOSE(19)
           GOTO 150

130	     CONTINUE !RSA error so close file and set istat to nonzero
           CLOSE(19)
           istat = 150150
           GOTO 150

150        CONTINUE
           print*,'RSA finished reading input file'
        ELSE
           print*,'RSA UNABLE TO OPEN OR READ TURBINE FILE'
           ! CALL wrf_error_fatal('RSA UNABLE TO OPEN OR READ TURBINE FILE')
        ENDIF
     ELSE 
        !RSA for hard coded turbine info
        PRINT*,"RSA NO TURBINE FILE SPECIFIED, USING HARD CODED TURBINE LOCATION."
        x_turbine(1)      = 1000.0 
        y_turbine(1)      = 1000.0
        hub_height(1)     = 80.0  
        rotor_diameter(1) = 77.0 
        blade_length(1)   = 38.5*(1.0-0.0444) !PSU1.5 turbine
        theta_turbine(1)  = 180.0
        cutinspd(1) = 3 !! Added by MSG
        cutoutspd(1) = 25 !! Added by MSG
        PRINT *,"Reading file turbines(MSG)",x_turbine(1),y_turbine(1),"locs"
     ENDIF !RSA end read file or not

     !RSA convert x_turbine and y_turbine from lat/lon to x/y if need be
     IF (turb_loc_in_latlon) THEN

#ifdef DM_PARALLEL
        !Get lat/lon of origin (southwest corner) of domain, which is on 
        !master proc (ie mytask=0). Then broadcast to all other procs.
        IF (mytask .EQ. master) THEN
           lat1 = xlat_v(its,jts)
           lon1 = xlong_u(its,jts)
        ENDIF

        CALL MPI_BCAST(lat1,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(lon1,1,MPI_REAL,master,MPI_COMM_WORLD,ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) !RSA not sure if this is needed, or if it should be before or after BCAST
#endif

        !RSA parameters for map_set. This function defines the structure 'proj'
        !which is used below in latlon_to_ij. proj_code should be set to match
        !the projection used in namelist.wps (see share/module_llxy.F for more info)
        knowni = 0.0 !RSA always 0 (origin)
        knownj = 0.0 !RSA always 0 (origin)
        CALL map_set(proj_code=PROJ_LC,proj=proj,lat1=lat1,lon1=lon1, &
                     knowni=knowni,knownj=knownj,dx=dx,stdlon=stand_lon, &
                     truelat1=truelat1,truelat2=truelat2)

        DO n = 1, n_turbines
           CALL latlon_to_ij(proj,y_turbine(n),x_turbine(n),inew,jnew)
           !RSA inew and jnew are REAL indices from the origin (knowni,knownj),
           !which is at the southwest corner of the domain.
           !Multiplying by dx and dy gives the distance on the grid from the
           !origin to the given point.
           x_turbine(n) = inew*dx
           y_turbine(n) = jnew*dy
        ENDDO

        !RSA for debugging
        ! print*,'y_turbine=',y_turbine(1)
        ! print*,'x_turbine=',x_turbine(1)
        ! print*,'lat1=',lat1
        ! print*,'lon1=',lon1
        ! print*,'dx=',dx
        ! print*,'stand_lon=',stand_lon
        ! print*,'truelat1=',truelat1
        ! print*,'truelat2=',truelat2

     ENDIF
  ENDIF !RSA end first time GAD is turned on

!print*,'CHECKING THAT VARIABLES ARE NON ZERO EVERY TIME A RESTART FILE INITIATES'
!print*,x_turbine(:),y_turbine(:),'locs',hub_height(:),'hub height',rotor_diameter(:),'diameter',blade_length(:),'blades',theta_turbine(:),'angle',cutinspd(:),'cutin',cutoutspd(:),'cutout'

!____________________________________________________________________________________________________________________________
!
! 1. Find x, y locations for u, v, and w nodes within the patch
!____________________________________________________________________________________________________________________________
 ! PRINT*,'Starting part 1'
  i_start = its
  i_end   = MIN(ite,ide-1)
  j_start = jts
  j_end   = MIN(jte,jde-1)

  dxo2 = dx/2.0
  dyo2 = dy/2.0

  ! RSA/JDM for debugging
  ! DO i = i_start, i_end !only for visualizing instantaneous position of multiple turbines
  ! DO j = j_start, j_end
  ! DO k = 1, 20
  !   wp_ts(i,k,j) = 0.0
  ! ENDDO
  ! ENDDO
  ! ENDDO

  DO i = i_start, i_end+1
    u_x(i) = (i-1)*dx
    v_x(i) = dxo2 + u_x(i)
    w_x(i) = v_x(i)
  ENDDO

  DO j = j_start, j_end+1
    v_y(j) = (j-1)*dy
    u_y(j) = dyo2 + v_y(j)
    w_y(j) = u_y(j)
  ENDDO

  !RSA determine if turbine(n) is on patch
  xmin_tile=u_x(i_start)
  xmax_tile=u_x(i_end+1)
  ymin_tile=v_y(j_start)
  ymax_tile=v_y(j_end+1)

  DO n = 1, n_turbines

    !Define if turbine falls within tile neighborhood
    turbine_on_patch(n) = 0

    xmin_turb = x_turbine(n) - rotor_diameter(n)
    xmax_turb = x_turbine(n) + rotor_diameter(n)
    ymin_turb = y_turbine(n) - rotor_diameter(n)
    ymax_turb = y_turbine(n) + rotor_diameter(n)

    !If turbine rectangle overlaps with tile do the following
    !the if statement actually checks if there is no overlap
    IF ( .NOT. ( (xmax_tile .LT. xmin_turb) .OR.              &  
                 (xmin_tile .GE. xmax_turb) .OR.              &
                 (ymax_tile .LT. ymin_turb) .OR.              &  
                 (ymin_tile .GE. ymax_turb) )  ) THEN

      turbine_on_patch(n) = 1
    ENDIF
  ENDDO
  
  !RSA calculate x_inflow and y_inflow positions
  DO n = 1, n_turbines

    !Draw ray upstream from turbine from hub to rotor_diameter using turbine yaw
    !RSA recall theta_turbine is measured counterclockwise starting at 0 deg. for a
    !west-facing turbine (ie. 90 facing south, 180 facing east, 270 facing
    !north, 360 back to facing west)
    x_inflow(n) = -rotor_diameter(n)*cos(theta_turbine(n)*pi/180.0) + x_turbine(n)
    y_inflow(n) = -rotor_diameter(n)*sin(theta_turbine(n)*pi/180.0) + y_turbine(n)

    !See if inflow location is on patch
    inflow_on_patch(n) = 0

    IF ( .NOT. ( (xmax_tile .LT. x_inflow(n)) .OR.              &  
                 (xmin_tile .GE. x_inflow(n)) .OR.              &
                 (ymax_tile .LT. y_inflow(n)) .OR.              &  
                 (ymin_tile .GE. y_inflow(n)) )  ) THEN

      inflow_on_patch(n) = 1
 !     print*,'YES, INFLOW LOCATION IS ON PATCH'
    ENDIF
  ENDDO
!PRINT*,'Ending part 1'
!____________________________________________________________________________________________________________________________
!
! 2. If some part of the turbine is on this patch, compute the height dependent gridpoint position variables
!____________________________________________________________________________________________________________________________
!PRINT*,'Starting part 2'
  DO n = 1, n_turbines

    IF (turbine_on_patch(n) .EQ. 1) THEN

      DO j = j_start-1, j_end+1   
        DO i = i_start-1, i_end+1
           w_z(i,kts,j)=0.0
        ENDDO
      ENDDO
  
      zmax_turb = hub_height(n) + rotor_diameter(n)*0.5                                      

      DO k = kts+1,kte             ! compute 1D column of z-values to determine k_end for the tile
        w_z(i_start,k,j_start) = w_z(i_start,k-1,j_start) + 1.0/rdzw(i_start,k-1,j_start)
      ENDDO
        
      DO k=kts+1,kte-1             ! k_end = 1 index greater to allow for nonuniform heights within the tile
        IF ( (w_z(i_start,k-1,j_start) .LT. zmax_turb) .AND.    &
             (w_z(i_start,k,j_start) .GE. zmax_turb)         ) THEN
          k_end=k+1
        ENDIF
      ENDDO 

      DO j = j_start, j_end        ! compute the rest of w_z up to k_end
        DO k = kts+1, k_end+1
          DO i = i_start, i_end
            w_z(i,k,j) = w_z(i,k-1,j) + 1.0/rdzw(i,k-1,j)
          ENDDO
        ENDDO
      ENDDO

      DO j = j_start, j_end         ! compute the i_start and i_end ghost nodes
        DO k = kts+1, k_end+1
          w_z(i_start-1,k,j)=w_z(i_start,k,j)
          w_z(i_end+1,k,j)=w_z(i_end,k,j)
        ENDDO
      ENDDO

      DO k = kts+1, k_end+1           ! compute the j_start and j_end ghost nodes
        DO i = i_start, i_end  
          w_z(i,k,j_start-1) = w_z(i,k,j_start)
          w_z(i,k,j_end+1) = w_z(i,k,j_end)
        ENDDO
      ENDDO

      DO k = kts+1, k_end+1           ! compute get the corner ghost nodes
        w_z(i_start-1,k,j_start-1)=0.5*(w_z(i_start-1,k,j_start) + &
                                        w_z(i_start,k,j_start-1) )
        w_z(i_start-1,k,j_end+1)=0.5*(w_z(i_start-1,k,j_end) +     &
                                      w_z(i_start,k,j_start+1) )
        w_z(i_end+1,k,j_start-1)=0.5*(w_z(i_end+1,k,j_start) +     &
                                      w_z(i_end,k,j_start-1) ) 
        w_z(i_end+1,k,j_end+1)=0.5*(w_z(i_end+1,k,j_end) +     &
                                    w_z(i_end,k,j_end+1) ) 
      ENDDO

      DO i = i_start-1, i_end+1
        DO k = kts, k_end
          DO j = j_start-1, j_end+1
            c_z(i,k,j) = 0.5 * ( w_z(i,k,j) + w_z(i,k+1,j) )
          ENDDO
        ENDDO
      ENDDO

      DO i = i_start, i_end
        DO k = kts, k_end
          DO j = j_start, j_end
            dz(i,k,j) = w_z(i,k+1,j) - w_z(i,k,j)
          ENDDO
        ENDDO
      ENDDO

      DO j = j_start, j_end           ! compute u_z and v_z up to k_end
        DO k = kts, k_end
          DO i = i_start, i_end+1
            ! 0.5 for averaging in horizontal, 0.5 for averaging in vertical
            u_z(i,k,j) = 0.5*(   0.5*( w_z(i-1,k+1,j) + w_z(i,k+1,j) )   &
                               + 0.5*( w_z(i-1,k,j)   + w_z(i,k,j)   ) )
          ENDDO
        ENDDO
      ENDDO

      DO j = j_start, j_end+1           ! compute u_z and v_z up to k_end
        DO k = kts, k_end
          DO i = i_start, i_end
            ! 0.5 for averaging in horizontal, 0.5 for averaging in vertical
            v_z(i,k,j) = 0.5*(   0.5*( w_z(i,k+1,j-1) + w_z(i,k+1,j) )   &
                               + 0.5*( w_z(i,k,j-1)   + w_z(i,k,j)   ) )
          ENDDO
        ENDDO
      ENDDO

    ENDIF ! turbine neighborhood is on patch
  ENDDO ! n_turbines
!PRINT*,'Ending part 2'
!____________________________________________________________________________________________________________________________
!
! 3. If inflow location is on patch, calculate the wind speed and direction at that location.
! These values are exchanged across all processors.
!____________________________________________________________________________________________________________________________

! shift wind speed for each turbine (column of v0t) back by one unit to make room
! for the current v0 value 
!PRINT*,'Starting part 3'

  DO n = 1, n_turbines
    DO k = n_timeseries-1,1,-1
      v0t(k+1,n)=v0t(k,n)
    ENDDO
  ENDDO
!PRINT*,'Starting part 3.1'
  IF (wp_opt .EQ. 3 ) THEN
    DO n = 1, n_turbines
      DO k = n_timeseries-1,1,-1
        d0t(k+1,n)=d0t(k,n)
      ENDDO
    ENDDO
  ENDIF
!PRINT*,'Starting part 3.2'
  DO n = 1, n_turbines
    IF (inflow_on_patch(n) .EQ. 1 ) THEN
 !     PRINT*,'inflow patch for turbine ',n
      theta =  theta_turbine(n)*pi/180.0

        !using w locations, will calculate cell-centered values below
        i_inflow = minloc(abs(w_x(i_start:i_end)-x_inflow(n)),1)+i_start-1
        j_inflow = minloc(abs(w_y(j_start:j_end)-y_inflow(n)),1)+j_start-1

        !RSA find k index of inflow point by following a grid line
        !Note that w_z is AGL, ie w_z(:,kts,:)=0.0
        DO k = kts+1, k_end
          inflow_z(k) = w_z(i_inflow,k,j_inflow) - hub_height(n)
        ENDDO
        k_inflow=minloc(abs(inflow_z(kts+1:k_end)),1)+kts
        IF (inflow_z(k_inflow) .GT. 0.0) k_inflow = k_inflow - 1    

        u_c = 0.5 * ( u(i_inflow,k_inflow,j_inflow) + u(i_inflow+1,k_inflow,j_inflow) )
        v_c = 0.5 * ( v(i_inflow,k_inflow,j_inflow) + v(i_inflow,k_inflow,j_inflow+1) )
        w_c = 0.5 * ( w(i_inflow,k_inflow,j_inflow) + w(i_inflow,k_inflow+1,j_inflow) )
        v0inst = u_c*cos(theta)*cos(tilt) + v_c*sin(theta)*cos(tilt) - w_c*sin(tilt) !Instantaneous velocity component normal to the rotor.

        IF (wp_opt .EQ. 3 ) THEN
          IF ( ( u_c .GE. 0.0 ) .AND. ( v_c .LT. 0.0 ) ) THEN   !Wind direction, independent of the
             d0inst = 360.0 - (atan2(u_c,v_c)*180.0/pi - 90.0)  !rotor angle. 0 and 360 are Westerly. 
          ELSE                                                  !Positive values are counterclockwise
             d0inst = 90.0 - atan2(u_c,v_c)*180.0/pi            !up to 360. No negative values.
          ENDIF
        ENDIF

#ifdef DM_PARALLEL
        tag=n
        CALL MPI_SEND(v0inst,1,MPI_REAL,master,tag,            &
                             MPI_COMM_WORLD,ierr )
        IF (wp_opt .EQ. 3 ) THEN
          CALL MPI_SEND(d0inst,1,MPI_REAL,master,tag,            &
                                 MPI_COMM_WORLD,ierr )
        ENDIF
#else
        v0n(n) = v0inst
        IF (wp_opt .EQ. 3 ) THEN
          d0n(n) = d0inst
        ENDIF
#endif

    ENDIF ! inflow location on patch
  ENDDO ! n_turbines
!PRINT*,'Starting part 3.3'
#ifdef DM_PARALLEL
  IF ( mytask .eq. master ) THEN
    DO n=1,n_turbines
      tag=n
      CALL MPI_RECV(v0inst,1,MPI_REAL,MPI_ANY_SOURCE,tag,                 &
                             MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )
      v0n(n)=v0inst
    ENDDO
!PRINT*,'Starting part 3.4'
    IF (wp_opt .EQ. 3 ) THEN
      DO n=1,n_turbines
        tag=n
        CALL MPI_RECV(d0inst,1,MPI_REAL,MPI_ANY_SOURCE,tag,                 &
                               MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr )
        d0n(n)=d0inst
      ENDDO
    ENDIF
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(v0n,n_turbines,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  IF (wp_opt .EQ. 3 ) THEN
    CALL MPI_BCAST(d0n,n_turbines,MPI_REAL,master,MPI_COMM_WORLD,ierr)
  ENDIF
#endif 
!____________________________________________________________________________________________________________________________
!
! 4. Add current v0 to timeseries in v0t, compute current average, and calculate pitch and rpm
!____________________________________________________________________________________________________________________________
  DO n = 1, n_turbines
  
      turb_power_local(n) = 0.0 !RSA initializing for power calculation

      operate = 0.0

      dum = 0.0                         ! Looking for problems...
      DO k = 1,n_timeseries
        dum = dum + v0t(k,n)
      ENDDO
      !If this is the first time through the GAD routine,
      !Assign all rows (times) of each column (turbine number) the current v0 and d0
      IF ( dum .eq. 0. ) THEN  
        DO k = 1,n_timeseries          ! Looking for problems...
          v0t(k,n) = v0n(n)
          d0t(k,n) = d0n(n)
        ENDDO
!        print*,'v0t in part 4 for first time running GAD should be equal to v0n for each turbine',v0t(1,n),v0n(n)
      ENDIF 

      v0t(1,n)=v0n(n)                  ! Assign the first element the current v0 
      d0t(1,n)=d0n(n)

      dum = 0.0                        ! Looking for problems...
      DO k = 1,n_timeseries
         dum = dum + v0t(k,n)
      ENDDO
      v0t_avg = dum/real(n_timeseries) !Running time-average

      dum = 0.0                        ! Looking for problems...
      DO k = 1,n_timeseries
         dum = dum + d0t(k,n)
      ENDDO
      d0t_avg = dum/real(n_timeseries) !Running time-average
      
      !!!!! Turbine yawing !!!!!!

      IF (wp_opt .EQ. 3) THEN

        IF ( abs(acc_yaw_err(n)) .GE. yaw_thresh ) THEN ! If acc_yaw_err(n) exceeds threshold, turbine n begins yawing.
                                                        ! During yawing, no new entries are stored in d0t(n) and 
                                                        ! acc_yaw_err(n) remains fixed.Turbine n rotates toward
                                                        ! the time average direction before yawing began,yaw_new. 
                                                        ! The turbine rotates at yaw_rate until it is less than 
                                                        ! yaw_rate*dt from yaw_new. At this point turbine_theta is
                                                        ! set to yaw_new, yawing is set to 0, and everything starts over.
          yawing(n) = 1
        ELSE
          yawing(n) = 0 
        ENDIF

        IF (yawing(n) .EQ. 0 ) THEN

          delta = theta_turbine(n) - d0t_avg
          IF ( delta .LT. -180.0 ) THEN 
            delta = theta_turbine(n) + ( 360.0 - d0t_avg )
          ELSE IF (delta .GE. 180.0 ) THEN 
            delta = theta_turbine(n) - ( 360.0 + d0t_avg )
          ENDIF

          sgn = 1.0
          IF ( delta .LT. 0.0 ) sgn = -1.0
          yaw_err = sgn*(delta)*(delta)*dt
          acc_yaw_err(n) = acc_yaw_err(n) + yaw_err

        ELSE ! turbine is yawing or ready to start
         ! PRINT *,"Turbine is yawing (MSG), theta turbine is",theta_turbine(n)
         yaw_new = d0t_avg !calculated above
         yaw_inc = yaw_rate*dt

         delta = theta_turbine(n) - yaw_new  
         IF ( delta .LT. -180.0 ) THEN 
           delta = theta_turbine(n) + ( 360.0 - d0t_avg )
         ELSE IF (delta .GE. 180.0 ) THEN 
           delta = theta_turbine(n) - ( 360.0 + d0t_avg )
         ENDIF

         sgn = 1.0
         IF ( delta .GT. 0.0 ) sgn = -1.0

         IF ( abs(yaw_inc) .GE. abs(delta) ) THEN 
           theta_turbine(n) = yaw_new
           yawing(n) = 0
           acc_yaw_err(n) = 0.0
         ELSE
           theta_turbine(n) = theta_turbine(n) + sgn*yaw_inc
         ENDIF

        ENDIF !yawing or not

        ! Keep turbine-theta between 0 and 360 degrees.
        IF ( theta_turbine(n) .GE. 360.0 ) theta_turbine(n) = theta_turbine(n) - 360.0
        IF ( theta_turbine(n) .LT. 0.0 ) theta_turbine(n) = 360.0 + theta_turbine(n)

      ENDIF !wp_opt .EQ. 3

    !!!!! Turbine operation !!!!!!
 !   print*,'Is the turbine on path?',turbine_on_patch(n),'1-> Yes, 0-> No'
    IF (turbine_on_patch(n) .EQ. 1 ) THEN !only do this if part of the turbine is on the patch

      theta =  theta_turbine(n)*pi/180.0
  !    print*,'current wind speed is:',v0t_avg,', cut in speed is:',cutinspd(n),', and cut out speed is:',cutoutspd(n)
      IF ( (v0t_avg .lt. cutinspd(n)) .or. (v0t_avg .gt. cutoutspd(n)) ) GOTO 200 
      operate = 1.0
      
      !Use time average velocity at hub location for rpm and pitch
   !   print*,'TURBINE SHOULD BE OPERATING'
      IF (turbine_opt .EQ. 15) THEN
         rpm   = calc_rpm_psu15( v0t_avg )
         pitchdeg = calc_pitch_psu15( v0t_avg )
      ELSEIF (turbine_opt .EQ. 16) THEN
         rpm   = calc_rpm_nrel5( v0t_avg )
         pitchdeg = calc_pitch_nrel5( v0t_avg )
      ENDIF
   
      omega(n) = 2.0*pi*rpm/60.0

      ! DO j = j_start, j_end
      !   DO i = i_start, i_end
      !      rotrate(i,j) = 2.0*pi*rpm/60.0
      !      thrust(i,j)  = 0.0
      !      torque(i,j)  = 0.0
      !      power(i,j)   = 0.0
      !   ENDDO
      ! ENDDO

      !____________________________________________________________________________________________________________________________
      !
      ! 5.  Compute distances from computational nodes to rotor place in order to 
      !     ii) Smear forces as a function of distance in the axial direction
      !____________________________________________________________________________________________________________________________
      !
      ! Compute distances from computational nodes to rotor place in order to 
      !  ii) Smear forces as a function of distance in the axial direction

      ! dau,dav,daw - Signed distance, in a horizontal plane, from the gridpoint to a line
      !               going through the turbine (nacelle) axis. 
      ! dru,drv,drw - signed distance, in a horizontal plane, from the vertical rotor plane
      CALL calc_distance( x_turbine(n), y_turbine(n), theta_turbine(n), &
                          its, ite, jts, jte, &
                          i_start, i_end, j_start, j_end,        &
                          u_x, u_y, v_x, v_y, w_x, w_y,          &
                          dau, dav, daw, dru, drv, drw )

      DO j = j_start, j_end
        DO k = kts, k_end
          DO i = i_start, i_end
            dzu(i,k,j) = u_z(i,k,j) - hub_height(n)
            dzv(i,k,j) = v_z(i,k,j) - hub_height(n)
            dzw(i,k,j) = w_z(i,k,j) - hub_height(n)
            dzc(i,k,j) = c_z(i,k,j) - hub_height(n)
          ENDDO
        ENDDO
      ENDDO

      DO j = j_start, j_end
        DO k = kts, k_end
          DO i = i_start, i_end
            ru(i,k,j) = sqrt(dau(i,j)**2+dzu(i,k,j)**2)
            rv(i,k,j) = sqrt(dav(i,j)**2+dzv(i,k,j)**2)
            rw(i,k,j) = sqrt(daw(i,j)**2+dzw(i,k,j)**2)
            rc(i,k,j) = sqrt(daw(i,j)**2+dzc(i,k,j)**2)
            wp_acc_u(i,k,j) = 0.0
            wp_acc_v(i,k,j) = 0.0
            wp_acc_w(i,k,j) = 0.0
          ENDDO
        ENDDO
      ENDDO

      rotor_radius = 0.5*rotor_diameter(n)

      cnst = 1.0/rho !JDM

!      sigma = sqrt(dx*dy)
!      adenom = sigma*sqrt(2.0*pi)
!      bdenom  = 2.0*sigma*sigma
      sigma   = (abs(dx*cos(theta))+abs(dy*sin(theta))) !NM so normal for rotor plane sqrt(dx*dy)	!!!This is the standard deviation
      adenom  = sigma*sqrt(2.0*pi)                  !want integral=1 for Gaussian ;accel will have units of 1/m
      bdenom  = 2.*sigma*sigma                      !!determines spread of gaussian

      !NM sigmar normal to rotor plane
      sigmar = abs(dx*cos(theta)*cos(tilt))+abs(dy*sin(theta)*cos(tilt)) !NM
      adenom2 = sigmar*sqrt(2.0*pi) !NM want integral to be 1 for gaussian
      bdenom2 = 2.*sigmar**2 !NM spread of 1-D gaussian 

      radmin = rotor_diameter(n)/2.0 - blade_length(n)
      !____________________________________________________________________________________________________________________________
      !
      ! 6. Compute local lift and drag accelerations using instantaneous local velocity, and apply to u,v, and w tendencies
      !____________________________________________________________________________________________________________________________
      DO j = j_start, j_end
        DO k = kts, k_end
          DO i = i_start, i_end

             v0 = u(i,k,j)*cos(theta)*cos(tilt) + v(i,k,j)*sin(theta)*cos(tilt) - w(i,k,j)*sin(tilt)  

             v0_u = u(i,k,j) * cos(theta)*cos(tilt) &
                    + 0.25 * ( v(i,k,j) + v(i-1,k,j) + v(i-1,k,j+1) + v(i,k,j+1) ) * sin(theta)*cos(tilt) &
                    - 0.25 * ( w(i,k,j) + w(i-1,k,j) + w(i-1,k+1,j) + w(i,k+1,j) ) * sin(tilt)

             v0_v =   0.25 * ( u(i,k,j) + u(i,k,j-1) + u(i+1,k,j-1) + u(i+1,k,j) ) * cos(theta)*cos(tilt) &
                    + v(i,k,j) * sin(theta)*cos(tilt) &
                    - 0.25 * ( w(i,k,j) + w(i,k,j-1) + w(i,k+1,j-1) + w(i,k+1,j) ) * sin(tilt)
             
             IF (k .GT. kts) THEN
                v0_w =   0.5 * ( fnm(k) * ( u(i,k,j) + u(i+1,k,j) ) + fnp(k) * ( u(i,k-1,j) + u(i+1,k-1,j) ) ) * cos(theta)*cos(tilt) &
                       + 0.5 * ( fnm(k) * ( v(i,k,j) + v(i,k,j+1) ) + fnp(k) * ( v(i,k-1,j) + v(i,k-1,j+1) ) ) * sin(theta)*cos(tilt) &
                       - w(i,k,j) * sin(tilt)
             ELSE
                v0_w = 0.0 !RSA-JDM This should be the kinematic surface BC
             ENDIF

             v0_c =  0.5 * ( u(i,k,j) + u(i+1,k,j) ) * cos(theta)*cos(tilt) &
                   + 0.5 * ( v(i,k,j) + v(i,k,j+1) ) * sin(theta)*cos(tilt) &
                   - 0.5 * ( w(i,k,j) + w(i,k+1,j) ) * sin(tilt)

             !Tower------------------------------------------------------	
             IF (cdt .GT. 0.0 ) THEN
                IF (dzu(i,k,j) .le. -(radmin)) then !NM tower below hub
                   Rtwr = abs(dzu(i,k,j))*(Rtwrb-Rtwrt)/hub_height(n)+Rtwrt !NM linear	 
                   IF ( abs( dau(i,j) ) .le. Rtwr) then !NM tower below hub; dau needs to be less than tower radius
                   !JDM-04-201 Note to self. This should be fixed so that the tower is not grid location dependent. Same goes for the hub!
                      fn = 0.5*rho*v0_u*v0_u*cdt !NM force per unit area for tower
                      ft = 0.0 !NM no tangential force
                      acc_u = (1.0/adenom)*exp(-(dru(i,j)**2)/bdenom) !NM 1-D !Gaussian
                      fx = fn*cos(theta) !NM to x direction
                      ru_tendf(i,k,j) = ru_tendf(i,k,j) - muu(i,j)*acc_u*cnst*fx !NM
                   ENDIF
                ENDIF

                IF (dzv(i,k,j) .le. -(radmin)) then !NM tower below hub
                   Rtwr = abs(dzv(i,k,j))*(Rtwrb-Rtwrt)/hub_height(n)+Rtwrt !NM linear	 
                   IF ( abs( dav(i,j) ) .le. Rtwr) then !NM tower below hub; dau needs to be less than tower radius
                   !JDM-04-201 Note to self. This should be fixed so that the tower is not grid location dependent. Same goes for the hub!
                      fn = 0.5*rho*v0_v*v0_v*cdt !NM force per unit area for tower
                      ft = 0.0 !NM no tangential force
                      acc_v = (1.0/adenom)*exp(-(drv(i,j)**2)/bdenom) !NM 1-D !Gaussian
                      fy = fn*sin(theta) !NM to x direction
                      rv_tendf(i,k,j) = rv_tendf(i,k,j) - muv(i,j)*acc_v*cnst*fy !NM
                  ENDIF
                ENDIF
             ENDIF !cdt > 0

!------ X ---- 

        ! IF ( (ru(i,k,j) .le. radmin) .and. (ru(i,k,j) .gt. -0.001) ) then !HUB-----------------------------------	
          IF (ru(i,k,j) .le. radmin) then !HUB----------------------------------- !RSA

           fn = 0.5*rho*v0_u*v0_u*cdh  !NM only need force per unit area

           ft = 0.0 !NM no tangential force

           !NM including off set of hub from tower hoff in + dru
           !direction or opposite of wind direction U

           acc_u = (1.0/adenom)*exp(-((dru(i,j)-hoff)**2)/bdenom) !NM 1-D gaussian, hub has no tilt.
 
                fx = fn*cos(theta) 
 
                ru_tendf(i,k,j) = ru_tendf(i,k,j) &
                                  - muu(i,j)*acc_u*cnst*fx

!JDM	        thrust(i,k,j) = thrust(i,k,j)-fn*cos(theta)*acc_u     		
                !!thrust everywhere is equal to current thrust plus force distribution due to gaussian times normal force

        ELSE !GAD-------------------------------------------------------------------------------------------------

          IF (tilt .eq. 0.0) THEN !No Tilt-----------------------------------------------------------------------

             drut = dru(i,j)-hoff !distance only depends on x value!!!if theta is zero the normal distance will be x0

             radius = sqrt( dau(i,j)**2 + dzu(i,k,j)**2 )

             zeta = atan2( dzu(i,k,j), dau(i,j) )

          ELSE !Tilt ------------------------------------------------------------------------------------------

! drut-----------------

             drut = calc_yprime( 1.0/tan(tilt), 1.0, 0.0, (dru(i,j)-hoff), dzu(i,k,j) ) 

! zd-----------------

             zd = calc_yprime( -tan(tilt), 1.0, 0.0, (dru(i,j)-hoff), dzu(i,k,j) ) 

             radius = sqrt( dau(i,j)**2 + zd**2 ) !NM new radial distance to tilted turbine:

             zeta = atan2( zd, dau(i,j) )

           ENDIF !Tilt or No Tilt-------------------------------------------------------------------------------

             acc_u = (1.0/adenom2)*exp(-(drut**2)/bdenom2) !NM 1-D

             !RSA adding logic to prevent fx from being NaN. See explanation in calc_fnft below.
             if ( v0_u .gt. 0.0 ) then
                if (turbine_opt .eq. 15) then 
     !              print*,'radius',radius,'blade_length',blade_length(n),'rotor_diameter',rotor_diameter(n),'dscale',dscale
                   fnft = calc_fnft_psu15( v0_u, rho, blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale)
     !              print*,'done u-component' 
                elseif (turbine_opt .eq. 16) then
                   fnft = calc_fnft_nrel5( v0_u, rho, blades,blade_length(n),rotor_diameter(n),   &
                                              radius, pitchdeg, rpm,dscale,cscale)
                endif

                fn = fnft(1)
                ft = fnft(2)
     
                fx = fn*cos(theta)*cos(tilt) + ft*sin(zeta)*sin(theta) - ft*cos(zeta)*sin(tilt)*cos(theta)
              else
                !RSA set fx to 0 for now, but turbine should eventually yaw into place.
                ! write(*,*) "RSA to prevent NaNs, skipping call to calc_fnft for fx: i,j,k=",i,j,k 
                fx = 0.0
              endif
                
!JDM TEST SIGN 12/2014     fx = fn*cos(theta)*cos(tilt) 
!JDM TEST SIGN 12/2014     fx = ft*sin(zeta)*sin(theta) 
!JDM TEST SIGN 12/2014     fx = -ft*cos(zeta)*sin(tilt)*cos(theta)

!print*,'i,j,k',k,j,i
!print*,'v0, pitchdeg, rpm',v0, pitchdeg, rpm
!print*,'fn,ft,fx,accu',fn,ft,fx,acc_u

              ru_tendf(i,k,j) = ru_tendf(i,k,j) &
                                 - muu(i,j)*acc_u*cnst*fx

              ! wp_ts(i,k,j) = -1.0*muu(i,j)*acc_u*cnst*fx

        ENDIF !Hub or GAD-----------------------------------------------------------------------------------------

!------- Y ----

        ! IF ( (rv(i,k,j) .le. radmin) .and. (rv(i,k,j) .gt. -0.001) ) then !HUB------------------------------------	
        IF (rv(i,k,j) .le. radmin) then !HUB--------------------------- !RSA

           fn = 0.5*rho*v0_v*v0_v*cdh 

           ft = 0.0

           acc_v = (1.0/adenom)*exp(-((drv(i,j)-hoff)**2)/bdenom) 

           fy = fn*sin(theta) 

           rv_tendf(i,k,j) = rv_tendf(i,k,j) &
                                   - muv(i,j)*acc_v*cnst*fy

        ELSE !GAD-------------------------------------------------------------------------------------------------	

           IF (tilt .eq. 0.0) THEN !No Tilt-----------------------------------------------------------------------

               drvt = drv(i,j)-hoff !distance only depends on x value !!!if theta is zero the normal distance will be x0

               radius = sqrt( dav(i,j)**2 + dzv(i,k,j)**2 )

               zeta = atan2( dzv(i,k,j), dav(i,j) )

           ELSE !Tilt---------------------------------------------------------------------------------------------

               drvt = calc_yprime( 1.0/tan(tilt), 1.0, 0.0, (drv(i,j)-hoff), dzv(i,k,j) )  

               zd = calc_yprime( -tan(tilt), 1.0, 0.0, (drv(i,j)-hoff), dzv(i,k,j) )

               radius = sqrt( dav(i,j)**2 + zd**2 )

               zeta = atan2( zd, dav(i,j) )

         ENDIF !Tilt or No Tilt-------------------------------------------------------------------------------

               acc_v = (1.0/adenom2)*exp(-(drvt**2)/bdenom2) !NM 1-D

               !RSA adding logic to prevent fy from being NaN. See explanation in calc_fnft below.
               if ( v0_v .gt. 0.0 ) then
                   if (turbine_opt .eq. 15) then
      !                print*,'radius',radius,'blade_length',blade_length(n),'rotor_diameter',rotor_diameter(n),'dscale',dscale
                      fnft = calc_fnft_psu15( v0_v, rho, blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale)
       !               print*,'done v-component'
                   elseif (turbine_opt .eq. 16) then
                      fnft = calc_fnft_nrel5( v0_v, rho,blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale)
                   endif

                   fn = fnft(1)
                   ft = fnft(2)

                   fy = fn*sin(theta)*cos(tilt) - ft*sin(zeta)*cos(theta) - ft*cos(zeta)*sin(tilt)*sin(theta)!NM 
                else
                   !RSA set fy to 0 for now, but turbine should eventually yaw into place.
                   ! write(*,*) "RSA to prevent NaNs, skipping call to calc_fnft for fy: i,j,k=",i,j,k 
                   fy = 0.0
                endif

!JDM TEST SIGN 12/2014      fy = fn*sin(theta)*cos(tilt)
!JDM TEST SIGN 12/2014      fy = -ft*sin(zeta)*cos(theta)
!JDM TEST SIGN 12/2014      fy = -ft*cos(zeta)*sin(tilt)*sin(theta)

                rv_tendf(i,k,j) = rv_tendf(i,k,j) &
                                  - muv(i,j)*acc_v*cnst*fy

!JDM            wp_ts(i,k,j) = -1.0*muv(i,j)*acc_v*cnst*fy

              ENDIF !Hub or GAD------------------------------------------------------------------------------------

!-------- Z ----

              IF (tilt .eq. 0.0) THEN !No Tilt --------------------------------------------------------------------

                 drwt = drw(i,j)-hoff !distance only depends on x value 	!!!if theta is zero the normal distance will be x0

                 radius = sqrt( daw(i,j)**2 + dzw(i,k,j)**2 )

                 zeta = atan2( dzw(i,k,j), daw(i,j) ) 

              ELSE !Tilt ------------------------------------------------------------------------------------------

                 drwt = calc_yprime( 1.0/tan(tilt), 1.0, 0.0, (drw(i,j)-hoff), dzw(i,k,j) ) 
 
                 zd = calc_yprime( -tan(tilt), 1.0, 0.0, (drw(i,j)-hoff), dzw(i,k,j) ) 

                 radius = sqrt( daw(i,j)**2 + zd**2 )

                 zeta = atan2( zd, daw(i,j) ) 

              ENDIF !Tilt or No Tilt-------------------------------------------------------------------------------

              acc_w = (1.0/adenom2)*exp(-(drwt**2)/bdenom2) !NM 1-D

              !RSA adding logic to prevent fz from being NaN. See explanation in calc_fnft below.
              if ( v0_w .gt. 0.0 ) then
                 if (turbine_opt .eq. 15) then
        !            print*,'radius',radius,'blade_length',blade_length(n),'rotor_diameter',rotor_diameter(n),'dscale',dscale
                    fnft = calc_fnft_psu15( v0_w, rho, blades, blade_length(n), rotor_diameter(n), &
                                            radius, pitchdeg, rpm, dscale, cscale )
         !           print*,'done w-component' 
                 elseif (turbine_opt .eq. 16) then
                      fnft = calc_fnft_nrel5( v0_w, rho,blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale)
                 endif

                 fn = fnft(1)
                 ft = fnft(2) 

                 fz = -fn*sin(tilt) - ft*cos(zeta)*cos(tilt)
                if (isnan(fz)) then
                        print*,'fz is NaN fn=',fn,'sin(tilt)=',sin(tilt),'ft=',ft,'cos(zeta)=',cos(zeta),'cos(tilt)=',cos(tilt)
                endif
              else
                 !RSA set fz to 0 for now, but turbine should eventually yaw into place.
                 ! write(*,*) "RSA to prevent NaNs, skipping call to calc_fnft for fz: i,j,k=",i,j,k 
                 fz = 0.0
              endif

!JDM TEST SIGN      fz = -fn*sin(tilt) 
!JDM TEST SIGN      fz = -ft*cos(zeta)*cos(tilt) 
              temp = rw_tendf(i,k,j)
              rw_tendf(i,k,j) = rw_tendf(i,k,j)  &
                              - mut(i,j)*acc_w*cnst*fz
                if (isnan(rw_tendf(i,k,j))) then
                        print*,'in phys/module_gen_act_disk mut(i,j)=',mut(i,j),'acc_w=',acc_w,'cnst=',cnst,'fz=',fz,'temp=',temp
                endif
!JDM          wp_ts(i,k,j) = -1.0*mut(i,j)*acc_w*cnst*fz
 
!------ C ---- 

        ! IF ( (ru(i,k,j) .le. radmin) .and. (ru(i,k,j) .gt. -0.001) ) then !HUB-----------------------------------	
          IF (rc(i,k,j) .le. radmin) then !HUB------------------------------ !RSA

           ft = 0.0 !NM no tangential force

          ELSE !GAD-------------------------------------------------------------------------------------------------

            IF (tilt .eq. 0.0) THEN !No Tilt-----------------------------------------------------------------------

                drct = drw(i,j)-hoff !distance only depends on x value!!!if theta is zero the normal distance will be x0

                radius = sqrt( daw(i,j)**2 + dzc(i,k,j)**2 )

                zeta = atan2( dzc(i,k,j), daw(i,j) )

            ELSE !Tilt ------------------------------------------------------------------------------------------

                drct = calc_yprime( 1.0/tan(tilt), 1.0, 0.0, (drw(i,j)-hoff), dzc(i,k,j) ) 

                zd = calc_yprime( -tan(tilt), 1.0, 0.0, (drw(i,j)-hoff), dzc(i,k,j) ) 

                radius = sqrt( daw(i,j)**2 + zd**2 ) !NM new radial distance to tilted turbine:

                zeta = atan2( zd, daw(i,j) )

             ENDIF !Tilt or No Tilt-------------------------------------------------------------------------------

             acc_c = (1.0/adenom2)*exp(-(drct**2)/bdenom2) !NM 1-D

             !RSA adding logic to prevent fx from being NaN. See explanation in calc_fnft below.
             if ( v0_c .gt. 0.0 ) then
                if (turbine_opt .eq. 15) then
!                   print*,'radius',radius,'blade_length',blade_length(n),'rotor_diameter',rotor_diameter(n),'dscale',dscale
                   fnft = calc_fnft_psu15( v0_c, rho, blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale) 
 !                  print*,'done c-component'
                elseif (turbine_opt .eq. 16) then
                      fnft = calc_fnft_nrel5( v0_c, rho,blades, blade_length(n), rotor_diameter(n),   &
                                              radius, pitchdeg, rpm, dscale, cscale)
                endif

                ft = fnft(2)

              else
                !RSA set fx to 0 for now, but turbine should eventually yaw into place.
                ! write(*,*) "RSA to prevent NaNs, skipping call to calc_fnft for fx: i,j,k=",i,j,k 
                ft = 0.0
              endif
                
        ENDIF !Hub or GAD-----------------------------------------------------------------------------------------

!!         Assuming turbine rotors are perpendicualar to x axis. Generalize to arbitrary orientation.  
!
!!            IF ( ru(i,k,j) .lt. rotor_radius ) THEN
!               thrust(i,j) = thrust(i,j) +                             &
!                             wp_acc_u(i,k,j)*fn
!               torque(i,j) = torque(i,j) +                             &
!                             wp_acc_u(i,k,j)*ru(i,k,j)*ft
!                power(i,j)  = power(i,j)  +                             &
!                             wp_acc_u(i,k,j)*rotrate(i,j)*ru(i,k,j)*ft
!!            ENDIF
!!

              !RSA sum up power for each turbine on local processor
              turb_power_local(n) = turb_power_local(n) + omega(n)*ft*acc_c*rc(i,k,j)*dx*dy*dz(i,k,j)

          ENDDO         ! j loop
        ENDDO         ! k loop
      ENDDO       ! i loop

    ENDIF !turbine_on_patch = 1

200 CONTINUE !time-averaged inflow is < cutin or > cutout speed
    
    if (operate .eq. 1.0) then
 !      print*,'Turbine did start operating'            
    elseif (operate .eq. 0.0) then
  !     print*,'Turbine did not start operating'
    endif

  ENDDO  !n loop over all turbines
!print*, 'finishing looping all turbines'
#ifdef DM_PARALLEL
  !RSA sum turbine power over all procs
!  print*, 'try setting turb_power to zero (allocate a space for it)'
  CALL MPI_REDUCE(turb_power_local,turb_power,n_turbines,MPI_REAL,MPI_SUM,master,MPI_COMM_WORLD,ierr)
#endif
!PRINT *,'GAD FINISHED, MSG MADE IT!!'

END SUBROUTINE gen_act_disk

!==============================================================================

SUBROUTINE calc_distance( x_turb, y_turb, theta_turb,     &
                          its, ite, jts, jte,             &
                          i_start, i_end, j_start, j_end, &
                          u_x, u_y, v_x, v_y, w_x, w_y,   &
                          dau, dav, daw, dru, drv, drw    )

!-----------------------------------------------------------------------------
!
! Purpose: To calculate distances from the axis and rotor plane in (x,y)
! plane. 
!
!  October, 2010: Branko Kosovic, NCAR
!  Summer, 2013: JDM. Modified for consistent sign relative to quadrant.
!
!=============================================================================

  REAL, INTENT( IN )                                 :: x_turb
  REAL, INTENT( IN )                                 :: y_turb
  REAL, INTENT( IN )                                 :: theta_turb
  INTEGER, INTENT( IN )                              :: its, ite, jts, jte
  INTEGER, INTENT( IN )                              :: i_start, i_end
  INTEGER, INTENT( IN )                              :: j_start, j_end

  REAL, DIMENSION( its:ite ), INTENT( IN )           :: u_x, v_x, w_x
  REAL, DIMENSION( jts:jte ), INTENT( IN )           :: u_y, v_y, w_y

  REAL, DIMENSION( its:ite, jts:jte ), INTENT( OUT ) :: dau, dav, daw !signed distance from each gridpoint to axis centerline
  REAL, DIMENSION( its:ite, jts:jte ), INTENT( OUT ) :: dru, drv, drw !signed distance from each gridpoint to rotor plane

! LOCAL VARIABLES

  REAL    :: aa, ba, ca
  REAL    :: ar, br, cr
  REAL    :: x0u, x0v, x0w
  REAL    :: y0u, y0v, y0w
  INTEGER :: i, j

  theta =  theta_turb*pi/180.0

! rp indicates rotor plane 

!  z                                               y    
!  |             rp                          rp    |      tan(theta) = y/x  
!                |                            \    |               .    
!                |                             \   |           . ^      
!                |                              \  |       .     |       
!                |___                            \ |   .   theta | 
!    .  .  .  .  |___|  .  .  .       ____________\|_____________|__ x
!        |       | ||                          .   |\
!    dau |       | ||                      .       | \
!        |       | ||                  .  \        |  \
!        p _ _ _ | ||                      \       | . \
!           dru    ||                   dau \    . |dru \
!                  ||                        p     |
!__________________||___ surface                            


!
! 1. Obtain coefficients for normal distance, dau, to axial plane, the plane 
! containing the dotted lines above, normal to the center of the rotor disk. 
! tan(theta) = y/x is the slope of this line in the x-y plane.
!
! If looking from a point p upwind of the turbine and facing it, positive
! distance is to the right and negative is to the left. Switching this
! sign convention along with the zd sign convention would be needed to get
! assumtion that ft rotates clockwise and not counterclockwise 	

  IF ( ( theta .LT. 0.5*pi ) .OR. (theta .GT. 1.5*pi) ) THEN   ! Quadrants 1 or 4 	
        aa = tan(theta)                                          
        ba = -1.0                      !neg because for ax + by + c = 0, slope m = -a/b                                        
     ELSEIF ( theta .EQ. 0.5*pi ) THEN                   ! INFINITE SLOPE
        aa = 1.0                       !x coordinate of point p 
        ba = 0.0                      
     ELSEIF ( theta .EQ. 1.5*pi ) THEN                   ! INFINITE SLOPE
        aa = -1.0                      !negative of x coordinate of point p 
        ba = 0.0
     ELSE                                                  ! Quadrants 2 or 3 	
        aa = -tan(theta)
        ba = 1.0
     ENDIF
     ca = 0.0                    !no slope intercept, center of rotor plane is the origin

! 2. Obtain coefficients for normal distance to rp, dru.
!    Negative downwind and positive (+) upwind

     IF ( theta .EQ. 0.5*pi ) THEN
        ar = 0.0
        br = -1.0                  !distance only depends on y value
     ELSEIF ( theta .EQ. 0.) THEN 
        ar = -1.0                  !distance only depends on x value, if theta is zero the normal distance will be x0
        br = 0.0  
     ELSEIF ( theta .LT. pi ) THEN 
        ar = -1.0/tan(theta)         !Since perpendicular, take negative reciprocal for slope
        br = -1.0                  !Must correspond to shortest line
     ELSEIF ( theta .EQ. 1.5*pi ) THEN
        ar = 0.0
        br = 1.0                   !distance only depends on y value
     ELSEIF ( theta .EQ. pi ) THEN 
        ar = 1.0                   !distance only depends on x value
        br = 0.0
     ELSE                          !theta between pi and 2pi
        ar = 1.0/tan(theta)
        br = 1.0
     ENDIF
     cr = 0.0

    DO j = j_start, j_end

      y0u= u_y(j) - y_turb
      y0v= v_y(j) - y_turb
      y0w= w_y(j) - y_turb

      DO i = i_start, i_end

         x0u= u_x(i) - x_turb
         x0v= v_x(i) - x_turb
         x0w= w_x(i) - x_turb

!JDM calc_distline returns unsigned distance.
!JDM ! calculate distance to the turbine axis in (x,y) plane
!JDM           dau(i,j)=calc_distline(aa,ba,ca,x0u,y0u)
!JDM           dav(i,j)=calc_distline(aa,ba,ca,x0v,y0v)
!JDM           daw(i,j)=calc_distline(aa,ba,ca,x0w,y0w)
!JDM 
!JDM ! calculate distance to the rotor plane of the turbine in (x,y) plane
!JDM            dru(i,j)=calc_distline(ar,br,cr,x0u,y0u)
!JDM            drv(i,j)=calc_distline(ar,br,cr,x0v,y0v)
!JDM            drw(i,j)=calc_distline(ar,br,cr,x0w,y0w)

!JDM calc_yprime returns signed distance. 

! calculate distance to the turbine axis in (x,y) plane

         dau(i,j)=calc_yprime(aa,ba,ca,x0u,y0u)
         dav(i,j)=calc_yprime(aa,ba,ca,x0v,y0v)
         daw(i,j)=calc_yprime(aa,ba,ca,x0w,y0w)
 
! calculate distance to the rotor plane of the turbine in (x,y) plane

         dru(i,j)=calc_yprime(ar,br,cr,x0u,y0u)
         drv(i,j)=calc_yprime(ar,br,cr,x0v,y0v)
         drw(i,j)=calc_yprime(ar,br,cr,x0w,y0w)

       ENDDO
    ENDDO

END SUBROUTINE calc_distance

!==============================================================================
!
!      Functions to compute distance from a line in 2d
!
!----------------------------------------------------------------------
FUNCTION calc_distline( a, b, c, x0, y0 ) result( distline ) 
!
    real, intent( in )   :: a, b, c
    real, intent( in )   :: x0, y0
    real                 :: distline
!   
    distline=abs(a*x0+b*y0+c)/sqrt(a*a+b*b)
!
END FUNCTION calc_distline

!==============================================================================
!
!      Function to compute signed distance from a line in 2D
!
!----------------------------------------------------------------------
FUNCTION calc_yprime( a, b, c, x0, y0 ) result( yprime ) 
!
    real, intent( in )   :: a, b, c
    real, intent( in )   :: x0, y0
    real                 :: yprime
!   
    yprime = (a*x0+b*y0+c)/sqrt(a*a+b*b)
!
END FUNCTION calc_yprime

!-------------------------------------------------------------

!==============================================================================
!
!      Functions to check a variable is NaN or Inf
!
!----------------------------------------------------------------------
FUNCTION check_naninf( rnaninf, snaninf ) result( naninf )
!
    real, intent( in )             :: rnaninf
    character(len=8), intent( in ) :: snaninf
    integer                        :: naninf
!
    naninf=0
!
    IF ( rnaninf*0 .NE. 0 ) THEN
      print *, ' Variable ',snaninf,' = ',rnaninf,' is Inf'
      naninf=1
    ENDIF
!   
    IF ( rnan .NE. rnan ) THEN
      print *, ' Variable ',srnaninf,' = ',rnaninf,' is NaN'
      naninf=2
    ENDIF
!_________________________________________________________________________
  
END FUNCTION check_naninf
!_________________________________________________________________________
  
END MODULE module_gen_act_disk

