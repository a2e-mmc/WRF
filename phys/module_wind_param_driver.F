
MODULE module_wind_param_driver

CONTAINS

!------------------------------------------------------------------

SUBROUTINE wind_param_driver( grid, config_flags,           &
                              ru_tendf, rv_tendf, rw_tendf, &
                              u_2, v_2, w_2, rdzw,          &
                              muu, muv, mut,                &
                              dx, dy, dt,                   & !JDM May 2, 2013
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
                              wp_acc_u, wp_acc_v, wp_acc_w, &
                              wp_ts,                        & !JDM May 2, 2013
                              thrust, torque,               &
                              rotrate, power,               &
                              n_turbines,                   &
                              n_timeseries,                 &
                              v0t, d0t, itimestep,          & !JDM May 2, 2013
                              xlat_v, xlong_u, hgt,         & !RSA
                              fnm, fnp,                     & !RSA
                              ids,ide,jds,jde,kds,kde,      & 
                              ims,ime,jms,jme,kms,kme,      &
                              ips,ipe,jps,jpe,kps,kpe,      &
                              its,ite,jts,jte,kts,kte       )

!------------------------------------------------------------------

! Driver layer modules
  USE module_domain
  USE module_configure
  USE module_tiles
  USE module_dm
  USE module_machine
  USE module_state_description
! Model layer modules
  USE module_bc

!! *** add new modules of schemes here

  USE module_gen_act_disk

  IMPLICIT NONE

! Input data.

  TYPE(domain) , TARGET          :: grid

  TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme), INTENT( INOUT ) &
  ::  ru_tendf,  & ! u tendency * mass           (kg s-4)                 
      rv_tendf,  & ! v                                                            
      rw_tendf     ! w

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme), INTENT( IN  ) &
  ::  u_2,       & ! u-velocity                  (m s-1)             
      v_2,       & ! v                                                              
      w_2,       & ! w                                                           
      rdzw         ! inverse vertical gid spacing at w-nodes (m-1)

  REAL, DIMENSION( ims:ime, jms:jme) ,         INTENT( IN  ) &
  ::  muu,       & ! column pressure at u nodes  (kg m-1 s-2)                     
      muv,       & !                    v                                                                
      mut,       & !                    w
      xlat_v,    &
      xlong_u,   &
      hgt

  REAL, DIMENSION( kms:kme) ,         INTENT( IN  ) &
  :: fnm, fnp

  REAL, DIMENSION( n_turbines ),                 INTENT( INOUT ) &
  :: x_turbine       & !                         (m)
   , y_turbine       & !                         (m)
   , x_inflow        & !                         (m)
   , y_inflow        & !                         (m)
   , hub_height      & !                         (m)
   , rotor_diameter  & !                         (m)
   , blade_length    & !                         (m)
   , theta_turbine   & !                         (deg)
   , cutinspd        & !                         (m/s)
   , cutoutspd       & !                         (m/s)
   , turb_power        !                         (kg*m2/s3)

 REAL, DIMENSION( n_turbines ),                 INTENT( INOUT ) &
  :: acc_yaw_err

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT( INOUT ) &
  :: wp_acc_u    & ! acceleration for u          (m s-2)
   , wp_acc_v    & !                  v          (m s-2)
   , wp_acc_w      !                  w          (m s-2)

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),  INTENT( INOUT ) & !JDM May 2, 2013
  :: wp_ts         ! 


  REAL, DIMENSION( ims:ime, jms:jme ),           INTENT( OUT ) &
  :: thrust          & !                         (N)
   , torque          & !                         (N m)
   , rotrate         & !                         (rad s-1)
   , power             !                         (W) 

  REAL, DIMENSION( n_timeseries, n_turbines ),   INTENT( INOUT ) &
  :: v0t         & !  
   , d0t
   
  REAL,                                           INTENT( IN  ) &
  :: dx          & !grid spacing in x-direction  (m)
   , dy            !                y

  REAL,                                           INTENT( IN  ) & !JDM May 2, 2013
  :: dt 

  INTEGER,                                        INTENT( IN  ) &
  :: n_turbines

  INTEGER,                                        INTENT( IN  ) &
  :: n_timeseries

  INTEGER, INTENT(IN) :: itimestep

  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde,               &
                         ims,ime,jms,jme,kms,kme,               &
                         ips,ipe,jps,jpe,kps,kpe,               &
                         its,ite,jts,jte,kts,kte

!local variables

  INTEGER :: ij, i, j, k

  IF ( config_flags%wp_opt .GE. 2 ) then

     CALL gen_act_disk ( ru_tendf, rv_tendf, rw_tendf,     & 
                         config_flags%wp_opt,              & 
                         config_flags%turbine_opt,         & !RSA
                         config_flags%windturb_spec,       & !RSA
                         config_flags%turb_loc_in_latlon,  & !RSA
                         config_flags%truelat1,            & !RSA
                         config_flags%truelat2,            & !RSA
                         config_flags%stand_lon,           & !RSA
                         n_turbines,                       &
                         n_timeseries,                     &
                         x_turbine,                        &
                         y_turbine,                        &
                         x_inflow,                         & !RSA
                         y_inflow,                         & !RSA
                         hub_height,                       &
                         rotor_diameter,                   &
                         blade_length,                     &
                         theta_turbine,                    &
                         cutinspd, cutoutspd,              & !RSA
                         turb_power,                       & !RSA
                         acc_yaw_err,                      & 
                         wp_acc_u, wp_acc_v, wp_acc_w,     & 
                         wp_ts,                            & !JDM May 2, 2013
                         thrust, torque,                   &
                         rotrate, power,                   &
                         u_2, v_2, w_2,                    & 
                         muu, muv, mut,                    & 
                         rdzw, dx, dy, dt,                 & !JDM May 2, 2013
                         v0t, d0t, itimestep,              & !JDM May 2, 2013
                         xlat_v, xlong_u, hgt,             & !RSA
                         fnm, fnp,                         & !RSA
                         ids, ide, jds, jde, kds, kde,     &
                         ims, ime, jms, jme, kms, kme,     &
                         ips, ipe, jps, jpe, kps, kpe,     &
                         its, ite, jts, jte, kts, kte      )

  ENDIF




RETURN

END SUBROUTINE wind_param_driver


END MODULE module_wind_param_driver
