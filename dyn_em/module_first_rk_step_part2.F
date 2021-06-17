!WRF:MEDIATION_LAYER:SOLVER

#define BENCH_START(A)
#define BENCH_END(A)

MODULE module_first_rk_step_part2

CONTAINS

  SUBROUTINE first_rk_step_part2 (   grid , config_flags        &
                             , moist , moist_old , moist_tend   &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , ru_tendf, rv_tendf               &
                             , rw_tendf, t_tendf                &
                             , ph_tendf, mu_tendf               &
                             , tke_tend                         &
                             , q_sq_prog_tend                   &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , gz1oz0 , chklowq   &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , t_phy           &
                             , dz8w , p8w , t8w                 &
                             , nba_mij, n_nba_mij               & !JDM
                             , nba_rij, n_nba_rij               & !JDM
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx,imex,jmsx,jmex,kmsx,kmex    &
                             , ipsx,ipex,jpsx,jpex,kpsx,kpex    &
                             , imsy,imey,jmsy,jmey,kmsy,kmey    &
                             , ipsy,ipey,jpsy,jpey,kpsy,kpey    &
                             , k_start , k_end                  &
                            )
    USE module_state_description
    USE module_model_constants
    USE module_domain, ONLY : domain
    USE module_configure, ONLY : grid_config_rec_type, model_config_rec
    USE module_big_step_utilities_em, ONLY : conv_t_tendf_to_moist
#ifdef DM_PARALLEL
    USE module_dm, ONLY : local_communicator, mytask, ntasks, ntasks_x, ntasks_y, local_communicator_periodic, &
                          wrf_dm_maxval, wrf_err_message, local_communicator_x, local_communicator_y
    USE module_comm_dm, ONLY : halo_em_tke_c_sub,halo_em_tke_d_sub,halo_em_tke_e_sub            &
            ,halo_em_phys_pbl_sub,halo_em_phys_shcu_sub &
            ,halo_em_fdda_sub,halo_em_phys_diffusion_sub,halo_em_tke_3_sub &
            ,halo_em_tke_5_sub,halo_obs_nudge_sub,period_bdy_em_a1_sub,period_bdy_em_phy_bc_sub &
            ,period_bdy_em_fdda_bc_sub,period_bdy_em_chem_sub,halo_em_phys_cu_sub,halo_em_helicity_sub &
            ,halo_em_phys_state_pbl3d_sub,halo_em_phys_substep_pbl3d_sub,halo_em_phys_diffusion_pbl3d_sub &
            ,halo_em_scalar_diffusion_pbl3d_sub,halo_em_pbl3d_3_sub,halo_em_pbl3d_5_sub &
            ,period_bdy_em_state_pbl3d_sub,period_bdy_em_substep_pbl3d_sub,period_bdy_em_phy_bc_pbl3d_sub &
            ,period_bdy_em_scalar_bc_pbl3d_sub
#endif

    USE module_driver_constants
    USE module_diffusion_em, ONLY : phy_bc, phy_bc_state_pbl3d, phy_bc_pbl3d, scalar_bc_pbl3d, cal_deform_and_div, compute_diff_metrics, &
                                    vertical_diffusion_2, horizontal_diffusion_2, calculate_km_kh, &
                                    tke_rhs, cal_helicity, &
                                    nonlocal_flux,meso_length_scale,free_atmos_length, &
                                    vertical_diffusion_implicit !XZ 
    USE module_em, ONLY : calculate_phy_tend
    USE module_fddaobs_driver, ONLY : fddaobs_driver
    USE module_bc, ONLY : set_physical_bc3d, set_physical_bc2d
    USE module_physics_addtendc, ONLY : update_phy_ten

    USE module_sfs_driver !JDM
    USE module_stoch, ONLY : update_stoch_ten, perturb_physics_tend,RAND_PERT_UPDATE
    USE module_scpm_jdm
    USE module_pbl3d, ONLY : Init_temp_arrays_for_substep, Update_wrf_tends_temp_state_and_zero_tends, Allocate_temp_tends, Deallocate_temp_tends, &
                             Calc_turb_fluxes_driver, Vertical_turb_mix, Horizontal_turb_mix


    IMPLICIT NONE

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           ips, ipe, jps, jpe, kps, kpe,     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,    &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                           imsy,imey,jmsy,jmey,kmsy,kmey,    &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey



    LOGICAL ,INTENT(IN)                        :: adapt_step_flag
    REAL, INTENT(IN)                           :: curr_secs

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(IN)      :: moist_old
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_fdda3d),INTENT(INOUT)  :: fdda3d
    REAL    ,DIMENSION(ims:ime,1:1,jms:jme,num_fdda2d),INTENT(INOUT)      :: fdda2d
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psim
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psih
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gz1oz0
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chklowq
    LOGICAL ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cu_act_flag
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: hol

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t8w

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: tke_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: q_sq_prog_tend

    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: mu_tendf

    INTEGER , INTENT(IN)                          ::  k_start, k_end

!JDM
  INTEGER, INTENT(  IN ) :: n_nba_mij, n_nba_rij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_mij) &
  :: nba_mij

  REAL ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,n_nba_rij) &
  :: nba_rij

! Local

    REAL, DIMENSION(:,:,:), ALLOCATABLE :: ru_tendf_tmp,   &
                                           rv_tendf_tmp,   &
                                           rw_tendf_tmp,   &
                                           t_tendf_tmp,    &
                                           moist_tend_tmp, &
                                           q_sq_tend_tmp

    REAL, DIMENSION( ims:ime, jms:jme ) :: ht_loc
    REAL :: scale_factor
    INTEGER, DIMENSION( ims:ime, jms:jme ) :: shadowmask 
    INTEGER                             :: ij, ss
    INTEGER  num_roof_layers
    INTEGER  num_wall_layers
    INTEGER  num_road_layers
    INTEGER  iswater
    INTEGER  rk_step 

#if ( WRF_DFI_RADAR == 1 )
    INTEGER  i_start,i_end,j_start,j_end,i,j,k
#endif

 ! initialize all tendencies to zero in order to update physics
 ! tendencies first (separate from dry dynamics).

   rk_step = 1

      IF ((grid%skebs_on==1).and.(grid%id .EQ. 1 )) then
          ! update and backtransform T
          CALL RAND_PERT_UPDATE(grid,'T',                                     &
                          grid%SPTFORCS,grid%SPTFORCC,                        &
                          grid%SPT_AMP,grid%ALPH_T,                           &
                          ips, ipe, jps, jpe, kps, kpe,                       &
                          ids, ide, jds, jde, kds, kde,                       &
                          ims, ime, jms, jme, kms, kme,                       &
                          k_start, k_end,                                     &
                          imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                          ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                          imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                          ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          config_flags%restart, grid%iseedarr_skebs,          &
                          config_flags%seed_dim,                              &
                          grid%DX,grid%DY,grid%skebs_vertstruc,               &
                          grid%rt_tendf_stoch,                                &
                          grid%stddev_cutoff_sppt,grid%gridpt_stddev_sppt, & 
                          grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPUV      )
          ! Update streamfunction, backtransform U
           CALL RAND_PERT_UPDATE(grid,'U',                                    &
                           grid%SPSTREAMFORCS,grid%SPSTREAMFORCC,             &
                           grid%SPSTREAM_AMP,grid%ALPH_PSI,                   &
                           ips, ipe, jps, jpe, kps, kpe,                      &
                           ids, ide, jds, jde, kds, kde,                      &
                           ims, ime, jms, jme, kms, kme,                      &
                           k_start, k_end,                                    &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                     &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                     &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                     &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           config_flags%restart, grid%iseedarr_skebs,         &
                           config_flags%seed_dim,                             &
                           grid%DX,grid%DY,grid%skebs_vertstruc,              &
                           grid%ru_tendf_stoch,                               &
                           grid%stddev_cutoff_sppt,grid%gridpt_stddev_sppt, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPUV     )
          ! Don't update streamfunction, backtransform V
           CALL RAND_PERT_UPDATE(grid,'V',                                    &
                           grid%SPSTREAMFORCS,grid%SPSTREAMFORCC,             &
                           grid%SPSTREAM_AMP,grid%ALPH_PSI,                   &
                           ips, ipe, jps, jpe, kps, kpe,                      &
                           ids, ide, jds, jde, kds, kde,                      &
                           ims, ime, jms, jme, kms, kme,                      &
                           k_start, k_end,                                    &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                     &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                     &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                     &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           grid% num_stoch_levels,grid% num_stoch_levels,     &
                           config_flags%restart, grid%iseedarr_skebs,         &
                           config_flags%seed_dim,                             &
                           grid%DX,grid%DY,grid%skebs_vertstruc,              &
                           grid%rv_tendf_stoch,                               &
                           grid%stddev_cutoff_sppt,grid%gridpt_stddev_sppt, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT      )
       ENDIF !skebs_on 

     if ((grid%sppt_on==1).and.(grid%id .EQ. 1 )) then
          CALL RAND_PERT_UPDATE(grid,'T',                                     &
                          grid%SPPTFORCS,grid%SPPTFORCC,                      &
                          grid%SPPT_AMP,grid%ALPH_SPPT,                       &
                          ips, ipe, jps, jpe, kps, kpe,                       &
                          ids, ide, jds, jde, kds, kde,                       &
                          ims, ime, jms, jme, kms, kme,                       &
                          k_start, k_end,                                     &
                          imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                          ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                          imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                          ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          grid%num_stoch_levels,grid%num_stoch_levels,        &
                          config_flags%restart, grid%iseedarr_sppt,           &
                          config_flags%seed_dim,                              &
                          grid%DX,grid%DY,grid%sppt_vertstruc,                &
                          grid%rstoch,                                        &
                          grid%stddev_cutoff_sppt,grid%gridpt_stddev_sppt, & 
                          grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
       ENDIF !sppt_on

      if ((grid%rand_perturb_on==1).and.(grid%id .EQ. 1 )) then
           CALL RAND_PERT_UPDATE(grid,'T',                                     &
                           grid%SPFORCS,grid%SPFORCC,                          &
                           grid%SP_AMP,grid%ALPH_RAND,                         &
                           ips, ipe, jps, jpe, kps, kpe,                       &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           k_start, k_end,                                     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           config_flags%restart, grid%iseedarr_rand_pert,      &
                           config_flags%seed_dim,                              &
                           grid%DX,grid%DY,grid%rand_pert_vertstruc,           &
                           grid%RAND_PERT,                                     &
                           grid%stddev_cutoff_rand_pert,grid%gridpt_stddev_rand_pert, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
      ENDIF !rand_perturb_on
      if ((grid%spp_conv==1).and.(grid%id .EQ. 1 )) then
           CALL RAND_PERT_UPDATE(grid,'T',                                     &
                           grid%SPFORCS2,grid%SPFORCC2,                        &
                           grid%SP_AMP2,grid%ALPH_RAND2,                       &
                           ips, ipe, jps, jpe, kps, kpe,                       &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           k_start, k_end,                                     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           config_flags%restart, grid%iseedarr_spp_conv,       &
                           config_flags%seed_dim,                              &
                           grid%DX,grid%DY,grid%vertstruc_spp_conv,            &
                           grid%pattern_spp_conv,                              & 
                           grid%stddev_cutoff_spp_conv,grid%gridpt_stddev_spp_conv, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
      ENDIF !spp_conv
      if ((grid%spp_pbl==1).and.(grid%id .EQ. 1 )) then
           CALL RAND_PERT_UPDATE(grid,'T',                                     &
                           grid%SPFORCS3,grid%SPFORCC3,                        &
                           grid%SP_AMP3,grid%ALPH_RAND3,                       &
                           ips, ipe, jps, jpe, kps, kpe,                       &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           k_start, k_end,                                     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           config_flags%restart, grid%iseedarr_spp_pbl,        &
                           config_flags%seed_dim,                              &
                           grid%DX,grid%DY,grid%vertstruc_spp_pbl,             &
                           grid%pattern_spp_pbl,                               & 
                           grid%stddev_cutoff_spp_pbl,grid%gridpt_stddev_spp_pbl, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
      ENDIF !spp_pbl
      if ((grid%spp_lsm==1).and.(grid%id .EQ. 1 )) then
           CALL RAND_PERT_UPDATE(grid,'T',                                     &
                           grid%SPFORCS4,grid%SPFORCC4,                        &
                           grid%SP_AMP4,grid%ALPH_RAND4,                       &
                           ips, ipe, jps, jpe, kps, kpe,                       &
                           ids, ide, jds, jde, kds, kde,                       &
                           ims, ime, jms, jme, kms, kme,                       &
                           k_start, k_end,                                     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,                      &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,                      &
                           imsy,imey,jmsy,jmey,kmsy,kmey,                      &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey,                      &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           grid%num_stoch_levels,grid%num_stoch_levels,        &
                           config_flags%restart, grid%iseedarr_spp_lsm,        &
                           config_flags%seed_dim,                              &
                           grid%DX,grid%DY,grid%vertstruc_spp_lsm,             &
                           grid%pattern_spp_lsm,                               & 
                           grid%stddev_cutoff_spp_lsm,grid%gridpt_stddev_spp_lsm, & 
                           grid%VERTSTRUCC,grid%VERTSTRUCS,grid%VERTAMPT       )
      ENDIF !spp_lsm

! calculate_phy_tend

BENCH_START(cal_phy_tend)
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )

      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call calculate_phy_tend' )
        CALL calculate_phy_tend (config_flags,grid%c1h,grid%c2h,       &
                     grid%mut,grid%muu,grid%muv,pi_phy,                &
                     grid%rthraten,                                    &
                     grid%rublten,grid%rvblten,grid%rthblten,          &
                     grid%rqvblten,grid%rqcblten,grid%rqiblten,        &
                     grid%rucuten,grid%rvcuten,grid%rthcuten,          &
                     grid%rqvcuten,grid%rqccuten,grid%rqrcuten,        &
                     grid%rqicuten,grid%rqscuten,                      &
                     grid%rushten,grid%rvshten,grid%rthshten,          &
                     grid%rqvshten,grid%rqcshten,grid%rqrshten,        &
                     grid%rqishten,grid%rqsshten,grid%rqgshten,        &
                     grid%RUNDGDTEN,grid%RVNDGDTEN,grid%RTHNDGDTEN,grid%RQVNDGDTEN, &
                     grid%RMUNDGDTEN,                                  &
                     scalar, scalar_tend, num_scalar,                  &
                     tracer, tracer_tend, num_tracer,                  &
                     ids,ide, jds,jde, kds,kde,                        &
                     ims,ime, jms,jme, kms,kme,                        &
                     grid%i_start(ij), min(grid%i_end(ij),ide-1),      &
                     grid%j_start(ij), min(grid%j_end(ij),jde-1),      &
                     k_start    , min(k_end,kde-1)                     )

      ENDDO
      !$OMP END PARALLEL DO
BENCH_END(cal_phy_tend)

! tke diffusion

       IF(config_flags%diff_opt .eq. 2 .OR. config_flags%diff_opt .eq. 1 .or. &
           ABS(config_flags%pbl3d_opt) > 0) THEN

BENCH_START(comp_diff_metrics_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL wrf_debug ( 200 , ' call compute_diff_metrics ' )
           CALL compute_diff_metrics ( config_flags, grid%ph_2, grid%phb, grid%z, grid%rdz, grid%rdzw, &
                                       grid%zx, grid%zy, grid%rdx, grid%rdy,                      &
                                       ids, ide, jds, jde, kds, kde,          &
                                       ims, ime, jms, jme, kms, kme,          &
                                       grid%i_start(ij), grid%i_end(ij),      &
                                       grid%j_start(ij), grid%j_end(ij),      &
                                       k_start    , k_end                    )
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(comp_diff_metrics_tim)

!JDM Perturbations


IF ( ( config_flags%m_pblh_opt .EQ. 1 ) .AND. ( config_flags%km_opt .EQ. 4 ) ) THEN   

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL force_down_meso_pblh( grid%m_pblh, grid%pblh,            &   
                                      ids, ide, jds, jde, kds, kde,      &
                                      ims, ime, jms, jme, kms, kme,      &
                                      grid%i_start(ij), grid%i_end(ij),  &
                                      grid%j_start(ij), grid%j_end(ij),  &
                                      k_start    , k_end                  )


        ENDDO
         !$OMP END PARALLEL DO

ENDIF


IF (  config_flags%scpm_jdm_opt .GT. 0 ) THEN   

   IF (  config_flags%scpm_jdm_opt .EQ. 1 ) THEN   
            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , grid%num_tiles

              CALL calc_scpm_jdm_t ( config_flags%scpm_jdm_opt,                  &
                                     config_flags%m_pblh_opt,                    &
                                     grid%prttms, grid%prtdt, grid%prtnk,        &
                                     grid%prtz, grid%prtseed, grid%pert_t,       &  
                                     grid%m_pblh,                                & 
                                     grid%t_2, grid%u_2, grid%v_2, grid%rdz,     &
                                     grid%dx, grid%dt,                           & 
                                     ids, ide, jds, jde, kds, kde,               &
                                     ims, ime, jms, jme, kms, kme,               &
                                     grid%i_start(ij), grid%i_end(ij),           &
                                     grid%j_start(ij), grid%j_end(ij),           &
                                     k_start    , k_end                          )

            ENDDO
            !$OMP END PARALLEL DO

   ENDIF

ENDIF


#ifdef DM_PARALLEL
#  include "HALO_EM_TKE_C.inc"
#  include "PERIOD_BDY_EM_A1.inc"
#endif

BENCH_START(tke_diff_bc_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call bc for diffusion_metrics ' )
           CALL set_physical_bc3d( grid%rdzw , 'w', config_flags,           &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%rdz , 'w', config_flags,            &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%z , 'w', config_flags,              &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%zx , 'e', config_flags,             &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc3d( grid%zy , 'f', config_flags,             &
                                   ids, ide, jds, jde, kds, kde,       &
                                   ims, ime, jms, jme, kms, kme,       &
                                   ips, ipe, jps, jpe, kps, kpe,       &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij),   &
                                   k_start    , k_end                 )
           CALL set_physical_bc2d( grid%ustm, 't', config_flags,            &
                                   ids, ide, jds, jde,                 &
                                   ims, ime, jms, jme,                 &
                                   ips, ipe, jps, jpe,                 &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij)   )
          CALL set_physical_bc2d( grid%ust, 't', config_flags,              &
                                   ids, ide, jds, jde,                 &
                                   ims, ime, jms, jme,                 &
                                   ips, ipe, jps, jpe,                 &
                                   grid%i_start(ij), grid%i_end(ij),   &
                                   grid%j_start(ij), grid%j_end(ij)   )
 
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(tke_diff_bc_tim)

         END IF

         IF (config_flags%diff_opt .eq. 2 .OR. config_flags%diff_opt .eq. 1) THEN

BENCH_START(deform_div_tim)

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call cal_deform_and_div' )
           CALL cal_deform_and_div ( config_flags,grid%u_2,grid%v_2,grid%w_2,grid%div,  &
                                     grid%defor11,grid%defor22,grid%defor33,            &
                                     grid%defor12,grid%defor13,grid%defor23,            &
                                     nba_rij, n_nba_rij,                                & !JDM
                                     grid%u_base, grid%v_base,grid%msfux,grid%msfuy,    &
                                     grid%msfvx,grid%msfvy,grid%msftx,grid%msfty,       &
                                     grid%rdx, grid%rdy, grid%dn, grid%dnw, grid%rdz,   &
                                     grid%rdzw,grid%fnm,grid%fnp,grid%cf1,grid%cf2,     &
                                     grid%cf3,grid%zx,grid%zy,            &
                                     ids, ide, jds, jde, kds, kde,        &
                                     ims, ime, jms, jme, kms, kme,        &
                                     grid%i_start(ij), grid%i_end(ij),    &
                                     grid%j_start(ij), grid%j_end(ij),    &
                                     k_start    , k_end                  )
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(deform_div_tim)

! Updraft helicity between output times

#ifdef DM_PARALLEL
#     include "HALO_EM_HELICITY.inc"
#endif

       IF ( ( config_flags%nwp_diagnostics .eq. 1 ) .OR. &
            ( ( config_flags%afwa_diag_opt .eq. 1 ) .AND. ( config_flags%afwa_severe_opt .EQ. 1 ) ) ) THEN
BENCH_START(helicity_tim)

       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

          CALL wrf_debug ( 200 , ' call cal_helicity' )
          CALL cal_helicity ( config_flags,grid%u_2,grid%v_2,grid%w_2,  &
                              grid%uh,                             &
                              grid%up_heli_max,                    &
                              grid%ph_2,grid%phb,                  &
                              grid%msfux,grid%msfuy,               &
                              grid%msfvx,grid%msfvy,               &
                              grid%ht,                             &
                              grid%rdx, grid%rdy, grid%dn, grid%dnw, grid%rdz, grid%rdzw,   &
                              grid%fnm,grid%fnp,grid%cf1,grid%cf2,grid%cf3,grid%zx,grid%zy, &
                              ids, ide, jds, jde, kds, kde,        &
                              ims, ime, jms, jme, kms, kme,        &
                              grid%i_start(ij), grid%i_end(ij),    &
                              grid%j_start(ij), grid%j_end(ij),    &
                              k_start    , k_end                  )
       ENDDO
       !$OMP END PARALLEL DO
BENCH_END(helicity_tim)
       ENDIF

#ifdef DM_PARALLEL
#     include "HALO_EM_TKE_D.inc"
#endif

! calculate tke, kmh, and kmv

BENCH_START(calc_tke_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

!XZ
         IF( config_flags%diff_opt .eq. 2 .and. config_flags%km_opt .eq. 5 ) THEN
           CALL nonlocal_flux(config_flags,grid%nlflux,grid%gamu,grid%gamv,      &
                              grid%pblh,grid%kpbl,grid%dx,grid%dy,grid%dt,       &
                              grid%ust,grid%hfx,grid%qfx,grid%br,ep_1,ep_2,      & 
                              karman,grid%u_phy,grid%v_phy,th_phy,grid%rho,      & 
                              moist,num_moist,                                   &
                              grid%msftx,grid%msfty,grid%rdzw,                   &
                              grid%u10,grid%v10,grid%wspd,                       &
                              ids, ide, jds, jde, kds, kde,                      &
                              ims, ime, jms, jme, kms, kme,                      &
                              grid%i_start(ij), grid%i_end(ij),                  &
                              grid%j_start(ij), grid%j_end(ij),                  &
                              k_start, k_end  )
 
           CALL free_atmos_length(config_flags,grid%dx,grid%dy,grid%rdzw,        &
                                  grid%rdz,grid%tke_2,th_phy,grid%elmin,         &
                                  grid%hfx,grid%qfx,moist,num_moist,             &
                                  ids, ide, jds, jde, kds, kde,                  &
                                  ims, ime, jms, jme, kms, kme,                  &
                                  grid%i_start(ij),grid%i_end(ij),               &
                                  grid%j_start(ij), grid%j_end(ij),              &
                                  k_start, k_end  )
  
           CALL meso_length_scale(config_flags,grid%dx,grid%dy,grid%rdzw,        &
                                  grid%rdz,grid%tke_2,p8w,t8w,th_phy,            &
                                  grid%dlk,grid%pblh,grid%elmin,                 &
                                  grid%rmol,grid%rho,grid%hfx,grid%qfx,          &
                                  moist,num_moist,                               &
                                  ids, ide, jds, jde, kds, kde,                  &
                                  ims, ime, jms, jme, kms, kme,                  &
                                  grid%i_start(ij), grid%i_end(ij),              &
                                  grid%j_start(ij), grid%j_end(ij),              &
                                  k_start, k_end  ) 
         ENDIF
!!

           CALL wrf_debug ( 200 , ' call calculate_km_kh' )
           CALL calculate_km_kh( config_flags,grid%dt,grid%dampcoef,grid%zdamp,         &
                                 config_flags%damp_opt,                                 &
                                 grid%xkmh,grid%xkmv,grid%xkhh,grid%xkhv,grid%bn2,      &
                                 grid%khdif,grid%kvdif,grid%div,                        &
                                 grid%defor11,grid%defor22,grid%defor33,grid%defor12,   &
                                 grid%defor13,grid%defor23,                             &
                                 grid%tke_2,p8w,t8w,th_phy,                             &
                                 t_phy,p_phy,moist,grid%dn,grid%dnw,                    &
                                 grid%dx,grid%dy,grid%rdz,grid%rdzw,                    &
                                 config_flags%mix_isotropic,num_moist,                  &
                                 grid%cf1, grid%cf2, grid%cf3, grid%warm_rain,          &
                                 grid%mix_upper_bound,                                  &
                                 grid%msftx, grid%msfty,                                &
                                 grid%zx, grid%zy,                                      &
                                 grid%pblh, grid%dlk, grid%xkmv_meso,                   &  
                                 ids,ide, jds,jde, kds,kde,                             &
                                 ims,ime, jms,jme, kms,kme,                             &
                                 grid%i_start(ij), grid%i_end(ij),                      &
                                 grid%j_start(ij), grid%j_end(ij),                      &
                                 k_start    , k_end                          )
         ENDDO
       !$OMP END PARALLEL DO
BENCH_END(calc_tke_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_TKE_E.inc"
#endif

       ENDIF

#ifdef DM_PARALLEL
#      include "PERIOD_BDY_EM_PHY_BC.inc"
       IF ( config_flags%grid_fdda .eq. 1) THEN
#      include "PERIOD_BDY_EM_FDDA_BC.inc"
       ENDIF
#      include "PERIOD_BDY_EM_CHEM.inc"
#endif

BENCH_START(phy_bc_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call phy_bc' )
         CALL phy_bc (config_flags,grid%div,grid%defor11,grid%defor22,grid%defor33,            &
                      grid%defor12,grid%defor13,grid%defor23,      &
                      grid%xkmh,grid%xkmv,grid%xkhh,grid%xkhv,     &
                      grid%tke_2,grid%rho,                         &
                      grid%rublten, grid%rvblten,                  &
                      grid%rucuten, grid%rvcuten,                  &
                      grid%rushten, grid%rvshten,                  &
                      grid%gamu, grid%gamv, grid%xkmv_meso,        &  ! XZ
                      ids, ide, jds, jde, kds, kde,                &
                      ims, ime, jms, jme, kms, kme,                &
                      ips, ipe, jps, jpe, kps, kpe,                &
                      grid%i_start(ij), grid%i_end(ij),            &
                      grid%j_start(ij), grid%j_end(ij),            &
                      k_start    , k_end                           )
       ENDDO
       !$OMP END PARALLEL DO
BENCH_END(phy_bc_tim)

!JDM
IF ( ( config_flags%sfs_opt .GT. 0 ) .AND. ( config_flags%diff_opt .eq. 2 ) ) THEN

 CALL sfs_driver( grid, config_flags,     &
                  nba_mij, n_nba_mij,     & 
                  nba_rij, n_nba_rij      ) 

ENDIF

#ifdef DM_PARALLEL
!-----------------------------------------------------------------------
!
! MPP for some physics tendency, km, kh, deformation, and divergence
!
!                                                         * * * * * * *
!                                            * * * * *    * * * * * * *
!               *                     *      * * * * *    * * * * * * *
!             * + *      * + *        +      * * + * *    * * * + * * *
!               *                     *      * * * * *    * * * * * * *
!                                            * * * * *    * * * * * * *
!                                                         * * * * * * *
!
! (for PBL)
! rublten                  x
! rvblten                             x
!
! (for Cumulus)
! rucuten                  x
! rvcuten                             x
!
! (for Shallow Cumulus)
! rushten                  x
! rvshten                             x
!
! (for FDDA)
! rundgdten     x
! rvndgdten     x
!
! (for TKE3)
! tke_2                                          x               
! (for TKE5)
! tke_2                                                         x
!
! (for diff_opt >= 1)
! defor11                  x
! defor22                             x
! defor12       x
! defor13                  x
! defor23                             x
! div           x
! xkmv          x
! xkmh          x
! xkhv          x
! xkhh          x
! tke           x
!
!-----------------------------------------------------------------------
       IF ( config_flags%bl_pbl_physics .ge. 1 ) THEN
#      include "HALO_EM_PHYS_PBL.inc"
       ENDIF
       IF ( config_flags%shcu_physics .gt. 1 ) THEN
#      include "HALO_EM_PHYS_SHCU.inc"
       ENDIF
       IF ( config_flags%cu_physics == SASSCHEME      .or.   &
            config_flags%cu_physics == TIEDTKESCHEME  .or.   &
            config_flags%cu_physics == NTIEDTKESCHEME .or.   &
            config_flags%cu_physics == MSKFSCHEME     .or.   &
            config_flags%cu_physics == CAMZMSCHEME    .or.   &
            config_flags%cu_physics == SCALESASSCHEME .or.   &
            config_flags%cu_physics == NSASSCHEME     .or.   &
            config_flags%cu_physics == KSASSCHEME ) THEN
#      include "HALO_EM_PHYS_CU.inc"
       ENDIF
       IF ( config_flags%grid_fdda .ge. 1) THEN
#      include "HALO_EM_FDDA.inc"
       ENDIF
       IF ( config_flags%diff_opt .ge. 1 ) THEN
#      include "HALO_EM_PHYS_DIFFUSION.inc"
       ENDIF

       IF      ( config_flags%h_sca_adv_order <= 4 ) THEN
#       include "HALO_EM_TKE_3.inc"
       ELSE IF ( config_flags%h_sca_adv_order <= 6 ) THEN
#       include "HALO_EM_TKE_5.inc"
       ELSE
         WRITE(wrf_err_message,*)'solve_em: invalid h_sca_adv_order = ',config_flags%h_sca_adv_order
         CALL wrf_error_fatal(TRIM(wrf_err_message))
       ENDIF

       IF ( (config_flags%bl_pbl_physics .eq. 0) .and. abs(config_flags%pbl3d_opt) .gt. 0 .and. (config_flags%pbl3d_prog .gt. 0) ) THEN
         IF      ( config_flags%h_pbl3d_adv_order <= 4 ) THEN
#         include "HALO_EM_PBL3D_3.inc"
         ELSE IF ( config_flags%h_pbl3d_adv_order <= 6 ) THEN
#         include "HALO_EM_PBL3D_5.inc"
         ELSE
           WRITE(wrf_err_message,*)'solve_em: invalid h_pbl3d_adv_order = ',config_flags%h_pbl3d_adv_order
           CALL wrf_error_fatal(TRIM(wrf_err_message))
         ENDIF
       END IF
#endif

BENCH_START(update_phy_ten_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call update_phy_ten' )
#if ( WRF_DFI_RADAR == 1 )
         if (config_flags%cu_physics .gt. 0) then
           i_start = grid%i_start(ij)
           i_end   = min( grid%i_end(ij),ide-1 )
           j_start = grid%j_start(ij)
           j_end   = min( grid%j_end(ij),jde-1 )
           if (grid%dfi_stage == DFI_FWD ) &
                 CALL wrf_debug ( 200 , ' Zero out cu_physics' )
           DO j = j_start, j_end
           DO k = k_start, min( k_end,kde-1 ) - 1
           DO i = i_start, i_end
             if (grid%dfi_stage ==DFI_FWD  &
             .and. grid%dfi_tten_rad(i,k,j) >= 1.0e-7 .and.  &
                   grid%dfi_tten_rad(i,k,j) <= 10.) then
! zero out cu-param temp tendency
                grid%rthcuten(i,k,j) = 0.0
             endif
           ENDDO
           ENDDO
           ENDDO
         ENDIF
#endif
         CALL update_phy_ten(ph_tendf,t_tendf, ru_tendf, rv_tendf,moist_tend ,&
                           scalar_tend, mu_tendf,                           &
                           grid%rthraten,grid%rthblten,grid%rthcuten,grid%rthshten, &
                           grid%rublten,grid%rucuten,grid%rushten,          &
                           grid%rvblten,grid%rvcuten,grid%rvshten,          &
                           grid%rqvblten,grid%rqcblten,grid%rqiblten,       &
                           grid%rqniblten,                                  & !CAMUWPBL scheme
                           grid%rqvcuten,grid%rqccuten,grid%rqrcuten,       &
                           grid%rqicuten,grid%rqscuten,                     &
                           grid%rqcncuten,grid%rqincuten,                   & !BSINGH - Added two CU tends
                           grid%rqvshten,grid%rqcshten,grid%rqrshten,       &
                           grid%rqishten,grid%rqsshten,grid%rqgshten,       &
                           grid%rqcnshten,grid%rqinshten,                   &!BSINGH - Added two SHCU tends
                           grid%RUNDGDTEN,                                  &
                           grid%RVNDGDTEN,grid%RTHNDGDTEN,grid%RPHNDGDTEN,  &
                           grid%RQVNDGDTEN,grid%RMUNDGDTEN,                 &
                           grid%rthfrten,grid%rqvfrten,                     &  ! fire
                           num_moist,num_scalar,config_flags,rk_step,       &
                           grid%adv_moist_cond,                             &
                           ids, ide, jds, jde, kds, kde,                    &
                           ims, ime, jms, jme, kms, kme,                    &
                           grid%i_start(ij), grid%i_end(ij),                &
                           grid%j_start(ij), grid%j_end(ij),                &
                           k_start, k_end                               )

       END DO
       !$OMP END PARALLEL DO
BENCH_END(update_phy_ten_tim)

      IF (grid%skebs_on==1) then
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles
               CALL wrf_debug ( 200 , ' call update_stoch_ten' )
                CALL update_stoch_ten(ru_tendf, rv_tendf, t_tendf,&
                               grid%ru_tendf_stoch,                          &
                               grid%rv_tendf_stoch,                          &
                               grid%rt_tendf_stoch,                          &
                               grid%mu_2, grid%mub, grid%c1h, grid%c2h,      &
                               ids, ide, jds, jde, kds, kde,                 &
                               ims, ime, jms, jme, kms, kme,                 &
                               grid%i_start(ij), grid%i_end(ij),             &
                               grid%j_start(ij), grid%j_end(ij),             &
                               k_start, k_end,                               &
                               grid%num_stoch_levels,grid%num_stoch_levels   )

           ENDDO
           !$OMP END PARALLEL DO
      ENDIF !skebs_on

      IF (grid%sppt_on==1) then
          !$OMP PARALLEL DO   &
          !$OMP PRIVATE ( ij )
          DO ij = 1 , grid%num_tiles
                 call perturb_physics_tend(grid%gridpt_stddev_sppt,          &
                        grid%stddev_cutoff_sppt,grid%rstoch,                 &
                        ru_tendf,rv_tendf,t_tendf,                           &
                        moist_tend(ims,kms,jms,p_qv),                        &
                        ids, ide, jds, jde, kds, kde,                        &
                        ims, ime, jms, jme, kms, kme,                        &
                        grid%i_start(ij), grid%i_end(ij),                    &
                        grid%j_start(ij), grid%j_end(ij),                    &
                        k_start, k_end,                                      &
                        grid%num_stoch_levels,grid%num_stoch_levels          )
           ENDDO
          !$OMP END PARALLEL DO
  ENDIF

#ifdef PLANET
       ! do rayleigh (and zonal-average newtonian) damping during
       ! first iteration of RK loop only

=============================================================
INTENTIONALLY PLACED HERE TO BREAK COMPILE
1. The DAMPTOP routine needs to have 1d column arrays c1 and c2
to correctly use the column pressures mut, muu, and muv.
2. The DAMPTOP routine needs to have the grid%t_2 evaluated
to see if perturbation of moist potential temperature is OK.
=============================================================

       IF ( (config_flags%damp_opt == 101) .OR. &
            (config_flags%damp_opt == 103)      ) THEN
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL damptop( grid%u_2, grid%v_2, grid%t_2, &
                         grid%mut, grid%muu, grid%muv, &
                         pi_phy,                                &
                         t_tendf, ru_tendf, rv_tendf, P2SI,     &
                         ids, ide, jds, jde, kds, kde,          &
                         ims, ime, jms, jme, kms, kme,          &
                         grid%i_start(ij), grid%i_end(ij),      &
                         grid%j_start(ij), grid%j_end(ij),      &
                         k_start, k_end                         )
         END DO
         !$OMP END PARALLEL DO
       END IF
#endif

       IF( config_flags%diff_opt .eq. 2 .and. ( config_flags%km_opt .eq. 2 .or. &
           config_flags%km_opt .eq. 5 ) ) THEN  ! XZ

BENCH_START(tke_rhs_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL tke_rhs  ( tke_tend,grid%bn2,                           &
                         config_flags,grid%defor11,grid%defor22,      &
                         grid%defor33,                                &
                         grid%defor12,grid%defor13,grid%defor23,      &
                         grid%u_2,grid%v_2,grid%w_2,grid%div,         &
                         grid%tke_2,grid%mut,grid%c1h,grid%c2h,       &
                         th_phy,p_phy,p8w,t8w,grid%z,grid%fnm,        & 
                         grid%fnp,grid%cf1,grid%cf2,grid%cf3,         &     
                         grid%msftx,grid%msfty,grid%xkmh,             &
                         grid%xkmv,grid%xkhv,grid%rdx,grid%rdy,       &
                         grid%dx,grid%dy,grid%dt,grid%zx,grid%zy,     &
                         grid%rdz,grid%rdzw,grid%dn,                  &
                         grid%dnw,config_flags%mix_isotropic,         &
                         grid%hfx, grid%qfx, moist(ims,kms,jms,P_QV), &
                         grid%ustm, grid%rho,                         &
                         grid%l_diss, grid%nlflux,                    & !XZ
                         grid%pblh, grid%dlk,                         & !XZ 
                         ids, ide, jds, jde, kds, kde,                &
                         ims, ime, jms, jme, kms, kme,                &
                         grid%i_start(ij), grid%i_end(ij),            &
                         grid%j_start(ij), grid%j_end(ij),            &
                         k_start    , k_end                           )

         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(tke_rhs_tim)

       ENDIF

       IF ( grid%obs_nudge_opt .EQ. 1 .AND. grid%xtime <= grid%fdda_end ) THEN
# ifdef DM_PARALLEL
#       include "HALO_OBS_NUDGE.inc"
#endif
!***********************************************************************
! This section for obs nudging
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )

         DO ij = 1 , grid%num_tiles

           CALL fddaobs_driver (grid%grid_id, model_config_rec%grid_id, &
                   model_config_rec%parent_id, config_flags%restart,    &
                   config_flags,                                        &
                   grid%obs_nudge_opt,                                  &
                   grid%obs_ipf_errob,                                  &
                   grid%obs_ipf_nudob,                                  &
                   grid%fdda_start,                                     &
                   grid%fdda_end,                                       &
                   grid%obs_nudge_wind,                                 &
                   grid%obs_nudge_temp,                                 &
                   grid%obs_nudge_mois,                                 &
                   grid%obs_nudge_pstr,                                 &
                   grid%obs_coef_wind,                                  &
                   grid%obs_coef_temp,                                  &
                   grid%obs_coef_mois,                                  &
                   grid%obs_coef_pstr,                                  &             
                   grid%obs_rinxy,                                      &
                   grid%obs_rinsig,                                     &
                   grid%obs_npfi,                                       &
                   grid%obs_ionf,                                       &
                   grid%obs_prt_max,                                    &
                   grid%obs_prt_freq,                                   &
                   grid%obs_idynin,                                     &
                   grid%obs_dtramp,                                     &
                   grid%parent_grid_ratio,                              &
                   grid%max_dom, grid%itimestep,                        &
                   grid%xtime,                                          &
                   grid%dt, grid%gmt, grid%julday, grid%fdob,           &
                   grid%max_obs,                                        &
                   model_config_rec%nobs_ndg_vars,                      &
                   model_config_rec%nobs_err_flds,                      &
                   grid%fdob%nstat, grid%fdob%varobs, grid%fdob%errf,   &
                   grid%dx, grid%KPBL,grid%HT,                          &
                   grid%mut, grid%muu, grid%muv, grid%c1h, grid%c2h,    &
                   grid%msftx, grid%msfty, grid%msfux, grid%msfuy, grid%msfvx, grid%msfvy, &
                   p_phy, t_tendf, t0,                                  &
                   grid%u_2, grid%v_2, grid%th_phy_m_t0,                &
                   moist(ims,kms,jms,P_QV),                             &
                   grid%pb, grid%p_top, grid%p, grid%phb, grid%ph_2,    &
                   grid%uratx, grid%vratx, grid%tratx,                  &
                   ru_tendf, rv_tendf,                                  &
                   moist_tend(ims,kms,jms,P_QV), grid%obs_savwt,        &
                   grid%regime, grid%pblh, grid%z_at_w, grid%z,         &
                   ids,ide, jds,jde, kds,kde,                           &
                   ims,ime, jms,jme, kms,kme,                           &
                   grid%i_start(ij), min(grid%i_end(ij),ide-1),         &
                   grid%j_start(ij), min(grid%j_end(ij),jde-1),         &
                   k_start    , min(k_end,kde-1)                     )
 
         ENDDO
         !$OMP END PARALLEL DO
       ENDIF  ! obs_nudge_opt .eq. 1

       IF ( config_flags%use_theta_m .EQ. 1 ) THEN
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL conv_t_tendf_to_moist ( grid%t_1 , moist_old(ims,kms,jms,P_qv) ,  &
                                        t_tendf  , moist_tend(ims,kms,jms,P_qv) , &
                                        ids, ide, jds, jde, kds, kde ,            &
                                        ims, ime, jms, jme, kms, kme ,            &
                                        grid%i_start(ij), grid%i_end(ij),         &
                                        grid%j_start(ij), grid%j_end(ij),         &
                                        k_start    , k_end                        )
         END DO 
         !$OMP END PARALLEL DO
       END IF ! use moist theta

! calculate vertical diffusion first and then horizontal
! (keep this order)

       IF(config_flags%diff_opt .eq. 2) THEN

         IF (config_flags%bl_pbl_physics .eq. 0 .AND. config_flags%pbl3d_opt .eq. 0) THEN

!XZ    implicit solver for vertical diffusion tendency
         IF( config_flags%km_opt .eq. 5 ) THEN
              !$OMP PARALLEL DO   &
              !$OMP PRIVATE ( ij )
              DO ij = 1 , grid%num_tiles
              CALL vertical_diffusion_implicit(ru_tendf, rv_tendf, rw_tendf,                              &
                                               t_tendf, tke_tend,                                         &
                                               moist_tend, num_moist,                                     &
                                               chem_tend, num_chem,                                       &
                                               scalar_tend, num_scalar,                                   &
                                               tracer_tend, num_tracer,                                   &
                                               grid%u_2, grid%v_2,grid%w_2,grid%dt,                       &
                                               grid%t_2,grid%u_base,grid%v_base,grid%t_base,grid%qv_base, &
                                               grid%mut,grid%tke_2,th_phy,config_flags,                   &
                                               moist, chem, scalar, tracer,                               &
                                               grid%xkmv, grid%xkhv, grid%xkmh, config_flags%km_opt,      &
                                               grid%fnm, grid%fnp, grid%dn, grid%dnw, grid%rdz, grid%rdzw,&
                                               grid%hfx, grid%qfx, grid%ustm, grid%rho,                   &
                                               grid%nlflux, grid%gamu, grid%gamv,                         & 
                                               grid%xkmv_meso, grid%l_diss,                               &
                                               grid%c1h, grid%c2h, grid%c1f, grid%c2f,                    &
                                               ids, ide, jds, jde, kds, kde,                              &
                                               ims, ime, jms, jme, kms, kme,                              &
                                               grid%i_start(ij), grid%i_end(ij),                          &
                                               grid%j_start(ij), grid%j_end(ij),                          &
                                               k_start, k_end      )
              ENDDO
              !$OMP END PARALLEL DO

         ELSE  ! original explicit 

BENCH_START(vert_diff_tim)
           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles

             CALL wrf_debug ( 200 , ' call vertical_diffusion_2 ' )
             CALL vertical_diffusion_2( ru_tendf, rv_tendf, rw_tendf,            &
                                      t_tendf, tke_tend,                         &
                                      moist_tend, num_moist,                      &
                                      chem_tend, num_chem,                       &
                                      scalar_tend, num_scalar,                     &
                                      tracer_tend, num_tracer,                     &
                                      grid%u_2, grid%v_2,                                  &
                                      grid%t_2,grid%u_base,grid%v_base,grid%t_base,grid%qv_base,          &
                                      grid%tke_2, th_phy, config_flags,              &
                                      grid%defor13,grid%defor23,grid%defor33,                   &
                                      nba_mij, num_nba_mij,          & !JDM
                                      grid%div, moist, chem, scalar,tracer,         &
                                      grid%xkmv, grid%xkhv, grid%xkmh, config_flags%km_opt,       & ! xkmh added Oct2013
                                      grid%fnm, grid%fnp, grid%dn, grid%dnw, grid%rdz, grid%rdzw, &
                                      grid%hfx, grid%qfx, grid%ustm, grid%rho,   &
                                      ids, ide, jds, jde, kds, kde,              &
                                      ims, ime, jms, jme, kms, kme,              &
                                      grid%i_start(ij), grid%i_end(ij),          &
                                      grid%j_start(ij), grid%j_end(ij),          &
                                      k_start, k_end                             )

           ENDDO
           !$OMP END PARALLEL DO
BENCH_END(vert_diff_tim)

         ENDIF
        ENDIF

        IF (config_flags%pbl3d_opt < 1) THEN
!
BENCH_START(hor_diff_tim)
         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles

           CALL wrf_debug ( 200 , ' call horizontal_diffusion_2' )
           CALL horizontal_diffusion_2( t_tendf, ru_tendf, rv_tendf, rw_tendf, &
                                      tke_tend,                              &
                                      moist_tend, num_moist,                  &
                                      chem_tend, num_chem,                   &
                                      scalar_tend, num_scalar,                 &
                                      tracer_tend, num_tracer,                 &
                                      grid%t_2, th_phy,                           &
                                      grid%tke_2, config_flags,              &
                                      grid%defor11, grid%defor22, grid%defor12,             &
                                      grid%defor13, grid%defor23,   &
                                      nba_mij, num_nba_mij,         & !JDM
                                      grid%div,                     &
                                      moist, chem, scalar,tracer,               &
                                      grid%msfux,grid%msfuy, grid%msfvx,grid%msfvy, grid%msftx,  &
                                      grid%msfty, grid%xkmh, grid%xkmv, grid%xkhh, config_flags%km_opt,     &
                                      grid%rdx, grid%rdy, grid%rdz, grid%rdzw,                   &
                                      grid%fnm, grid%fnp, grid%cf1, grid%cf2, grid%cf3,          &
                                      grid%zx, grid%zy, grid%dn, grid%dnw, grid%rho,             &
                                      ids, ide, jds, jde, kds, kde,          &
                                      ims, ime, jms, jme, kms, kme,          &
                                      grid%i_start(ij), grid%i_end(ij),      &
                                      grid%j_start(ij), grid%j_end(ij),      &
                                      k_start    , k_end                    )
         ENDDO
         !$OMP END PARALLEL DO
BENCH_END(hor_diff_tim)
        ENDIF
       ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PAJ/TWJ: Adding three dimensional turbulent mixing

       IF (config_flags%bl_pbl_physics .eq. 0 .and. ABS(config_flags%pbl3d_opt) .gt. 0) THEN

#ifdef DM_PARALLEL
#        include "HALO_EM_PHYS_STATE_PBL3D.inc"
#        include "PERIOD_BDY_EM_STATE_PBL3D.inc"
#endif

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL wrf_debug ( 200 , ' call phy_bc_state_pbl3d' )
           CALL phy_bc_state_pbl3d (config_flags,              &
                grid%u_2, grid%v_2, grid%w_2, grid%t_2,        &
                grid%thetav, moist(:,:,:,P_QV), grid%rho,      &
                grid%q_sq_prog_2, grid%q_sq, grid%l_master,    &
                ids, ide, jds, jde, kds, kde,                  &
                ims, ime, jms, jme, kms, kme,                  &
                ips, ipe, jps, jpe, kps, kpe,                  &
                grid%i_start(ij), grid%i_end(ij),              &
                grid%j_start(ij), grid%j_end(ij),              &
                k_start    , k_end                             )
         ENDDO
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL Init_temp_arrays_for_substep (u_2 = grid%u_2, v_2 = grid%v_2, &
                w_2 = grid%w_2, t_2 = grid%t_2, moist = moist(:,:,:,P_QV), rho = grid%rho, &
                q_sq_2 = grid%q_sq_prog_2, &
                u_tmp = grid%u_tmp, v_tmp = grid%v_tmp, w_tmp = grid%w_tmp, &
                th_tmp = grid%th_tmp, qv_tmp = grid%qv_tmp, rho_tmp = grid%rho_tmp, q_sq_tmp = grid%q_sq_tmp, &
                ru_tendf_tmp = ru_tendf_tmp, rv_tendf_tmp = rv_tendf_tmp, &
                rw_tendf_tmp = rw_tendf_tmp, t_tendf_tmp = t_tendf_tmp, &
                moist_tend_tmp = moist_tend_tmp, q_sq_tend_tmp = q_sq_tend_tmp, &
                ide = ide, jde = jde, kde = kde, ims = ims, ime = ime, jms = jms, &
                jme = jme, kms = kms, kme = kme, its = grid%i_start(ij), &
                ite = grid%i_end(ij), jts = grid%j_start(ij), jte = grid%j_end(ij), &
                kts = k_start, kte = k_end)
         END DO
         !$OMP END PARALLEL DO

           ! Sub-stepping for 3DPBL
         DO ss = 1, config_flags%pbl3d_nsteps

#ifdef DM_PARALLEL
#          include "HALO_EM_PHYS_SUBSTEP_PBL3D.inc"
#          include "PERIOD_BDY_EM_SUBSTEP_PBL3D.inc"
#endif

           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles
             CALL wrf_debug ( 200 , ' call phy_bc_state_pbl3d' )
             CALL phy_bc_state_pbl3d (config_flags,                &
                  grid%u_tmp, grid%v_tmp, grid%w_tmp, grid%th_tmp, &
                  grid%thetav, grid%qv_tmp, grid%rho_tmp,          &
                  grid%q_sq_tmp, grid%q_sq, grid%l_master,         &
                  ids, ide, jds, jde, kds, kde,                    &
                  ims, ime, jms, jme, kms, kme,                    &
                  ips, ipe, jps, jpe, kps, kpe,                    &
                  grid%i_start(ij), grid%i_end(ij),                &
                  grid%j_start(ij), grid%j_end(ij),                &
                  k_start    , k_end                               )
           ENDDO
           !$OMP END PARALLEL DO

           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles
             CALL wrf_debug (200 , ' call Calc_turb_fluxes_driver')
             CALL Calc_turb_fluxes_driver (config_flags = config_flags, &
             dz = dz8w, z_at_mass = grid%z, z_at_w = grid%z_at_w, u = grid%u_tmp, v = grid%v_tmp, w = grid%w_tmp, &
             th = th_phy, th2 = grid%th_tmp, qv = grid%qv_tmp, qc = moist(:,:,:,P_QC), &
             qr = moist(:,:,:,P_QR), qi = moist(:,:,:,P_QI), qs = moist(:,:,:,P_QS), qg = moist(:,:,:,P_QG), &
             qnc = scalar(:,:,:,P_QNC), qnr = scalar(:,:,:,P_QNR), qni = scalar(:,:,:,P_QNI), &
             qnwfa = scalar(:,:,:,P_QNWFA), qnifa = scalar(:,:,:,P_QNIFA), &
             thetav = grid%thetav, tsk = grid%tsk, t2=grid%t2, q2 = grid%q2, psfc = grid%psfc, &
             rho = grid%rho_tmp, ust = grid%ust, hfx = grid%hfx, qfx = grid%qfx, rmol = grid%rmol, &
             itimestep = grid%itimestep, ss = ss, n_tracer = num_tracer, &
             q_sq = grid%q_sq, q_sq_prog = grid%q_sq_tmp, q_sq_tend = q_sq_tend_tmp, &
             q_sq_vdiff_tend = grid%q_sq_vdiff, q_sq_hdiff_tend = grid%q_sq_hdiff, q_sq_shear_tend = grid%q_sq_shear, &
             q_sq_buoyancy_tend = grid%q_sq_buoyancy, q_sq_dissip_tend = grid%q_sq_dissip, &
             khx = grid%tf_khx, khy = grid%tf_khy, khz = grid%tf_khz, &
             l_master = grid%l_master, l_master_at_mass = grid%l_master_at_mass, l_boulac = grid%l_boulac, & ! ok
             msfux = grid%msfux, msfuy = grid%msfuy, msfvx = grid%msfvx, &
             msfvy = grid%msfvy, msftx = grid%msftx, msfty = grid%msfty, &
             dx = grid%dx, dy = grid%dy, &
             rdx = grid%rdx, rdy = grid%rdy, dn = grid%dn, dnw = grid%dnw, rdz = grid%rdz, dt = grid%dt, id = grid%id,  &
             rdzw = grid%rdzw, fnm = grid%fnm, fnp = grid%fnp, cf1 = grid%cf1,cf2 = grid%cf2,  &
             cf3 = grid%cf3, zx = grid%zx, zy = grid%zy, &
             mut = grid%mut, c1h = grid%c1h, c2h = grid%c2h, &
             turb_flux_u2 = grid%turb_flux_u2, turb_flux_v2 = grid%turb_flux_v2,& ! OK
             turb_flux_w2 = grid%turb_flux_w2, turb_flux_uv = grid%turb_flux_uv,& ! OK
             turb_flux_uw = grid%turb_flux_uw, turb_flux_vw = grid%turb_flux_vw,& ! OK
             turb_flux_utheta_v = grid%turb_flux_utheta_v, & ! OK
             turb_flux_vtheta_v = grid%turb_flux_vtheta_v, & ! OK
             turb_flux_wtheta_v = grid%turb_flux_wtheta_v, & ! OK
             turb_flux_theta2_v = grid%turb_flux_theta2_v, & ! OK
             turb_flux_uqv = grid%turb_flux_uqv, turb_flux_vqv = grid%turb_flux_vqv, turb_flux_wqv = grid%turb_flux_wqv, &
             turb_flux_uqc = grid%turb_flux_uqc, turb_flux_vqc = grid%turb_flux_vqc, turb_flux_wqc = grid%turb_flux_wqc, &
             turb_flux_uqr = grid%turb_flux_uqr, turb_flux_vqr = grid%turb_flux_vqr, turb_flux_wqr = grid%turb_flux_wqr, &
             turb_flux_uqi = grid%turb_flux_uqi, turb_flux_vqi = grid%turb_flux_vqi, turb_flux_wqi = grid%turb_flux_wqi, &
             turb_flux_uqs = grid%turb_flux_uqs, turb_flux_vqs = grid%turb_flux_vqs, turb_flux_wqs = grid%turb_flux_wqs, &
             turb_flux_uqg = grid%turb_flux_uqg, turb_flux_vqg = grid%turb_flux_vqg, turb_flux_wqg = grid%turb_flux_wqg, &
             turb_flux_uqnc = grid%turb_flux_uqnc, turb_flux_vqnc = grid%turb_flux_vqnc, turb_flux_wqnc = grid%turb_flux_wqnc, &
             turb_flux_uqnr = grid%turb_flux_uqnr, turb_flux_vqnr = grid%turb_flux_vqnr, turb_flux_wqnr = grid%turb_flux_wqnr, &
             turb_flux_uqni = grid%turb_flux_uqni, turb_flux_vqni = grid%turb_flux_vqni, turb_flux_wqni = grid%turb_flux_wqni, &
             turb_flux_uqnwfa = grid%turb_flux_uqnwfa, turb_flux_vqnwfa = grid%turb_flux_vqnwfa, turb_flux_wqnwfa = grid%turb_flux_wqnwfa, &
             turb_flux_uqnifa = grid%turb_flux_uqnifa, turb_flux_vqnifa = grid%turb_flux_vqnifa, turb_flux_wqnifa = grid%turb_flux_wqnifa, &
             turb_flux_utheta = grid%turb_flux_utheta, &
             turb_flux_vtheta = grid%turb_flux_vtheta, &
             turb_flux_wtheta = grid%turb_flux_wtheta, &
             turb_flux_u2_mass = grid%turb_flux_u2_mass, turb_flux_v2_mass = grid%turb_flux_v2_mass,& ! OK
             turb_flux_w2_mass = grid%turb_flux_w2_mass, turb_flux_uv_mass = grid%turb_flux_uv_mass,& ! OK
             turb_flux_uw_mass = grid%turb_flux_uw_mass, turb_flux_vw_mass = grid%turb_flux_vw_mass,& ! OK
             turb_flux_wtheta_v_mass = grid%turb_flux_wtheta_v_mass, &
             ids = ids, ide = ide, jds = jds, jde = jde, kds = kds, kde = kde, & ! ok
             ims = ims, ime = ime, jms = jms, jme = jme, kms = kms, kme = kme, & ! ok
             its = grid%i_start(ij), ite = grid%i_end(ij), & ! ok
             jts = grid%j_start(ij), jte = grid%j_end(ij), & ! ok
             kts = k_start, kte = k_end                    & ! ok
             )
           END DO
           !$OMP END PARALLEL DO

#ifdef DM_PARALLEL
#        include "HALO_EM_PHYS_DIFFUSION_PBL3D.inc"
#        include "PERIOD_BDY_EM_PHY_BC_PBL3D.inc"
         IF ( config_flags%pbl3d_scalar_mix .eq. 1 .and. ss .eq. config_flags%pbl3d_nsteps ) THEN
#          include "HALO_EM_SCALAR_DIFFUSION_PBL3D.inc"
#          include "PERIOD_BDY_EM_SCALAR_BC_PBL3D.inc"
         END IF
#endif

           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles
             CALL wrf_debug ( 200 , ' call phy_bc_pbl3d' )
             CALL phy_bc_pbl3d (config_flags,                    &
                  grid%turb_flux_u2, grid%turb_flux_v2,          &
                  grid%turb_flux_w2, grid%turb_flux_uv,          &
                  grid%turb_flux_uw, grid%turb_flux_vw,          &
                  grid%turb_flux_utheta, grid%turb_flux_vtheta,  &
                  grid%turb_flux_wtheta, grid%turb_flux_uqv,     &
                  grid%turb_flux_vqv, grid%turb_flux_wqv,        &
                  ids, ide, jds, jde, kds, kde,                  &
                  ims, ime, jms, jme, kms, kme,                  &
                  ips, ipe, jps, jpe, kps, kpe,                  &
                  grid%i_start(ij), grid%i_end(ij),              &
                  grid%j_start(ij), grid%j_end(ij),              &
                  k_start    , k_end                             )
           ENDDO
           !$OMP END PARALLEL DO

           IF ( config_flags%pbl3d_scalar_mix .eq. 1 .and. ss .eq. config_flags%pbl3d_nsteps ) THEN
             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij )
             DO ij = 1 , grid%num_tiles
               CALL wrf_debug ( 200 , ' call scalar_bc_pbl3d' )
               CALL scalar_bc_pbl3d (config_flags,                                  &
                    grid%turb_flux_uqc, grid%turb_flux_vqc, grid%turb_flux_uqc,     &
                    grid%turb_flux_uqr, grid%turb_flux_vqr, grid%turb_flux_uqr,     &
                    grid%turb_flux_uqi, grid%turb_flux_vqi, grid%turb_flux_uqi,     &
                    grid%turb_flux_uqs, grid%turb_flux_vqs, grid%turb_flux_uqs,     &
                    grid%turb_flux_uqg, grid%turb_flux_vqg, grid%turb_flux_uqg,     &
                    grid%turb_flux_uqnc, grid%turb_flux_vqnc, grid%turb_flux_wqnc,  &
                    grid%turb_flux_uqnr, grid%turb_flux_vqnr, grid%turb_flux_wqnr,  &
                    grid%turb_flux_uqni, grid%turb_flux_vqni, grid%turb_flux_wqni,  &
                    grid%turb_flux_uqnwfa, grid%turb_flux_vqnwfa, grid%turb_flux_wqnwfa, &
                    grid%turb_flux_uqnifa, grid%turb_flux_vqnifa, grid%turb_flux_wqnifa, &
                    ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme,     &
                    ips, ipe, jps, jpe, kps, kpe, grid%i_start(ij), grid%i_end(ij), &
                    grid%j_start(ij), grid%j_end(ij), k_start, k_end                )
             ENDDO
           END IF
           !$OMP END PARALLEL DO

           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles
             CALL wrf_debug (200 , ' call Vertical_turb_mix')
             CALL Vertical_turb_mix (dth = t_tendf_tmp, du = ru_tendf_tmp, &
                  dv = rv_tendf_tmp, dw = rw_tendf_tmp, dqv = moist_tend_tmp, dqc = moist_tend(:,:,:,P_QC), &
                  dqr = moist_tend(:,:,:,P_QR), dqi = moist_tend(:,:,:,P_QI), dqs = moist_tend(:,:,:,P_QS), &
                  dqg = moist_tend(:,:,:,P_QG), dqnc = scalar_tend(:,:,:,P_QNC), dqnr = scalar_tend(:,:,:,P_QNR), &
                  dqni = scalar_tend(:,:,:,P_QNI), dqnwfa = scalar_tend(:,:,:,P_QNWFA), dqnifa = scalar_tend(:,:,:,P_QNIFA), &
                  khz = grid%tf_khz, tracer = tracer, tracer_tend = tracer_tend, &
                  pbl3d_opt = config_flags%pbl3d_opt, ss = ss, &
                  config_flags = config_flags, &
                  turb_flux_wtheta = grid%turb_flux_wtheta, turb_flux_uw = grid%turb_flux_uw, turb_flux_vw = grid%turb_flux_vw, &
                  turb_flux_w2 = grid%turb_flux_w2, &
                  turb_flux_wqv = grid%turb_flux_wqv, turb_flux_wqc = grid%turb_flux_wqc, turb_flux_wqr = grid%turb_flux_wqr, &
                  turb_flux_wqi = grid%turb_flux_wqi, turb_flux_wqs = grid%turb_flux_wqs, turb_flux_wqg = grid%turb_flux_wqg, &
                  turb_flux_wqnc = grid%turb_flux_wqnc, turb_flux_wqnr = grid%turb_flux_wqnr, turb_flux_wqni = grid%turb_flux_wqni, &
                  turb_flux_wqnwfa = grid%turb_flux_wqnwfa, turb_flux_wqnifa = grid%turb_flux_wqnifa, &
                  rdz = grid%rdz, dnw = grid%dnw, dn = grid%dn,  rho = grid%rho_tmp, &
                  fnm = grid%fnm, fnp = grid%fnp,  &
                  n_tracer = num_tracer, &
                  ids = ids, ide = ide, jds = jds, jde = jde, kds = kds, kde = kde, &
                  ims = ims, ime = ime, jms = jms, jme = jme, kms = kms, kme = kme, &
                  its = grid%i_start(ij), ite = grid%i_end(ij), &
                  jts = grid%j_start(ij), jte = grid%j_end(ij), &
                  kts = k_start, kte = k_end                    )
           END DO
           !$OMP END PARALLEL DO

           IF (config_flags%pbl3d_opt > 0) THEN
             !$OMP PARALLEL DO   &
             !$OMP PRIVATE ( ij )
             DO ij = 1 , grid%num_tiles
               CALL wrf_debug (200 , ' call Horizontal_turb_mix')
               CALL Horizontal_turb_mix (rt_tendf = t_tendf_tmp, ru_tendf = ru_tendf_tmp, rv_tendf = rv_tendf_tmp, &
                    rw_tendf = rw_tendf_tmp, moist_tendf = moist_tend_tmp, qc_tendf = moist_tend(:,:,:,P_QC), &
                    qr_tendf = moist_tend(:,:,:,P_QR), qi_tendf = moist_tend(:,:,:,P_QI), &
                    qs_tendf = moist_tend(:,:,:,P_QS), qg_tendf = moist_tend(:,:,:,P_QG), &
                    qnc_tendf = scalar_tend(:,:,:,P_QNC), qnr_tendf = scalar_tend(:,:,:,P_QNR), &
                    qni_tendf = scalar_tend(:,:,:,P_QNI), qnwfa_tendf = scalar_tend(:,:,:,P_QNWFA), qnifa_tendf = scalar_tend(:,:,:,P_QNIFA), &
                    tracer_tend = tracer_tend, tracer = tracer, n_tracer = num_tracer, khx = grid%tf_khx, khy = grid%tf_khy, &
                    ss = ss, config_flags = config_flags, &
                    turb_flux_u2 = grid%turb_flux_u2, turb_flux_v2 = grid%turb_flux_v2, &
                    turb_flux_uv = grid%turb_flux_uv, turb_flux_uw = grid%turb_flux_uw, &
                    turb_flux_vw = grid%turb_flux_vw, turb_flux_w2 = grid%turb_flux_w2, &
                    turb_flux_utheta = grid%turb_flux_utheta, turb_flux_vtheta = grid%turb_flux_vtheta, &
                    turb_flux_uqv = grid%turb_flux_uqv, turb_flux_vqv = grid%turb_flux_vqv, &
                    turb_flux_uqc = grid%turb_flux_uqc, turb_flux_vqc = grid%turb_flux_vqc, &
                    turb_flux_uqr = grid%turb_flux_uqr, turb_flux_vqr = grid%turb_flux_vqr, &
                    turb_flux_uqi = grid%turb_flux_uqi, turb_flux_vqi = grid%turb_flux_vqi, &
                    turb_flux_uqs = grid%turb_flux_uqs, turb_flux_vqs = grid%turb_flux_vqs, &
                    turb_flux_uqg = grid%turb_flux_uqg, turb_flux_vqg = grid%turb_flux_vqg, &
                    turb_flux_uqnc = grid%turb_flux_uqnc, turb_flux_vqnc = grid%turb_flux_vqnc, &
                    turb_flux_uqnr = grid%turb_flux_uqnr, turb_flux_vqnr = grid%turb_flux_vqnr, &
                    turb_flux_uqni = grid%turb_flux_uqni, turb_flux_vqni = grid%turb_flux_vqni, &
                    turb_flux_uqnwfa = grid%turb_flux_uqnwfa, turb_flux_vqnwfa = grid%turb_flux_vqnwfa, &
                    turb_flux_uqnifa = grid%turb_flux_uqnifa, turb_flux_vqnifa = grid%turb_flux_vqnifa, &
                    msfux = grid%msfux, msfuy = grid%msfuy, msfvx = grid%msfvx, msfvy = grid%msfvy, &
                    msftx = grid%msftx, msfty = grid%msfty, dx = grid%dx, dy = grid%dy, &
                    rdx = grid%rdx, rdy = grid%rdy, rdz = grid%rdz, rdzw = grid%rdzw, &
                    fnm = grid%fnm, fnp = grid%fnp, cf1 = grid%cf1, cf2 = grid%cf2, cf3 = grid%cf3,          &
                    zx = grid%zx, zy = grid%zy, dn = grid%dn, dnw = grid%dnw, rho = grid%rho_tmp, &
                    ids = ids, ide = ide, jds = jds, jde = jde, kds = kds, kde = kde, &
                    ims = ims, ime = ime, jms = jms, jme = jme, kms = kms, kme = kme, &
                    its = grid%i_start(ij), ite = grid%i_end(ij), &
                    jts = grid%j_start(ij), jte = grid%j_end(ij), &
                    kts = k_start, kte = k_end                    )
             END DO
             !$OMP END PARALLEL DO
           END IF

           !$OMP PARALLEL DO   &
           !$OMP PRIVATE ( ij )
           DO ij = 1 , grid%num_tiles
             CALL Update_wrf_tends_temp_state_and_zero_tends (ru_tendf = ru_tendf, rv_tendf = rv_tendf, &
                  rw_tendf = rw_tendf, t_tendf = t_tendf, moist_tend = moist_tend(:,:,:,P_QV), &
                  q_sq_tend = q_sq_prog_tend, &
                  ru_tendf_tmp = ru_tendf_tmp, rv_tendf_tmp = rv_tendf_tmp, &
                  rw_tendf_tmp = rw_tendf_tmp, t_tendf_tmp = t_tendf_tmp, &
                  moist_tend_tmp = moist_tend_tmp, q_sq_tend_tmp = q_sq_tend_tmp, &
                  u_tmp = grid%u_tmp, v_tmp = grid%v_tmp, &
                  w_tmp = grid%w_tmp, th_tmp = grid%th_tmp, qv_tmp = grid%qv_tmp, rho_tmp = grid%rho_tmp, &
                  q_sq_tmp = grid%q_sq_tmp, &
                  dt = grid%dt, c1h = grid%c1h, c2h = grid%c2h, c1f = grid%c1f, &
                  c2f = grid%c2f, muu = grid%muu, muv = grid%muv, mut = grid%mut, &
                  ide = ide, jde = jde, kde = kde, ims = ims, ime = ime, jms = jms, &
                  jme = jme, kms = kms, kme = kme, its = grid%i_start(ij), &
                  ite = grid%i_end(ij), jts = grid%j_start(ij), jte = grid%j_end(ij), &
                  kts = k_start, kte = k_end)
           END DO
           !$OMP END PARALLEL DO


         END DO ! Sub-stepping for 3DPBL

         !$OMP PARALLEL DO   &
         !$OMP PRIVATE ( ij )
         DO ij = 1 , grid%num_tiles
           CALL Deallocate_temp_tends (ru_tendf_tmp = ru_tendf_tmp, &
                rv_tendf_tmp = rv_tendf_tmp, rw_tendf_tmp = rw_tendf_tmp, &
                t_tendf_tmp = t_tendf_tmp, moist_tend_tmp = moist_tend_tmp, &
                q_sq_tend_tmp = q_sq_tend_tmp, &
                its = grid%i_start(ij), ite = grid%i_end(ij), jts = grid%j_start(ij), &
                jte = grid%j_end(ij), kts = k_start, kte = k_end)
         END DO
         !$OMP END PARALLEL DO

       END IF

! 
!***********************************************************************

  END SUBROUTINE first_rk_step_part2

END MODULE module_first_rk_step_part2

