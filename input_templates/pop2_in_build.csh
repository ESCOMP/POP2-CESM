#!/bin/csh -f

set verbose

#==========================================================================
# This script constructs a CCSM pop2_in file. It is invoked from within 
# pop2.buildnml_prestage.csh, the CCSM pop2 build script
#
# Except for $rundir, variables in UPPER CASE are set externally, variables
# in lower case are set internally
#==========================================================================


#--------------------------------------------------------------------------
#  Define directory locations
#    $rundir, $CASE, $LID, $INPUT, and $POP2_NMLFILE are environment variables 
#       defined in pop2.buildnml_prestage.csh
#--------------------------------------------------------------------------

set output_L = $rundir
set output_r = $rundir/rest/$CASE.pop.r
set output_h = $rundir/hist/$CASE.pop.h
set output_d = $rundir/hist/$CASE.pop.d

set log_filename =  ocn.log.$LID
set pop2_pointer = $rundir/rpointer.ocn


#--------------------------------------------------------------------------
#  domain_nml
#     $NPROCS_CLINIC and $NPROCS_TROPIC are defined in 
#      pop2.buildnml_prestage.csh
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&domain_nml
  nprocs_clinic = $NPROCS_CLINIC
  nprocs_tropic = $NPROCS_TROPIC
  clinic_distribution_type = 'balanced'
  tropic_distribution_type = 'cartesian'
  ew_boundary_type = 'cyclic'
  ns_boundary_type = 'closed'
/

EOF
 
#--------------------------------------------------------------------------
#  io_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&io_nml
  num_iotasks        = 1 
  lredirect_stdout   = .true. 
  log_filename       = '$output_L/$log_filename'
  luse_pointer_files = .true.
  pointer_filename   = '$pop2_pointer' 
/

EOF

#--------------------------------------------------------------------------
#  time_manager_nml
#    $IYEAR0, $IMONTH0, and $IDAY0 are defined in pop2.buildnml_prestage.csh
#--------------------------------------------------------------------------

if      ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
  setenv DT_COUNT 12
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
  setenv DT_COUNT 23
endif

cat >> $POP2_NMLFILE << EOF
&time_manager_nml
  runid             = '$CASE'
  time_mix_opt      = 'avgfit'
  time_mix_freq     = 17
  dt_option         = 'steps_per_day'
  dt_count          = $DT_COUNT
  impcor            = .true.
  laccel            = .false.
  accel_file        = '$INPUT/depth_accel'
  dtuxcel           = 1.0 
  allow_leapyear    = .false.
  iyear0            = $IYEAR0
  imonth0           = $IMONTH0
  iday0             = $IDAY0
  ihour0            = 0
  iminute0          = 0
  isecond0          = 0
  date_separator    = '-'
  stop_option       = 'nyear'
  stop_count        =  1000
  fit_freq          = 1
/

EOF

#--------------------------------------------------------------------------
#  grid_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&grid_nml
   horiz_grid_opt       = 'file'
   horiz_grid_file      = '$INPUT/horiz_grid'
   vert_grid_opt        = 'file'
   vert_grid_file       = '$INPUT/vert_grid'
   topography_opt       = 'file'
   topography_file      = '$INPUT/topography'
   partial_bottom_cells = .false.
   bottom_cell_file     = 'unknown_bottom_cell'
   topo_smooth          = .false.
   flat_bottom          = .false.
   lremove_points       = .false.
   region_mask_file     = '$INPUT/region_mask'
   region_info_file     = '$INPUT/region_ids'
   sfc_layer_opt        = 'varthick'
/

EOF

#--------------------------------------------------------------------------
#  init_ts_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&init_ts_nml
   init_ts_option    = '$runtype'
   init_ts_suboption = 'null'
   init_ts_file      = '$INPUT/ts'
   init_ts_file_fmt  = 'bin'
/

EOF

#--------------------------------------------------------------------------
#  diagnostics_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&diagnostics_nml
   diag_global_freq_opt   = 'nmonth'
   diag_global_freq       = 1 
   diag_cfl_freq_opt      = 'nmonth'
   diag_cfl_freq          = 1
   diag_transp_freq_opt   = 'nmonth'
   diag_transp_freq       = 1
   diag_transport_file    = '$INPUT/transport_contents'
   diag_outfile           = '${output_d}d'
   diag_transport_outfile = '${output_d}t'
   cfl_all_levels         = .false.
   diag_all_levels        = .false.
   diag_velocity_outfile  = '${output_d}v'
/

EOF

#--------------------------------------------------------------------------
#  budget_diagnostics_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&budget_diagnostics_nml
   ldiag_global_tracer_budgets = .true.
/

EOF

#--------------------------------------------------------------------------
#  restart_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&restart_nml
   restart_freq_opt    = 'nyear' 
   restart_freq        = 100000
   restart_outfile     = '$output_r'
   restart_fmt         = 'bin'
   leven_odd_on        = .false. 
   even_odd_freq       = 100000
   pressure_correction = .false.
/

EOF

#--------------------------------------------------------------------------
#  tavg_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&tavg_nml
   tavg_freq_opt  = 'nmonth'
   tavg_freq      = 1
   tavg_start_opt = 'nstep'
   tavg_start     = 0
   tavg_infile    = '${output_h}restart.end'
   tavg_fmt_in    = 'nc'
   tavg_outfile   = '$output_h'
   tavg_fmt_out   = 'nc'
   tavg_contents  = '$INPUT/tavg_contents'
/

EOF

#--------------------------------------------------------------------------
#  history_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&history_nml
   history_freq_opt  = 'never'
   history_freq      = 1
   history_outfile   = '${output_h}s'
   history_fmt       = 'nc'
   history_contents  = '$INPUT/history_contents'
/

EOF

#--------------------------------------------------------------------------
#  movie_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&movie_nml
   movie_freq_opt = 'never'
   movie_freq     = 1
   movie_outfile  = '${output_h}m'
   movie_fmt      = 'nc'
   movie_contents = '$INPUT/movie_contents'
/

EOF

#--------------------------------------------------------------------------
#  solver_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&solver_nml
   solv_convrg    = 1.0e-13 
   solv_max_iters = 1000
   solv_ncheck    = 10
   lprecond       = .false.
   solv_type      = 'cg1'
   precond_file   = 'empty'
/

EOF

#--------------------------------------------------------------------------
#  vertical_mix_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&vertical_mix_nml
   vmix_choice           = 'kpp'
   aidif                 = 1.0
   implicit_vertical_mix = .true.
   convection_type       = 'diffusion'
   nconvad               = 2
   convect_diff          = 1000.0
   convect_visc          = 1000.0
   bottom_drag           = 1.0e-3
   bottom_heat_flx       = 0.0
   bottom_heat_flx_depth = 1000.0e2
/

EOF

#--------------------------------------------------------------------------
#  vmix_const_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&vmix_const_nml
   const_vvc = 0.25
   const_vdc = 0.25
/

EOF

#--------------------------------------------------------------------------
#  vmix_rich_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&vmix_rich_nml
   bckgrnd_vvc = 1.0
   bckgrnd_vdc = 0.1
   rich_mix    = 50.0
/

EOF

#--------------------------------------------------------------------------
#  vmix_kpp_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v6 ) then
 set bckgrnd_vdc1 = 0.1
 set bckgrnd_vdc2 = 0.0
else if ( ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 set bckgrnd_vdc1 = 0.1
 set bckgrnd_vdc2 = 0.0
else
 set bckgrnd_vdc1 = 0.524
 set bckgrnd_vdc2 = 0.313
endif

cat >> $POP2_NMLFILE << EOF
&vmix_kpp_nml
   bckgrnd_vdc1    = $bckgrnd_vdc1
   bckgrnd_vdc2    = $bckgrnd_vdc2
   bckgrnd_vdc_dpth= 1000.0e02 
   bckgrnd_vdc_linv= 4.5e-05
   Prandtl         = 10.0
   rich_mix        = 50.0
   lrich           = .true.
   ldbl_diff       = .true.
   lshort_wave     = .true.
   lcheckekmo      = .false.
   num_v_smooth_Ri = 1
/

EOF

#--------------------------------------------------------------------------
#  advect_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&advect_nml
   tadvect_ctype = 'upwind3'
/

EOF

#--------------------------------------------------------------------------
#  hmix_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&hmix_nml
   hmix_momentum_choice = 'anis'
   hmix_tracer_choice   = 'gent'
/

EOF

#--------------------------------------------------------------------------
#  hmix_del2u_nml 
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
  set am_del2_value = 3.0e9
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
  set am_del2_value = 0.5e8
endif

cat >> $POP2_NMLFILE << EOF
&hmix_del2u_nml
   lauto_hmix           = .false. 
   lvariable_hmix       = .false. 
   am                   = $am_del2_value
/

EOF

#--------------------------------------------------------------------------
# hmix_del2t_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
  set ah_del2_value = 1.0e7
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
  set ah_del2_value = 0.6e7
endif

cat >> $POP2_NMLFILE << EOF
&hmix_del2t_nml
   lauto_hmix           = .false. 
   lvariable_hmix       = .false. 
   ah                   = $ah_del2_value
/

EOF

#--------------------------------------------------------------------------
# hmix_del4u_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&hmix_del4u_nml
   lauto_hmix           = .false. 
   lvariable_hmix       = .false. 
   am                   = -0.6e20
/

EOF

#--------------------------------------------------------------------------
# hmix_del4t_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&hmix_del4t_nml
   lauto_hmix           = .false. 
   lvariable_hmix       = .false. 
   ah                   = -0.2e20
/

EOF

#--------------------------------------------------------------------------
#  hmix_gm_nml
#--------------------------------------------------------------------------

set kappa_isop_choice = bfre
set kappa_thic_choice = bfre

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.8e7
   set ah_bolus_value = 0.8e7
   set ah_bkg_srfbl   = 0.8e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 4.0e7
   set ah_bolus_value = 4.0e7
   set ah_bkg_srfbl   = 4.0e7
 endif
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.6e7
   set ah_bolus_value = 0.6e7
   set ah_bkg_srfbl   = 0.6e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 3.0e7
   set ah_bolus_value = 3.0e7
   set ah_bkg_srfbl   = 3.0e7
 endif
endif


cat >> $POP2_NMLFILE << EOF
&hmix_gm_nml
   kappa_isop_choice      = '$kappa_isop_choice'
   kappa_thic_choice      = '$kappa_thic_choice'
   kappa_freq_choice      = 'once_a_day'
   slope_control_choice   = 'notanh'
   kappa_depth_1          = 1.0
   kappa_depth_2          = 0.0
   kappa_depth_scale      = 150000.0
   ah                     = $ah_gm_value
   ah_bolus               = $ah_bolus_value
   use_const_ah_bkg_srfbl = .true.
   ah_bkg_srfbl           = $ah_bkg_srfbl
   ah_bkg_bottom          = 0.0
   slm_r                  = 0.3
   slm_b                  = 0.3
   diag_gm_bolus          = .true.
   transition_layer_on    = .true.
   read_n2_data           = .false.
   buoyancy_freq_filename = '$INPUT/buoyancy_freq'
   buoyancy_freq_fmt      = 'nc'
/

EOF

#--------------------------------------------------------------------------
#  hmix_aniso_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
 set hmix_alignment_choice =  grid
 set lvariable_hmix_aniso  =  .true.
 set lsmag_aniso           =  .false.
 set visc_para =  1.0
 set visc_perp =  1.0
 set c_para    =  0.0
 set c_perp    =  0.0
 set u_para    =  0.0
 set u_perp    =  0.0
 set vconst_1  =  1.0e7
 set vconst_2  = 24.5
 set vconst_3  =  0.2
 set vconst_4  =  1.0e-8
 set vconst_5  =  3
 set vconst_6  = 1.0e7
 set vconst_7  = 90.0
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 set hmix_alignment_choice =  east
 set lvariable_hmix_aniso =  .true.
 set lsmag_aniso          =  .false.
 set visc_para = 50.0e7
 set visc_perp = 50.0e7
 set c_para    =  8.0
 set c_perp    =  8.0
 set u_para    =  5.0
 set u_perp    =  5.0
 set vconst_1  =  0.6e7
 set vconst_2  =  0.5
 set vconst_3  =  0.16
 set vconst_4  =  2.e-8
 set vconst_5  =  3
 set vconst_6  =  0.6e7
 set vconst_7  =  45.0
endif


cat >> $POP2_NMLFILE << EOF
&hmix_aniso_nml
   hmix_alignment_choice = '$hmix_alignment_choice'
   lvariable_hmix_aniso  = $lvariable_hmix_aniso
   lsmag_aniso    = $lsmag_aniso
   visc_para      = $visc_para
   visc_perp      = $visc_perp
   c_para         = $c_para
   c_perp         = $c_perp
   u_para         = $u_para
   u_perp         = $u_perp
   vconst_1       = $vconst_1
   vconst_2       = $vconst_2
   vconst_3       = $vconst_3
   vconst_4       = $vconst_4
   vconst_5       = $vconst_5
   vconst_6       = $vconst_6
   vconst_7       = $vconst_7
   smag_lat       = 20.0
   smag_lat_fact  = 0.98
   smag_lat_gauss = 98.0
   var_viscosity_infile      = 'ccsm-internal'
   var_viscosity_infile_fmt  = 'bin'
   var_viscosity_outfile     = '${output_h}v'
   var_viscosity_outfile_fmt = 'nc'
/

EOF

#--------------------------------------------------------------------------
# state_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&state_nml
   state_choice     = 'mwjf'
   state_file       = 'internal'
   state_range_opt  = 'enforce'
   state_range_freq = 100000
/

EOF

#--------------------------------------------------------------------------
# baroclinic_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&baroclinic_nml
   reset_to_freezing = .false.
/

EOF

#--------------------------------------------------------------------------
#  ice_nml
#   OCN_ICE_FORCING is a CCSM environment variable; see $case/env_conf
#--------------------------------------------------------------------------

if ($OCN_ICE_FORCING =~ *inactive*) then
  set lactive_ice = .false.
else
  set lactive_ice = .true.
endif

cat >> $POP2_NMLFILE << EOF
&ice_nml
   ice_freq_opt     = 'coupled'
   ice_freq         =  100000
   kmxice           = 1
   lactive_ice      = $lactive_ice
/

EOF

#--------------------------------------------------------------------------
# pressure_grad_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&pressure_grad_nml
   lpressure_avg  = .true.
   lbouss_correct = .false.
/

EOF

#--------------------------------------------------------------------------
# topostress_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&topostress_nml
   ltopostress = .false.
   nsmooth_topo = 0
/

EOF

#--------------------------------------------------------------------------
# xdisplay_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&xdisplay_nml
   lxdisplay      = .false.
   nstep_xdisplay = 1
/

EOF

#--------------------------------------------------------------------------
# forcing_ws_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&forcing_ws_nml
   ws_data_type      = 'none'
   ws_data_inc       = 24.
   ws_interp_freq    = 'every-timestep'
   ws_interp_type    = 'linear'
   ws_interp_inc     = 72.
   ws_filename       = 'unknown-ws'
   ws_file_fmt       = 'bin'
   ws_data_renorm(1) = 10.
/

EOF

#--------------------------------------------------------------------------
#  forcing_shf_nml and forcing_sfwf_nml
#    OCN_COUPLING is a CCSM environment variable; see $case/env_conf
#--------------------------------------------------------------------------

if ($OCN_COUPLING  =~ *partial*) then
  #... forcing_shf_nml and forcing_sfwf_nml
  set formulation       = partially-coupled
  set data_type         = monthly

  #... forcing_shf_nml
  set luse_cpl_ifrac    = .true.

  #...forcing_sfwf_nml
  set ladjust_precip    = .true.
  set lms_balance       = .false.
  set lsend_precip_fact = .true.
else
  #... forcing_shf_nml and forcing_sfwf_nml
  set formulation       = restoring
  set data_type         = none

  #... forcing_shf_nml
  set luse_cpl_ifrac    = .false.

  #... forcing_sfwf_nml
  set ladjust_precip    = .false.
  set lms_balance       = .true.
  set lsend_precip_fact = .false.
endif

#--------------------------------------------------------------------------
#  forcing_shf_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&forcing_shf_nml
   shf_formulation      = '$formulation'
   shf_data_type        = '$data_type'
   shf_data_inc         = 24.
   shf_interp_freq      = 'every-timestep'
   shf_interp_type      = 'linear'
   shf_interp_inc       = 72.
   shf_restore_tau      = 30.
   shf_filename         = '$INPUT/shf'
   shf_file_fmt         = 'bin'
   shf_data_renorm(3)   = 0.94
   shf_weak_restore     = 0.
   shf_strong_restore   = 0.0
   luse_cpl_ifrac       = $luse_cpl_ifrac
   shf_strong_restore_ms= 92.64
/

EOF

#--------------------------------------------------------------------------
#  forcing_sfwf_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
 set sfwf_weak_restore = 0.092
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 set sfwf_weak_restore = 0.0115
endif

cat >> $POP2_NMLFILE << EOF
&forcing_sfwf_nml
   sfwf_formulation       = '$formulation'
   sfwf_data_type         = '$data_type'
   sfwf_data_inc          = 24.
   sfwf_interp_freq       = 'every-timestep'
   sfwf_interp_type       = 'linear'
   sfwf_interp_inc        = 72.
   sfwf_restore_tau       = 30.
   sfwf_filename          = '$INPUT/sfwf'
   sfwf_file_fmt          = 'bin'
   sfwf_data_renorm(1)    = 0.001
   sfwf_weak_restore      = $sfwf_weak_restore
   sfwf_strong_restore    = 0.0
   sfwf_strong_restore_ms = 0.6648
   ladjust_precip         = $ladjust_precip
   lms_balance            = $lms_balance
   lfw_as_salt_flx        = .true.
   lsend_precip_fact      = $lsend_precip_fact
/

EOF

#--------------------------------------------------------------------------
# forcing_pt_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&forcing_pt_interior_nml
   pt_interior_data_type         = 'none'
   pt_interior_data_inc          = 24.
   pt_interior_interp_freq       = 'every-timestep'
   pt_interior_interp_type       = 'linear'
   pt_interior_interp_inc        = 72.
   pt_interior_restore_tau       = 365.
   pt_interior_filename          = 'unknown-pt_interior'
   pt_interior_file_fmt          = 'bin'
   pt_interior_restore_max_level = 0 
   pt_interior_formulation       = 'restoring'
   pt_interior_data_renorm(1)    = 1.
   pt_interior_variable_restore  = .false.
   pt_interior_restore_filename  = 'unknown-pt_interior_restore'
   pt_interior_restore_file_fmt  = 'bin'
/

EOF

#--------------------------------------------------------------------------
# forcing_s_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&forcing_s_interior_nml
   s_interior_data_type         = 'none'
   s_interior_data_inc          = 24.
   s_interior_interp_freq       = 'every-timestep'
   s_interior_interp_type       = 'linear'
   s_interior_interp_inc        = 72.
   s_interior_restore_tau       = 365.
   s_interior_filename          = 'unknown-s_interior'
   s_interior_file_fmt          = 'bin'
   s_interior_restore_max_level = 0 
   s_interior_formulation       = 'restoring'
   s_interior_data_renorm(1)    = 1.
   s_interior_variable_restore  = .false.
   s_interior_restore_filename  = 'unknown-s_interior_restore'
   s_interior_restore_file_fmt  = 'bin'
/

EOF

#--------------------------------------------------------------------------
# forcing_ap_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2_NMLFILE << EOF
&forcing_ap_nml
   ap_data_type      = 'none'
   ap_data_inc       = 1.e20
   ap_interp_freq    = 'never'
   ap_interp_type    = 'nearest'
   ap_interp_inc     = 1.e20
   ap_filename       = 'unknown-ap'
   ap_file_fmt       = 'bin'
   ap_data_renorm    = 1.
/

EOF

#--------------------------------------------------------------------------
#  coupled_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then

#... qsw_diurnal_cycle is always false for gx3
 set qsw_diurnal_cycle = .false.

else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
  if ($OCN_COUPLING  =~ *partial*) then
    set qsw_diurnal_cycle = .true.
  else if ($OCN_COUPLING  =~ *full*) then
    set qsw_diurnal_cycle = .true.
  endif
endif 

cat >> $POP2_NMLFILE << EOF
&coupled_nml
   coupled_freq_opt  = 'nday'
   coupled_freq      = 1
   qsw_diurnal_cycle = $qsw_diurnal_cycle
/

EOF

#--------------------------------------------------------------------------
#  sw_absorption_nml
#   OCN_CHL_TYPE is a CCSM environment variable; see $case/env_conf
#--------------------------------------------------------------------------

if ($OCN_CHL_TYPE == prognostic) then
 setenv CHL_OPTION model
else if ($OCN_CHL_TYPE == diagnostic) then
 setenv CHL_OPTION file
endif

cat >> $POP2_NMLFILE << EOF
&sw_absorption_nml
   sw_absorption_type = 'chlorophyll'
   chl_option         = '$CHL_OPTION'
   chl_filename       = '$INPUT/chl_data'
   chl_file_fmt       = 'bin'
   jerlov_water_type  = 3
/

EOF

#--------------------------------------------------------------------------
#  transports_nml
#--------------------------------------------------------------------------

if      ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
 set transport_reg2_names = ("'Atlantic Ocean'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay '")
 set moc = .true.
 set n_heat_trans = .true.
 set n_salt_trans = .true.
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 set transport_reg2_names = ("'Atlantic Ocean'","'Mediterranean Sea'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay'")
 set moc = .true.
 set n_heat_trans = .true.
 set n_salt_trans = .true.
else # turn off transport diags with high-resolution grid
 set moc = .false.
 set n_heat_trans = .false.
 set n_salt_trans = .false.
endif

cat >> $POP2_NMLFILE << EOF
The present code makes assumptions about the region boundaries, so
DO NOT change transport_reg2_names unless you know exactly what you are doing.
&transports_nml
  lat_aux_grid_type    = 'southern'
  lat_aux_begin        = -90.0
  lat_aux_end          =  90.0
  n_lat_aux_grid       = 180 
  moc                  = $moc
  n_heat_trans         = $n_heat_trans
  n_salt_trans         = $n_salt_trans
  transport_reg2_names = $transport_reg2_names
  n_transport_reg      = 2
/

EOF

#--------------------------------------------------------------------------
#  tidal_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx3v6 ) then
 set ltidal_mixing = .true.
else if ( ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
 set ltidal_mixing = .true.
else
 set ltidal_mixing = .false.
endif

cat >> $POP2_NMLFILE << EOF
Only the gx3v6, gx1v4, and gx1v5 versions of tidal_energy exist.  
For all other resolutions, set ltidal_mixing false.

&tidal_nml
  ltidal_mixing          = $ltidal_mixing
  local_mixing_fraction  = 0.33
  mixing_efficiency      = 0.2
  vertical_decay_scale   = 500.0e02
  tidal_mix_max          = 100.0
  tidal_energy_file      = '$INPUT/tidal_energy'
  tidal_energy_file_fmt  = 'bin'
/
EOF

cat >> $POP2_NMLFILE << EOF
&context_nml
   lcoupled          = .true.
   lccsm             = .true.
/

EOF
