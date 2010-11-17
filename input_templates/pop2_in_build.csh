#! /bin/csh -f

#=========================================================================
# create pop2_in by echoing semi-resolved namelists into POP2_IN 
#=========================================================================

#--------------------------------------------------------------------------
#  First, test for supported OCN_GRID resolution  
#--------------------------------------------------------------------------

if      ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6 || ${OCN_GRID} == gx3v7 || ${OCN_GRID} == gx1v5 || ${OCN_GRID} == gx1v5a || ${OCN_GRID} == gx1v5b || ${OCN_GRID} == gx1v6 ) then   
# supported dipole resolutions
else if ( ${OCN_GRID} == tx0.1v2 || ${OCN_GRID} == tx1v1 ) then   
# tripole resolutions
else
   echo " "
   echo "   =============================================================="
   echo "   FATAL ERROR detected building pop2_in:                        "
   echo "     ${OCN_GRID} is not a supported grid.               "
   echo "     Supported grids are: gx3v5, gx3v7, gx1v5, gx1v5a, gx1v5b,   "
   echo "                          gx1v6 and tx0.1v2                      "
  #echo "     Experimental grid: gx3v6                                    "
  #echo "     Testing grids:   tx1v1                                      "
   echo "   =============================================================="
   echo " "
   exit -99
endif



#--------------------------------------------------------------------------
#  domain_nml
#     $NPROCS_CLINIC and $NPROCS_TROPIC are defined in 
#      pop2.buildnml.csh
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx* ) then
  set ns_boundary_type = closed
else if ( ${OCN_GRID} =~ tx* ) then
  set ns_boundary_type = tripole
endif




cat >> $POP2BLDSCRIPT << EOF2

#--------------------------------------------------------------------------
# define variables needed by the pop2_in file
#--------------------------------------------------------------------------
set output_L = \$rundir
set output_d = \$rundir/\$CASE.pop.d

set output_r     = ./\$CASE.pop.r
set output_h     = ./\$CASE.pop.h
set pop2_pointer = ./rpointer.ocn


if ( \$POP_DECOMPTYPE == spacecurve) then
  set clinic_distribution_type = spacecurve
  set tropic_distribution_type = spacecurve
else
  set clinic_distribution_type = cartesian
  set tropic_distribution_type = cartesian
endif

if (\$CALENDAR == GREGORIAN) then
 set allow_leapyear = .true.
else
 set allow_leapyear = .false.
endif

cat >> \$POP2_IN << EOF1

#==========================================================================
#  Begin pop2_in namelist build
#==========================================================================

&domain_nml
  nprocs_clinic            = \$NPROCS_CLINIC
  nprocs_tropic            = \$NPROCS_TROPIC
  clinic_distribution_type = '\$clinic_distribution_type'
  tropic_distribution_type = '\$tropic_distribution_type'
  ew_boundary_type         = 'cyclic'
  ns_boundary_type         = '$ns_boundary_type'
/

EOF2

 
#--------------------------------------------------------------------------
#  io_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&io_nml
  num_iotasks          = 1 
  lredirect_stdout     = .true. 
  log_filename         = '\$output_L/ocn.log.\$LID'
  luse_pointer_files   = .true.
  pointer_filename     = './rpointer.ocn'
  luse_nf_64bit_offset = .true.
/

EOF2

#--------------------------------------------------------------------------
#  time_manager_nml
#    $IYEAR0, $IMONTH0, $IDAY0, and $IHOUR0 are defined in pop2.buildnml_prestage.csh
#--------------------------------------------------------------------------

if      ( ${OCN_GRID} =~ gx3* ) then
  setenv DT_COUNT 12
else if ( ${OCN_GRID} =~ gx1* ) then
  setenv DT_COUNT 23
else if ( ${OCN_GRID} == tx1v1 ) then
  setenv DT_COUNT 23
else if ( ${OCN_GRID} == tx0.1v2 ) then 
  if ($OCN_COUPLING  =~ *partial*) then
    setenv DT_COUNT 500
  else
    setenv DT_COUNT 300
  endif
endif

cat >> $POP2BLDSCRIPT << EOF2

##########################################################
WARNING: DO NOT CHANGE iyear0, imonth0, iday0, ihour0 
##########################################################

&time_manager_nml
  runid             = '\$CASE'
  time_mix_opt      = 'avgfit'
  time_mix_freq     = 17
  dt_option         = 'steps_per_day'
  dt_count          = $DT_COUNT
  impcor            = .true.
  laccel            = .false.
  accel_file        = '\$depth_accel_filename'
  dtuxcel           = 1.0 
  allow_leapyear    = \$allow_leapyear
  iyear0            = $IYEAR0            
  imonth0           = $IMONTH0      
  iday0             = $IDAY0       
  ihour0            = $IHOUR0    
  iminute0          = 0         
  isecond0          = 0        
  date_separator    = '-'
  stop_option       = 'nyear'
  stop_count        =  1000
  fit_freq          = $OCN_NCPL
/

EOF2


#--------------------------------------------------------------------------
#  grid_nml -- see error checking below
#--------------------------------------------------------------------------

set partial_bottom_cells = .false.

if (${OCN_GRID} == tx0.1v2 ) then
  set partial_bottom_cells = .true.
  set bottom_cell_file = ''$INPUT'/dzbc'
endif

set  topography_opt = file

#--------------------------------------------------------------------------
#  some logic to determine proper settings in pop2_in
#--------------------------------------------------------------------------

# Note: topography_opt = bathymetry is a nonstandard option that requires
# the user to provide nonstandard files in the users' <CASE>/SourceMods/src.pop2 directory

set  topography_opt = $topography_opt
if ($topography_opt == 'bathymetry') then
 set lremove_points = .true.
else
 set lremove_points = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&grid_nml
   horiz_grid_opt       = 'file'
   horiz_grid_file      = '\$horiz_grid_filename'
   vert_grid_opt        = 'file'
   vert_grid_file       = '\$vert_grid_filename'
   topography_opt       = '$topography_opt'
   kmt_kmin             = 3
   topography_file      = '\$topography_filename'
   topography_outfile   = '\${output_h}.topography_bathymetry.ieeer8'
   bathymetry_file      = '\$bathymetry_filename'
   partial_bottom_cells =  $partial_bottom_cells
   bottom_cell_file     = '\$bottom_cell_filename'
   n_topo_smooth        = 0
   flat_bottom          = .false.
   lremove_points       =  $lremove_points
   region_mask_file     = '\$regionmask_filename'
   region_info_file     = '\$region_ids_filename'
   sfc_layer_opt        = 'varthick'
/
EOF2


#--------------------------------------------------------------------------
#  init_ts_nml
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
#  some logic to determine proper settings in pop2_in
#--------------------------------------------------------------------------

set init_ts_option = ccsm_\$runtype

if ($runtype == 'startup' && $topography_opt == 'bathymetry') then
   set init_ts_option = PHC
   set init_ts_file = 'ts_PHC2_jan_ic_resindpt'
   set init_ts_file_fmt = nc
else
   set init_ts_file = 'ts'
   set init_ts_file_fmt = \$RESTART_INPUT_TS_FMT
endif

set init_ts_suboption = null
set init_ts_outfile_fmt = nc

if ( ${OCN_GRID} == tx0.1v2) then 
   set init_ts_suboption = spunup
   set init_ts_outfile_fmt = bin
endif


cat >> $POP2BLDSCRIPT << EOF2
&init_ts_nml
   init_ts_option      = '$init_ts_option'
   init_ts_suboption   = '$init_ts_suboption'
EOF2
if ($runtype == 'startup' && $topography_opt == 'bathymetry') then
cat >> $POP2BLDSCRIPT << EOF2
   init_ts_file        = '$init_ts_file'
EOF2
else
cat >> $POP2BLDSCRIPT << EOF2
   init_ts_file        = '\$init_ts_filename'
EOF2
endif
cat >> $POP2BLDSCRIPT << EOF2
   init_ts_file_fmt    = '$init_ts_file_fmt'
   init_ts_outfile     = '\${output_h}.ts_ic'
   init_ts_outfile_fmt = '$init_ts_outfile_fmt'
/

EOF1
EOF2

#--------------------------------------------------------------------------
#  diagnostics_nml
#--------------------------------------------------------------------------

if( ${OCN_GRID} == tx0.1v2) then
   set diag_freq_opt = nday
else
   set diag_freq_opt = nmonth
endif

if( ${OCN_GRID} == tx0.1v2 || ${OCN_GRID} =~ tx1v*) then
   set ldiag_velocity       = .false.
else
   set ldiag_velocity       = .true.
endif

cat >> $POP2BLDSCRIPT << EOF2

if ( \$INFO_DBUG == 2 || \$INFO_DBUG == 3 ) then
  set diag_freq_opt = nstep
else
  set diag_freq_opt = $diag_freq_opt
endif

cat >> \$POP2_IN << EOF1
&diagnostics_nml
   diag_global_freq_opt   = '\$diag_freq_opt'
   diag_global_freq       = 1
   diag_cfl_freq_opt      = '\$diag_freq_opt'
   diag_cfl_freq          = 1
   diag_transp_freq_opt   = '\$diag_freq_opt'
   diag_transp_freq       = 1
   diag_transport_file    = '\$transport_contents_filename'
   diag_outfile           = '\${output_d}d'
   diag_transport_outfile = '\${output_d}t'
   cfl_all_levels         = .false.
   diag_all_levels        = .false.
   diag_velocity_outfile  = '\${output_d}v'
   ldiag_velocity         = $ldiag_velocity
/

EOF1
EOF2


#--------------------------------------------------------------------------
#  budget_diagnostics_nml
#--------------------------------------------------------------------------

  set ldiag_global_tracer_budgets = .true.

if (${OCN_GRID} =~ tx0.1*  || ${OCN_GRID} =~ tx1* ) then
  set ldiag_global_tracer_budgets = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2

if ( \$OCN_TAVG_HIFREQ == FALSE ) then
  set ldiag_global_tracer_budgets = $ldiag_global_tracer_budgets
else
  set ldiag_global_tracer_budgets = .false.
endif

cat >> \$POP2_IN << EOF1
&budget_diagnostics_nml
   ldiag_global_tracer_budgets = \$ldiag_global_tracer_budgets
/

EOF2

#--------------------------------------------------------------------------
#  bsf_diagnostic_nml
#--------------------------------------------------------------------------

  set ldiag_bsf = .true.

if (${OCN_GRID} =~ tx0.1* || ${OCN_GRID} =~ tx1* ) then
  set ldiag_bsf = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&bsf_diagnostic_nml
   ldiag_bsf = $ldiag_bsf
/

EOF2

#--------------------------------------------------------------------------
#  restart_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&restart_nml
   restart_freq_opt    = 'nyear' 
   restart_freq        = 100000
   restart_start_opt   = 'nstep'
   restart_start       =  0
   restart_outfile     = '\$output_r'
   restart_fmt         = 'bin'
   leven_odd_on        = .false. 
   even_odd_freq       = 100000
   pressure_correction = .false.
/

EOF1
EOF2


#--------------------------------------------------------------------------
#  tavg_nml 
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#   first, define the tavg_nml settings for the base model
#--------------------------------------------------------------------------
if ( ${OCN_GRID} =~ gx* ) then
  #------------------- shut off time-invariant stream until vertical grid issues are resolved
  set n_tavg_streams               = 3
  set ltavg_streams_index_present  = .true.
  set tavg_freq_opt_values         = ("'nmonth'" "'nday'"   "'once'" )
  set tavg_freq_values             = (    1         1           1    )
  set tavg_stream_filestrings      = ("'nmonth1'" "'nday1'" "'once'" )
  set tavg_file_freq_opt           = ("'nmonth'" "'nmonth'" "'once'" )
  set tavg_file_freq_values        = (    1         1           1    )
  set tavg_start_opt_values        = ("'nstep'"   "'nstep'" "'nstep'")
  set tavg_start_values            = (    0         0           0    )
  set tavg_fmt_in_values           = ("'nc'"      "'nc'"      "'nc'" )
  set tavg_fmt_out_values          = ("'nc'"      "'nc'"      "'nc'" )
  set ltavg_has_offset_date_values = (.false.    .false.      .false.)
  set tavg_offset_year_values      = (    1         1           1    )
  set tavg_offset_month_values     = (    1         1           1    )
  set tavg_offset_day_values       = (    2         2           2    )
  set ltavg_one_time_header        = (.false.    .false.      .false.)
  set ltavg_nino_diags_requested   = .true. 
else if (${OCN_GRID} =~ tx0.1* || ${OCN_GRID} =~ tx1* ) then
  set n_tavg_streams               = 2
  set ltavg_streams_index_present  = .true.
  set tavg_freq_opt_values         = ("'nmonth'"  "'nday'"  )
  set tavg_freq_values             = (   1           1      )
  set tavg_stream_filestrings      = ("'nmonth1'" "'nday1'" )
  set tavg_file_freq_opt           = ("'nmonth'"  "'nmonth'")
  set tavg_file_freq_values        = (   1           1      )
  set tavg_start_opt_values        = ("'nstep'"   "'nstep'" )
  set tavg_start_values            = (   0           0      )
  set tavg_fmt_in_values           = ( "'nc'"      "'nc'"   )
  set tavg_fmt_out_values          = ( "'nc'"      "'nc'"   )
  set ltavg_has_offset_date_values = (.false.      .false.  )
  set tavg_offset_year_values      = (   1           1      )
  set tavg_offset_month_values     = (   1           1      )
  set tavg_offset_day_values       = (   2           2      )
  set ltavg_one_time_header        = (.false.      .false.  )
  set ltavg_nino_diags_requested   = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2

#############################################################################
# The following setting are for the base model only.
#   Change base-model   tavg_nml setting here.
#   Change tracer-model tavg_nml settings in the ocn.*.setup.csh scripts.
#   The final tavg_nml (base model + extra-tracer models) is constructed below.
#############################################################################

cat >&! \$POP2_TAVG_NML_BASE << EOF1
&tavg_nml
   n_tavg_streams              = $n_tavg_streams
   ltavg_streams_index_present = $ltavg_streams_index_present
   tavg_freq_opt               = $tavg_freq_opt_values
   tavg_freq                   = $tavg_freq_values
   tavg_file_freq_opt          = $tavg_file_freq_opt
   tavg_file_freq              = $tavg_file_freq_values
   tavg_stream_filestrings     = $tavg_stream_filestrings
   tavg_start_opt              = $tavg_start_opt_values
   tavg_start                  = $tavg_start_values
   tavg_fmt_in                 = $tavg_fmt_in_values
   tavg_fmt_out                = $tavg_fmt_out_values
   tavg_contents               ='\$tavg_contents_filename'
   ltavg_nino_diags_requested  = $ltavg_nino_diags_requested
   tavg_infile                 ='\${output_h}restart.end'
   tavg_outfile                ='\$output_h'
   ltavg_has_offset_date       = $ltavg_has_offset_date_values
   tavg_offset_years           = $tavg_offset_year_values
   tavg_offset_months          = $tavg_offset_month_values
   tavg_offset_days            = $tavg_offset_day_values
   ltavg_one_time_header       = $ltavg_one_time_header
/

EOF1
EOF2

#--------------------------------------------------------------------------
#  history_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
cat >> \$POP2_IN << EOF1
&history_nml
   history_freq_opt  = 'never'
   history_freq      = 1
   history_outfile   = '\${output_h}s'
   history_fmt       = 'nc'
   history_contents  = '\$history_contents_filename'
/

EOF2


#--------------------------------------------------------------------------
#  movie_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&movie_nml
   movie_freq_opt = 'never'
   movie_freq     = 1
   movie_outfile  = '\${output_h}m'
   movie_fmt      = 'nc'
   movie_contents = '\$movie_contents_filename'
/

EOF2


#--------------------------------------------------------------------------
#  solver_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&solvers
   solverChoice         = 'ChronGear'
   convergenceCriterion = 1.0e-13 
   maxIterations        = 1000
   convergenceCheckFreq = 10
   preconditionerChoice = 'diagonal'
   preconditionerFile   = 'unknownPrecondFile'
/

EOF2


#--------------------------------------------------------------------------
#  vertical_mix_nml 
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&vertical_mix_nml
   vmix_choice           = 'kpp'
   aidif                 = 1.0
   implicit_vertical_mix = .true.
   convection_type       = 'diffusion'
   nconvad               = 2
   convect_diff          = 10000.0
   convect_visc          = 10000.0
   bottom_drag           = 1.0e-3
   bottom_heat_flx       = 0.0
   bottom_heat_flx_depth = 1000.0e2
/

EOF2


#--------------------------------------------------------------------------
#  vmix_const_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&vmix_const_nml
   const_vvc = 0.25
   const_vdc = 0.25
/

EOF2


#--------------------------------------------------------------------------
#  vmix_rich_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&vmix_rich_nml
   bckgrnd_vvc = 1.0
   bckgrnd_vdc = 0.1
   rich_mix    = 50.0
/

EOF2


#--------------------------------------------------------------------------
#  tidal_nml -- must be defined prior to vmix_kpp_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx3* ) then
 set ltidal_mixing = .true.
else if ( ${OCN_GRID} =~ gx1* ) then
 set ltidal_mixing = .true.
else
 # tidal mixing forcing files exist only for gx1 and gx3 resolutions
 set ltidal_mixing = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&tidal_nml
  ltidal_mixing          = $ltidal_mixing
  local_mixing_fraction  = 0.33
  mixing_efficiency      = 0.2
  vertical_decay_scale   = 500.0e02
  tidal_mix_max          = 100.0
  tidal_energy_file      = '\$tidal_mixing_filename'
  tidal_energy_file_fmt  = 'bin'
/

EOF2


#--------------------------------------------------------------------------
#  vmix_kpp_nml
#--------------------------------------------------------------------------

set llangmuir = .false.
set linertial = .false.
set bckgrnd_vdc_dpth = 1000.0e02

if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6 ) then
 set lhoriz_varying_bckgrnd = .false.
 if ($lhoriz_varying_bckgrnd == .true.) then
   echo "   =============================================================="
   echo "   FATAL ERROR detected building pop2_in:                        "
   echo "     you cannot set lhoriz_varying_bckgrnd = .true. at this     "
   echo "     resolution: ${OCN_GRID}                            "
   echo "   =============================================================="
   echo " "
   exit -99
 else
   if ($ltidal_mixing == .true.) then
     set bckgrnd_vdc1       = 0.2
     set bckgrnd_vdc2       = 0.0
   else
     set bckgrnd_vdc1       = 0.524
     set bckgrnd_vdc2       = 0.313
   endif
 
   set bckgrnd_vdc_eq     = 0.0
   set bckgrnd_vdc_psim   = 0.0
   set bckgrnd_vdc_ban    = 0.0
 endif
else if ( ${OCN_GRID} =~ gx1* || ${OCN_GRID} == gx3v7) then
 set lhoriz_varying_bckgrnd = .true.
 if ($lhoriz_varying_bckgrnd == .true.) then
   if ($ltidal_mixing == .true.) then
     set bckgrnd_vdc1       = 0.16
     set bckgrnd_vdc2       = 0.0
   else
      echo "   =============================================================="
      echo "   FATAL ERROR detected building pop2_in:                        "
      echo "     do not set lhoriz_varying_bckgrnd = .true. and              "
      echo "     ltidal_mixing = .false.  Proper settings for bckgrnd_vdc    "
      echo "     variables are unknown at this time.                         "
      echo "   =============================================================="
      echo " "
      exit -99
   endif
   set bckgrnd_vdc_eq     = 0.01
   set bckgrnd_vdc_psim   = 0.13
   set bckgrnd_vdc_ban    = 1.0
 else
   if ($ltidal_mixing == .true.) then
     set bckgrnd_vdc1       = 0.16
     set bckgrnd_vdc2       = 0.0
   else
     set bckgrnd_vdc1       = 0.524
     set bckgrnd_vdc2       = 0.313
   endif
   set bckgrnd_vdc_eq     = 0.0
   set bckgrnd_vdc_psim   = 0.0
   set bckgrnd_vdc_ban    = 0.0
 endif
else if ( ${OCN_GRID} =~ tx* ) then
 set lhoriz_varying_bckgrnd = .false.
 if ($lhoriz_varying_bckgrnd == .true.) then
   echo "   =============================================================="
   echo "   FATAL ERROR detected building pop2_in:                        "
   echo "     you cannot set lhoriz_varying_bckgrnd = .true. at this     "
   echo "     resolution: ${OCN_GRID}                            "
   echo "   =============================================================="
   echo " "
   exit -99
 endif
 set llangmuir = .false.
 set linertial = .false.
 set bckgrnd_vdc1 = 0.55
 set bckgrnd_vdc2 = 0.303615
 set bckgrnd_vdc_dpth = 2500.0e02
 set bckgrnd_vdc_eq     = 0.0
 set bckgrnd_vdc_psim   = 0.0
 set bckgrnd_vdc_ban    = 0.0
endif

cat >> $POP2BLDSCRIPT << EOF2
&vmix_kpp_nml
   bckgrnd_vdc1           = $bckgrnd_vdc1
   bckgrnd_vdc2           = $bckgrnd_vdc2
   bckgrnd_vdc_eq         = $bckgrnd_vdc_eq
   bckgrnd_vdc_psim       = $bckgrnd_vdc_psim
   bckgrnd_vdc_ban        = $bckgrnd_vdc_ban
   bckgrnd_vdc_dpth       = $bckgrnd_vdc_dpth
   bckgrnd_vdc_linv       = 4.5e-05
   Prandtl                = 10.0
   rich_mix               = 50.0
   lrich                  = .true.
   ldbl_diff              = .true.
   lshort_wave            = .true.
   lcheckekmo             = .false.
   num_v_smooth_Ri        = 1
   lhoriz_varying_bckgrnd = $lhoriz_varying_bckgrnd
   llangmuir              = $llangmuir
   linertial              = $linertial
/

EOF2


#--------------------------------------------------------------------------
#  advect_nml
#--------------------------------------------------------------------------

if( ${OCN_GRID} =~ gx3* || ${OCN_GRID} =~ gx1* ) then 
  set advect_type = 'upwind3'
else if (${OCN_GRID} =~ tx1* ) then
  set advect_type = 'centered'
else if (${OCN_GRID} =~ tx0.1* ) then
  set advect_type = 'centered'
endif

cat >> $POP2BLDSCRIPT << EOF2
&advect_nml
   tadvect_ctype = '$advect_type'
/

EOF2


#--------------------------------------------------------------------------
#  hmix_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx* ) then
  set hmix_momentum_choice = 'anis'
  set hmix_tracer_choice   = 'gent'
  set lsubmesoscale_mixing = .true.
else if (${OCN_GRID} =~ tx1* ) then
  set hmix_momentum_choice = 'del4'
  set hmix_tracer_choice   = 'del4'
  set lsubmesoscale_mixing = .false.
else if (${OCN_GRID} =~ tx0.1* ) then
  set hmix_momentum_choice = 'del4'
  set hmix_tracer_choice   = 'del4'
  set lsubmesoscale_mixing = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&hmix_nml
   hmix_momentum_choice = '$hmix_momentum_choice'
   hmix_tracer_choice   = '$hmix_tracer_choice'
   lsubmesoscale_mixing = $lsubmesoscale_mixing
/

EOF2


#--------------------------------------------------------------------------
#  hmix_del2u_nml 
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx3* ) then
  set lauto_hmix      = .false.
  set lvariable_hmix  = .false.
  set am_del2_value   = 3.0e9
else if ( ${OCN_GRID} =~ gx1* ) then
  set lauto_hmix      = .false.
  set lvariable_hmix  = .false.
  set am_del2_value   = 0.5e8
else if ( ${OCN_GRID} == tx1v1 ) then 
  set lauto_hmix      = .true.
  set lvariable_hmix  = .false.
  set am_del2_value   = 0.5e8
else if ( ${OCN_GRID} == tx0.1v2 ) then 
  set lauto_hmix      = .true.
  set lvariable_hmix  = .false.
  set am_del2_value   = 1.e8
endif

cat >> $POP2BLDSCRIPT << EOF2
&hmix_del2u_nml
   lauto_hmix          = $lauto_hmix 
   lvariable_hmix      = $lvariable_hmix 
   am                  = $am_del2_value
/

EOF2


#--------------------------------------------------------------------------
# hmix_del2t_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx3* ) then
  set lauto_hmix      = .false.
  set lvariable_hmix  = .false.
  set ah_del2_value   = 1.0e7
else if ( ${OCN_GRID} =~ gx1* ) then
  set lauto_hmix      = .false.
  set lvariable_hmix  = .false.
  set ah_del2_value   = 0.6e7
else if ( ${OCN_GRID} == tx1v1 ) then 
  set lauto_hmix      = .false.
  set lvariable_hmix  = .true.
  set ah_del2_value   = 0.6e7
else if ( ${OCN_GRID} == tx0.1v2 ) then 
  set lauto_hmix      = .false.
  set lvariable_hmix  = .true.
  set ah_del2_value   = 1.e8
endif

cat >> $POP2BLDSCRIPT << EOF2
&hmix_del2t_nml
   lauto_hmix          = $lauto_hmix
   lvariable_hmix      = $lvariable_hmix
   ah                  = $ah_del2_value
/

EOF2


#--------------------------------------------------------------------------
# hmix_del4u_nml
#--------------------------------------------------------------------------

  set lauto_hmix       = .false.
  set lvariable_hmix   = .false.
  set am_del4_value    = -0.6e20

if ( ${OCN_GRID} == tx0.1v2 ) then 
  set lauto_hmix       = .false.
  set lvariable_hmix   = .true.
  set am_del4_value    = -27.0e17
endif

cat >> $POP2BLDSCRIPT << EOF2
&hmix_del4u_nml
   lauto_hmix          = $lauto_hmix 
   lvariable_hmix      = $lvariable_hmix
   am                  = $am_del4_value
/

EOF2


#--------------------------------------------------------------------------
# hmix_del4t_nml
#--------------------------------------------------------------------------

  set lauto_hmix       = .false.
  set lvariable_hmix   = .false.
  set ah_del4_value    = -0.2e20

if ( ${OCN_GRID} == tx0.1v2 ) then 
  set lauto_hmix       = .false.
  set lvariable_hmix   = .true.
  set ah_del4_value    = -3.0e17
endif

cat >> $POP2BLDSCRIPT << EOF2
&hmix_del4t_nml
   lauto_hmix          = $lauto_hmix
   lvariable_hmix      = $lvariable_hmix
   ah                  = $ah_del4_value
/

EOF2


#--------------------------------------------------------------------------
#  hmix_gm_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx1* ) then
   set kappa_isop_choice = bfre
   set kappa_thic_choice = bfre
else
   set kappa_isop_choice = bfre
   set kappa_thic_choice = bfre
endif

if ($kappa_isop_choice == 'edgr' && $kappa_thic_choice == 'edgr') then
   # in this instance, the following values are irrelevent
   set ah_gm_value    = 1.0
   set ah_bolus       = 1.0
   set ah_bkg_srfbl   = 1.0
endif

if ($kappa_isop_choice == 'edgr' && $kappa_thic_choice == 'edgr') then
   set use_const_ah_bkg_srfbl = .false.
else
   set use_const_ah_bkg_srfbl = .true.
endif

if ( ${OCN_GRID} =~ gx3* ) then
   set diag_gm_bolus = .true.
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.8e7
   set ah_bolus       = 0.8e7
   set ah_bkg_srfbl   = 0.8e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 4.0e7
   set ah_bolus       = 4.0e7
   set ah_bkg_srfbl   = 4.0e7
 endif
else if ( ${OCN_GRID} =~ gx1* ) then
   set diag_gm_bolus = .true.
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.6e7
   set ah_bolus       = 0.6e7
   set ah_bkg_srfbl   = 0.6e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 3.0e7
   set ah_bolus       = 3.0e7
   set ah_bkg_srfbl   = 3.0e7
 endif
else if ( ${OCN_GRID} == tx1v1 ) then
   set diag_gm_bolus = .false.
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.6e7
   set ah_bolus       = 0.6e7
   set ah_bkg_srfbl   = 0.6e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 3.0e7
   set ah_bolus       = 3.0e7
   set ah_bkg_srfbl   = 3.0e7
 endif
else if ( ${OCN_GRID} == tx0.1v2 ) then
   set diag_gm_bolus = .false.
 if ($kappa_isop_choice == 'constant' && $kappa_thic_choice == 'constant') then 
   set ah_gm_value    = 0.6e7
   set ah_bolus       = 0.6e7
   set ah_bkg_srfbl   = 0.6e7
 else if ($kappa_isop_choice == 'bfre' && $kappa_thic_choice == 'bfre') then 
   set ah_gm_value    = 3.0e7
   set ah_bolus       = 3.0e7
   set ah_bkg_srfbl   = 3.0e7
 endif
endif


cat >> $POP2BLDSCRIPT << EOF2
&hmix_gm_nml
   kappa_isop_choice      = '$kappa_isop_choice'
   kappa_thic_choice      = '$kappa_thic_choice'
   kappa_freq_choice      = 'once_a_day'
   slope_control_choice   = 'notanh'
   kappa_depth_1          = 1.0
   kappa_depth_2          = 0.0
   kappa_depth_scale      = 150000.0
   ah                     = $ah_gm_value
   ah_bolus               = $ah_bolus
   use_const_ah_bkg_srfbl = $use_const_ah_bkg_srfbl
   ah_bkg_srfbl           = $ah_bkg_srfbl
   ah_bkg_bottom          = 0.0
   slm_r                  = 0.3
   slm_b                  = 0.3
   diag_gm_bolus          = $diag_gm_bolus
   transition_layer_on    = .true.
   read_n2_data           = .false.
   buoyancy_freq_filename = '$INPUT/buoyancy_freq'
   buoyancy_freq_fmt      = 'nc'
   const_eg               = 1.2
   gamma_eg               = 500.0
   kappa_min_eg           = 0.35e7
   kappa_max_eg           = 2.0e7
/

EOF2

#--------------------------------------------------------------------------
#  mix_submeso_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&mix_submeso_nml
   efficiency_factor          = 0.07
   time_scale_constant        = 8.64e4
   luse_const_horiz_len_scale = .false.
   hor_length_scale           = 5.0e5
/

EOF2

#--------------------------------------------------------------------------
#  hmix_aniso_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx3* ) then
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
else if ( ${OCN_GRID} =~ gx1* ) then
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
else if ( ${OCN_GRID} == tx1v1 ) then
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
else if ( ${OCN_GRID} == tx0.1v2 ) then
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

cat >> $POP2BLDSCRIPT << EOF2
&hmix_aniso_nml
   hmix_alignment_choice     = '$hmix_alignment_choice'
   lvariable_hmix_aniso      = $lvariable_hmix_aniso
   lsmag_aniso               = $lsmag_aniso
   visc_para                 = $visc_para
   visc_perp                 = $visc_perp
   c_para                    = $c_para
   c_perp                    = $c_perp
   u_para                    = $u_para
   u_perp                    = $u_perp
   vconst_1                  = $vconst_1
   vconst_2                  = $vconst_2
   vconst_3                  = $vconst_3
   vconst_4                  = $vconst_4
   vconst_5                  = $vconst_5
   vconst_6                  = $vconst_6
   vconst_7                  = $vconst_7
   smag_lat                  = 20.0
   smag_lat_fact             = 0.98
   smag_lat_gauss            = 98.0
   var_viscosity_infile      = 'ccsm-internal'
   var_viscosity_infile_fmt  = 'bin'
   var_viscosity_outfile     = '\${output_h}v'
   var_viscosity_outfile_fmt = 'nc'
/

EOF2


#--------------------------------------------------------------------------
# state_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&state_nml
   state_choice     = 'mwjf'
   state_file       = 'internal'
   state_range_opt  = 'enforce'
   state_range_freq = 100000
/

EOF2


#--------------------------------------------------------------------------
# baroclinic_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&baroclinic_nml
   reset_to_freezing = .false.
/

EOF2


#--------------------------------------------------------------------------
#  ice_nml
#   OCN_ICE_FORCING is a CCSM environment variable; see $case/env_conf
#--------------------------------------------------------------------------

if ($OCN_ICE_FORCING =~ *inactive*) then
  set lactive_ice = .false.
else
  set lactive_ice = .true.
endif

cat >> $POP2BLDSCRIPT << EOF2
&ice_nml
   ice_freq_opt     = 'coupled'
   ice_freq         =  100000
   kmxice           = 1
   lactive_ice      = $lactive_ice
/

EOF2


#--------------------------------------------------------------------------
# pressure_grad_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&pressure_grad_nml
   lpressure_avg  = .true.
   lbouss_correct = .false.
/

EOF2


#--------------------------------------------------------------------------
# topostress_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
&topostress_nml
   ltopostress  = .false.
   nsmooth_topo = 0
/

EOF2


#--------------------------------------------------------------------------
# forcing_ws_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
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

EOF2


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

cat >> $POP2BLDSCRIPT << EOF2
&forcing_shf_nml
   shf_formulation      = '$formulation'
   shf_data_type        = '$data_type'
   shf_data_inc         = 24.
   shf_interp_freq      = 'every-timestep'
   shf_interp_type      = 'linear'
   shf_interp_inc       = 72.
   shf_restore_tau      = 30.
   shf_filename         = '\$shf_filename'
   shf_file_fmt         = 'bin'
   shf_data_renorm(3)   = 0.94
   shf_weak_restore     = 0.
   shf_strong_restore   = 0.0
   luse_cpl_ifrac       = $luse_cpl_ifrac
   shf_strong_restore_ms= 92.64
/

EOF2


#--------------------------------------------------------------------------
#  forcing_sfwf_nml
#--------------------------------------------------------------------------

if ( ${OCN_GRID} =~ gx3* ) then
 set sfwf_weak_restore = 0.0115
else if ( ${OCN_GRID} =~ gx1* ) then
 set sfwf_weak_restore = 0.0115
else if ( ${OCN_GRID} == tx1v1 ) then
 set sfwf_weak_restore = 0.0115
else if ( ${OCN_GRID} == tx0.1v2 ) then
 set sfwf_weak_restore = 0.0115
endif

cat >> $POP2BLDSCRIPT << EOF2
&forcing_sfwf_nml
   sfwf_formulation       = '$formulation'
   sfwf_data_type         = '$data_type'
   sfwf_data_inc          = 24.
   sfwf_interp_freq       = 'every-timestep'
   sfwf_interp_type       = 'linear'
   sfwf_interp_inc        = 72.
   sfwf_restore_tau       = 30.
   sfwf_filename          = '\$sfwf_filename'
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

EOF2


#--------------------------------------------------------------------------
# forcing_pt_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
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

EOF2


#--------------------------------------------------------------------------
# forcing_s_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
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

EOF2


#--------------------------------------------------------------------------
# forcing_ap_interior_nml
#--------------------------------------------------------------------------

cat >> $POP2BLDSCRIPT << EOF2
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

EOF2


#--------------------------------------------------------------------------
#  coupled_nml
#--------------------------------------------------------------------------

  set coupled_freq_opt = nhour
  @ coupled_freq = 24 / $OCN_NCPL

if ( ${OCN_GRID} =~ gx3* ) then
 set qsw_distrb_opt = cosz
else if ( ${OCN_GRID} =~ gx1* ) then
  if ($OCN_COUPLING  =~ *partial*) then
    set qsw_distrb_opt = cosz
  else if ($OCN_COUPLING  =~ *full*) then
    set qsw_distrb_opt = cosz
  endif
else if ( ${OCN_GRID} == tx1v1 ) then
    set qsw_distrb_opt = const
else if ( ${OCN_GRID} == tx0.1v2 ) then
    set qsw_distrb_opt = const
endif 

#@ ny = \$ntask / \$NX; setenv NY \$ny

cat >> $POP2BLDSCRIPT << EOF2
&coupled_nml
   coupled_freq_opt  = '$coupled_freq_opt'
   coupled_freq      =  $coupled_freq
   qsw_distrb_opt    = '$qsw_distrb_opt'
/

EOF2


#--------------------------------------------------------------------------
#  sw_absorption_nml
#   OCN_CHL_TYPE is a CCSM environment variable; see $case/env_conf
#--------------------------------------------------------------------------

if ($OCN_CHL_TYPE == prognostic) then
 setenv CHL_OPTION model
else if ($OCN_CHL_TYPE == diagnostic) then
 setenv CHL_OPTION file
endif

cat >> $POP2BLDSCRIPT << EOF2
&sw_absorption_nml
   sw_absorption_type = 'chlorophyll'
   chl_option         = '$CHL_OPTION'
   chl_filename       = '\$chl_filename'
   chl_file_fmt       = 'bin'
   jerlov_water_type  = 3
/

EOF2


#--------------------------------------------------------------------------
#  transports_nml
#--------------------------------------------------------------------------

if      ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
 set transport_reg2_names = ("'Atlantic Ocean'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay '")
 set moc = .true.
 set n_heat_trans = .true.
 set n_salt_trans = .true.
else if ( ${OCN_GRID} =~ gx1* || ${OCN_GRID} == gx3v7 ) then
 set transport_reg2_names = ("'Atlantic Ocean'","'Mediterranean Sea'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay'")
 set moc = .true.
 set n_heat_trans = .true.
 set n_salt_trans = .true.
else # turn off transport diagnostics 
 set transport_reg2_names = ("'Atlantic Ocean'","'Mediterranean Sea'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay'")
 set moc = .false.
 set n_heat_trans = .false.
 set n_salt_trans = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
The present code makes assumptions about the region boundaries, so
DO NOT change transport_reg2_names unless you know exactly what you are doing.
&transports_nml
  lat_aux_grid_type      = 'southern'
  lat_aux_begin          = -90.0
  lat_aux_end            =  90.0
  n_lat_aux_grid         = 180 
  moc_requested          = $moc
  n_heat_trans_requested = $n_heat_trans
  n_salt_trans_requested = $n_salt_trans
  transport_reg2_names   = $transport_reg2_names
  n_transport_reg        = 2
/

EOF2

if ( ${OCN_GRID} == gx1v6 ) then
 set lccsm_control_compatible = .false.
else
 set lccsm_control_compatible = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&context_nml
   lcoupled                 = .true.
   lccsm                    = .true.
   b4b_flag                 = .false.
   lccsm_control_compatible = $lccsm_control_compatible
/

EOF2

#--------------------------------------------------------------------------
#  overflows -- see error checking below
#--------------------------------------------------------------------------

if ( ${OCN_GRID} == gx1v5a || ${OCN_GRID} == gx1v5b  || ${OCN_GRID} =~ gx1v6* ) then
 # ONLY gx1v5a, gx1v5b, or gx1v6 -- not gx1v5
 set overflows_on = .true.
 set overflows_interactive = .true.
else if ( ${OCN_GRID} == gx3v7) then
 set overflows_on = .true.
 set overflows_interactive = .true.
else
 # for all other resolutions
 set overflows_on = .false.
 set overflows_interactive = .false.
endif

cat >> $POP2BLDSCRIPT << EOF2
&overflows_nml
   overflows_on           = $overflows_on
   overflows_interactive  = $overflows_interactive
   overflows_infile       = '\$overflow_filename'
   overflows_diag_outfile = '\${output_d}o'
   overflows_restart_type = 'ccsm_\$runtype'
   overflows_restfile     = '\${output_r}o'
/

EOF1

EOF2

#--------------------------------------------------------------------------
#  Finally, test for conflicting options:
#--------------------------------------------------------------------------

if ($topography_opt == 'bathymetry' && $overflows_on == .true.) then
   echo " "
   echo "   =============================================================="
   echo "   FATAL ERROR detected building pop2_in:                        "
   echo "     You cannot select $topography_opt = 'bathymetry' when       "
   echo "     overflows are active                                        "
   echo "   =============================================================="
   echo " "
   exit -99
endif


if ($ltidal_mixing == .true.) then
  if !($bckgrnd_vdc2 == 0.0) then
   echo "   =============================================================="
   echo "   FATAL ERROR detected building pop2_in:                        "
   echo "     when ltidal_mixing == .true. , you must set bckgrnd_vdc2 = 0.0    "
   echo "     ltidal_mixing = $ltidal_mixing                              "
   echo "     bckgrnd_vdc2  = $bckgrnd_vdc2                               "
   echo "   =============================================================="
   echo " "
   exit -99
  endif
endif

