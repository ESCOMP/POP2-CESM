#!/bin/csh -f

#===============================================================================
#  ocn.ecosys.setup.csh : perform setup tasks for ecosys module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    set_tavg_nml
#    tavg_contents  tavg_contents_filename my_stream
#    prestage       res_dpt_dir res_indpt_dir
#    document       pop2_document_files
#    ccsm_prestage  ccsm_prestage_file
#
#  CVS:$Id: ocn.ecosys.setup.csh 2448 2006-11-14 21:01:38Z njn01 $ 
#  CVS:$Name$ 
#
#===============================================================================

if ($#argv < 1) then
   echo ocn.ecosys.setup.csh : command argument missing
   exit 1
endif

#===============================================================================
#  set IC and forcing files
#===============================================================================

if (($OCN_TRANSIENT != unset) && ($OCN_TRANSIENT != 1850-2000)) then
   echo OCN_TRANSIENT=$OCN_TRANSIENT not supported by ecosystem module
   exit 2
endif

set relpath = ocn/pop/$OCN_GRID
if ($OCN_GRID == gx3v7) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx3v7_20100514.nc
   set IC_fmt   = nc
   set ALK_scale_factor = 1.025
   set DIC_scale_factor = 1.025
   set O2_scale_factor  = 44.66
   set DST_file   = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx3v7_20100305.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_CCSM4gx3v7_orgc_etopo2v2_20100307.nc
   if ($OCN_TRANSIENT == unset) then
     set ndep_file  = $DIN_LOC_ROOT/$relpath/forcing/ndep_ocn_1850_gx3v7_c100428.nc
   endif
   if ($OCN_TRANSIENT == 1850-2000) then
     set ndep_file = $DIN_LOC_ROOT/$relpath/forcing/ndep_ocn_1850-2005_gx3v7_c100612.nc
   endif

else if ($OCN_GRID == gx1v6) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx1v6_20100514.nc
   set IC_fmt   = nc
   set ALK_scale_factor = 1.025
   set DIC_scale_factor = 1.025
   set O2_scale_factor  = 44.66
   set DST_file   = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx1v6_090416.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_CCSM4gx1v6_orgc_etopo2v2_20090802.nc
   if ($OCN_TRANSIENT == unset) then
     set ndep_file  = $DIN_LOC_ROOT/$relpath/forcing/ndep_ocn_1850_gx1v6_c100428.nc
   endif
   if ($OCN_TRANSIENT == 1850-2000) then
     set ndep_file = $DIN_LOC_ROOT/$relpath/forcing/ndep_ocn_1850-2005_gx1v6_c100612.nc
   endif

else

   echo OCN_GRID=$OCN_GRID not supported by ecosystem module
   exit 2

endif

if ($runtype == startup) then
  if ($OCN_PRESTAGE == TRUE) then
   set IC_file_nml = $INPUT/$IC_file:t
  else
   set IC_file_nml = $IC_file
  endif
else
   set IC_file_nml = same_as_TS
endif

if ($OCN_PRESTAGE == TRUE) then
  set DST_file_nml   = $INPUT/$DST_file:t
  set fesed_file_nml = $INPUT/$fesed_file:t
  set ndep_file_nml  = $INPUT/$ndep_file:t
else
  set DST_file_nml   = $DST_file
  set fesed_file_nml = $fesed_file
  set ndep_file_nml  = $ndep_file
endif

#===============================================================================
#  set module name, which is required for tavg_nml
#  module name must match the name of this setup script
#===============================================================================
set module = ecosys

#===============================================================================

set command = $1

if ($command == set_nt) then

   echo ocn.ecosys.setup.csh : setting nt                                        >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo nt_filename argument missing
      exit 3
   endif

   set nt_filename = $2

   if !(-f $nt_filename) then
      echo nt_filename = $nt_filename does not exist
      exit 3
   endif

   @ nt_in = `cat $nt_filename`
   @ nt_in += 24
   echo $nt_in > $nt_filename

else if ($command == namelist) then

   echo ocn.ecosys.setup.csh : setting namelist options                          >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo pop_in_filename argument missing
      exit 4
   endif

   set pop_in_filename = $2

   if !(-f $pop_in_filename) then
      echo pop_in_filename = $pop_in_filename does not exist
      exit 4
   endif

   if ($OCN_CO2_TYPE == constant) then
      set atm_co2_opt = const
   else if ($OCN_CO2_TYPE == prognostic) then
      set atm_co2_opt = drv_prog
   else if ($OCN_CO2_TYPE == diagnostic) then
      set atm_co2_opt = drv_diag
   else
      echo error specifying atm_co2_opt in ocn.ecosys.setup.csh
      echo unknown OCN_CO2_TYPE:  $OCN_CO2_TYPE
      exit 4
   endif

   if ($OCN_CO2_FLUX_OCMIP_BUG_FIX == TRUE) then
      set locmip_k1_k2_bug_fix = .true.
   else
      set locmip_k1_k2_bug_fix = .false.
   endif

   if ($runtype == startup) then
      set use_nml_surf_vals = .true.
   else
      set use_nml_surf_vals = .false.
   endif

   # note: $runtype is the ccsm script-level variable;
   #       ccsm_$runtype is the pop2 model-level variable
   set init_ecosys_option = ccsm_$runtype

   cat >> $pop_in_filename << EOF

&ecosys_nml
   init_ecosys_option              = '$init_ecosys_option'
   init_ecosys_init_file           = '$IC_file_nml'
   init_ecosys_init_file_fmt       = '$IC_fmt'
   tracer_init_ext(1)%mod_varname = 'ALK'
   tracer_init_ext(1)%scale_factor = $ALK_scale_factor
   tracer_init_ext(2)%mod_varname = 'DIC'
   tracer_init_ext(2)%scale_factor = $DIC_scale_factor
   tracer_init_ext(3)%mod_varname = 'O2'
   tracer_init_ext(3)%scale_factor = $O2_scale_factor
   lflux_gas_o2                    = .true.
   lflux_gas_co2                   = .true.
   locmip_k1_k2_bug_fix            = $locmip_k1_k2_bug_fix
   atm_co2_opt                     = '$atm_co2_opt'
   atm_co2_const                   = $CCSM_CO2_PPMV
   ecosys_tadvect_ctype            = 'base_model'
   gas_flux_forcing_opt            = 'drv'
   lmarginal_seas                  = .true.
   lsource_sink                    = .true.
   comp_surf_avg_freq_opt          = 'never'
   comp_surf_avg_freq              = 1
   use_nml_surf_vals               = $use_nml_surf_vals
   surf_avg_dic_const              = 1944.0
   surf_avg_alk_const              = 2225.0
   ecosys_qsw_distrb_const         = .false.
!  iron_dust_flx_data_type         = 'monthly'
   dust_flux_input%filename        = '$DST_file_nml'
   dust_flux_input%file_fmt        = 'nc'
   dust_flux_input%file_varname    = 'DSTSF'
   dust_flux_input%scale_factor    = 1.0e-1    ! kg/m^2/sec -> g/cm^2/sec
   iron_flux_input%filename        = '$DST_file_nml'
   iron_flux_input%file_fmt        = 'nc'
   iron_flux_input%file_varname    = 'DSTSF'
   iron_flux_input%scale_factor    = 6.2668e4  ! kg/m^2/sec -> nmol/cm^2/sec, 3.5% iron by weight
   fesedflux_input%filename       = '$fesed_file_nml'
   fesedflux_input%file_varname   = 'FESEDFLUXIN'
   fesedflux_input%file_fmt        = 'nc'
   fesedflux_input%scale_factor   = 1.1574e-6 ! umolFe/m2/day -> nmolFe/cm2/s
EOF

if ($OCN_TRANSIENT == unset) then
   cat >> $pop_in_filename << EOF
   ndep_data_type                         = 'monthly-calendar'
   nox_flux_monthly_input%filename        = '$ndep_file_nml'
   nox_flux_monthly_input%file_fmt        = 'nc'
   nox_flux_monthly_input%file_varname    = 'NOy_deposition'
   nox_flux_monthly_input%scale_factor    = 7.1429e+06  ! kgN/m^2/sec -> nmolN/cm^2/sec
   nhy_flux_monthly_input%filename        = '$ndep_file_nml'
   nhy_flux_monthly_input%file_fmt        = 'nc'
   nhy_flux_monthly_input%file_varname    = 'NHx_deposition'
   nhy_flux_monthly_input%scale_factor    = 7.1429e+06  ! kgN/m^2/sec -> nmolN/cm^2/sec
EOF
endif

if ($OCN_TRANSIENT == 1850-2000) then
   cat >> $pop_in_filename << EOF
   ndep_data_type                         = 'shr_stream'
   ndep_shr_stream_year_first             = 1849
   ndep_shr_stream_year_last              = 2006
   ndep_shr_stream_year_align             = 1849
   ndep_shr_stream_file                   = '$ndep_file_nml'
   ndep_shr_stream_scale_factor           = 7.1429e+06  ! kgN/m^2/sec -> nmolN/cm^2/sec
EOF
endif

   cat >> $pop_in_filename << EOF
   lecovars_full_depth_tavg = .false.
/

&ecosys_parms_nml
/
EOF

else if ($command == set_tavg_nml) then

  #-------------------------------------------------------------------------------------
  # if there is no module-related tavg output, set n_tavg_streams_tracer = 0
  #-------------------------------------------------------------------------------------
    set n_tavg_streams_tracer = 2
cat >&! $POP2_DOCDIR/$module.tavg << EOF
n_tavg_streams_tracer =  $n_tavg_streams_tracer
EOF

  #-------------------------------------------------------------------------------------
  # optional: the following definitions are necessary only if n_tavg_streams_tracer > 0 
  #           number of settings must agree with the value of n_tavg_streams_tracer
  #-------------------------------------------------------------------------------------
  if ($n_tavg_streams_tracer > 0) then
    cat >> $POP2_DOCDIR/$module.tavg << EOF
tavg_freq_opt             = 'nday'          'nyear'
tavg_freq                 =  1              1
tavg_stream_filestrings   = 'ecosys.nday1'  'ecosys.nyear1'
tavg_file_freq_opt        = 'nmonth'        'nyear'
tavg_file_freq            =  1              1
tavg_start_opt            = 'nstep'         'nstep'
tavg_start                =  0              0
tavg_fmt_in               = 'nc'            'nc'
tavg_fmt_out              = 'nc'            'nc'
ltavg_has_offset_date     = .false.         .false.
tavg_offset_years         =  1              1
tavg_offset_months        =  1              1
tavg_offset_days          =  2              2
ltavg_one_time_header     = .false.         .false.
EOF
  endif  #n_tavg_streams_tracer

else if ($command == tavg_contents) then

   echo ocn.ecosys.setup.csh : setting tavg_contents variables                   >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 3) then
      echo tavg_contents_filename argument missing
      exit 5
   endif

   set tavg_contents_filename = $2

   if !(-f $tavg_contents_filename) then
      echo tavg_contents_filename = $tavg_contents_filename does not exist
      exit 5
   endif

   @ my_stream = $3
   if ($my_stream < 1) then
      echo invalid my_stream number  ($my_stream)
      exit 5
   endif

   #------------------------------------------------------------------------------------
   # For now, set streams manually. You must only set as many streams as are declared
   #  in the tavg_nml section. For example, if there are three streams:
   #  @ s1 = $my_stream
   #  @ s2 = $s1 + 1
   #  @ s3 = $s2 + 1
   #------------------------------------------------------------------------------------

      @ s1 = 1             # use the base-model stream 1
      @ s2 = $my_stream    # use an ecosystem-defined stream
      @ s3 = $s2 + 1       # use an ecosystem-defined stream

   cat >> $tavg_contents_filename << EOF
$s1  ECOSYS_ATM_PRESS
$s1  ECOSYS_IFRAC
$s1  ECOSYS_XKW
$s1  SCHMIDT_O2
$s1  SCHMIDT_CO2
$s1  IRON_FLUX
$s1  NOx_FLUX
$s1  NHy_FLUX
$s1  PH
$s1  O2SAT
$s1  STF_O2
$s1  CO2STAR
$s1  DCO2STAR
$s1  pCO2SURF
$s1  DpCO2
$s1  FG_CO2
$s1  ATM_CO2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  PO4
$s1  NO3
$s1  SiO3
$s1  NH4
$s1  Fe
$s1  O2
$s1  O2_ZMIN
$s1  O2_ZMIN_DEPTH
$s1  O2_PRODUCTION
$s1  O2_CONSUMPTION
$s1  AOU
$s1  DIC
$s1  J_DIC
$s1  ALK
$s1  H2CO3
$s1  HCO3
$s1  CO3
$s1  pH_3D
$s1  co3_sat_calc
$s1  zsatcalc
$s1  co3_sat_arag
$s1  zsatarag
$s1  DOC
$s1  DOC_prod
$s1  DOC_remin
$s1  spC
$s1  spChl
$s1  spCaCO3
$s1  diatC
$s1  diatChl
$s1  zooC
$s1  spFe
$s1  diatSi
$s1  diatFe
$s1  diazC
$s1  diazChl
$s1  diazFe
$s1  DON
$s1  DOFe
$s1  DOP
$s1  graze_sp
$s1  graze_diat
$s1  graze_diaz
$s1  sp_agg
$s1  diat_agg
$s1  photoC_sp
$s1  CaCO3_form
$s1  photoC_diat
$s1  photoC_diaz
$s1  photoC_NO3_sp
$s1  photoC_NO3_diat
$s1  photoC_NO3_diaz
$s1  Fe_scavenge
$s1  Fe_scavenge_rate
$s1  diaz_Nfix
$s1  bSi_form
$s1  NITRIF
$s1  DENITRIF
$s1  POC_PROD
$s1  CaCO3_PROD
$s1  SiO2_PROD
$s1  P_iron_PROD
$s1  POC_FLUX_IN
$s1  CaCO3_FLUX_IN
$s1  SiO2_FLUX_IN
$s1  P_iron_FLUX_IN
$s1  PAR_avg
$s1  sp_Fe_lim
$s1  diat_Fe_lim
$s1  diaz_Fe_lim
$s1  sp_N_lim
$s1  diat_N_lim
$s1  sp_PO4_lim
$s1  diat_PO4_lim
$s1  diaz_P_lim
$s1  diat_SiO3_lim
$s1  sp_light_lim
$s1  diat_light_lim
$s1  diaz_light_lim
$s1  DON_prod
$s1  DOFe_prod
$s1  DOP_prod
$s1  sp_loss
$s1  diat_loss
$s1  zoo_loss
$s1  diaz_loss
$s1  Jint_100m_DIC
$s1  Jint_100m_NO3
$s1  Jint_100m_NH4
$s1  Jint_100m_PO4
$s1  Jint_100m_Fe
$s1  Jint_100m_SiO3
$s1  Jint_100m_ALK
$s1  Jint_100m_O2
$s1  Jint_100m_DOC
$s1  tend_zint_100m_DIC
$s1  tend_zint_100m_NO3
$s1  tend_zint_100m_NH4
$s1  tend_zint_100m_PO4
$s1  tend_zint_100m_Fe
$s1  tend_zint_100m_SiO3
$s1  tend_zint_100m_ALK
$s1  tend_zint_100m_O2
$s1  tend_zint_100m_DOC
$s2  photoC_sp_zint
$s2  CaCO3_form_zint
$s2  photoC_diaz_zint
$s2  photoC_diat_zint
$s1  photoC_NO3_sp_zint
$s1  photoC_NO3_diat_zint
$s1  photoC_NO3_diaz_zint
$s2  ECOSYS_IFRAC_2
$s2  ECOSYS_XKW_2
$s2  DpCO2_2
$s2  FG_CO2_2
$s2  STF_O2_2
$s2  spC_zint_100m
$s2  spCaCO3_zint_100m
$s2  diazC_zint_100m
$s2  diatC_zint_100m
$s2  zooC_zint_100m
$s2  spChl_SURF
$s2  diazChl_SURF
$s2  diatChl_SURF
$s3  J_NO3
$s3  J_NH4
$s3  J_PO4
$s3  J_Fe
$s3  J_SiO3
$s3  J_ALK
$s3  UE_O2
$s3  VN_O2
$s3  WT_O2
$s3  KPP_SRC_O2
$s3  DIA_IMPVF_O2
$s3  HDIFE_O2
$s3  HDIFN_O2
$s3  HDIFB_O2
$s3  UE_DOC
$s3  VN_DOC
$s3  WT_DOC
$s3  DIA_IMPVF_DOC
$s3  HDIFE_DOC
$s3  HDIFN_DOC
$s3  HDIFB_DOC
$s3  UE_DIC
$s3  VN_DIC
$s3  WT_DIC
$s3  KPP_SRC_DIC
$s3  DIA_IMPVF_DIC
$s3  HDIFE_DIC
$s3  HDIFN_DIC
$s3  HDIFB_DIC
$s3  UE_Fe
$s3  VN_Fe
$s3  WT_Fe
$s3  KPP_SRC_Fe
$s3  DIA_IMPVF_Fe
$s3  HDIFE_Fe
$s3  HDIFN_Fe
$s3  HDIFB_Fe
EOF

#1  dust_FLUX_IN
#1   DON_remin
#1   DOFe_remin
#1   DOP_remin
#1   photoFe_diaz
#1   photoFe_diat
#1   photoFe_sp
#1  Jint_PO4
#1  Jint_NO3
#1  Jint_SiO3
#1  Jint_NH4
#1  Jint_Fe
#1  Jint_O2
#1  Jint_DIC
#1  Jint_ALK
#1  Jint_DOC
#1  Jint_spC
#1  Jint_spChl
#1  Jint_spCaCO3
#1  Jint_diatC
#1  Jint_diatChl
#1  Jint_zooC

else if ($command == prestage) then

 if ($OCN_PRESTAGE == TRUE) then

   echo ocn.ecosys.setup.csh : prestaging data files                             >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 1) then
      echo input data directories missing
      exit 6
   endif

   if ($runtype != continue) then
      \cp -f $IC_file    $IC_file_nml  || exit 6
   endif

   \cp -f $DST_file $DST_file_nml || exit 6
   \cp -f $fesed_file $fesed_file_nml || exit 6
   \cp -f $ndep_file $ndep_file_nml || exit 6

 else

   echo ocn.ecosys.setup.csh : data files are not being prestaged                >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

 endif

else if ($command == document) then

   echo ocn.ecosys.setup.csh : documenting inputdata files                       >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo error documenting input files
      exit 7
   endif

   set pop2_document_files   = $2

  if (-e $IC_file)    \ls -la $IC_file    >> $pop2_document_files
  if (-e $DST_file)   \ls -la $DST_file   >> $pop2_document_files
  if (-e $fesed_file) \ls -la $fesed_file >> $pop2_document_files
  if (-e $ndep_file)  \ls -la $ndep_file  >> $pop2_document_files

else if ($command == ccsm_prestage) then

  #echo ocn.ecosys.setup.csh : writing ccsm prestaging information               >> $POP2_BLDNML
  #echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo error documenting input files
      exit 8
   endif

   set ccsm_prestage_file = $2

   #----------------------------------------------------------------------
   # Note: the following info is used only by the ccsm check_input_data
   # data-checking utility. Everything to the right of the equals sign
   # is ignored by check_input_data.
   #----------------------------------------------------------------------
   if (-f $ccsm_prestage_file) then
     if ($IC_file:h == $DIN_LOC_ROOT/$relpath/ic) then
       echo "    set ecosys_IC_file    = "\$DIN_LOC_ROOT/$relpath/ic/$IC_file:t         >> $ccsm_prestage_file || exit 8
     endif
     if ($DST_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set ecosys_DST_file   = "\$DIN_LOC_ROOT/$relpath/forcing/$DST_file:t   >> $ccsm_prestage_file || exit 8
     endif
     if ($fesed_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set ecosys_fesed_file = "\$DIN_LOC_ROOT/$relpath/forcing/$fesed_file:t >> $ccsm_prestage_file || exit 8
     endif
     if ($ndep_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set ecosys_ndep_file = "\$DIN_LOC_ROOT/$relpath/forcing/$ndep_file:t  >> $ccsm_prestage_file || exit 8
     endif
   endif

else

   echo ocn.ecosys.setup.csh : unrecognized command $command                     >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   exit 9

endif
