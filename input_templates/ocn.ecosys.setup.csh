#!/bin/csh -f

#===============================================================================
#  ocn.ecosys.setup.csh : perform setup tasks for ecosys module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    tavg_contents  tavg_contents_filename
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

set relpath = ocn/pop/$OCN_GRID_INTERNAL
if ($OCN_GRID_INTERNAL == gx3v5) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx3v5_20051214.nc
   set ALK_scale_factor = 1.0
   set DIC_scale_factor = 1.0
   set O2_scale_factor  = 1.0
   set DST_file = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx3v5_20040426.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_gx3v5_20070521.nc

else if ($OCN_GRID_INTERNAL == gx1v5) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx1v5_20070529.nc
   set ALK_scale_factor = 1.025
   set DIC_scale_factor = 1.025
   set O2_scale_factor  = 44.66
   set DST_file = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx1v5_20070329.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_gx1v5_20070518.nc

else

   echo $OCN_GRID_INTERNAL not supported by ecosystem module
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
else
  set DST_file_nml   = $DST_file
  set fesed_file_nml = $fesed_file
endif

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
      set atm_co2_opt = model
   else
      echo error specifying atm_co2_opt in ocn.ecosys.setup.csh
      echo unknown OCN_CO2_TYPE:  $OCN_CO2_TYPE
      exit 4
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
   init_ecosys_init_file_fmt       = 'nc'
   tracer_init_ext(1)%mod_varname = 'ALK'
   tracer_init_ext(1)%scale_factor = $ALK_scale_factor
   tracer_init_ext(2)%mod_varname = 'DIC'
   tracer_init_ext(2)%scale_factor = $DIC_scale_factor
   tracer_init_ext(3)%mod_varname = 'O2'
   tracer_init_ext(3)%scale_factor = $O2_scale_factor
   lflux_gas_o2                    = .true.
   lflux_gas_co2                   = .true.
   atm_co2_opt                     = '$atm_co2_opt'
   atm_co2_const                   = $CCSM_CO2_PPMV
   ecosys_tadvect_ctype            = 'base_model'
   gas_flux_forcing_opt            = 'model'
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
!   dustflx_daily_input%interp_type    = 'linear'
!   dustflx_daily_input%filename       = '/fiji/home/ivan/data/Mahowald/daily/dstgnx'
!   dustflx_daily_input%data_renorm(1) = 1.e-1    ! kg/m^2/sec -> g/cm^2/sec
!   ironflx_daily_input%interp_type    = 'linear'
!   ironflx_daily_input%filename       = '/fiji/home/ivan/data/Mahowald/daily/dstgnx'
!   ironflx_daily_input%data_renorm(1) = 6.2668e4 ! kg/m^2/sec->nmol/cm^2/sec, 3.5% Fe
/

&ecosys_parms_nml
/
EOF

else if ($command == tavg_contents) then

   echo ocn.ecosys.setup.csh : setting tavg_contents variables                   >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo tavg_contents_filename argument missing
      exit 5
   endif

   set tavg_contents_filename = $2

   if !(-f $tavg_contents_filename) then
      echo tavg_contents_filename = $tavg_contents_filename does not exist
      exit 5
   endif

   cat >> $tavg_contents_filename << EOF
1  ECOSYS_ATM_PRESS
1  ECOSYS_IFRAC
1  ECOSYS_XKW
1  SCHMIDT_O2
1  SCHMIDT_CO2
1  IRON_FLUX
1  PH
1  O2SAT
1  STF_O2
1  CO2STAR
1  DCO2STAR
1  pCO2SURF
1  DpCO2
1  FG_CO2
1  ATM_CO2
1  FvPER_DIC
1  FvICE_DIC
1  FvPER_ALK
1  FvICE_ALK
1  Jint_PO4
1  Jint_NO3
1  Jint_SiO3
1  Jint_NH4
1  Jint_Fe
1  Jint_O2
1  Jint_DIC
1  Jint_ALK
1  Jint_DOC
1  Jint_spC
1  Jint_spChl
1  Jint_spCaCO3
1  Jint_diatC
1  Jint_diatChl
1  Jint_zooC
1  PO4
1  NO3
1  SiO3
1  NH4
1  Fe
1  O2
1  DIC
1  ALK
1  DOC
1  spC
1  spChl
1  spCaCO3
1  diatC
1  diatChl
1  zooC
1  spFe
1  diatSi
1  diatFe
1  diazC
1  diazChl
1  diazFe
1  DON
1  DOFe
1  DOP
1  graze_sp
1  graze_diat
1  graze_diaz
1  sp_agg
1  diat_agg
1  photoC_sp
1  photoC_diat
1  photoC_diaz
1  Fe_scavenge
1  Fe_scavenge_rate
1  CaCO3_form
1  diaz_Nfix
1  bSi_form
1  NITRIF
1  DENITRIF
1  POC_PROD
1  CaCO3_PROD
1  SiO2_PROD
1  P_iron_PROD
1  POC_FLUX_IN
1  CaCO3_FLUX_IN
1  SiO2_FLUX_IN
1  P_iron_FLUX_IN
1  PAR_avg
1  sp_Fe_lim
1  diat_Fe_lim
1  diaz_Fe_lim
1  sp_N_lim
1  diat_N_lim
1  sp_PO4_lim
1  diat_PO4_lim
1  diaz_P_lim
1  diat_SiO3_lim
1  sp_light_lim
1  diat_light_lim
1  diaz_light_lim
1  DOC_prod
1  DON_prod
1  DOFe_prod
1  DOP_prod
1  sp_loss
1  diat_loss
1  zoo_loss
1  diaz_loss
EOF

#1  dust_FLUX_IN
#1   DOC_remin
#1   DON_remin
#1   DOFe_remin
#1   DOP_remin
#1   photoFe_diaz
#1   photoFe_diat
#1   photoFe_sp

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

  \ls -la $IC_file    >> $pop2_document_files || exit 7
  \ls -la $DST_file   >> $pop2_document_files || exit 7
  \ls -la $fesed_file >> $pop2_document_files || exit 7

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
   endif

else

   echo ocn.ecosys.setup.csh : unrecognized command $command                     >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   exit 9

endif
