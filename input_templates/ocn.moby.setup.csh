#!/bin/csh -f
### verbose

#===============================================================================
#  ocn.moby.setup.csh : perform setup tasks for moby module
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
#===============================================================================

if ($#argv < 1) then
   echo ocn.moby.setup.csh : command argument missing
   exit 1
endif

#===============================================================================
#  set IC and forcing files
#===============================================================================

#Note: $DIN_LOC_ROOT references the CESM "inputdata" directory

set relpath = ocn/pop/$OCN_GRID

if ($OCN_GRID == gx3v7) then
   set DST_file   = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx3v7_20100305.nc
   set ALK_scale_factor = 1.025     # to convert ALK to units of mmol C/m^3
   set DIC_scale_factor = 1.025     # to convert DIC to units of mmol C/m^3
else if ($OCN_GRID == gx1v6) then
   set DST_file   = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx1v6_090416.nc
   set ALK_scale_factor = 1.025
   set DIC_scale_factor = 1.025
else
   echo OCN_GRID=$OCN_GRID not supported by moby module
   exit 2
endif

if ($runtype == startup) then
   set IC_file     = 'unknown'
   set IC_file_nml = 'unknown'
else
   set IC_file     = 'unknown'
   set IC_file_nml = same_as_TS
endif
   set IC_fmt = 'nc'

set DST_file_nml   = $DST_file


#===============================================================================
#  set variables whose values depend on the compset vars
#    DEBUG: check that this is the correct CESM flag to key off from
#===============================================================================
if ($CCSM_BGC == CO2C) then
   set lflux_gas_co2 = .true.
else
   set lflux_gas_co2 = .false.
endif

#===============================================================================
#  set module name, which is required for tavg_nml
#  module name must match the name of this setup script
#===============================================================================
set module = moby
set my_model = darwin

#===============================================================================

set command = $1


if ($command == set_nt) then

   echo ocn.moby.setup.csh : setting nt                                          >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo nt_filename argument missing
      exit 3
   endif

   if ($#argv < 3) then
      echo my_path argument missing
      exit 3
   endif

   set nt_filename = $2
   set my_path     = $3

   if !(-f $nt_filename) then
      echo nt_filename = $nt_filename does not exist
      exit 3
   endif

   echo ocn.moby.setup.csh : determine nt from  ${OCN_GRID}_data.ptracers        >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML


   if (-e $my_path/${OCN_GRID}_data.ptracers) then
    set search_dir = $my_path
   else if (-e $OCN_MOBY/$my_model/input/${OCN_GRID}_data.ptracers) then
     set search_dir = $OCN_MOBY/$my_model/input
   else
     exit 31
   endif

   set nt_moby = `grep PTRACERS_numInUse $search_dir/${OCN_GRID}_data.ptracers | cut -f 2 -d = | cut -f 1 -d","`
   if ($status != 0) exit 32


   @ nt_in = `cat $nt_filename`
   @ nt_in += $nt_moby
   echo $nt_in > $nt_filename


else if ($command == namelist) then
#============================================================================
#           The following settings are for the POP ecosystem model           
#           and must be changed to support POP darwin/quota/etc
#     
# Presently, the moby.cpl7.template script performs this funciton; later,
# it should be moved here.
#============================================================================

   echo ocn.moby.setup.csh : setting namelist options                            >> $POP2_BLDNML
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

   if ($runtype == startup) then
      set use_nml_surf_vals = .true.
   else
      set use_nml_surf_vals = .false.
   endif

   #------------------------------------------------------
   # note: $runtype is the ccsm script-level variable;
   #       ccsm_$runtype is the pop2 model-level variable
   #------------------------------------------------------
   set init_moby_option = ccsm_$runtype
   set moby_log_filename = $rundir/moby.log.$LID

   cat >> $pop_in_filename << EOF

&moby_nml
   moby_log_filename             = '$moby_log_filename'
   ldarwin                       = .true.
   init_moby_option              = '$init_moby_option'
   init_moby_init_file           = '$IC_file_nml'
   init_moby_init_file_fmt       = '$IC_fmt'
   lflux_gas_co2                 = $lflux_gas_co2
   moby_tadvect_ctype            = 'base_model'
   lmarginal_seas                = .true.
   comp_surf_avg_freq_opt        = 'never'
   comp_surf_avg_freq            = 1
   moby_qsw_distrb_const         = .false.
   use_nml_surf_vals             = $use_nml_surf_vals
   surf_avg_dic_const            = 1944.0
   surf_avg_alk_const            = 2225.0
   iron_flux_input%filename      = '$DST_file_nml'
   iron_flux_input%file_fmt      = 'nc'
   iron_flux_input%file_varname  = 'DSTSF'
   iron_flux_input%scale_factor  = 0.62668 ! kg/m^2/sec -> mol/m^2/sec, 3.5% iron by weight
   tracer_init_ext(1)%mod_varname = 'ALK'
   tracer_init_ext(1)%scale_factor = $ALK_scale_factor
   tracer_init_ext(2)%mod_varname = 'DIC'
   tracer_init_ext(2)%scale_factor = $DIC_scale_factor
   lecovars_full_depth_tavg      = .false.
/

&moby_parms_nml
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
tavg_stream_filestrings   = '$my_model.nday1'  '$my_model.nyear1'
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

   echo ocn.moby.setup.csh : setting tavg_contents variables                     >> $POP2_BLDNML
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
   # For now, set streams manually. You must only set AS MANY STREAMS AS ARE DEFINED
   #  in the tavg_nml section. For example, if there are three streams, set:
   #  @ s1 = $my_stream
   #  @ s2 = $s1 + 1
   #  @ s3 = $s2 + 1
   #------------------------------------------------------------------------------------

      @ s1 = 1             # use the base-model stream 1
      @ s2 = $my_stream    # use a moby-defined stream
      @ s3 = $s2 + 1       # use a moby-defined stream

   cat >> $tavg_contents_filename << EOF
$s1  PO4
$s1  NO3
$s1  FeT
$s1  SiO2
$s1  DOP
$s1  DON
$s1  DOFe
$s1  ZOO1P
$s1  ZOO1N
$s1  ZOO1Fe
$s1  ZOO1Si
$s1  ZOO2P
$s1  ZOO2N
$s1  ZOO2Fe
$s1  ZOO2Si
$s1  POP
$s1  PON
$s1  POFe
$s1  POSi
$s1  NH4
$s1  NO2
$s1  Phy01
$s1  Phy02
$s1  Phy03
$s1  Phy04
$s1  Phy05
$s1  Phy06
$s1  Phy07
$s1  Phy08
$s1  Phy09
$s1  Chl01
$s1  Chl02
$s1  Chl03
$s1  Chl04
$s1  Chl05
$s1  Chl06
$s1  Chl07
$s1  Chl08
$s1  Chl09
$s1  DIC
$s1  DOC
$s1  POC
$s1  PIC
$s1  ALK
$s1  O2
$s1  ZOOC1
$s1  ZOOC2
#s1  MOBY_CO2_FLUX # same as FG_CO2
$s1  FG_CO2
$s1  STF_O2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  J_DIC
$s1  Jint_100m_PO4
$s1  Jint_100m_NO3
$s1  Jint_100m_FeT
$s1  Jint_100m_SiO2
$s1  Jint_100m_NH4
$s1  Jint_100m_NO2
$s1  Jint_100m_DIC
$s1  Jint_100m_DOC
$s1  Jint_100m_ALK
$s1  Jint_100m_O2
$s1  Jint_100m_ZOOC1
$s1  Jint_100m_ZOOC2
$s1  tend_zint_100m_PO4
$s1  tend_zint_100m_NO3
$s1  tend_zint_100m_FeT
$s1  tend_zint_100m_SiO2
$s1  tend_zint_100m_NH4
$s1  tend_zint_100m_NO2
$s1  tend_zint_100m_DIC
$s1  tend_zint_100m_DOC
$s1  tend_zint_100m_POC
$s1  tend_zint_100m_PIC
$s1  tend_zint_100m_ALK
$s1  tend_zint_100m_O2
$s2  Phy01_zint_100m
$s2  Phy02_zint_100m
$s2  Phy03_zint_100m
$s2  Phy04_zint_100m
$s2  Phy05_zint_100m
$s2  Phy06_zint_100m
$s2  Phy07_zint_100m
$s2  Phy08_zint_100m
$s2  Phy09_zint_100m
$s2  ZOO1P_zint_100m
$s2  ZOO1N_zint_100m
$s2  ZOO1Fe_zint_100m
$s2  ZOO1Si_zint_100m
$s2  ZOO2P_zint_100m
$s2  ZOO2N_zint_100m
$s2  ZOO2Fe_zint_100m
$s2  ZOO2Si_zint_100m
$s2  Chl01_SURF
$s2  Chl02_SURF
$s2  Chl03_SURF
$s2  Chl04_SURF
$s2  Chl05_SURF
$s2  Chl06_SURF
$s2  Chl07_SURF
$s2  Chl08_SURF
$s2  Chl09_SURF
$s3  J_PO4
$s3  J_NO3
$s3  J_FeT
$s3  J_SiO2
$s3  J_NH4
$s3  J_NO2
$s3  J_ALK
$s3  UE_FeT
$s3  VN_FeT
$s3  WT_FeT
$s3  UE_DIC
$s3  VN_DIC
$s3  WT_DIC
$s3  KPP_SRC_DIC
$s3  DIA_IMPVF_DIC
$s3  HDIFE_DIC
$s3  HDIFN_DIC
$s3  HDIFB_DIC
$s3  UE_DOC
$s3  VN_DOC
$s3  WT_DOC
$s3  DIA_IMPVF_DOC
$s3  HDIFE_DOC
$s3  HDIFN_DOC
$s3  HDIFB_DOC
$s3  UE_POC
$s3  VN_POC
$s3  WT_POC
$s3  DIA_IMPVF_POC
$s3  HDIFE_POC
$s3  HDIFN_POC
$s3  HDIFB_POC
$s3  UE_PIC
$s3  VN_PIC
$s3  WT_PIC
$s3  DIA_IMPVF_PIC
$s3  HDIFE_PIC
$s3  HDIFN_PIC
$s3  HDIFB_PIC
$s3  KPP_SRC_FeT
$s3  DIA_IMPVF_FeT
$s3  HDIFE_FeT
$s3  HDIFN_FeT
$s3  HDIFB_FeT
$s3  UE_O2
$s3  VN_O2
$s3  WT_O2
$s3  KPP_SRC_O2
$s3  DIA_IMPVF_O2
$s3  HDIFE_O2
$s3  HDIFN_O2
$s3  HDIFB_O2
# NOTE: the following fields are commented out; to uncomment, delete the leading #
#$s1  Jint_PO4
#$s1  Jint_NO3
#$s1  Jint_Fe
#$s1  Jint_SiO2
#$s1  Jint_ZOO1P
#$s1  Jint_ZOO1N
#$s1  Jint_ZOO1Fe
#$s1  Jint_ZOO1Si
#$s1  Jint_ZOO2P
#$s1  Jint_ZOO2N
#$s1  Jint_ZOO2Fe
#$s1  Jint_ZOO2Si
#$s1  Jint_NH4
#$s1  Jint_NO2
#$s1  Jint_Phy01
#$s1  Jint_Phy02
#$s1  Jint_Phy03
#$s1  Jint_Phy04
#$s1  Jint_Phy05
#$s1  Jint_Phy06
#$s1  Jint_Phy07
#$s1  Jint_Phy08
#$s1  Jint_Phy09
#$s1  Jint_Chl01
#$s1  Jint_Chl02
#$s1  Jint_Chl03
#$s1  Jint_Chl04
#$s1  Jint_Chl05
#$s1  Jint_Chl06
#$s1  Jint_Chl07
#$s1  Jint_Chl08
#$s1  Jint_Chl09
#$s1  Jint_DIC
#$s1  Jint_DOC
#$s1  Jint_POC
#$s1  Jint_PIC
#$s1  Jint_ALK
#$s1  Jint_O2
#$s1  Jint_ZOOC1
#$s1  Jint_ZOOC2
EOF

#1  Jint_zooC

else if ($command == prestage) then

 if ($OCN_PRESTAGE == TRUE) then

   echo ocn.moby.setup.csh : prestaging data files is NOT supported; fatal error >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   exit 6

 else

   echo ocn.moby.setup.csh : data files are not being prestaged                  >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

 endif

else if ($command == document) then

   echo ocn.moby.setup.csh : documenting inputdata files                         >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo error documenting input files
      exit 7
   endif

   set pop2_document_files   = $2

  if (-e $IC_file)    \ls -la $IC_file    >> $pop2_document_files
  if (-e $DST_file)   \ls -la $DST_file   >> $pop2_document_files

else if ($command == ccsm_prestage) then

  #echo ocn.moby.setup.csh : writing ccsm prestaging information                 >> $POP2_BLDNML
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
       echo "    set moby_IC_file    = "\$DIN_LOC_ROOT/$relpath/ic/$IC_file:t         >> $ccsm_prestage_file || exit 8
     endif
     if ($DST_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set moby_DST_file   = "\$DIN_LOC_ROOT/$relpath/forcing/$DST_file:t   >> $ccsm_prestage_file || exit 8
     endif
   endif

else

   echo ocn.moby.setup.csh : unrecognized command $command                       >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   exit 9

endif
