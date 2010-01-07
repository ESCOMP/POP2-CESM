#!/bin/csh -f

#===============================================================================
#  ocn.cfc.setup.csh : perform setup tasks for cfc module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    set_tavg_nml
#    tavg_contents  tavg_contents_filename my_stream
#    prestage       res_dpt_dir res_indpt_dir
#    document       pop2_document
#    ccsm_prestage  ccsm_prestage_file
#
#  SVN:$Id$
#
#===============================================================================

if ($#argv < 1) then
   echo ocn.cfc.setup.csh : command argument missing
   exit 1
endif

#===============================================================================
#  set forcing files
#===============================================================================

set relpath = ocn/pop/res_indpt
set pcfc_file = $DIN_LOC_ROOT/$relpath/forcing/pcfc1112_atm_20091117.nc

if ($OCN_PRESTAGE == TRUE) then
  set pcfc_file_nml = $INPUT/$pcfc_file:t
else
  set pcfc_file_nml = $pcfc_file
endif

#===============================================================================
#  set module name, which is required for tavg_nml
#  module name must match the name of this setup script
#===============================================================================
set module = cfc

#===============================================================================

set command = $1

if ($command == set_nt) then

   echo ocn.cfc.setup.csh : setting nt                                           >> $POP2_BLDNML
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
   @ nt_in += 2
   echo $nt_in >! $nt_filename

else if ($command == namelist) then

   echo ocn.cfc.setup.csh : setting namelist options                             >> $POP2_BLDNML
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

#===============================================================================
#  some namelist options are different for the
#  Climate of the 20th Century IPCC run
#===============================================================================

   set init_cfc_option = ccsm_$runtype

   if ($OCN_TRANSIENT == 1850-2000) then
      if ($runtype != continue) set init_cfc_option = zero
      set model_year = 1850
      set data_year  = 1850
   else
      set model_year = 1
      set data_year  = 1981
   endif

   cat >> $pop_in_filename << EOF

&cfc_nml
   init_cfc_option    = '$init_cfc_option'
   init_cfc_init_file = 'same_as_TS'
   model_year         = $model_year
   data_year          = $data_year
   pcfc_file          = '$pcfc_file_nml'
   cfc_formulation    = 'model'
/
EOF

else if ($command == set_tavg_nml) then

  #-------------------------------------------------------------------------------------
  # if there is no module-related tavg output, or all tavg output is sent to the base 
  # model output streams, then set n_tavg_streams_tracer = 0
  #-------------------------------------------------------------------------------------
    set n_tavg_streams_tracer = 0
cat >&! $module.tavg << EOF
n_tavg_streams_tracer =  $n_tavg_streams_tracer
EOF

else if ($command == tavg_contents) then

   echo ocn.cfc.setup.csh : setting tavg_contents variables                      >> $POP2_BLDNML
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

   @ s1 = 1   # use base-model stream 1

   cat >> $tavg_contents_filename << EOF
$s1  CFC_IFRAC
$s1  CFC_XKW
$s1  CFC_ATM_PRESS
$s1  STF_CFC11
$s1  STF_CFC12
$s1  CFC11
$s1  CFC12
EOF

#===============================================================================
# The following are fields computed by the CFC modules that are not placed in
# the tavg file by default.
#
#1  pCFC11
#1  pCFC12
#1  CFC11_SCHMIDT
#1  CFC12_SCHMIDT
#1  CFC11_PV
#1  CFC11_surf_sat
#1  CFC12_PV
#1  CFC12_surf_sat
#===============================================================================

else if ($command == prestage) then

   echo ocn.cfc.setup.csh : prestaging data files                                >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 3) then
      echo input data directories missing
      exit 6
   endif

   set res_dpt_dir   = $2
   set res_indpt_dir = $3

if ($OCN_PRESTAGE == TRUE) then
   \cp -f $pcfc_file $pcfc_file_nml || exit 6
endif

else if ($command == document) then

   echo ocn.cfc.setup.csh : documenting inputdata files                          >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   if ($#argv < 2) then
      echo error documenting input files
      exit 7
   endif

   set pop2_document   = $2

  \ls -la $pcfc_file    >> $pop2_document || exit 7


else if ($command == ccsm_prestage) then

  #echo ocn.cfc.setup.csh : writing ccsm prestaging information                  >> $POP2_BLDNML
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
    if ($pcfc_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
      echo "    set cfc_pcfc_file     = "\$DIN_LOC_ROOT/$relpath/forcing/$pcfc_file:t >> $ccsm_prestage_file || exit 8
    endif
   endif


else

   echo ocn.cfc.setup.csh : unrecognized command $command                        >> $POP2_BLDNML
   echo ------------------------------------------------------------------------ >> $POP2_BLDNML

   exit 2

endif
