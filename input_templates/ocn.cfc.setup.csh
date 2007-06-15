#!/bin/csh -f

#===============================================================================
#  ocn.cfc.setup.csh : perform setup tasks for cfc module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    tavg_contents  tavg_contents_filename
#    prestage       res_dpt_dir res_indpt_dir
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

set pcfc_file = pcfc1112_atm_20070125.nc

#===============================================================================

set command = $1

if ($command == set_nt) then

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : setting nt
   echo -----------------------------------------------------------------

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

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : setting namelist options
   echo -----------------------------------------------------------------

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

   set init_cfc_option = $RUN_TYPE
   if ($CONTINUE_RUN == TRUE) set init_cfc_option = continue

   if ($CCSM_IPCC == 1870_to_2000_control) then
      if ($CONTINUE_RUN == FALSE) set init_cfc_option = zero
      set model_year = 1870
      set data_year  = 1870
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
   pcfc_file          = '$INPUT/$pcfc_file'
   cfc_formulation    = 'model'
/
EOF

else if ($command == tavg_contents) then

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : setting tavg_contents variables
   echo -----------------------------------------------------------------

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
CFC_IFRAC
CFC_XKW
CFC_ATM_PRESS
STF_CFC11
STF_CFC12
CFC11
CFC12
EOF

#===============================================================================
# The following are fields computed by the CFC modules that are not placed in
# the tavg file by default.
#
# pCFC11
# pCFC12
# CFC11_SCHMIDT
# CFC12_SCHMIDT
# CFC11_PV
# CFC11_surf_sat
# CFC12_PV
# CFC12_surf_sat
#===============================================================================

else if ($command == prestage) then

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : prestaging data files
   echo -----------------------------------------------------------------

   if ($#argv < 3) then
      echo input data directories missing
      exit 6
   endif

   set res_dpt_dir   = $2
   set res_indpt_dir = $3

   \cp -f $res_indpt_dir/forcing/$pcfc_file . || exit 6

else

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : unrecognized command $command
   echo -----------------------------------------------------------------

   exit 2

endif
