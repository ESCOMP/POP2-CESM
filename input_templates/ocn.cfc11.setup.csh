#!/bin/csh -fv

#===============================================================================
#  ocn.cfc.setup.csh : perform setup tasks for cfc module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    tavg_contents  tavg_contents_filename
#    prestage       res_dpt_dir res_indpt_dir
#
#  CVS:$Id$
#  CVS:$Name$
#
#===============================================================================

if ($#argv < 1) then
   echo ocn.cfc.setup.csh : command argument missing
   exit 1
endif

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
   @ nt_in += 1
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

   set init_cfc_option = RUNTYPE
   if ($IPCC_MODE == 1870_TO_PRESENT) then
      if ($CONTINUE_RUN == FALSE) set init_cfc_option = zero
      set model_year = 1870
      set data_year  = 1870
   else
      set model_year = 1
      set data_year  = 1981
   endif

   cat >> $pop_in_filename << EOF

&cfc_nml
   init_cfc_option = '$init_cfc_option'
   model_year      = $model_year
   data_year       = $data_year
   pcfc_file       = 'INPUT/pcfc1112_atm.nc'
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

   if ! { grep ATM_PRESS $tavg_contents_filename } then
      echo "2110 'ATM_PRESS' 'Atmospheric Pressure' 'dyne/cm^2'" >> \
         $tavg_contents_filename
   endif

   if ! { grep IFRAC $tavg_contents_filename } then
      echo "2110 'IFRAC' 'Ice Fraction' 'cm^2/cm^2'" >> $tavg_contents_filename
   endif

   if ! { grep U10_SQR $tavg_contents_filename } then
      echo "2110 'U10_SQR' '10m Wind Speed Squared' 'cm^2/s^2'" >> $tavg_contents_filename
   endif

   cat >> $tavg_contents_filename << EOF
2110  'STF_CFC11'  'Surface Flux of CFC-11'  'fmol/cm^2/s'
3111  'CFC11'      'CFC-11'                  'fmol/cm^3'
EOF

#===============================================================================
# The following are fields computed by the CFC modules that by default are not
# placed in the tavg file.
#
# 2110  'pCFC11'         'Atmospheric CFC-11 surface partial pressure'  'pmol/mol'
# 2110  'pCFC12'         'Atmospheric CFC-12 surface partial pressure'  'pmol/mol'
# 2110  'surf_sat_CFC11' 'surface saturation concentration of CFC11'    'fmol/cm^3'
# 2110  'surf_sat_CFC12' 'surface saturation concentration of CFC12'    'fmol/cm^3'
# 2110  'schmidt_CFC11'  'surface Schmidt number for CFC11'             '1'
# 2110  'schmidt_CFC12'  'surface Schmidt number for CFC12'             '1'
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

   \cp -f $res_indpt_dir/forcing/pcfc1112_atm.nc pcfc1112_atm.nc || exit 6

else

   echo -----------------------------------------------------------------
   echo ocn.cfc.setup.csh : unrecognized command $command
   echo -----------------------------------------------------------------

   exit 2

endif
