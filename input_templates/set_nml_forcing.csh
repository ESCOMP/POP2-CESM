#! /bin/csh -f

cat >> $POP2BLDSCRIPT << EOF3

#-----------------------------------------------------------------------------
#  namelist:  the following settings are selected based upon 
#             coupling to active/inactive atm/ice models
#-----------------------------------------------------------------------------

EOF3
 

if ($OCN_COUPLING  =~ *partial*) then
 cat >> $POP2BLDSCRIPT << EOF3
#... forcing_shf_nml and forcing_sfwf_nml
  set formulation       = partially-coupled
  set data_type         = monthly

#... forcing_sfwf_nml
  set ladjust_precip    = .true.
  set lms_balance       = .false.
  set lsend_precip_fact = .true.

EOF3
else
 cat >> $POP2BLDSCRIPT << EOF3
#... forcing_shf_nml and forcing_sfwf_nml
  set formulation       = restoring
  set data_type         = none

#... forcing_sfwf_nml
  set ladjust_precip    = .false.
  set lms_balance       = .true.
  set lsend_precip_fact = .false.
EOF3
endif

if ($OCN_ICE_FORCING =~ *inactive*) then
 cat >> $POP2BLDSCRIPT << EOF3
#... ice_nml
  set lactive_ice = .false.
EOF3
else
 cat >> $POP2BLDSCRIPT << EOF3
#... ice_nml
  set lactive_ice = .true.
EOF3
endif

