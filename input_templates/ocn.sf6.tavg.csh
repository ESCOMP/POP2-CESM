#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1   # use base-model stream 1

cat >! $CASEROOT/Buildconf/popconf/sf6_tavg_contents << EOF
$s1  SF6_IFRAC
$s1  SF6_XKW
$s1  SF6_ATM_PRESS
$s1  STF_SF6
$s1  SF6
EOF

if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $CASEROOT/Buildconf/popconf/sf6_tavg_contents << EOF
$s1  KPP_SRC_SF6
$s1  DIA_IMPVF_SF6
$s1  HDIFE_SF6
$s1  HDIFN_SF6
$s1  HDIFB_SF6
$s1  UE_SF6
$s1  VN_SF6
$s1  WT_SF6
EOF
endif

#===============================================================================
# The following are fields computed by the SF6 modules that are not placed in
# the tavg file by default.
#
#1  pSF6
#1  SF6_SCHMIDT
#1  SF6_PV
#1  SF6_surf_sat
#===============================================================================
