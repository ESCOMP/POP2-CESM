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
   echo invalid my_stream number $my_stream
   exit 5
endif

@ s1 = 1   # use base-model stream 1

set irf_tracer_file = $2
set irf_tracer_file_ind_start = $3
set IRF_NT = $4
set IRF_MODE = $5

#
# all modes require VDC_GM
#
echo "1  VDC_GM" >! $CASEROOT/Buildconf/popconf/IRF_tavg_contents

#
# only NK_precond mode requires positive and negative velocity components
#
if ($IRF_MODE == NK_precond) then
   echo "1  UTE_POS" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   echo "1  UTE_NEG" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   echo "1  VTN_POS" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   echo "1  VTN_NEG" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   echo "1  WTK_POS" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   echo "1  WTK_NEG" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
endif

#
# extract impulse variable names from variable var_names in impulse variable file
#
set impulse_dim = `ncdump -h $irf_tracer_file | grep 'impulse_dim =' | head -1 | awk '{print $3}'`
@ impulse_dim_p1 = $impulse_dim + 1
set var_names = ( `ncdump -v var_names $irf_tracer_file | tail -n $impulse_dim_p1 | head -n $impulse_dim | cut -f2 -d\"` )

#
# output lateral mixing tendencies for specified tracers, independent of mode
# output adevective tendencies for offline_transport mode
#
@ count = 1
@ irf_tracer_ind = $irf_tracer_file_ind_start
while ( $count <= $IRF_NT )
   set var_name = $var_names[$irf_tracer_ind]
   echo "1  HDIF_EXPLICIT_3D_$var_name" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   if ($IRF_MODE == offline_transport) then
      echo "1  ADV_3D_$var_name" >> $CASEROOT/Buildconf/popconf/IRF_tavg_contents
   endif
   @ count++
   @ irf_tracer_ind++
end
