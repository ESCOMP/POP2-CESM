#! /bin/csh -f


#-----------------------------------------------------------------------------
#  namelist:  the following default settings are selected based upon 
#             coupling to active/inactive atm/ice models
#             *AND* ocean-model resolution
#-----------------------------------------------------------------------------


#======================================================
if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then

#             qsw_diurnal_cycle is always false
#======================================================

cat >> $POP2BLDSCRIPT << EOF3

#..... coupled_nml
 set qsw_diurnal_cycle = .false.

EOF3

#============================================================
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4) then

#             qsw_diurnal_cycle is false unless fully coupled
#============================================================

  #------------------------------------
  if ($OCN_COUPLING  =~ *partial*) then
  #------------------------------------

cat >> $POP2BLDSCRIPT << EOF3

#..... coupled_nml
 set qsw_diurnal_cycle = .true.

EOF3

  #--------------------------------------
  else if ($OCN_COUPLING  =~ *full*) then
  #--------------------------------------

cat >> $POP2BLDSCRIPT << EOF3

#..... coupled_nml
 set qsw_diurnal_cycle = .true.

EOF3

  #----
  endif
  #----

endif 
#====
