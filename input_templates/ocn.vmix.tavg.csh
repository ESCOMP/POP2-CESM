#!/bin/csh -f
#set verbose

#------------------------------------------------------------
#  This script adds mixing-related variables to the tavg_contents
#  file, depending on various input namelist settings.
#
#  Variables related to tidal mixing, geothermal heat flux,
#  and near-interial wave mixing are all controlled in this
#  script.
# 
#------------------------------------------------------------

#------------------------------------------------------------
#  Define base-model stream number (s1) and stream number
#    assigned to this module
#------------------------------------------------------------
@ s1          = 1         # base-model stream number 1
@ my_stream   = $argv[1]  # stream number assigned by scripts to this module
@ once_stream = $argv[2]  # stream number of the first "once" stream
if ($once_stream < 1) then
  set sO = '#'
else
  @ sO = $once_stream  # short-hand for once_stream
endif

#------------------------------------------------------------
#  Read script arguments
#------------------------------------------------------------
set ltidal_mixing              = $argv[3]
set ltidal_lunar_cycle         = $argv[4]
set tidal_mixing_method_choice = $argv[5]
set tidal_energy_choice        = $argv[6]
set ltidal_melet_plot          = $argv[7]
set lniw_mixing                = $argv[8]
set geoheatflux_choice         = $argv[9]
set lcvmix                     = $argv[10]

#------------------------------------------------------------
#  Error checking
#------------------------------------------------------------
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 1
endif

#------------------------------------------------------------
#  Clear the temporary work file, vmix_tavg_contents
#------------------------------------------------------------
rm -f $CASEROOT/Buildconf/popconf/vmix_tavg_contents
touch $CASEROOT/Buildconf/popconf/vmix_tavg_contents

#------------------------------------------------------------
#  Append tidal_mixing variables to vmix_tavg_contents file
#------------------------------------------------------------
if ($ltidal_mixing =~ ".true.") then
  #--------------------------------
  #  General tidal-mixing variables
  #--------------------------------
  echo "#..Tidal Mixing"     >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  TIDAL_DIFF"     >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$sO  VERTICAL_FUNC"  >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$sO  REGION_BOX3D"  >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents

  #... these variables are also supported, but probably not of general interest:
  # echo "$s1  TIDAL_N2"       >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  # echo "$s1  TIDAL_N2_eps"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  # echo "$s1  TIDAL_GAMMA_EPS">> $CASEROOT/Buildconf/popconf/vmix_tavg_contents

  #-------------------------
  #  Polzin method variables
  #-------------------------
  if ($tidal_mixing_method_choice == "polzin") then
    echo "#..Tidal Mixing -- Polzin" \
                               >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    echo "$sO  H2_P"           >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    echo "$sO  U_P"            >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    if ($lcvmix =~ ".false.") then 
    # can only accumulate in 2D mode, not column mode
    echo "$sO  POLZIN_EQ2"     >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    endif
  endif

  #------------------------------------------------------
  #  Variables used to recreate plots in Melet et al 2013
  #------------------------------------------------------
  if ($ltidal_melet_plot =~ ".true." && $OCN_GRID == "gx1v6") then
      echo "#..Tidal Mixing -- Melet plot variables" \
                              >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP1_2km"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP1_1p5km" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP1_3km"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP2_3km"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP2_4km"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP3p5_6km" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    # echo "$s1  TEMP4_6km"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  endif

  #---------------------------------------------
  #  Time-varying tidal coefficients 2D methods
  #---------------------------------------------
  if ($tidal_mixing_method_choice == "polzin" || $tidal_mixing_method_choice == "jayne") then
    if ($ltidal_lunar_cycle =~ ".true.") then
      echo "$s1  TIDAL_QE_2D"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    else
      echo "$sO  TIDAL_QE_2D"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    endif
  endif

#------------- for plotting
      echo "$sO  TIDAL_QE_2D"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
#------------- for plotting

  #------------------------------
  #  3D-method tidal coefficients
  #------------------------------
  if ($tidal_mixing_method_choice == "Schmittner" ) then
    if ($ltidal_lunar_cycle =~ ".true.") then
      echo "$s1  TIDAL_COEF_3D" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    else
      echo "$sO  TIDAL_COEF_3D" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    endif
  endif

  #------------------------
  #  3D tidal-energy fields
  #------------------------
  if ($tidal_energy_choice == "ER03" || $tidal_energy_choice == "GN13"     || \
      $tidal_energy_choice == "LGM0" || $tidal_energy_choice == "LGMi5g21" || \
      $tidal_energy_choice == "LGMi6g21") then
    echo "#..3D Tidal Energy"  >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    if ($ltidal_lunar_cycle =~ ".true.") then
      echo "$s1  TIDAL_QE_3D"  >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    else
      echo "$sO  TIDAL_QE_3D"  >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    endif
    echo "$sO  TCM2"           >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    echo "$sO  TCS2"           >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    echo "$sO  TCK1"           >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
    echo "$sO  TCO1"           >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  endif

endif # ltidal_mixing

#------------------------------------------------------------
#  Append geoheatflux variables to vmix_tavg_contents file
#------------------------------------------------------------
if ($geoheatflux_choice == "spatial") then
  echo "#..Spatially varying geothermal heat flux" \
                       >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$sO  GEOHFLUX" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
endif

#------------------------------------------------------------
#  Append niw_mixing variables to vmix_tavg_contents file
#------------------------------------------------------------
if ($lniw_mixing == ".true.") then
  echo "#..Near-inertial wave mixing" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  KVNIW"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  N2"      >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  BFNIW"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  KVNIW_M" >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  KE_BL"   >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
  echo "$s1  En"      >> $CASEROOT/Buildconf/popconf/vmix_tavg_contents
endif

exit
