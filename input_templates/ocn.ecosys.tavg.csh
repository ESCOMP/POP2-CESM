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

@ s1 = 1             # use the base-model stream 1
@ s2 = $my_stream    # use an ecosystem-defined stream
@ s3 = $s2 + 1       # use an ecosystem-defined stream

set lecosys_debug = $2

if ($lecosys_debug == ".false.") then

  cat >! $CASEBUILD/popconf/ecosys.tavg.nml << EOF
tavg_freq_opt             = 'nday'   'nyear'
tavg_freq                 =  1       1
tavg_file_freq_opt        = 'nmonth' 'nyear'
tavg_file_freq            =  1       1
tavg_start_opt            = 'nstep'  'nstep'
tavg_start                =  0       0
tavg_fmt_in               = 'nc'     'nc'
tavg_fmt_out              = 'nc'     'nc'
ltavg_has_offset_date     = .false.  .false.
tavg_offset_years         =  1       1
tavg_offset_months        =  1       1
tavg_offset_days          =  2       2
ltavg_one_time_header     = .false.  .false.
tavg_stream_filestrings   = 'ecosys.nday1' 'ecosys.nyear1'
EOF

else
  rm -f $CASEBUILD/popconf/ecosys.tavg.nml
  touch $CASEBUILD/popconf/ecosys.tavg.nml
endif

if ($lecosys_debug == ".true.") then

  cat >! $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
#  GENERAL INTERIOR DIAGNOSTICS
1  zsatcalc
1  zsatarag
1  O2_ZMIN
1  O2_ZMIN_DEPTH
1  photoC_TOT_zint
1  photoC_NO3_TOT_zint
1  Jint_Ctot
1  Jint_100m_Ctot
1  Jint_Ntot
1  Jint_100m_Ntot
1  Jint_Ptot
1  Jint_100m_Ptot
1  Jint_Sitot
1  Jint_100m_Sitot
1  Jint_Fetot
1  Jint_100m_Fetot
1  CO3
1  HCO3
1  H2CO3
1  pH_3D
1  CO3_ALT_CO2
1  HCO3_ALT_CO2
1  H2CO3_ALT_CO2
1  pH_3D_ALT_CO2
1  co3_sat_calc
1  co3_sat_arag
1  NITRIF
1  DENITRIF
1  O2_PRODUCTION
1  O2_CONSUMPTION
1  AOU
1  PAR_avg
1  graze_auto_TOT
1  photoC_TOT
1  photoC_NO3_TOT
1  DOC_prod
1  DOC_remin
1  DOCr_remin
1  DON_prod
1  DON_remin
1  DONr_remin
1  DOP_prod
1  DOP_remin
1  DOPr_remin
1  Fe_scavenge
1  Fe_scavenge_rate
#  PARTICULATE
1  calcToSed
1  pocToSed
1  ponToSed
1  SedDenitrif
1  OtherRemin
1  popToSed
1  bsiToSed
1  dustToSed
1  pfeToSed
1  POC_FLUX_IN
1  POC_PROD
1  POC_REMIN
1  POC_REMIN_DIC
1  PON_REMIN_NH4
1  POP_REMIN_PO4
1  CaCO3_FLUX_IN
1  CaCO3_PROD
1  CaCO3_REMIN
1  SiO2_FLUX_IN
1  SiO2_PROD
1  SiO2_REMIN
1  dust_FLUX_IN
1  dust_REMIN
1  P_iron_FLUX_IN
1  P_iron_PROD
1  P_iron_REMIN
#  FORCING FIELDS
1  ECOSYS_IFRAC
1  ECOSYS_XKW
1  ECOSYS_ATM_PRESS
1  PV_O2
1  SCHMIDT_O2
1  O2SAT
1  CO2STAR
1  DCO2STAR
1  pCO2SURF
1  DpCO2
1  PV_CO2
1  SCHMIDT_CO2
1  FG_CO2
1  PH
1  ATM_CO2
1  CO2STAR_ALT_CO2
1  DCO2STAR_ALT_CO2
1  pCO2SURF_ALT_CO2
1  DpCO2_ALT_CO2
1  FG_ALT_CO2
1  PH_ALT_CO2
1  ATM_ALT_CO2
1  IRON_FLUX
1  DUST_FLUX
1  DUST_FLUX_CPL
1  BLACK_CARBON_FLUX_CPL
1  NOx_FLUX
1  NHy_FLUX
1  DIN_RIV_FLUX
1  DIP_RIV_FLUX
1  DON_RIV_FLUX
1  DONr_RIV_FLUX
1  DOP_RIV_FLUX
1  DOPr_RIV_FLUX
1  DSI_RIV_FLUX
1  DFE_RIV_FLUX
1  DIC_RIV_FLUX
1  ALK_RIV_FLUX
1  DOC_RIV_FLUX
1  DOCr_RIV_FLUX
#  DUPLICATE TAVG VARS
1  ECOSYS_IFRAC_2
1  ECOSYS_XKW_2
1  DpCO2_2
1  FG_CO2_2
1  STF_O2_2
EOF

  # interior autotroph fields
  foreach autotroph ( sp diat diaz )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
#  AUTOTROPHS (${autotroph})
1  photoC_${autotroph}_zint
1  photoC_NO3_${autotroph}_zint
1  ${autotroph}_N_lim
1  ${autotroph}_P_lim
1  ${autotroph}_Fe_lim
1  ${autotroph}_light_lim
1  photoC_${autotroph}
1  photoC_NO3_${autotroph}
1  photoFe_${autotroph}
1  photoNO3_${autotroph}
1  photoNH4_${autotroph}
1  DOP_${autotroph}_uptake
1  PO4_${autotroph}_uptake
1  graze_${autotroph}
1  graze_${autotroph}_poc
1  graze_${autotroph}_doc
1  graze_${autotroph}_zoo
1  ${autotroph}_loss
1  ${autotroph}_loss_poc
1  ${autotroph}_loss_doc
1  ${autotroph}_agg
EOF
  end
  cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
1  sp_CaCO3_form_zint
1  sp_CaCO3_form
1  diat_SiO3_lim
1  diatbSi_form
1  diaz_Nfix
1  CaCO3_form_zint
1  bSi_form
1  CaCO3_form
1  Nfix
EOF

  foreach zooplankton ( zoo )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
#  ZOOPLANKTON (${zooplankton})
1  ${zooplankton}_loss
1  ${zooplankton}_loss_poc
1  ${zooplankton}_loss_doc
1  graze_${zooplankton}
1  graze_${zooplankton}_poc
1  graze_${zooplankton}_doc
1  graze_${zooplankton}_zoo
1  x_graze_${zooplankton}
EOF
  end

  foreach tracer ( PO4 NO3 SiO3 NH4 Fe O2 DIC DIC_ALT_CO2 ALK DOC DON DOCr    \
                   DOP DOPr DONr zooC spChl spC spFe spCaCO3 diatChl diatC    \
                   diatFe diatSi diazChl diazC diazFe )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
#  TRACERS (${tracer})
1  STF_${tracer}
1  J_${tracer}
1  ${tracer}_RESTORE
EOF
  end
  cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
1  FvPER_DIC
1  FvPER_ALK
1  FvICE_DIC
1  FvICE_ALK
EOF

else # ecosys not in debug mode
  cat >! $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ECOSYS_ATM_PRESS
$s1  ECOSYS_IFRAC
$s1  ECOSYS_XKW
$s1  SCHMIDT_O2
$s1  SCHMIDT_CO2
$s1  IRON_FLUX
$s1  NOx_FLUX
$s1  NHy_FLUX
$s1  STF_ALK
$s1  PH
$s1  O2SAT
$s1  STF_O2
$s1  CO2STAR
$s1  DCO2STAR
$s1  pCO2SURF
$s1  DpCO2
$s1  FG_CO2
$s1  ATM_CO2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  PO4
$s1  NO3
$s1  SiO3
$s1  NH4
$s1  Fe
$s1  O2
$s1  O2_ZMIN
$s1  O2_ZMIN_DEPTH
$s1  O2_PRODUCTION
$s1  O2_CONSUMPTION
$s1  AOU
$s1  DIC
$s1  J_DIC
$s1  ALK
$s1  H2CO3
$s1  HCO3
$s1  CO3
$s1  pH_3D
$s1  co3_sat_calc
$s1  zsatcalc
$s1  co3_sat_arag
$s1  zsatarag
$s1  DOC
$s1  DOC_prod
$s1  DOC_remin
$s1  DOCr_remin
$s1  zooC
$s1  DON
$s1  DON_remin
$s1  DONr_remin
$s1  DOP
$s1  DOP_remin
$s1  DOPr_remin
$s1  DONr
$s1  DOPr
$s1  DOCr
$s1  DIN_RIV_FLUX
$s1  DIP_RIV_FLUX
$s1  DON_RIV_FLUX
$s1  DONr_RIV_FLUX
$s1  DOP_RIV_FLUX
$s1  DOPr_RIV_FLUX
$s1  DOC_RIV_FLUX
$s1  DOCr_RIV_FLUX
$s1  DSI_RIV_FLUX
$s1  DFE_RIV_FLUX
$s1  DIC_RIV_FLUX
$s1  ALK_RIV_FLUX
$s1  calcToSed
$s1  pocToSed
$s1  ponToSed
$s1  popToSed
$s1  pfeToSed
$s1  dustToSed
$s1  SedDenitrif
$s1  bsiToSed
$s1  CaCO3_form
$s1  Fe_scavenge
$s1  Fe_scavenge_rate
$s1  bSi_form
$s1  NITRIF
$s1  DENITRIF
$s1  POC_PROD
$s1  POC_REMIN
$s1  POC_REMIN_DIC
$s1  PON_REMIN_NH4
$s1  POP_REMIN_PO4
$s1  CaCO3_PROD
$s1  SiO2_PROD
$s1  P_iron_PROD
$s1  POC_FLUX_IN
$s1  CaCO3_FLUX_IN
$s1  SiO2_FLUX_IN
$s1  P_iron_FLUX_IN
$s1  dust_FLUX_IN
$s1  PAR_avg
$s1  DON_prod
$s1  DOP_prod
$s1  zoo_loss
$s1  Jint_100m_DIC
$s1  Jint_100m_NO3
$s1  Jint_100m_NH4
$s1  Jint_100m_PO4
$s1  Jint_100m_Fe
$s1  Jint_100m_SiO3
$s1  Jint_100m_ALK
$s1  Jint_100m_O2
$s1  Jint_100m_DOC
$s1  tend_zint_100m_DIC
$s1  tend_zint_100m_NO3
$s1  tend_zint_100m_NH4
$s1  tend_zint_100m_PO4
$s1  tend_zint_100m_Fe
$s1  tend_zint_100m_SiO3
$s1  tend_zint_100m_ALK
$s1  tend_zint_100m_O2
$s1  tend_zint_100m_DOC
$s2  CaCO3_form_zint
$s2  ECOSYS_IFRAC_2
$s2  ECOSYS_XKW_2
$s2  DpCO2_2
$s2  FG_CO2_2
$s2  STF_O2_2
$s2  zooC_zint_100m
$s3  J_NO3
$s3  J_NH4
$s3  J_PO4
$s3  J_Fe
$s3  J_SiO3
$s3  J_ALK
$s3  UE_O2
$s3  VN_O2
$s3  WT_O2
$s3  KPP_SRC_O2
$s3  DIA_IMPVF_O2
$s3  HDIFE_O2
$s3  HDIFN_O2
$s3  HDIFB_O2
$s3  UE_DOC
$s3  VN_DOC
$s3  WT_DOC
$s3  DIA_IMPVF_DOC
$s3  HDIFE_DOC
$s3  HDIFN_DOC
$s3  HDIFB_DOC
$s3  UE_DIC
$s3  VN_DIC
$s3  WT_DIC
$s3  KPP_SRC_DIC
$s3  DIA_IMPVF_DIC
$s3  HDIFE_DIC
$s3  HDIFN_DIC
$s3  HDIFB_DIC
$s3  UE_Fe
$s3  VN_Fe
$s3  WT_Fe
$s3  KPP_SRC_Fe
$s3  DIA_IMPVF_Fe
$s3  HDIFE_Fe
$s3  HDIFN_Fe
$s3  HDIFB_Fe
EOF

  # generic autotroph fields
  # skip N_lim for diaz
  foreach autotroph ( sp diat diaz )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}Chl
$s1  ${autotroph}C
$s1  ${autotroph}Fe
$s1  graze_${autotroph}
$s1  ${autotroph}_agg
$s1  photoC_${autotroph}
$s1  photoC_NO3_${autotroph}
$s1  photoNO3_${autotroph}
$s1  photoNH4_${autotroph}
$s1  photoFe_${autotroph}
$s1  DOP_${autotroph}_uptake
$s1  PO4_${autotroph}_uptake
$s1  ${autotroph}_Fe_lim
$s1  ${autotroph}_P_lim
$s1  ${autotroph}_light_lim
$s1  ${autotroph}_loss
$s2  photoC_${autotroph}_zint
$s1  photoC_NO3_${autotroph}_zint
$s2  ${autotroph}C_zint_100m
$s2  ${autotroph}Chl_SURF
EOF
    if !($autotroph == diaz) then
      cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}_N_lim
EOF
    endif
  end

  # Nfix terms from N fixers 
  foreach autotroph ( diaz )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}_Nfix
EOF
  end

  # CaCO3 terms from calcifiers 
  foreach autotroph ( sp )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}CaCO3
$s2  ${autotroph}CaCO3_zint_100m
EOF
  end

  # Si terms from silicifiers
  foreach autotroph ( diat )
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}Si
$s1  ${autotroph}_SiO3_lim
EOF
  end

  setenv OCN_TAVG_DIC_ALT_CO2 FALSE

  if ($OCN_TAVG_DIC_ALT_CO2 == TRUE) then
cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  PH_ALT_CO2
$s1  DCO2STAR_ALT_CO2
$s1  DpCO2_ALT_CO2
$s1  FG_ALT_CO2
$s1  ATM_ALT_CO2
$s1  DIC_ALT_CO2
!  CO3_ALT_CO2
!  pH_3D_ALT_CO2
$s1  tend_zint_100m_DIC_ALT_CO2
$s3  UE_DIC_ALT_CO2
$s3  VN_DIC_ALT_CO2
$s3  WT_DIC_ALT_CO2
$s3  KPP_SRC_DIC_ALT_CO2
$s3  DIA_IMPVF_DIC_ALT_CO2
$s3  HDIFE_DIC_ALT_CO2
$s3  HDIFN_DIC_ALT_CO2
$s3  HDIFB_DIC_ALT_CO2
EOF
  endif

  # include these budget check fields when doing development
  if ( 1 ) then
    cat >> $CASEROOT/Buildconf/popconf/ecosys_tavg_contents << EOF
$s1  Jint_Ctot
$s1  Jint_100m_Ctot
$s1  Jint_Ntot
$s1  Jint_100m_Ntot
$s1  Jint_Ptot
$s1  Jint_100m_Ptot
$s1  Jint_Sitot
$s1  Jint_100m_Sitot
$s1  Jint_Fetot
$s1  Jint_100m_Fetot
EOF
  endif

#1  Jint_PO4
#1  Jint_NO3
#1  Jint_SiO3
#1  Jint_NH4
#1  Jint_Fe
#1  Jint_O2
#1  Jint_DIC
#1  Jint_ALK
#1  Jint_DOC
#1  Jint_spC
#1  Jint_spChl
#1  Jint_spCaCO3
#1  Jint_diatC
#1  Jint_diatChl
#1  Jint_zooC
endif
