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

set lecosys_tavg_all     = $2
set lecosys_tavg_alt_co2 = $3
set ladjust_bury_coeff   = $4
set lvariable_PtoC       = $5
set tracer_restore_vars  = "PO4 NO3 SiO3 ALK"

set MARBL_args = "-i $CASEBUILD/popconf/marbl_diagnostics -t $CASEBUILD/popconf/ecosys_tavg_contents -o $CASEBUILD/popconf/marbl_diagnostics_operators --low_frequency_stream $s3 --medium_frequency_stream $s1 --high_frequency_stream $s2 --append True"

if ($lecosys_tavg_all == ".true.") then
  set MARBL_args = "$MARBL_args --lMARBL_tavg_all True --lMARBL_tavg_alt_co2 True"
else if ($lecosys_tavg_alt_co2 == ".true.") then
  set MARBL_args = "$MARBL_args --lMARBL_tavg_alt_co2 True"
endif

if ($lecosys_tavg_all == ".false.") then

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

if ($lecosys_tavg_all == ".true.") then

  cat >! $CASEBUILD/popconf/ecosys_tavg_contents << EOF
#  GENERAL INTERIOR DIAGNOSTICS
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
1  Lig_prod
1  Lig_loss
1  Lig_scavenge
1  Fefree
1  Lig_photochem
1  Lig_deg
#  PARTICULATE
1  calcToSed
1  calcToSed_ALT_CO2
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
1  POP_FLUX_IN
1  POP_PROD
1  POP_REMIN
1  POP_REMIN_PO4
1  PON_REMIN_NH4
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
1  FINE_DUST_FLUX_CPL
1  COARSE_DUST_FLUX_CPL
1  BLACK_CARBON_FLUX_CPL
#  DUPLICATE TAVG VAR
1  STF_O2_2
EOF

  cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
#  bury coefficients (only available if ladjust_bury_coeff==.true.)
EOF
if ($ladjust_bury_coeff == ".true.") then
  cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
1  MARBL_rmean_glo_scalar_POC_bury_coeff
1  MARBL_rmean_glo_scalar_POP_bury_coeff
1  MARBL_rmean_glo_scalar_bSi_bury_coeff
EOF
endif

echo "#  River Fluxes" >> $CASEBUILD/popconf/ecosys_tavg_contents
  set tracer_list = ( NO3 PO4 DON DONr DOP DOPr SiO3 Fe DIC ALK DOC DOCr DIC_ALT_CO2 ALK_ALT_CO2 )
  foreach tracer ( $tracer_list )
    echo "1  ${tracer}_RIV_FLUX" >> $CASEBUILD/popconf/ecosys_tavg_contents
  end

  # interior autotroph fields
echo "#  AUTOTROPH TRACER FIELDS" >> $CASEBUILD/popconf/ecosys_tavg_contents
  foreach autotroph ( sp diat diaz )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
#  ${autotroph} TRACER FIELDS
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
if ($lvariable_PtoC == ".true.") then
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
1  ${autotroph}_Qp
EOF
endif
  end
  cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
1  sp_CaCO3_form_zint
1  sp_CaCO3_form
1  diat_SiO3_lim
1  diat_bSi_form
1  diaz_Nfix
1  CaCO3_form_zint
1  bSi_form
1  CaCO3_form
1  Nfix
EOF

echo "#  ZOOPLANKTON TRACER FIELDS" >> $CASEBUILD/popconf/ecosys_tavg_contents
  foreach zooplankton ( zoo )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
#  ${zooplankton} TRACER FIELDS
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

echo "#  TRACER FIELDS" >> $CASEBUILD/popconf/ecosys_tavg_contents
  set tracer_list = ( PO4 NO3 SiO3 NH4 Fe Lig O2 DIC DIC_ALT_CO2 ALK ALK_ALT_CO2 DOC DON DOCr \
                   DOP DOPr DONr zooC spChl spC spFe spCaCO3 diatChl diatC \
                   diatFe diatSi diazChl diazC diazFe )
  if ($lvariable_PtoC == ".true.") then
    set tracer_list = ( $tracer_list spP diatP diazP )
  endif
  foreach tracer ( $tracer_list )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
#  ${tracer} TRACER FIELDS
1  ${tracer}
1  STF_${tracer}
1  J_${tracer}
1  ${tracer}_RESTORE_TEND
EOF
  end
  cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
1  FvPER_DIC
1  FvPER_ALK
1  FvPER_DIC_ALT_CO2
1  FvPER_ALK_ALT_CO2
1  FvICE_DIC
1  FvICE_ALK
EOF

else # ecosys not in lecosys_tavg_all mode
  cat >! $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  FINE_DUST_FLUX_CPL
$s1  COARSE_DUST_FLUX_CPL
$s1  BLACK_CARBON_FLUX_CPL
$s1  STF_ALK
$s1  STF_O2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  PO4
$s1  NO3
$s1  SiO3
$s1  NH4
$s1  Fe
$s1  Lig
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
$s1  Lig_prod
$s1  Lig_loss
$s1  Lig_scavenge
$s1  Fefree
$s1  Lig_photochem
$s1  Lig_deg
$s1  bSi_form
$s1  NITRIF
$s1  DENITRIF
$s1  CaCO3_FLUX_IN
$s1  CaCO3_PROD
$s1  CaCO3_REMIN
$s1  dust_FLUX_IN
$s1  dust_REMIN
$s1  P_iron_FLUX_IN
$s1  P_iron_PROD
$s1  P_iron_REMIN
$s1  POC_FLUX_IN
$s1  POC_PROD
$s1  POC_REMIN
$s1  POC_REMIN_DIC
$s1  PON_REMIN_NH4
$s1  POP_FLUX_IN
$s1  POP_PROD
$s1  POP_REMIN
$s1  POP_REMIN_PO4
$s1  SiO2_FLUX_IN
$s1  SiO2_PROD
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

  if ($ladjust_bury_coeff == ".true.") then
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  MARBL_rmean_glo_scalar_POC_bury_coeff
$s1  MARBL_rmean_glo_scalar_POP_bury_coeff
$s1  MARBL_rmean_glo_scalar_bSi_bury_coeff
EOF
  endif

  # River Fluxes
  set tracer_list = ( NO3 PO4 DON DONr DOP DOPr SiO3 Fe DIC ALK DOC DOCr )
  if ($lecosys_tavg_alt_co2 == ".true.") then
    set tracer_list = ( $tracer_list DIC_ALT_CO2 ALK_ALT_CO2 )
  endif
  foreach tracer ( $tracer_list )
    echo "$s1  ${tracer}_RIV_FLUX" >> $CASEBUILD/popconf/ecosys_tavg_contents
  end

  # generic autotroph fields
  # skip N_lim for diaz
  foreach autotroph ( sp diat diaz )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
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
if ($lvariable_PtoC == ".true.") then
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}P
$s1  ${autotroph}_Qp
EOF
endif
    if !($autotroph == diaz) then
      cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}_N_lim
EOF
    endif
  end

  # Nfix terms from N fixers 
  foreach autotroph ( diaz )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}_Nfix
EOF
  end

  # CaCO3 terms from calcifiers 
  foreach autotroph ( sp )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}CaCO3
$s2  ${autotroph}CaCO3_zint_100m
EOF
  end

  # Si terms from silicifiers
  foreach autotroph ( diat )
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}Si
$s1  ${autotroph}_SiO3_lim
EOF
  end

  # restoring terms for tracers that have restoring enabled
  foreach tracer_restore_var ( `echo $tracer_restore_vars | tr ',' ' '` )
    echo "$s1  ${tracer_restore_var}_RESTORE_TEND" >> $CASEBUILD/popconf/ecosys_tavg_contents
  end

  if ($lecosys_tavg_alt_co2 == ".true.") then
cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  DIC_ALT_CO2
$s1  ALK_ALT_CO2
$s1  calcToSed_ALT_CO2
$s1  CaCO3_ALT_CO2_REMIN
$s1  CaCO3_ALT_CO2_FLUX_IN
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
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
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

# Add MARBL diagnostics to tavg_contents
$CASEBUILD/popconf/MARBL_diags_to_tavg.py $MARBL_args

if ($lecosys_tavg_all == ".true.") then
  cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
1  FG_CO2_2
EOF
endif
