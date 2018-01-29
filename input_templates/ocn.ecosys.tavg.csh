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

if ( -f $CASEROOT/SourceMods/src.pop/marbl_diagnostics ) then
  set MARBL_args = "-i $CASEROOT/SourceMods/src.pop/marbl_diagnostics"
else
  set MARBL_args = "-i $CASEBUILD/popconf/marbl_diagnostics"
endif
set MARBL_args = "$MARBL_args -t $CASEBUILD/popconf/ecosys_tavg_contents -o $CASEBUILD/popconf/marbl_diagnostics_operators --low_frequency_stream $s3 --medium_frequency_stream $s1 --high_frequency_stream $s2 --append True"

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
$s1  DIC
$s1  J_DIC
$s1  ALK
$s1  DOC
$s1  zooC
$s1  DON
$s1  DOP
$s1  DONr
$s1  DOPr
$s1  DOCr
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
$s2  ${autotroph}C_zint_100m
$s2  ${autotroph}Chl_SURF
EOF
if ($lvariable_PtoC == ".true.") then
    cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  ${autotroph}P
EOF
endif
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
EOF
  end

  if ($lecosys_tavg_alt_co2 == ".true.") then
cat >> $CASEBUILD/popconf/ecosys_tavg_contents << EOF
$s1  DIC_ALT_CO2
$s1  ALK_ALT_CO2
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

endif

# Add MARBL diagnostics to tavg_contents
$CASEBUILD/popconf/MARBL_diags_to_tavg.py $MARBL_args
if ($status != 0) then
  echo "ERROR in MARBL_diags_to_tavg.py"
  exit 1
endif
