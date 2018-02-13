#!/bin/csh -f

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#
# To change what ecosystem diagnostics are output in tavg, copy
#    Buildconf/popconf/ecosys_diagnostics
# to SourceMods/src.pop and make any changes desired.
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

rm -f $CASEBUILD/popconf/ecosys_tavg_contents

# Command line arguments for call to MARBL_diags_to_tavg.py
# Which will read ecosys_diagnostics and produce ecosys_tavg_contents
set MARBL_args = "-t $CASEBUILD/popconf/ecosys_tavg_contents"
set MARBL_args = "$MARBL_args -d $CASEBUILD/popconf/marbl_diagnostics_list"
set MARBL_args = "$MARBL_args -o $CASEBUILD/popconf/marbl_diagnostics_operators"
# MARBL recommends one of three temporal averaging frequencies: low, medium, or high
# By default, POP uses low -> annual, medium -> monthly, and high -> daily
# You can change the low and high frequency streams in the block above
# To change the medium frequency stream, edit ocn.base.tavg.csh
set MARBL_args = "$MARBL_args --low_frequency_stream $s3 "
set MARBL_args = "$MARBL_args --medium_frequency_stream $s1"
set MARBL_args = "$MARBL_args --high_frequency_stream $s2"

if ($lecosys_tavg_all == ".true.") then
  set MARBL_args = "$MARBL_args --lMARBL_tavg_all True --lMARBL_tavg_alt_co2 True"
else if ($lecosys_tavg_alt_co2 == ".true.") then
  set MARBL_args = "$MARBL_args --lMARBL_tavg_alt_co2 True"
endif

# Add POP-based ecosys diagnostics to tavg_contents
if ( -f $CASEROOT/SourceMods/src.pop/ecosys_diagnostics ) then
  set MARBL_args_filename = "-i $CASEROOT/SourceMods/src.pop/ecosys_diagnostics"
else
  set MARBL_args_filename = "-i $CASEBUILD/popconf/ecosys_diagnostics"
endif
$POPROOT/MARBL_scripts/MARBL_diags_to_tavg.py $MARBL_args $MARBL_args_filename
if ($status != 0) then
  echo "ERROR in MARBL_diags_to_tavg.py"
  exit 1
endif
