#! /bin/csh -f

#=========================================================================
# Purpose: create tavg_nml for inclusion in the file pop2_in
#          by combining base-model and extra-tracer-model information
#=========================================================================

#----------------------------------------------------------------------------
#  initialize
#----------------------------------------------------------------------------
  @ ntracer_stream                         = 1
  set ltavg_ignore_extra_streams = .false.
  set tavg_freq_opt_values                 = ( )
  set tavg_freq_values                     = ( )
  set tavg_stream_filestrings              = ( )
  set tavg_file_freq_opt                   = ( )
  set tavg_file_freq_values                = ( )
  set tavg_start_opt_values                = ( )
  set tavg_start_values                    = ( )
  set tavg_fmt_in_values                   = ( )
  set tavg_fmt_out_values                  = ( )
  set ltavg_has_offset_date_values         = ( )
  set tavg_offset_year_values              = ( )
  set tavg_offset_month_values             = ( )
  set tavg_offset_day_values               = ( )
  set ltavg_one_time_header                = ( )

#----------------------------------------------------------------------------
# parse the base-model tavg_nml, located in $POP2_TAVG_NML_BASE, and 
#  extract base-model values
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#  extract single-valued variable settings
#----------------------------------------------------------------------------
@ n_tavg_streams                = `grep "^ *n_tavg_streams"              $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `

set ltavg_temp  = `grep "^ *ltavg_ignore_extra_streams"   $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `

if ($ltavg_temp == .true. || $ltavg_temp == .false.) then
set ltavg_ignore_extra_streams  = $ltavg_temp
endif

set ltavg_streams_index_present = `grep "^ *ltavg_streams_index_present" $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `

set ltavg_nino_diags_requested  = `grep "^ *ltavg_nino_diags_requested"  $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `

set tavg_contents_filename      = `grep "^ *tavg_contents"  $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `
set tavg_infile                 = `grep "^ *tavg_infile"    $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `
set tavg_outfile                = `grep "^ *tavg_outfile"   $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `



#----------------------------------------------------------------------------
#  extract vector-valued variable settings
#----------------------------------------------------------------------------
set VARS = (tavg_freq_opt tavg_freq tavg_stream_filestrings tavg_file_freq_opt tavg_file_freq tavg_start_opt tavg_start tavg_fmt_in tavg_fmt_out ltavg_has_offset_date tavg_offset_years tavg_offset_months tavg_offset_days ltavg_one_time_header )

foreach VAR ($VARS)
  set opts = `grep "^ *$VAR " $POP2_TAVG_NML_BASE | awk -F= '{print $2}' `
  if ($VAR == tavg_freq_opt) then
                                             set tavg_freq_opt_values    = ($tavg_freq_opt_values  $opts)
  else if ($VAR == tavg_freq) then
                                             set tavg_freq_values        = ($tavg_freq_values  $opts)
  else if ($VAR == tavg_stream_filestrings) then
                                             set tavg_stream_filestrings = ($tavg_stream_filestrings  $opts)
  else if ($VAR == tavg_file_freq_opt) then
                                             set tavg_file_freq_opt      = ($tavg_file_freq_opt  $opts)
  else if ($VAR == tavg_file_freq) then
                                             set tavg_file_freq_values   = ($tavg_file_freq_values  $opts)
  else if ($VAR == tavg_start_opt) then
                                             set tavg_start_opt_values   = ($tavg_start_opt_values  $opts)
  else if ($VAR == tavg_start) then
                                             set tavg_start_values       = ($tavg_start_values  $opts)
  else if ($VAR == tavg_fmt_in) then
                                             set tavg_fmt_in_values      = ($tavg_fmt_in_values  $opts)
  else if ($VAR == tavg_fmt_out) then
                                             set tavg_fmt_out_values     = ($tavg_fmt_out_values  $opts)
  else if ($VAR == ltavg_has_offset_date) then
                                             set ltavg_has_offset_date_values = ($ltavg_has_offset_date_values  $opts)
  else if ($VAR == tavg_offset_years) then
                                             set tavg_offset_year_values  = ($tavg_offset_year_values  $opts)
  else if ($VAR == tavg_offset_months) then
                                             set tavg_offset_month_values = ($tavg_offset_month_values  $opts)
  else if ($VAR == tavg_offset_days) then
                                             set tavg_offset_day_values   = ($tavg_offset_day_values  $opts)
  else if ($VAR == ltavg_one_time_header) then
                                             set ltavg_one_time_header    = ($ltavg_one_time_header  $opts)
  else
    echo pop2_in_build.csh failure in $module loop 
    echo    undefined VAR = $VAR
    exit -999
  endif
end # VAR

#----------------------------------------------------------------------------
#   now, define the tavg_nml settings for all of the extra-tracer modules.
#   This involves the creation and subsequent reading of a tracer-module file 
#----------------------------------------------------------------------------
set ocn_tracers = (`echo $OCN_TRACER_MODULES`)
set tracer_stream_numbers = ( )
foreach module ($ocn_tracers)
    if (-f ${MY_PATH}/ocn.${module}.setup.csh) then
       set setup_path = ${MY_PATH}
    else if (-f $SRCDIR/input_templates/ocn.${module}.setup.csh) then
       set setup_path = $SRCDIR/input_templates
    else
       echo error in pop2_in_build.csh unknown tracer: $module
       exit -3
    endif

    ${setup_path}/ocn.${module}.setup.csh set_tavg_nml || exit $status

    @ ntracer_stream = `grep "^ *n_tavg_streams_tracer" $POP2_DOCDIR/$module.tavg | awk -F= '{print $2}' `

  if ($ntracer_stream == 0) then
    @ module_stream  = $n_tavg_streams
  else
    #------------------------------------------------------------------------------------
    # increment number of streams and identify the starting stream for each tracer module
    #------------------------------------------------------------------------------------
    @ module_stream  = $n_tavg_streams + 1
    @ n_tavg_streams = $n_tavg_streams + $ntracer_stream



#------------------------------------------------------------------------------------
# parse the information in the $POP2_DOCDIR/$module.tavg file (created by the setup.csh scripts)
#------------------------------------------------------------------------------------
    foreach VAR ($VARS)
      set opts = `grep "^ *$VAR " $POP2_DOCDIR/$module.tavg | awk -F= '{print $2}' `
      if ($VAR == tavg_freq_opt) then
        set tavg_freq_opt_values = ($tavg_freq_opt_values  $opts)
      else if ($VAR == tavg_freq) then
        set tavg_freq_values = ($tavg_freq_values  $opts)
      else if ($VAR == tavg_stream_filestrings) then
        set tavg_stream_filestrings = ($tavg_stream_filestrings  $opts)
      else if ($VAR == tavg_file_freq_opt) then
        set tavg_file_freq_opt = ($tavg_file_freq_opt  $opts)
      else if ($VAR == tavg_file_freq) then
        set tavg_file_freq_values = ($tavg_file_freq_values  $opts)
      else if ($VAR == tavg_start_opt) then
        set tavg_start_opt_values = ($tavg_start_opt_values  $opts)
      else if ($VAR == tavg_start) then
        set tavg_start_values = ($tavg_start_values  $opts)
      else if ($VAR == tavg_fmt_in) then
        set tavg_fmt_in_values = ($tavg_fmt_in_values  $opts)
      else if ($VAR == tavg_fmt_out) then
        set tavg_fmt_out_values = ($tavg_fmt_out_values  $opts)
      else if ($VAR == ltavg_has_offset_date) then
        set ltavg_has_offset_date_values = ($ltavg_has_offset_date_values  $opts)
      else if ($VAR == tavg_offset_years) then
        set tavg_offset_year_values = ($tavg_offset_year_values  $opts)
      else if ($VAR == tavg_offset_months) then
        set tavg_offset_month_values = ($tavg_offset_month_values  $opts)
      else if ($VAR == tavg_offset_days) then
        set tavg_offset_day_values = ($tavg_offset_day_values  $opts)
      else if ($VAR == ltavg_one_time_header) then
        set ltavg_one_time_header = ($ltavg_one_time_header  $opts)
      else
        echo pop2_in_build.csh failure in $module loop 
        echo    undefined VAR = $VAR
        exit -999
      endif
    end # VAR
   endif
    #------------------------------------------------------------------------------------
    # store information about stream numbers for later use by ocn.*.setup.csh scripts
    # see section below, after the complete specification of pop2_in
    #------------------------------------------------------------------------------------
   set tracer_stream_numbers = ($tracer_stream_numbers $module_stream)
   #---------
   # clean up
   #---------
   ### NEVER rm $POP2_DOCDIR/$module.tavg
end # loop over tracer setup scripts

#---------
# clean up
#---------
### NEVER rm $POP2_TAVG_NML_BASE
cat >&! $POP2_TAVG_NML << EOF2

##########################################################################
 Combined tavg stream information from base model and extra-tracer models
##########################################################################

&tavg_nml
   n_tavg_streams                       = $n_tavg_streams
   ltavg_ignore_extra_streams           = $ltavg_ignore_extra_streams
   ltavg_streams_index_present          = $ltavg_streams_index_present
   tavg_freq_opt                        = $tavg_freq_opt_values
   tavg_freq                            = $tavg_freq_values
   tavg_file_freq_opt                   = $tavg_file_freq_opt
   tavg_file_freq                       = $tavg_file_freq_values
   tavg_stream_filestrings              = $tavg_stream_filestrings
   tavg_start_opt                       = $tavg_start_opt_values
   tavg_start                           = $tavg_start_values
   tavg_fmt_in                          = $tavg_fmt_in_values
   tavg_fmt_out                         = $tavg_fmt_out_values
   tavg_contents                        = $tavg_contents_filename
   ltavg_nino_diags_requested           = $ltavg_nino_diags_requested
   tavg_infile                          = $tavg_infile 
   tavg_outfile                         = $tavg_outfile 
   ltavg_has_offset_date                = $ltavg_has_offset_date_values
   tavg_offset_years                    = $tavg_offset_year_values
   tavg_offset_months                   = $tavg_offset_month_values
   tavg_offset_days                     = $tavg_offset_day_values
   ltavg_one_time_header                = $ltavg_one_time_header
/

EOF2

