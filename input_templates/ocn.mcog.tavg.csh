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

set lmcog_debug = $2
set mcog_ncols  = $3
set mcog_nbins  = $4

if ($lmcog_debug == ".true.") then
   set stream_normal_var = "2"
   set stream_debug_var  = "2"
else
   set stream_normal_var = "1"
   set stream_debug_var  = "!"
endif

rm -f $CASEROOT/Buildconf/popconf/mcog_tavg_contents
touch $CASEROOT/Buildconf/popconf/mcog_tavg_contents

@ nbin = 1
while ($nbin <= $mcog_nbins)
   set nn = `printf "%02d" $nbin`
   echo "$stream_debug_var  FRAC_BIN_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   echo "$stream_normal_var  FRACR_BIN_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   echo "$stream_debug_var  QSW_RAW_BIN_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   echo "$stream_normal_var  QSW_BIN_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   @ nbin++
end

@ ncol = 1
while ($ncol <= $mcog_ncols)
   set nn = `printf "%02d" $ncol`
   echo "$stream_debug_var  FRAC_COL_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   echo "$stream_debug_var  FRACR_COL_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   echo "$stream_debug_var  QSW_RAW_COL_$nn" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
   @ ncol++
end

echo "$stream_debug_var  QSW_RAW_COL_DAGG" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
echo "$stream_debug_var  QSW_RAW_BIN_DAGG" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
echo "$stream_debug_var  FRAC_ADJUST_FACT" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
echo "$stream_debug_var  FRACR_ADJUST_FACT" >> $CASEROOT/Buildconf/popconf/mcog_tavg_contents
