#!/bin/csh -f

#===============================================================================
#  ocn.ecosys.setup.csh : perform setup tasks for ecosys module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       pop_in_filename
#    tavg_contents  tavg_contents_filename
#    prestage       res_dpt_dir res_indpt_dir
#    document       pop2_document_files
#    ccsm_prestage  ccsm_prestage_file
#
#  CVS:$Id: ocn.ecosys.setup.csh 2448 2006-11-14 21:01:38Z njn01 $ 
#  CVS:$Name$ 
#
#===============================================================================

if ($#argv < 1) then
   echo ocn.ecosys.setup.csh : command argument missing
   exit 1
endif

#===============================================================================
#  set IC and forcing files
#===============================================================================

set relpath = ocn/pop/$OCN_GRID_INTERNAL
if ($OCN_GRID_INTERNAL == gx3v5) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx3v5_20051214.nc
   set ALK_scale_factor = 1.0
   set DIC_scale_factor = 1.0
   set O2_scale_factor  = 1.0
   set DST_file = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx3v5_20040426.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_gx3v5_20070521.nc

else if ($OCN_GRID_INTERNAL == gx1v5) then

   set IC_file  = $DIN_LOC_ROOT/$relpath/ic/ecosys_jan_IC_gx1v5_20070529.nc
   set ALK_scale_factor = 1.025
   set DIC_scale_factor = 1.025
   set O2_scale_factor  = 44.66
   set DST_file = $DIN_LOC_ROOT/$relpath/forcing/dst79gnx_gx1v5_20070329.nc
   set fesed_file = $DIN_LOC_ROOT/$relpath/forcing/fesedflux_gx1v5_20070518.nc

else

   echo $OCN_GRID_INTERNAL not supported by ecosystem module
   exit 2

endif

if (($RUN_TYPE == startup) && ($CONTINUE_RUN == FALSE)) then
  if ($OCN_PRESTAGE == TRUE) then
   set IC_file_nml = $INPUT/$IC_file:t
  else
   set IC_file_nml = $IC_file
  endif
else
   set IC_file_nml = same_as_TS
endif

if ($OCN_PRESTAGE == TRUE) then
  set DST_file_nml   = $INPUT/$DST_file:t
  set fesed_file_nml = $INPUT/$fesed_file:t
else
  set DST_file_nml   = $DST_file
  set fesed_file_nml = $fesed_file
endif

#===============================================================================

set command = $1

if ($command == set_nt) then

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : setting nt
   echo -----------------------------------------------------------------

   if ($#argv < 2) then
      echo nt_filename argument missing
      exit 3
   endif

   set nt_filename = $2

   if !(-f $nt_filename) then
      echo nt_filename = $nt_filename does not exist
      exit 3
   endif

   @ nt_in = `cat $nt_filename`
   @ nt_in += 24
   echo $nt_in > $nt_filename

else if ($command == namelist) then

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : setting namelist options
   echo -----------------------------------------------------------------

   if ($#argv < 2) then
      echo pop_in_filename argument missing
      exit 4
   endif

   set pop_in_filename = $2

   if !(-f $pop_in_filename) then
      echo pop_in_filename = $pop_in_filename does not exist
      exit 4
   endif

   if ($OCN_CO2_TYPE == constant) then
      set atm_co2_opt = const
   else if ($OCN_CO2_TYPE == prognostic) then
      set atm_co2_opt = model
   else
      echo error specifying atm_co2_opt in ocn.ecosys.setup.csh
      echo unknown OCN_CO2_TYPE:  $OCN_CO2_TYPE
      exit 4
   endif

   if (($RUN_TYPE == startup) && ($CONTINUE_RUN == FALSE)) then
      set use_nml_surf_vals = .true.
   else
      set use_nml_surf_vals = .false.
   endif

   set init_ecosys_option = $RUN_TYPE
   if ($CONTINUE_RUN == TRUE) set init_ecosys_option = continue

   cat >> $pop_in_filename << EOF

&ecosys_nml
   init_ecosys_option              = '$init_ecosys_option'
   init_ecosys_init_file           = '$IC_file_nml'
   init_ecosys_init_file_fmt       = 'nc'
   tracer_init_ext(1)%mod_varname = 'ALK'
   tracer_init_ext(1)%scale_factor = $ALK_scale_factor
   tracer_init_ext(2)%mod_varname = 'DIC'
   tracer_init_ext(2)%scale_factor = $DIC_scale_factor
   tracer_init_ext(3)%mod_varname = 'O2'
   tracer_init_ext(3)%scale_factor = $O2_scale_factor
   lflux_gas_o2                    = .true.
   lflux_gas_co2                   = .true.
   atm_co2_opt                     = '$atm_co2_opt'
   atm_co2_const                   = $CCSM_CO2_PPMV
   ecosys_tadvect_ctype            = 'base_model'
   gas_flux_forcing_opt            = 'model'
   lmarginal_seas                  = .true.
   lsource_sink                    = .true.
   comp_surf_avg_freq_opt          = 'never'
   comp_surf_avg_freq              = 1
   use_nml_surf_vals               = $use_nml_surf_vals
   surf_avg_dic_const              = 1944.0
   surf_avg_alk_const              = 2225.0
   ecosys_qsw_distrb_const         = .false.
!  iron_dust_flx_data_type         = 'monthly'
   dust_flux_input%filename        = '$DST_file_nml'
   dust_flux_input%file_fmt        = 'nc'
   dust_flux_input%file_varname    = 'DSTSF'
   dust_flux_input%scale_factor    = 1.0e-1    ! kg/m^2/sec -> g/cm^2/sec
   iron_flux_input%filename        = '$DST_file_nml'
   iron_flux_input%file_fmt        = 'nc'
   iron_flux_input%file_varname    = 'DSTSF'
   iron_flux_input%scale_factor    = 6.2668e4  ! kg/m^2/sec -> nmol/cm^2/sec, 3.5% iron by weight
   fesedflux_input%filename       = '$fesed_file_nml'
   fesedflux_input%file_varname   = 'FESEDFLUXIN'
   fesedflux_input%file_fmt        = 'nc'
   fesedflux_input%scale_factor   = 1.1574e-6 ! umolFe/m2/day -> nmolFe/cm2/s
!   dustflx_daily_input%interp_type    = 'linear'
!   dustflx_daily_input%filename       = '/fiji/home/ivan/data/Mahowald/daily/dstgnx'
!   dustflx_daily_input%data_renorm(1) = 1.e-1    ! kg/m^2/sec -> g/cm^2/sec
!   ironflx_daily_input%interp_type    = 'linear'
!   ironflx_daily_input%filename       = '/fiji/home/ivan/data/Mahowald/daily/dstgnx'
!   ironflx_daily_input%data_renorm(1) = 6.2668e4 ! kg/m^2/sec->nmol/cm^2/sec, 3.5% Fe
/

&ecosys_parms_nml
/
EOF

else if ($command == tavg_contents) then

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : setting tavg_contents variables
   echo -----------------------------------------------------------------

   if ($#argv < 2) then
      echo tavg_contents_filename argument missing
      exit 5
   endif

   set tavg_contents_filename = $2

   if !(-f $tavg_contents_filename) then
      echo tavg_contents_filename = $tavg_contents_filename does not exist
      exit 5
   endif

   cat >> $tavg_contents_filename << EOF
ECOSYS_ATM_PRESS
ECOSYS_IFRAC
ECOSYS_XKW
SCHMIDT_O2
SCHMIDT_CO2
IRON_FLUX
PH
O2SAT
STF_O2
CO2STAR
DCO2STAR
pCO2SURF
DpCO2
FG_CO2
ATM_CO2
FvPER_DIC
FvICE_DIC
FvPER_ALK
FvICE_ALK
Jint_PO4
Jint_NO3
Jint_SiO3
Jint_NH4
Jint_Fe
Jint_O2
Jint_DIC
Jint_ALK
Jint_DOC
Jint_spC
Jint_spChl
Jint_spCaCO3
Jint_diatC
Jint_diatChl
Jint_zooC
PO4
NO3
SiO3
NH4
Fe
O2
DIC
ALK
DOC
spC
spChl
spCaCO3
diatC
diatChl
zooC
spFe
diatSi
diatFe
diazC
diazChl
diazFe
DON
DOFe
DOP
graze_sp
graze_diat
graze_diaz
sp_agg
diat_agg
photoC_sp
photoC_diat
photoC_diaz
Fe_scavenge
Fe_scavenge_rate
CaCO3_form
diaz_Nfix
bSi_form
NITRIF
DENITRIF
POC_PROD
CaCO3_PROD
SiO2_PROD
P_iron_PROD
POC_FLUX_IN
CaCO3_FLUX_IN
SiO2_FLUX_IN
P_iron_FLUX_IN
PAR_avg
sp_Fe_lim
diat_Fe_lim
diaz_Fe_lim
sp_N_lim
diat_N_lim
sp_PO4_lim
diat_PO4_lim
diaz_P_lim
diat_SiO3_lim
sp_light_lim
diat_light_lim
diaz_light_lim
DOC_prod
DON_prod
DOFe_prod
DOP_prod
sp_loss
diat_loss
zoo_loss
diaz_loss
EOF

# dust_FLUX_IN
# DOC_remin
# DON_remin
# DOFe_remin
# DOP_remin
# photoFe_diaz
# photoFe_diat
# photoFe_sp

else if ($command == prestage) then

 if ($OCN_PRESTAGE == TRUE) then

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : prestaging data files
   echo -----------------------------------------------------------------

   if ($#argv < 1) then
      echo input data directories missing
      exit 6
   endif

   if ($CONTINUE_RUN == FALSE) then
      \cp -f $IC_file    $IC_file_nml  || exit 6
   endif

   \cp -f $DST_file $DST_file_nml || exit 6
   \cp -f $fesed_file $fesed_file_nml || exit 6

 else

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : data files are not being prestaged
   echo -----------------------------------------------------------------

 endif

else if ($command == document) then

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : documenting inputdata files
   echo -----------------------------------------------------------------

   if ($#argv < 2) then
      echo error documenting input files
      exit 7
   endif

   set pop2_document_files   = $2

  \ls -la $IC_file    >> $pop2_document_files || exit 7
  \ls -la $DST_file   >> $pop2_document_files || exit 7
  \ls -la $fesed_file >> $pop2_document_files || exit 7

else if ($command == ccsm_prestage) then

  #echo -----------------------------------------------------------------
  #echo ocn.ecosys.setup.csh : writing ccsm prestaging information
  #echo -----------------------------------------------------------------

   if ($#argv < 2) then
      echo error documenting input files
      exit 8
   endif

   set ccsm_prestage_file = $2

   #----------------------------------------------------------------------
   # Note: the following info is used only by the ccsm check_input_data
   # data-checking utility. Everything to the right of the equals sign
   # is ignored by check_input_data.
   #----------------------------------------------------------------------
   if (-f $ccsm_prestage_file) then
     if ($IC_file:h == $DIN_LOC_ROOT/$relpath/ic) then
       echo "    set ecosys_IC_file    = "\$DIN_LOC_ROOT/$relpath/ic/$IC_file:t         >> $ccsm_prestage_file || exit 8
     endif
     if ($DST_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set ecosys_DST_file   = "\$DIN_LOC_ROOT/$relpath/forcing/$DST_file:t   >> $ccsm_prestage_file || exit 8
     endif
     if ($fesed_file:h == $DIN_LOC_ROOT/$relpath/forcing) then
       echo "    set ecosys_fesed_file = "\$DIN_LOC_ROOT/$relpath/forcing/$fesed_file:t >> $ccsm_prestage_file || exit 8
     endif
   endif

else

   echo -----------------------------------------------------------------
   echo ocn.ecosys.setup.csh : unrecognized command $command
   echo -----------------------------------------------------------------

   exit 9

endif
