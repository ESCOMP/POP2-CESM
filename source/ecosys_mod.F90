! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_mod

  !BOP
  ! !MODULE: ecosys_mod
  !
  ! !DESCRIPTION:
  !
  !  Multispecies ecosystem based on Doney et al. 1996, Moore et al., 2002
  !  Based on POP Global NCAR Nitrogen Ecosystem Model
  !  version 0.0 (June 15th, 1998) from S.C. Doney.
  !  Based on Doney et al., 1996 model.
  !  Climate and Global Dynamics, NCAR
  !  (doney@whoi.edu)
  !
  !  Version 1.0
  !  Multispecies, multiple limiting nutrient version of ecosystem
  !  based on mixed layer model of Moore et al.(2002).  Implemented here with
  !  fixed elemental ratios and including only the diatoms and small
  !  phytoplankton, with a parameterization of calcification,
  !  by Keith Lindsay and Keith Moore, Fall 2001 - Spring 2002.
  !  Calcification parameterization based on Moore et al. 2002.
  !
  !  Version 2.0, January 2003
  !    Adds diazotrophs as a phytoplankton group, (based on Moore et al., 2002a)
  !    Allows for variable fe/C for all phytoplankton groups
  !     Allows for variable si/C for the diatoms
  !     Adds explicit tracers for DON, DOP, DOFe
  !     variable remin length scale for detrital soft POM and bSi f(temperature)
  !     Extensive modifications to iron scavenging parameterization
  !     Addition of a sedimentary dissolved iron source,
  !        (implemented in ballast code as excess remin in bottom cell)
  !        coded by J.K. Moore, (jkmoore@uci.edu)
  !
  !   Version 2.01. March 2003
  !     corrected O2 bug
  !     corrected grazing parameter z_grz bug at depth
  !     dust dissolution at depth releases iron,
  !     increased length scale for dust diss., increased hard fraction dust
  !     no deep ocean reduction in scavenging rates,
  !     increase bSi OC/ballast ratio 0.3 -> 0.35,
  !     corrected bug in diazotroph photoadaptation, and diat and sp adapatation
  !
  !   Version 2.02.
  !     corrected bug in Fe_scavenge (units for dust), May 2003
  !     changed C/N/P ratios to 117/16/1 (Anderson & Sarmiento, 1994)
  !
  !   Version 2.03., July 2003
  !     Remin of DOM no longer temperature dependent,
  !     new iron scavenging parameterization added,
  !     some dissolution of hard fraction of ballast materials added
  !
  !   Version 2.1, September 2003
  !     modfied iron scavenging and dust dissolution at depth
  !
  !   Version 2.11, March 2004
  !     fixed bug in iron scavenging code, replace dust and POC flux_in w/ flux_out
  !
  !   Version 2.12, April 2004 - Final version for GBC paper revision,
  !     (Questions/comments, Keith Moore - jkmoore@uci.edu
  !
  !   References
  !   Doney, S.C., Glover, D.M., Najjar, R.G., 1996. A new coupled, one-dimensional
  !   biological-physical model for the upper ocean: applications to the JGOFS
  !   Bermuda Time-Series Study (BATS) site. Deep-Sea Res. II, 43: 591-624.
  !
  !   Moore, JK, Doney, SC, Kleypas, JA, Glover, DM, Fung, IY, 2002. An intermediate
  !   complexity marine ecosystem model for the global domain. Deep-Sea Res. II, 49:
  !   403-462.
  !
  !   Moore, JK, Doney, SC, Glover, DM, Fung, IY, 2002. Iron cycling and nutrient
  !   limitation patterns in surface waters of the world ocean. Deep-Sea Res. II,
  !   49: 463-507.

  ! !REVISION HISTORY:

  !  SVN:$Id:  $

  !-----------------------------------------------------------------------
  !  variables/subroutines/function used from other modules
  !  The following are used extensively in this ecosys, so are used at
  !  the module level. The use statements for variables that are only needed
  !  locally are located at the module subprogram level.
  !-----------------------------------------------------------------------

  ! !USES:

  use POP_KindsMod, only : POP_i4 ! pop, need similar in marbl
  use POP_KindsMod, only : POP_r8 ! pop, need similar in marbl

  use POP_ErrorMod, only : POP_ErrorSet ! pop, need similar in marbl
  use POP_ErrorMod, only : POP_Success  ! pop, need similar in marbl

  use POP_CommMod, only : POP_communicator ! pop

  use POP_GridHorzMod, only : POP_gridHorzLocCenter ! pop

  use POP_FieldMod, only : POP_fieldKindScalar

  use POP_HaloMod, only : POP_HaloUpdate ! pop

  use kinds_mod, only : log_kind
  use kinds_mod, only : int_kind
  use kinds_mod, only : r8
  use kinds_mod, only : char_len

  use constants, only : p5
  use constants, only : c0
  use constants, only : c1
  use constants, only : c2
  use constants, only : c10
  use constants, only : char_blank
  use constants, only : field_loc_center
  use constants, only : field_type_scalar
  use constants, only : mpercm
  use constants, only : hflux_factor
  use constants, only : blank_fmt
  use constants, only : delim_fmt
  use constants, only : ndelim_fmt
  use constants, only : xkw_coeff

  use communicate, only : master_task
  use communicate, only : my_task

  !   use global_reductions

  use blocks, only : block
  use blocks, only : nx_block
  use blocks, only : ny_block
  use blocks, only : get_block

  use domain_size, only : km
  use domain_size, only : nx_global
  use domain_size, only : ny_global
  use domain_size, only : max_blocks_clinic

  use domain, only : distrb_clinic
  use domain, only : nblocks_clinic
  use domain, only : POP_haloClinic
  use domain, only : blocks_clinic

  use exit_mod, only : exit_POP
  use exit_mod, only : sigAbort

  use prognostic, only : tracer_field
  use prognostic, only : curtime
  use prognostic, only : oldtime

  use grid, only : KMT
  use grid, only : DZT
  use grid, only : zt
  use grid, only : partial_bottom_cells
  use grid, only : n_topo_smooth
  use grid, only : REGION_MASK
  use grid, only : dz
  use grid, only : zw
  use grid, only : fill_points

  use io, only : datafile, data_set

  use io_types, only : io_dim
  use io_types, only : io_field_desc
  use io_types, only : nml_filename
  use io_types, only : nml_in
  use io_types, only : stdout
  use io_types, only : construct_io_dim
  use io_types, only : construct_io_field
  use io_types, only : add_attrib_file

  use io_tools, only : document

  use timers, only : timer_start
  use timers, only : timer_stop
  use timers, only : get_timer

  use passive_tracer_tools, only : ind_name_pair
  use passive_tracer_tools, only : tracer_read
  use passive_tracer_tools, only : field_exists_in_file
  use passive_tracer_tools, only : forcing_monthly_every_ts
  use passive_tracer_tools, only : comp_surf_avg
  use passive_tracer_tools, only : extract_surf_avg
  use passive_tracer_tools, only : file_read_tracer_block
  use passive_tracer_tools, only : init_forcing_monthly_every_ts
  use passive_tracer_tools, only : read_field
  use passive_tracer_tools, only : rest_read_tracer_block

  use named_field_mod, only : named_field_get
  use named_field_mod, only : named_field_register
  use named_field_mod, only : named_field_set

  use forcing_tools, only : find_forcing_times
  use forcing_tools, only : interpolate_forcing
  use forcing_tools, only : update_forcing_data

  use time_management, only : isecond
  use time_management, only : iminute
  use time_management, only : ihour
  use time_management, only : iday
  use time_management, only : imonth
  use time_management, only : iyear
  use time_management, only : thour00
  use time_management, only : freq_opt_never
  use time_management, only : freq_opt_nmonth
  use time_management, only : freq_opt_nyear
  use time_management, only : check_time_flag
  use time_management, only : init_time_flag
  use time_management, only : eval_time_flag

  use ecosys_constants, only : ecosys_tracer_cnt
  ! tavg_forcing diagnostics
  use ecosys_diagnostics_mod, only : forcing_diag_cnt
  use ecosys_diagnostics_mod, only : ECOSYS_IFRAC_diag_ind
  use ecosys_diagnostics_mod, only : ECOSYS_XKW_diag_ind
  use ecosys_diagnostics_mod, only : ECOSYS_ATM_PRESS_diag_ind
  use ecosys_diagnostics_mod, only : PV_O2_diag_ind
  use ecosys_diagnostics_mod, only : SCHMIDT_O2_diag_ind
  use ecosys_diagnostics_mod, only : O2SAT_diag_ind
  use ecosys_diagnostics_mod, only : O2_GAS_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : CO2STAR_diag_ind
  use ecosys_diagnostics_mod, only : DCO2STAR_diag_ind
  use ecosys_diagnostics_mod, only : pCO2SURF_diag_ind
  use ecosys_diagnostics_mod, only : DpCO2_diag_ind
  use ecosys_diagnostics_mod, only : PV_CO2_diag_ind
  use ecosys_diagnostics_mod, only : SCHMIDT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DIC_GAS_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : PH_diag_ind
  use ecosys_diagnostics_mod, only : ATM_CO2_diag_ind
  use ecosys_diagnostics_mod, only : CO2STAR_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DCO2STAR_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : pCO2SURF_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DpCO2_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : DIC_GAS_FLUX_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : PH_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : ATM_ALT_CO2_diag_ind
  use ecosys_diagnostics_mod, only : IRON_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DUST_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : NOx_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : NHy_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIN_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIP_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DoN_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DoNr_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOP_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOPr_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DSI_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DFE_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DIC_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : ALK_RIV_FLUX_diag_ind
  use ecosys_diagnostics_mod, only : DOC_RIV_FLUX_diag_ind
  ! Diagnostic subroutines
  use ecosys_diagnostics_mod, only : store_diagnostics_carbonate
  use ecosys_diagnostics_mod, only : store_diagnostics_nitrification
  use ecosys_diagnostics_mod, only : store_diagnostics_autotrophs
  use ecosys_diagnostics_mod, only : store_diagnostics_autotroph_sums
  use ecosys_diagnostics_mod, only : store_diagnostics_particulates
  use ecosys_diagnostics_mod, only : store_diagnostics_oxygen
  use ecosys_diagnostics_mod, only : store_diagnostics_photosynthetically_available_radiation
  use ecosys_diagnostics_mod, only : store_diagnostics_misc
  use ecosys_diagnostics_mod, only : store_diagnostics_zooplankton
  use ecosys_diagnostics_mod, only : store_diagnostics_dissolved_organic_matter
  use ecosys_diagnostics_mod, only : store_diagnostics_carbon_fluxes
  use ecosys_diagnostics_mod, only : store_diagnostics_nitrogen_fluxes
  use ecosys_diagnostics_mod, only : store_diagnostics_phosphorus_fluxes
  use ecosys_diagnostics_mod, only : store_diagnostics_silicon_fluxes

  use marbl_parms, only : marbl_params_init, marbl_params_print
  use marbl_parms, only : grz_fnc_michaelis_menten
  use marbl_parms, only : grz_fnc_sigmoidal
  use marbl_parms, only : f_qsw_par
  use marbl_parms, only : parm_Fe_bioavail
  use marbl_parms, only : dust_to_Fe
  use marbl_parms, only : parm_BSIbury
  use marbl_parms, only : parm_POMbury
  use marbl_parms, only : denitrif_C_N
  use marbl_parms, only : parm_Red_Fe_C
  use marbl_parms, only : Q
  use marbl_parms, only : Qp_zoo_pom
  use marbl_parms, only : spd
  use marbl_parms, only : parm_CaCO3_diss
  use marbl_parms, only : parm_POC_diss
  use marbl_parms, only : parm_SiO2_diss
  use marbl_parms, only : parm_scalelen_z
  use marbl_parms, only : parm_scalelen_vals
  use marbl_parms, only : caco3_poc_min
  use marbl_parms, only : CaCO3_sp_thres
  use marbl_parms, only : CaCO3_temp_thres1
  use marbl_parms, only : CaCO3_temp_thres2
  use marbl_parms, only : DOC_reminR
  use marbl_parms, only : DOFe_reminR
  use marbl_parms, only : DON_reminR
  use marbl_parms, only : DONr_reminR
  use marbl_parms, only : DONrefract
  use marbl_parms, only : DOP_reminR
  use marbl_parms, only : DOPr_reminR
  use marbl_parms, only : DOPrefract
  use marbl_parms, only : dps
  use marbl_parms, only : dust_fescav_scale
  use marbl_parms, only : f_graze_CaCO3_REMIN
  use marbl_parms, only : f_graze_si_remin
  use marbl_parms, only : f_graze_sp_poc_lim
  use marbl_parms, only : f_photosp_CaCO3
  use marbl_parms, only : fe_max_scale2
  use marbl_parms, only : Fe_scavenge_thres1
  use marbl_parms, only : parm_f_prod_sp_CaCO3
  use marbl_parms, only : parm_Fe_scavenge_rate0
  use marbl_parms, only : parm_kappa_nitrif
  use marbl_parms, only : parm_labile_ratio
  use marbl_parms, only : parm_nitrif_par_lim
  use marbl_parms, only : parm_o2_min
  use marbl_parms, only : parm_o2_min_delta
  use marbl_parms, only : parm_red_d_c_o2
  use marbl_parms, only : parm_red_d_c_o2_diaz
  use marbl_parms, only : parm_Remin_D_C_O2
  use marbl_parms, only : q_10
  use marbl_parms, only : QCaCO3_max
  use marbl_parms, only : Qfe_zoo
  use marbl_parms, only : r_Nfix_photo
  use marbl_parms, only : spc_poc_fac
  use marbl_parms, only : thres_z1
  use marbl_parms, only : thres_z2
  use marbl_parms, only : yps

  use registry, only : registry_match

  use named_field_mod, only : named_field_get_index

  use marbl_share_mod, only : autotrophs, autotroph_cnt
  use marbl_share_mod, only : zooplankton, zooplankton_cnt
  use marbl_share_mod, only : grazing, grazer_prey_cnt

  use marbl_interface_types, only : carbonate_type
  use marbl_interface_types, only : zooplankton_secondary_species_type
  use marbl_interface_types, only : autotroph_secondary_species_type
  use marbl_interface_types, only : photosynthetically_available_radiation_type
  use marbl_interface_types, only : dissolved_organic_matter_type

  use ecosys_restore_mod, only : ecosys_restore_type
#ifdef CCSMCOUPLED
  use POP_MCT_vars_mod, only : pop_mct_dom_o
  use POP_MCT_vars_mod, only : POP_MCT_gsMap_o
  use POP_MCT_vars_mod, only : POP_MCT_OCNID
#endif
  use strdata_interface_mod, only : strdata_input_type

  ! !INPUT PARAMETERS:
  !-----------------------------------------------------------------------
  !  include ecosystem parameters
  !  all variables from this modules have a parm_ prefix
  !-----------------------------------------------------------------------

  implicit none
  save
  private

  !-----------------------------------------------------------------------
  !  public/private member procedure declarations
  !-----------------------------------------------------------------------

  public :: &
       ecosys_init,                  &
       ecosys_tracer_ref_val,        &
       ecosys_set_sflux,             &
       ecosys_tavg_forcing,          &
       ecosys_set_interior,          &
       ecosys_write_restart

  private :: &
       ecosys_init_non_autotroph_tracer_metadata, &
       ecosys_init_forcing_monthly_every_ts, &
       ecosys_init_tavg, &
       init_particulate_terms, &
       update_particulate_terms_from_prior_level, &
       update_sinking_particle_from_prior_level, &
       compute_particulate_terms, &
       ecosys_init_sflux, &
       ecosys_init_interior_restore, &
       check_ecosys_tracer_count_consistency, &
       initialize_zooplankton_tracer_metadata, &
       initialize_autotroph_tracer_metadata, &
       setup_local_tracers, &
       setup_local_zooplankton, &
       setup_local_autotrophs, &
       autotroph_consistency_check, &
       compute_autotroph_elemental_ratios, &
       compute_photosynthetically_available_radiation, &
       compute_carbonate_chemistry, &
       compute_function_scaling, &
       compute_Pprime, &
       compute_Zprime

  type, private :: zooplankton_local_type
     real (r8) :: C  ! local copy of model zooplankton C
  end type zooplankton_local_type

  type, private :: autotroph_local_type
     real (r8) :: Chl   ! local copy of model autotroph Chl
     real (r8) :: C     ! local copy of model autotroph C
     real (r8) :: Fe    ! local copy of model autotroph Fe
     real (r8) :: Si    ! local copy of model autotroph Si
     real (r8) :: CaCO3 ! local copy of model autotroph CaCO3
  end type autotroph_local_type

  !-----------------------------------------------------------------------
  !  flags controlling which portion of code are executed
  !  usefull for debugging
  !-----------------------------------------------------------------------

  logical (log_kind) ::  lsource_sink
  logical (log_kind) ::  lflux_gas_o2
  logical (log_kind) ::  lflux_gas_co2
  logical (log_kind) ::  locmip_k1_k2_bug_fix

  !-----------------------------------------------------------------------
  !  non-autotroph relative tracer indices
  !  autotroph relative tracer indices are in autotroph derived type and are determined at run time
  !-----------------------------------------------------------------------

  integer (int_kind), parameter :: &
       po4_ind         =  1,  & ! dissolved inorganic phosphate
       no3_ind         =  2,  & ! dissolved inorganic nitrate
       sio3_ind        =  3,  & ! dissolved inorganic silicate
       nh4_ind         =  4,  & ! dissolved ammonia
       fe_ind          =  5,  & ! dissolved inorganic iron
       o2_ind          =  6,  & ! dissolved oxygen
       dic_ind         =  7,  & ! dissolved inorganic carbon
       dic_alt_co2_ind =  8,  & ! dissolved inorganic carbon with alternative CO2
       alk_ind         =  9,  & ! alkalinity
       doc_ind         = 10,  & ! dissolved organic carbon
       don_ind         = 11,  & ! dissolved organic nitrogen
       dofe_ind        = 12,  & ! dissolved organic iron
       dop_ind         = 13,  & ! dissolved organic phosphorus
       dopr_ind        = 14,  & ! refractory DOP
       donr_ind        = 15     ! refractory DON

  !-----------------------------------------------------------------------
  !  derived type & parameter for tracer index lookup
  !-----------------------------------------------------------------------

  type(ind_name_pair), dimension(ecosys_tracer_cnt) :: &
       ind_name_table

  !-----------------------------------------------------------------------
  !  options for forcing of gas fluxes
  !-----------------------------------------------------------------------

  integer (int_kind), parameter :: &
       gas_flux_forcing_iopt_drv   = 1,   &
       gas_flux_forcing_iopt_file  = 2,   &
       atm_co2_iopt_const          = 1,   &
       atm_co2_iopt_drv_prog       = 2,   &
       atm_co2_iopt_drv_diag       = 3

  integer (int_kind) :: &
       gas_flux_forcing_iopt,             &
       atm_co2_iopt,                      &
       atm_alt_co2_iopt

  real (r8)       :: &
       atm_co2_const,                     &  ! value of atmospheric co2 (ppm, dry-air, 1 atm)
       atm_alt_co2_const                     ! value of atmospheric alternative co2 (ppm, dry-air, 1 atm)

  character(char_len) :: &
       gas_flux_forcing_file    ! file containing gas flux forcing fields

  !-----------------------------------------------------------------------

  type(tracer_read) :: &
       gas_flux_fice,       & ! ice fraction for gas fluxes
       gas_flux_ws,         & ! wind speed for gas fluxes
       gas_flux_ap,         & ! atmospheric pressure for gas fluxes
       fesedflux_input        ! namelist input for iron_flux

  !-----------------------------------------------------------------------
  !  module variables related to ph computations
  !-----------------------------------------------------------------------

  real (r8), dimension(:, :, :)   , allocatable, target :: PH_PREV            ! computed ph from previous time step
  real (r8), dimension(:, :, :)   , allocatable, target :: PH_PREV_ALT_CO2    ! computed ph from previous time step, alternative CO2
  real (r8), dimension(:, :, :)   , allocatable, target :: IRON_PATCH_FLUX    ! localized iron patch flux

  !-----------------------------------------------------------------------
  !  restoring climatologies for nutrients
  !-----------------------------------------------------------------------

  real (r8), dimension(:, :, :, :), allocatable, target :: &
       FESEDFLUX      !  sedimentary Fe inputs

  character(char_len) :: &
       nutr_rest_file               ! file containing nutrient fields

  !maltrud variable restoring
  logical (log_kind) :: &
       lnutr_variable_restore       ! geographically varying nutrient restoring

  character(char_len) :: &
       nutr_variable_rest_file,   & ! file containing variable restoring info
       nutr_variable_rest_file_fmt  ! format of file containing variable restoring info

  real (r8), dimension(:, :, :), allocatable, target :: &
       NUTR_RESTORE_RTAU            ! inverse restoring timescale for variable
  ! interior restoring

  integer (int_kind), dimension(:, :, :), allocatable :: &
       NUTR_RESTORE_MAX_LEVEL       ! maximum level for applying variable
  ! interior restoring

  real (r8), dimension(:, :, :, :), allocatable :: &
       INTERP_WORK                  ! temp array for interpolate_forcing output

  type(forcing_monthly_every_ts) :: &
       dust_flux,                 & ! surface dust flux
       iron_flux,                 & ! iron component of surface dust flux
       fice_file,                 & ! ice fraction, if read from file
       xkw_file,                  & ! a * wind-speed ** 2, if read from file
       ap_file                      ! atmoshperic pressure, if read from file

  character(char_len) :: &
       ndep_data_type               ! type of ndep forcing

  type(forcing_monthly_every_ts) :: &
       nox_flux_monthly,          & ! surface NOx species flux, added to nitrate pool
       nhy_flux_monthly             ! surface NHy species flux, added to ammonium pool

  integer (int_kind) :: &
       ndep_shr_stream_year_first, & ! first year in stream to use
       ndep_shr_stream_year_last,  & ! last year in stream to use
       ndep_shr_stream_year_align    ! align ndep_shr_stream_year_first with this model year

  integer (int_kind), parameter :: &
       ndep_shr_stream_var_cnt = 2, & ! number of variables in ndep shr_stream
       ndep_shr_stream_no_ind  = 1, & ! index for NO forcing
       ndep_shr_stream_nh_ind  = 2    ! index for NH forcing

  character(char_len) :: &
       ndep_shr_stream_file          ! file containing domain and input data

  real (r8) :: &
       ndep_shr_stream_scale_factor  ! unit conversion factor

  type(strdata_input_type) :: ndep_inputlist

  type(forcing_monthly_every_ts) :: &
       din_riv_flux,              & ! river DIN species flux, added to nitrate pool
       dip_riv_flux,              & ! river DIP species flux, added to phosphate pool
       don_riv_flux,              & ! river DON flux, added to semi-lab don pool
       dop_riv_flux,              & ! river DOP flux, added to semi-lab dop pool
       dsi_riv_flux,              & ! river DSI flux, added to dsi pool
       dfe_riv_flux,              & ! river dfe flux, added to dfe pool
       dic_riv_flux,              & ! river dic flux, added to dic pool
       alk_riv_flux,              & ! river alk flux, added to alk pool
       doc_riv_flux                 ! river doc flux, added to semi-labile DOC

  !-----------------------------------------------------------------------
  !  buffer indices (into ECO_SFLUX_TAVG) for 2d fields related to surface fluxes
  !  duplicates, which are used for placing fields into multiple tavg streams,
  !  do not need separate buffer indices
  !  fields that are recoverable from the STF field do not need separate buffer indices
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       buf_ind_ECOSYS_IFRAC,          &! ice fraction
       buf_ind_ECOSYS_XKW,            &! xkw
       buf_ind_ECOSYS_ATM_PRESS,      &! atmospheric pressure
       buf_ind_PV_O2,                 &! o2 piston velocity
       buf_ind_SCHMIDT_O2,            &! O2 schmidt number
       buf_ind_O2SAT,                 &! O2 saturation
       buf_ind_CO2STAR,               &! co2star
       buf_ind_DCO2STAR,              &! dco2star
       buf_ind_pCO2SURF,              &! surface pco2
       buf_ind_DpCO2,                 &! delta pco2
       buf_ind_PV_CO2,                &! co2 piston velocity
       buf_ind_SCHMIDT_CO2,           &! co2 schmidt number
       buf_ind_DIC_GAS_FLUX,          &! dic flux
       buf_ind_PH,                    &! surface pH
       buf_ind_ATM_CO2,               &! atmospheric CO2
       buf_ind_CO2STAR_ALT_CO2,       &! co2star alternative CO2
       buf_ind_DCO2STAR_ALT_CO2,      &! dco2star alternative CO2
       buf_ind_pCO2SURF_ALT_CO2,      &! surface pco2 alternative CO2
       buf_ind_DpCO2_ALT_CO2,         &! delta pco2 alternative CO2
       buf_ind_DIC_GAS_FLUX_ALT_CO2,  &! dic flux alternative CO2
       buf_ind_PH_ALT_CO2,            &! surface pH alternative CO2
       buf_ind_ATM_ALT_CO2,           &! atmospheric alternative CO2
       buf_ind_IRON_FLUX,             &! iron flux
       buf_ind_NOx_FLUX,              &! nox flux
       buf_ind_DIN_RIV_FLUX,          &! din river flux
       buf_ind_DFE_RIV_FLUX,          &! dfe river flux
       buf_ind_DIC_RIV_FLUX,          &! dic river flux
       buf_ind_ALK_RIV_FLUX            ! alk river flux

  !-----------------------------------------------------------------------
  !  define array for holding flux-related quantities that need to be time-averaged
  !  this is necessary since the forcing routines are called before tavg flags
  !-----------------------------------------------------------------------

  real (r8), dimension(:, :, :, :), allocatable :: ECO_SFLUX_TAVG

  !-----------------------------------------------------------------------
  !  average surface tracer value related variables
  !  used as reference value for virtual flux computations
  !-----------------------------------------------------------------------

  logical (log_kind) , dimension(ecosys_tracer_cnt) :: vflux_flag ! which tracers get virtual fluxes applied
  real (r8)          , dimension(ecosys_tracer_cnt) :: surf_avg   ! average surface tracer values
  integer (int_kind)          :: comp_surf_avg_flag               ! time flag id for computing average surface tracer values
  logical (log_kind) , public :: ecosys_qsw_distrb_const

  !-----------------------------------------------------------------------
  !  iron patch fertilization
  !-----------------------------------------------------------------------

  logical (log_kind)  :: liron_patch               ! flag for iron patch fertilization
  character(char_len) :: iron_patch_flux_filename  ! file containing name of iron patch file
  integer (int_kind)  :: iron_patch_month          !  integer month to add patch flux

  !-----------------------------------------------------------------------
  !  timers
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       ecosys_shr_strdata_advance_timer,        &
       ecosys_comp_CO3terms_timer,              &
       ecosys_sflux_timer

  !-----------------------------------------------------------------------
  !  named field indices
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
       totChl_surf_nf_ind = 0,    & ! total chlorophyll in surface layer
       sflux_co2_nf_ind   = 0,    & ! air-sea co2 gas flux
       atm_co2_nf_ind     = 0       ! atmospheric co2

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------

  real (r8), parameter :: &
       phlo_surf_init = 7.0_r8, & ! low bound for surface ph for no prev soln
       phhi_surf_init = 9.0_r8, & ! high bound for surface ph for no prev soln
       phlo_3d_init = 6.0_r8,   & ! low bound for subsurface ph for no prev soln
       phhi_3d_init = 9.0_r8,   & ! high bound for subsurface ph for no prev soln
       del_ph = 0.20_r8           ! delta-ph for prev soln


  !-----------------------------------------------------------------------

  !EOC
  !*****************************************************************************


contains

  !*****************************************************************************
  !BOP
  ! !IROUTINE: ecosys_init
  ! !INTERFACE:

  subroutine ecosys_init(nl_buffer, init_ts_file_fmt, read_restart_filename, &
       tracer_d_module, TRACER_MODULE, &
       lmarginal_seas, ecosys_restore, &
       saved_state, &
       errorCode, marbl_status)

    ! !DESCRIPTION:
    !  Initialize ecosys tracer module. This involves setting metadata, reading
    !  the module namelist, setting initial conditions, setting up forcing,
    !  and defining additional tavg variables.
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_interface_constants, only : marbl_nl_buffer_size, marbl_status_ok, marbl_status_could_not_read_namelist
    use marbl_interface_types, only : marbl_status_type
    use marbl_interface_types, only : marbl_saved_state_type

    implicit none

    ! !INPUT PARAMETERS:

    character(marbl_nl_buffer_size), intent(in) :: nl_buffer

    character (*), intent(in) :: &
         init_ts_file_fmt,     & ! format (bin or nc) for input file
         read_restart_filename   ! file name for restart file

    logical (kind=log_kind), intent(in) :: &
         lmarginal_seas               ! Is ecosystem active in marginal seas ?

    ! !INPUT/OUTPUT PARAMETERS:

    type (tracer_field), dimension(:), intent(inout) :: &
         tracer_d_module   ! descriptors for each tracer

    real (r8), dimension(:, :, :, :, :, :), &
         intent(inout) :: TRACER_MODULE

    type(ecosys_restore_type), intent(inout) :: ecosys_restore

    type(marbl_saved_state_type), intent(inout) :: saved_state

    ! !OUTPUT PARAMETERS:

    integer (POP_i4), intent(out) :: &
         errorCode

    type(marbl_status_type), intent(out) :: marbl_status

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_init'

    character(char_len) :: &
         init_ecosys_option,        & ! option for initialization of bgc
         init_ecosys_init_file,     & ! filename for option 'file'
         init_ecosys_init_file_fmt, & ! file format for option 'file'
         comp_surf_avg_freq_opt,    & ! choice for freq of comp_surf_avg
         gas_flux_forcing_opt,      & ! option for forcing gas fluxes
         atm_co2_opt,               & ! option for atmospheric co2 concentration
         atm_alt_co2_opt              ! option for atmospheric alternative CO2

    type(tracer_read), dimension(ecosys_tracer_cnt) :: &
         tracer_init_ext              ! namelist variable for initializing tracers

    type(tracer_read) :: &
         dust_flux_input,           & ! namelist input for dust_flux
         iron_flux_input,           & ! namelist input for iron_flux
         nox_flux_monthly_input,    & ! namelist input for nox_flux_monthly
         nhy_flux_monthly_input,    & ! namelist input for nhy_flux_monthly
         din_riv_flux_input,        & ! namelist input for din_riv_flux
         dip_riv_flux_input,        & ! namelist input for dip_riv_flux
         don_riv_flux_input,        & ! namelist input for don_riv_flux
         dop_riv_flux_input,        & ! namelist input for dop_riv_flux
         dsi_riv_flux_input,        & ! namelist input for dsi_riv_flux
         dfe_riv_flux_input,        & ! namelist input for dfe_riv_flux
         dic_riv_flux_input,        & ! namelist input for dic_riv_flux
         alk_riv_flux_input,        & ! namelist input for alk_riv_flux
         doc_riv_flux_input           ! namelist input for doc_riv_flux

    integer (int_kind) :: &
         non_living_biomass_ecosys_tracer_cnt, & ! number of non-autotroph ecosystem tracers
         auto_ind,                  & ! autotroph functional group index
         n,                         & ! index for looping over tracers
         k,                         & ! index for looping over depth levels
         iblock,                    & ! index for looping over blocks
         nml_error                    ! namelist i/o error flag

    integer (int_kind) :: &
         zoo_ind                    ! zooplankton functional group index


    integer (int_kind) :: &
         comp_surf_avg_freq_iopt,   & ! choice for freq of comp_surf_avg
         comp_surf_avg_freq           ! choice for freq of comp_surf_avg

    logical (log_kind) :: &
         use_nml_surf_vals            ! do namelist surf values override values from restart file

    logical (log_kind) :: &
         lecovars_full_depth_tavg     ! should ecosystem vars be written full depth

    !-----------------------------------------------------------------------
    !  values to be used when comp_surf_avg_freq_opt==never
    !-----------------------------------------------------------------------

    real (r8) :: &
         surf_avg_dic_const, surf_avg_alk_const

    namelist /ecosys_nml/ &
         init_ecosys_option, init_ecosys_init_file, tracer_init_ext, &
         init_ecosys_init_file_fmt, &
         dust_flux_input, iron_flux_input, fesedflux_input, &
         ndep_data_type, nox_flux_monthly_input, nhy_flux_monthly_input, &
         ndep_shr_stream_year_first, ndep_shr_stream_year_last, &
         ndep_shr_stream_year_align, ndep_shr_stream_file, &
         ndep_shr_stream_scale_factor, &
         din_riv_flux_input, dip_riv_flux_input, don_riv_flux_input, &
         dop_riv_flux_input, dsi_riv_flux_input, dfe_riv_flux_input, &
         dic_riv_flux_input, alk_riv_flux_input, doc_riv_flux_input, &
         gas_flux_forcing_opt, gas_flux_forcing_file, &
         gas_flux_fice, gas_flux_ws, gas_flux_ap, &
         nutr_rest_file, &
         comp_surf_avg_freq_opt, comp_surf_avg_freq,  &
         use_nml_surf_vals, surf_avg_dic_const, surf_avg_alk_const, &
         ecosys_qsw_distrb_const, &
         lsource_sink, lflux_gas_o2, lflux_gas_co2, locmip_k1_k2_bug_fix, &
         lnutr_variable_restore, nutr_variable_rest_file,  &
         nutr_variable_rest_file_fmt, atm_co2_opt, atm_co2_const, &
         atm_alt_co2_opt, atm_alt_co2_const, &
         liron_patch, iron_patch_flux_filename, iron_patch_month, &
         lecovars_full_depth_tavg

    character (char_len) :: &
         ecosys_restart_filename  ! modified file name for restart file

    real (r8), dimension (nx_block, ny_block) :: WORK
    marbl_status%status = marbl_status_ok
    marbl_status%message = ''
    errorCode = POP_Success

    call ecosys_init_forcing_monthly_every_ts()

    call marbl_params_init(nl_buffer, marbl_status)
    if (marbl_status%status /= marbl_status_ok) then
       return
    end if
    ! FIXME(bja, 2015-01) need to have a flag if params should be printed!
    ! if (stdout > 0) then
    !    call marbl_params_print(stdout)
    ! end if

    call ecosys_init_non_autotroph_tracer_metadata(tracer_d_module, non_living_biomass_ecosys_tracer_cnt)

    call check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt)

    call initialize_zooplankton_tracer_metadata(tracer_d_module, &
         non_living_biomass_ecosys_tracer_cnt, n)

    call initialize_autotroph_tracer_metadata(tracer_d_module, n)
    !-----------------------------------------------------------------------
    !  initialize ind_name_table
    !-----------------------------------------------------------------------

    do n = 1, ecosys_tracer_cnt
       ind_name_table(n) = ind_name_pair(n, tracer_d_module(n)%short_name)
    end do

    !-----------------------------------------------------------------------
    !  default namelist settings
    !-----------------------------------------------------------------------

    init_ecosys_option = 'unknown'
    init_ecosys_init_file = 'unknown'
    init_ecosys_init_file_fmt = 'bin'

    gas_flux_forcing_opt  = 'drv'
    gas_flux_forcing_file = 'unknown'

    gas_flux_fice%filename     = 'unknown'
    gas_flux_fice%file_varname = 'FICE'
    gas_flux_fice%scale_factor = c1
    gas_flux_fice%default_val  = c0
    gas_flux_fice%file_fmt     = 'bin'

    gas_flux_ws%filename     = 'unknown'
    gas_flux_ws%file_varname = 'XKW'
    gas_flux_ws%scale_factor = c1
    gas_flux_ws%default_val  = c0
    gas_flux_ws%file_fmt     = 'bin'

    gas_flux_ap%filename     = 'unknown'
    gas_flux_ap%file_varname = 'P'
    gas_flux_ap%scale_factor = c1
    gas_flux_ap%default_val  = c0
    gas_flux_ap%file_fmt     = 'bin'

    nutr_rest_file = 'unknown'

    !maltrud variable restoring
    lnutr_variable_restore      = .false.
    nutr_variable_rest_file     = 'unknown'
    nutr_variable_rest_file_fmt = 'bin'

    dust_flux_input%filename     = 'unknown'
    dust_flux_input%file_varname = 'dust_flux'
    dust_flux_input%scale_factor = c1
    dust_flux_input%default_val  = c0
    dust_flux_input%file_fmt     = 'bin'

    iron_flux_input%filename     = 'unknown'
    iron_flux_input%file_varname = 'iron_flux'
    iron_flux_input%scale_factor = c1
    iron_flux_input%default_val  = c0
    iron_flux_input%file_fmt     = 'bin'

    fesedflux_input%filename     = 'unknown'
    fesedflux_input%file_varname = 'FESEDFLUXIN'
    fesedflux_input%scale_factor = c1
    fesedflux_input%default_val  = c0
    fesedflux_input%file_fmt     = 'bin'

    ndep_data_type              = 'monthly-calendar'

    nox_flux_monthly_input%filename     = 'unknown'
    nox_flux_monthly_input%file_varname = 'nox_flux'
    nox_flux_monthly_input%scale_factor = c1
    nox_flux_monthly_input%default_val  = c0
    nox_flux_monthly_input%file_fmt     = 'bin'

    nhy_flux_monthly_input%filename     = 'unknown'
    nhy_flux_monthly_input%file_varname = 'nhy_flux'
    nhy_flux_monthly_input%scale_factor = c1
    nhy_flux_monthly_input%default_val  = c0
    nhy_flux_monthly_input%file_fmt     = 'bin'

    ndep_shr_stream_year_first = 1
    ndep_shr_stream_year_last  = 1
    ndep_shr_stream_year_align = 1
    ndep_shr_stream_file       = 'unknown'
    ndep_shr_stream_scale_factor = c1

    din_riv_flux_input%filename     = 'unknown'
    din_riv_flux_input%file_varname = 'din_riv_flux'
    din_riv_flux_input%scale_factor = c1
    din_riv_flux_input%default_val  = c0
    din_riv_flux_input%file_fmt     = 'nc'

    dip_riv_flux_input%filename     = 'unknown'
    dip_riv_flux_input%file_varname = 'dip_riv_flux'
    dip_riv_flux_input%scale_factor = c1
    dip_riv_flux_input%default_val  = c0
    dip_riv_flux_input%file_fmt     = 'nc'

    don_riv_flux_input%filename     = 'unknown'
    don_riv_flux_input%file_varname = 'don_riv_flux'
    don_riv_flux_input%scale_factor = c1
    don_riv_flux_input%default_val  = c0
    don_riv_flux_input%file_fmt     = 'nc'

    dop_riv_flux_input%filename     = 'unknown'
    dop_riv_flux_input%file_varname = 'dop_riv_flux'
    dop_riv_flux_input%scale_factor = c1
    dop_riv_flux_input%default_val  = c0
    dop_riv_flux_input%file_fmt     = 'nc'

    dsi_riv_flux_input%filename     = 'unknown'
    dsi_riv_flux_input%file_varname = 'dsi_riv_flux'
    dsi_riv_flux_input%scale_factor = c1
    dsi_riv_flux_input%default_val  = c0
    dsi_riv_flux_input%file_fmt     = 'nc'

    dfe_riv_flux_input%filename     = 'unknown'
    dfe_riv_flux_input%file_varname = 'dfe_riv_flux'
    dfe_riv_flux_input%scale_factor = c1
    dfe_riv_flux_input%default_val  = c0
    dfe_riv_flux_input%file_fmt     = 'nc'

    dic_riv_flux_input%filename     = 'unknown'
    dic_riv_flux_input%file_varname = 'dic_riv_flux'
    dic_riv_flux_input%scale_factor = c1
    dic_riv_flux_input%default_val  = c0
    dic_riv_flux_input%file_fmt     = 'nc'

    alk_riv_flux_input%filename     = 'unknown'
    alk_riv_flux_input%file_varname = 'alk_riv_flux'
    alk_riv_flux_input%scale_factor = c1
    alk_riv_flux_input%default_val  = c0
    alk_riv_flux_input%file_fmt     = 'nc'

    doc_riv_flux_input%filename     = 'unknown'
    doc_riv_flux_input%file_varname = 'doc_riv_flux'
    doc_riv_flux_input%scale_factor = c1
    doc_riv_flux_input%default_val  = c0
    doc_riv_flux_input%file_fmt     = 'nc'

    do n = 1, ecosys_tracer_cnt
       tracer_init_ext(n)%mod_varname  = 'unknown'
       tracer_init_ext(n)%filename     = 'unknown'
       tracer_init_ext(n)%file_varname = 'unknown'
       tracer_init_ext(n)%scale_factor = c1
       tracer_init_ext(n)%default_val  = c0
       tracer_init_ext(n)%file_fmt     = 'bin'
    end do

    lsource_sink          = .true.
    lflux_gas_o2          = .true.
    lflux_gas_co2         = .true.
    locmip_k1_k2_bug_fix  = .true.

    comp_surf_avg_freq_opt        = 'never'
    comp_surf_avg_freq            = 1
    use_nml_surf_vals             = .false.
    surf_avg_dic_const            = 1944.0_r8
    surf_avg_alk_const            = 2225.0_r8

    ecosys_qsw_distrb_const  = .true.

    liron_patch              = .false.
    iron_patch_flux_filename = 'unknown_iron_patch_filename'
    iron_patch_month         = 1

    atm_co2_opt   = 'const'
    atm_co2_const = 280.0_r8

    atm_alt_co2_opt   = 'const'
    atm_alt_co2_const = 280.0_r8

    lecovars_full_depth_tavg = .false.

    ! read the namelist buffer on every processor
    read(nl_buffer, nml=ecosys_nml, iostat=nml_error)
    if (nml_error /= 0) then
       marbl_status%status = marbl_status_could_not_read_namelist
       marbl_status%message = "ERROR: "//subname//"(): could not read ecosys_nml."
       return
    end if

    if (my_task == master_task) then
       write(stdout, blank_fmt)
       write(stdout, ndelim_fmt)
       write(stdout, blank_fmt)
       write(stdout, *) ' ecosys:'
       write(stdout, blank_fmt)
       write(stdout, *) ' ecosys_nml namelist settings:'
       write(stdout, blank_fmt)
       write(stdout, ecosys_nml)
       write(stdout, blank_fmt)
       write(stdout, delim_fmt)
    endif

    !-----------------------------------------------------------------------
    ! reassign values temporary input values to correct arrays
    !-----------------------------------------------------------------------

    if (trim(gas_flux_forcing_opt) == 'drv') then
       gas_flux_forcing_iopt = gas_flux_forcing_iopt_drv
    else if (trim(gas_flux_forcing_opt) == 'file') then
       gas_flux_forcing_iopt = gas_flux_forcing_iopt_file
    else
       call document(subname, 'gas_flux_forcing_opt', gas_flux_forcing_opt)
       call exit_POP(sigAbort, 'unknown gas_flux_forcing_opt')
    endif

    fice_file%input = gas_flux_fice

    xkw_file%input = gas_flux_ws

    ap_file%input = gas_flux_ap

    dust_flux%input = dust_flux_input

    iron_flux%input = iron_flux_input

    nox_flux_monthly%input = nox_flux_monthly_input

    nhy_flux_monthly%input = nhy_flux_monthly_input

    din_riv_flux%input = din_riv_flux_input

    dip_riv_flux%input = dip_riv_flux_input

    don_riv_flux%input = don_riv_flux_input

    dop_riv_flux%input = dop_riv_flux_input

    dsi_riv_flux%input = dsi_riv_flux_input

    dfe_riv_flux%input = dfe_riv_flux_input

    dic_riv_flux%input = dic_riv_flux_input

    alk_riv_flux%input = alk_riv_flux_input

    doc_riv_flux%input = doc_riv_flux_input


    !-----------------------------------------------------------------------
    !  initialize modules with dependancies on above variables
    !-----------------------------------------------------------------------

    call ecosys_restore%Init(nml_filename, nml_in, &
         ind_name_table)

    !-----------------------------------------------------------------------
    !  set variables immediately dependent on namelist variables
    !-----------------------------------------------------------------------

    select case (comp_surf_avg_freq_opt)
    case ('never')
       comp_surf_avg_freq_iopt = freq_opt_never
    case ('nyear')
       comp_surf_avg_freq_iopt = freq_opt_nyear
    case ('nmonth')
       comp_surf_avg_freq_iopt = freq_opt_nmonth
    case default
       call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
       call exit_POP(sigAbort, 'unknown comp_surf_avg_freq_opt')
    end select

    call init_time_flag('ecosys_comp_surf_avg', comp_surf_avg_flag, &
         default=.false., freq_opt=comp_surf_avg_freq_iopt,  &
         freq=comp_surf_avg_freq, owner='ecosys_init')

    select case (atm_co2_opt)
    case ('const')
       atm_co2_iopt = atm_co2_iopt_const
    case ('drv_prog')
       atm_co2_iopt = atm_co2_iopt_drv_prog
    case ('drv_diag')
       atm_co2_iopt = atm_co2_iopt_drv_diag
    case default
       call document(subname, 'atm_co2_opt', atm_co2_opt)
       call exit_POP(sigAbort, 'unknown atm_co2_opt')
    end select

    select case (atm_alt_co2_opt)
    case ('const')
       atm_alt_co2_iopt = atm_co2_iopt_const
    case default
       call document(subname, 'atm_alt_co2_opt', atm_alt_co2_opt)
       call exit_POP(sigAbort, 'unknown atm_alt_co2_opt')
    end select

    !-----------------------------------------------------------------------
    !  namelist consistency checking
    !-----------------------------------------------------------------------

    if (use_nml_surf_vals .and. comp_surf_avg_freq_iopt /= freq_opt_never) then
       call document(subname, 'use_nml_surf_vals', use_nml_surf_vals)
       call document(subname, 'comp_surf_avg_freq_opt', comp_surf_avg_freq_opt)
       call exit_POP(sigAbort, 'use_nml_surf_vals can only be .true. if ' /&
            &/ ' comp_surf_avg_freq_opt is never')
    endif

    !-----------------------------------------------------------------------
    !  initialize virtual flux flag array
    !-----------------------------------------------------------------------

    vflux_flag = .false.
    vflux_flag(dic_ind) = .true.
    vflux_flag(alk_ind) = .true.
    vflux_flag(dic_alt_co2_ind) = .true.

    !-----------------------------------------------------------------------
    !  allocate various ecosys allocatable module variables
    !-----------------------------------------------------------------------

    allocate( PH_PREV(nx_block, ny_block, max_blocks_clinic) )
    allocate( PH_PREV_ALT_CO2(nx_block, ny_block, max_blocks_clinic) )

    !-----------------------------------------------------------------------
    !  allocate and initialize LAND_MASK
    !-----------------------------------------------------------------------

    if (lmarginal_seas) then
       saved_state%land_mask = REGION_MASK /= 0
    else
       saved_state%land_mask = REGION_MASK > 0
    endif

    !-----------------------------------------------------------------------
    !  initialize tracers
    !-----------------------------------------------------------------------

    select case (init_ecosys_option)

    case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid')

       ecosys_restart_filename = char_blank

       if (init_ecosys_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call document(subname, 'no restart file to read ecosys from')
             call exit_POP(sigAbort, 'stopping in ' /&
                  &/ subname)
          endif
          ecosys_restart_filename = read_restart_filename
          init_ecosys_init_file_fmt = init_ts_file_fmt

       else  ! do not read from TS restart file

          ecosys_restart_filename = trim(init_ecosys_init_file)

       endif

       call rest_read_tracer_block(init_ecosys_init_file_fmt, &
            ecosys_restart_filename,   &
            tracer_d_module,           &
            TRACER_MODULE)

       if (field_exists_in_file(init_ecosys_init_file_fmt, &
            ecosys_restart_filename, &
            'PH_SURF')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_SURF', PH_PREV)
       else
          call document(subname, 'PH_SURF does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV to 0')
          PH_PREV = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, &
            ecosys_restart_filename, &
            'PH_SURF_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_SURF_ALT_CO2', PH_PREV_ALT_CO2)
       else
          call document(subname, 'PH_SURF_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2 to 0')
          PH_PREV_ALT_CO2 = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, &
            ecosys_restart_filename, &
            'PH_3D')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D', saved_state%PH_PREV_3D)
       else
          call document(subname, 'PH_3D does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_3D to 0')
          saved_state%PH_PREV_3D  = c0
       endif

       if (field_exists_in_file(init_ecosys_init_file_fmt, &
            ecosys_restart_filename, &
            'PH_3D_ALT_CO2')) then
          call read_field(init_ecosys_init_file_fmt, &
               ecosys_restart_filename,   &
               'PH_3D_ALT_CO2', saved_state%PH_PREV_ALT_CO2_3D)
       else
          call document(subname, 'PH_3D_ALT_CO2 does not exist in ' /&
               &/ trim(ecosys_restart_filename) /&
               &/ ', setting PH_PREV_ALT_CO2_3D to 0')
          saved_state%PH_PREV_ALT_CO2_3D = c0
       endif

       if (use_nml_surf_vals) then
          surf_avg = c0
          surf_avg(dic_ind) = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind) = surf_avg_alk_const
       else
          call extract_surf_avg(init_ecosys_init_file_fmt,     &
               ecosys_restart_filename,       &
               ecosys_tracer_cnt, vflux_flag, &
               ind_name_table, surf_avg)
       endif

       call eval_time_flag(comp_surf_avg_flag) ! evaluates time_flag(comp_surf_avg_flag)%value via time_to_do

       if (check_time_flag(comp_surf_avg_flag)) &
            call comp_surf_avg(TRACER_MODULE(:, :, 1, :, oldtime, :), &
            TRACER_MODULE(:, :, 1, :, curtime, :), &
            ecosys_tracer_cnt, vflux_flag, surf_avg)

    case ('file', 'ccsm_startup')
       call document(subname, 'ecosystem vars being read from separate files')

       call file_read_tracer_block(init_ecosys_init_file_fmt, &
            init_ecosys_init_file,     &
            tracer_d_module,           &
            ind_name_table,            &
            tracer_init_ext,           &
            TRACER_MODULE)

       if (n_topo_smooth > 0) then
          do n = 1, ecosys_tracer_cnt
             do k=1, km
                call fill_points(k, TRACER_MODULE(:, :, k, n, oldtime, :), &
                     errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(oldtime)')
                   return
                endif

                call fill_points(k, TRACER_MODULE(:, :, k, n, curtime, :), &
                     errorCode)

                if (errorCode /= POP_Success) then
                   call POP_ErrorSet(errorCode, &
                        'ecosys_init: error in fill points for tracers(newtime)')
                   return
                endif

             enddo
          enddo
       endif

       PH_PREV = c0
       PH_PREV_ALT_CO2 = c0
       saved_state%PH_PREV_3D = c0
       saved_state%PH_PREV_ALT_CO2_3D = c0

       if (use_nml_surf_vals) then
          surf_avg = c0
          surf_avg(dic_ind) = surf_avg_dic_const
          surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
          surf_avg(alk_ind) = surf_avg_alk_const
       else
          call comp_surf_avg(TRACER_MODULE(:, :, 1, :, oldtime, :), &
               TRACER_MODULE(:, :, 1, :, curtime, :), &
               ecosys_tracer_cnt, vflux_flag, surf_avg)
       endif

    case default
       call document(subname, 'init_ecosys_option', init_ecosys_option)
       call exit_POP(sigAbort, 'unknown init_ecosys_option')

    end select

    !-----------------------------------------------------------------------
    !  register Chl field for short-wave absorption
    !  apply land mask to tracers
    !  set Chl field for short-wave absorption
    !-----------------------------------------------------------------------

    call named_field_register('model_chlorophyll', totChl_surf_nf_ind)

    !$OMP PARALLEL DO PRIVATE(iblock, n, k, WORK)
    do iblock=1, nblocks_clinic
       do n = 1, ecosys_tracer_cnt
          do k = 1, km
             where (.not. saved_state%land_mask(:, :, iblock) .or. k > KMT(:, :, iblock))
                TRACER_MODULE(:, :, k, n, curtime, iblock) = c0
                TRACER_MODULE(:, :, k, n, oldtime, iblock) = c0
             end where
          end do
       end do

       WORK = c0
       do auto_ind = 1, autotroph_cnt
          n = autotrophs(auto_ind)%Chl_ind
          WORK = WORK + max(c0, p5*(TRACER_MODULE(:, :, 1, n, oldtime, iblock) + &
               TRACER_MODULE(:, :, 1, n, curtime, iblock)))
       end do
       call named_field_set(totChl_surf_nf_ind, iblock, WORK)
    enddo
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !  timer init
    !-----------------------------------------------------------------------

    call get_timer(ecosys_comp_CO3terms_timer, 'comp_CO3terms', &
         nblocks_clinic, distrb_clinic%nprocs)
    call get_timer(ecosys_sflux_timer, 'ECOSYS_SFLUX', 1, &
         distrb_clinic%nprocs)
    if (ndep_data_type == 'shr_stream') then
       call get_timer(ecosys_shr_strdata_advance_timer, &
            'ecosys_shr_strdata_advance', 1, distrb_clinic%nprocs)
    endif

    !-----------------------------------------------------------------------
    !  call other initialization subroutines
    !-----------------------------------------------------------------------

    call ecosys_init_tavg()
    call ecosys_init_sflux(saved_state)
    call ecosys_init_interior_restore(saved_state, ecosys_restore)

    !-----------------------------------------------------------------------
    !  set lfull_depth_tavg flag for short-lived ecosystem tracers
    !-----------------------------------------------------------------------
    do zoo_ind = 1, zooplankton_cnt
       n = zooplankton(zoo_ind)%C_ind
       tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg
    end do

    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%Chl_ind
       tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%C_ind
       tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%Fe_ind
       tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%Si_ind
       if (n > 0) then
          tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg
       endif

       n = autotrophs(auto_ind)%CaCO3_ind
       if (n > 0) then
          tracer_d_module(n)%lfull_depth_tavg = lecovars_full_depth_tavg
       endif
    end do

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_init

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_init_tavg
  ! !INTERFACE:

  subroutine ecosys_init_tavg

    ! !DESCRIPTION:
    !  call define_tavg_field for nonstandard tavg fields
    !
    ! !REVISION HISTORY:
    !  same as module
    !

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: &
         buf_len           ! how many surface flux fields are stored in ECO_SFLUX_TAVG

    !-----------------------------------------------------------------------
    !  2D fields related to surface fluxes
    !-----------------------------------------------------------------------

    buf_len = 0

    buf_len = buf_len+1
    buf_ind_ECOSYS_IFRAC = buf_len

    buf_len = buf_len+1
    buf_ind_ECOSYS_XKW = buf_len

    buf_len = buf_len+1
    buf_ind_ECOSYS_ATM_PRESS = buf_len

    buf_len = buf_len+1
    buf_ind_PV_O2 = buf_len

    buf_len = buf_len+1
    buf_ind_SCHMIDT_O2 = buf_len

    buf_len = buf_len+1
    buf_ind_O2SAT = buf_len

    buf_len = buf_len+1
    buf_ind_CO2STAR = buf_len

    buf_len = buf_len+1
    buf_ind_DCO2STAR = buf_len

    buf_len = buf_len+1
    buf_ind_pCO2SURF = buf_len

    buf_len = buf_len+1
    buf_ind_DpCO2 = buf_len

    buf_len = buf_len+1
    buf_ind_PV_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_SCHMIDT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_DIC_GAS_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_PH = buf_len

    buf_len = buf_len+1
    buf_ind_ATM_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_CO2STAR_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_DCO2STAR_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_pCO2SURF_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_DpCO2_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_DIC_GAS_FLUX_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_PH_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_ATM_ALT_CO2 = buf_len

    buf_len = buf_len+1
    buf_ind_IRON_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_NOx_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_DIN_RIV_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_DFE_RIV_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_DIC_RIV_FLUX = buf_len

    buf_len = buf_len+1
    buf_ind_ALK_RIV_FLUX = buf_len

    !-----------------------------------------------------------------------

    allocate(ECO_SFLUX_TAVG(nx_block, ny_block, buf_len, max_blocks_clinic))
    ECO_SFLUX_TAVG = c0

  end subroutine ecosys_init_tavg

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_set_interior
  ! !INTERFACE:

  subroutine ecosys_set_interior(i, c, num_columns, domain, &
       marbl_diagnostics, &
       saved_state, ecosys_restore, &
       ecosys_interior_share, ecosys_zooplankton_share, &
       ecosys_autotroph_share, ecosys_particulate_share, &
       tracer_module, &
       dtracer, lexport_shared_vars,         &
       bid)

    ! !DESCRIPTION:
    !  Compute time derivatives for ecosystem state variables
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_share_mod, only : sinking_particle
    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_share_mod, only : column_sinking_particle_to_slab_sinking_particle
    use marbl_share_mod, only : slab_sinking_particle_to_column_sinking_particle
    use marbl_share_mod, only : ecosys_interior_share_type
    use marbl_share_mod, only : ecosys_autotroph_share_type
    use marbl_share_mod, only : ecosys_zooplankton_share_type
    use marbl_share_mod, only : ecosys_particulate_share_type
    use marbl_share_mod, only : marbl_interior_share_type
    use marbl_share_mod, only : marbl_autotroph_share_type
    use marbl_share_mod, only : marbl_zooplankton_share_type
    use marbl_share_mod, only : marbl_particulate_share_type
    use marbl_share_mod, only : column_interior_share_to_slab_interior_share
    use marbl_share_mod, only : column_zooplankton_share_to_slab_zooplankton_share
    use marbl_share_mod, only : column_autotroph_share_to_slab_autotroph_share
    use marbl_share_mod, only : column_particulate_share_to_slab_particulate_share

    use marbl_interface_types, only : marbl_diagnostics_type
    use marbl_interface_types, only : marbl_column_domain_type
    use marbl_interface_types, only : marbl_saved_state_type

    use marbl_parms, only : epsC, epsTinv

    integer (int_kind), intent(in) :: i         ! index for looping over nx_block dimension
    integer (int_kind), intent(in) :: c         ! column index
    integer (int_kind), intent(in) :: num_columns
    ! FIXME(bja, 2015-08) domain will become intent(in) once the loops
    ! are moved!
    type(marbl_column_domain_type), intent(inout) :: domain

    real (r8), intent(in)  :: tracer_module(ecosys_tracer_cnt, km)       ! tracer values
    logical (log_kind), intent(in)  :: lexport_shared_vars ! flag to save shared_vars or not
    integer (int_kind), intent(in) :: bid       ! local_block id
    real (r8), intent(out) :: dtracer(ecosys_tracer_cnt, km)      ! computed source/sink terms

    type(marbl_diagnostics_type), intent(inout) :: marbl_diagnostics(km)

    type(marbl_saved_state_type), intent(inout) :: saved_state
    type(ecosys_restore_type), intent(inout) :: ecosys_restore

    type(ecosys_interior_share_type)   , intent(inout) :: ecosys_interior_share(:)
    type(ecosys_zooplankton_share_type)   , intent(inout) :: ecosys_zooplankton_share(:)
    type(ecosys_autotroph_share_type)   , intent(inout) :: ecosys_autotroph_share(:)
    type(ecosys_particulate_share_type), intent(inout) :: ecosys_particulate_share(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: &
         subname = 'ecosys_mod:ecosys_set_interior'

    real (r8) :: f_loss_thres ! fraction of grazing loss reduction at depth

    real (r8) :: nitrif(km)    ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
    real (r8) :: denitrif(km)  ! WC nitrification (NO3 -> N2) (mmol N/m^3/sec)

    real (r8) :: O2_production(km)  ! O2 production
    real (r8) :: O2_consumption(km) ! O2 consumption

    integer (int_kind) :: auto_ind  ! autotroph functional group index
    integer (int_kind) :: auto_ind2 ! autotroph functional group index
    integer (int_kind) :: zoo_ind   ! zooplankton functional group index
    integer (int_kind) :: zoo_ind2  ! zooplankton functional group index
    integer (int_kind) :: prey_ind  ! grazee group index
    integer (int_kind) :: pred_ind  ! grazer group index
    integer (int_kind) :: kk        ! index for looping over k levels
    integer (int_kind) :: d         ! diag index index
    integer (int_kind) :: n         ! tracer index
    integer (int_kind) :: k         ! vertical level index

    real (r8) :: Tfunc(km)
    real (r8) :: Fe_scavenge_rate(km) ! annual scavenging rate of iron as % of ambient
    real (r8) :: Fe_scavenge(km)      ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)

    real(r8) :: tracer_local(ecosys_tracer_cnt, km)
    ! FIXME(bja, 2014-10) size should be
    ! (non_living_biomass_ecosys_tracer_cnt, km) but
    ! non-living-biomass-tracer-count isn't global and I'm reluctant
    ! to make it right now.
    real(r8) :: restore_local(ecosys_tracer_cnt, km) ! local restoring terms for nutrients (mmol ./m^3/sec)

    type(zooplankton_local_type) :: zooplankton_local(zooplankton_cnt, km)
    type(autotroph_local_type) :: autotroph_local(autotroph_cnt, km)

    real(r8) :: QA_dust_def(km)
    real(r8) :: PAR_out ! photosynthetically available radiation (W/m^2)
    real(r8) :: dust_flux_in

    type(autotroph_secondary_species_type) :: autotroph_secondary_species(autotroph_cnt, km)
    type(zooplankton_secondary_species_type) :: zooplankton_secondary_species(zooplankton_cnt, km)

    type(photosynthetically_available_radiation_type) :: PAR(km)

    type(dissolved_organic_matter_type) :: dissolved_organic_matter(km)

    type(carbonate_type) :: carbonate(km)

    type(marbl_interior_share_type) :: marbl_interior_share(km)
    type(marbl_zooplankton_share_type) :: marbl_zooplankton_share(zooplankton_cnt, km)
    type(marbl_autotroph_share_type) :: marbl_autotroph_share(autotroph_cnt, km)
    type(marbl_particulate_share_type) :: marbl_particulate_share

    real(r8) :: ph_prev_3d(km)
    real(r8) :: ph_prev_alt_co2_3d(km)

    real(r8) :: zsat_calcite(km) ! Calcite Saturation Depth
    real(r8) :: zsat_aragonite(km) ! Aragonite Saturation Depth

    real(r8) :: co3_calc_anom(km) ! CO3 concentration above calcite saturation at k-1
    real(r8) :: co3_arag_anom(km) ! CO3 concentration above aragonite saturation at k-1

    real(r8) :: sed_denitrif(km) ! sedimentary denitrification (nmol N/cm^3/sec)
    real(r8) :: other_remin(km)  ! organic C remin not due oxic or denitrif (nmolC/cm^3/sec)
    ! NOTE(bja, 2015-07) vectorization: arrays that are (n, k, c, i)
    ! probably can not be vectorized reasonably over c without memory
    ! copies. If we break up the main k loop, some of the (k, c) loops
    ! can probably be vectorized over k and / or c!

    !-----------------------------------------------------------------------

          dust_flux_in = saved_state%dust_FLUX_IN(i, c, bid)
          PAR_out = saved_state%PAR_out(i, c, bid)
          ph_prev_3d(:) = saved_state%PH_PREV_3D(i, c, :, bid)
          ph_prev_alt_co2_3d(:) = saved_state%PH_PREV_ALT_CO2_3D(i, c, :, bid)

          ! NOTE(bja, 2015-07) dTracer=0 must come before the "not
          ! lsource_sink check to ensure correct answer when not doing
          ! computations.
          dtracer(:, :) = c0

          if (.not. lsource_sink) then
             !-----------------------------------------------------------------------
             !  exit immediately if computations are not to be performed
             !-----------------------------------------------------------------------
             return
          endif


       !! NOTE(bja, 2015-07) this comment got mis-placed during
       !! refactoring, not sure where it should actually go.
       !! MNL NOTE: removing logical flags triggered by accumulate_tavg_now
       !!           I realize this greatly increases the amount of work done
       !!           in cases where we aren't accumulating, but we need a better
       !!           way to handle this logic and my first pass ignores that fact

       !-----------------------------------------------------------------------
       !  create local copies of model tracers
       !-----------------------------------------------------------------------
          associate(&
               POC => marbl_particulate_share%POC, &
               P_CaCO3 => marbl_particulate_share%P_CaCO3, &
               P_SiO2 => marbl_particulate_share%P_SiO2, &
               dust => marbl_particulate_share%dust, &
               P_iron => marbl_particulate_share%P_iron &
               )

          do k = 1, domain%km
             !write(*, *) 'set_interior loop: ', k, i, c
             call setup_local_tracers(k, domain%land_mask, domain%kmt, &
                  tracer_module(:, k), tracer_local(:, k))

             call setup_local_zooplankton(k, domain%land_mask, domain%kmt, &
                  tracer_module(:, k), zooplankton_cnt, zooplankton, zooplankton_local(:, k))

             call setup_local_autotrophs(k, domain%land_mask, domain%kmt, tracer_module(:, k), &
                  autotroph_cnt, autotrophs, autotroph_local(:, k))

             !  set tracer restore fields
             call ecosys_restore%restore_tracers(ecosys_tracer_cnt, &
                  vert_level=k, x_index=i, y_index=c, block_id=bid, local_data=tracer_local(:, k), &
                  restore_data=restore_local(:, k))

          enddo
          call init_particulate_terms(1, POC, P_CaCO3, P_SiO2, &
               dust, P_iron, &
               QA_dust_def(:), dust_flux_in)

          call compute_carbonate_chemistry(domain%km, bid, &
               domain%land_mask, domain%kmt, &
               domain%temperature(:), domain%salinity(:), &
               tracer_local(:, :), &
               carbonate(:), &
               ph_prev_3d(:), ph_prev_alt_co2_3d(:), &
               zsat_calcite(:), zsat_aragonite(:), &
               co3_calc_anom(:), co3_arag_anom(:))

          do k = 1, domain%km



             call autotroph_consistency_check(autotroph_cnt, autotrophs, autotroph_local(:, k))

             call compute_autotroph_elemental_ratios(k, autotroph_cnt, autotrophs, autotroph_local(:, k), &
                  tracer_local(:, k), autotroph_secondary_species(:, k))

             call compute_photosynthetically_available_radiation(k, autotroph_cnt, &
                  autotroph_local(:, k), &
                  domain%land_mask, domain%kmt, domain%dzt(k), domain%dz(k), &
                  PAR_out, PAR(k))


             call compute_function_scaling(domain%temperature(k), Tfunc(k))

             call compute_Pprime(k, autotroph_cnt, autotrophs, autotroph_local(:, k), &
                  domain%temperature(k), autotroph_secondary_species(:, k))

             call compute_autotroph_uptake(autotroph_cnt, autotrophs, &
                  tracer_local(:, k), &
                  autotroph_secondary_species(:, k))

             call compute_autotroph_photosynthesis(autotroph_cnt, autotrophs, &
                  autotroph_local(:, k), domain%temperature(k), Tfunc(k), &
                  PAR(k)%avg, autotroph_secondary_species(:, k))

             call compute_autotroph_phyto_diatoms (autotroph_cnt, autotrophs, &
                  autotroph_local(:, k), PAR(k)%avg, autotroph_secondary_species(:, k))

             call compute_autotroph_calcification(autotroph_cnt, autotrophs, &
                  autotroph_local(:, k),  domain%temperature(k), autotroph_secondary_species(:, k))

             call compute_autotroph_nfixation(autotroph_cnt, autotrophs, &
                  autotroph_secondary_species(:, k))

             call compute_autotroph_loss(autotroph_cnt, autotrophs, &
                  Tfunc(k), autotroph_secondary_species(:, k))

             call compute_Zprime(k, zooplankton_cnt, zooplankton, zooplankton_local(:, k)%C, &
                  Tfunc(k), zooplankton_secondary_species(:, k))

             call compute_grazing (autotroph_cnt, zooplankton_cnt, grazer_prey_cnt, autotrophs, &
                  Tfunc(k), zooplankton_local(:, k), &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

             call compute_routing (autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

             call compute_dissolved_organic_matter (autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k), &
                  PAR(k)%avg, tracer_local(:, k), &
                  dissolved_organic_matter(k))

             call compute_large_detritus(k, autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k), tracer_local(fe_ind, k), &
                  POC, P_CaCO3, P_SiO2, dust, P_iron, &
                  Fe_scavenge(k), Fe_scavenge_rate(k))

             call compute_particulate_terms(i, c, k, &
                  domain%land_mask, domain%kmt, domain%dzt(k), domain%dz(k), &
                  marbl_particulate_share, &
                  POC, P_CaCO3, P_SiO2, dust, P_iron, &
                  QA_dust_def(k), domain%temperature(k), tracer_local(:, k), &
                  sed_denitrif(k), other_remin(k), lexport_shared_vars, &
                  bid)

             call compute_nitrif(PAR_out, PAR(k)%in, PAR(k)%KPARdz, &
                  tracer_local(nh4_ind, k), nitrif(k))

             call compute_denitrif(tracer_local(o2_ind, k), tracer_local(no3_ind, k), &
                  dissolved_organic_matter(k)%DOC_remin, &
                  POC%remin(k), other_remin(k), sed_denitrif(k), denitrif(k))

             call compute_dtracer_local (autotroph_cnt, zooplankton_cnt, autotrophs, zooplankton, &
                  autotroph_secondary_species(:, k), &
                  zooplankton_secondary_species(:, k), &
                  dissolved_organic_matter(k), &
                  nitrif(k), denitrif(k), sed_denitrif(k), &
                  Fe_scavenge(k) , Fe_scavenge_rate(k), &
                  P_iron%remin(k), POC%remin(k), &
                  P_SiO2%remin(k), P_CaCO3%remin(k), other_remin(k), &
                  restore_local(:, k), &
                  tracer_local(o2_ind, k), &
                  o2_production(k), o2_consumption(k), &
                  dtracer(:, k) )


             ! store_diagnostics_restore
             marbl_diagnostics(k)%restore_diags(:) = restore_local(:, k)

             call store_diagnostics_particulates(k, domain%dz(k), POC, &
                  P_CaCO3, P_SiO2, &
                  dust,  P_iron, &
                  sed_denitrif(k), other_remin(k), &
                  marbl_diagnostics(k)%part_diags(:))

             call store_diagnostics_zooplankton(zooplankton_cnt, &
                  zooplankton_secondary_species(:, k), &
                  marbl_diagnostics(k)%zoo_diags(:, :))

             call store_diagnostics_dissolved_organic_matter(&
                  dissolved_organic_matter(k), &
                  fe_scavenge(k), fe_scavenge_rate(k), &
                  marbl_diagnostics(k)%diags_3d(:))

             call store_diagnostics_carbon_fluxes(&
                  k, domain%kmt, domain%dzt(k), domain%dz(k), zw, &
                  POC, P_CaCO3, &
                  autotroph_cnt, autotrophs, zooplankton_cnt, zooplankton, &
                  dtracer(:, k), marbl_diagnostics(k)%diags_3d(:))

             call store_diagnostics_nitrogen_fluxes(&
                  k, domain%kmt, domain%dzt(k), domain%dz(k), zw, &
                  POC, denitrif(k), sed_denitrif(k), &
                  autotroph_cnt, autotrophs, autotroph_secondary_species(:, k), &
                  zooplankton_cnt, zooplankton, &
                  dtracer(:, k), marbl_diagnostics(k)%diags_3d(:))

             call store_diagnostics_phosphorus_fluxes(&
                  k, domain%kmt, domain%dzt(k), domain%dz(k), zw, &
                  POC, &
                  autotroph_cnt, autotrophs, &
                  zooplankton_cnt, zooplankton, &
                  dtracer(:, k), marbl_diagnostics(k)%diags_3d(:))

             call store_diagnostics_silicon_fluxes(&
                  k, domain%kmt, domain%dzt(k), domain%dz(k), zw, &
                  P_SiO2, &
                  autotroph_cnt, autotrophs, &
                  dtracer(:, k), marbl_diagnostics(k)%diags_3d(:))

             if (lexport_shared_vars) then
                call export_interior_shared_variables(tracer_local(:, k), &
                     carbonate(k), dissolved_organic_matter(k), &
                     QA_dust_def(k), &
                     marbl_interior_share(k))

                call export_zooplankton_shared_variables(zooplankton_cnt, &
                     zooplankton_local(:, k), &
                     zooplankton_secondary_species(:, k), &
                     marbl_zooplankton_share(:, k))

                call export_autotroph_shared_variables(autotroph_cnt, &
                     autotroph_local(:, k), &
                     autotroph_secondary_species(:, k), &
                     marbl_autotroph_share(:, k))

                ! FIXME(bja, 2015-08) need to pull particulate share
                ! out of compute_particulate_terms!
             end if

       ! FIXME(bja, 2015-07) copy column ordered local tracers back to
       ! slab ordered for remaining computations. This will go a way
       ! once the slab --> column reordering is complete.

             saved_state%PAR_out(i, c, bid)               = PAR_out
             saved_state%PH_PREV_3D(i, c, k, bid)         = ph_prev_3d(k)
             saved_state%PH_PREV_ALT_CO2_3D(i, c, k, bid) = ph_prev_alt_co2_3d(k)

       ! FIXME(bja, 2015-07) copy column diags back to slab eventually
       ! gets moved into ecosys_driver!

             if (lexport_shared_vars) then
                ! FIXME(bja, 2015-08) copy column shared data back to slab
                ! eventually gets moved into ecosys_driver!
                call column_interior_share_to_slab_interior_share(i, c, k, bid, &
                     marbl_interior_share(k), ecosys_interior_share(k))

                call column_zooplankton_share_to_slab_zooplankton_share(i, c, k, bid, &
                     marbl_zooplankton_share, ecosys_zooplankton_share(k))


                call column_autotroph_share_to_slab_autotroph_share(i, c, k, bid, &
                     marbl_autotroph_share, ecosys_autotroph_share(k))

                call column_particulate_share_to_slab_particulate_share(i, c, k, bid, &
                     marbl_particulate_share, ecosys_particulate_share(k))
             end if
             if(k<domain%km) THEN
                call update_particulate_terms_from_prior_level(k+1, POC, P_CaCO3, &
                     P_SiO2, dust, P_iron, QA_dust_def(:))
             endif

          end do ! k

          call store_diagnostics_carbonate(carbonate, zsat_calcite,           &
                                           zsat_aragonite, marbl_diagnostics)

          call store_diagnostics_autotrophs(domain, autotrophs,               &
                                            autotroph_secondary_species,      &
                                            marbl_diagnostics)

          call store_diagnostics_autotroph_sums(domain,                       &
                                                autotroph_secondary_species,  &
                                                marbl_diagnostics)

          call store_diagnostics_nitrification(nitrif, denitrif,              &
                                               marbl_diagnostics)

          call store_diagnostics_oxygen(domain, zt, tracer_module(o2_ind, :), &
                                        o2_production, o2_consumption,        &
                                        marbl_diagnostics)

          call store_diagnostics_photosynthetically_available_radiation(PAR,  &
                                                            marbl_diagnostics)

          end associate

  end subroutine ecosys_set_interior

  !***********************************************************************
  !BOP
  ! !IROUTINE: init_particulate_terms
  ! !INTERFACE:

  subroutine init_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, NET_DUST_IN)

    ! !DESCRIPTION:
    !  Set incoming fluxes (put into outgoing flux for first level usage).
    !  Set dissolution length, production fraction and mass terms.
    !
    !  The first 6 arguments are intent(inout) in
    !  order to preserve contents on other blocks.

    use marbl_share_mod, only : column_sinking_particle_type

    ! !INPUT/OUTPUT PARAMETERS:

    integer(int_kind), intent(in) :: k

    type(column_sinking_particle_type), intent(inout) :: &
         POC,          & ! base units = nmol C
         P_CaCO3,      & ! base units = nmol CaCO3
         P_SiO2,       & ! base units = nmol SiO2
         dust,         & ! base units = g
         P_iron          ! base units = nmol Fe

    real (r8), intent(inout) :: QA_dust_def(km) ! incoming deficit in the QA(dust) POC flux

    ! !INPUT PARAMETERS:

    real (r8), intent(in) :: net_dust_in ! dust flux

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  parameters, from Armstrong et al. 2000
    !
    !  July 2002, length scale for excess POC and bSI modified by temperature
    !  Value given here is at Tref of 30 deg. C, JKM
    !-----------------------------------------------------------------------

    POC%diss      = parm_POC_diss   ! diss. length (cm), modified by TEMP
    POC%gamma     = c0              ! not used
    POC%mass      = 12.01_r8        ! molecular weight of POC
    POC%rho       = c0              ! not used

    P_CaCO3%diss  = parm_CaCO3_diss ! diss. length (cm)
    P_CaCO3%gamma = 0.30_r8         ! prod frac -> hard subclass
    P_CaCO3%mass  = 100.09_r8       ! molecular weight of CaCO
    P_CaCO3%rho   = 0.05_r8 * P_CaCO3%mass / POC%mass ! QA mass ratio for CaCO3

    P_SiO2%diss   = parm_SiO2_diss  ! diss. length (cm), modified by TEMP
    P_SiO2%gamma  = 0.030_r8        ! prod frac -> hard subclass
    P_SiO2%mass   = 60.08_r8        ! molecular weight of SiO2
    P_SiO2%rho    = 0.05_r8 * P_SiO2%mass / POC%mass ! QA mass ratio for SiO2

    dust%diss     = 20000.0_r8      ! diss. length (cm)
    dust%gamma    = 0.97_r8         ! prod frac -> hard subclass
    dust%mass     = 1.0e9_r8        ! base units are already grams
    dust%rho      = 0.05_r8 * dust%mass / POC%mass ! QA mass ratio for dust

    P_iron%diss   = 60000.0_r8      ! diss. length (cm) - not used
    P_iron%gamma  = c0              ! prod frac -> hard subclass - not used
    P_iron%mass   = c0              ! not used
    P_iron%rho    = c0              ! not used

    !-----------------------------------------------------------------------
    !  Set incoming fluxes
    !-----------------------------------------------------------------------

    P_CaCO3%sflux_out(k) = c0
    P_CaCO3%hflux_out(k) = c0
    P_CaCO3%sflux_in(k) = P_CaCO3%sflux_out(k)
    P_CaCO3%hflux_in(k) = P_CaCO3%hflux_out(k)

    P_SiO2%sflux_out(k) = c0
    P_SiO2%hflux_out(k) = c0
    P_SiO2%sflux_in(k) = P_SiO2%sflux_out(k)
    P_SiO2%hflux_in(k) = P_SiO2%hflux_out(k)

    if (dust_flux%has_data) then
       dust%sflux_out(k) = (c1 - dust%gamma) * net_dust_in
       dust%hflux_out(k) = dust%gamma * net_dust_in
    else
       dust%sflux_out(k) = c0
       dust%hflux_out(k) = c0
    endif
    dust%sflux_in(k) = dust%sflux_out(k)
    dust%hflux_in(k) = dust%hflux_out(k)

    P_iron%sflux_out(k) = c0
    P_iron%hflux_out(k) = c0
    P_iron%sflux_in(k) = P_iron%sflux_out(k)
    P_iron%hflux_in(k) = P_iron%hflux_out(k)

    !-----------------------------------------------------------------------
    !  Hard POC is QA flux and soft POC is excess POC.
    !-----------------------------------------------------------------------

    POC%sflux_out(k) = c0
    POC%hflux_out(k) = c0
    POC%sflux_in(k) = POC%sflux_out(k)
    POC%hflux_in(k) = POC%hflux_out(k)

    !-----------------------------------------------------------------------
    !  Compute initial QA(dust) POC flux deficit.
    !-----------------------------------------------------------------------

    QA_dust_def(k) = dust%rho * &
         (dust%sflux_out(k) + dust%hflux_out(k))

    !-----------------------------------------------------------------------
    !EOC

  end subroutine init_particulate_terms

  !***********************************************************************

  subroutine update_particulate_terms_from_prior_level(k, POC, P_CaCO3, P_SiO2, dust, P_iron, QA_dust_def)

    use marbl_share_mod, only : column_sinking_particle_type

    integer (int_kind), intent(in) :: k ! vertical model level
    type(column_sinking_particle_type), intent(inout) :: POC, P_CaCO3, P_SiO2, dust, P_iron
    real(r8), intent(inout) :: QA_dust_def(km)

    ! NOTE(bja, 2015-04) assume that k == 1 condition was handled by
    ! call to init_particulate_terms()
    if (k > 1) then
       !-----------------------------------------------------------------------
       ! NOTE: incoming fluxes are outgoing fluxes from previous level
       !
       ! initialize loss to sediments = 0
       !-----------------------------------------------------------------------
       call update_sinking_particle_from_prior_level(k, P_CaCO3)

       call update_sinking_particle_from_prior_level(k, P_SiO2)

       call update_sinking_particle_from_prior_level(k, dust)

       call update_sinking_particle_from_prior_level(k, POC)

       call update_sinking_particle_from_prior_level(k, P_iron)

       QA_dust_def(k) = QA_dust_def(k-1)
    end if

  end subroutine update_particulate_terms_from_prior_level

  !***********************************************************************

  subroutine update_sinking_particle_from_prior_level(k, sinking_particle)

    use marbl_share_mod, only : column_sinking_particle_type

    integer (int_kind), intent(in) :: k
    type(column_sinking_particle_type), intent(inout) :: sinking_particle

    ! NOTE(bja, 201504) level k influx is equal to the level k-1 outflux.
    sinking_particle%sflux_out(k) = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_out(k) = sinking_particle%hflux_out(k-1)
    sinking_particle%sflux_in(k)  = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_in(k)  = sinking_particle%hflux_out(k-1)

  end subroutine update_sinking_particle_from_prior_level

  !***********************************************************************
  !BOP
  ! !IROUTINE: compute_particulate_terms
  ! !INTERFACE:

  subroutine compute_particulate_terms(i, j, k, &
       column_land_mask, column_kmt, column_dzt, column_dz, &
       marbl_particulate_share, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, temperature, tracer_local, sed_denitrif, other_remin, &
       lexport_shared_vars, bid)

    ! !DESCRIPTION:
    !  Compute outgoing fluxes and remineralization terms. Assumes that
    !  production terms have been set. Incoming fluxes are assumed to be the
    !  outgoing fluxes from the previous level.
    !
    !  It is assumed that there is no production of dust.
    !
    !  Instantaneous remineralization in the bottom cell is implemented by
    !  setting the outgoing flux to zero.
    !
    !  For POC, the hard subclass is the POC flux qualitatively associated
    !  with the ballast flux. The soft subclass is the excess POC flux.
    !
    !  Remineralization for the non-iron particulate pools is computing
    !  by first computing the outgoing flux and then computing the
    !  remineralization from conservation, i.e.
    !     flux_in - flux_out + prod * dz - remin * dz == 0.
    !
    !  For iron, remineralization is first computed from POC remineralization
    !  and then flux_out is computed from conservation. If the resulting
    !  flux_out is negative or should be zero because of the sea floor, the
    !  remineralization is adjusted.
    !  Note: all the sinking iron is in the P_iron%sflux pool, hflux Fe not
    !        explicitly tracked, it is assumed that total iron remin is
    !        proportional to total POC remin.
    !
    !  Based upon Armstrong et al. 2000
    !
    !  July 2002, added temperature effect on remin length scale of
    !  excess POC (all soft POM& Iron) and on SiO2.
    !  new variable passed into ballast, Tfunc, main Temperature function
    !  computed in ecosystem routine.  scaling factor for dissolution
    !  of excess POC, Fe, and Bsi now varies with location (f(temperature)).
    !
    !  Added diffusive iron flux from sediments at depths < 1100m,
    !  based on Johnson et al., 1999, value of 5 umolFe/m2/day,
    !      this value too high, using 2 umolFe/m2/day here
    !
    !  Allow hard fraction of ballast to remin with long length scale 40, 000m
    !     thus ~ 10% of hard ballast remins over 4000m water column.
    !
    !  Sinking dust flux is decreased by assumed instant solubility/dissolution
    !     at ocean surface from the parm_Fe_bioavail.
    !
    !  Modified to allow different Q10 factors for soft POM and bSI remin,
    !  water TEMP is now passed in instead of Tfunc (1/2005, JKM)

    ! !USES:

    use marbl_share_mod, only : sinking_particle
    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_share_mod, only : marbl_particulate_share_type
    use marbl_parms, only : Tref
    use constants, only : T0_Kelvin

#ifdef CCSMCOUPLED
    use shr_sys_mod, only: shr_sys_abort
#endif

    ! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: i, j, k ! vertical model level
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt
    real (r8), intent(in) :: column_dzt, column_dz

    real (r8), intent(in) :: temperature ! temperature for scaling functions bsi%diss

    real (r8), dimension(ecosys_tracer_cnt), intent(in) :: tracer_local ! local copies of model tracer concentrations

    logical (log_kind), intent(in) :: &
         lexport_shared_vars ! flag to save shared_vars or not

    integer (int_kind), intent(in) :: bid ! block info for the current block

    ! !INPUT/OUTPUT PARAMETERS:

    type(column_sinking_particle_type), intent(inout) :: POC ! base units = nmol C

    type(column_sinking_particle_type), intent(inout) :: P_CaCO3 ! base units = nmol CaCO3
    type(column_sinking_particle_type), intent(inout) :: P_SiO2 ! base units = nmol SiO2
    type(column_sinking_particle_type), intent(inout) :: dust ! base units = g
    type(column_sinking_particle_type), intent(inout) :: P_iron ! base units = nmol Fe

    real (r8), intent(inout) :: QA_dust_def     ! incoming deficit in the QA(dust) POC flux

    real (r8), intent(out) :: sed_denitrif    ! sedimentary denitrification (umolN/cm^2/s)
    real (r8), intent(out) :: other_remin     ! sedimentary remin not due to oxic or denitrification

    type(marbl_particulate_share_type), intent(inout) :: marbl_particulate_share

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    real (r8) :: poc_diss, & ! diss. length used (cm)
         sio2_diss, & ! diss. length varies spatially with O2
         caco3_diss, &
         dust_diss

    character(*), parameter :: &
         subname = 'ecosys_mod:compute_particulate_terms'

    real (r8) :: TfuncS  ! temperature scaling from soft POM remin

    real (r8) :: &
         DECAY_Hard,         & ! scaling factor for dissolution of Hard Ballast
         DECAY_HardDust        ! scaling factor for dissolution of Hard dust

    real (r8) :: &
         decay_POC_E,        & ! scaling factor for dissolution of excess POC
         decay_SiO2,         & ! scaling factor for dissolution of SiO2
         decay_CaCO3,        & ! scaling factor for dissolution of CaCO3
         decay_dust,         & ! scaling factor for dissolution of dust
         POC_PROD_avail,     & ! POC production available for excess POC flux
         new_QA_dust_def,    & ! outgoing deficit in the QA(dust) POC flux
         scalelength,        & ! used to scale dissolution length scales
         flux, flux_alt,     & ! temp variables used to update sinking flux
         dz_loc, dzr_loc       ! dz, dzr at a particular i, j location

    integer (int_kind) :: n     ! loop indices

    logical (log_kind) :: poc_error   ! POC error flag


          associate(                                                                          &
               O2_loc   => tracer_local(o2_ind),       &
               NO3_loc => tracer_local(no3_ind),     &
               POC_PROD_avail_fields    => marbl_particulate_share%POC_PROD_avail_fields,    & ! IN/OUT
               decay_POC_E_fields       => marbl_particulate_share%decay_POC_E_fields,       & ! IN/OUT
               decay_CaCO3_fields       => marbl_particulate_share%decay_CaCO3_fields,       & ! IN/OUT
               poc_diss_fields          => marbl_particulate_share%poc_diss_fields,          & ! IN/OUT
               caco3_diss_fields        => marbl_particulate_share%caco3_diss_fields,        & ! IN/OUT
               P_CaCO3_sflux_out_fields => marbl_particulate_share%P_CaCO3_sflux_out_fields, & ! IN/OUT
               P_CaCO3_hflux_out_fields => marbl_particulate_share%P_CaCO3_hflux_out_fields, & ! IN/OUT
               POC_sflux_out_fields     => marbl_particulate_share%POC_sflux_out_fields,     & ! IN/OUT
               POC_hflux_out_fields     => marbl_particulate_share%POC_hflux_out_fields,     & ! IN/OUT
               POC_remin_fields         => marbl_particulate_share%POC_remin_fields,         & ! IN/OUT
               P_CaCO3_remin_fields     => marbl_particulate_share%P_CaCO3_remin_fields,     & ! IN/OUT
               DECAY_Hard_fields        => marbl_particulate_share%DECAY_Hard_fields         & ! IN/OUT
               )

          !-----------------------------------------------------------------------
          !  initialize local copy of percent sed
          !-----------------------------------------------------------------------
          sed_denitrif = c0
          other_remin = c0

          !-----------------------------------------------------------------------
          !  compute scalelength and decay factors
          !-----------------------------------------------------------------------

          if (zw(k) < parm_scalelen_z(1)) then
             scalelength = parm_scalelen_vals(1)
          else if (zw(k) >= parm_scalelen_z(size(parm_scalelen_z))) then
             scalelength = parm_scalelen_vals(size(parm_scalelen_z))
          else
             do n = 2, size(parm_scalelen_z)
                if (zw(k) < parm_scalelen_z(n)) then
                   scalelength = parm_scalelen_vals(n-1) &
                        + (parm_scalelen_vals(n) - parm_scalelen_vals(n-1)) &
                        * (zw(k) - parm_scalelen_z(n-1))/(parm_scalelen_z(n) - parm_scalelen_z(n-1))
                   exit
                endif
             end do
          endif

          if (partial_bottom_cells) then
             DECAY_Hard     = exp(-column_dzt / 4.0e6_r8)
             DECAY_HardDust = exp(-column_dzt / 1.2e7_r8)
          else
             DECAY_Hard     = exp(-column_dz / 4.0e6_r8)
             DECAY_HardDust = exp(-column_dz / 1.2e7_r8)
          endif

          !----------------------------------------------------------------------
          !   Tref = 30.0 reference temperature (deg. C)
          !-----------------------------------------------------------------------
          TfuncS = 1.5_r8**(((temperature + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

          poc_error = .false.
          dz_loc = column_dz

          if (column_land_mask .and. k <= column_kmt) then

             if (partial_bottom_cells) then
                dz_loc = column_dzt
             endif
             dzr_loc = c1 / dz_loc

             poc_diss = POC%diss
             sio2_diss = P_SiO2%diss
             caco3_diss = P_CaCO3%diss
             dust_diss = dust%diss

             !-----------------------------------------------------------------------
             !  increase POC diss length scale where O2 concentrations are low
             !-----------------------------------------------------------------------

             if ((O2_loc >= 5.0_r8) .and. (O2_loc < 40.0_r8)) then
                poc_diss = POC%diss*(c1+(3.3_r8-c1)*(40.0_r8 - O2_loc)/35.0_r8)
             else if (O2_loc < 5.0_r8) then
                poc_diss = POC%diss * 3.3_r8
             endif

             !-----------------------------------------------------------------------
             !  apply scalelength factor to length scales
             !-----------------------------------------------------------------------

             poc_diss = scalelength * poc_diss
             sio2_diss = scalelength * sio2_diss
             caco3_diss = scalelength * caco3_diss
             dust_diss = scalelength * dust_diss

             !-----------------------------------------------------------------------
             !  apply temperature dependence to sio2_diss length scale
             !-----------------------------------------------------------------------

             sio2_diss = sio2_diss / TfuncS

             !-----------------------------------------------------------------------
             !  decay_POC_E and decay_SiO2 set locally, modified by O2
             !-----------------------------------------------------------------------

             decay_POC_E = exp(-dz_loc / poc_diss)
             decay_SiO2  = exp(-dz_loc / sio2_diss)
             decay_CaCO3 = exp(-dz_loc / caco3_diss)
             decay_dust  = exp(-dz_loc / dust_diss)

             !-----------------------------------------------------------------------
             !  Set outgoing fluxes for non-iron pools.
             !  The outoing fluxes for ballast materials are from the
             !  solution of the coresponding continuous ODE across the model
             !  level. The ODE has a constant source term and linear decay.
             !  It is assumed that there is no sub-surface dust production.
             !-----------------------------------------------------------------------

             P_CaCO3%sflux_out(k) = P_CaCO3%sflux_in(k) * decay_CaCO3 + &
                  P_CaCO3%prod(k) * ((c1 - P_CaCO3%gamma) * (c1 - decay_CaCO3) &
                  * caco3_diss)

             P_CaCO3%hflux_out(k) = P_CaCO3%hflux_in(k) * DECAY_Hard + &
                  P_CaCO3%prod(k) * (P_CaCO3%gamma * dz_loc)

             P_SiO2%sflux_out(k) = P_SiO2%sflux_in(k) * decay_SiO2 + &
                  P_SiO2%prod(k) * ((c1 - P_SiO2%gamma) * (c1 - decay_SiO2) &
                  * sio2_diss)

             P_SiO2%hflux_out(k) = P_SiO2%hflux_in(k) * DECAY_Hard + &
                  P_SiO2%prod(k) * (P_SiO2%gamma * dz_loc)

             dust%sflux_out(k) = dust%sflux_in(k) * decay_dust

             dust%hflux_out(k) = dust%hflux_in(k) * DECAY_HardDust

             !-----------------------------------------------------------------------
             !  Compute how much POC_PROD is available for deficit reduction
             !  and excess POC flux after subtracting off fraction of non-dust
             !  ballast production from net POC_PROD.
             !-----------------------------------------------------------------------

             POC_PROD_avail = POC%prod(k) - &
                  P_CaCO3%rho * P_CaCO3%prod(k) - &
                  P_SiO2%rho * P_SiO2%prod(k)

             !-----------------------------------------------------------------------
             !  Check for POC production bounds violations
             !-----------------------------------------------------------------------

             if (POC_PROD_avail < c0) then
                poc_error = .true.
             endif

             !-----------------------------------------------------------------------
             !  Compute 1st approximation to new QA_dust_def, the QA_dust
             !  deficit leaving the cell. Ignore POC_PROD_avail at this stage.
             !-----------------------------------------------------------------------

             if (QA_dust_def > 0) then
                new_QA_dust_def = QA_dust_def * &
                     (dust%sflux_out(k) + dust%hflux_out(k)) / &
                     (dust%sflux_in(k) + dust%hflux_in(k))
             else
                new_QA_dust_def = c0
             endif

             !-----------------------------------------------------------------------
             !  Use POC_PROD_avail to reduce new_QA_dust_def.
             !-----------------------------------------------------------------------

             if (new_QA_dust_def > c0) then
                new_QA_dust_def = new_QA_dust_def - POC_PROD_avail * dz_loc
                if (new_QA_dust_def < c0) then
                   POC_PROD_avail = -new_QA_dust_def * dzr_loc
                   new_QA_dust_def = c0
                else
                   POC_PROD_avail = c0
                endif
             endif

             QA_dust_def = new_QA_dust_def

             ! Save certain fields for use by other modules
             if (lexport_shared_vars) then
                POC_PROD_avail_fields(k) = POC_PROD_avail
                decay_POC_E_fields(k)    = decay_POC_E
                decay_CaCO3_fields(k)    = decay_CaCO3
                poc_diss_fields(k)       = poc_diss
                caco3_diss_fields(k)     = caco3_diss
             endif

             !-----------------------------------------------------------------------
             !  Compute outgoing POC fluxes. QA POC flux is computing using
             !  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
             !  and no production occured, then no QA POC flux goes out. This
             !  shortcut is present to avoid roundoff cancellation errors from
             !  the dust%rho * dust_flux_out - QA_dust_def computation.
             !  Any POC_PROD_avail still remaining goes into excess POC flux.
             !-----------------------------------------------------------------------

             if (POC%hflux_in(k) == c0 .and. POC%prod(k) == c0) then
                POC%hflux_out(k) = c0
             else
                POC%hflux_out(k) = P_CaCO3%rho * &
                     (P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)) + &
                     P_SiO2%rho * &
                     (P_SiO2%sflux_out(k) + P_SiO2%hflux_out(k)) + &
                     dust%rho * &
                     (dust%sflux_out(k) + dust%hflux_out(k)) - &
                     new_QA_dust_def
                POC%hflux_out(k) = max(POC%hflux_out(k), c0)
             endif

             POC%sflux_out(k) = POC%sflux_in(k) * decay_POC_E + &
                  POC_PROD_avail *((c1 - decay_POC_E) * &
                  poc_diss)

             !-----------------------------------------------------------------------
             !  Compute remineralization terms. It is assumed that there is no
             !  sub-surface dust production.
             !-----------------------------------------------------------------------

             P_CaCO3%remin(k) = P_CaCO3%prod(k) + &
                  ((P_CaCO3%sflux_in(k) - P_CaCO3%sflux_out(k)) + &
                  (P_CaCO3%hflux_in(k) - P_CaCO3%hflux_out(k))) * dzr_loc

             P_SiO2%remin(k) = P_SiO2%prod(k) + &
                  ((P_SiO2%sflux_in(k) - P_SiO2%sflux_out(k)) + &
                  (P_SiO2%hflux_in(k) - P_SiO2%hflux_out(k))) * dzr_loc

             POC%remin(k) = POC%prod(k) + &
                  ((POC%sflux_in(k) - POC%sflux_out(k)) + &
                  (POC%hflux_in(k) - POC%hflux_out(k))) * dzr_loc

             dust%remin(k) = &
                  ((dust%sflux_in(k) - dust%sflux_out(k)) + &
                  (dust%hflux_in(k) - dust%hflux_out(k))) * dzr_loc

             !-----------------------------------------------------------------------
             !  Compute iron remineralization and flux out.
             !-----------------------------------------------------------------------

             if (POC%sflux_in(k) + POC%hflux_in(k) == c0) then
                P_iron%remin(k) = (POC%remin(k) * parm_Red_Fe_C)
             else
                P_iron%remin(k) = (POC%remin(k) * &
                     (P_iron%sflux_in(k) + P_iron%hflux_in(k)) / &
                     (POC%sflux_in(k) + POC%hflux_in(k)))
             endif
             P_iron%remin(k) = P_iron%remin(k) +                &
                  (P_iron%sflux_in(k) * 1.5e-5_r8)

             P_iron%sflux_out(k) = P_iron%sflux_in(k) + dz_loc * &
                  ((c1 - P_iron%gamma) * P_iron%prod(k) - P_iron%remin(k))

             if (P_iron%sflux_out(k) < c0) then
                P_iron%sflux_out(k) = c0
                P_iron%remin(k) = P_iron%sflux_in(k) * dzr_loc + &
                     (c1 - P_iron%gamma) * P_iron%prod(k)
             endif

             !-----------------------------------------------------------------------
             !  Compute iron release from dust remin/dissolution
             !
             !  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
             !                      gFe     molFe     nmolFe
             !  Also add in Fe source from sediments if applicable to this cell.
             !-----------------------------------------------------------------------


             P_iron%remin(k) = P_iron%remin(k) &
                  + dust%remin(k) * dust_to_Fe &
                  + (FESEDFLUX(i, j, k, bid) * dzr_loc)

             P_iron%hflux_out(k) = P_iron%hflux_in(k)

          else
             P_CaCO3%sflux_out(k) = c0
             P_CaCO3%hflux_out(k) = c0
             P_CaCO3%remin(k) = c0

             P_SiO2%sflux_out(k) = c0
             P_SiO2%hflux_out(k) = c0
             P_SiO2%remin(k) = c0

             dust%sflux_out(k) = c0
             dust%hflux_out(k) = c0
             dust%remin(k) = c0

             POC%sflux_out(k) = c0
             POC%hflux_out(k) = c0
             POC%remin(k) = c0

             P_iron%sflux_out(k) = c0
             P_iron%hflux_out(k) = c0
             P_iron%remin(k) = c0
          endif

          ! Save some fields for use by other modules before setting outgoing fluxes to 0.0 in bottom cell below
          if (lexport_shared_vars) then
             P_CaCO3_sflux_out_fields(k) = P_CaCO3%sflux_out(k)
             P_CaCO3_hflux_out_fields(k) = P_CaCO3%hflux_out(k)
             POC_sflux_out_fields(k)     = POC%sflux_out(k)
             POC_hflux_out_fields(k)     = POC%hflux_out(k)
             POC_remin_fields(k)         = POC%remin(k)
             P_CaCO3_remin_fields(k)     = P_CaCO3%remin(k)
             DECAY_Hard_fields(k)        = DECAY_Hard
          endif

          !-----------------------------------------------------------------------
          !  Bottom Sediments Cell?
          !  If so compute sedimentary burial and denitrification N losses.
          !  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
          !  other_remin estimates organic matter remineralized in the sediments
          !      by the processes other than oxic remin and denitrification (SO4 and CO2,
          !      etc..)
          !      based on Soetaert et al., 1996, varies between 10% and 50%
          !      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
          !      other_remin is 50% above this high flux value,
          !      In special case where bottom O2 has been depleted to < 1.0 uM,
          !               all sedimentary remin is due to DENITRIFICATION + other_remin
          !  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
          !  Bsi preservation in sediments based on
          !     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
          !  Calcite is preserved in sediments above the lysocline, dissolves below.
          !       Here a constant depth is used for lysocline.
          !-----------------------------------------------------------------------

          POC%sed_loss(k)     = c0
          P_SiO2%sed_loss(k)  = c0
          P_CaCO3%sed_loss(k) = c0
          P_iron%sed_loss(k)  = c0
          dust%sed_loss(k)    = c0

          if (column_land_mask .and. (k == column_kmt)) then

             flux = POC%sflux_out(k) + POC%hflux_out(k)

             if (flux > c0) then
                flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day

                POC%sed_loss(k) = flux * min(0.8_r8, parm_POMbury &
                     * (0.013_r8 + 0.53_r8 * flux_alt*flux_alt / (7.0_r8 + flux_alt)**2))

                sed_denitrif = dzr_loc * flux &
                     * (0.06_r8 + 0.19_r8 * 0.99_r8**(O2_loc-NO3_loc))

                flux_alt = flux*1.0e-6_r8*spd*365.0_r8 ! convert to mmol/cm^2/year
                other_remin = dzr_loc &
                     * min(min(0.1_r8 + flux_alt, 0.5_r8) * (flux - POC%sed_loss(k)), &
                     (flux - POC%sed_loss(k) - (sed_denitrif*dz_loc*denitrif_C_N)))

                !----------------------------------------------------------------------------------
                !              if bottom water O2 is depleted, assume all remin is denitrif + other
                !----------------------------------------------------------------------------------

                if (O2_loc < c1) then
                   other_remin = dzr_loc * &
                        (flux - POC%sed_loss(k) - (sed_denitrif*dz_loc*denitrif_C_N))
                endif

             endif

             flux = P_SiO2%sflux_out(k) + P_SiO2%hflux_out(k)
             flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day
             ! first compute burial efficiency, then compute loss to sediments
             if (flux_alt > c2) then
                P_SiO2%sed_loss(k) = 0.2_r8
             else
                P_SiO2%sed_loss(k) = 0.04_r8
             endif
             P_SiO2%sed_loss(k) = flux * parm_BSIbury * P_SiO2%sed_loss(k)

             if (zw(k) < 3300.0e2_r8) then
                flux = P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)
                P_CaCO3%sed_loss(k) = flux
             endif

             !----------------------------------------------------------------------------------
             !  Update sinking fluxes and remin fluxes, accounting for sediments.
             !  flux used to hold sinking fluxes before update.
             !----------------------------------------------------------------------------------

             flux = P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)
             if (flux > c0) then
                P_CaCO3%remin(k) = P_CaCO3%remin(k) &
                     + ((flux - P_CaCO3%sed_loss(k)) * dzr_loc)
             endif

             flux = P_SiO2%sflux_out(k) + P_SiO2%hflux_out(k)
             if (flux > c0) then
                P_SiO2%remin(k) = P_SiO2%remin(k) &
                     + ((flux - P_SiO2%sed_loss(k)) * dzr_loc)
             endif

             flux = POC%sflux_out(k) + POC%hflux_out(k)
             if (flux > c0) then
                POC%remin(k) = POC%remin(k) &
                     + ((flux - POC%sed_loss(k)) * dzr_loc)
             endif

             !-----------------------------------------------------------------------
             !   Remove all Piron and dust that hits bottom, sedimentary Fe source
             !        accounted for by FESEDFLUX elsewhere.
             !-----------------------------------------------------------------------

             flux = (P_iron%sflux_out(k) + P_iron%hflux_out(k))
             if (flux > c0) then
                P_iron%sed_loss(k) = flux
             endif

             dust%sed_loss(k) = dust%sflux_out(k) + dust%hflux_out(k)

             !-----------------------------------------------------------------------
             !   Set all outgoing fluxes to 0.0
             !-----------------------------------------------------------------------

             if (k == column_kmt) then
                P_CaCO3%sflux_out(k) = c0
                P_CaCO3%hflux_out(k) = c0

                P_SiO2%sflux_out(k) = c0
                P_SiO2%hflux_out(k) = c0

                dust%sflux_out(k) = c0
                dust%hflux_out(k) = c0

                POC%sflux_out(k) = c0
                POC%hflux_out(k) = c0

                P_iron%sflux_out(k) = c0
                P_iron%hflux_out(k) = c0
             endif

          endif

          if (poc_error) then
             call shr_sys_abort(subname /&
                  &/ ': mass ratio of ballast ' /&
                  &/ 'production exceeds POC production')
          endif

          end associate

  end subroutine compute_particulate_terms

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_init_sflux
  ! !INTERFACE:

  subroutine ecosys_init_sflux(saved_state)

    ! !DESCRIPTION:
    !  Initialize surface flux computations for ecosys tracer module.
    !
    ! !REVISION HISTORY:
    !  same as module
    use marbl_interface_types, only : marbl_saved_state_type
    type(marbl_saved_state_type), intent(inout) :: saved_state

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_init_sflux'

    logical (log_kind) :: &
         luse_INTERP_WORK     ! does INTERP_WORK need to be allocated

    integer (int_kind) :: &
         n,                 & ! index for looping over tracers
         iblock               ! index for looping over blocks

    real (r8), dimension (nx_block, ny_block) :: WORK

    real (r8), dimension (nx_block, ny_block, 12, max_blocks_clinic), target :: &
         WORK_READ            ! temporary space to read in fields

    !-----------------------------------------------------------------------

    luse_INTERP_WORK = .false.

    !-----------------------------------------------------------------------
    !  read gas flux forcing (if required)
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
         gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

       luse_INTERP_WORK = .true.

       !-----------------------------------------------------------------------
       !  first, read ice file
       !-----------------------------------------------------------------------

       allocate(fice_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(fice_file%input%filename) == 'unknown') &
            fice_file%input%filename = gas_flux_forcing_file

       call read_field(fice_file%input%file_fmt, &
            fice_file%input%filename, &
            fice_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             fice_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  fice_file%DATA(:, :, iblock, 1, n) = c0
             fice_file%DATA(:, :, iblock, 1, n) = &
                  fice_file%DATA(:, :, iblock, 1, n) * fice_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(fice_file%data_time, &
            fice_file%data_inc, fice_file%interp_type, &
            fice_file%data_next, fice_file%data_time_min_loc, &
            fice_file%data_update, fice_file%data_type)

       !-----------------------------------------------------------------------
       !  next, read piston velocity file
       !-----------------------------------------------------------------------

       allocate(xkw_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(xkw_file%input%filename) == 'unknown') &
            xkw_file%input%filename = gas_flux_forcing_file

       call read_field(xkw_file%input%file_fmt, &
            xkw_file%input%filename, &
            xkw_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             xkw_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  xkw_file%DATA(:, :, iblock, 1, n) = c0
             xkw_file%DATA(:, :, iblock, 1, n) = &
                  xkw_file%DATA(:, :, iblock, 1, n) * xkw_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(xkw_file%data_time, &
            xkw_file%data_inc, xkw_file%interp_type, &
            xkw_file%data_next, xkw_file%data_time_min_loc, &
            xkw_file%data_update, xkw_file%data_type)

       !-----------------------------------------------------------------------
       !  last, read atmospheric pressure file
       !-----------------------------------------------------------------------

       allocate(ap_file%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(ap_file%input%filename) == 'unknown') &
            ap_file%input%filename = gas_flux_forcing_file

       call read_field(ap_file%input%file_fmt, &
            ap_file%input%filename, &
            ap_file%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             ap_file%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  ap_file%DATA(:, :, iblock, 1, n) = c0
             ap_file%DATA(:, :, iblock, 1, n) = &
                  ap_file%DATA(:, :, iblock, 1, n) * ap_file%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(ap_file%data_time, &
            ap_file%data_inc, ap_file%interp_type, &
            ap_file%data_next, ap_file%data_time_min_loc, &
            ap_file%data_update, ap_file%data_type)

    endif

    !-----------------------------------------------------------------------
    !  load dust flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(dust_flux%input%filename) /= 'none' .and. &
         trim(dust_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dust_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(dust_flux%input%filename) == 'unknown') &
            dust_flux%input%filename = gas_flux_forcing_file

       call read_field(dust_flux%input%file_fmt, &
            dust_flux%input%filename, &
            dust_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dust_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dust_flux%DATA(:, :, iblock, 1, n) = c0
             dust_flux%DATA(:, :, iblock, 1, n) = &
                  dust_flux%DATA(:, :, iblock, 1, n) * dust_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dust_flux%data_time, &
            dust_flux%data_inc, dust_flux%interp_type, &
            dust_flux%data_next, dust_flux%data_time_min_loc, &
            dust_flux%data_update, dust_flux%data_type)

       dust_flux%has_data = .true.
    else
       dust_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  load iron flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(iron_flux%input%filename) /= 'none' .and. &
         trim(iron_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(iron_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(iron_flux%input%filename) == 'unknown') &
            iron_flux%input%filename = gas_flux_forcing_file

       call read_field(iron_flux%input%file_fmt, &
            iron_flux%input%filename, &
            iron_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             iron_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  iron_flux%DATA(:, :, iblock, 1, n) = c0
             iron_flux%DATA(:, :, iblock, 1, n) = &
                  iron_flux%DATA(:, :, iblock, 1, n) * iron_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(iron_flux%data_time, &
            iron_flux%data_inc, iron_flux%interp_type, &
            iron_flux%data_next, iron_flux%data_time_min_loc, &
            iron_flux%data_update, iron_flux%data_type)

       iron_flux%has_data = .true.
    else
       iron_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  load nox & noy flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(ndep_data_type) /= 'none' .and. &
         trim(ndep_data_type) /= 'monthly-calendar' .and. &
         trim(ndep_data_type) /= 'shr_stream') then
       call document(subname, 'ndep_data_type', ndep_data_type)
       call exit_POP(sigAbort, 'unknown ndep_data_type')
    endif

    if (trim(ndep_data_type) == 'monthly-calendar' .and. &
         trim(nox_flux_monthly%input%filename) /= 'none' .and. &
         trim(nox_flux_monthly%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(nox_flux_monthly%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(nox_flux_monthly%input%filename) == 'unknown') &
            nox_flux_monthly%input%filename = gas_flux_forcing_file

       call read_field(nox_flux_monthly%input%file_fmt, &
            nox_flux_monthly%input%filename, &
            nox_flux_monthly%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             nox_flux_monthly%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  nox_flux_monthly%DATA(:, :, iblock, 1, n) = c0
             nox_flux_monthly%DATA(:, :, iblock, 1, n) = &
                  nox_flux_monthly%DATA(:, :, iblock, 1, n) * nox_flux_monthly%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(nox_flux_monthly%data_time, &
            nox_flux_monthly%data_inc, nox_flux_monthly%interp_type, &
            nox_flux_monthly%data_next, nox_flux_monthly%data_time_min_loc, &
            nox_flux_monthly%data_update, nox_flux_monthly%data_type)

       nox_flux_monthly%has_data = .true.
    else
       nox_flux_monthly%has_data = .false.
    endif

    if (trim(ndep_data_type) == 'monthly-calendar' .and. &
         trim(nhy_flux_monthly%input%filename) /= 'none' .and. &
         trim(nhy_flux_monthly%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(nhy_flux_monthly%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))
       if (trim(nhy_flux_monthly%input%filename) == 'unknown') &
            nhy_flux_monthly%input%filename = gas_flux_forcing_file

       call read_field(nhy_flux_monthly%input%file_fmt, &
            nhy_flux_monthly%input%filename, &
            nhy_flux_monthly%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             nhy_flux_monthly%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  nhy_flux_monthly%DATA(:, :, iblock, 1, n) = c0
             nhy_flux_monthly%DATA(:, :, iblock, 1, n) = &
                  nhy_flux_monthly%DATA(:, :, iblock, 1, n) * nhy_flux_monthly%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(nhy_flux_monthly%data_time, &
            nhy_flux_monthly%data_inc, nhy_flux_monthly%interp_type, &
            nhy_flux_monthly%data_next, nhy_flux_monthly%data_time_min_loc, &
            nhy_flux_monthly%data_update, nhy_flux_monthly%data_type)

       nhy_flux_monthly%has_data = .true.
    else
       nhy_flux_monthly%has_data = .false.
    endif

    !-----------------------------------------------------------------------

#ifndef CCSMCOUPLED
    if (trim(ndep_data_type) == 'shr_stream') then
       call document(subname, 'ndep_data_type', ndep_data_type)
       call exit_POP(sigAbort, &
            'shr_stream option only supported when CCSMCOUPLED is defined')
    endif
#endif

    !-----------------------------------------------------------------------
    !  load river nutrient flux fields (if required)
    !-----------------------------------------------------------------------

    if (trim(din_riv_flux%input%filename) /= 'none' .and. &
         trim(din_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(din_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(din_riv_flux%input%file_fmt, &
            din_riv_flux%input%filename, &
            din_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             din_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  din_riv_flux%DATA(:, :, iblock, 1, n) = c0
             din_riv_flux%DATA(:, :, iblock, 1, n) = &
                  din_riv_flux%DATA(:, :, iblock, 1, n) * din_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(din_riv_flux%data_time, &
            din_riv_flux%data_inc, din_riv_flux%interp_type, &
            din_riv_flux%data_next, din_riv_flux%data_time_min_loc, &
            din_riv_flux%data_update, din_riv_flux%data_type)

       din_riv_flux%has_data = .true.
    else
       din_riv_flux%has_data = .false.
    endif


    if (trim(dip_riv_flux%input%filename) /= 'none' .and. &
         trim(dip_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dip_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dip_riv_flux%input%file_fmt, &
            dip_riv_flux%input%filename, &
            dip_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dip_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dip_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dip_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dip_riv_flux%DATA(:, :, iblock, 1, n) * dip_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dip_riv_flux%data_time, &
            dip_riv_flux%data_inc, dip_riv_flux%interp_type, &
            dip_riv_flux%data_next, dip_riv_flux%data_time_min_loc, &
            dip_riv_flux%data_update, dip_riv_flux%data_type)

       dip_riv_flux%has_data = .true.
    else
       dip_riv_flux%has_data = .false.
    endif


    if (trim(don_riv_flux%input%filename) /= 'none' .and. &
         trim(don_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(don_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(don_riv_flux%input%file_fmt, &
            don_riv_flux%input%filename, &
            don_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             don_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  don_riv_flux%DATA(:, :, iblock, 1, n) = c0
             don_riv_flux%DATA(:, :, iblock, 1, n) = &
                  don_riv_flux%DATA(:, :, iblock, 1, n) * don_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(don_riv_flux%data_time, &
            don_riv_flux%data_inc, don_riv_flux%interp_type, &
            don_riv_flux%data_next, don_riv_flux%data_time_min_loc, &
            don_riv_flux%data_update, don_riv_flux%data_type)

       don_riv_flux%has_data = .true.
    else
       don_riv_flux%has_data = .false.
    endif


    if (trim(dop_riv_flux%input%filename) /= 'none' .and. &
         trim(dop_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dop_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dop_riv_flux%input%file_fmt, &
            dop_riv_flux%input%filename, &
            dop_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dop_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dop_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dop_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dop_riv_flux%DATA(:, :, iblock, 1, n) * dop_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dop_riv_flux%data_time, &
            dop_riv_flux%data_inc, dop_riv_flux%interp_type, &
            dop_riv_flux%data_next, dop_riv_flux%data_time_min_loc, &
            dop_riv_flux%data_update, dop_riv_flux%data_type)

       dop_riv_flux%has_data = .true.
    else
       dop_riv_flux%has_data = .false.
    endif

    if (trim(dsi_riv_flux%input%filename) /= 'none' .and. &
         trim(dsi_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dsi_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dsi_riv_flux%input%file_fmt, &
            dsi_riv_flux%input%filename, &
            dsi_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dsi_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dsi_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dsi_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dsi_riv_flux%DATA(:, :, iblock, 1, n) * dsi_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dsi_riv_flux%data_time, &
            dsi_riv_flux%data_inc, dsi_riv_flux%interp_type, &
            dsi_riv_flux%data_next, dsi_riv_flux%data_time_min_loc, &
            dsi_riv_flux%data_update, dsi_riv_flux%data_type)

       dsi_riv_flux%has_data = .true.
    else
       dsi_riv_flux%has_data = .false.
    endif


    if (trim(dfe_riv_flux%input%filename) /= 'none' .and. &
         trim(dfe_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dfe_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dfe_riv_flux%input%file_fmt, &
            dfe_riv_flux%input%filename, &
            dfe_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dfe_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dfe_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dfe_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dfe_riv_flux%DATA(:, :, iblock, 1, n) * dfe_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dfe_riv_flux%data_time, &
            dfe_riv_flux%data_inc, dfe_riv_flux%interp_type, &
            dfe_riv_flux%data_next, dfe_riv_flux%data_time_min_loc, &
            dfe_riv_flux%data_update, dfe_riv_flux%data_type)

       dfe_riv_flux%has_data = .true.
    else
       dfe_riv_flux%has_data = .false.
    endif


    if (trim(dic_riv_flux%input%filename) /= 'none' .and. &
         trim(dic_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(dic_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(dic_riv_flux%input%file_fmt, &
            dic_riv_flux%input%filename, &
            dic_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             dic_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  dic_riv_flux%DATA(:, :, iblock, 1, n) = c0
             dic_riv_flux%DATA(:, :, iblock, 1, n) = &
                  dic_riv_flux%DATA(:, :, iblock, 1, n) * dic_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(dic_riv_flux%data_time, &
            dic_riv_flux%data_inc, dic_riv_flux%interp_type, &
            dic_riv_flux%data_next, dic_riv_flux%data_time_min_loc, &
            dic_riv_flux%data_update, dic_riv_flux%data_type)

       dic_riv_flux%has_data = .true.
    else
       dic_riv_flux%has_data = .false.
    endif


    if (trim(alk_riv_flux%input%filename) /= 'none' .and. &
         trim(alk_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(alk_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(alk_riv_flux%input%file_fmt, &
            alk_riv_flux%input%filename, &
            alk_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             alk_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  alk_riv_flux%DATA(:, :, iblock, 1, n) = c0
             alk_riv_flux%DATA(:, :, iblock, 1, n) = &
                  alk_riv_flux%DATA(:, :, iblock, 1, n) * alk_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(alk_riv_flux%data_time, &
            alk_riv_flux%data_inc, alk_riv_flux%interp_type, &
            alk_riv_flux%data_next, alk_riv_flux%data_time_min_loc, &
            alk_riv_flux%data_update, alk_riv_flux%data_type)

       alk_riv_flux%has_data = .true.
    else
       alk_riv_flux%has_data = .false.
    endif

    if (trim(doc_riv_flux%input%filename) /= 'none' .and. &
         trim(doc_riv_flux%input%filename) /= 'unknown') then

       luse_INTERP_WORK = .true.

       allocate(doc_riv_flux%DATA(nx_block, ny_block, max_blocks_clinic, 1, 12))

       call read_field(doc_riv_flux%input%file_fmt, &
            doc_riv_flux%input%filename, &
            doc_riv_flux%input%file_varname, &
            WORK_READ)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             doc_riv_flux%DATA(:, :, iblock, 1, n) = WORK_READ(:, :, n, iblock)
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  doc_riv_flux%DATA(:, :, iblock, 1, n) = c0
             doc_riv_flux%DATA(:, :, iblock, 1, n) = &
                  doc_riv_flux%DATA(:, :, iblock, 1, n) * doc_riv_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

       call find_forcing_times(doc_riv_flux%data_time, &
            doc_riv_flux%data_inc, doc_riv_flux%interp_type, &
            doc_riv_flux%data_next, doc_riv_flux%data_time_min_loc, &
            doc_riv_flux%data_update, doc_riv_flux%data_type)

       doc_riv_flux%has_data = .true.
    else
       doc_riv_flux%has_data = .false.
    endif

    !-----------------------------------------------------------------------
    !  allocate space for interpolate_forcing
    !-----------------------------------------------------------------------

    if (luse_INTERP_WORK) &
         allocate(INTERP_WORK(nx_block, ny_block, max_blocks_clinic, 1))

    !-----------------------------------------------------------------------
    !  load iron PATCH flux fields (if required)
    !-----------------------------------------------------------------------

    if (liron_patch) then

       !maltrud iron patch
       !  assume patch file has same normalization and format as deposition file

       allocate(IRON_PATCH_FLUX(nx_block, ny_block, max_blocks_clinic))

       if (trim(iron_flux%input%filename) == 'unknown') &
            iron_flux%input%filename = gas_flux_forcing_file

       call read_field(iron_flux%input%file_fmt, &
            iron_flux%input%filename, &
            iron_patch_flux_filename, &
            IRON_PATCH_FLUX)

       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1, nblocks_clinic
          do n=1, 12
             where (.not. saved_state%land_mask(:, :, iblock)) &
                  IRON_PATCH_FLUX(:, :, iblock) = c0
             iron_flux%DATA(:, :, iblock, 1, n) = &
                  IRON_PATCH_FLUX(:, :, iblock) * iron_flux%input%scale_factor
          end do
       end do
       !$OMP END PARALLEL DO

    endif

    !-----------------------------------------------------------------------
    !  register and set
    !     fco2, the air-sea co2 gas flux
    !-----------------------------------------------------------------------

    call named_field_register('SFLUX_CO2', sflux_co2_nf_ind)
    !$OMP PARALLEL DO PRIVATE(iblock, WORK)
    do iblock=1, nblocks_clinic
       WORK = c0
       call named_field_set(sflux_co2_nf_ind, iblock, WORK)
    end do
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !  verify running coupled if gas fluxes use coupler forcing
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
         (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv .or. &
         atm_co2_iopt == atm_co2_iopt_drv_prog .or. &
         atm_co2_iopt == atm_co2_iopt_drv_diag) .and. &
         .not. registry_match('lcoupled')) then
       call exit_POP(sigAbort, 'ecosys_init: ecosys module requires the ' /&
            &/ 'flux coupler when gas_flux_forcing_opt=drv')
    endif

    !-----------------------------------------------------------------------
    !  get named field index for atmospheric CO2, if required
    !-----------------------------------------------------------------------

    if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_prog) then
       call named_field_get_index('ATM_CO2_PROG', atm_co2_nf_ind, &
            exit_on_err=.false.)
       if (atm_co2_nf_ind == 0) then
          call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
               &/ 'atmopsheric CO2 from the flux coupler ' /&
               &/ 'and it is not present')
       endif
    endif

    if (lflux_gas_co2 .and. atm_co2_iopt == atm_co2_iopt_drv_diag) then
       call named_field_get_index('ATM_CO2_DIAG', atm_co2_nf_ind, &
            exit_on_err=.false.)
       if (atm_co2_nf_ind == 0) then
          call exit_POP(sigAbort, 'ecosys_init: ecosys module requires ' /&
               &/ 'atmopsheric CO2 from the flux coupler ' /&
               &/ 'and it is not present')
       endif
    endif

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_init_sflux

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_init_interior_restore
  ! !INTERFACE:

  subroutine ecosys_init_interior_restore(saved_state, ecosys_restore)

    ! !DESCRIPTION:
    !  Initialize interior restoring computations for ecosys tracer module.
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_interface_types, only : marbl_saved_state_type

    type(marbl_saved_state_type), intent(in) :: saved_state
    type(ecosys_restore_type), intent(inout) :: ecosys_restore
    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: &
         k,                   & ! index for looping over levels
         i, j,                 & ! index for looping over horiz. dims.
         iblock                 ! index for looping over blocks

    real (r8) :: &
         subsurf_fesed          ! sum of subsurface fesed values

    !-----------------------------------------------------------------------
    !  initialize restoring timescale (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%initialize_restoring_timescale(nml_filename, nml_in, zt)

    !-----------------------------------------------------------------------
    !  load restoring fields (if required)
    !-----------------------------------------------------------------------

    call ecosys_restore%read_restoring_fields(saved_state%land_mask)

    !-----------------------------------------------------------------------
    !  load fesedflux
    !  add subsurface positives to 1 level shallower, to accomodate overflow pop-ups
    !-----------------------------------------------------------------------

    allocate(FESEDFLUX(nx_block, ny_block, km, max_blocks_clinic))

    call read_field(fesedflux_input%file_fmt, &
         fesedflux_input%filename, &
         fesedflux_input%file_varname, &
         FESEDFLUX)

    do iblock=1, nblocks_clinic
       do j=1, ny_block
          do i=1, nx_block
             if (KMT(i, j, iblock) > 0 .and. KMT(i, j, iblock) < km) then
                subsurf_fesed = c0
                do k=KMT(i, j, iblock)+1, km
                   subsurf_fesed = subsurf_fesed + FESEDFLUX(i, j, k, iblock)
                enddo
                FESEDFLUX(i, j, KMT(i, j, iblock), iblock) = FESEDFLUX(i, j, KMT(i, j, iblock), iblock) + subsurf_fesed
             endif
          enddo
       enddo

       do k = 1, km
          where (.not. saved_state%land_mask(:, :, iblock) .or. k > KMT(:, :, iblock)) &
               FESEDFLUX(:, :, k, iblock) = c0
          FESEDFLUX(:, :, k, iblock) = FESEDFLUX(:, :, k, iblock) * fesedflux_input%scale_factor
       enddo
    end do

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_init_interior_restore

  !***********************************************************************
  !BOP
  ! !IROUTINE: ecosys_set_sflux
  ! !INTERFACE:

  subroutine ecosys_set_sflux( &
       saved_state, &
       marbl_surface_share, &
       SHF_QSW_RAW, SHF_QSW, &
       U10_SQR, IFRAC, PRESS, SST, SSS, &
       SURF_VALS_OLD, SURF_VALS_CUR, STF_MODULE, &
       lexport_shared_vars)

    use co2calc, only : co2calc_row
    use schmidt_number, only : SCHMIDT_CO2
    use marbl_oxygen, only : schmidt_o2
    use marbl_oxygen, only : o2sat
    use marbl_share_mod, only : ecosys_surface_share_type
    use marbl_interface_types, only : marbl_saved_state_type
    use strdata_interface_mod, only : POP_strdata_create
    use strdata_interface_mod, only : POP_strdata_advance

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:

    type(marbl_saved_state_type), intent(inout) :: saved_state

    real (r8), dimension(nx_block, ny_block, max_blocks_clinic), intent(in) :: &
         SHF_QSW_RAW,  &! penetrative solar heat flux, from coupler (degC*cm/s)
         SHF_QSW,      &! SHF_QSW used by physics, may have diurnal cylce imposed (degC*cm/s)
         U10_SQR,      &! 10m wind speed squared (cm/s)**2
         IFRAC,        &! sea ice fraction (non-dimensional)
         PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
         SST,          &! sea surface temperature (C)
         SSS            ! sea surface salinity (psu)

    real (r8), dimension(:, :, :, :), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

    logical (log_kind), intent(in) :: &
         lexport_shared_vars ! flag to save shared_vars or not

    ! !INPUT/OUTPUT PARAMETERS:

    real (r8), dimension(:, :, :, :), &
         intent(inout) :: STF_MODULE

    type(ecosys_surface_share_type), intent(inout) :: marbl_surface_share
    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_set_sflux'

    logical (log_kind) :: first_call = .true.

    type (block) :: &
         this_block      ! block info for the current block

    integer (int_kind) :: &
         i, j, iblock, n, & ! loop indices
         auto_ind,        & ! autotroph functional group index
         errorCode          ! errorCode from HaloUpdate call

    real (r8), dimension(nx_block, ny_block, max_blocks_clinic) :: &
         IFRAC_USED,   & ! used ice fraction (non-dimensional)
         XKW_USED,     & ! portion of piston velocity (cm/s)
         AP_USED,      & ! used atm pressure (converted from dyne/cm**2 to atm)
         IRON_FLUX_IN    ! iron flux

    real (r8), dimension(nx_block, ny_block, max_blocks_clinic) :: &
         SHR_STREAM_WORK

    real (r8), dimension(nx_block, ny_block) :: &
         XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
         SCHMIDT_USED, & ! used Schmidt number
         PV,           & ! piston velocity (cm/s)
         O2SAT_1atm,   & ! O2 saturation @ 1 atm (mmol/m^3)
         O2SAT_USED,   & ! used O2 saturation (mmol/m^3)
         XCO2,         & ! atmospheric co2 conc. (dry-air, 1 atm)
         XCO2_ALT_CO2, & ! atmospheric alternative CO2 (dry-air, 1 atm)
         FLUX,         & ! tracer flux (nmol/cm^2/s)
         FLUX_ALT_CO2    ! tracer flux alternative CO2 (nmol/cm^2/s)


    real (r8), dimension(nx_block) :: &
         PHLO,         & ! lower bound for ph in solver
         PHHI,         & ! upper bound for ph in solver
         PH_NEW,       & ! computed PH from solver
         DIC_ROW,      & ! row of DIC values for solver
         ALK_ROW,      & ! row of ALK values for solver
         PO4_ROW,      & ! row of PO4 values for solver
         SiO3_ROW,     & ! row of SiO3 values for solver
         CO2STAR_ROW,  & ! CO2STAR from solver
         DCO2STAR_ROW, & ! DCO2STAR from solver
         pCO2SURF_ROW, & ! pCO2SURF from solver
         DpCO2_ROW,    & ! DpCO2 from solver
         CO3_ROW

    character (char_len) :: &
         tracer_data_label          ! label for what is being updated

    character (char_len), dimension(1) :: &
         tracer_data_names          ! short names for input data fields

    integer (int_kind), dimension(1) :: &
         tracer_bndy_loc,         & ! location and field type for ghost
         tracer_bndy_type           !    cell updates

    real (r8), dimension(nx_block, ny_block) :: &
         WORK1 ! temporaries for averages

    real (r8) :: scalar_temp

    associate(                                                              &
         PV_SURF_fields       => marbl_surface_share%PV_SURF_fields,       & ! IN/OUT
         DIC_SURF_fields      => marbl_surface_share%DIC_SURF_fields,      & ! IN/OUT
         CO2STAR_SURF_fields  => marbl_surface_share%CO2STAR_SURF_fields,  & ! IN/OUT
         DCO2STAR_SURF_fields => marbl_surface_share%DCO2STAR_SURF_fields, & ! IN/OUT
         CO3_SURF_fields      => marbl_surface_share%CO3_SURF_fields,      & ! IN/OUT
         dic_riv_flux_fields  => marbl_surface_share%dic_riv_flux_fields,  & ! IN/OUT
         doc_riv_flux_fields  => marbl_surface_share%doc_riv_flux_fields   & ! IN/OUT
         )

    !-----------------------------------------------------------------------

    call timer_start(ecosys_sflux_timer)

    !-----------------------------------------------------------------------

    if (check_time_flag(comp_surf_avg_flag))  &
         call comp_surf_avg(SURF_VALS_OLD, SURF_VALS_CUR, &
         ecosys_tracer_cnt, vflux_flag, surf_avg)

    !-----------------------------------------------------------------------
    !  fluxes initially set to 0
    !  set Chl field for short-wave absorption
    !  store incoming shortwave in PAR_out field, converting to W/m^2
    !-----------------------------------------------------------------------

    scalar_temp = f_qsw_par / hflux_factor

    !$OMP PARALLEL DO PRIVATE(iblock, WORK1, auto_ind, n)
    do iblock = 1, nblocks_clinic
       STF_MODULE(:, :, :, iblock) = c0

       WORK1 = c0
       do auto_ind = 1, autotroph_cnt
          n = autotrophs(auto_ind)%Chl_ind
          WORK1 = WORK1 + max(c0, p5*(SURF_VALS_OLD(:, :, n, iblock) + &
               SURF_VALS_CUR(:, :, n, iblock)))
       end do
       call named_field_set(totChl_surf_nf_ind, iblock, WORK1)

       if (ecosys_qsw_distrb_const) then
          saved_state%PAR_out(:, :, iblock) = SHF_QSW_RAW(:, :, iblock)
       else
          saved_state%PAR_out(:, :, iblock) = SHF_QSW(:, :, iblock)
       endif

       where (saved_state%land_mask(:, :, iblock))
          saved_state%PAR_out(:, :, iblock) = max(c0, scalar_temp * saved_state%PAR_out(:, :, iblock))
       elsewhere
          saved_state%PAR_out(:, :, iblock) = c0
       end where
    enddo
    !$OMP END PARALLEL DO

    !-----------------------------------------------------------------------
    !  Interpolate gas flux forcing data if necessary
    !-----------------------------------------------------------------------

    if ((lflux_gas_o2 .or. lflux_gas_co2) .and. &
         gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(fice_file%data_time,      &
               fice_file%data_time_min_loc,  fice_file%interp_type,    &
               fice_file%data_next,          fice_file%data_update,    &
               fice_file%data_type,          fice_file%data_inc,       &
               fice_file%DATA(:, :, :, :, 1:12), fice_file%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:, :, :, :, 1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_USED = INTERP_WORK(:, :, :, 1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(xkw_file%data_time,      &
               xkw_file%data_time_min_loc,  xkw_file%interp_type,    &
               xkw_file%data_next,          xkw_file%data_update,    &
               xkw_file%data_type,          xkw_file%data_inc,       &
               xkw_file%DATA(:, :, :, :, 1:12), xkw_file%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            xkw_file%DATA(:, :, :, :, 1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_USED = INTERP_WORK(:, :, :, 1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(ap_file%data_time,    &
               ap_file%data_time_min_loc,  ap_file%interp_type,    &
               ap_file%data_next,          ap_file%data_update,    &
               ap_file%data_type,          ap_file%data_inc,       &
               ap_file%DATA(:, :, :, :, 1:12), ap_file%data_renorm,    &
               tracer_data_label,          tracer_data_names,      &
               tracer_bndy_loc,            tracer_bndy_type,       &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:, :, :, :, 1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = INTERP_WORK(:, :, :, 1)

    endif

    !-----------------------------------------------------------------------
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

    if (lflux_gas_o2 .or. lflux_gas_co2) then

       !$OMP PARALLEL DO PRIVATE(iblock, j, XKW_ICE, SCHMIDT_USED, PV, O2SAT_USED, &
       !$OMP                     O2SAT_1atm, FLUX, FLUX_ALT_CO2, XCO2, XCO2_ALT_CO2, &
       !$OMP                     PHLO, PHHI, DIC_ROW, ALK_ROW, &
       !$OMP                     PO4_ROW, SiO3_ROW, PH_NEW, CO2STAR_ROW, &
       !$OMP                     DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, &
       !$OMP                     CO3_ROW)

       do iblock = 1, nblocks_clinic

          !-----------------------------------------------------------------------
          !  Apply OCMIP ice fraction mask when input is from a file.
          !-----------------------------------------------------------------------

          if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then
             where (IFRAC_USED(:, :, iblock) < 0.2000_r8) &
                  IFRAC_USED(:, :, iblock) = 0.2000_r8
             where (IFRAC_USED(:, :, iblock) > 0.9999_r8) &
                  IFRAC_USED(:, :, iblock) = 0.9999_r8
          endif

          if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv) then
             IFRAC_USED(:, :, iblock) = IFRAC(:, :, iblock)
             where (IFRAC_USED(:, :, iblock) < c0) IFRAC_USED(:, :, iblock) = c0
             where (IFRAC_USED(:, :, iblock) > c1) IFRAC_USED(:, :, iblock) = c1
             XKW_USED(:, :, iblock) = xkw_coeff * U10_SQR(:, :, iblock)
             AP_USED(:, :, iblock) = PRESS(:, :, iblock)
          endif

          !-----------------------------------------------------------------------
          !  assume PRESS is in cgs units (dyne/cm**2) since that is what is
          !    required for pressure forcing in barotropic
          !  want units to be atmospheres
          !  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
          !  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
          !-----------------------------------------------------------------------

          AP_USED(:, :, iblock) = PRESS(:, :, iblock) / 101.325e+4_r8

          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_IFRAC, iblock) = IFRAC_USED(:, :, iblock)
          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_XKW, iblock) = XKW_USED(:, :, iblock)
          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_ATM_PRESS, iblock) = AP_USED(:, :, iblock)

          !-----------------------------------------------------------------------
          !  Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
          !-----------------------------------------------------------------------

          XKW_ICE = (c1 - IFRAC_USED(:, :, iblock)) * XKW_USED(:, :, iblock)

          !-----------------------------------------------------------------------
          !  compute O2 flux
          !-----------------------------------------------------------------------

          if (lflux_gas_o2) then
             SCHMIDT_USED = SCHMIDT_O2(nx_block, ny_block, SST(:, :, iblock), saved_state%land_mask(:, :, iblock))

             O2SAT_1atm = O2SAT(nx_block, ny_block, SST(:, :, iblock), SSS(:, :, iblock),   &
                  saved_state%land_mask(:, :, iblock))

             where (saved_state%land_mask(:, :, iblock))
                PV = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED)
                O2SAT_USED = AP_USED(:, :, iblock) * O2SAT_1atm
                FLUX = PV * (O2SAT_USED - p5*(SURF_VALS_OLD(:, :, o2_ind, iblock) + &
                     SURF_VALS_CUR(:, :, o2_ind, iblock)))
                STF_MODULE(:, :, o2_ind, iblock) = STF_MODULE(:, :, o2_ind, iblock) + FLUX
             elsewhere
                O2SAT_USED = c0
             end where

             ECO_SFLUX_TAVG(:, :, buf_ind_PV_O2, iblock) = PV
             ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_O2, iblock) = SCHMIDT_USED
             ECO_SFLUX_TAVG(:, :, buf_ind_O2SAT, iblock) = O2SAT_USED

          endif  ! lflux_gas_o2

          !-----------------------------------------------------------------------
          !  compute CO2 flux, computing disequilibrium one row at a time
          !-----------------------------------------------------------------------

          if (lflux_gas_co2) then

             SCHMIDT_USED = SCHMIDT_CO2(SST(:, :, iblock), saved_state%land_mask(:, :, iblock))

             where (saved_state%land_mask(:, :, iblock))
                PV = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED)
             elsewhere
                PV = c0
             end where

             ! Save surface field of PV for use in other modules
             if (lexport_shared_vars) PV_SURF_fields(:, :, iblock) = PV

             !-----------------------------------------------------------------------
             !  Set XCO2
             !-----------------------------------------------------------------------

             select case (atm_co2_iopt)
             case (atm_co2_iopt_const)
                XCO2 = atm_co2_const
             case (atm_co2_iopt_drv_prog, atm_co2_iopt_drv_diag)
                call named_field_get(atm_co2_nf_ind, iblock, XCO2)
             end select

             select case (atm_alt_co2_iopt)
             case (atm_co2_iopt_const)
                XCO2_ALT_CO2 = atm_alt_co2_const
             end select

             do j = 1, ny_block
                where (PH_PREV(:, j, iblock) /= c0)
                   PHLO = PH_PREV(:, j, iblock) - del_ph
                   PHHI = PH_PREV(:, j, iblock) + del_ph
                elsewhere
                   PHLO = phlo_surf_init
                   PHHI = phhi_surf_init
                end where

                DIC_ROW = p5*(SURF_VALS_OLD(:, j, dic_ind, iblock) + &
                     SURF_VALS_CUR(:, j, dic_ind, iblock))
                ALK_ROW = p5*(SURF_VALS_OLD(:, j, alk_ind, iblock) + &
                     SURF_VALS_CUR(:, j, alk_ind, iblock))
                PO4_ROW = p5*(SURF_VALS_OLD(:, j, po4_ind, iblock) + &
                     SURF_VALS_CUR(:, j, po4_ind, iblock))
                SiO3_ROW = p5*(SURF_VALS_OLD(:, j, sio3_ind, iblock) + &
                     SURF_VALS_CUR(:, j, sio3_ind, iblock))

                call co2calc_row(iblock, j, saved_state%land_mask(:, j, iblock), &
                     locmip_k1_k2_bug_fix, .true., &
                     SST(:, j, iblock), SSS(:, j, iblock), &
                     DIC_ROW, ALK_ROW, PO4_ROW, SiO3_ROW, &
                     PHLO, PHHI, PH_NEW, XCO2(:, j), &
                     AP_USED(:, j, iblock), CO2STAR_ROW, &
                     DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, &
                     CO3_ROW)

                PH_PREV(:, j, iblock) = PH_NEW

                FLUX(:, j) = PV(:, j) * DCO2STAR_ROW

                ECO_SFLUX_TAVG(:, j, buf_ind_CO2STAR, iblock)  = CO2STAR_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_DCO2STAR, iblock) = DCO2STAR_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_pCO2SURF, iblock) = pCO2SURF_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_DpCO2, iblock)    = DpCO2_ROW

                !-------------------------------------------------------------------
                !  The following variables need to be shared with other modules, and
                !  are now defined in marbl_share as targets.
                !-------------------------------------------------------------------
                if (lexport_shared_vars) then
                   DIC_SURF_fields(:, j, iblock)      = DIC_ROW
                   CO2STAR_SURF_fields(:, j, iblock)  = CO2STAR_ROW
                   DCO2STAR_SURF_fields(:, j, iblock) = DCO2STAR_ROW
                   CO3_SURF_fields(:, j, iblock)      = CO3_ROW
                endif

                where (PH_PREV_ALT_CO2(:, j, iblock) /= c0)
                   PHLO = PH_PREV_ALT_CO2(:, j, iblock) - del_ph
                   PHHI = PH_PREV_ALT_CO2(:, j, iblock) + del_ph
                elsewhere
                   PHLO = phlo_surf_init
                   PHHI = phhi_surf_init
                end where

                DIC_ROW = p5*(SURF_VALS_OLD(:, j, dic_alt_co2_ind, iblock) + &
                     SURF_VALS_CUR(:, j, dic_alt_co2_ind, iblock))

                call co2calc_row(iblock, j, saved_state%land_mask(:, j, iblock), &
                     locmip_k1_k2_bug_fix, .false., &
                     SST(:, j, iblock), SSS(:, j, iblock), &
                     DIC_ROW, ALK_ROW, PO4_ROW, SiO3_ROW, &
                     PHLO, PHHI, PH_NEW, XCO2_ALT_CO2(:, j), &
                     AP_USED(:, j, iblock), CO2STAR_ROW, &
                     DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, &
                     CO3_ROW)

                PH_PREV_ALT_CO2(:, j, iblock) = PH_NEW

                FLUX_ALT_CO2(:, j) = PV(:, j) * DCO2STAR_ROW

                ECO_SFLUX_TAVG(:, j, buf_ind_CO2STAR_ALT_CO2, iblock) = CO2STAR_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_DCO2STAR_ALT_CO2, iblock) = DCO2STAR_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_pCO2SURF_ALT_CO2, iblock) = pCO2SURF_ROW
                ECO_SFLUX_TAVG(:, j, buf_ind_DpCO2_ALT_CO2, iblock) = DpCO2_ROW

             end do

             !-----------------------------------------------------------------------
             !  set air-sea co2 gas flux named field, converting units from
             !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
             !-----------------------------------------------------------------------

             call named_field_set(sflux_co2_nf_ind, iblock, 44.0e-8_r8 * FLUX)

             STF_MODULE(:, :, dic_ind, iblock) = STF_MODULE(:, :, dic_ind, iblock) + FLUX

             ECO_SFLUX_TAVG(:, :, buf_ind_PV_CO2, iblock) = PV
             ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_CO2, iblock) = SCHMIDT_USED
             ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX, iblock) = FLUX
             ECO_SFLUX_TAVG(:, :, buf_ind_PH, iblock) = PH_PREV(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_ATM_CO2, iblock) = XCO2

             STF_MODULE(:, :, dic_alt_co2_ind, iblock) = STF_MODULE(:, :, dic_alt_co2_ind, iblock) + FLUX_ALT_CO2

             ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX_ALT_CO2, iblock) = FLUX_ALT_CO2
             ECO_SFLUX_TAVG(:, :, buf_ind_PH_ALT_CO2, iblock) = PH_PREV_ALT_CO2(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_ATM_ALT_CO2, iblock) = XCO2_ALT_CO2

          endif  !  lflux_gas_co2

       enddo
       !$OMP END PARALLEL DO

    endif  ! lflux_gas_o2 .or. lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    if (iron_flux%has_data) then
       if (thour00 >= iron_flux%data_update) then
          tracer_data_names = iron_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Iron Flux'
          call update_forcing_data(iron_flux%data_time,    &
               iron_flux%data_time_min_loc,  iron_flux%interp_type,    &
               iron_flux%data_next,          iron_flux%data_update,    &
               iron_flux%data_type,          iron_flux%data_inc,       &
               iron_flux%DATA(:, :, :, :, 1:12), iron_flux%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               iron_flux%filename,           iron_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            iron_flux%DATA(:, :, :, :, 1:12), &
            iron_flux%data_time,         iron_flux%interp_type, &
            iron_flux%data_time_min_loc, iron_flux%interp_freq, &
            iron_flux%interp_inc,        iron_flux%interp_next, &
            iron_flux%interp_last,       0)
       if (liron_patch .and. imonth == iron_patch_month) then
          IRON_FLUX_IN = INTERP_WORK(:, :, :, 1) + IRON_PATCH_FLUX
       else
          IRON_FLUX_IN = INTERP_WORK(:, :, :, 1)
       endif
    else
       IRON_FLUX_IN = c0
    endif

    IRON_FLUX_IN = IRON_FLUX_IN * parm_Fe_bioavail

    STF_MODULE(:, :, fe_ind, :) = STF_MODULE(:, :, fe_ind, :) + IRON_FLUX_IN
    ECO_SFLUX_TAVG(:, :, buf_ind_IRON_FLUX, :) = IRON_FLUX_IN

    if (dust_flux%has_data) then
       if (thour00 >= dust_flux%data_update) then
          tracer_data_names = dust_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Dust Flux'
          call update_forcing_data(dust_flux%data_time,    &
               dust_flux%data_time_min_loc,  dust_flux%interp_type,    &
               dust_flux%data_next,          dust_flux%data_update,    &
               dust_flux%data_type,          dust_flux%data_inc,       &
               dust_flux%DATA(:, :, :, :, 1:12), dust_flux%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               dust_flux%filename,           dust_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            dust_flux%DATA(:, :, :, :, 1:12),    &
            dust_flux%data_time,         dust_flux%interp_type, &
            dust_flux%data_time_min_loc, dust_flux%interp_freq, &
            dust_flux%interp_inc,        dust_flux%interp_next, &
            dust_flux%interp_last,       0)
       saved_state%dust_FLUX_IN = INTERP_WORK(:, :, :, 1)

       !-----------------------------------------------------------------------
       !  Reduce surface dust flux due to assumed instant surface dissolution
       !  Can't use parm_fe_bioavail when using solFe input files
       !-----------------------------------------------------------------------

       !     dust_FLUX_IN = dust_FLUX_IN * (c1 - parm_Fe_bioavail)
       saved_state%dust_FLUX_IN = saved_state%dust_FLUX_IN * 0.98_r8
    else
       saved_state%dust_FLUX_IN = c0
    endif

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------

    if (nox_flux_monthly%has_data) then
       if (thour00 >= nox_flux_monthly%data_update) then
          tracer_data_names = nox_flux_monthly%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'NOx Flux'
          call update_forcing_data(nox_flux_monthly%data_time,    &
               nox_flux_monthly%data_time_min_loc,  nox_flux_monthly%interp_type, &
               nox_flux_monthly%data_next,          nox_flux_monthly%data_update, &
               nox_flux_monthly%data_type,          nox_flux_monthly%data_inc,    &
               nox_flux_monthly%DATA(:, :, :, :, 1:12), nox_flux_monthly%data_renorm, &
               tracer_data_label,                   tracer_data_names,            &
               tracer_bndy_loc,                     tracer_bndy_type,             &
               nox_flux_monthly%filename,           nox_flux_monthly%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            nox_flux_monthly%DATA(:, :, :, :, 1:12), &
            nox_flux_monthly%data_time,         nox_flux_monthly%interp_type, &
            nox_flux_monthly%data_time_min_loc, nox_flux_monthly%interp_freq, &
            nox_flux_monthly%interp_inc,        nox_flux_monthly%interp_next, &
            nox_flux_monthly%interp_last,       0)
       STF_MODULE(:, :, no3_ind, :) = STF_MODULE(:, :, no3_ind, :) + INTERP_WORK(:, :, :, 1)
       ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, :) = INTERP_WORK(:, :, :, 1)
    endif

    if (nhy_flux_monthly%has_data) then
       if (thour00 >= nhy_flux_monthly%data_update) then
          tracer_data_names = nhy_flux_monthly%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'NHy Flux'
          call update_forcing_data(nhy_flux_monthly%data_time,    &
               nhy_flux_monthly%data_time_min_loc,  nhy_flux_monthly%interp_type, &
               nhy_flux_monthly%data_next,          nhy_flux_monthly%data_update, &
               nhy_flux_monthly%data_type,          nhy_flux_monthly%data_inc,    &
               nhy_flux_monthly%DATA(:, :, :, :, 1:12), nhy_flux_monthly%data_renorm, &
               tracer_data_label,                   tracer_data_names,            &
               tracer_bndy_loc,                     tracer_bndy_type,             &
               nhy_flux_monthly%filename,           nhy_flux_monthly%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            nhy_flux_monthly%DATA(:, :, :, :, 1:12), &
            nhy_flux_monthly%data_time,         nhy_flux_monthly%interp_type, &
            nhy_flux_monthly%data_time_min_loc, nhy_flux_monthly%interp_freq, &
            nhy_flux_monthly%interp_inc,        nhy_flux_monthly%interp_next, &
            nhy_flux_monthly%interp_last,       0)
       STF_MODULE(:, :, nh4_ind, :) = STF_MODULE(:, :, nh4_ind, :) + INTERP_WORK(:, :, :, 1)
    endif

#ifdef CCSMCOUPLED
    if (trim(ndep_data_type) == 'shr_stream') then
       if (first_call) then

          ndep_inputlist%field_name = 'ndep data'
          ndep_inputlist%short_name = 'ndep'
          ndep_inputlist%year_first = ndep_shr_stream_year_first
          ndep_inputlist%year_last = ndep_shr_stream_year_last
          ndep_inputlist%year_align = ndep_shr_stream_year_align
          ndep_inputlist%file_name = ndep_shr_stream_file
          ndep_inputlist%field_list = ' '
          do n = 1, ndep_shr_stream_var_cnt
             if (n == ndep_shr_stream_no_ind) &
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) /&
                  &/ 'NOy_deposition'
             if (n == ndep_shr_stream_nh_ind) &
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) /&
                  &/ 'NHx_deposition'
             if (n < ndep_shr_stream_var_cnt) &
                  ndep_inputlist%field_list = trim(ndep_inputlist%field_list) /&
                  &/ ':'
          end do

          call POP_strdata_create(ndep_inputlist)
          first_call = .false.
       endif

       ndep_inputlist%date = iyear*10000 + imonth*100 + iday
       ndep_inputlist%time = isecond + 60 * (iminute + 60 * ihour)

       call timer_start(ecosys_shr_strdata_advance_timer)
       call POP_strdata_advance(ndep_inputlist)
       call timer_stop(ecosys_shr_strdata_advance_timer)

       !
       ! process NO3 flux, store results in SHR_STREAM_WORK array
       ! instead of directly into STF_MODULE
       ! to avoid argument copies in HaloUpdate calls
       !
       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock), iblock)
          do j=this_block%jb, this_block%je
             do i=this_block%ib, this_block%ie
                n = n + 1
                SHR_STREAM_WORK(i, j, iblock) = &
                     ndep_inputlist%sdat%avs(1)%rAttr(ndep_shr_stream_no_ind, n)
             enddo
          enddo
       enddo

       call POP_HaloUpdate(SHR_STREAM_WORK, POP_haloClinic, &
            POP_gridHorzLocCenter,          &
            POP_fieldKindScalar, errorCode, &
            fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call exit_POP(sigAbort, subname /&
               &/ ': error updating halo for Ndep fields')
       endif

       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1, nblocks_clinic
          where (saved_state%land_mask(:, :, iblock))
             STF_MODULE(:, :, no3_ind, iblock) = STF_MODULE(:, :, no3_ind, iblock) &
                  + ndep_shr_stream_scale_factor * SHR_STREAM_WORK(:, :, iblock)
          endwhere
          ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, iblock) = &
               ndep_shr_stream_scale_factor * SHR_STREAM_WORK(:, :, iblock)
       enddo
       !$OMP END PARALLEL DO

       !
       ! process NH4 flux, store results in SHR_STREAM_WORK array
       ! instead of directly into STF_MODULE
       ! to avoid argument copies in HaloUpdate calls
       !
       n = 0
       do iblock = 1, nblocks_clinic
          this_block = get_block(blocks_clinic(iblock), iblock)
          do j=this_block%jb, this_block%je
             do i=this_block%ib, this_block%ie
                n = n + 1
                SHR_STREAM_WORK(i, j, iblock) = &
                     ndep_inputlist%sdat%avs(1)%rAttr(ndep_shr_stream_nh_ind, n)
             enddo
          enddo
       enddo

       call POP_HaloUpdate(SHR_STREAM_WORK, POP_haloClinic, &
            POP_gridHorzLocCenter,          &
            POP_fieldKindScalar, errorCode, &
            fillValue = 0.0_POP_r8)
       if (errorCode /= POP_Success) then
          call exit_POP(sigAbort, subname /&
               &/ ': error updating halo for Ndep fields')
       endif

       !$OMP PARALLEL DO PRIVATE(iblock)
       do iblock = 1, nblocks_clinic
          where (saved_state%land_mask(:, :, iblock))
             STF_MODULE(:, :, nh4_ind, iblock) = STF_MODULE(:, :, nh4_ind, iblock) &
                  + ndep_shr_stream_scale_factor * SHR_STREAM_WORK(:, :, iblock)
          endwhere
       enddo
       !$OMP END PARALLEL DO

    endif
#endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------

    if (din_riv_flux%has_data) then
       if (thour00 >= din_riv_flux%data_update) then
          tracer_data_names = din_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DIN River Flux'
          call update_forcing_data(din_riv_flux%data_time,    &
               din_riv_flux%data_time_min_loc,  din_riv_flux%interp_type,    &
               din_riv_flux%data_next,          din_riv_flux%data_update,    &
               din_riv_flux%data_type,          din_riv_flux%data_inc,       &
               din_riv_flux%DATA(:, :, :, :, 1:12), din_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               din_riv_flux%filename,           din_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            din_riv_flux%DATA(:, :, :, :, 1:12), &
            din_riv_flux%data_time,         din_riv_flux%interp_type, &
            din_riv_flux%data_time_min_loc, din_riv_flux%interp_freq, &
            din_riv_flux%interp_inc,        din_riv_flux%interp_next, &
            din_riv_flux%interp_last,       0)
       STF_MODULE(:, :, no3_ind, :) = STF_MODULE(:, :, no3_ind, :) + INTERP_WORK(:, :, :, 1)
       ECO_SFLUX_TAVG(:, :, buf_ind_DIN_RIV_FLUX, :) = INTERP_WORK(:, :, :, 1)
    endif

    if (dip_riv_flux%has_data) then
       if (thour00 >= dip_riv_flux%data_update) then
          tracer_data_names = dip_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DIP River Flux'
          call update_forcing_data(dip_riv_flux%data_time,    &
               dip_riv_flux%data_time_min_loc,  dip_riv_flux%interp_type,    &
               dip_riv_flux%data_next,          dip_riv_flux%data_update,    &
               dip_riv_flux%data_type,          dip_riv_flux%data_inc,       &
               dip_riv_flux%DATA(:, :, :, :, 1:12), dip_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dip_riv_flux%filename,           dip_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dip_riv_flux%DATA(:, :, :, :, 1:12), &
            dip_riv_flux%data_time,         dip_riv_flux%interp_type, &
            dip_riv_flux%data_time_min_loc, dip_riv_flux%interp_freq, &
            dip_riv_flux%interp_inc,        dip_riv_flux%interp_next, &
            dip_riv_flux%interp_last,       0)
       STF_MODULE(:, :, po4_ind, :) = STF_MODULE(:, :, po4_ind, :) + INTERP_WORK(:, :, :, 1)
    endif

    if (don_riv_flux%has_data) then
       if (thour00 >= don_riv_flux%data_update) then
          tracer_data_names = don_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DON River Flux'
          call update_forcing_data(don_riv_flux%data_time,    &
               don_riv_flux%data_time_min_loc,  don_riv_flux%interp_type,    &
               don_riv_flux%data_next,          don_riv_flux%data_update,    &
               don_riv_flux%data_type,          don_riv_flux%data_inc,       &
               don_riv_flux%DATA(:, :, :, :, 1:12), don_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               don_riv_flux%filename,           don_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            don_riv_flux%DATA(:, :, :, :, 1:12), &
            don_riv_flux%data_time,         don_riv_flux%interp_type, &
            don_riv_flux%data_time_min_loc, don_riv_flux%interp_freq, &
            don_riv_flux%interp_inc,        don_riv_flux%interp_next, &
            don_riv_flux%interp_last,       0)
       STF_MODULE(:, :, don_ind, :) = STF_MODULE(:, :, don_ind, :) + (INTERP_WORK(:, :, :, 1) * 0.9_r8)
       STF_MODULE(:, :, donr_ind, :) = STF_MODULE(:, :, donr_ind, :) + (INTERP_WORK(:, :, :, 1) * 0.1_r8)
    endif

    if (dop_riv_flux%has_data) then
       if (thour00 >= dop_riv_flux%data_update) then
          tracer_data_names = dop_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(dop_riv_flux%data_time,    &
               dop_riv_flux%data_time_min_loc,  dop_riv_flux%interp_type,    &
               dop_riv_flux%data_next,          dop_riv_flux%data_update,    &
               dop_riv_flux%data_type,          dop_riv_flux%data_inc,       &
               dop_riv_flux%DATA(:, :, :, :, 1:12), dop_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dop_riv_flux%filename,           dop_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dop_riv_flux%DATA(:, :, :, :, 1:12), &
            dop_riv_flux%data_time,         dop_riv_flux%interp_type, &
            dop_riv_flux%data_time_min_loc, dop_riv_flux%interp_freq, &
            dop_riv_flux%interp_inc,        dop_riv_flux%interp_next, &
            dop_riv_flux%interp_last,       0)
       STF_MODULE(:, :, dop_ind, :) = STF_MODULE(:, :, dop_ind, :) + (INTERP_WORK(:, :, :, 1) * 0.975_r8)
       STF_MODULE(:, :, dopr_ind, :) = STF_MODULE(:, :, dopr_ind, :) + (INTERP_WORK(:, :, :, 1) * 0.025_r8)
    endif

    if (dsi_riv_flux%has_data) then
       if (thour00 >= dsi_riv_flux%data_update) then
          tracer_data_names = dsi_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(dsi_riv_flux%data_time,    &
               dsi_riv_flux%data_time_min_loc,  dsi_riv_flux%interp_type,    &
               dsi_riv_flux%data_next,          dsi_riv_flux%data_update,    &
               dsi_riv_flux%data_type,          dsi_riv_flux%data_inc,       &
               dsi_riv_flux%DATA(:, :, :, :, 1:12), dsi_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dsi_riv_flux%filename,           dsi_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dsi_riv_flux%DATA(:, :, :, :, 1:12), &
            dsi_riv_flux%data_time,         dsi_riv_flux%interp_type, &
            dsi_riv_flux%data_time_min_loc, dsi_riv_flux%interp_freq, &
            dsi_riv_flux%interp_inc,        dsi_riv_flux%interp_next, &
            dsi_riv_flux%interp_last,       0)
       STF_MODULE(:, :, sio3_ind, :) = STF_MODULE(:, :, sio3_ind, :) + INTERP_WORK(:, :, :, 1)
    endif

    if (dfe_riv_flux%has_data) then
       if (thour00 >= dfe_riv_flux%data_update) then
          tracer_data_names = dfe_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(dfe_riv_flux%data_time,    &
               dfe_riv_flux%data_time_min_loc,  dfe_riv_flux%interp_type,    &
               dfe_riv_flux%data_next,          dfe_riv_flux%data_update,    &
               dfe_riv_flux%data_type,          dfe_riv_flux%data_inc,       &
               dfe_riv_flux%DATA(:, :, :, :, 1:12), dfe_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dfe_riv_flux%filename,           dfe_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dfe_riv_flux%DATA(:, :, :, :, 1:12), &
            dfe_riv_flux%data_time,         dfe_riv_flux%interp_type, &
            dfe_riv_flux%data_time_min_loc, dfe_riv_flux%interp_freq, &
            dfe_riv_flux%interp_inc,        dfe_riv_flux%interp_next, &
            dfe_riv_flux%interp_last,       0)
       STF_MODULE(:, :, fe_ind, :) = STF_MODULE(:, :, fe_ind, :) + INTERP_WORK(:, :, :, 1)
       ECO_SFLUX_TAVG(:, :, buf_ind_DFE_RIV_FLUX, :) = INTERP_WORK(:, :, :, 1)
    endif

    if (dic_riv_flux%has_data) then
       if (thour00 >= dic_riv_flux%data_update) then
          tracer_data_names = dic_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(dic_riv_flux%data_time,    &
               dic_riv_flux%data_time_min_loc,  dic_riv_flux%interp_type,    &
               dic_riv_flux%data_next,          dic_riv_flux%data_update,    &
               dic_riv_flux%data_type,          dic_riv_flux%data_inc,       &
               dic_riv_flux%DATA(:, :, :, :, 1:12), dic_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               dic_riv_flux%filename,           dic_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            dic_riv_flux%DATA(:, :, :, :, 1:12), &
            dic_riv_flux%data_time,         dic_riv_flux%interp_type, &
            dic_riv_flux%data_time_min_loc, dic_riv_flux%interp_freq, &
            dic_riv_flux%interp_inc,        dic_riv_flux%interp_next, &
            dic_riv_flux%interp_last,       0)
       STF_MODULE(:, :, dic_ind, :) = STF_MODULE(:, :, dic_ind, :) + INTERP_WORK(:, :, :, 1)
       STF_MODULE(:, :, dic_alt_co2_ind, :) = STF_MODULE(:, :, dic_alt_co2_ind, :) + INTERP_WORK(:, :, :, 1)
       ECO_SFLUX_TAVG(:, :, buf_ind_DIC_RIV_FLUX, :) = INTERP_WORK(:, :, :, 1)
       if (lexport_shared_vars) dic_riv_flux_fields=INTERP_WORK(:, :, :, 1)
    endif

    if (alk_riv_flux%has_data) then
       if (thour00 >= alk_riv_flux%data_update) then
          tracer_data_names = alk_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'DOP River Flux'
          call update_forcing_data(alk_riv_flux%data_time,    &
               alk_riv_flux%data_time_min_loc,  alk_riv_flux%interp_type,    &
               alk_riv_flux%data_next,          alk_riv_flux%data_update,    &
               alk_riv_flux%data_type,          alk_riv_flux%data_inc,       &
               alk_riv_flux%DATA(:, :, :, :, 1:12), alk_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               alk_riv_flux%filename,           alk_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            alk_riv_flux%DATA(:, :, :, :, 1:12), &
            alk_riv_flux%data_time,         alk_riv_flux%interp_type, &
            alk_riv_flux%data_time_min_loc, alk_riv_flux%interp_freq, &
            alk_riv_flux%interp_inc,        alk_riv_flux%interp_next, &
            alk_riv_flux%interp_last,       0)
       STF_MODULE(:, :, alk_ind, :) = STF_MODULE(:, :, alk_ind, :) + INTERP_WORK(:, :, :, 1)
       ECO_SFLUX_TAVG(:, :, buf_ind_ALK_RIV_FLUX, :) = INTERP_WORK(:, :, :, 1)
    endif

    if (doc_riv_flux%has_data) then
       if (thour00 >= doc_riv_flux%data_update) then
          tracer_data_names = doc_riv_flux%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'PP River Flux'
          call update_forcing_data(doc_riv_flux%data_time,    &
               doc_riv_flux%data_time_min_loc,  doc_riv_flux%interp_type,    &
               doc_riv_flux%data_next,          doc_riv_flux%data_update,    &
               doc_riv_flux%data_type,          doc_riv_flux%data_inc,       &
               doc_riv_flux%DATA(:, :, :, :, 1:12), doc_riv_flux%data_renorm,    &
               tracer_data_label,           tracer_data_names,       &
               tracer_bndy_loc,             tracer_bndy_type,        &
               doc_riv_flux%filename,           doc_riv_flux%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK,     &
            doc_riv_flux%DATA(:, :, :, :, 1:12), &
            doc_riv_flux%data_time,         doc_riv_flux%interp_type, &
            doc_riv_flux%data_time_min_loc, doc_riv_flux%interp_freq, &
            doc_riv_flux%interp_inc,        doc_riv_flux%interp_next, &
            doc_riv_flux%interp_last,       0)
       STF_MODULE(:, :, doc_ind, :) = STF_MODULE(:, :, doc_ind, :) + INTERP_WORK(:, :, :, 1)
       if (lexport_shared_vars) doc_riv_flux_fields=INTERP_WORK(:, :, :, 1)
    endif

    !-----------------------------------------------------------------------
    !  Apply NO & NH fluxes to alkalinity
    !-----------------------------------------------------------------------

    STF_MODULE(:, :, alk_ind, :) = STF_MODULE(:, :, alk_ind, :) &
         + STF_MODULE(:, :, nh4_ind, :) - STF_MODULE(:, :, no3_ind, :)

    !-----------------------------------------------------------------------

    call timer_stop(ecosys_sflux_timer)

    end associate
    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_set_sflux

  !*****************************************************************************
  !BOP
  ! !IROUTINE: ecosys_tavg_forcing
  ! !INTERFACE:

  subroutine ecosys_tavg_forcing(saved_state, STF_MODULE, FLUX_DIAGS)

    ! !DESCRIPTION:
    !  Accumulate non-standard forcing related tavg variables.
    !
    ! !REVISION HISTORY:
    !  same as module

    use marbl_interface_types, only : marbl_saved_state_type

    type(marbl_saved_state_type), intent(in) :: saved_state

    ! !INPUT PARAMETERS:

    real (r8), dimension(:, :, :, :), &
         intent(in) :: STF_MODULE

    ! !OUTPUT PARAMETERS:
    real (r8), dimension(nx_block, ny_block, forcing_diag_cnt, nblocks_clinic), &
         intent(out) :: &
         FLUX_DIAGS                ! Computed diagnostics for surface fluxes

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! accumulate surface flux related fields in the order in which they are declared
    !
    !  multiply IRON, DUST fluxes by mpercm (.01) to convert from model
    !    units (cm/s)(mmol/m^3) to mmol/s/m^2
    !-----------------------------------------------------------------------

    FLUX_DIAGS(:, :, ECOSYS_IFRAC_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_IFRAC, :)
    FLUX_DIAGS(:, :, ECOSYS_XKW_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_XKW, :)
    FLUX_DIAGS(:, :, ECOSYS_ATM_PRESS_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_ATM_PRESS, :)
    FLUX_DIAGS(:, :, PV_O2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_PV_O2, :)
    FLUX_DIAGS(:, :, SCHMIDT_O2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_O2, :)
    FLUX_DIAGS(:, :, O2SAT_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_O2SAT, :)
    FLUX_DIAGS(:, :, O2_GAS_FLUX_diag_ind, :) = STF_MODULE(:, :, o2_ind, :)
    FLUX_DIAGS(:, :, CO2STAR_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR, :)
    FLUX_DIAGS(:, :, DCO2STAR_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR, :)
    FLUX_DIAGS(:, :, pCO2SURF_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF, :)
    FLUX_DIAGS(:, :, DpCO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2, :)
    FLUX_DIAGS(:, :, PV_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_PV_CO2, :)
    FLUX_DIAGS(:, :, SCHMIDT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_CO2, :)
    FLUX_DIAGS(:, :, DIC_GAS_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX, :)
    FLUX_DIAGS(:, :, PH_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_PH, :)
    FLUX_DIAGS(:, :, ATM_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ATM_CO2, :)
    FLUX_DIAGS(:, :, CO2STAR_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR_ALT_CO2, :)
    FLUX_DIAGS(:, :, DCO2STAR_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR_ALT_CO2, :)
    FLUX_DIAGS(:, :, pCO2SURF_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF_ALT_CO2, :)
    FLUX_DIAGS(:, :, DpCO2_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2_ALT_CO2, :)
    FLUX_DIAGS(:, :, DIC_GAS_FLUX_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX_ALT_CO2, :)
    FLUX_DIAGS(:, :, PH_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_PH_ALT_CO2, :)
    FLUX_DIAGS(:, :, ATM_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ATM_ALT_CO2, :)
    FLUX_DIAGS(:, :, IRON_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_IRON_FLUX, :)*mpercm
    FLUX_DIAGS(:, :, DUST_FLUX_diag_ind, :) = saved_state%dust_FLUX_IN(:, :, :)*mpercm
    FLUX_DIAGS(:, :, NOx_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, :)
    FLUX_DIAGS(:, :, NHy_FLUX_diag_ind, :) = STF_MODULE(:, :, nh4_ind, :)
    FLUX_DIAGS(:, :, DIN_RIV_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DIN_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DIP_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, po4_ind, :)
    FLUX_DIAGS(:, :, DON_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, don_ind, :)
    FLUX_DIAGS(:, :, DONr_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, donr_ind, :)
    FLUX_DIAGS(:, :, DOP_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, dop_ind, :)
    FLUX_DIAGS(:, :, DOPr_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, dopr_ind, :)
    FLUX_DIAGS(:, :, DSI_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, sio3_ind, :)
    FLUX_DIAGS(:, :, DFE_RIV_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DFE_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DIC_RIV_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_RIV_FLUX, :)
    FLUX_DIAGS(:, :, ALK_RIV_FLUX_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_ALK_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DOC_RIV_FLUX_diag_ind, :) = STF_MODULE(:, :, doc_ind, :)

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_tavg_forcing

  !*****************************************************************************
  !BOP
  ! !IROUTINE: ecosys_write_restart
  ! !INTERFACE:

  subroutine ecosys_write_restart(saved_state, restart_file, action)

    ! !DESCRIPTION:
    !  write auxiliary fields & scalars to restart files
    !
    ! !REVISION HISTORY:
    !  same as module
    use marbl_interface_types, only : marbl_saved_state_type


    ! !INPUT PARAMETERS:
    type(marbl_saved_state_type), intent(in) :: saved_state
    character(*), intent(in) :: action

    ! !INPUT/OUTPUT PARAMETERS:

    type (datafile), intent (inout)  :: restart_file

    !EOP
    !BOC
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character (char_len) :: &
         short_name  ! tracer name temporaries

    type (io_dim) :: &
         i_dim, j_dim, & ! dimension descriptors
         k_dim           ! dimension descriptor for vertical levels

    integer (int_kind) :: n

    type (io_field_desc), save :: PH_SURF, PH_SURF_ALT_CO2, &
         PH_3D_ALT_CO2,  PH_3D

    !-----------------------------------------------------------------------

    if (trim(action) == 'add_attrib_file') then
       short_name = char_blank
       do n=1, ecosys_tracer_cnt
          if (vflux_flag(n)) then
             short_name = 'surf_avg_' /&
                  &/ ind_name_table(n)%name
             call add_attrib_file(restart_file, trim(short_name), surf_avg(n))
          endif
       end do
    endif

    if (trim(action) == 'define') then
       i_dim = construct_io_dim('i', nx_global)
       j_dim = construct_io_dim('j', ny_global)
       k_dim = construct_io_dim('k', km)

       PH_SURF = construct_io_field('PH_SURF', i_dim, j_dim,     &
            long_name='surface pH at current time',      &
            units='pH', grid_loc='2110',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d2d_array = PH_PREV(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_SURF)

       PH_SURF_ALT_CO2 = construct_io_field('PH_SURF_ALT_CO2', i_dim, j_dim, &
            long_name='surface pH, alternate CO2, at current time', &
            units='pH', grid_loc='2110',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d2d_array = PH_PREV_ALT_CO2(:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_SURF_ALT_CO2)

       PH_3D_ALT_CO2 = construct_io_field('PH_3D_ALT_CO2', i_dim, j_dim, k_dim, &
            long_name='3D pH, alternate CO2, at current time', &
            units='pH', grid_loc='3111',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d3d_array = saved_state%PH_PREV_ALT_CO2_3D(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_3D_ALT_CO2)

       PH_3D = construct_io_field('PH_3D', i_dim, j_dim, k_dim, &
            long_name='3D pH at current time', &
            units='pH', grid_loc='3111',            &
            field_loc = field_loc_center,                &
            field_type = field_type_scalar,              &
            d3d_array = saved_state%PH_PREV_3D(:,:,:,1:nblocks_clinic))
       call data_set (restart_file, 'define', PH_3D)

    endif

    if (trim(action) == 'write') then
       call data_set (restart_file, 'write', PH_SURF)
       call data_set (restart_file, 'write', PH_SURF_ALT_CO2)
       call data_set (restart_file, 'write', PH_3D)
       call data_set (restart_file, 'write', PH_3D_ALT_CO2)
    endif

    !-----------------------------------------------------------------------
    !EOC

  end subroutine ecosys_write_restart

  !*****************************************************************************

  !*****************************************************************************
  !BOP
  ! !IROUTINE: ecosys_tracer_ref_val
  ! !INTERFACE:

  function ecosys_tracer_ref_val(ind)

    ! !DESCRIPTION:
    !  return reference value for tracers using virtual fluxes
    !
    ! !REVISION HISTORY:
    !  same as module

    ! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: ind

    ! !OUTPUT PARAMETERS:

    real (r8) :: ecosys_tracer_ref_val

    !EOP
    !BOC
    !-----------------------------------------------------------------------

    if (vflux_flag(ind)) then
       ecosys_tracer_ref_val = surf_avg(ind)
    else
       ecosys_tracer_ref_val = c0
    endif

    !-----------------------------------------------------------------------
    !EOC

  end function ecosys_tracer_ref_val

  !***********************************************************************

  subroutine ecosys_init_forcing_monthly_every_ts()
    !
    ! initialize
    !
    implicit none

    call init_forcing_monthly_every_ts(dust_flux)
    call init_forcing_monthly_every_ts(iron_flux)
    call init_forcing_monthly_every_ts(fice_file)
    call init_forcing_monthly_every_ts(xkw_file)
    call init_forcing_monthly_every_ts(ap_file)
    call init_forcing_monthly_every_ts(nox_flux_monthly)
    call init_forcing_monthly_every_ts(nhy_flux_monthly)
    call init_forcing_monthly_every_ts(din_riv_flux)
    call init_forcing_monthly_every_ts(dip_riv_flux)
    call init_forcing_monthly_every_ts(don_riv_flux)
    call init_forcing_monthly_every_ts(dop_riv_flux)
    call init_forcing_monthly_every_ts(dsi_riv_flux)
    call init_forcing_monthly_every_ts(dfe_riv_flux)
    call init_forcing_monthly_every_ts(dic_riv_flux)
    call init_forcing_monthly_every_ts(alk_riv_flux)
    call init_forcing_monthly_every_ts(doc_riv_flux)

  end subroutine ecosys_init_forcing_monthly_every_ts

  !***********************************************************************

  subroutine ecosys_init_non_autotroph_tracer_metadata(tracer_d_module, &
       non_living_biomass_ecosys_tracer_cnt)
    !-----------------------------------------------------------------------
    !  initialize non-autotroph tracer_d values and accumulate
    !  non_living_biomass_ecosys_tracer_cnt
    !-----------------------------------------------------------------------
    implicit none

    type (tracer_field), dimension(:), intent(inout) :: &
         tracer_d_module   ! descriptors for each tracer

    integer (int_kind), intent(out) :: &
         non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers

    integer(int_kind) :: n

    non_living_biomass_ecosys_tracer_cnt = 0

    tracer_d_module(po4_ind)%short_name='PO4'
    tracer_d_module(po4_ind)%long_name='Dissolved Inorganic Phosphate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(no3_ind)%short_name='NO3'
    tracer_d_module(no3_ind)%long_name='Dissolved Inorganic Nitrate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(sio3_ind)%short_name='SiO3'
    tracer_d_module(sio3_ind)%long_name='Dissolved Inorganic Silicate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(nh4_ind)%short_name='NH4'
    tracer_d_module(nh4_ind)%long_name='Dissolved Ammonia'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(fe_ind)%short_name='Fe'
    tracer_d_module(fe_ind)%long_name='Dissolved Inorganic Iron'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(o2_ind)%short_name='O2'
    tracer_d_module(o2_ind)%long_name='Dissolved Oxygen'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(dic_ind)%short_name='DIC'
    tracer_d_module(dic_ind)%long_name='Dissolved Inorganic Carbon'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(dic_alt_co2_ind)%short_name='DIC_ALT_CO2'
    tracer_d_module(dic_alt_co2_ind)%long_name='Dissolved Inorganic Carbon, Alternative CO2'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(alk_ind)%short_name='ALK'
    tracer_d_module(alk_ind)%long_name='Alkalinity'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(doc_ind)%short_name='DOC'
    tracer_d_module(doc_ind)%long_name='Dissolved Organic Carbon'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(don_ind)%short_name='DON'
    tracer_d_module(don_ind)%long_name='Dissolved Organic Nitrogen'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(dofe_ind)%short_name='DOFe'
    tracer_d_module(dofe_ind)%long_name='Dissolved Organic Iron'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(dop_ind)%short_name='DOP'
    tracer_d_module(dop_ind)%long_name='Dissolved Organic Phosphorus'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(dopr_ind)%short_name='DOPr'
    tracer_d_module(dopr_ind)%long_name='Refractory DOP'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    tracer_d_module(donr_ind)%short_name='DONr'
    tracer_d_module(donr_ind)%long_name='Refractory DON'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    do n = 1, non_living_biomass_ecosys_tracer_cnt
       if (n == alk_ind) then
          tracer_d_module(n)%units      = 'meq/m^3'
          tracer_d_module(n)%tend_units = 'meq/m^3/s'
          tracer_d_module(n)%flux_units = 'meq/m^3 cm/s'
       else
          tracer_d_module(n)%units      = 'mmol/m^3'
          tracer_d_module(n)%tend_units = 'mmol/m^3/s'
          tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       endif
    end do

  end subroutine ecosys_init_non_autotroph_tracer_metadata

  !***********************************************************************

  subroutine check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt)

    implicit none

    integer (int_kind), intent(in) :: &
         non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers

    integer (int_kind) :: &
         n,               &
         auto_ind,        & ! autotroph functional group index
         zoo_ind            ! zooplankton functional group index

    character(*), parameter :: subname = 'ecosys_mod:check_ecosys_tracer_count_consistency'

    !-----------------------------------------------------------------------
    !  confirm that ecosys_tracer_cnt is consistent with autotroph declarations
    !-----------------------------------------------------------------------
    n = non_living_biomass_ecosys_tracer_cnt
    ! Do we really need a loop here, or would simple addition work?!
    do zoo_ind = 1, zooplankton_cnt
       n = n + 1 ! C
    end do

    do auto_ind = 1, autotroph_cnt
       n = n + 3 ! Chl, C, Fe tracers
       if (autotrophs(auto_ind)%kSiO3 > c0) n = n + 1 ! Si tracer
       if (autotrophs(auto_ind)%imp_calcifier .or. &
            autotrophs(auto_ind)%exp_calcifier) n = n + 1 ! CaCO3 tracer
    end do

    if (ecosys_tracer_cnt /= n) then
       call document(subname, 'actual ecosys_tracer_cnt', ecosys_tracer_cnt)
       call document(subname, 'computed ecosys_tracer_cnt', n)
       call exit_POP(sigAbort, 'inconsistency between actual ecosys_tracer_cnt and computed ecosys_tracer_cnt')
    endif

  end subroutine check_ecosys_tracer_count_consistency

  !***********************************************************************

  subroutine initialize_zooplankton_tracer_metadata(tracer_d_module, &
       non_living_biomass_ecosys_tracer_cnt, n)
    !-----------------------------------------------------------------------
    !  initialize zooplankton tracer_d values and tracer indices
    !-----------------------------------------------------------------------
    type (tracer_field), dimension(:), intent(inout) :: &
         tracer_d_module   ! descriptors for each tracer

    integer (int_kind), intent(in) :: &
         non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers

    integer (int_kind), intent(inout) :: n

    integer (int_kind) :: &
         zoo_ind            ! zooplankton functional group index

    n = non_living_biomass_ecosys_tracer_cnt + 1

    do zoo_ind = 1, zooplankton_cnt
       tracer_d_module(n)%short_name = trim(zooplankton(zoo_ind)%sname) // 'C'
       tracer_d_module(n)%long_name  = trim(zooplankton(zoo_ind)%lname) // ' Carbon'
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       zooplankton(zoo_ind)%C_ind = n
       n = n + 1
    end do

    if (my_task == master_task) THEN
       write (stdout, *) '----- zooplankton tracer indices -----'
       do zoo_ind = 1, zooplankton_cnt
          write (stdout, *) 'C_ind(', trim(zooplankton(zoo_ind)%sname), ') = ', zooplankton(zoo_ind)%C_ind
       end do
       write (stdout, *) '------------------------------------'
    endif

  end subroutine initialize_zooplankton_tracer_metadata

  !***********************************************************************

  subroutine initialize_autotroph_tracer_metadata(tracer_d_module, n)
    !-----------------------------------------------------------------------
    !  initialize autotroph tracer_d values and tracer indices
    !-----------------------------------------------------------------------
    type (tracer_field), dimension(:), intent(inout) :: &
         tracer_d_module   ! descriptors for each tracer

    integer(int_kind), intent(inout) :: n

    integer (int_kind) :: &
         auto_ind            ! zooplankton functional group index

    do auto_ind = 1, autotroph_cnt
       tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Chl'
       tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Chlorophyll'
       tracer_d_module(n)%units      = 'mg/m^3'
       tracer_d_module(n)%tend_units = 'mg/m^3/s'
       tracer_d_module(n)%flux_units = 'mg/m^3 cm/s'
       autotrophs(auto_ind)%Chl_ind = n
       n = n + 1

       tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'C'
       tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon'
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%C_ind = n
       n = n + 1

       tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Fe'
       tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Iron'
       tracer_d_module(n)%units      = 'mmol/m^3'
       tracer_d_module(n)%tend_units = 'mmol/m^3/s'
       tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%Fe_ind = n
       n = n + 1

       if (autotrophs(auto_ind)%kSiO3 > c0) then
          tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Si'
          tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Silicon'
          tracer_d_module(n)%units      = 'mmol/m^3'
          tracer_d_module(n)%tend_units = 'mmol/m^3/s'
          tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%Si_ind = n
          n = n + 1
       else
          autotrophs(auto_ind)%Si_ind = 0
       endif

       if (autotrophs(auto_ind)%imp_calcifier .or. &
            autotrophs(auto_ind)%exp_calcifier) then
          tracer_d_module(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'CaCO3'
          tracer_d_module(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' CaCO3'
          tracer_d_module(n)%units      = 'mmol/m^3'
          tracer_d_module(n)%tend_units = 'mmol/m^3/s'
          tracer_d_module(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%CaCO3_ind = n
          n = n + 1
       else
          autotrophs(auto_ind)%CaCO3_ind = 0
       endif
    end do

    if (my_task == master_task) THEN
       write (stdout, *) '----- autotroph tracer indices -----'
       do auto_ind = 1, autotroph_cnt
          write (stdout, *) 'Chl_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Chl_ind
          write (stdout, *) 'C_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%C_ind
          write (stdout, *) 'Fe_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Fe_ind
          write (stdout, *) 'Si_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%Si_ind
          write (stdout, *) 'CaCO3_ind(', trim(autotrophs(auto_ind)%sname), ') = ', autotrophs(auto_ind)%CaCO3_ind
       end do
       write (stdout, *) '------------------------------------'
    endif

  end subroutine initialize_autotroph_tracer_metadata

  !***********************************************************************

  subroutine setup_local_tracers(k, column_land_mask, column_kmt, tracer_module, tracer_local)
    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------

    integer(int_kind), intent(in) :: k
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt

    real (r8), dimension(ecosys_tracer_cnt), intent(in) :: tracer_module ! tracer values

    real (r8), dimension(ecosys_tracer_cnt), intent(out) :: tracer_local ! local copies of model tracer concentrations

    integer (int_kind) :: &
         n ! tracer index


    ! FIXME(bja, 2015-06) only need to loop over
    ! non-living-biomass-ecosys-tracer-cnt. Does it actually need to be
    ! a loop?
    do n = 1, ecosys_tracer_cnt
       tracer_local(n) = max(c0, tracer_module(n))
       if (.not. column_land_mask .or. k > column_kmt) then
          tracer_local(n) = c0
       end if
    end do

  end subroutine setup_local_tracers

  !***********************************************************************

  subroutine setup_local_zooplankton(k, column_land_mask, column_kmt, tracer_module, zoo_cnt, zoo, zooplankton_local)
    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------
    use marbl_share_mod, only : zooplankton_type

    integer (int_kind), intent(in) :: k
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt
    real (r8), dimension(:), intent(in) :: tracer_module ! tracer values

    ! FIXME(bja, 2015-07) shortening zooplankton to zoo to avoid
    ! a namespace collision with the global imported into the
    ! module. Need to fix after global is removed.
    integer(int_kind), intent(in) :: zoo_cnt
    type(zooplankton_type), dimension(zoo_cnt), intent(in) :: zoo

    type(zooplankton_local_type), dimension(zoo_cnt), intent(out) :: zooplankton_local


    integer (int_kind) :: &
         zoo_ind, n ! tracer index


    do zoo_ind = 1, zoo_cnt
       n = zoo(zoo_ind)%C_ind
       zooplankton_local(zoo_ind)%C = max(c0, tracer_module(n))
       if (.not. column_land_mask .or. k > column_kmt) then
          zooplankton_local(zoo_ind)%C = c0
       end if
    end do

  end subroutine setup_local_zooplankton


  !***********************************************************************

  subroutine setup_local_autotrophs(k, column_land_mask, column_kmt, tracer_module, auto_cnt, auto_meta, autotroph_loc)
    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------
    use marbl_share_mod, only : autotroph_type

    use constants, only : T0_Kelvin

    ! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: k
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt
    real (r8), dimension(:), intent(in) :: tracer_module ! tracer values

    ! FIXME(bja, 2015-07) autotroph --> auto are horrible names, but
    ! can't use full name until it is removed from the global
    ! namespace!
    integer(int_kind), intent(in) :: auto_cnt ! autotroph_cnt
    type(autotroph_type), dimension(auto_cnt), intent(in) :: auto_meta ! autotrophs

    type(autotroph_local_type), dimension(auto_cnt), intent(out) :: autotroph_loc

    integer (int_kind) :: &
         auto_ind, tracer_ind ! tracer index


    do auto_ind = 1, autotroph_cnt

       tracer_ind = autotrophs(auto_ind)%Chl_ind
       autotroph_loc(auto_ind)%Chl = max(c0, tracer_module(tracer_ind))

       tracer_ind = autotrophs(auto_ind)%C_ind
       autotroph_loc(auto_ind)%C = max(c0, tracer_module(tracer_ind))

       tracer_ind = autotrophs(auto_ind)%Fe_ind
       autotroph_loc(auto_ind)%Fe = max(c0, tracer_module(tracer_ind))

       tracer_ind = autotrophs(auto_ind)%Si_ind
       if (tracer_ind > 0) then
          autotroph_loc(auto_ind)%Si = max(c0, tracer_module(tracer_ind))
       endif

       tracer_ind = autotrophs(auto_ind)%CaCO3_ind
       if (tracer_ind > 0) then
          autotroph_loc(auto_ind)%CaCO3 = max(c0, tracer_module(tracer_ind))
       endif

       if (.not. column_land_mask .or. k > column_kmt) then
          autotroph_loc(auto_ind)%Chl = c0
          autotroph_loc(auto_ind)%C = c0
          autotroph_loc(auto_ind)%Fe = c0
          autotroph_loc(auto_ind)%Si = c0
          autotroph_loc(auto_ind)%CaCO3 = c0
       end if
    end do

  end subroutine setup_local_autotrophs

  !***********************************************************************

  subroutine autotroph_consistency_check(auto_cnt, auto_meta, &
       autotroph_local)
    !-----------------------------------------------------------------------
    !  If any phyto box are zero, set others to zeros.
    !-----------------------------------------------------------------------
    use marbl_share_mod, only : autotroph_type

    ! FIXME(bja, 2015-07) autotroph --> auto are horrible names, but
    ! can't use full name until it is removed from the global
    ! namespace!
    integer(int_kind), intent(in) :: auto_cnt ! autotroph_cnt
    type(autotroph_type), dimension(auto_cnt), intent(in) :: auto_meta ! autotrophs

    type(autotroph_local_type), dimension(auto_cnt), intent(inout) :: autotroph_local


    integer(int_kind) :: auto_ind
    logical (log_kind) :: zero_mask

    do auto_ind = 1, auto_cnt
       zero_mask = (autotroph_local(auto_ind)%Chl == c0 .or. &
            autotroph_local(auto_ind)%C == c0 .or. &
            autotroph_local(auto_ind)%Fe == c0)
       if (auto_meta(auto_ind)%Si_ind > 0) then
          zero_mask = zero_mask .or. autotroph_local(auto_ind)%Si == c0
       end if
       if (zero_mask) then
          autotroph_local(auto_ind)%Chl = c0
          autotroph_local(auto_ind)%C = c0
          autotroph_local(auto_ind)%Fe = c0
       end if
       if (auto_meta(auto_ind)%Si_ind > 0) then
          if (zero_mask) then
             autotroph_local(auto_ind)%Si = c0
          end if
       end if
       if (auto_meta(auto_ind)%CaCO3_ind > 0) then
          if (zero_mask) then
             autotroph_local(auto_ind)%CaCO3 = c0
          end if
       end if
    end do

  end subroutine autotroph_consistency_check

  !***********************************************************************

  subroutine compute_autotroph_elemental_ratios(k, auto_cnt, auto_meta, autotroph_local, tracer_local, autotroph_secondary_species)

    use marbl_share_mod, only : autotroph_type

    use marbl_parms, only : epsC
    use marbl_parms, only : gQsi_0
    use marbl_parms, only : gQsi_max
    use marbl_parms, only : gQsi_min

    integer (int_kind), intent(in) :: k
    integer (int_kind), intent(in) :: auto_cnt
    type(autotroph_type), dimension(auto_cnt), intent(in) :: auto_meta ! autotrophs

    type(autotroph_local_type), intent(in) :: autotroph_local(auto_cnt)
    real (r8), intent(in) :: tracer_local(ecosys_tracer_cnt) ! local copies of model tracer concentrations
    type(autotroph_secondary_species_type), intent(inout) :: autotroph_secondary_species(auto_cnt)

    real :: &
         cks,            & ! constant used in Fe quota modification
         cksi              ! constant used in Si quota modification

    integer(int_kind) :: auto_ind

    associate(                                          &
         Fe_loc   => tracer_local(fe_ind),       &
         SiO3_loc => tracer_local(sio3_ind),     &
         auto_C   => autotroph_local(:)%C,      &
         auto_Chl => autotroph_local(:)%Chl,      &
         auto_Fe  => autotroph_local(:)%Fe,      &
         auto_Si  => autotroph_local(:)%Si,      &
         auto_CaCO3  => autotroph_local(:)%CaCO3, &

         thetaC   => autotroph_secondary_species(:)%thetaC , & ! local Chl/C ratio (mg Chl/mmol C)
         QCaCO3   => autotroph_secondary_species(:)%QCaCO3 , & ! CaCO3/C ratio (mmol CaCO3/mmol C)
         Qfe      => autotroph_secondary_species(:)%Qfe,     & ! init fe/C ratio (mmolFe/mmolC)
         gQfe     => autotroph_secondary_species(:)%gQfe,    & ! fe/C for growth
         Qsi      => autotroph_secondary_species(:)%Qsi,     & ! initial Si/C ratio (mmol Si/mmol C)
         gQsi     => autotroph_secondary_species(:)%gQsi     & ! diatom Si/C ratio for growth (new biomass)
         )

    !-----------------------------------------------------------------------
    !  set local variables, with incoming ratios
    !-----------------------------------------------------------------------

    do auto_ind = 1, autotroph_cnt
       thetaC(auto_ind) = auto_Chl(auto_ind) / (auto_C(auto_ind) + epsC)
       Qfe(auto_ind) = auto_Fe(auto_ind) / (auto_C(auto_ind) + epsC)
       if (autotrophs(auto_ind)%Si_ind > 0) then
          Qsi(auto_ind) = min(auto_Si(auto_ind) / (auto_C(auto_ind) + epsC), gQsi_max)
       endif
    end do

    !-----------------------------------------------------------------------
    !  DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)
    !  Modify these initial ratios under low ambient iron conditions
    !  Modify the initial si/C ratio under low ambient Si conditions
    !-----------------------------------------------------------------------

    cks = 9._r8
    cksi = 5._r8

    do auto_ind = 1, autotroph_cnt
       gQfe(auto_ind) = autotrophs(auto_ind)%gQfe_0
       if (Fe_loc < cks * autotrophs(auto_ind)%kFe) then
          gQfe(auto_ind) = &
               max((gQfe(auto_ind) * Fe_loc / (cks * autotrophs(auto_ind)%kFe)), &
               autotrophs(auto_ind)%gQfe_min)
       end if

       if (autotrophs(auto_ind)%Si_ind > 0) then
          gQsi(auto_ind) = gQsi_0
          if ((Fe_loc < cksi * autotrophs(auto_ind)%kFe) .and. &
               (Fe_loc > c0) .and. &
               (SiO3_loc > (cksi * autotrophs(auto_ind)%kSiO3))) then
             gQsi(auto_ind) = min((gQsi(auto_ind) * cksi * autotrophs(auto_ind)%kFe / Fe_loc), gQsi_max)
          end if

          if (Fe_loc == c0) then
             gQsi(auto_ind) = gQsi_max
          end if

          if (SiO3_loc < (cksi * autotrophs(auto_ind)%kSiO3)) then
             gQsi(auto_ind) = max((gQsi(auto_ind) * SiO3_loc / (cksi * autotrophs(auto_ind)%kSiO3)), &
                  gQsi_min)
          end if
       endif

       !-----------------------------------------------------------------------
       !  QCaCO3 is the percentage of sp organic matter which is associated
       !  with coccolithophores
       !-----------------------------------------------------------------------

       if (autotrophs(auto_ind)%CaCO3_ind > 0) then
          QCaCO3(auto_ind) = auto_CaCO3(auto_ind) / (auto_C(auto_ind) + epsC)
          if (QCaCO3(auto_ind) > QCaCO3_max) then
             QCaCO3(auto_ind) = QCaCO3_max
          end if
       end if
    end do
    end associate
  end subroutine compute_autotroph_elemental_ratios

  !***********************************************************************

  subroutine compute_photosynthetically_available_radiation(k, auto_cnt, &
       autotroph_local, &
       column_land_mask, column_kmt, column_dzt, column_dz, &
       PAR_out, PAR)
    !-----------------------------------------------------------------------
    !  compute PAR related quantities
    !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
    !  0.45   fraction of incoming SW -> PAR (non-dim)
    !-----------------------------------------------------------------------

    integer(int_kind), intent(in) :: k, auto_cnt
    type(autotroph_local_type), intent(in) :: autotroph_local(auto_cnt)
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt
    real (r8), intent(in) :: &
         column_dzt, &
         column_dz

    real (r8), intent(inout) :: PAR_out

    type(photosynthetically_available_radiation_type), intent(out) :: PAR


    real (r8) :: WORK1

    PAR%in = PAR_out
    if (.not. column_land_mask .or. k > column_kmt) then
       PAR%in = c0
    end if

    ! TODO(bja, 2015-07) move calculation outside and just pass in this
    ! work array as autotroph_Chl instead of passing in all of
    ! autotroph_loc?
    WORK1 = max(sum(autotroph_local(:)%Chl, dim=1), 0.02_r8)
    if (WORK1 < 0.13224_r8) then
       PAR%KPARdz = 0.000919_r8*(WORK1**0.3536_r8)
    else
       PAR%KPARdz = 0.001131_r8*(WORK1**0.4562_r8)
    end if

    if (partial_bottom_cells) then
       PAR%KPARdz = PAR%KPARdz * column_dzt
    else
       PAR%KPARdz = PAR%KPARdz * column_dz
    endif

    PAR_out = PAR%in * exp(-PAR%KPARdz)
    PAR%avg = PAR%in * (c1 - exp(-PAR%KPARdz)) / PAR%KPARdz


  end subroutine compute_photosynthetically_available_radiation

  !***********************************************************************


  subroutine compute_carbonate_chemistry(dkm, bid, &
       column_land_mask, column_kmt, &
       temperature, salinity, &
       tracer_local, carbonate, &
       ph_prev_3d, ph_prev_alt_co2_3d, &
       zsat_calcite, zsat_aragonite, &
       co3_calcite_anom, co3_aragonite_anom)

    use co2calc_column, only : comp_co3terms
    use co2calc_column, only : comp_co3_sat_vals
    use co2calc_column, only : thermodynamic_coefficients_type
    USE state_mod,      ONLY : ref_pressure

    integer (int_kind), intent(in) :: dkm
    integer (int_kind), intent(in) :: bid
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt

    real (r8), intent(in) :: &
         temperature(dkm),   &! old potential temperature (C)
         salinity(dkm)        ! current salinity (msu)

    real (r8), intent(in) :: tracer_local(ecosys_tracer_cnt,dkm) ! local copies of model tracer concentrations

    type(carbonate_type), intent(out) :: carbonate(dkm)

    real(r8), intent(inout) :: ph_prev_3d(dkm)
    real(r8), intent(inout) :: ph_prev_alt_co2_3d(dkm)

    real(r8), intent(inout) :: zsat_calcite(km) ! Calcite Saturation Depth
    real(r8), intent(inout) :: zsat_aragonite(km) ! Aragonite Saturation Depth

    real(r8), intent(inout) :: co3_calcite_anom(km) ! CO3 concentration above calcite saturation at k-1
    real(r8), intent(inout) :: co3_aragonite_anom(km) ! CO3 concentration above aragonite saturation at k-1

    integer :: k
    logical(log_kind), dimension(dkm) :: mask
    logical(log_kind), dimension(dkm) :: pressure_correct
    real(r8), dimension(dkm) :: ph_lower_bound
    real(r8), dimension(dkm) :: ph_upper_bound
    real(r8), dimension(dkm) :: press_bar           ! pressure at level (bars)
    type(thermodynamic_coefficients_type), dimension(dkm) :: co3_coeffs

    call timer_start(ecosys_comp_CO3terms_timer, block_id=bid)

    associate( &
         DIC_loc => tracer_local(dic_ind,:), &
         DIC_ALT_CO2_loc => tracer_local(dic_alt_co2_ind,:),  &
         ALK_loc => tracer_local(alk_ind,:), &
         PO4_loc => tracer_local(po4_ind,:), &
         SiO3_loc => tracer_local(sio3_ind,:), &
         CO3 => carbonate%CO3, &
         HCO3 => carbonate%HCO3, &
         H2CO3 => carbonate%H2CO3, &
         pH => carbonate%pH, &
         CO3_sat_calcite => carbonate%CO3_sat_calcite, &
         CO3_sat_aragonite => carbonate%CO3_sat_aragonite, &
         CO3_ALT_CO2 => carbonate%CO3_ALT_CO2, &
         HCO3_ALT_CO2 => carbonate%HCO3_ALT_CO2, &
         H2CO3_ALT_CO2 => carbonate%H2CO3_ALT_CO2, &
         pH_ALT_CO2 => carbonate%pH_ALT_CO2 &
         )


    pressure_correct = .TRUE.
    pressure_correct(1) = .FALSE.
    do k=1,dkm

      mask(k) = column_land_mask .and. k <= column_kmt
      press_bar(k) = ref_pressure(k)

       ! -------------------
       if (ph_prev_3d(k)  /= c0) then
          ph_lower_bound(k) = ph_prev_3d(k) - del_ph
          ph_upper_bound(k) = ph_prev_3d(k) + del_ph
       else
          ph_lower_bound(k) = phlo_3d_init
          ph_upper_bound(k) = phhi_3d_init
       end if

    enddo

    call comp_CO3terms(dkm, mask, pressure_correct, .true., co3_coeffs, temperature, &
                       salinity, press_bar, DIC_loc, ALK_loc, PO4_loc, SiO3_loc, &
                       ph_lower_bound, ph_upper_bound, pH, H2CO3, HCO3, CO3)
       
    do k=1,dkm

       ph_prev_3d(k) = pH(k)
       
       ! -------------------
       if (ph_prev_alt_co2_3d(k) /= c0) then
          ph_lower_bound(k) = ph_prev_alt_co2_3d(k) - del_ph
          ph_upper_bound(k) = ph_prev_alt_co2_3d(k) + del_ph
       else
          ph_lower_bound(k) = phlo_3d_init
          ph_upper_bound(k) = phhi_3d_init
       end if

    enddo

    call comp_CO3terms(dkm, mask, pressure_correct, .false., co3_coeffs, temperature,    &
                       salinity, press_bar, DIC_ALT_CO2_loc, ALK_loc, PO4_loc, SiO3_loc, &
                       ph_lower_bound, ph_upper_bound, pH_ALT_CO2, H2CO3_ALT_CO2,        &
                       HCO3_ALT_CO2, CO3_ALT_CO2)
       
    ph_prev_alt_co2_3d = pH_ALT_CO2

    call timer_stop(ecosys_comp_CO3terms_timer, block_id=bid)

    call comp_co3_sat_vals(dkm, mask, pressure_correct, temperature, salinity, &
                           press_bar, CO3_sat_calcite, CO3_sat_aragonite)
       
    co3_calcite_anom   = CO3 - CO3_sat_calcite
    co3_aragonite_anom = CO3 - CO3_sat_aragonite

    do k=1,dkm

       if (k == 1) then
          ! set to -1, i.e. depth not found yet,
          ! if mask == .true. and surface supersaturated to -1
          zsat_calcite(k) = merge(-c1, c0, column_land_mask .and. CO3(k) > CO3_sat_calcite(k))
       else
          zsat_calcite(k) = zsat_calcite(k-1)
          if (zsat_calcite(k) == -c1 .and. CO3(k) <= CO3_sat_calcite(k)) then
             zsat_calcite(k) = zt(k-1) + (zt(k) - zt(k-1)) * &
                  co3_calcite_anom(k-1) / (co3_calcite_anom(k-1) - co3_calcite_anom(k))
          end if
          if (zsat_calcite(k) == -c1 .and. column_kmt == k) then
             zsat_calcite(k) = zw(k)
          end if
       end if

       ! -------------------
       if (k == 1) then
          ! set to -1, i.e. depth not found yet,
          ! if mask == .true. and surface supersaturated to -1
          zsat_aragonite(k) = merge(-c1, c0, column_land_mask .and. CO3(k) > CO3_sat_aragonite(k))
       else
          zsat_aragonite(k) = zsat_aragonite(k-1)
          if (zsat_aragonite(k) == -c1 .and. CO3(k) <= CO3_sat_aragonite(k)) then
             zsat_aragonite(k) = zt(k-1) + (zt(k) - zt(k-1)) * &
                  co3_aragonite_anom(k-1) / (co3_aragonite_anom(k-1) - co3_aragonite_anom(k))
          end if
          if (zsat_aragonite(k) == -c1 .and. column_kmt == k) then
             zsat_aragonite(k) = zw(k)
          end if
       end if
    enddo
    end associate

  end subroutine compute_carbonate_chemistry

  !***********************************************************************

  subroutine compute_function_scaling(column_temperature, Tfunc )
    !-----------------------------------------------------------------------
    !  Tref = 30.0 reference temperature (deg. C)
    !  Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
    !  growth, mort and grazing rates scaled by Tfunc where they are computed
    !-----------------------------------------------------------------------
    use marbl_parms, only : Q_10
    use marbl_parms, only : Tref
    use constants  , only : T0_Kelvin
    use constants  , only : c10

    real(r8), intent(in)  :: column_temperature
    real(r8), intent(out) :: Tfunc

    Tfunc = Q_10**(((column_temperature + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

  end subroutine compute_function_scaling

  !***********************************************************************

  subroutine compute_Pprime(k, auto_cnt, auto_meta, &
       autotroph_local, column_temperature, autotroph_secondary_species)

    use marbl_share_mod, only : autotroph_type
    use grid           , only : zt
    use marbl_parms    , only : thres_z1
    use marbl_parms    , only : thres_z2

    integer(int_kind)                          , intent(in)  :: k
    integer(int_kind)                          , intent(in)  :: auto_cnt
    type(autotroph_type)                       , intent(in)  :: auto_meta(auto_cnt)
    type(autotroph_local_type)          , intent(in)  :: autotroph_local(auto_cnt)
    real(r8)                                   , intent(in)  :: column_temperature
    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind
    real(r8) :: f_loss_thres
    real(r8) :: C_loss_thres

    associate(                                                       &
         Pprime => autotroph_secondary_species(:)%Pprime     & ! output
         )

    !  calculate the loss threshold interpolation factor
    if (zt(k) > thres_z1) then
       if (zt(k) < thres_z2) then
          f_loss_thres = (thres_z2 - zt(k))/(thres_z2 - thres_z1)
       else
          f_loss_thres = c0
       endif
    else
       f_loss_thres = c1
    endif

    !  Compute Pprime for all autotrophs, used for loss terms
    do auto_ind = 1, auto_cnt
       if (column_temperature < auto_meta(auto_ind)%temp_thres) then
          C_loss_thres = f_loss_thres * auto_meta(auto_ind)%loss_thres2
       else
          C_loss_thres = f_loss_thres * auto_meta(auto_ind)%loss_thres
       end if
       Pprime(auto_ind) = max(autotroph_local(auto_ind)%C - C_loss_thres, c0)
    end do

    end associate
  end subroutine compute_Pprime

  !***********************************************************************

  subroutine compute_Zprime(k, zoo_cnt, zoo_meta, zooC, &
       Tfunc, zooplankton_secondary_species)

    use constants       , only : c1, c0
    use marbl_share_mod , only : zooplankton_type
    use grid            , only : zt
    use marbl_parms     , only : thres_z1
    use marbl_parms     , only : thres_z2

    integer(int_kind)                            , intent(in)    :: k
    integer(int_kind)                            , intent(in)    :: zoo_cnt
    type(zooplankton_type)                       , intent(in)    :: zoo_meta(zoo_cnt)
    real(r8)                                     , intent(in)    :: zooC(zoo_cnt)
    real(r8)                                     , intent(in)    :: Tfunc
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)

    integer  :: zoo_ind
    real(r8) :: f_loss_thres
    real(r8) :: C_loss_thres

    associate(                                                       &
         zoo_loss => zooplankton_secondary_species(:)%zoo_loss,  & ! output
         Zprime   => zooplankton_secondary_species(:)%Zprime     & ! output
         )

    !  calculate the loss threshold interpolation factor
    if (zt(k) > thres_z1) then
       if (zt(k) < thres_z2) then
          f_loss_thres = (thres_z2 - zt(k))/(thres_z2 - thres_z1)
       else
          f_loss_thres = c0
       endif
    else
       f_loss_thres = c1
    endif

    do zoo_ind = 1, zoo_cnt
       C_loss_thres = f_loss_thres * zoo_meta(zoo_ind)%loss_thres
       Zprime(zoo_ind) = max(zooC(zoo_ind) - C_loss_thres, c0)

       zoo_loss(zoo_ind) = ( zoo_meta(zoo_ind)%z_mort2_0 * Zprime(zoo_ind)**1.5_r8 + &
            zoo_meta(zoo_ind)%z_mort_0  * Zprime(zoo_ind)) * Tfunc
    end do

    end associate
  end subroutine compute_Zprime

  !***********************************************************************

  subroutine compute_autotroph_uptake (auto_cnt, auto_meta, &
       tracer_local, &
       autotroph_secondary_species)

    use constants       , only : c1
    use marbl_share_mod , only : autotroph_type

    integer(int_kind)                 , intent(in)  :: auto_cnt
    type(autotroph_type)              , intent(in)  :: auto_meta(auto_cnt)
    real(r8), intent(in) :: tracer_local(ecosys_tracer_cnt)

    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind

    associate(                                               &
         DOP_loc => tracer_local(dop_ind), &
         NO3_loc => tracer_local(no3_ind), &
         NH4_loc => tracer_local(nh4_ind), &
         PO4_loc => tracer_local(po4_ind), &
         Fe_loc   => tracer_local(fe_ind),       &
         SiO3_loc => tracer_local(sio3_ind),     &
         VNO3  => autotroph_secondary_species(:)%VNO3  , & ! output
         VNH4  => autotroph_secondary_species(:)%VNH4  , & ! output
         VNtot => autotroph_secondary_species(:)%VNtot , & ! output
         VFe   => autotroph_secondary_species(:)%VFe   , & ! output
         f_nut => autotroph_secondary_species(:)%f_nut , & ! output
         VDOP  => autotroph_secondary_species(:)%VDOP  , & ! output
         VPO4  => autotroph_secondary_species(:)%VPO4  , & ! output
         VPtot => autotroph_secondary_species(:)%VPtot , & ! output
         VSiO3 => autotroph_secondary_species(:)%VSiO3   & ! output
         )

    !-----------------------------------------------------------------------
    !  Get relative nutrient uptake rates for autotrophs,
    !  min. relative uptake rate modifies C fixation in the manner
    !  that the min. cell quota does in GD98.
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt

       VNO3(auto_ind) = (NO3_loc / auto_meta(auto_ind)%kNO3) / (c1 + (NO3_loc / auto_meta(auto_ind)%kNO3) + (NH4_loc / auto_meta(auto_ind)%kNH4))
       VNH4(auto_ind) = (NH4_loc / auto_meta(auto_ind)%kNH4) / (c1 + (NO3_loc / auto_meta(auto_ind)%kNO3) + (NH4_loc / auto_meta(auto_ind)%kNH4))
       VNtot(auto_ind) = VNO3(auto_ind) + VNH4(auto_ind)
       if (auto_meta(auto_ind)%Nfixer) then
          VNtot(auto_ind) = c1
       end if

       VFe(auto_ind) = Fe_loc / (Fe_loc + auto_meta(auto_ind)%kFe)

       VPO4(auto_ind) = (PO4_loc / auto_meta(auto_ind)%kPO4) / (c1 + (PO4_loc / auto_meta(auto_ind)%kPO4) + (DOP_loc / auto_meta(auto_ind)%kDOP))
       VDOP(auto_ind) = (DOP_loc / auto_meta(auto_ind)%kDOP) / (c1 + (PO4_loc / auto_meta(auto_ind)%kPO4) + (DOP_loc / auto_meta(auto_ind)%kDOP))
       VPtot(auto_ind) = VPO4(auto_ind) + VDOP(auto_ind)

       if (auto_meta(auto_ind)%kSiO3 > c0) then
          VSiO3(auto_ind) = SiO3_loc / (SiO3_loc + auto_meta(auto_ind)%kSiO3)
       endif

       f_nut(auto_ind) = min(VNtot(auto_ind), VFe(auto_ind))
       f_nut(auto_ind) = min(f_nut(auto_ind), VPO4(auto_ind))
       if (auto_meta(auto_ind)%kSiO3 > c0) then
          f_nut(auto_ind) = min(f_nut(auto_ind), VSiO3(auto_ind))
       endif

    end do

    end associate

  end subroutine compute_autotroph_uptake

  !***********************************************************************

  subroutine compute_autotroph_photosynthesis (auto_cnt, auto_meta, &
       autotroph_loc, temperature, Tfunc, PAR_avg, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !     get photosynth. rate, phyto C biomass change, photoadapt
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : autotroph_type
    use constants       , only : c0, c1
    use marbl_parms     , only : epsTinv

    integer(int_kind)                          , intent(in)    :: auto_cnt
    type(autotroph_type)                       , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)          , intent(in)    :: autotroph_loc(auto_cnt)
    real(r8)                                   , intent(in)    :: temperature
    real(r8)                                   , intent(in)    :: Tfunc
    real(r8)                                   , intent(in)    :: PAR_avg
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer   :: auto_ind
    real (r8) :: PCmax  ! max value of PCphoto at temperature TEMP (1/sec)

    associate(                                               &
         thetaC    => autotroph_secondary_species(:)%thetaC,    & ! local Chl/C ratio (mg Chl/mmol C)
         f_nut     => autotroph_secondary_species(:)%f_nut,     & ! input
         light_lim => autotroph_secondary_species(:)%light_lim, & ! output
         PCPhoto   => autotroph_secondary_species(:)%PCPhoto,   & ! output
         photoC    => autotroph_secondary_species(:)%photoC     & ! output
         )

    do auto_ind = 1, auto_cnt

       PCmax = auto_meta(auto_ind)%PCref * f_nut(auto_ind) * Tfunc
       if (temperature < autotrophs(auto_ind)%temp_thres) then
          PCmax = c0
       end if
       light_lim(auto_ind) = (c1 - exp((-c1 * auto_meta(auto_ind)%alphaPI * thetaC(auto_ind) * PAR_avg) / (PCmax + epsTinv)))
       PCphoto(auto_ind) = PCmax * light_lim(auto_ind)
       photoC(auto_ind) = PCphoto(auto_ind) * autotroph_loc(auto_ind)%C

    end do

    end associate

  end subroutine compute_autotroph_photosynthesis

  !***********************************************************************

  subroutine compute_autotroph_phyto_diatoms (auto_cnt, auto_meta, &
       autotroph_loc, PAR_avg, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !  Get nutrient uptakes by small phyto based on calculated C fixation
    !  total N uptake VNC is used in photoadaption
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : autotroph_type
    use constants       , only : c0
    use marbl_parms     , only : Q

    integer(int_kind)                          , intent(in)    :: auto_cnt
    type(autotroph_type)                       , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)          , intent(in)    :: autotroph_loc(auto_cnt)
    real(r8)                                   , intent(in)    :: PAR_avg
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind
    real(r8) :: WORK1 ! temporary
    real(r8) :: pChl             ! Chl synth. regulation term (mg Chl/mmol N)

    associate(                                            &
         thetaC   => autotroph_secondary_species(:)%thetaC , & ! local Chl/C ratio (mg Chl/mmol C)
         gQfe     => autotroph_secondary_species(:)%gQfe,    & ! fe/C for growth
         gQsi     => autotroph_secondary_species(:)%gQsi,    & ! diatom Si/C ratio for growth (new biomass)
         VNO3     => autotroph_secondary_species(:)%VNO3,    & ! input
         VNH4     => autotroph_secondary_species(:)%VNH4,    & ! input
         VNtot    => autotroph_secondary_species(:)%VNtot,   & ! input
         VPO4     => autotroph_secondary_species(:)%VPO4,    & ! input
         VDOP     => autotroph_secondary_species(:)%VDOP,    & ! input
         VPtot    => autotroph_secondary_species(:)%VPtot,   & ! input
         PCphoto  => autotroph_secondary_species(:)%PCphoto, & ! input
         photoC   => autotroph_secondary_species(:)%photoC,  & ! input
         NO3_V    => autotroph_secondary_species(:)%NO3_V,   & ! output
         NH4_V    => autotroph_secondary_species(:)%NH4_V,   & ! output
         VNC      => autotroph_secondary_species(:)%VNC,     & ! output
         PO4_V    => autotroph_secondary_species(:)%PO4_V,   & ! output
         DOP_V    => autotroph_secondary_species(:)%DOP_V,   & ! output
         photoFe  => autotroph_secondary_species(:)%photoFe, & ! output
         photoSi  => autotroph_secondary_species(:)%photoSi, & ! output
         photoacc => autotroph_secondary_species(:)%photoacc & ! output
         )

    do auto_ind = 1, auto_cnt

       if (VNtot(auto_ind) > c0) then
          NO3_V(auto_ind) = (VNO3(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
          NH4_V(auto_ind) = (VNH4(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
          VNC(auto_ind) = PCphoto(auto_ind) * Q
       else
          NO3_V(auto_ind) = c0
          NH4_V(auto_ind) = c0
          VNC(auto_ind)   = c0
       end if

       if (VPtot(auto_ind) > c0) then
          PO4_V(auto_ind) = (VPO4(auto_ind) / VPtot(auto_ind)) * photoC(auto_ind) * auto_meta(auto_ind)%Qp
          DOP_V(auto_ind) = (VDOP(auto_ind) / VPtot(auto_ind)) * photoC(auto_ind) * auto_meta(auto_ind)%Qp
       else
          PO4_V(auto_ind) = c0
          DOP_V(auto_ind) = c0
       end if

       !-----------------------------------------------------------------------
       !  Get nutrient uptake by diatoms based on C fixation
       !-----------------------------------------------------------------------

       photoFe(auto_ind) = photoC(auto_ind) * gQfe(auto_ind)

       if (autotrophs(auto_ind)%Si_ind > 0) then
          photoSi(auto_ind) = photoC(auto_ind) * gQsi(auto_ind)
       endif

       !-----------------------------------------------------------------------
       !  calculate pChl, (used in photoadapt., GD98)
       !  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
       !  GD 98 Chl. synth. term
       !-----------------------------------------------------------------------

       WORK1 = autotrophs(auto_ind)%alphaPI * thetaC(auto_ind) * PAR_avg
       if (WORK1 > c0) then
          pChl = autotrophs(auto_ind)%thetaN_max * PCphoto(auto_ind) / WORK1
          photoacc(auto_ind) = (pChl * VNC(auto_ind) / thetaC(auto_ind)) * autotroph_loc(auto_ind)%Chl
       else
          photoacc(auto_ind) = c0
       end if

    end do

    end associate

  end subroutine compute_autotroph_phyto_diatoms

  !***********************************************************************

  subroutine compute_autotroph_calcification (auto_cnt, auto_meta, &
       autotroph_loc, temperature, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !  CaCO3 Production, parameterized as function of small phyto production
    !  decrease CaCO3 as function of nutrient limitation decrease CaCO3 prod
    !  at low temperatures increase CaCO3 prod under bloom conditions
    !  maximum calcification rate is 40% of primary production
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : autotroph_type
    use marbl_parms     , only : parm_f_prod_sp_CaCO3
    use marbl_parms     , only : CaCO3_sp_thres
    use marbl_parms     , only : CaCO3_temp_thres1
    use marbl_parms     , only : CaCO3_temp_thres2
    use marbl_parms     , only : f_photosp_CaCO3

    integer(int_kind)                          , intent(in)    :: auto_cnt
    type(autotroph_type)                       , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)          , intent(in)    :: autotroph_loc(auto_cnt)
    real(r8)                                   , intent(in)    :: temperature
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind

    associate(                                                       &
         f_nut      => autotroph_secondary_species(:)%f_nut,     & ! input
         photoC     => autotroph_secondary_species(:)%photoC,    & ! input
         CaCO3_PROD => autotroph_secondary_species(:)%CaCO3_PROD & ! output
         )

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%imp_calcifier) then
          CaCO3_PROD(auto_ind) = parm_f_prod_sp_CaCO3 * photoC(auto_ind)
          CaCO3_PROD(auto_ind) = CaCO3_PROD(auto_ind) * f_nut(auto_ind)

          if (temperature < CaCO3_temp_thres1)  then
             CaCO3_PROD(auto_ind) = CaCO3_PROD(auto_ind) * max((temperature - CaCO3_temp_thres2), c0) / &
                  (CaCO3_temp_thres1-CaCO3_temp_thres2)
          end if

          if (autotroph_loc(auto_ind)%C > CaCO3_sp_thres) then
             CaCO3_PROD(auto_ind) = min((CaCO3_PROD(auto_ind) * autotroph_loc(auto_ind)%C / CaCO3_sp_thres), &
                  (f_photosp_CaCO3 * photoC(auto_ind)))
          end if
       end if
    end do

    end associate
  end subroutine compute_autotroph_calcification

  !***********************************************************************

  subroutine compute_autotroph_nfixation (auto_cnt, auto_meta, &
       autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !  Get N fixation by diazotrophs based on C fixation,
    !  Diazotrophs fix more than they need then 20% is excreted
    !-----------------------------------------------------------------------

    use marbl_share_mod , only : autotroph_type
    use marbl_parms     , only : Q
    use marbl_parms     , only : r_Nfix_photo

    integer(int_kind)                          , intent(in)  :: auto_cnt
    type(autotroph_type)                       , intent(in)  :: auto_meta(auto_cnt)
    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind
    real(r8) :: work1

    associate(                                                   &
         photoC   => autotroph_secondary_species(:)%photoC,  & ! input
         NO3_V    => autotroph_secondary_species(:)%NO3_V ,  & ! input
         NH4_V    => autotroph_secondary_species(:)%NH4_V ,  & ! input
         Nfix     => autotroph_secondary_species(:)%Nfix  ,  & ! output total Nitrogen fixation (mmol N/m^3/sec)
         Nexcrete => autotroph_secondary_species(:)%Nexcrete & ! output fixed N excretion
         )

    do auto_ind = 1, autotroph_cnt

       if (auto_meta(auto_ind)%Nfixer) then
          work1 = photoC(auto_ind) * Q
          Nfix(auto_ind) = (work1 * r_Nfix_photo) - NO3_V(auto_ind) - NH4_V(auto_ind)
          Nexcrete(auto_ind) = Nfix(auto_ind) + NO3_V(auto_ind) + NH4_V(auto_ind) - work1
       endif
    end do

    end associate
  end subroutine compute_autotroph_nfixation

  !***********************************************************************

  subroutine compute_autotroph_loss (auto_cnt, auto_meta, &
       Tfunc, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    ! Compute autotroph-loss, autotroph aggregation loss and routine of
    ! loss terms
    !-----------------------------------------------------------------------

    use marbl_parms     , only : dps
    use marbl_share_mod , only : autotroph_type

    integer(int_kind)                          , intent(in)    :: auto_cnt
    type(autotroph_type)                       , intent(in)    :: auto_meta(auto_cnt)
    real(r8)                                   , intent(in)    :: Tfunc
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind

    associate(                                                               &
         QCaCO3        => autotroph_secondary_species(:)%QCaCO3        , & ! input
         Pprime        => autotroph_secondary_species(:)%Pprime        , & ! input
         auto_loss     => autotroph_secondary_species(:)%auto_loss     , & ! output
         auto_loss_poc => autotroph_secondary_species(:)%auto_loss_poc , & ! output
         auto_loss_dic => autotroph_secondary_species(:)%auto_loss_dic , & ! output
         auto_loss_doc => autotroph_secondary_species(:)%auto_loss_doc , & ! output
         auto_agg      => autotroph_secondary_species(:)%auto_agg        & ! output
         )

    do auto_ind = 1, autotroph_cnt
       !-----------------------------------------------------------------------
       !  get autotroph loss (in C units)
       !  autotroph agg loss
       !-----------------------------------------------------------------------

       auto_loss(auto_ind) = auto_meta(auto_ind)%mort * Pprime(auto_ind) * Tfunc

       auto_agg(auto_ind) = min((auto_meta(auto_ind)%agg_rate_max * dps) * Pprime(auto_ind), &
            auto_meta(auto_ind)%mort2 * Pprime(auto_ind) * Pprime(auto_ind))
       auto_agg(auto_ind) = max((auto_meta(auto_ind)%agg_rate_min * dps) * Pprime(auto_ind), auto_agg(auto_ind))

       !-----------------------------------------------------------------------
       !  routing of loss terms
       !  all aggregation goes to POC
       !  min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
       !-----------------------------------------------------------------------

       if (auto_meta(auto_ind)%imp_calcifier) then
          auto_loss_poc(auto_ind) = QCaCO3(auto_ind) * auto_loss(auto_ind)
       else
          auto_loss_poc(auto_ind) = auto_meta(auto_ind)%loss_poc * auto_loss(auto_ind)
       endif
       auto_loss_doc(auto_ind) = (c1 - parm_labile_ratio) * (auto_loss(auto_ind) - auto_loss_poc(auto_ind))
       auto_loss_dic(auto_ind) = parm_labile_ratio * (auto_loss(auto_ind) - auto_loss_poc(auto_ind))
    end do  ! auto_ind = 1, autotroph_cnt

    end associate
  end subroutine compute_autotroph_loss

  !***********************************************************************

  subroutine compute_grazing (auto_cnt, zoo_cnt, grazer_prey_cnt, auto_meta, &
       Tfunc, zooplankton_loc, &
       zooplankton_secondary_species, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !  CALCULATE GRAZING
    !
    !  Autotroph prey
    !  routing of grazing terms
    !  all aggregation goes to POC
    !  currently assumes that 33% of grazed caco3 is remineralized
    !  if autotrophs(sp_ind)%graze_zoo ever changes, coefficients on routing grazed sp must change!
    !  min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
    !  NOTE: if autotrophs(diat_ind)%graze_zoo is changed, coeff.s for poc, doc and dic must change!
    !-----------------------------------------------------------------------

    use marbl_parms     , only : epsC
    use marbl_parms     , only : epsTinv
    use marbl_share_mod , only : autotroph_type
    use marbl_share_mod , only : zooplankton_type
    use marbl_share_mod , only : grazing
    use marbl_parms     , only : grz_fnc_michaelis_menten
    use marbl_parms     , only : grz_fnc_sigmoidal
    use constants       , only : c0

    integer(int_kind)                            , intent(in)    :: auto_cnt
    integer(int_kind)                            , intent(in)    :: zoo_cnt
    integer(int_kind)                            , intent(in)    :: grazer_prey_cnt
    type(autotroph_type)                         , intent(in)    :: auto_meta(auto_cnt)
    real(r8)                                     , intent(in)    :: Tfunc
    type(zooplankton_local_type)          , intent(in)    :: zooplankton_loc(zoo_cnt)
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind, auto_ind2
    integer  :: zoo_ind, zoo_ind2
    integer  :: pred_ind
    integer  :: prey_ind
    real(r8) :: work1, work2, work3, work4
    real(r8) :: graze_rate

    associate(                                                                  &
         Pprime         => autotroph_secondary_species(:)%Pprime          , & ! input
         QCaCO3         => autotroph_secondary_species(:)%QCaCO3          , & ! input
         Zprime         => zooplankton_secondary_species(:)%Zprime        , & ! input
         auto_graze     => autotroph_secondary_species(:)%auto_graze      , & ! output
         auto_graze_poc => autotroph_secondary_species(:)%auto_graze_poc  , & ! output
         auto_graze_dic => autotroph_secondary_species(:)%auto_graze_dic  , & ! output
         auto_graze_doc => autotroph_secondary_species(:)%auto_graze_doc  , & ! output
         auto_graze_zoo => autotroph_secondary_species(:)%auto_graze_zoo  , & ! output
         zoo_graze      => zooplankton_secondary_species(:)%zoo_graze     , & ! output
         zoo_graze_poc  => zooplankton_secondary_species(:)%zoo_graze_poc , & ! output
         zoo_graze_dic  => zooplankton_secondary_species(:)%zoo_graze_dic , & ! output
         zoo_graze_doc  => zooplankton_secondary_species(:)%zoo_graze_doc , & ! output
         zoo_graze_zoo  => zooplankton_secondary_species(:)%zoo_graze_zoo , & ! output
         x_graze_zoo    => zooplankton_secondary_species(:)%x_graze_zoo   , & ! output
         f_zoo_detr     => zooplankton_secondary_species(:)%f_zoo_detr      & ! output
         )

    auto_graze(:)     = c0 ! total grazing losses from autotroph pool at auto_ind
    auto_graze_zoo(:) = c0 ! autotroph grazing losses routed to zooplankton at auto_ind
    auto_graze_poc(:) = c0 ! autotroph grazing losses routed to poc
    auto_graze_doc(:) = c0 ! autotroph grazing losses routed to doc
    auto_graze_dic(:) = c0 ! autotroph grazing losses routed to dic (computed by residual)

    zoo_graze(:)     = c0 ! total grazing losses from zooplankton pool at zoo_ind
    zoo_graze_zoo(:) = c0 ! zooplankton grazing losses routed to zooplankton at zoo_ind
    zoo_graze_poc(:) = c0 ! zooplankton grazing losses routed to poc
    zoo_graze_doc(:) = c0 ! zooplankton grazing losses routed to doc
    zoo_graze_dic(:) = c0 ! zooplankton grazing losses routed to dic (computed by residual)

    x_graze_zoo(:)   = c0 ! grazing gains by zooplankton at zoo_ind

    do pred_ind = 1, zoo_cnt

       work3 = c0
       work4 = c0

       do prey_ind = 1, grazer_prey_cnt

          !-----------------------------------------------------------------------
          !  compute sum of carbon in the grazee class, both autotrophs and zoop
          !-----------------------------------------------------------------------
          work1 = c0 ! biomass in prey class prey_ind
          do auto_ind2 = 1, grazing(prey_ind, pred_ind)%auto_ind_cnt
             auto_ind = grazing(prey_ind, pred_ind)%auto_ind(auto_ind2)
             work1 = work1 + Pprime(auto_ind)
          end do

          do zoo_ind2 = 1, grazing(prey_ind, pred_ind)%zoo_ind_cnt
             zoo_ind = grazing(prey_ind, pred_ind)%zoo_ind(zoo_ind2)
             work1 = work1 + Zprime(zoo_ind)
          end do

          ! compute grazing rate
          select case (grazing(prey_ind, pred_ind)%grazing_function)

          case (grz_fnc_michaelis_menten)

             if (work1 > c0) then
                graze_rate = grazing(prey_ind, pred_ind)%z_umax_0 * Tfunc * zooplankton_loc(pred_ind)%C &
                     * ( work1 / (work1 + grazing(prey_ind, pred_ind)%z_grz) )
             else
                graze_rate = c0
             end if

          case (grz_fnc_sigmoidal)

             if (work1 > c0) then
                graze_rate = grazing(prey_ind, pred_ind)%z_umax_0 * Tfunc * zooplankton_loc(pred_ind)%C &
                     * ( work1**2 / (work1**2 + grazing(prey_ind, pred_ind)%z_grz**2) )
             else
                graze_rate = c0
             end if

          end select

          !-----------------------------------------------------------------------
          !  autotroph prey
          !-----------------------------------------------------------------------

          do auto_ind2 = 1, grazing(prey_ind, pred_ind)%auto_ind_cnt
             auto_ind = grazing(prey_ind, pred_ind)%auto_ind(auto_ind2)

             ! scale by biomass from autotroph pool
             if (work1 > c0) then
                work2 = (Pprime(auto_ind) / work1) * graze_rate ! total grazing loss from auto_ind
             else
                work2 = c0
             end if
             auto_graze(auto_ind) = auto_graze(auto_ind) + work2

             ! routed to zooplankton
             auto_graze_zoo(auto_ind) = auto_graze_zoo(auto_ind) + grazing(prey_ind, pred_ind)%graze_zoo * work2
             x_graze_zoo(pred_ind)    = x_graze_zoo(pred_ind)    + grazing(prey_ind, pred_ind)%graze_zoo * work2

             ! routed to POC
             if (auto_meta(auto_ind)%imp_calcifier) then
                auto_graze_poc(auto_ind) = auto_graze_poc(auto_ind) &
                     + work2 * max((caco3_poc_min * QCaCO3(auto_ind)),  &
                     min(spc_poc_fac * max(1.0_r8, Pprime(auto_ind)),    &
                     f_graze_sp_poc_lim))
             else
                auto_graze_poc(auto_ind) = auto_graze_poc(auto_ind) + grazing(prey_ind, pred_ind)%graze_poc * work2
             endif

             ! routed to DOC
             auto_graze_doc(auto_ind) = auto_graze_doc(auto_ind) + grazing(prey_ind, pred_ind)%graze_doc * work2

             !  get fractional factor for routing of zoo losses, based on food supply
             work3 = work3 + grazing(prey_ind, pred_ind)%f_zoo_detr * (work2 + epsC * epsTinv)
             work4 = work4 + (work2 + epsC * epsTinv)

          end do

          !-----------------------------------------------------------------------
          !  Zooplankton prey
          !-----------------------------------------------------------------------
          do zoo_ind2 = 1, grazing(prey_ind, pred_ind)%zoo_ind_cnt
             zoo_ind = grazing(prey_ind, pred_ind)%zoo_ind(zoo_ind2)

             ! scale by biomass from zooplankton pool
             if (work1 > c0) then
                work2 = (Zprime(zoo_ind) / work1) * graze_rate
             else
                work2 = c0
             end if

             ! grazing loss from zooplankton prey pool
             zoo_graze(zoo_ind) = zoo_graze(zoo_ind) + work2

             ! routed to zooplankton
             zoo_graze_zoo(zoo_ind) = zoo_graze_zoo(zoo_ind) + grazing(prey_ind, pred_ind)%graze_zoo * work2
             x_graze_zoo(pred_ind) = x_graze_zoo(pred_ind)   + grazing(prey_ind, pred_ind)%graze_zoo * work2

             ! routed to POC/DOC
             zoo_graze_poc(zoo_ind) = zoo_graze_poc(zoo_ind) + grazing(prey_ind, pred_ind)%graze_poc * work2
             zoo_graze_doc(zoo_ind) = zoo_graze_doc(zoo_ind) + grazing(prey_ind, pred_ind)%graze_doc * work2

             !  get fractional factor for routing of zoo losses, based on food supply
             work3 = work3 + grazing(prey_ind, pred_ind)%f_zoo_detr * (work2 + epsC * epsTinv)
             work4 = work4 + (work2 + epsC * epsTinv)

          end do
       end do

       f_zoo_detr(pred_ind) = work3 / work4
    end do

    end associate

  end subroutine compute_grazing

  !***********************************************************************

  subroutine compute_routing (auto_cnt, zoo_cnt,  auto_meta, &
       zooplankton_secondary_species, autotroph_secondary_species)

    use constants       , only : c1
    use marbl_parms     , only : Qp_zoo_pom
    use marbl_parms     , only : parm_labile_ratio
    use marbl_share_mod , only : autotroph_type

    integer(int_kind)                            , intent(in)    :: auto_cnt
    integer(int_kind)                            , intent(in)    :: zoo_cnt
    type(autotroph_type)                         , intent(in)    :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(inout) :: autotroph_secondary_species(auto_cnt)

    integer  :: auto_ind, zoo_ind
    real(r8) :: remaining_P      ! used in routing P from autotrophs w/ Qp different from Qp_zoo_pom

    associate(                                                                   &
         auto_graze      => autotroph_secondary_species(:)%auto_graze      , & ! input
         auto_graze_zoo  => autotroph_secondary_species(:)%auto_graze_zoo  , & ! input
         auto_graze_poc  => autotroph_secondary_species(:)%auto_graze_poc  , & ! input
         auto_graze_doc  => autotroph_secondary_species(:)%auto_graze_doc  , & ! input
         auto_loss       => autotroph_secondary_species(:)%auto_loss       , & ! input
         auto_loss_poc   => autotroph_secondary_species(:)%auto_loss_poc   , & ! input
         auto_agg        => autotroph_secondary_species(:)%auto_agg        , & ! input
         zoo_graze       => zooplankton_secondary_species(:)%zoo_graze     , & ! input
         zoo_graze_poc   => zooplankton_secondary_species(:)%zoo_graze_poc , & ! input
         zoo_graze_doc   => zooplankton_secondary_species(:)%zoo_graze_doc , & ! input
         zoo_graze_zoo   => zooplankton_secondary_species(:)%zoo_graze_zoo , & ! input
         zoo_loss        => zooplankton_secondary_species(:)%zoo_loss      , & ! input
         f_zoo_detr      => zooplankton_secondary_species(:)%f_zoo_detr    , & ! input

         auto_graze_dic  => autotroph_secondary_species(:)%auto_graze_dic  , & ! output
         remaining_P_dop => autotroph_secondary_species(:)%remaining_P_dop , & ! output
         remaining_P_dip => autotroph_secondary_species(:)%remaining_P_dip , & ! output
         zoo_graze_dic   => zooplankton_secondary_species(:)%zoo_graze_dic , & ! output
         zoo_loss_poc    => zooplankton_secondary_species(:)%zoo_loss_poc  , & ! output
         zoo_loss_doc    => zooplankton_secondary_species(:)%zoo_loss_doc  , & ! output
         zoo_loss_dic    => zooplankton_secondary_species(:)%zoo_loss_dic    & ! output
         )

    !-----------------------------------------------------------------------
    ! compute routing to dic of grazed material
    ! call this and the one below compute_routing
    !-----------------------------------------------------------------------
    do auto_ind = 1, auto_cnt
       auto_graze_dic(auto_ind) = auto_graze(auto_ind) &
            - (auto_graze_zoo(auto_ind) + auto_graze_poc(auto_ind) + auto_graze_doc(auto_ind))
    end do
    do zoo_ind = 1, zoo_cnt
       zoo_graze_dic(zoo_ind) = zoo_graze(zoo_ind)  &
            - (zoo_graze_zoo(zoo_ind) + zoo_graze_poc(zoo_ind) + zoo_graze_doc(zoo_ind))
    end do

    !-----------------------------------------------------------------------
    ! compute zooplankton loss routing
    ! call this compute_routing_zooplankton_loss
    !-----------------------------------------------------------------------
    do zoo_ind = 1, zoo_cnt
       zoo_loss_poc(zoo_ind) = f_zoo_detr(zoo_ind) * zoo_loss(zoo_ind)
       zoo_loss_doc(zoo_ind) = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr(zoo_ind)) * zoo_loss(zoo_ind)
       zoo_loss_dic(zoo_ind) = parm_labile_ratio * (c1 - f_zoo_detr(zoo_ind)) * zoo_loss(zoo_ind)
    end do

    !-----------------------------------------------------------------------
    ! P from some autotrophs w/ Qp different from Qp_zoo_pom must be routed differently than other
    ! elements to ensure that sinking detritus and zooplankton pools get their fixed P/C ratios.
    ! The remaining P is split evenly between DOP and PO4.
    !-----------------------------------------------------------------------
    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%Qp /= Qp_zoo_pom) then
          remaining_P = ((auto_graze(auto_ind) + auto_loss(auto_ind) + auto_agg(auto_ind)) * auto_meta(auto_ind)%Qp) &
               - ((auto_graze_zoo(auto_ind)) * Qp_zoo_pom) &
               - ((auto_graze_poc(auto_ind) + auto_loss_poc(auto_ind) + auto_agg(auto_ind)) * Qp_zoo_pom)
          remaining_P_dop(auto_ind) = (c1 - parm_labile_ratio) * remaining_P
          remaining_P_dip(auto_ind) = parm_labile_ratio * remaining_P
       endif
    end do

    end associate

  end subroutine compute_routing

  !***********************************************************************

  subroutine compute_dissolved_organic_matter (auto_cnt, zoo_cnt, auto_meta, &
       zooplankton_secondary_species, autotroph_secondary_species, &
       PAR_avg, tracer_local, &
       dissolved_organic_matter)

    use marbl_share_mod , only : autotroph_type
    use marbl_parms     , only : DOC_reminR
    use marbl_parms     , only : DOFe_reminR
    use marbl_parms     , only : DON_reminR
    use marbl_parms     , only : DONr_reminR
    use marbl_parms     , only : DONrefract
    use marbl_parms     , only : DOP_reminR
    use marbl_parms     , only : DOPr_reminR
    use marbl_parms     , only : Qfe_zoo
    use marbl_parms     , only : Qp_zoo_pom
    use marbl_parms     , only : Q

    integer                                      , intent(in)  :: auto_cnt
    integer                                      , intent(in)  :: zoo_cnt
    type(autotroph_type)                         , intent(in)  :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type) , intent(in)  :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(in)  :: autotroph_secondary_species(auto_cnt)
    real(r8)                                     , intent(in)  :: PAR_avg
    real(r8), intent(in) :: tracer_local(ecosys_tracer_cnt)

    type(dissolved_organic_matter_type)   , intent(out) :: dissolved_organic_matter

    integer :: auto_ind

    associate(                                                                   &
         DOC_loc  => tracer_local(doc_ind), &
         DON_loc  => tracer_local(don_ind), &
         DOFe_loc => tracer_local(dofe_ind), &
         DOP_loc  => tracer_local(dop_ind), &
         DONr_loc => tracer_local(donr_ind), &
         DOPr_loc => tracer_local(dopr_ind), &
         Qfe             => autotroph_secondary_species(:)%Qfe             , & ! input
         remaining_P_dop => autotroph_secondary_species(:)%remaining_P_dop , & ! input
         auto_loss_doc   => autotroph_secondary_species(:)%auto_loss_doc   , & ! input
         auto_graze_doc  => autotroph_secondary_species(:)%auto_graze_doc  , & ! inpug
         zoo_loss_doc    => zooplankton_secondary_species(:)%zoo_loss_doc  , & ! input
         zoo_graze_doc   => zooplankton_secondary_species(:)%zoo_graze_doc , & ! input

         DOC_prod        => dissolved_organic_matter%DOC_prod           , & ! output production of DOC (mmol C/m^3/sec)
         DOC_remin       => dissolved_organic_matter%DOC_remin          , & ! output remineralization of DOC (mmol C/m^3/sec)
         DON_prod        => dissolved_organic_matter%DON_prod           , & ! output production of dissolved organic N
         DON_remin       => dissolved_organic_matter%DON_remin          , & ! output portion of DON remineralized
         DOFe_prod       => dissolved_organic_matter%DOFe_prod          , & ! output produciton of dissolved organic Fe
         DOFe_remin      => dissolved_organic_matter%DOFe_remin         , & ! output portion of DOFe remineralized
         DOP_prod        => dissolved_organic_matter%DOP_prod           , & ! output production of dissolved organic P
         DOP_remin       => dissolved_organic_matter%DOP_remin          , & ! output portion of DOP remineralized
         DONr_remin      => dissolved_organic_matter%DONr_remin         , & ! output portion of refractory DON remineralized
         DOPr_remin      => dissolved_organic_matter%DOPr_remin           & ! output portion of refractory DOP remineralized
         )

      !-----------------------------------------------------------------------
      !  compute terms for DOM
      !-----------------------------------------------------------------------

      DOC_prod = sum(zoo_loss_doc(:)) + sum(auto_loss_doc(:)) + sum(auto_graze_doc(:)) + sum(zoo_graze_doc(:))
      DON_prod = Q * DOC_prod
      DOP_prod = Qp_zoo_pom * ( sum(zoo_loss_doc(:)) + sum(zoo_graze_doc(:)) )
      do auto_ind = 1, auto_cnt
         if (auto_meta(auto_ind)%Qp == Qp_zoo_pom) then
            DOP_prod = DOP_prod + auto_meta(auto_ind)%Qp * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind))
         else
            DOP_prod = DOP_prod + remaining_P_dop(auto_ind)
         endif
      end do
      DOFe_prod = Qfe_zoo * ( sum(zoo_loss_doc(:)) + sum(zoo_graze_doc(:)) )
      do auto_ind = 1, auto_cnt
         DOFe_prod = DOFe_prod + Qfe(auto_ind) * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind))
      end do

      DOC_remin  = DOC_loc  * DOC_reminR
      DON_remin  = DON_loc  * DON_reminR
      DOFe_remin = DOFe_loc * DOFe_reminR
      DOP_remin  = DOP_loc  * DOP_reminR

      !  Refractory remin rate due to photochemistry
      !  below euphotic zone remin rate sharply decrease

      if (PAR_avg > 1.0_r8) then
         DONr_remin = DONr_loc * DONr_reminR
         DOPr_remin = DOPr_loc * DOPr_reminR
      else
         DONr_remin = DONr_loc * (c1/(365.0_r8*670.0_r8)) * dps  ! 1/670 yrs
         DOPr_remin = DOPr_loc * (c1/(365.0_r8*460.0_r8)) * dps  ! 1/460 yrs
         DOC_remin  = DOC_remin * 0.0685_r8
         DON_remin  = DON_remin * 0.1_r8
         DOFe_remin = DOFe_remin * 0.05_r8
         DOP_remin  = DOP_remin * 0.05_r8
      end if

    end associate
  end subroutine compute_dissolved_organic_matter

  !***********************************************************************

  subroutine compute_large_detritus(k, auto_cnt, zoo_cnt, auto_meta, &
       zooplankton_secondary_species, autotroph_secondary_species, Fe_loc, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, &
       Fe_scavenge, Fe_scavenge_rate)

    use marbl_share_mod , only : autotroph_type
    use marbl_share_mod , only : column_sinking_particle_type
    use marbl_parms     , only : f_graze_CaCO3_remin
    use marbl_parms     , only : f_graze_si_remin
    use marbl_parms     , only : Qfe_zoo
    use marbl_parms     , only : parm_Fe_scavenge_rate0
    use marbl_parms     , only : dust_fescav_scale
    use marbl_parms     , only : Fe_scavenge_thres1
    use marbl_parms     , only : fe_max_scale2
    use marbl_parms     , only : yps

    integer                                      , intent(in)  :: k
    integer                                      , intent(in)  :: auto_cnt
    integer                                      , intent(in)  :: zoo_cnt
    type(autotroph_type)                         , intent(in)  :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type) , intent(in)  :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(in)  :: autotroph_secondary_species(auto_cnt)
    real(r8)                                     , intent(in)  :: Fe_loc
    type(column_sinking_particle_type)           , intent(out) :: POC
    type(column_sinking_particle_type)           , intent(out) :: P_CaCO3
    type(column_sinking_particle_type)           , intent(out) :: P_SiO2
    type(column_sinking_particle_type)           , intent(out) :: dust
    type(column_sinking_particle_type)           , intent(out) :: P_iron
    real(r8)                                     , intent(out) :: Fe_scavenge
    real(r8)                                     , intent(out) :: Fe_scavenge_rate

    integer :: auto_ind

    associate(                                                                 &
         QCaCO3         => autotroph_secondary_species(:)%QCaCO3         , & ! input
         Qsi            => autotroph_secondary_species(:)%Qsi            , & ! input
         Qfe            => autotroph_secondary_species(:)%Qfe            , & ! input
         auto_graze     => autotroph_secondary_species(:)%auto_graze     , & ! input
         auto_graze_poc => autotroph_secondary_species(:)%auto_graze_poc , & ! input
         auto_agg       => autotroph_secondary_species(:)%auto_agg       , & ! input
         auto_loss      => autotroph_secondary_species(:)%auto_loss      , & ! input
         auto_loss_poc  => autotroph_secondary_species(:)%auto_loss_poc  , & ! input
         zoo_loss_poc   => zooplankton_secondary_species(:)%zoo_loss_poc , & ! input
         zoo_graze_poc  => zooplankton_secondary_species(:)%zoo_graze_poc  & ! input
         )

    !-----------------------------------------------------------------------
    !  large detritus C
    !-----------------------------------------------------------------------

    POC%prod(k) = sum(zoo_loss_poc(:)) + sum(auto_graze_poc(:)) + sum(zoo_graze_poc(:)) &
         + sum(auto_agg(:)) + sum(auto_loss_poc(:))

    !-----------------------------------------------------------------------
    !  large detrital CaCO3
    !  33% of CaCO3 is remin when phyto are grazed
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%CaCO3_ind > 0) then
          P_CaCO3%prod(k) = ((c1 - f_graze_CaCO3_remin) * auto_graze(auto_ind) + &
               auto_loss(auto_ind) + auto_agg(auto_ind)) * QCaCO3(auto_ind)
       endif
    end do

    !-----------------------------------------------------------------------
    !  large detritus SiO2
    !  grazed diatom SiO2, 60% is remineralized
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%Si_ind > 0) then
          P_SiO2%prod(k) = Qsi(auto_ind) * ((c1 - f_graze_si_remin) * auto_graze(auto_ind) + auto_agg(auto_ind) &
               + auto_meta(auto_ind)%loss_poc * auto_loss(auto_ind))
       endif
    end do

    !-----------------------------------------------------------------------
    ! Dust
    !-----------------------------------------------------------------------
    dust%prod(k) = c0

    !-----------------------------------------------------------------------
    !  Compute iron scavenging :
    !  1) compute in terms of loss per year per unit iron (%/year/fe)
    !  2) scale by sinking POMx10 + Dust + bSi + CaCO3 flux
    !  3) increase scavenging at higher iron (>0.6nM)
    !  4) convert to net loss per second
    !-----------------------------------------------------------------------
    Fe_scavenge_rate = parm_Fe_scavenge_rate0

    Fe_scavenge_rate = Fe_scavenge_rate * &
         ((POC%sflux_out(k)     + POC%hflux_out(k)    ) * 120.1_r8 + &
          (P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)) * P_CaCO3%mass + &
          (P_SiO2%sflux_out(k)  + P_SiO2%hflux_out(k) ) * P_SiO2%mass + &
          (dust%sflux_out(k)    + dust%hflux_out(k)   ) * dust_fescav_scale)

    if (Fe_loc > Fe_scavenge_thres1) then
       Fe_scavenge_rate = Fe_scavenge_rate + (Fe_loc - Fe_scavenge_thres1) * fe_max_scale2
    end if

    Fe_scavenge = yps * Fe_loc * Fe_scavenge_rate

    P_iron%prod(k) = (sum(zoo_loss_poc(:)) + sum(zoo_graze_poc(:))) * Qfe_zoo + Fe_scavenge

    do auto_ind = 1, autotroph_cnt
       P_iron%prod(k) = P_iron%prod(k) + Qfe(auto_ind) * &
            (auto_agg(auto_ind) + auto_graze_poc(auto_ind) + auto_loss_poc(auto_ind))
    end do

    end associate
  end subroutine compute_large_detritus

  !***********************************************************************

  subroutine compute_nitrif(PAR_out, PAR_in, KPARdz, NH4_loc, nitrif)

    !-----------------------------------------------------------------------
    !  nitrate & ammonium
    !  nitrification in low light
    !  use exponential decay of PAR across model level to compute taper factor
    !-----------------------------------------------------------------------

    use marbl_parms, only : parm_nitrif_par_lim
    use marbl_parms, only : parm_kappa_nitrif

    real(r8), intent(in)  :: PAR_out
    real(r8), intent(in)  :: PAR_in
    real(r8), intent(in)  :: kPARdz
    real(r8), intent(in)  :: NH4_loc
    real(r8), intent(out) :: nitrif

    if (PAR_out < parm_nitrif_par_lim) then
       nitrif = parm_kappa_nitrif * NH4_loc
       if (PAR_in > parm_nitrif_par_lim) then
          nitrif = nitrif * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
       end if
    else
       nitrif = c0
    end if

  end subroutine compute_nitrif

  !***********************************************************************

  subroutine compute_denitrif(O2_loc, NO3_loc, DOC_remin, POC_remin, other_remin, sed_denitrif, &
       denitrif)

    !-----------------------------------------------------------------------
    !  Compute denitrification under low O2 conditions
    !-----------------------------------------------------------------------

    real(r8), intent(in)  :: O2_loc
    real(r8), intent(in)  :: NO3_loc
    real(r8), intent(in)  :: DOC_remin
    real(r8), intent(in)  :: POC_remin
    real(r8), intent(in)  :: OTHER_REMIN
    real(r8), intent(in)  :: SED_DENITRIF
    real(r8), intent(out) :: denitrif

    real(r8) :: work

    work = ((parm_o2_min + parm_o2_min_delta) - O2_loc) / parm_o2_min_delta
    work = min(max(work, c0), c1)
    work = merge(c0, work, NO3_loc == c0)
    denitrif = work * ((DOC_remin + POC_remin - other_remin) / denitrif_C_N  - sed_denitrif)

  end subroutine compute_denitrif

  !***********************************************************************

  subroutine compute_dtracer_local (auto_cnt, zoo_cnt, auto_meta, zoo_meta, &
       autotroph_secondary_species, &
       zooplankton_secondary_species, &
       dissolved_organic_matter, &
       nitrif, denitrif, sed_denitrif, &
       Fe_scavenge, Fe_scavenge_rate, &
       P_iron_remin, POC_remin, &
       P_SiO2_remin, P_CaCO3_remin, other_remin, &
       restore_local, &
       O2_loc, o2_production, o2_consumption, &
       dtracer)

    use marbl_share_mod, only : autotroph_type
    use marbl_share_mod, only : zooplankton_type

    integer                                      , intent(in)  :: auto_cnt
    integer                                      , intent(in)  :: zoo_cnt
    type(autotroph_type)                         , intent(in)  :: auto_meta(auto_cnt)
    type(zooplankton_type)                       , intent(in)  :: zoo_meta(zoo_cnt)
    type(zooplankton_secondary_species_type) , intent(in)  :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(in)  :: autotroph_secondary_species(auto_cnt)
    type(dissolved_organic_matter_type)   , intent(in)  :: dissolved_organic_matter
    real(r8)                                     , intent(in)  :: nitrif
    real(r8)                                     , intent(in)  :: denitrif
    real(r8)                                     , intent(in)  :: sed_denitrif
    real(r8)                                     , intent(in)  :: Fe_scavenge
    real(r8)                                     , intent(in)  :: Fe_scavenge_rate
    real(r8)                                     , intent(in)  :: P_iron_remin
    real(r8)                                     , intent(in)  :: POC_remin
    real(r8)                                     , intent(in)  :: P_SiO2_remin
    real(r8)                                     , intent(in)  :: P_CaCO3_remin
    real(r8)                                     , intent(in)  :: other_remin
    real(r8)                                     , intent(in)  :: restore_local(ecosys_tracer_cnt)
    real(r8)                                     , intent(in)  :: O2_loc
    real(r8)                                     , intent(out) :: o2_production
    real(r8)                                     , intent(out) :: o2_consumption
    real(r8)                                     , intent(out) :: dtracer(ecosys_tracer_cnt)

    integer  :: auto_ind, zoo_ind, n
    real(r8) :: auto_sum

    associate(                                                                &
         thetaC          => autotroph_secondary_species%thetaC          , & ! local Chl/C ratio (mg Chl/mmol C)
         QCaCO3          => autotroph_secondary_species%QCaCO3          , & ! CaCO3/C ratio (mmol CaCO3/mmol C)
         Qfe             => autotroph_secondary_species%Qfe             , & ! init fe/C ratio (mmolFe/mmolC)
         Qsi             => autotroph_secondary_species%Qsi             , & ! initial Si/C ratio (mmol Si/mmol C)
         NO3_V           => autotroph_secondary_species%NO3_V           , & ! nitrate uptake (mmol NO3/m^3/sec)
         NH4_V           => autotroph_secondary_species%NH4_V           , & ! ammonium uptake (mmol NH4/m^3/sec)
         PO4_V           => autotroph_secondary_species%PO4_V           , & ! PO4 uptake (mmol PO4/m^3/sec)
         DOP_V           => autotroph_secondary_species%DOP_V           , & ! DOP uptake (mmol DOP/m^3/sec)
         photoC          => autotroph_secondary_species%photoC          , & ! C-fixation (mmol C/m^3/sec)
         photoFe         => autotroph_secondary_species%photoFe         , & ! iron uptake
         photoSi         => autotroph_secondary_species%photoSi         , & ! silicon uptake (mmol Si/m^3/sec)
         photoacc        => autotroph_secondary_species%photoacc        , & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
         auto_loss       => autotroph_secondary_species%auto_loss       , & ! autotroph non-grazing mort (mmol C/m^3/sec)
         auto_loss_dic   => autotroph_secondary_species%auto_loss_dic   , & ! auto_loss routed to dic (mmol C/m^3/sec)
         auto_agg        => autotroph_secondary_species%auto_agg        , & ! autotroph aggregation (mmol C/m^3/sec)
         auto_graze      => autotroph_secondary_species%auto_graze      , & ! autotroph grazing rate (mmol C/m^3/sec)
         auto_graze_zoo  => autotroph_secondary_species%auto_graze_zoo  , & ! auto_graze routed to zoo (mmol C/m^3/sec)
         auto_graze_dic  => autotroph_secondary_species%auto_graze_dic  , & ! auto_graze routed to dic (mmol C/m^3/sec)
         CaCO3_PROD      => autotroph_secondary_species%CaCO3_PROD      , & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
         Nfix            => autotroph_secondary_species%Nfix            , & ! total Nitrogen fixation (mmol N/m^3/sec)
         Nexcrete        => autotroph_secondary_species%Nexcrete        , & ! fixed N excretion
         remaining_P_dip => autotroph_secondary_species%remaining_P_dip , & ! remaining_P from mort routed to remin

         f_zoo_detr      => zooplankton_secondary_species%f_zoo_detr    , & ! frac of zoo losses into large detrital pool (non-dim)
         x_graze_zoo     => zooplankton_secondary_species%x_graze_zoo   , & ! {auto, zoo}_graze routed to zoo (mmol C/m^3/sec)
         zoo_graze       => zooplankton_secondary_species%zoo_graze     , & ! zooplankton losses due to grazing (mmol C/m^3/sec)
         zoo_graze_zoo   => zooplankton_secondary_species%zoo_graze_zoo , & ! grazing of zooplankton routed to zoo (mmol C/m^3/sec)
         zoo_graze_dic   => zooplankton_secondary_species%zoo_graze_dic , & ! grazing of zooplankton routed to dic (mmol C/m^3/sec)
         zoo_loss        => zooplankton_secondary_species%zoo_loss      , & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
         zoo_loss_dic    => zooplankton_secondary_species%zoo_loss_dic  , & ! zoo_loss routed to dic (mmol C/m^3/sec)

         DOC_prod        => dissolved_organic_matter%DOC_prod        , & ! production of DOC (mmol C/m^3/sec)
         DOC_remin       => dissolved_organic_matter%DOC_remin       , & ! remineralization of DOC (mmol C/m^3/sec)
         DON_prod        => dissolved_organic_matter%DON_prod        , & ! production of dissolved organic N
         DON_remin       => dissolved_organic_matter%DON_remin       , & ! portion of DON remineralized
         DOFe_prod       => dissolved_organic_matter%DOFe_prod       , & ! produciton of dissolved organic Fe
         DOFe_remin      => dissolved_organic_matter%DOFe_remin      , & ! portion of DOFe remineralized
         DOP_prod        => dissolved_organic_matter%DOP_prod        , & ! production of dissolved organic P
         DOP_remin       => dissolved_organic_matter%DOP_remin       , & ! portion of DOP remineralized
         DONr_remin      => dissolved_organic_matter%DONr_remin      , & ! portion of refractory DON remineralized
         DOPr_remin      => dissolved_organic_matter%DOPr_remin        & ! portion of refractory DOP remineralized
         )

    !-----------------------------------------------------------------------
    !  nitrate & ammonium
    !-----------------------------------------------------------------------

    dtracer(no3_ind) = restore_local(no3_ind) + nitrif - denitrif - sed_denitrif - sum(NO3_V(:))

    dtracer(nh4_ind) = -sum(NH4_V(:)) - nitrif + DON_remin + DONr_remin  &
         + Q * (sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) + sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) &
         + POC_remin * (c1 - DONrefract) )

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%Nfixer) then
          dtracer(nh4_ind) = dtracer(nh4_ind) + Nexcrete(auto_ind)
       end if
    end do

    !-----------------------------------------------------------------------
    !  dissolved iron
    !-----------------------------------------------------------------------

    dtracer(fe_ind) = P_iron_remin + (Qfe_zoo * ( sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) )) &
         + DOFe_remin - sum(photofe(:)) - Fe_scavenge

    do auto_ind = 1, autotroph_cnt
       dtracer(fe_ind) = dtracer(fe_ind) &
            + (Qfe(auto_ind) * (auto_loss_dic(auto_ind) + auto_graze_dic(auto_ind))) &
            + auto_graze_zoo(auto_ind) * (Qfe(auto_ind) - Qfe_zoo)
    end do

    !-----------------------------------------------------------------------
    !  dissolved SiO3
    !-----------------------------------------------------------------------

    dtracer(sio3_ind) = restore_local(sio3_ind) + P_SiO2_remin

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%Si_ind > 0) then
          dtracer(sio3_ind) = dtracer(sio3_ind) &
               - photoSi(auto_ind) + Qsi(auto_ind) * (f_graze_si_remin * auto_graze(auto_ind) &
               + (c1 - auto_meta(auto_ind)%loss_poc) * auto_loss(auto_ind))
       endif
    end do

    !-----------------------------------------------------------------------
    !  phosphate
    !-----------------------------------------------------------------------

    dtracer(po4_ind) = restore_local(po4_ind) + DOP_remin + DOPr_remin - sum(PO4_V(:)) &
         + Qp_zoo_pom * ( (c1 - DOPrefract) * POC_remin + sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) )

    do auto_ind = 1, autotroph_cnt
       if (auto_meta(auto_ind)%Qp == Qp_zoo_pom) then
          dtracer(po4_ind) = dtracer(po4_ind) &
               + auto_meta(auto_ind)%Qp * (auto_loss_dic(auto_ind) + auto_graze_dic(auto_ind))
       else
          dtracer(po4_ind) = dtracer(po4_ind) &
               + remaining_P_dip(auto_ind)
       endif
    end do

    !-----------------------------------------------------------------------
    !  zoo Carbon
    !-----------------------------------------------------------------------
    do zoo_ind = 1, zoo_cnt
       n = zoo_meta(zoo_ind)%C_ind
       dtracer(n) = x_graze_zoo(zoo_ind) - zoo_graze(zoo_ind) - zoo_loss(zoo_ind)
    end do

    !-----------------------------------------------------------------------
    !  autotroph Carbon
    !  autotroph Chlorophyll
    !  autotroph Fe
    !  autotroph Si
    !  autotroph CaCO3
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt
       auto_sum = auto_graze(auto_ind) + auto_loss(auto_ind) + auto_agg(auto_ind)

       n = autotrophs(auto_ind)%C_ind
       dtracer(n) = photoC(auto_ind) - auto_sum

       n = autotrophs(auto_ind)%Chl_ind
       dtracer(n) = photoacc(auto_ind) - thetaC(auto_ind) * auto_sum

       n = autotrophs(auto_ind)%Fe_ind
       dtracer(n) =  photoFe(auto_ind) - Qfe(auto_ind) * auto_sum

       n = autotrophs(auto_ind)%Si_ind
       if (n > 0) then
          dtracer(n) =  photoSi(auto_ind) - Qsi(auto_ind) * auto_sum
       endif

       n = autotrophs(auto_ind)%CaCO3_ind
       if (n > 0) then
          dtracer(n) = CaCO3_PROD(auto_ind) - QCaCO3(auto_ind) * auto_sum
       endif
    end do


    !-----------------------------------------------------------------------
    !  dissolved organic Matter
    !  from sinking remin small fraction to refractory pool
    !-----------------------------------------------------------------------

    dtracer(doc_ind) = DOC_prod - DOC_remin

    dtracer(don_ind) = (DON_prod * (c1 - DONrefract)) - DON_remin

    dtracer(donr_ind) = (DON_prod * DONrefract) - DONr_remin + (POC_remin * DONrefract * Q)

    dtracer(dop_ind) = (DOP_prod * (c1 - DOPrefract)) - DOP_remin - sum(DOP_V(:))

    dtracer(dopr_ind) = (DOP_prod * DOPrefract) - DOPr_remin + (POC_remin * DOPrefract * Qp_zoo_pom)

    dtracer(dofe_ind) = DOFe_prod - DOFe_remin


    !-----------------------------------------------------------------------
    !  dissolved inorganic Carbon
    !-----------------------------------------------------------------------

    dtracer(dic_ind) = &
         sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) - sum(photoC(:)) &
            + DOC_remin + POC_remin + sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) + P_CaCO3_remin

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%CaCO3_ind > 0) then
          dtracer(dic_ind) = dtracer(dic_ind) &
               + f_graze_CaCO3_REMIN * auto_graze(auto_ind) * QCaCO3(auto_ind) - CaCO3_PROD(auto_ind)
       end if
    end do

    dtracer(dic_alt_co2_ind) = dtracer(dic_ind)


    !-----------------------------------------------------------------------
    !  alkalinity
    !-----------------------------------------------------------------------

    dtracer(alk_ind) = -dtracer(no3_ind) + dtracer(nh4_ind) + c2 * P_CaCO3_remin

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%CaCO3_ind > 0) then
          dtracer(alk_ind) = dtracer(alk_ind) &
               + c2 * (f_graze_CaCO3_REMIN * auto_graze(auto_ind) * QCaCO3(auto_ind) - CaCO3_PROD(auto_ind))
       end if
    end do

    !-----------------------------------------------------------------------
    !  oxygen
    !-----------------------------------------------------------------------

    o2_production = c0
    do auto_ind = 1, auto_cnt
       if (.not. auto_meta(auto_ind)%Nfixer) then
          if (photoC(auto_ind) > c0) then
             o2_production = o2_production + photoC(auto_ind) * &
                  ((NO3_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind))) / parm_Red_D_C_O2 + &
                   (NH4_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind))) / parm_Remin_D_C_O2)
          end if
       else
          if (photoC(auto_ind) > c0) then
             o2_production = o2_production + photoC(auto_ind) * &
                  ((NO3_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Red_D_C_O2 + &
                   (NH4_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Remin_D_C_O2 + &
                   (Nfix(auto_ind)  / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Red_D_C_O2_diaz)
          end if
       endif
    end do

    o2_consumption = (O2_loc - parm_o2_min) / parm_o2_min_delta
    o2_consumption = min(max(o2_consumption, c0), c1)
    o2_consumption = o2_consumption * ( (POC_remin + DOC_remin - (sed_denitrif * denitrif_C_N) - other_remin &
         + sum(zoo_loss_dic(:))  +  sum(zoo_graze_dic(:)) + sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) ) &
         / parm_Remin_D_C_O2 + (c2 * nitrif))

    dtracer(o2_ind) = o2_production - o2_consumption

    end associate
  end subroutine compute_dtracer_local

  !-----------------------------------------------------------------------

  subroutine export_interior_shared_variables (&
       tracer_local, &
       carbonate, &
       dissolved_organic_matter, &
       QA_dust_def, &
       marbl_interior_share)

    use marbl_share_mod, only : marbl_interior_share_type

    real(r8), intent(in) :: tracer_local(ecosys_tracer_cnt)
    type(carbonate_type), intent(in) :: carbonate
    type(dissolved_organic_matter_type), intent(in) :: dissolved_organic_matter
    real(r8), intent(in) :: QA_dust_def
    type(marbl_interior_share_type), intent(inout) :: marbl_interior_share

    associate( &
         share => marbl_interior_share &
         )

    share%QA_dust_def = QA_dust_def
    share%DIC_loc_fields = tracer_local(DIC_ind)
    share%DOC_loc_fields = tracer_local(DOC_ind)
    share%O2_loc_fields = tracer_local(O2_ind)
    share%NO3_loc_fields = tracer_local(NO3_ind)


    share%CO3_fields = carbonate%CO3
    share%HCO3_fields = carbonate%HCO3
    share%H2CO3_fields = carbonate%H2CO3

    share%DOC_remin_fields = dissolved_organic_matter%DOC_remin

    end associate
  end subroutine export_interior_shared_variables

  !-----------------------------------------------------------------------

  subroutine export_zooplankton_shared_variables (&
       zoo_cnt, &
       zooplankton_local, &
       zooplankton_secondary_species, &
       marbl_zooplankton_share)

    use marbl_share_mod, only : marbl_zooplankton_share_type

    integer(int_kind) :: zoo_cnt
    type(zooplankton_local_type), intent(in) :: zooplankton_local(zoo_cnt)
    type(zooplankton_secondary_species_type), intent(in) :: zooplankton_secondary_species(zoo_cnt)
    type(marbl_zooplankton_share_type), intent(inout) :: marbl_zooplankton_share(zoo_cnt)

    integer(int_kind) :: n

    associate( &
         share => marbl_zooplankton_share(:) &
         )

    do n = 1, zoo_cnt
       share(n)%zooC_loc_fields     = zooplankton_local(n)%C
       share(n)%zoo_loss_fields     = zooplankton_secondary_species(n)%zoo_loss
       share(n)%zoo_loss_poc_fields = zooplankton_secondary_species(n)%zoo_loss_poc
       share(n)%zoo_loss_doc_fields = zooplankton_secondary_species(n)%zoo_loss_doc
       share(n)%zoo_loss_dic_fields = zooplankton_secondary_species(n)%zoo_loss_dic
    end do

    end associate
  end subroutine export_zooplankton_shared_variables

  !-----------------------------------------------------------------------

  subroutine export_autotroph_shared_variables (&
       auto_cnt, &
       autotroph_local, &
       autotroph_secondary_species, &
       marbl_autotroph_share)

    use marbl_share_mod, only : marbl_autotroph_share_type

    integer(int_kind) :: auto_cnt
    type(autotroph_local_type), intent(in) :: autotroph_local(auto_cnt)
    type(autotroph_secondary_species_type), intent(in) :: autotroph_secondary_species(auto_cnt)
    type(marbl_autotroph_share_type), intent(inout) :: marbl_autotroph_share(auto_cnt)

    integer(int_kind) :: n

    associate( &
         share => marbl_autotroph_share(:) &
         )

    do n = 1, auto_cnt
       share(n)%autotrophChl_loc_fields = autotroph_local(n)%Chl
       share(n)%autotrophC_loc_fields = autotroph_local(n)%C
       share(n)%autotrophFe_loc_fields = autotroph_local(n)%Fe

       if (autotrophs(n)%Si_ind > 0) then
          share(n)%autotrophSi_loc_fields = autotroph_local(n)%Si
       end if

       if (autotrophs(n)%CaCO3_ind > 0) then
          share(n)%autotrophCaCO3_loc_fields = autotroph_local(n)%CaCO3
       end if

       share(n)%QCaCO3_fields         = autotroph_secondary_species(n)%QCaCO3
       share(n)%auto_graze_fields     = autotroph_secondary_species(n)%auto_graze
       share(n)%auto_graze_zoo_fields = autotroph_secondary_species(n)%auto_graze_zoo
       share(n)%auto_graze_poc_fields = autotroph_secondary_species(n)%auto_graze_poc
       share(n)%auto_graze_doc_fields = autotroph_secondary_species(n)%auto_graze_doc
       share(n)%auto_graze_dic_fields = autotroph_secondary_species(n)%auto_graze_dic
       share(n)%auto_loss_fields      = autotroph_secondary_species(n)%auto_loss
       share(n)%auto_loss_poc_fields  = autotroph_secondary_species(n)%auto_loss_poc
       share(n)%auto_loss_doc_fields  = autotroph_secondary_species(n)%auto_loss_doc
       share(n)%auto_loss_dic_fields  = autotroph_secondary_species(n)%auto_loss_dic
       share(n)%auto_agg_fields       = autotroph_secondary_species(n)%auto_agg
       share(n)%photoC_fields         = autotroph_secondary_species(n)%photoC
       share(n)%CaCO3_PROD_fields     = autotroph_secondary_species(n)%CaCO3_PROD
       share(n)%PCphoto_fields        = autotroph_secondary_species(n)%PCphoto
    end do
    end associate
  end subroutine export_autotroph_shared_variables

end module ecosys_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
