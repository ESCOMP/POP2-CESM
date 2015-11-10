! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_mod

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

  !-----------------------------------------------------------------------
  !  variables/subroutines/function used from other modules
  !  The following are used extensively in this ecosys, so are used at
  !  the module level. The use statements for variables that are only needed
  !  locally are located at the module subprogram level.
  !-----------------------------------------------------------------------

  ! !USES:

  use marbl_kinds_mod, only : log_kind
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : char_len

  use POP_ErrorMod         , only : POP_Success  
  use constants            , only : c0
  use constants            , only : c1
  use constants            , only : c2
  use constants            , only : c10
  use constants            , only : mpercm
  use constants            , only : blank_fmt
  use constants            , only : delim_fmt
  use constants            , only : ndelim_fmt
  use constants            , only : xkw_coeff
  use communicate          , only : master_task
  use communicate          , only : my_task
  use blocks               , only : nx_block
  use blocks               , only : ny_block
  use domain_size          , only : km
  use domain_size          , only : max_blocks_clinic
  use domain               , only : nblocks_clinic
  use exit_mod             , only : exit_POP
  use exit_mod             , only : sigAbort
  use grid                 , only : partial_bottom_cells
  use grid                 , only : zt
  use grid                 , only : zw
  use io_types             , only : stdout
  use io_tools             , only : document
  use passive_tracer_tools , only : tracer_read
  use time_management      , only : freq_opt_never
  use time_management      , only : freq_opt_nmonth
  use time_management      , only : freq_opt_nyear
  use time_management      , only : check_time_flag
  use time_management      , only : init_time_flag
  use time_management      , only : eval_time_flag
  use ecosys_constants     , only : ecosys_tracer_cnt

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

  !
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
  use marbl_parms, only : po4_ind 
  use marbl_parms, only : no3_ind         
  use marbl_parms, only : sio3_ind        
  use marbl_parms, only : nh4_ind         
  use marbl_parms, only : fe_ind          
  use marbl_parms, only : o2_ind          
  use marbl_parms, only : dic_ind         
  use marbl_parms, only : dic_alt_co2_ind 
  use marbl_parms, only : alk_ind         
  use marbl_parms, only : doc_ind         
  use marbl_parms, only : don_ind         
  use marbl_parms, only : dofe_ind        
  use marbl_parms, only : dop_ind         
  use marbl_parms, only : dopr_ind        
  use marbl_parms, only : donr_ind        

  use marbl_share_mod, only : autotrophs, autotroph_cnt
  use marbl_share_mod, only : zooplankton, zooplankton_cnt
  use marbl_share_mod, only : grazing, grazer_prey_cnt

  use marbl_interface_types, only : carbonate_type
  use marbl_interface_types, only : zooplankton_secondary_species_type
  use marbl_interface_types, only : autotroph_secondary_species_type
  use marbl_interface_types, only : photosynthetically_available_radiation_type
  use marbl_interface_types, only : dissolved_organic_matter_type

  implicit none
  private

  !-----------------------------------------------------------------------
  !  public/private member procedure declarations
  !-----------------------------------------------------------------------

  ! initialization routines
  public  :: marbl_ecosys_init_nml
  public  :: marbl_ecosys_init_tracer_metadata
  private :: marbl_ecosys_init_non_autotroph_tracer_metadata
  private :: marbl_ecosys_init_forcing_metadata
  public  :: marbl_ecosys_init_tavg

  ! set_interior routines
  public  :: marbl_ecosys_set_interior
  private :: marbl_init_particulate_terms
  private :: marbl_update_particulate_terms_from_prior_level      
  private :: marbl_update_sinking_particle_from_prior_level       
  private :: marbl_compute_particulate_terms                      
  private :: marbl_check_ecosys_tracer_count_consistency          
  private :: marbl_initialize_zooplankton_tracer_metadata         
  private :: marbl_initialize_autotroph_tracer_metadata           
  private :: marbl_setup_local_tracers                            
  private :: marbl_setup_local_zooplankton                        
  private :: marbl_setup_local_autotrophs                         
  private :: marbl_autotroph_consistency_check                    
  private :: marbl_compute_autotroph_elemental_ratios             
  private :: marbl_compute_photosynthetically_available_radiation 
  private :: marbl_compute_carbonate_chemistry                    
  private :: marbl_compute_function_scaling                       
  private :: marbl_compute_Pprime                                 
  private :: marbl_compute_Zprime
       
  ! set surface co2 flux
  public  :: marbl_ecosys_set_sflux

  ! tavg
  public  :: marbl_ecosys_tavg_forcing

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

  type(tracer_read) :: &
       gas_flux_fice,       & ! ice fraction for gas fluxes
       gas_flux_ws,         & ! wind speed for gas fluxes
       gas_flux_ap            ! atmospheric pressure for gas fluxes

  !-----------------------------------------------------------------------
  !  restoring climatologies for nutrients
  !-----------------------------------------------------------------------

  character(char_len) :: &
       nutr_rest_file               ! file containing nutrient fields

  !maltrud variable restoring
  logical (log_kind) :: &
       lnutr_variable_restore       ! geographically varying nutrient restoring

  character(char_len) :: &
       nutr_variable_rest_file,   & ! file containing variable restoring info
       nutr_variable_rest_file_fmt  ! format of file containing variable restoring info

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
  !  iron patch fertilization
  !-----------------------------------------------------------------------

  logical (log_kind)  :: liron_patch               ! flag for iron patch fertilization
  character(char_len) :: iron_patch_flux_filename  ! file containing name of iron patch file
  integer (int_kind)  :: iron_patch_month          !  integer month to add patch flux

  !-----------------------------------------------------------------------
  !  named field indices
  !-----------------------------------------------------------------------

  integer (int_kind) :: &
!       sflux_co2_nf_ind   = 0,    & ! air-sea co2 gas flux TEMPORARY
       atm_co2_nf_ind     = 0       ! atmospheric co2

  !-----------------------------------------------------------------------

  real (r8), parameter :: &
       phlo_surf_init = 7.0_r8, & ! low bound for surface ph for no prev soln
       phhi_surf_init = 9.0_r8, & ! high bound for surface ph for no prev soln
       phlo_3d_init = 6.0_r8,   & ! low bound for subsurface ph for no prev soln
       phhi_3d_init = 9.0_r8,   & ! high bound for subsurface ph for no prev soln
       del_ph = 0.20_r8           ! delta-ph for prev soln

  !-----------------------------------------------------------------------

  logical (log_kind)  :: lecovars_full_depth_tavg ! should ecosystem vars be written full depth

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine marbl_ecosys_init_nml(nl_buffer, marbl_status)

    ! !DESCRIPTION:
    !  Initialize ecosys tracer module. This involves setting metadata, reading
    !  the module namelist, setting initial conditions, setting up forcing,
    !  and defining additional tavg variables.
    !
    use marbl_interface_constants , only : marbl_nl_buffer_size
    use marbl_interface_constants , only : marbl_status_ok
    use marbl_interface_constants , only : marbl_status_could_not_read_namelist
    use marbl_interface_types     , only : marbl_status_type
    use marbl_interface_types     , only : marbl_tracer_read_type
    use marbl_share_mod           , only : surf_avg_dic_const, surf_avg_alk_const
    use marbl_share_mod           , only : use_nml_surf_vals         
    use marbl_share_mod           , only : init_ecosys_option        
    use marbl_share_mod           , only : init_ecosys_init_file     
    use marbl_share_mod           , only : init_ecosys_init_file_fmt 
    use marbl_share_mod           , only : tracer_init_ext
    use marbl_share_mod           , only : ndep_data_type 
    use marbl_share_mod           , only : ndep_shr_stream_year_first
    use marbl_share_mod           , only : ndep_shr_stream_year_last
    use marbl_share_mod           , only : ndep_shr_stream_year_align
    use marbl_share_mod           , only : ndep_shr_stream_file
    use marbl_share_mod           , only : ndep_shr_stream_scale_factor
    use marbl_share_mod           , only : lflux_gas_co2
    use marbl_share_mod           , only : lflux_gas_o2
    use marbl_share_mod           , only : gas_flux_forcing_iopt_drv
    use marbl_share_mod           , only : gas_flux_forcing_iopt_file
    use marbl_share_mod           , only : gas_flux_forcing_iopt
    use marbl_share_mod           , only : gas_flux_forcing_file
    use marbl_share_mod           , only : atm_co2_const
    use marbl_share_mod           , only : atm_alt_co2_const
    use marbl_share_mod           , only : atm_co2_iopt
    use marbl_share_mod           , only : atm_co2_iopt_const
    use marbl_share_mod           , only : atm_co2_iopt_drv_prog
    use marbl_share_mod           , only : atm_co2_iopt_drv_diag
    use marbl_share_mod           , only : atm_alt_co2_iopt
    use marbl_share_mod           , only : dust_flux        
    use marbl_share_mod           , only : iron_flux        
    use marbl_share_mod           , only : fice_file        
    use marbl_share_mod           , only : xkw_file         
    use marbl_share_mod           , only : ap_file          
    use marbl_share_mod           , only : nox_flux_monthly 
    use marbl_share_mod           , only : nhy_flux_monthly 
    use marbl_share_mod           , only : din_riv_flux     
    use marbl_share_mod           , only : dip_riv_flux     
    use marbl_share_mod           , only : don_riv_flux     
    use marbl_share_mod           , only : dop_riv_flux     
    use marbl_share_mod           , only : dsi_riv_flux     
    use marbl_share_mod           , only : dfe_riv_flux     
    use marbl_share_mod           , only : dic_riv_flux     
    use marbl_share_mod           , only : alk_riv_flux     
    use marbl_share_mod           , only : doc_riv_flux     
    use marbl_share_mod           , only : liron_patch  
    use marbl_share_mod           , only : iron_patch_flux_filename  
    use marbl_share_mod           , only : iron_patch_month  
    use marbl_share_mod           , only : ecosys_qsw_distrb_const
    use marbl_share_mod           , only : fesedflux_input 

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use marbl_share_mod           , only : comp_surf_avg_flag 

    use ecosys_constants          , only : ecosys_tracer_cnt

    implicit none

    ! !INPUT PARAMETERS:
    character(marbl_nl_buffer_size) , intent(in) :: nl_buffer

    ! !OUTPUT PARAMETERS:
    type(marbl_status_type) , intent(out) :: marbl_status

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_init_nml'

    integer (int_kind)  :: n                        ! index for looping over tracers
    character(char_len) :: comp_surf_avg_freq_opt   ! choice for freq of comp_surf_avg
    character(char_len) :: gas_flux_forcing_opt     ! option for forcing gas fluxes
    character(char_len) :: atm_co2_opt              ! option for atmospheric co2 concentration
    character(char_len) :: atm_alt_co2_opt          ! option for atmospheric alternative CO2
    type(tracer_read)   :: dust_flux_input          ! namelist input for dust_flux
    type(tracer_read)   :: iron_flux_input          ! namelist input for iron_flux
    type(tracer_read)   :: nox_flux_monthly_input   ! namelist input for nox_flux_monthly
    type(tracer_read)   :: nhy_flux_monthly_input   ! namelist input for nhy_flux_monthly
    type(tracer_read)   :: din_riv_flux_input       ! namelist input for din_riv_flux
    type(tracer_read)   :: dip_riv_flux_input       ! namelist input for dip_riv_flux
    type(tracer_read)   :: don_riv_flux_input       ! namelist input for don_riv_flux
    type(tracer_read)   :: dop_riv_flux_input       ! namelist input for dop_riv_flux
    type(tracer_read)   :: dsi_riv_flux_input       ! namelist input for dsi_riv_flux
    type(tracer_read)   :: dfe_riv_flux_input       ! namelist input for dfe_riv_flux
    type(tracer_read)   :: dic_riv_flux_input       ! namelist input for dic_riv_flux
    type(tracer_read)   :: alk_riv_flux_input       ! namelist input for alk_riv_flux
    type(tracer_read)   :: doc_riv_flux_input       ! namelist input for doc_riv_flux
    integer (int_kind)  :: nml_error                ! namelist i/o error flag
    integer (int_kind)  :: zoo_ind                  ! zooplankton functional group index
    integer (int_kind)  :: comp_surf_avg_freq_iopt  ! choice for freq of comp_surf_avg
    integer (int_kind)  :: comp_surf_avg_freq       ! choice for freq of comp_surf_avg

    !-----------------------------------------------------------------------
    !  values to be used when comp_surf_avg_freq_opt==never
    !-----------------------------------------------------------------------

    namelist /ecosys_nml/                                                 &
         init_ecosys_option, init_ecosys_init_file, tracer_init_ext,      &
         init_ecosys_init_file_fmt,                                       &
         dust_flux_input, iron_flux_input, fesedflux_input,               &
         ndep_data_type, nox_flux_monthly_input, nhy_flux_monthly_input,  &
         ndep_shr_stream_year_first, ndep_shr_stream_year_last,           &
         ndep_shr_stream_year_align, ndep_shr_stream_file,                &
         ndep_shr_stream_scale_factor,                                    &
         din_riv_flux_input, dip_riv_flux_input, don_riv_flux_input,      &
         dop_riv_flux_input, dsi_riv_flux_input, dfe_riv_flux_input,      &
         dic_riv_flux_input, alk_riv_flux_input, doc_riv_flux_input,      &
         gas_flux_forcing_opt, gas_flux_forcing_file,                     &
         gas_flux_fice, gas_flux_ws, gas_flux_ap,                         &
         nutr_rest_file,                                                  &
         comp_surf_avg_freq_opt, comp_surf_avg_freq,                      &
         use_nml_surf_vals, surf_avg_dic_const, surf_avg_alk_const,       &
         ecosys_qsw_distrb_const,                                         &
         lsource_sink, lflux_gas_o2, lflux_gas_co2, locmip_k1_k2_bug_fix, &
         lnutr_variable_restore, nutr_variable_rest_file,                 &
         nutr_variable_rest_file_fmt, atm_co2_opt, atm_co2_const,         &
         atm_alt_co2_opt, atm_alt_co2_const,                              &
         liron_patch, iron_patch_flux_filename, iron_patch_month,         &
         lecovars_full_depth_tavg

    marbl_status%status = marbl_status_ok
    marbl_status%message = ''

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

    !-----------------------------------------------------------------------
    ! read the namelist buffer on every processor
    !-----------------------------------------------------------------------

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

    fice_file%input        = gas_flux_fice
    xkw_file%input         = gas_flux_ws
    ap_file%input          = gas_flux_ap
    dust_flux%input        = dust_flux_input
    iron_flux%input        = iron_flux_input
    nox_flux_monthly%input = nox_flux_monthly_input
    nhy_flux_monthly%input = nhy_flux_monthly_input
    din_riv_flux%input     = din_riv_flux_input
    dip_riv_flux%input     = dip_riv_flux_input
    don_riv_flux%input     = don_riv_flux_input
    dop_riv_flux%input     = dop_riv_flux_input
    dsi_riv_flux%input     = dsi_riv_flux_input
    dfe_riv_flux%input     = dfe_riv_flux_input
    dic_riv_flux%input     = dic_riv_flux_input
    alk_riv_flux%input     = alk_riv_flux_input
    doc_riv_flux%input     = doc_riv_flux_input

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
    !  read ecosys_parms_nml namelist
    !-----------------------------------------------------------------------

    call marbl_params_init(nl_buffer, marbl_status)
    if (marbl_status%status /= marbl_status_ok) then
       return
    end if
    ! FIXME(bja, 2015-01) need to have a flag if params should be printed!
    ! if (stdout > 0) then
    !    call marbl_params_print(stdout)
    ! end if

  end subroutine marbl_ecosys_init_nml

  !*****************************************************************************
  
  subroutine marbl_ecosys_init_tracer_metadata(tracer_d_module)

    ! !DESCRIPTION:
    !  Set tracer and forcing metadata

    use marbl_interface_constants , only : marbl_status_ok
    use marbl_interface_constants , only : marbl_status_could_not_read_namelist
    use marbl_interface_types     , only : marbl_status_type
    use marbl_interface_types     , only : marbl_tracer_metadata_type

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type (marbl_tracer_metadata_type), intent(inout) :: tracer_d_module(:)   ! descriptors for each tracer

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_ecosys_init_tracer_metadata'

    integer (int_kind) :: non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers
    integer (int_kind) :: n        ! index for looping over tracers
    integer (int_kind) :: zoo_ind  ! zooplankton functional group index
    integer (int_kind) :: auto_ind ! autotroph functional group index

    !-----------------------------------------------------------------------
    ! initialize tracer metatdata
    !-----------------------------------------------------------------------

    call marbl_ecosys_init_forcing_metadata()

    call marbl_ecosys_init_non_autotroph_tracer_metadata(tracer_d_module, non_living_biomass_ecosys_tracer_cnt)

    call marbl_check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt)

    call marbl_initialize_zooplankton_tracer_metadata(tracer_d_module, non_living_biomass_ecosys_tracer_cnt, n)

    call marbl_initialize_autotroph_tracer_metadata(tracer_d_module, n)

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

  end subroutine marbl_ecosys_init_tracer_metadata

  !***********************************************************************

  subroutine marbl_ecosys_init_tavg

    ! !DESCRIPTION:
    !  call define_tavg_field for nonstandard tavg fields

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: buf_len  ! how many surface flux fields are stored in ECO_SFLUX_TAVG

    !-----------------------------------------------------------------------
    !  2D fields related to surface fluxes
    !-----------------------------------------------------------------------

    buf_len = 0
    buf_len = buf_len+1;     buf_ind_ECOSYS_IFRAC = buf_len
    buf_len = buf_len+1;     buf_ind_ECOSYS_XKW = buf_len
    buf_len = buf_len+1;     buf_ind_ECOSYS_ATM_PRESS = buf_len
    buf_len = buf_len+1;     buf_ind_PV_O2 = buf_len
    buf_len = buf_len+1;     buf_ind_SCHMIDT_O2 = buf_len
    buf_len = buf_len+1;     buf_ind_O2SAT = buf_len
    buf_len = buf_len+1;     buf_ind_CO2STAR = buf_len
    buf_len = buf_len+1;     buf_ind_DCO2STAR = buf_len
    buf_len = buf_len+1;     buf_ind_pCO2SURF = buf_len
    buf_len = buf_len+1;     buf_ind_DpCO2 = buf_len
    buf_len = buf_len+1;     buf_ind_PV_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_SCHMIDT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_DIC_GAS_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_PH = buf_len
    buf_len = buf_len+1;     buf_ind_ATM_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_CO2STAR_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_DCO2STAR_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_pCO2SURF_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_DpCO2_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_DIC_GAS_FLUX_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_PH_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_ATM_ALT_CO2 = buf_len
    buf_len = buf_len+1;     buf_ind_IRON_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_NOx_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_DIN_RIV_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_DFE_RIV_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_DIC_RIV_FLUX = buf_len
    buf_len = buf_len+1;     buf_ind_ALK_RIV_FLUX = buf_len

    allocate(ECO_SFLUX_TAVG(nx_block, ny_block, buf_len, max_blocks_clinic))
    ECO_SFLUX_TAVG(:,:,:,:) = c0

  end subroutine marbl_ecosys_init_tavg

  !***********************************************************************

  subroutine marbl_ecosys_set_interior( &
       lexport_shared_vars,             &
       domain,                          &
       gcm_state,                       &
       marbl_diagnostics,               &
       restore_local,                   &
       marbl_interior_share,            &
       marbl_zooplankton_share,         &
       marbl_autotroph_share,           &
       marbl_particulate_share,         &
       tracer_module,                   &
       dtracer,                         &
       fesedflux, &
       dust_flux_in, PAR_out, &
       ph_prev_3d, ph_prev_alt_co2_3d)
    
    ! !DESCRIPTION:
    !  Compute time derivatives for ecosystem state variables

    use marbl_interface_types , only : marbl_diagnostics_type
    use marbl_interface_types , only : marbl_column_domain_type
    use marbl_interface_types , only : marbl_gcm_state_type
    use marbl_share_mod       , only : marbl_interior_share_type
    use marbl_share_mod       , only : marbl_autotroph_share_type
    use marbl_share_mod       , only : marbl_zooplankton_share_type
    use marbl_share_mod       , only : marbl_particulate_share_type

    logical (log_kind)                 , intent(in)    :: lexport_shared_vars                  ! flag to save shared_vars or not
    real(r8)                           , intent(in)    :: fesedflux(:)
    real(r8)                           , intent(in)    :: restore_local(ecosys_tracer_cnt, km) ! local restoring terms for nutrients (mmol ./m^3/sec)
    real(r8)                           , intent(in)    :: tracer_module(ecosys_tracer_cnt, km) ! tracer values
    type(marbl_column_domain_type)     , intent(inout) :: domain                               ! FIXME(bja, 2015-08) domain will become intent(in) once the loops are moved!
    type(marbl_gcm_state_type)         , intent(inout) :: gcm_state
    real (r8)                          , intent(out)   :: dtracer(ecosys_tracer_cnt, km)       ! computed source/sink terms
    type(marbl_diagnostics_type)       , intent(inout) :: marbl_diagnostics
    type(marbl_interior_share_type)    , intent(out)   :: marbl_interior_share(km)
    type(marbl_zooplankton_share_type) , intent(out)   :: marbl_zooplankton_share(zooplankton_cnt, km)
    type(marbl_autotroph_share_type)   , intent(out)   :: marbl_autotroph_share(autotroph_cnt, km)
    type(marbl_particulate_share_type) , intent(out)   :: marbl_particulate_share
    real(r8)                           , intent(in)    :: dust_flux_in
    real(r8)                           , intent(inout) :: PAR_out
    real(r8)                           , intent(inout) :: ph_prev_3d(km)        
    real(r8)                           , intent(inout) :: ph_prev_alt_co2_3d(km)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_ecosys_set_interior'

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

    real(r8) :: QA_dust_def(km)

    type(zooplankton_local_type) :: zooplankton_local(zooplankton_cnt, km)
    type(autotroph_local_type) :: autotroph_local(autotroph_cnt, km)

    type(autotroph_secondary_species_type)   :: autotroph_secondary_species(autotroph_cnt, km)
    type(zooplankton_secondary_species_type) :: zooplankton_secondary_species(zooplankton_cnt, km)

    type(photosynthetically_available_radiation_type) :: PAR(km)

    type(dissolved_organic_matter_type) :: dissolved_organic_matter(km)

    type(carbonate_type) :: carbonate(km)

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
         POC     => marbl_particulate_share%POC,     &
         P_CaCO3 => marbl_particulate_share%P_CaCO3, &
         P_SiO2  => marbl_particulate_share%P_SiO2,  &
         dust    => marbl_particulate_share%dust,    &
         P_iron  => marbl_particulate_share%P_iron   &
         )

          do k = 1, domain%km
             !write(*, *) 'set_interior loop: ', k, i, c
             call marbl_setup_local_tracers(k, domain%land_mask, domain%kmt, &
                  tracer_module(:, k), tracer_local(:, k))

             call marbl_setup_local_zooplankton(k, domain%land_mask, domain%kmt, &
                  tracer_module(:, k), zooplankton_cnt, zooplankton, zooplankton_local(:, k))

             call marbl_setup_local_autotrophs(k, domain%land_mask, domain%kmt, tracer_module(:, k), &
                  autotroph_cnt, autotrophs, autotroph_local(:, k))

          enddo
          call marbl_init_particulate_terms(1, POC, P_CaCO3, P_SiO2, &
               dust, P_iron, &
               QA_dust_def(:), dust_flux_in)

          !FIXME (mvertens, 2015-11), new marbl timers need to be implemented to turn on timers here
          ! around this subroutine call
          call marbl_compute_carbonate_chemistry(domain%km,  &
               domain%land_mask, domain%kmt, &
               gcm_state%temperature(:), gcm_state%salinity(:), &
               tracer_local(:, :), &
               carbonate(:), &
               ph_prev_3d(:), ph_prev_alt_co2_3d(:), &
               zsat_calcite(:), zsat_aragonite(:), &
               co3_calc_anom(:), co3_arag_anom(:))

          do k = 1, domain%km

             call marbl_autotroph_consistency_check(autotroph_cnt, autotrophs, autotroph_local(:, k))

             call marbl_compute_autotroph_elemental_ratios(k, autotroph_cnt, autotrophs, autotroph_local(:, k), &
                  tracer_local(:, k), autotroph_secondary_species(:, k))

             call marbl_compute_photosynthetically_available_radiation(k, autotroph_cnt, &
                  autotroph_local(:, k), &
                  domain%land_mask, domain%kmt, domain%dzt(k), domain%dz(k), &
                  PAR_out, PAR(k))

             call marbl_compute_function_scaling(gcm_state%temperature(k), Tfunc(k))

             call marbl_compute_Pprime(k, autotroph_cnt, autotrophs, autotroph_local(:, k), &
                  gcm_state%temperature(k), autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_uptake(autotroph_cnt, autotrophs, &
                  tracer_local(:, k), &
                  autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_photosynthesis(autotroph_cnt, autotrophs, &
                  autotroph_local(:, k), gcm_state%temperature(k), Tfunc(k), &
                  PAR(k)%avg, autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_phyto_diatoms (autotroph_cnt, autotrophs, &
                  autotroph_local(:, k), PAR(k)%avg, autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_calcification(autotroph_cnt, autotrophs, &
                  autotroph_local(:, k),  gcm_state%temperature(k), autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_nfixation(autotroph_cnt, autotrophs, &
                  autotroph_secondary_species(:, k))

             call marbl_compute_autotroph_loss(autotroph_cnt, autotrophs, &
                  Tfunc(k), autotroph_secondary_species(:, k))

             call marbl_compute_Zprime(k, zooplankton_cnt, zooplankton, zooplankton_local(:, k)%C, &
                  Tfunc(k), zooplankton_secondary_species(:, k))

             call marbl_compute_grazing (autotroph_cnt, zooplankton_cnt, grazer_prey_cnt, autotrophs, &
                  Tfunc(k), zooplankton_local(:, k), &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

             call marbl_compute_routing (autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

             call marbl_compute_dissolved_organic_matter (autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k), &
                  PAR(k)%avg, tracer_local(:, k), &
                  dissolved_organic_matter(k))

             call marbl_compute_large_detritus(k, autotroph_cnt, zooplankton_cnt, autotrophs, &
                  zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k), tracer_local(fe_ind, k), &
                  POC, P_CaCO3, P_SiO2, dust, P_iron, &
                  Fe_scavenge(k), Fe_scavenge_rate(k))

             call marbl_compute_particulate_terms(k, &
                  domain%land_mask, domain%kmt, domain%dzt(k), domain%dz(k), &
                  marbl_particulate_share, &
                  POC, P_CaCO3, P_SiO2, dust, P_iron, &
                  QA_dust_def(k), gcm_state%temperature(k), tracer_local(:, k), &
                  sed_denitrif(k), other_remin(k), fesedflux(k), lexport_shared_vars)

             call marbl_compute_nitrif(PAR_out, PAR(k)%in, PAR(k)%KPARdz, &
                  tracer_local(nh4_ind, k), nitrif(k))

             call marbl_compute_denitrif(tracer_local(o2_ind, k), tracer_local(no3_ind, k), &
                  dissolved_organic_matter(k)%DOC_remin, &
                  POC%remin(k), other_remin(k), sed_denitrif(k), denitrif(k))

             call marbl_compute_dtracer_local (autotroph_cnt, zooplankton_cnt, autotrophs, zooplankton, &
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

             if (lexport_shared_vars) then
                call marbl_export_interior_shared_variables(tracer_local(:, k), &
                     carbonate(k), dissolved_organic_matter(k), &
                     QA_dust_def(k), &
                     marbl_interior_share(k))

                call marbl_export_zooplankton_shared_variables(zooplankton_cnt, &
                     zooplankton_local(:, k), &
                     zooplankton_secondary_species(:, k), &
                     marbl_zooplankton_share(:, k))

                call marbl_export_autotroph_shared_variables(autotroph_cnt, &
                     autotroph_local(:, k), &
                     autotroph_secondary_species(:, k), &
                     marbl_autotroph_share(:, k))

                ! FIXME(bja, 2015-08) need to pull particulate share
                ! out of compute_particulate_terms!
             end if

             if  (k<domain%km) then
                call marbl_update_particulate_terms_from_prior_level(k+1, POC, P_CaCO3, &
                     P_SiO2, dust, P_iron, QA_dust_def(:))
             endif

          end do ! k

          call store_diagnostics_carbonate(carbonate, zsat_calcite,           &
                                           zsat_aragonite, marbl_diagnostics)

          call store_diagnostics_autotrophs(domain,                           &
                                            autotroph_secondary_species,      &
                                            marbl_diagnostics)

          call store_diagnostics_particulates(domain, POC, P_CaCO3, P_SiO2,   &
                                              dust,  P_iron, sed_denitrif,    &
                                              other_remin, marbl_diagnostics)

          call store_diagnostics_autotroph_sums(domain,                       &
                                                autotroph_secondary_species,  &
                                                marbl_diagnostics)

          call store_diagnostics_nitrification(nitrif, denitrif,              &
                                               marbl_diagnostics)

          call store_diagnostics_oxygen(domain, gcm_state, zt, tracer_module(o2_ind, :), &
                                        o2_production, o2_consumption,        &
                                        marbl_diagnostics)

          call store_diagnostics_photosynthetically_available_radiation(PAR,  &
                                                            marbl_diagnostics)

          call store_diagnostics_zooplankton(zooplankton_secondary_species,   &
                                             marbl_diagnostics)

          call store_diagnostics_dissolved_organic_matter(                    &
                                               dissolved_organic_matter,      &
                                               fe_scavenge, fe_scavenge_rate, &
                                               marbl_diagnostics)

          call store_diagnostics_carbon_fluxes(domain, zw, POC, P_CaCO3,      &
                                               dtracer, marbl_diagnostics)

          call store_diagnostics_nitrogen_fluxes(domain, zw, POC, denitrif,   &
                                   sed_denitrif, autotroph_secondary_species, &
                                   dtracer, marbl_diagnostics)

          call store_diagnostics_phosphorus_fluxes(domain, zw, POC, dtracer,  &
                                                   marbl_diagnostics)

          call store_diagnostics_silicon_fluxes(domain, zw, P_SiO2, dtracer,  &
                                                marbl_diagnostics)

          do k = 1, domain%km
             ! store_diagnostics_restore (transpose restore_local!)
             marbl_diagnostics%restore_diags(k, :) = restore_local(:, k)
          end do

          end associate

  end subroutine marbl_ecosys_set_interior

  !***********************************************************************

  subroutine marbl_init_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, NET_DUST_IN)

    ! !DESCRIPTION:
    !  Set incoming fluxes (put into outgoing flux for first level usage).
    !  Set dissolution length, production fraction and mass terms.
    !
    !  The first 6 arguments are intent(inout) in
    !  order to preserve contents on other blocks.

    use marbl_share_mod, only : column_sinking_particle_type
    use marbl_share_mod, only : dust_flux        

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

  end subroutine marbl_init_particulate_terms

  !***********************************************************************

  subroutine marbl_update_particulate_terms_from_prior_level(k, POC, P_CaCO3, P_SiO2, dust, P_iron, QA_dust_def)

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
       call marbl_update_sinking_particle_from_prior_level(k, P_CaCO3)

       call marbl_update_sinking_particle_from_prior_level(k, P_SiO2)

       call marbl_update_sinking_particle_from_prior_level(k, dust)

       call marbl_update_sinking_particle_from_prior_level(k, POC)

       call marbl_update_sinking_particle_from_prior_level(k, P_iron)

       QA_dust_def(k) = QA_dust_def(k-1)
    end if

  end subroutine marbl_update_particulate_terms_from_prior_level

  !***********************************************************************

  subroutine marbl_update_sinking_particle_from_prior_level(k, sinking_particle)

    use marbl_share_mod, only : column_sinking_particle_type

    integer (int_kind), intent(in) :: k
    type(column_sinking_particle_type), intent(inout) :: sinking_particle

    ! NOTE(bja, 201504) level k influx is equal to the level k-1 outflux.
    sinking_particle%sflux_out(k) = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_out(k) = sinking_particle%hflux_out(k-1)
    sinking_particle%sflux_in(k)  = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_in(k)  = sinking_particle%hflux_out(k-1)

  end subroutine marbl_update_sinking_particle_from_prior_level

  !***********************************************************************

  subroutine marbl_compute_particulate_terms(k, &
       column_land_mask, column_kmt, column_dzt, column_dz, &
       marbl_particulate_share, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, temperature, tracer_local, sed_denitrif, other_remin, &
       fesedflux, lexport_shared_vars)

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

    integer (int_kind), intent(in) :: k ! vertical model level
    logical(log_kind), intent(in) :: column_land_mask
    integer(int_kind), intent(in) :: column_kmt
    real (r8), intent(in) :: column_dzt, column_dz

    real (r8), intent(in) :: temperature ! temperature for scaling functions bsi%diss

    real (r8), dimension(ecosys_tracer_cnt), intent(in) :: tracer_local ! local copies of model tracer concentrations

    logical (log_kind), intent(in) :: lexport_shared_vars ! flag to save shared_vars or not

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

    real(r8), intent(in) :: fesedflux  ! sedimentary Fe input
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
                  + (fesedflux * dzr_loc)

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
             !        accounted for by fesedflux elsewhere.
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

  end subroutine marbl_compute_particulate_terms

  !***********************************************************************

  subroutine marbl_ecosys_set_sflux(                       &
       saved_state,                                        &
       marbl_surface_share,                                &
       U10_SQR, IFRAC, PRESS, SST, SSS,                    &
       SURF_VALS, MARBL_STF,                               &
       XCO2, XCO2_ALT_CO2, &
       STF_MODULE,                                         &
       IFRAC_USED, XKW_USED, AP_USED, IRON_FLUX_IN,        &
       lexport_shared_vars, PH_PREV, PH_PREV_ALT_CO2, FLUX) 

    use co2calc               , only : co2calc_row
    use schmidt_number        , only : SCHMIDT_CO2
    use marbl_oxygen          , only : schmidt_o2
    use marbl_oxygen          , only : o2sat
    use marbl_share_mod       , only : ecosys_surface_share_type
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : ndep_shr_stream_scale_factor
    use marbl_share_mod       , only : gas_flux_forcing_iopt_drv
    use marbl_share_mod       , only : gas_flux_forcing_iopt_file
    use marbl_share_mod       , only : gas_flux_forcing_iopt
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
    use marbl_share_mod       , only : ecosys_qsw_distrb_const
    use marbl_share_mod       , only : ndep_shr_stream_year_first
    use marbl_share_mod       , only : ndep_shr_stream_year_last
    use marbl_share_mod       , only : ndep_shr_stream_year_align
    use marbl_share_mod       , only : ndep_shr_stream_file
    use marbl_share_mod       , only : ndep_shr_stream_var_cnt
    use marbl_share_mod       , only : ndep_shr_stream_no_ind
    use marbl_share_mod       , only : ndep_shr_stream_nh_ind
    use marbl_share_mod       , only : dust_flux        
    use marbl_share_mod       , only : iron_flux        
    use marbl_share_mod       , only : fice_file        
    use marbl_share_mod       , only : xkw_file         
    use marbl_share_mod       , only : ap_file          
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
    use marbl_share_mod       , only : din_riv_flux     
    use marbl_share_mod       , only : dip_riv_flux     
    use marbl_share_mod       , only : don_riv_flux     
    use marbl_share_mod       , only : dop_riv_flux     
    use marbl_share_mod       , only : dsi_riv_flux     
    use marbl_share_mod       , only : dfe_riv_flux     
    use marbl_share_mod       , only : dic_riv_flux     
    use marbl_share_mod       , only : alk_riv_flux     
    use marbl_share_mod       , only : doc_riv_flux     
    use marbl_parms           , only : ind_nox_flux
    use marbl_parms           , only : ind_nhy_flux
    use marbl_parms           , only : ind_no3_flux
    use marbl_parms           , only : ind_nh4_flux
    use marbl_parms           , only : ind_din_riv_flux
    use marbl_parms           , only : ind_dip_riv_flux
    use marbl_parms           , only : ind_don_riv_flux
    use marbl_parms           , only : ind_dop_riv_flux
    use marbl_parms           , only : ind_dsi_riv_flux
    use marbl_parms           , only : ind_dfe_riv_flux
    use marbl_parms           , only : ind_dic_riv_flux
    use marbl_parms           , only : ind_alk_riv_flux
    use marbl_parms           , only : ind_doc_riv_flux
    use marbl_interface_types , only : marbl_saved_state_type

    !FIXME (mvertens, 2015-11) where does this variable belong?
    use marbl_share_mod           , only : comp_surf_avg_flag 

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.
    !
    ! !INPUT PARAMETERS:

    type(marbl_saved_state_type), intent(in) :: saved_state
    real (r8)          , intent(in) :: U10_SQR(:, :, :)      ! 10m wind speed squared (cm/s)**2
    real (r8)          , intent(in) :: IFRAC  (:, :, :)      ! sea ice fraction (non-dimensional)
    real (r8)          , intent(in) :: PRESS  (:, :, :)      ! sea level atmospheric pressure (dyne/cm**2)
    real (r8)          , intent(in) :: SST    (:, :, :)      ! sea surface temperature (C)
    real (r8)          , intent(in) :: SSS    (:, :, :)      ! sea surface salinity (psu)
    real (r8)          , intent(in) :: SURF_VALS(:, :, :, :) ! module tracers
    real (r8)          , intent(in) :: MARBL_STF(:, :, :, :) 
    real (r8)          , intent(in) :: XCO2(:, :, :)         ! atmospheric co2 conc. (dry-air, 1 atm)
    real (r8)          , intent(in) :: XCO2_ALT_CO2(:, :, :) ! atmospheric alternative CO2 (dry-air, 1 atm)
    logical (log_kind) , intent(in) :: lexport_shared_vars   ! flag to save shared_vars or not

    ! !INPUT/OUTPUT PARAMETERS:
    type(ecosys_surface_share_type) , intent(inout) :: marbl_surface_share
    real (r8), dimension(:, :, :)   , intent(inout) :: IFRAC_USED      ! used ice fraction (non-dimensional)
    real (r8), dimension(:, :, :)   , intent(inout) :: XKW_USED        ! portion of piston velocity (cm/s)
    real (r8), dimension(:, :, :)   , intent(inout) :: AP_USED         ! used atm pressure (converted from dyne/cm**2 to atm)
    real (r8), dimension(:, :, :)   , intent(inout) :: IRON_FLUX_IN    ! iron flux
    real (r8), dimension(:, :, :)   , intent(inout) :: PH_PREV         ! computed ph from previous time step
    real (r8), dimension(:, :, :)   , intent(inout) :: PH_PREV_ALT_CO2 ! computed ph from previous time step

    ! !OUTPUT PARAMETERS:
    real (r8), dimension(:, :, :, :) , intent(out) :: STF_MODULE
    real (r8), dimension(:, :, :)    , intent(out) :: FLUX              ! tracer flux (nmol/cm^2/s)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_set_sflux'

    logical (log_kind) :: first_call = .true.

    integer (int_kind) :: &
         i, j, iblock, n, & ! loop indices
         auto_ind           ! autotroph functional group index

    ! the following arrays need to be used by both set_sflux and store_sflux
    ! in the original routine, they were dimensioned ny_block and reused
    ! could also go into a marbl datatype
    real (r8), dimension(nx_block, ny_block, max_blocks_clinic) :: &
         CO2STAR,      &
         DCO2STAR,     &
         pCO2SURF,     &
         DpCO2,        &
         CO2STAR_ALT,  &
         DCO2STAR_ALT, &
         pCO2SURF_ALT, &
         DpCO2_ALT

    real (r8), dimension(nx_block, ny_block, max_blocks_clinic) :: &
         SCHMIDT_USED_CO2, & ! used Schmidt number
         SCHMIDT_USED_O2,  & ! used Schmidt number
         PV,               & ! piston velocity (cm/s)
         O2SAT_USED,       & ! used O2 saturation (mmol/m^3)
         FLUX_ALT_CO2        ! tracer flux alternative CO2 (nmol/cm^2/s)

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

    real (r8), dimension(nx_block, ny_block) :: &
         XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
         O2SAT_1atm      ! O2 saturation @ 1 atm (mmol/m^3)

    character (char_len) :: &
         tracer_data_label,       & ! label for what is being updated
         ndep_shr_stream_fldList

    character (char_len), dimension(1) :: &
         tracer_data_names          ! short names for input data fields

    integer (int_kind), dimension(1) :: &
         tracer_bndy_loc,         & ! location and field type for ghost
         tracer_bndy_type           !    cell updates

    !-----------------------------------------------------------------------

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
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  fluxes initially set to 0
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic
       STF_MODULE(:, :, :, iblock) = c0
    enddo

    !-----------------------------------------------------------------------
    !  compute CO2 flux, computing disequilibrium one row at a time
    !-----------------------------------------------------------------------

    if (lflux_gas_o2 .or. lflux_gas_co2) then

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

          !-----------------------------------------------------------------------
          !  Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
          !-----------------------------------------------------------------------

          XKW_ICE = (c1 - IFRAC_USED(:, :, iblock)) * XKW_USED(:, :, iblock)

          !-----------------------------------------------------------------------
          !  compute O2 flux
          !-----------------------------------------------------------------------

          if (lflux_gas_o2) then
             SCHMIDT_USED_O2(:, :, iblock) = SCHMIDT_O2(nx_block, ny_block, SST(:, :, iblock), saved_state%land_mask(:, :, iblock))

             O2SAT_1atm = O2SAT(nx_block, ny_block, SST(:, :, iblock), SSS(:, :, iblock), saved_state%land_mask(:, :, iblock))

             where (saved_state%land_mask(:, :, iblock))
                PV(:, :, iblock)         = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED_O2(:, :, iblock))
                O2SAT_USED(:, :, iblock) = AP_USED(:, :, iblock) * O2SAT_1atm
                FLUX(:, :, iblock)       = PV(:, :, iblock) * (O2SAT_USED(:, :, iblock) &
                                         - SURF_VALS(:, :, o2_ind, iblock))
                STF_MODULE(:, :, o2_ind, iblock) = STF_MODULE(:, :, o2_ind, iblock) + FLUX(:, :, iblock)
             elsewhere
                O2SAT_USED(:, :, iblock) = c0
             end where

          endif  ! lflux_gas_o2

          !-----------------------------------------------------------------------
          !  compute CO2 flux, computing disequilibrium one row at a time
          !-----------------------------------------------------------------------

          if (lflux_gas_co2) then

             SCHMIDT_USED_CO2(:, :, iblock) = SCHMIDT_CO2(SST(:, :, iblock), saved_state%land_mask(:, :, iblock))

             where (saved_state%land_mask(:, :, iblock))
                PV(:, :, iblock) = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED_CO2(:, :, iblock))
             elsewhere
                PV(:, :, iblock) = c0
             end where

             ! Save surface field of PV for use in other modules
             if (lexport_shared_vars) PV_SURF_fields(:, :, iblock) = PV(:, :, iblock)

             !-----------------------------------------------------------------------
             !  Set XCO2
             !-----------------------------------------------------------------------

             do j = 1, ny_block
                where (PH_PREV(:, j, iblock) /= c0)
                   PHLO = PH_PREV(:, j, iblock) - del_ph
                   PHHI = PH_PREV(:, j, iblock) + del_ph
                elsewhere
                   PHLO = phlo_surf_init
                   PHHI = phhi_surf_init
                end where

                DIC_ROW  = SURF_VALS(:, j,  dic_ind, iblock)
                ALK_ROW  = SURF_VALS(:, j,  alk_ind, iblock)
                PO4_ROW  = SURF_VALS(:, j,  po4_ind, iblock)
                SiO3_ROW = SURF_VALS(:, j, sio3_ind, iblock)

                call co2calc_row(iblock, j, saved_state%land_mask(:, j, iblock), &
                     locmip_k1_k2_bug_fix, .true., &
                     SST(:, j, iblock), SSS(:, j, iblock), &
                     DIC_ROW, ALK_ROW, PO4_ROW, SiO3_ROW, &
                     PHLO, PHHI, PH_NEW, XCO2(:, j, iblock), &
                     AP_USED(:, j, iblock), CO2STAR_ROW, &
                     DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, &
                     CO3_ROW)

                CO2STAR(:, j, iblock)  = CO2STAR_ROW
                DCO2STAR(:, j, iblock) = DCO2STAR_ROW
                pCO2SURF(:, j, iblock) = pCO2SURF_ROW
                DpCO2(:, j, iblock)    = DpCO2_ROW
                PH_PREV(:, j, iblock)  = PH_NEW

                FLUX(:, j, iblock) = PV(:, j, iblock) * DCO2STAR_ROW
 
                !-------------------------------------------------------------------
                !  The following variables need to be shared with other modules,
                !  and are now defined in marbl_share as targets.
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

                DIC_ROW = SURF_VALS(:, j, dic_alt_co2_ind, iblock)

                call co2calc_row(iblock, j, saved_state%land_mask(:, j, iblock), &
                     locmip_k1_k2_bug_fix, .false., &
                     SST(:, j, iblock), SSS(:, j, iblock), &
                     DIC_ROW, ALK_ROW, PO4_ROW, SiO3_ROW, &
                     PHLO, PHHI, PH_NEW, XCO2_ALT_CO2(:, j, iblock), &
                     AP_USED(:, j, iblock), CO2STAR_ROW, &
                     DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW, &
                     CO3_ROW)

                CO2STAR_ALT(:, j, iblock)  = CO2STAR_ROW
                DCO2STAR_ALT(:, j, iblock) = DCO2STAR_ROW
                pCO2SURF_ALT(:, j, iblock) = pCO2SURF_ROW
                DpCO2_ALT(:, j, iblock)    = DpCO2_ROW
                PH_PREV_ALT_CO2(:, j, iblock) = PH_NEW

                FLUX_ALT_CO2(:, j, iblock) = PV(:, j, iblock) * DCO2STAR_ROW

             end do

             !-----------------------------------------------------------------------
             !  set air-sea co2 gas flux named field, converting units from
             !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
             !-----------------------------------------------------------------------

             STF_MODULE(:, :, dic_ind, iblock) = STF_MODULE(:, :, dic_ind, iblock) + FLUX(:, :, iblock)

             STF_MODULE(:, :, dic_alt_co2_ind, iblock) = STF_MODULE(:, :, dic_alt_co2_ind, iblock) + FLUX_ALT_CO2(:, :, iblock)

          endif  !  lflux_gas_co2

       enddo

    endif  ! lflux_gas_o2 .or. lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    IRON_FLUX_IN = IRON_FLUX_IN * parm_Fe_bioavail
   
    STF_MODULE(:, :, fe_ind, :) = STF_MODULE(:, :, fe_ind, :) + IRON_FLUX_IN

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------

    if (nox_flux_monthly%has_data) then
       STF_MODULE(:, :, no3_ind, :) = STF_MODULE(:, :, no3_ind, :) + MARBL_STF(:, :, ind_nox_flux, :)
    endif

    if (nhy_flux_monthly%has_data) then
       STF_MODULE(:, :, nh4_ind, :) = STF_MODULE(:, :, nh4_ind, :) + MARBL_STF(:, :, ind_nhy_flux, :)
    endif

    if (trim(ndep_data_type) == 'shr_stream') then
       do iblock = 1, nblocks_clinic
          where (saved_state%land_mask(:, :, iblock))
             STF_MODULE(:, :, no3_ind, iblock) = STF_MODULE(:, :, no3_ind, iblock) &
                  + ndep_shr_stream_scale_factor * MARBL_STF(:, :, ind_no3_flux, iblock)
          endwhere
       enddo

       do iblock = 1, nblocks_clinic
          where (saved_state%land_mask(:, :, iblock))
             STF_MODULE(:, :, nh4_ind, iblock) = STF_MODULE(:, :, nh4_ind, iblock) &
                  + ndep_shr_stream_scale_factor * MARBL_STF(:, :, ind_nh4_flux, iblock)
          endwhere
       enddo
    endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------

    if (din_riv_flux%has_data) then
       STF_MODULE(:, :, no3_ind, :) = STF_MODULE(:, :, no3_ind, :) + MARBL_STF(:, :, ind_din_riv_flux, :)
    endif

    if (dip_riv_flux%has_data) then
       STF_MODULE(:, :, po4_ind, :) = STF_MODULE(:, :, po4_ind, :) + MARBL_STF(:, :, ind_dip_riv_flux, :)
    endif

    if (don_riv_flux%has_data) then
       STF_MODULE(:, :, don_ind, :)  = STF_MODULE(:, :, don_ind, :)  + (MARBL_STF(:, :, ind_don_riv_flux, :) * 0.9_r8)
       STF_MODULE(:, :, donr_ind, :) = STF_MODULE(:, :, donr_ind, :) + (MARBL_STF(:, :, ind_don_riv_flux, :) * 0.1_r8)
    endif

    if (dop_riv_flux%has_data) then
       STF_MODULE(:, :, dop_ind, :)  = STF_MODULE(:, :, dop_ind, :)  + (MARBL_STF(:, :, ind_dop_riv_flux, :) * 0.975_r8)
       STF_MODULE(:, :, dopr_ind, :) = STF_MODULE(:, :, dopr_ind, :) + (MARBL_STF(:, :, ind_dop_riv_flux, :) * 0.025_r8)
    endif

    if (dsi_riv_flux%has_data) then
       STF_MODULE(:, :, sio3_ind, :) = STF_MODULE(:, :, sio3_ind, :) + MARBL_STF(:, :, ind_dsi_riv_flux, :)
    endif

    if (dfe_riv_flux%has_data) then
       STF_MODULE(:, :, fe_ind, :) = STF_MODULE(:, :, fe_ind, :) + MARBL_STF(:, :, ind_dfe_riv_flux, :)
    endif

    if (dic_riv_flux%has_data) then
       STF_MODULE(:, :, dic_ind, :) = STF_MODULE(:, :, dic_ind, :) + MARBL_STF(:, :, ind_dic_riv_flux, :)
       STF_MODULE(:, :, dic_alt_co2_ind, :) = STF_MODULE(:, :, dic_alt_co2_ind, :) + MARBL_STF(:, :, ind_dic_riv_flux, :)
       if (lexport_shared_vars) dic_riv_flux_fields = MARBL_STF(:, :, ind_dic_riv_flux, :)
    endif

    if (alk_riv_flux%has_data) then
       STF_MODULE(:, :, alk_ind, :) = STF_MODULE(:, :, alk_ind, :) + MARBL_STF(:, :, ind_alk_riv_flux, :)
    endif

    if (doc_riv_flux%has_data) then
       STF_MODULE(:, :, doc_ind, :) = STF_MODULE(:, :, doc_ind, :) + MARBL_STF(:, :, ind_doc_riv_flux, :)
       !JW change INTERP_WORK to MARBL_STF
       if (lexport_shared_vars) doc_riv_flux_fields=MARBL_STF(:, :, ind_doc_riv_flux, :)
    endif

    !-----------------------------------------------------------------------
    !  Apply NO & NH fluxes to alkalinity
    !-----------------------------------------------------------------------

    STF_MODULE(:, :, alk_ind, :) = STF_MODULE(:, :, alk_ind, :) &
                                 + STF_MODULE(:, :, nh4_ind, :) &
                                 - STF_MODULE(:, :, no3_ind, :)

    !-----------------------------------------------------------------------

    end associate

    call marbl_ecosys_store_sflux(                           &
         STF_MODULE,                                         &
         SURF_VALS, MARBL_STF,                               &
         SCHMIDT_USED_CO2, SCHMIDT_USED_O2, PV, O2SAT_USED,  &
         XCO2, XCO2_ALT_CO2, FLUX, FLUX_ALT_CO2, IFRAC_USED, &
         XKW_USED, AP_USED, IRON_FLUX_IN,                    &
         CO2STAR, DCO2STAR, pCO2SURF, DpCO2,                 &
         CO2STAR_ALT, DCO2STAR_ALT, pCO2SURF_ALT, DpCO2_ALT, &
         PH_PREV, PH_PREV_ALT_CO2,                           &
         lexport_shared_vars)

  end subroutine marbl_ecosys_set_sflux

  !***********************************************************************

  subroutine marbl_ecosys_store_sflux(                           &
       STF_MODULE,                                         &
       SURF_VALS, MARBL_STF,                               &
       SCHMIDT_USED_CO2, SCHMIDT_USED_O2, PV, O2SAT_USED,  & 
       XCO2, XCO2_ALT_CO2, FLUX, FLUX_ALT_CO2, IFRAC_USED, & 
       XKW_USED, AP_USED, IRON_FLUX_IN,                    &
       CO2STAR, DCO2STAR, pCO2SURF, DpCO2,                 &
       CO2STAR_ALT, DCO2STAR_ALT, pCO2SURF_ALT, DpCO2_ALT, &
       PH_PREV, PH_PREV_ALT_CO2,                           &
       lexport_shared_vars)

    use marbl_share_mod       , only : ecosys_surface_share_type
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : ndep_shr_stream_scale_factor
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
    use marbl_share_mod       , only : iron_flux     
    use marbl_share_mod       , only : din_riv_flux     
    use marbl_share_mod       , only : dip_riv_flux     
    use marbl_share_mod       , only : don_riv_flux     
    use marbl_share_mod       , only : dop_riv_flux     
    use marbl_share_mod       , only : dsi_riv_flux     
    use marbl_share_mod       , only : dfe_riv_flux     
    use marbl_share_mod       , only : dic_riv_flux     
    use marbl_share_mod       , only : alk_riv_flux     
    use marbl_parms           , only : ind_nox_flux
    use marbl_parms           , only : ind_no3_flux
    use marbl_parms           , only : ind_din_riv_flux
    use marbl_parms           , only : ind_dfe_riv_flux
    use marbl_parms           , only : ind_dic_riv_flux
    use marbl_parms           , only : ind_alk_riv_flux

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.

    ! !INPUT PARAMETERS:
    logical (log_kind), intent(in) :: lexport_shared_vars ! flag to save shared_vars or not

    real (r8), dimension(:, :, :)    , intent(in) :: SCHMIDT_USED_CO2 ! used Schmidt number
    real (r8), dimension(:, :, :)    , intent(in) :: SCHMIDT_USED_O2  ! used Schmidt number
    real (r8), dimension(:, :, :)    , intent(in) :: PV               ! piston velocity (cm/s)
    real (r8), dimension(:, :, :)    , intent(in) :: O2SAT_USED       ! used O2 saturation (mmol/m^3)
    real (r8), dimension(:, :, :)    , intent(in) :: XCO2             ! atmospheric co2 conc. (dry-air, 1 atm)
    real (r8), dimension(:, :, :)    , intent(in) :: XCO2_ALT_CO2     ! atmospheric alternative CO2 (dry-air, 1 atm)
    real (r8), dimension(:, :, :)    , intent(in) :: FLUX             ! tracer flux (nmol/cm^2/s)
    real (r8), dimension(:, :, :)    , intent(in) :: FLUX_ALT_CO2     ! tracer flux alternative CO2 (nmol/cm^2/s)
    real (r8), dimension(:, :, :)    , intent(in) :: IFRAC_USED       ! used ice fraction (non-dimensional)
    real (r8), dimension(:, :, :)    , intent(in) :: XKW_USED         ! portion of piston velocity (cm/s)
    real (r8), dimension(:, :, :)    , intent(in) :: AP_USED          ! used atm pressure (converted from dyne/cm**2 to atm)
    real (r8), dimension(:, :, :)    , intent(in) :: IRON_FLUX_IN     ! iron flux! 
    real (r8), dimension(:, :, :)    , intent(in) :: CO2STAR
    real (r8), dimension(:, :, :)    , intent(in) :: DCO2STAR
    real (r8), dimension(:, :, :)    , intent(in) :: pCO2SURF
    real (r8), dimension(:, :, :)    , intent(in) :: DpCO2
    real (r8), dimension(:, :, :)    , intent(in) :: CO2STAR_ALT
    real (r8), dimension(:, :, :)    , intent(in) :: DCO2STAR_ALT
    real (r8), dimension(:, :, :)    , intent(in) :: pCO2SURF_ALT
    real (r8), dimension(:, :, :)    , intent(in) :: DpCO2_ALT
    real (r8), dimension(:, :, :, :) , intent(in) :: STF_MODULE
    real (r8), dimension(:, :, :, :) , intent(in) :: MARBL_STF
    real (r8), dimension(:, :, :, :) , intent(in) :: SURF_VALS
    real (r8), dimension(:, :, :)    , intent(in) :: PH_PREV         
    real (r8), dimension(:, :, :)    , intent(in) :: PH_PREV_ALT_CO2 

    !JW    type(ecosys_surface_share_type), intent(inout) :: marbl_surface_share

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:ecosys_store_sflux'

    logical (log_kind) :: first_call = .true.

    integer (int_kind) :: &
         i, j, iblock, n, & ! loop indices
         auto_ind           ! autotroph functional group index

    real (r8), dimension(nx_block, ny_block, max_blocks_clinic) :: &
         SHR_STREAM_WORK

    character (char_len) :: &
         tracer_data_label,       & ! label for what is being updated
         ndep_shr_stream_fldList

    character (char_len), dimension(1) :: &
         tracer_data_names          ! short names for input data fields

    integer (int_kind), dimension(1) :: &
         tracer_bndy_loc,         & ! location and field type for ghost
         tracer_bndy_type           !    cell updates

    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

   if (lflux_gas_o2 .or. lflux_gas_co2) then

       do iblock = 1, nblocks_clinic

          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_IFRAC, iblock)     = IFRAC_USED(:, :, iblock)
          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_XKW, iblock)       = XKW_USED(:, :, iblock)
          ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_ATM_PRESS, iblock) = AP_USED(:, :, iblock)

          if (lflux_gas_o2) then

             !JW this could be in post, but would require returning  SCHMIDT and O2SAT_USED
             ECO_SFLUX_TAVG(:, :, buf_ind_PV_O2, iblock)      = PV(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_O2, iblock) = SCHMIDT_USED_O2(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_O2SAT, iblock)      = O2SAT_USED(:, :, iblock)

          endif  ! lflux_gas_o2

          !-----------------------------------------------------------------------
          !  compute CO2 flux, computing disequilibrium one row at a time
          !-----------------------------------------------------------------------

          if (lflux_gas_co2) then

             ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR, iblock)  = CO2STAR(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR, iblock) = DCO2STAR(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF, iblock) = pCO2SURF(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2, iblock)    = DpCO2(:, :, iblock)

             ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR_ALT_CO2, iblock)  = CO2STAR_ALT(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR_ALT_CO2, iblock) = DCO2STAR_ALT(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF_ALT_CO2, iblock) = pCO2SURF_ALT(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2_ALT_CO2, iblock)    = DpCO2_ALT(:, :, iblock)

   
             ECO_SFLUX_TAVG(:, :, buf_ind_PV_CO2,       iblock) = PV(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_CO2,  iblock) = SCHMIDT_USED_CO2(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX, iblock) = FLUX(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_PH,           iblock) = PH_PREV(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_ATM_CO2,      iblock) = XCO2(:, :, iblock)

             ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX_ALT_CO2, iblock) = FLUX_ALT_CO2(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_PH_ALT_CO2,           iblock) = PH_PREV_ALT_CO2(:, :, iblock)
             ECO_SFLUX_TAVG(:, :, buf_ind_ATM_ALT_CO2,          iblock) = XCO2_ALT_CO2(:, :, iblock)

          endif  !  lflux_gas_co2

       enddo

    endif  ! lflux_gas_o2 .or. lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    if (iron_flux%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_IRON_FLUX, :) = IRON_FLUX_IN
    endif

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------

    if (nox_flux_monthly%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, :) = MARBL_STF(:, :, ind_nox_flux, :)
    endif

    if (trim(ndep_data_type) == 'shr_stream') then
       do iblock = 1, nblocks_clinic
          ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, iblock) = &
                 ndep_shr_stream_scale_factor * MARBL_STF(:, :, ind_no3_flux, iblock)
       enddo
    endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------
 
    if (din_riv_flux%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_DIN_RIV_FLUX, :) = MARBL_STF(:, :, ind_din_riv_flux, :)
    endif

    if (dfe_riv_flux%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_DFE_RIV_FLUX, :) = MARBL_STF(:, :, ind_dfe_riv_flux, :)
    endif

    if (dic_riv_flux%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_DIC_RIV_FLUX, :) = MARBL_STF(:, :, ind_dic_riv_flux, :)
    endif

    if (alk_riv_flux%has_data) then
       ECO_SFLUX_TAVG(:, :, buf_ind_ALK_RIV_FLUX, :) = MARBL_STF(:, :, ind_alk_riv_flux, :)
    endif

  end subroutine marbl_ecosys_store_sflux

  !*****************************************************************************

  subroutine marbl_ecosys_tavg_forcing(saved_state, STF_MODULE, FLUX_DIAGS)

    ! !DESCRIPTION:
    !  Accumulate non-standard forcing related tavg variables.

    use marbl_interface_types, only : marbl_saved_state_type

    type(marbl_saved_state_type), intent(in) :: saved_state

    ! !INPUT PARAMETERS:
    real (r8), intent(in) :: STF_MODULE(:, :, :, :)

    ! !OUTPUT PARAMETERS:
    real (r8),  intent(out) :: FLUX_DIAGS(nx_block, ny_block, forcing_diag_cnt, nblocks_clinic)  ! Computed diagnostics for surface fluxes

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! accumulate surface flux related fields in the order in which they are declared
    !
    !  multiply IRON, DUST fluxes by mpercm (.01) to convert from model
    !    units (cm/s)(mmol/m^3) to mmol/s/m^2
    !-----------------------------------------------------------------------

    FLUX_DIAGS(:, :, ECOSYS_IFRAC_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_IFRAC, :)
    FLUX_DIAGS(:, :, ECOSYS_XKW_diag_ind, :)           = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_XKW, :)
    FLUX_DIAGS(:, :, ECOSYS_ATM_PRESS_diag_ind, :)     = ECO_SFLUX_TAVG(:, :, buf_ind_ECOSYS_ATM_PRESS, :)
    FLUX_DIAGS(:, :, PV_O2_diag_ind, :)                = ECO_SFLUX_TAVG(:, :, buf_ind_PV_O2, :)
    FLUX_DIAGS(:, :, SCHMIDT_O2_diag_ind, :)           = ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_O2, :)
    FLUX_DIAGS(:, :, O2SAT_diag_ind, :)                = ECO_SFLUX_TAVG(:, :, buf_ind_O2SAT, :)
    FLUX_DIAGS(:, :, O2_GAS_FLUX_diag_ind, :)          = STF_MODULE(:, :, o2_ind, :)
    FLUX_DIAGS(:, :, CO2STAR_diag_ind, :)              = ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR, :)
    FLUX_DIAGS(:, :, DCO2STAR_diag_ind, :)             = ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR, :)
    FLUX_DIAGS(:, :, pCO2SURF_diag_ind, :)             = ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF, :)
    FLUX_DIAGS(:, :, DpCO2_diag_ind, :)                = ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2, :)
    FLUX_DIAGS(:, :, PV_CO2_diag_ind, :)               = ECO_SFLUX_TAVG(:, :, buf_ind_PV_CO2, :)
    FLUX_DIAGS(:, :, SCHMIDT_CO2_diag_ind, :)          = ECO_SFLUX_TAVG(:, :, buf_ind_SCHMIDT_CO2, :)
    FLUX_DIAGS(:, :, DIC_GAS_FLUX_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX, :)
    FLUX_DIAGS(:, :, PH_diag_ind, :)                   = ECO_SFLUX_TAVG(:, :, buf_ind_PH, :)
    FLUX_DIAGS(:, :, ATM_CO2_diag_ind, :)              = ECO_SFLUX_TAVG(:, :, buf_ind_ATM_CO2, :)
    FLUX_DIAGS(:, :, CO2STAR_ALT_CO2_diag_ind, :)      = ECO_SFLUX_TAVG(:, :, buf_ind_CO2STAR_ALT_CO2, :)
    FLUX_DIAGS(:, :, DCO2STAR_ALT_CO2_diag_ind, :)     = ECO_SFLUX_TAVG(:, :, buf_ind_DCO2STAR_ALT_CO2, :)
    FLUX_DIAGS(:, :, pCO2SURF_ALT_CO2_diag_ind, :)     = ECO_SFLUX_TAVG(:, :, buf_ind_pCO2SURF_ALT_CO2, :)
    FLUX_DIAGS(:, :, DpCO2_ALT_CO2_diag_ind, :)        = ECO_SFLUX_TAVG(:, :, buf_ind_DpCO2_ALT_CO2, :)
    FLUX_DIAGS(:, :, DIC_GAS_FLUX_ALT_CO2_diag_ind, :) = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_GAS_FLUX_ALT_CO2, :)
    FLUX_DIAGS(:, :, PH_ALT_CO2_diag_ind, :)           = ECO_SFLUX_TAVG(:, :, buf_ind_PH_ALT_CO2, :)
    FLUX_DIAGS(:, :, ATM_ALT_CO2_diag_ind, :)          = ECO_SFLUX_TAVG(:, :, buf_ind_ATM_ALT_CO2, :)
    FLUX_DIAGS(:, :, IRON_FLUX_diag_ind, :)            = ECO_SFLUX_TAVG(:, :, buf_ind_IRON_FLUX, :)*mpercm
    FLUX_DIAGS(:, :, DUST_FLUX_diag_ind, :)            = saved_state%dust_FLUX_IN(:, :, :)*mpercm
    FLUX_DIAGS(:, :, NOx_FLUX_diag_ind, :)             = ECO_SFLUX_TAVG(:, :, buf_ind_NOx_FLUX, :)
    FLUX_DIAGS(:, :, NHy_FLUX_diag_ind, :)             = STF_MODULE(:, :, nh4_ind, :)
    FLUX_DIAGS(:, :, DIN_RIV_FLUX_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_DIN_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DIP_RIV_FLUX_diag_ind, :)         = STF_MODULE(:, :, po4_ind, :)
    FLUX_DIAGS(:, :, DON_RIV_FLUX_diag_ind, :)         = STF_MODULE(:, :, don_ind, :)
    FLUX_DIAGS(:, :, DONr_RIV_FLUX_diag_ind, :)        = STF_MODULE(:, :, donr_ind, :)
    FLUX_DIAGS(:, :, DOP_RIV_FLUX_diag_ind, :)         = STF_MODULE(:, :, dop_ind, :)
    FLUX_DIAGS(:, :, DOPr_RIV_FLUX_diag_ind, :)        = STF_MODULE(:, :, dopr_ind, :)
    FLUX_DIAGS(:, :, DSI_RIV_FLUX_diag_ind, :)         = STF_MODULE(:, :, sio3_ind, :)
    FLUX_DIAGS(:, :, DFE_RIV_FLUX_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_DFE_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DIC_RIV_FLUX_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_DIC_RIV_FLUX, :)
    FLUX_DIAGS(:, :, ALK_RIV_FLUX_diag_ind, :)         = ECO_SFLUX_TAVG(:, :, buf_ind_ALK_RIV_FLUX, :)
    FLUX_DIAGS(:, :, DOC_RIV_FLUX_diag_ind, :)         = STF_MODULE(:, :, doc_ind, :)

  end subroutine marbl_ecosys_tavg_forcing

  !***********************************************************************

  subroutine marbl_ecosys_init_forcing_metadata()

    !-----------------------------------------------------------------------
    ! initialize surface forcing metadata
    !-----------------------------------------------------------------------

    use passive_tracer_tools , only : init_forcing_monthly_every_ts
    use marbl_share_mod      , only : dust_flux        
    use marbl_share_mod      , only : iron_flux        
    use marbl_share_mod      , only : fice_file        
    use marbl_share_mod      , only : xkw_file         
    use marbl_share_mod      , only : ap_file          
    use marbl_share_mod      , only : nox_flux_monthly 
    use marbl_share_mod      , only : nhy_flux_monthly 
    use marbl_share_mod      , only : din_riv_flux     
    use marbl_share_mod      , only : dip_riv_flux     
    use marbl_share_mod      , only : don_riv_flux     
    use marbl_share_mod      , only : dop_riv_flux     
    use marbl_share_mod      , only : dsi_riv_flux     
    use marbl_share_mod      , only : dfe_riv_flux     
    use marbl_share_mod      , only : dic_riv_flux     
    use marbl_share_mod      , only : alk_riv_flux     
    use marbl_share_mod      , only : doc_riv_flux     

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

  end subroutine marbl_ecosys_init_forcing_metadata

  !***********************************************************************

  subroutine marbl_ecosys_init_non_autotroph_tracer_metadata(tracer_d_module, &
       non_living_biomass_ecosys_tracer_cnt)

    !-----------------------------------------------------------------------
    !  initialize non-autotroph tracer_d values and accumulate
    !  non_living_biomass_ecosys_tracer_cnt
    !-----------------------------------------------------------------------

    use marbl_interface_types     , only : marbl_tracer_metadata_type

    implicit none

    type (marbl_tracer_metadata_type), dimension(:), intent(inout) :: tracer_d_module   ! descriptors for each tracer

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

  end subroutine marbl_ecosys_init_non_autotroph_tracer_metadata

  !***********************************************************************

  subroutine marbl_check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt)

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

  end subroutine marbl_check_ecosys_tracer_count_consistency

  !***********************************************************************

  subroutine marbl_initialize_zooplankton_tracer_metadata(tracer_d_module, &
       non_living_biomass_ecosys_tracer_cnt, n)

    use marbl_interface_types     , only : marbl_tracer_metadata_type

    !-----------------------------------------------------------------------
    !  initialize zooplankton tracer_d values and tracer indices
    !-----------------------------------------------------------------------

    type (marbl_tracer_metadata_type), dimension(:), intent(inout) :: tracer_d_module   ! descriptors for each tracer

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

  end subroutine marbl_initialize_zooplankton_tracer_metadata

  !***********************************************************************

  subroutine marbl_initialize_autotroph_tracer_metadata(tracer_d_module, n)

    !-----------------------------------------------------------------------
    !  initialize autotroph tracer_d values and tracer indices
    !-----------------------------------------------------------------------

    use marbl_interface_types     , only : marbl_tracer_metadata_type

    type (marbl_tracer_metadata_type), dimension(:), intent(inout) :: tracer_d_module   ! descriptors for each tracer
    integer(int_kind), intent(inout) :: n

    integer (int_kind) :: auto_ind ! zooplankton functional group index

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

  end subroutine marbl_initialize_autotroph_tracer_metadata

  !***********************************************************************

  subroutine marbl_setup_local_tracers(k, column_land_mask, column_kmt, tracer_module, tracer_local)
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

  end subroutine marbl_setup_local_tracers

  !***********************************************************************

  subroutine marbl_setup_local_zooplankton(k, column_land_mask, column_kmt, &
       tracer_module, zoo_cnt, zoo, zooplankton_local)

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

  end subroutine marbl_setup_local_zooplankton


  !***********************************************************************

  subroutine marbl_setup_local_autotrophs(k, column_land_mask, column_kmt, &
       tracer_module, auto_cnt, auto_meta, autotroph_loc)

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

  end subroutine marbl_setup_local_autotrophs

  !***********************************************************************

  subroutine marbl_autotroph_consistency_check(auto_cnt, auto_meta, &
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

  end subroutine marbl_autotroph_consistency_check

  !***********************************************************************

  subroutine marbl_compute_autotroph_elemental_ratios(k, auto_cnt, auto_meta, autotroph_local, &
       tracer_local, autotroph_secondary_species)

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
  end subroutine marbl_compute_autotroph_elemental_ratios

  !***********************************************************************

  subroutine marbl_compute_photosynthetically_available_radiation(k, auto_cnt, &
       autotroph_local, &
       column_land_mask, column_kmt, column_dzt, column_dz, &
       PAR_out, PAR)
    !-----------------------------------------------------------------------
    !  compute PAR related quantities
    !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
    !  0.45   fraction of incoming SW -> PAR (non-dim)
    !-----------------------------------------------------------------------

    integer(int_kind)          , intent(in) :: k, auto_cnt
    type(autotroph_local_type) , intent(in) :: autotroph_local(auto_cnt)
    logical(log_kind)          , intent(in) :: column_land_mask
    integer(int_kind)          , intent(in) :: column_kmt
    real (r8)                  , intent(in) :: column_dzt
    real (r8)                  , intent(in) :: column_dz

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


  end subroutine marbl_compute_photosynthetically_available_radiation

  !***********************************************************************

  subroutine marbl_compute_carbonate_chemistry(dkm, &
       column_land_mask, column_kmt, &
       temperature, salinity, &
       tracer_local, carbonate, &
       ph_prev_3d, ph_prev_alt_co2_3d, &
       zsat_calcite, zsat_aragonite, &
       co3_calcite_anom, co3_aragonite_anom)

    use co2calc_column, only : comp_co3terms
    use co2calc_column, only : comp_co3_sat_vals
    use co2calc_column, only : thermodynamic_coefficients_type
    use state_mod,      only : ref_pressure

    integer (int_kind), intent(in) :: dkm
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

  end subroutine marbl_compute_carbonate_chemistry

  !***********************************************************************

  subroutine marbl_compute_function_scaling(column_temperature, Tfunc )

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

  end subroutine marbl_compute_function_scaling

  !***********************************************************************

  subroutine marbl_compute_Pprime(k, auto_cnt, auto_meta, &
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
  end subroutine marbl_compute_Pprime

  !***********************************************************************

  subroutine marbl_compute_Zprime(k, zoo_cnt, zoo_meta, zooC, &
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
  end subroutine marbl_compute_Zprime

  !***********************************************************************

  subroutine marbl_compute_autotroph_uptake (auto_cnt, auto_meta, &
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

  end subroutine marbl_compute_autotroph_uptake

  !***********************************************************************

  subroutine marbl_compute_autotroph_photosynthesis (auto_cnt, auto_meta, &
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

  end subroutine marbl_compute_autotroph_photosynthesis

  !***********************************************************************

  subroutine marbl_compute_autotroph_phyto_diatoms (auto_cnt, auto_meta, &
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

  end subroutine marbl_compute_autotroph_phyto_diatoms

  !***********************************************************************

  subroutine marbl_compute_autotroph_calcification (auto_cnt, auto_meta, &
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
  end subroutine marbl_compute_autotroph_calcification

  !***********************************************************************

  subroutine marbl_compute_autotroph_nfixation (auto_cnt, auto_meta, &
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
  end subroutine marbl_compute_autotroph_nfixation

  !***********************************************************************

  subroutine marbl_compute_autotroph_loss (auto_cnt, auto_meta, &
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
  end subroutine marbl_compute_autotroph_loss

  !***********************************************************************

  subroutine marbl_compute_grazing (auto_cnt, zoo_cnt, grazer_prey_cnt, auto_meta, &
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

  end subroutine marbl_compute_grazing

  !***********************************************************************

  subroutine marbl_compute_routing (auto_cnt, zoo_cnt,  auto_meta, &
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

  end subroutine marbl_compute_routing

  !***********************************************************************

  subroutine marbl_compute_dissolved_organic_matter (auto_cnt, zoo_cnt, auto_meta, &
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
  end subroutine marbl_compute_dissolved_organic_matter

  !***********************************************************************

  subroutine marbl_compute_large_detritus(k, auto_cnt, zoo_cnt, auto_meta, &
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
  end subroutine marbl_compute_large_detritus

  !***********************************************************************

  subroutine marbl_compute_nitrif(PAR_out, PAR_in, KPARdz, NH4_loc, nitrif)

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

  end subroutine marbl_compute_nitrif

  !***********************************************************************

  subroutine marbl_compute_denitrif(O2_loc, NO3_loc, DOC_remin, &
       POC_remin, other_remin, sed_denitrif, denitrif)

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

  end subroutine marbl_compute_denitrif

  !***********************************************************************

  subroutine marbl_compute_dtracer_local (auto_cnt, zoo_cnt, auto_meta, zoo_meta, &
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
  end subroutine marbl_compute_dtracer_local

  !-----------------------------------------------------------------------

  subroutine marbl_export_interior_shared_variables (&
       tracer_local, &
       carbonate, &
       dissolved_organic_matter, &
       QA_dust_def, &
       marbl_interior_share)

    use marbl_share_mod, only : marbl_interior_share_type

    real(r8)                            , intent(in)    :: tracer_local(ecosys_tracer_cnt)
    type(carbonate_type)                , intent(in)    :: carbonate
    type(dissolved_organic_matter_type) , intent(in)    :: dissolved_organic_matter
    real(r8)                            , intent(in)    :: QA_dust_def
    type(marbl_interior_share_type)     , intent(inout) :: marbl_interior_share

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
  end subroutine marbl_export_interior_shared_variables

  !-----------------------------------------------------------------------

  subroutine marbl_export_zooplankton_shared_variables (&
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
  end subroutine marbl_export_zooplankton_shared_variables

  !-----------------------------------------------------------------------

  subroutine marbl_export_autotroph_shared_variables (&
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
  end subroutine marbl_export_autotroph_shared_variables

end module ecosys_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
