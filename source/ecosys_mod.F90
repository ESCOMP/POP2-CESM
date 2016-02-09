! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module ecosys_mod

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
  !  The following are used extensively in this ecosys, so are used at
  !  the module level. The use statements for variables that are only needed
  !  locally are located at the module subprogram level.
  !-----------------------------------------------------------------------

  ! !USES:

  use shr_sys_mod          , only : shr_sys_abort        !FIXME
  use communicate          , only : master_task, my_task !FIXME
  use io_types             , only : stdout               !FIXME
  use constants            , only : T0_Kelvin            !FIXME
  use state_mod            , only : ref_pressure         !FIXME

  use marbl_kinds_mod, only : log_kind
  use marbl_kinds_mod, only : int_kind
  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : char_len
  use co2calc_column , only : thermodynamic_coefficients_type

  use marbl_parms, only : c0
  use marbl_parms, only : c1
  use marbl_parms, only : c2
  use marbl_parms, only : c10
  use marbl_parms, only : mpercm
  use marbl_parms, only : blank_fmt
  use marbl_parms, only : delim_fmt
  use marbl_parms, only : ndelim_fmt
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
  use marbl_parms, only : DOC_reminR_light
  use marbl_parms, only : DON_reminR_light
  use marbl_parms, only : DOP_reminR_light
  use marbl_parms, only : DOC_reminR_dark
  use marbl_parms, only : DON_reminR_dark
  use marbl_parms, only : DOP_reminR_dark
  use marbl_parms, only : DOCr_reminR0
  use marbl_parms, only : DONr_reminR0
  use marbl_parms, only : DOPr_reminR0
  use marbl_parms, only : DOCprod_refract
  use marbl_parms, only : DONprod_refract
  use marbl_parms, only : DOPprod_refract
  use marbl_parms, only : POCremin_refract
  use marbl_parms, only : PONremin_refract
  use marbl_parms, only : POPremin_refract
  use marbl_parms, only : DOCriv_refract
  use marbl_parms, only : DONriv_refract
  use marbl_parms, only : DOPriv_refract
  use marbl_parms, only : f_toDON
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
  use marbl_parms, only : dop_ind
  use marbl_parms, only : dopr_ind
  use marbl_parms, only : donr_ind
  use marbl_parms, only : docr_ind

  use marbl_sizes, only : ecosys_tracer_cnt    
  use marbl_sizes, only : autotroph_cnt
  use marbl_sizes, only : zooplankton_cnt
  use marbl_sizes, only : grazer_prey_cnt

  use marbl_parms, only : grazing  
  use marbl_parms, only : autotrophs
  use marbl_parms, only : zooplankton

  use marbl_share_mod, only : marbl_forcing_ind
  use marbl_share_mod, only : max_forcing_fields

  use marbl_internal_types, only : carbonate_type
  use marbl_internal_types, only : zooplankton_type
  use marbl_internal_types, only : autotroph_type
  use marbl_internal_types, only : zooplankton_secondary_species_type
  use marbl_internal_types, only : autotroph_secondary_species_type
  use marbl_internal_types, only : dissolved_organic_matter_type
  use marbl_internal_types, only : column_sinking_particle_type
  use marbl_internal_types, only : marbl_PAR_type
  use marbl_internal_types, only : marbl_particulate_share_type
  use marbl_internal_types, only : marbl_interior_share_type
  use marbl_internal_types, only : marbl_autotroph_share_type
  use marbl_internal_types, only : marbl_zooplankton_share_type
  use marbl_internal_types, only : marbl_forcing_share_type

  use marbl_interface_types, only : marbl_domain_type
  use marbl_interface_types, only : marbl_gcm_state_type
  use marbl_interface_types, only : marbl_tracer_metadata_type
  use marbl_interface_types, only : marbl_tracer_read_type
  use marbl_interface_types, only : marbl_forcing_input_type
  use marbl_interface_types, only : marbl_forcing_output_type
  use marbl_interface_types, only : marbl_forcing_fields_type
  use marbl_interface_types, only : marbl_diagnostics_type
  use marbl_interface_types, only : forcing_monthly_every_ts

  use marbl_logging             , only : marbl_log_type

  implicit none
  private

  !-----------------------------------------------------------------------
  !  public/private member procedure declarations
  !-----------------------------------------------------------------------

  public  :: marbl_init_nml
  public  :: marbl_sflux_forcing_fields_init
  public  :: marbl_init_tracer_metadata
  public  :: marbl_set_interior_forcing
  public  :: marbl_set_surface_forcing
  public  :: marbl_compute_totalChl

  private :: marbl_init_non_autotroph_tracer_metadata
  private :: marbl_init_forcing_metadata
  private :: marbl_init_particulate_terms
  private :: marbl_init_monthly_forcing_metadata
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
  private :: marbl_compute_PAR
  private :: marbl_compute_carbonate_chemistry                    
  private :: marbl_compute_function_scaling                       
  private :: marbl_compute_Pprime                                 
  private :: marbl_compute_Zprime

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

  type(marbl_tracer_read_type) :: &
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
  !  iron patch fertilization
  !-----------------------------------------------------------------------

  logical (log_kind)  :: liron_patch               ! flag for iron patch fertilization
  character(char_len) :: iron_patch_flux_filename  ! file containing name of iron patch file
  integer (int_kind)  :: iron_patch_month          !  integer month to add patch flux

  !-----------------------------------------------------------------------
  !  bury to sediment options
  !-----------------------------------------------------------------------

  character(char_len) :: caco3_bury_thres_opt    ! option of threshold of caco3 burial ['fixed_depth', 'omega_calc']
  integer (int_kind)  :: caco3_bury_thres_iopt   ! integer version of caco3_bury_thres_opt
  integer (int_kind), parameter :: caco3_bury_thres_iopt_fixed_depth = 1
  integer (int_kind), parameter :: caco3_bury_thres_iopt_omega_calc  = 2
  real (r8)           :: caco3_bury_thres_depth  ! threshold depth for caco3_bury_thres_opt='fixed_depth'

  real (r8) :: PON_bury_coeff ! PON_sed_loss = PON_bury_coeff * Q * POC_sed_loss
                              ! factor is used to avoid overburying PON like POC
                              ! is when total C burial is matched to C riverine input

  real (r8) :: POP_bury_coeff ! POP_sed_loss = POP_bury_coeff * Qp_zoo_pom * POC_sed_loss
                              ! factor is used to enable forced closure of the P cycle
                              ! i.e. POP_sed_loss = P inputs (riverine + atm dep)

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

  type(forcing_monthly_every_ts), target :: dust_flux_loc
  type(forcing_monthly_every_ts), target :: iron_flux_loc
  type(forcing_monthly_every_ts), target :: fice_file_loc
  type(forcing_monthly_every_ts), target :: xkw_file_loc
  type(forcing_monthly_every_ts), target :: ap_file_loc
  type(forcing_monthly_every_ts), target :: nox_flux_monthly_loc
  type(forcing_monthly_every_ts), target :: nhy_flux_monthly_loc
  type(forcing_monthly_every_ts), target :: din_riv_flux_loc
  type(forcing_monthly_every_ts), target :: dip_riv_flux_loc
  type(forcing_monthly_every_ts), target :: don_riv_flux_loc
  type(forcing_monthly_every_ts), target :: dop_riv_flux_loc
  type(forcing_monthly_every_ts), target :: dsi_riv_flux_loc
  type(forcing_monthly_every_ts), target :: dfe_riv_flux_loc
  type(forcing_monthly_every_ts), target :: dic_riv_flux_loc
  type(forcing_monthly_every_ts), target :: alk_riv_flux_loc
  type(forcing_monthly_every_ts), target :: doc_riv_flux_loc

  !-----------------------------------------------------------------------
  !  private module string for storing error messages
  !-----------------------------------------------------------------------

  character(len=char_len), private :: error_msg

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine marbl_init_nml(nl_buffer, marbl_status_log)

    ! !DESCRIPTION:
    !  Initialize ecosys tracer module. This involves setting metadata, reading
    !  the module namelist, setting initial conditions, setting up forcing,
    !  and defining additional tavg variables.
    !
    use marbl_namelist_mod        , only : marbl_nl_cnt
    use marbl_namelist_mod        , only : marbl_nl_buffer_size
    use marbl_namelist_mod        , only : marbl_namelist
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
    use marbl_share_mod           , only : fesedflux_input 
    use marbl_share_mod           , only : marbl_freq_opt_never  
    use marbl_share_mod           , only : marbl_freq_opt_nmonth 
    use marbl_share_mod           , only : marbl_freq_opt_nyear  

    implicit none

    character(marbl_nl_buffer_size), intent(in)  :: nl_buffer(marbl_nl_cnt)
    type(marbl_log_type)           , intent(inout) :: marbl_status_log

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_init_nml'
    character(len=marbl_nl_buffer_size) :: tmp_nl_buffer

    integer (int_kind)           :: n                        ! index for looping over tracers
    character(char_len)          :: comp_surf_avg_freq_opt   ! choice for freq of comp_surf_avg
    character(char_len)          :: gas_flux_forcing_opt     ! option for forcing gas fluxes
    character(char_len)          :: atm_co2_opt              ! option for atmospheric co2 concentration
    character(char_len)          :: atm_alt_co2_opt          ! option for atmospheric alternative CO2
    type(marbl_tracer_read_type) :: dust_flux_input          ! namelist input for dust_flux
    type(marbl_tracer_read_type) :: iron_flux_input          ! namelist input for iron_flux
    type(marbl_tracer_read_type) :: nox_flux_monthly_input   ! namelist input for nox_flux_monthly
    type(marbl_tracer_read_type) :: nhy_flux_monthly_input   ! namelist input for nhy_flux_monthly
    type(marbl_tracer_read_type) :: din_riv_flux_input       ! namelist input for din_riv_flux
    type(marbl_tracer_read_type) :: dip_riv_flux_input       ! namelist input for dip_riv_flux
    type(marbl_tracer_read_type) :: don_riv_flux_input       ! namelist input for don_riv_flux
    type(marbl_tracer_read_type) :: dop_riv_flux_input       ! namelist input for dop_riv_flux
    type(marbl_tracer_read_type) :: dsi_riv_flux_input       ! namelist input for dsi_riv_flux
    type(marbl_tracer_read_type) :: dfe_riv_flux_input       ! namelist input for dfe_riv_flux
    type(marbl_tracer_read_type) :: dic_riv_flux_input       ! namelist input for dic_riv_flux
    type(marbl_tracer_read_type) :: alk_riv_flux_input       ! namelist input for alk_riv_flux
    type(marbl_tracer_read_type) :: doc_riv_flux_input       ! namelist input for doc_riv_flux
    integer (int_kind)           :: nml_error                ! namelist i/o error flag
    integer (int_kind)           :: zoo_ind                  ! zooplankton functional group index
    integer (int_kind)           :: comp_surf_avg_freq_iopt  ! choice for freq of comp_surf_avg
    integer (int_kind)           :: comp_surf_avg_freq       ! choice for freq of comp_surf_avg

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
         lsource_sink, lflux_gas_o2, lflux_gas_co2, locmip_k1_k2_bug_fix, &
         lnutr_variable_restore, nutr_variable_rest_file,                 &
         nutr_variable_rest_file_fmt, atm_co2_opt, atm_co2_const,         &
         atm_alt_co2_opt, atm_alt_co2_const,                              &
         liron_patch, iron_patch_flux_filename, iron_patch_month,         &
         caco3_bury_thres_opt, caco3_bury_thres_depth, &
         PON_bury_coeff, &
         lecovars_full_depth_tavg

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

    liron_patch              = .false.
    iron_patch_flux_filename = 'unknown_iron_patch_filename'
    iron_patch_month         = 1

    atm_co2_opt   = 'const'
    atm_co2_const = 280.0_r8

    atm_alt_co2_opt   = 'const'
    atm_alt_co2_const = 280.0_r8

    caco3_bury_thres_opt = 'omega_calc'
    caco3_bury_thres_depth = 3000.0e2

    PON_bury_coeff = 0.5_r8
    POP_bury_coeff = 1.0_r8

    lecovars_full_depth_tavg = .false.

    !-----------------------------------------------------------------------
    ! read the namelist buffer on every processor
    !-----------------------------------------------------------------------

    tmp_nl_buffer = marbl_namelist(nl_buffer, 'ecosys_nml')
    read(tmp_nl_buffer, nml=ecosys_nml, iostat=nml_error)
    if (nml_error /= 0) then
      call marbl_status_log%log_error("error reading &ecosys_nml", "ecosys_mod::marbl_init_nml()")
       return
    else
      ! FIXME(mnl,2016-02): this is printing contents of pop_in, not the entire
      !                     ecosys_nml
      call marbl_status_log%log_namelist('ecosys_nml', tmp_nl_buffer, 'ecosys_mod::marbl_init_nml')
    end if

    !-----------------------------------------------------------------------
    ! reassign values temporary input values to correct arrays
    !-----------------------------------------------------------------------

    if (trim(gas_flux_forcing_opt) == 'drv') then
       gas_flux_forcing_iopt = gas_flux_forcing_iopt_drv
    else if (trim(gas_flux_forcing_opt) == 'file') then
       gas_flux_forcing_iopt = gas_flux_forcing_iopt_file
    else
       write(error_msg, "(2A)"), "unknown gas_flux_forcing_opt: ", trim(gas_flux_forcing_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    endif

    fice_file_loc%input        = gas_flux_fice
    xkw_file_loc%input         = gas_flux_ws
    ap_file_loc%input          = gas_flux_ap
    dust_flux_loc%input        = dust_flux_input
    iron_flux_loc%input        = iron_flux_input
    nox_flux_monthly_loc%input = nox_flux_monthly_input
    nhy_flux_monthly_loc%input = nhy_flux_monthly_input
    din_riv_flux_loc%input     = din_riv_flux_input
    dip_riv_flux_loc%input     = dip_riv_flux_input
    don_riv_flux_loc%input     = don_riv_flux_input
    dop_riv_flux_loc%input     = dop_riv_flux_input
    dsi_riv_flux_loc%input     = dsi_riv_flux_input
    dfe_riv_flux_loc%input     = dfe_riv_flux_input
    dic_riv_flux_loc%input     = dic_riv_flux_input
    alk_riv_flux_loc%input     = alk_riv_flux_input
    doc_riv_flux_loc%input     = doc_riv_flux_input

    !-----------------------------------------------------------------------
    !  set variables immediately dependent on namelist variables
    !-----------------------------------------------------------------------

    select case (comp_surf_avg_freq_opt)
    case ('never')
       comp_surf_avg_freq_iopt = marbl_freq_opt_never
    case ('nyear')
       comp_surf_avg_freq_iopt = marbl_freq_opt_nyear
    case ('nmonth')
       comp_surf_avg_freq_iopt = marbl_freq_opt_nmonth
    case default
       write(error_msg, "(2A)"), "unknown comp_surf_avg_freq_opt: ", trim(comp_surf_avg_freq_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    end select

    select case (atm_co2_opt)
    case ('const')
       atm_co2_iopt = atm_co2_iopt_const
    case ('drv_prog')
       atm_co2_iopt = atm_co2_iopt_drv_prog
    case ('drv_diag')
       atm_co2_iopt = atm_co2_iopt_drv_diag
    case default
       write(error_msg, "(2A)"), "unknown atm_co2_opt: ", trim(atm_co2_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    end select

    select case (atm_alt_co2_opt)
    case ('const')
       atm_alt_co2_iopt = atm_co2_iopt_const
    case default
       write(error_msg, "(2A)"), "unknown atm_alt_co2_opt: ", trim(atm_alt_co2_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    end select

    select case (caco3_bury_thres_opt)
    case ('fixed_depth')
       caco3_bury_thres_iopt = caco3_bury_thres_iopt_fixed_depth
    case ('omega_calc')
       caco3_bury_thres_iopt = caco3_bury_thres_iopt_omega_calc
    case default
       write(error_msg, "(2A)"), "unknown caco3_bury_thres_opt: ", trim(caco3_bury_thres_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    end select

    !-----------------------------------------------------------------------
    !  namelist consistency checking
    !-----------------------------------------------------------------------

    if (use_nml_surf_vals .and. comp_surf_avg_freq_iopt /= marbl_freq_opt_never) then
       write(error_msg, "(4A)"), "use_nml_surf_vals can only be .true. if ", &
                                 "comp_surf_avg_freq_opt is 'never', but",   &
                                 "comp_surf_avg_freq_opt = ", trim(comp_surf_avg_freq_opt)
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_nml()")
       return
    endif

    !-----------------------------------------------------------------------
    !  read ecosys_parms_nml namelist
    !-----------------------------------------------------------------------

    ! FIXME(mnl, 2016-01): eliminate marbl_parms!
    call marbl_params_init(nl_buffer, marbl_status_log)
    if (marbl_status_log%labort_marbl) then
       return
    end if

    dust_flux        => dust_flux_loc
    iron_flux        => iron_flux_loc
    fice_file        => fice_file_loc
    xkw_file         => xkw_file_loc
    ap_file          => ap_file_loc
    nox_flux_monthly => nox_flux_monthly_loc
    nhy_flux_monthly => nhy_flux_monthly_loc
    din_riv_flux     => din_riv_flux_loc
    dip_riv_flux     => dip_riv_flux_loc
    don_riv_flux     => don_riv_flux_loc
    dop_riv_flux     => dop_riv_flux_loc
    dsi_riv_flux     => dsi_riv_flux_loc
    dfe_riv_flux     => dfe_riv_flux_loc
    dic_riv_flux     => dic_riv_flux_loc
    alk_riv_flux     => alk_riv_flux_loc
    doc_riv_flux     => doc_riv_flux_loc

  end subroutine marbl_init_nml

  !*****************************************************************************

  subroutine marbl_sflux_forcing_fields_init(num_elements, marbl_forcing_fields)

    ! !DESCRIPTION:
    !  Initialize the sflux forcing_fields datatype with information from the
    !  namelist read
    !
    use marbl_share_mod, only : gas_flux_forcing_iopt_drv
    use marbl_share_mod, only : gas_flux_forcing_iopt_file
    use marbl_share_mod, only : gas_flux_forcing_iopt
    use marbl_share_mod, only : gas_flux_forcing_file
    use marbl_share_mod, only : atm_alt_co2_const
    use marbl_share_mod, only : dust_flux
    use marbl_share_mod, only : iron_flux
    use marbl_share_mod, only : fice_file
    use marbl_share_mod, only : xkw_file
    use marbl_share_mod, only : ap_file
    use marbl_share_mod, only : nox_flux_monthly
    use marbl_share_mod, only : nhy_flux_monthly
    use marbl_share_mod, only : din_riv_flux
    use marbl_share_mod, only : dip_riv_flux
    use marbl_share_mod, only : don_riv_flux
    use marbl_share_mod, only : dop_riv_flux
    use marbl_share_mod, only : dsi_riv_flux
    use marbl_share_mod, only : dfe_riv_flux
    use marbl_share_mod, only : dic_riv_flux
    use marbl_share_mod, only : alk_riv_flux
    use marbl_share_mod, only : doc_riv_flux

    implicit none

    ! !INPUT PARAMETERS:
    integer (KIND=int_kind),         intent(in)   :: num_elements

    ! !OUTPUT PARAMETERS:
    type(marbl_forcing_fields_type), intent(out) :: marbl_forcing_fields

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_sflux_forcing_fields_init'

    character(char_len) :: fsource                  ! JW TODO
    character(char_len) :: drivername               ! JW TODO
    character(char_len) :: filename                 ! JW TODO
    character(char_len) :: varname                  ! JW TODO
    character(char_len) :: units                    ! JW TODO

    real (KIND=r8)      :: constant

    !-----------------------------------------------------------------------
    !  load namelist output into forcing field type
    !-----------------------------------------------------------------------

     ! Allocate memory for surface forcing fields
    call marbl_forcing_fields%construct(num_elements, max_forcing_fields)

    fsource    = 'driver'
    varname    = 'u10_sqr'
    drivername = 'U10_SQR'
    units      = 'unknown'
    call marbl_forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=marbl_forcing_ind%u10_sqr_id)

    associate(                                   &
         ind            => marbl_forcing_ind,    &
         forcing_fields => marbl_forcing_fields  &
         )

    fsource    = 'driver'
    varname    = 'ifrac'
    drivername = 'IFRAC'
    units      = 'dimensionless'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%ifrac_id)

    fsource    = 'driver'
    varname    = 'sst'
    drivername = 'SST'
    units      = 'Temperature (C)'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%sst_id)

    fsource    = 'driver'
    varname    = 'sss'
    drivername = 'SSS'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%sss_id)

    fsource    = 'driver'
    varname    = 'xco2'
    drivername = 'XCO2'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%xco2_id)

    fsource    = 'constant'
    varname    = 'xco2_alt_co2'
    constant   = atm_alt_co2_const
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, field_constant=constant, &
                                          id=ind%xco2_alt_co2_id)

    fsource    = 'driver'
    varname    = 'ph_prev'
    drivername = 'PH_PREV'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%ph_prev_id)

    fsource    = 'driver'
    varname    = 'ph_prev_alt_co2'
    drivername = 'PH_PREV_ALT_CO2'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                          id=ind%ph_prev_alt_co2_id)

    if (gas_flux_forcing_iopt == gas_flux_forcing_iopt_drv) then

       fsource    = 'driver'
       varname    = 'fice'
       drivername = 'FICE_USED'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                             id=ind%fice_id)

       fsource    = 'driver'
       varname    = 'xkw'
       drivername = 'XKW_ICE'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                             id=ind%xkw_id)

       fsource    = 'driver'
       varname    = 'atm_pressure'
       drivername = 'AP_FILE_INPUT'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_driver_varname=drivername, &
                                             id=ind%atm_pressure_id)

    elseif (gas_flux_forcing_iopt == gas_flux_forcing_iopt_file) then

       fsource    = 'POP monthly calendar'
       varname    = 'fice'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=fice_file, &
                                             id=ind%fice_id)

       fsource    = 'POP monthly calendar'
       varname    = 'xkw'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=xkw_file, &
                                             id=ind%xkw_id)

       fsource    = 'POP monthly calendar'
       varname    = 'atm_pressure'
       units      = 'unknown'
       call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=ap_file, &
                                             id=ind%atm_pressure_id)

    endif 

    fsource    = 'POP monthly calendar'
    varname    = 'dust flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dust_flux, &
                                          id=ind%dust_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'iron flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=iron_flux, &
                                          id=ind%iron_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'nox flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=nox_flux_monthly, &
                                          id=ind%nox_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'nhy flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=nhy_flux_monthly, &
                                          id=ind%nhy_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'din river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=din_riv_flux, &
                                          id=ind%din_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'dip river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dip_riv_flux, &
                                          id=ind%dip_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'don river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=don_riv_flux, &
                                          id=ind%don_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'dop river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dop_riv_flux, &
                                          id=ind%dop_riv_flux_id)


    fsource    = 'POP monthly calendar'
    varname    = 'dsi river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dsi_riv_flux, &
                                          id=ind%dsi_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'dfe river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dfe_riv_flux, &
                                          id=ind%dfe_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'dic river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=dic_riv_flux, &
                                          id=ind%dic_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'alk river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=alk_riv_flux, &
                                          id=ind%alk_riv_flux_id)

    fsource    = 'POP monthly calendar'
    varname    = 'doc river flux'
    units      = 'unknown'
    call forcing_fields%add_forcing_field(fsource, varname, units, marbl_forcing_calendar_name=doc_riv_flux, &
                                          id=ind%doc_riv_flux_id)

    end associate

  end subroutine marbl_sflux_forcing_fields_init

  !*****************************************************************************
  
  subroutine marbl_init_tracer_metadata(marbl_tracer_metadata, marbl_status_log)

    ! !DESCRIPTION:
    !  Set tracer and forcing metadata

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type (marbl_tracer_metadata_type), intent(inout) :: marbl_tracer_metadata(:)   ! descriptors for each tracer
    type(marbl_log_type)           , intent(inout) :: marbl_status_log

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_init_tracer_metadata'

    integer (int_kind) :: non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers
    integer (int_kind) :: n        ! index for looping over tracers
    integer (int_kind) :: zoo_ind  ! zooplankton functional group index
    integer (int_kind) :: auto_ind ! autotroph functional group index

    !-----------------------------------------------------------------------
    ! initialize tracer metatdata
    !-----------------------------------------------------------------------

    call marbl_init_forcing_metadata()

    call marbl_init_non_autotroph_tracer_metadata(marbl_tracer_metadata, non_living_biomass_ecosys_tracer_cnt)

    call marbl_check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt, marbl_status_log)
    if (marbl_status_log%labort_marbl) then
      error_msg = "error code returned from marbl_check_ecosys_tracer_count_consistency"
      call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_init_tracer_metadata")
                                    
      return
    end if

    call marbl_initialize_zooplankton_tracer_metadata(marbl_tracer_metadata, non_living_biomass_ecosys_tracer_cnt, n)

    call marbl_initialize_autotroph_tracer_metadata(marbl_tracer_metadata, n)

    !-----------------------------------------------------------------------
    !  set lfull_depth_tavg flag for short-lived ecosystem tracers
    !-----------------------------------------------------------------------

    ! Should be done in ecosys_diagnostics, and without the _tavg name
    do zoo_ind = 1, zooplankton_cnt
       n = zooplankton(zoo_ind)%C_ind
       marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg
    end do

    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%Chl_ind
       marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%C_ind
       marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%Fe_ind
       marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg

       n = autotrophs(auto_ind)%Si_ind
       if (n > 0) then
          marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg
       endif

       n = autotrophs(auto_ind)%CaCO3_ind
       if (n > 0) then
          marbl_tracer_metadata(n)%lfull_depth_tavg = lecovars_full_depth_tavg
       endif
    end do

  end subroutine marbl_init_tracer_metadata

  !***********************************************************************

  subroutine marbl_set_interior_forcing( &
       ciso_on,                          &
       domain,                           &
       gcm_state,                        &
       restore_local ,                   &
       dust_flux_in,                     &
       fesedflux,                        &
       tracer_module,                    &
       marbl_interior_diags,             &
       marbl_restore_diags,              &
       ph_prev_3d,                       &
       ph_prev_alt_co2_3d,               &
       dtracer,                          &
       marbl_interior_share,             &
       marbl_zooplankton_share,          &
       marbl_autotroph_share,            &
       marbl_particulate_share)
    
    ! !DESCRIPTION:
    !  Compute time derivatives for ecosystem state variables

    use ecosys_diagnostics_mod, only : store_diagnostics_carbonate
    use ecosys_diagnostics_mod, only : store_diagnostics_nitrification
    use ecosys_diagnostics_mod, only : store_diagnostics_autotrophs
    use ecosys_diagnostics_mod, only : store_diagnostics_autotroph_sums
    use ecosys_diagnostics_mod, only : store_diagnostics_particulates
    use ecosys_diagnostics_mod, only : store_diagnostics_oxygen
    use ecosys_diagnostics_mod, only : store_diagnostics_PAR
    use ecosys_diagnostics_mod, only : store_diagnostics_misc
    use ecosys_diagnostics_mod, only : store_diagnostics_zooplankton
    use ecosys_diagnostics_mod, only : store_diagnostics_dissolved_organic_matter
    use ecosys_diagnostics_mod, only : store_diagnostics_carbon_fluxes
    use ecosys_diagnostics_mod, only : store_diagnostics_nitrogen_fluxes
    use ecosys_diagnostics_mod, only : store_diagnostics_phosphorus_fluxes
    use ecosys_diagnostics_mod, only : store_diagnostics_silicon_fluxes
    use ecosys_diagnostics_mod, only : store_diagnostics_iron_fluxes

    logical (log_kind)                     , intent(in)    :: ciso_on   ! flag to turn on carbon isotope calculations
    type    (marbl_domain_type)            , intent(in)    :: domain                                
    type    (marbl_gcm_state_type)         , intent(in)    :: gcm_state
    real    (r8)                           , intent(in)    :: restore_local(:,:)    ! (ecosys_used_tracer_cnt, km) local restoring terms for nutrients (mmol ./m^3/sec) 
    real    (r8)                           , intent(in)    :: dust_flux_in
    real    (r8)                           , intent(in)    :: fesedflux(:)
    real    (r8)                           , intent(in)    :: tracer_module(:,: )   ! (ecosys_used_tracer_cnt, km) tracer values 
    real    (r8)                           , intent(inout) :: ph_prev_3d(:)         ! (km)
    real    (r8)                           , intent(inout) :: ph_prev_alt_co2_3d(:) ! (km)
    type    (marbl_diagnostics_type)       , intent(inout) :: marbl_interior_diags
    type    (marbl_diagnostics_type)       , intent(inout) :: marbl_restore_diags
    real    (r8)                           , intent(out)   :: dtracer(:,:)          ! (ecosys_used_tracer_cnt, km) computed source/sink terms
    type    (marbl_interior_share_type)    , intent(inout) :: marbl_interior_share(domain%km)  !FIXME - intent is inout due to DIC_Loc
    type    (marbl_zooplankton_share_type) , intent(inout) :: marbl_zooplankton_share(zooplankton_cnt, domain%km)
    type    (marbl_autotroph_share_type)   , intent(inout) :: marbl_autotroph_share(autotroph_cnt, domain%km)
    type    (marbl_particulate_share_type) , intent(inout) :: marbl_particulate_share

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    character(*), parameter :: subname = 'ecosys_mod:marbl_set_interior_forcing'

    real (r8) :: nitrif(domain%km)    ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
    real (r8) :: denitrif(domain%km)  ! WC nitrification (NO3 -> N2) (mmol N/m^3/sec)

    real (r8) :: O2_production(domain%km)  ! O2 production
    real (r8) :: O2_consumption(domain%km) ! O2 consumption

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

    real (r8) :: Tfunc(domain%km)
    real (r8) :: Fe_scavenge_rate(domain%km) ! annual scavenging rate of iron as % of ambient
    real (r8) :: Fe_scavenge(domain%km)      ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
    real (r8) :: tracer_local(ecosys_tracer_cnt, domain%km)
    real (r8) :: QA_dust_def(domain%km)
    real (r8) :: zsat_calcite(domain%km)     ! Calcite Saturation Depth
    real (r8) :: zsat_aragonite(domain%km)   ! Aragonite Saturation Depth
    real (r8) :: sed_denitrif(domain%km)     ! sedimentary denitrification (nmol N/cm^3/sec)
    real (r8) :: other_remin(domain%km)      ! organic C remin not due oxic or denitrif (nmolC/cm^3/sec)
    real (r8) :: PON_remin(domain%km)        ! remin of PON
    real (r8) :: PON_sed_loss(domain%km)     ! loss of PON to sediments
    real (r8) :: POP_remin(domain%km)        ! remin of POP
    real (r8) :: POP_sed_loss(domain%km)     ! loss of POP to sediments

    type(zooplankton_local_type)             :: zooplankton_local(zooplankton_cnt, domain%km)
    type(autotroph_local_type)               :: autotroph_local(autotroph_cnt, domain%km)
    type(autotroph_secondary_species_type)   :: autotroph_secondary_species(autotroph_cnt, domain%km)
    type(zooplankton_secondary_species_type) :: zooplankton_secondary_species(zooplankton_cnt, domain%km)
    type(dissolved_organic_matter_type)      :: dissolved_organic_matter(domain%km)
    type(carbonate_type)                     :: carbonate(domain%km)

    ! NOTE(bja, 2015-07) vectorization: arrays that are (n, k, c, i)
    ! probably can not be vectorized reasonably over c without memory
    ! copies. If we break up the main k loop, some of the (k, c) loops
    ! can probably be vectorized over k and / or c!

    type(marbl_PAR_type) :: PAR
    !-----------------------------------------------------------------------

    call PAR%construct(num_levels=domain%km, num_PAR_subcols=domain%num_PAR_subcols)

    ! NOTE(bja, 2015-07) dTracer=0 must come before the "not
    ! lsource_sink check to ensure correct answer when not doing
    ! computations.
    ! NOTE(mvertens, 2015-12) the following includes carbon isotopes if 
    ! ciso_on is true

    dtracer(:, :) = c0

    if (.not. lsource_sink) then
       !-----------------------------------------------------------------------
       !  exit immediately if computations are not to be performed
       !-----------------------------------------------------------------------
       return
    endif

    associate(&
         POC     => marbl_particulate_share%POC,     &
         P_CaCO3 => marbl_particulate_share%P_CaCO3, &
         P_SiO2  => marbl_particulate_share%P_SiO2,  &
         dust    => marbl_particulate_share%dust,    &
         P_iron  => marbl_particulate_share%P_iron   &
         )

    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !-----------------------------------------------------------------------

    do k = 1, domain%km
       !write(*, *) 'set_interior forcing loop: ', k, i, c
       call marbl_setup_local_tracers(k, domain%kmt, &
            tracer_module(:, k), tracer_local(:, k))

       call marbl_setup_local_zooplankton(k, domain%kmt, &
            tracer_module(:, k), zooplankton_cnt, zooplankton, zooplankton_local(:, k))

       call marbl_setup_local_autotrophs(k, domain%kmt, &
            tracer_module(:, k), autotroph_cnt, autotrophs, autotroph_local(:, k))

    enddo

    call marbl_init_particulate_terms(1, &
         POC, P_CaCO3, P_SiO2, dust, P_iron, QA_dust_def(:), dust_flux_in)

    !FIXME (mvertens, 2015-11), new marbl timers need to be implemented to turn on timers here
    ! around this subroutine call
    call marbl_compute_carbonate_chemistry(domain, &
         gcm_state%temperature(:), gcm_state%salinity(:), &
         tracer_local(:, :), carbonate(:), &
         ph_prev_3d(:), ph_prev_alt_co2_3d(:), &
         zsat_calcite(:), zsat_aragonite(:))

    call marbl_autotroph_consistency_check(autotroph_cnt, &
         domain%kmt, autotrophs, autotroph_local(:,1:domain%kmt))

    call marbl_compute_PAR(domain, gcm_state, autotroph_cnt, autotroph_local, PAR)

    do k = 1, domain%km

       call marbl_compute_autotroph_elemental_ratios( autotroph_cnt,    &
            autotrophs, autotroph_local(:, k), tracer_local(:, k),      &
            autotroph_secondary_species(:, k))

       call marbl_compute_function_scaling(gcm_state%temperature(k), Tfunc(k))

       call marbl_compute_Pprime(k, domain, autotroph_cnt, autotrophs, &
            autotroph_local(:, k), gcm_state%temperature(k), autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_uptake(autotroph_cnt, autotrophs, &
            tracer_local(:, k), &
            autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_photosynthesis(autotroph_cnt,      &
            domain%num_PAR_subcols, autotrophs, autotroph_local(:, k), &
            gcm_state%temperature(k), Tfunc(k), PAR%col_frac(:),       &
            PAR%avg(k,:), autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_phyto_diatoms (autotroph_cnt, &
            autotrophs, autotroph_local(:, k),                    &
            autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_calcification(autotroph_cnt, autotrophs, &
            autotroph_local(:, k),  gcm_state%temperature(k), autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_nfixation(autotroph_cnt, autotrophs, &
            autotroph_secondary_species(:, k))

       call marbl_compute_autotroph_loss(autotroph_cnt, autotrophs, &
            Tfunc(k), autotroph_secondary_species(:, k))

       call marbl_compute_Zprime(k, domain, &
            zooplankton_cnt, zooplankton, zooplankton_local(:, k)%C, &
            Tfunc(k), zooplankton_secondary_species(:, k))

       call marbl_compute_grazing (autotroph_cnt, zooplankton_cnt, grazer_prey_cnt, autotrophs, &
            Tfunc(k), zooplankton_local(:, k), &
            zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

       call marbl_compute_routing (autotroph_cnt, zooplankton_cnt, autotrophs, &
            zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k))

       call marbl_compute_dissolved_organic_matter (k, autotroph_cnt, zooplankton_cnt, &
            domain%num_PAR_subcols, autotrophs,        &
            zooplankton_secondary_species(:, k),                        &
            autotroph_secondary_species(:, k),                          &
            PAR%col_frac(:), PAR%interface(k-1,:), PAR%avg(k,:),        &
            domain%delta_z(1), tracer_local(:, k),                      &
            dissolved_organic_matter(k))

       call marbl_compute_large_detritus(k, autotroph_cnt, zooplankton_cnt, autotrophs, &
            zooplankton_secondary_species(:, k), autotroph_secondary_species(:, k), tracer_local(fe_ind, k), &
            POC, P_CaCO3, P_SiO2, dust, P_iron, &
            Fe_scavenge(k), Fe_scavenge_rate(k))

       ! FIXME(bja, 2015-08) need to pull particulate share out of compute_particulate_terms!
       call marbl_compute_particulate_terms(k, domain,                 &
            marbl_particulate_share, POC, P_CaCO3, P_SiO2, dust,       &
            P_iron, PON_remin(k), PON_sed_loss(k), POP_remin(k),       &
            POP_sed_loss(k), QA_dust_def(k), gcm_state%temperature(k), &
            tracer_local(:, k), carbonate(k), sed_denitrif(k),         &
            other_remin(k), fesedflux(k), ciso_on)

       call marbl_compute_nitrif(k, domain%num_PAR_subcols, domain%kmt, &
            PAR%col_frac(:), PAR%interface(k-1,:), PAR%interface(k,:),  &
            PAR%KPARdz(k), tracer_local(nh4_ind, k), nitrif(k))

       call marbl_compute_denitrif(tracer_local(o2_ind, k), tracer_local(no3_ind, k), &
            dissolved_organic_matter(k)%DOC_remin, &
            dissolved_organic_matter(k)%DOCr_remin, &
            POC%remin(k), other_remin(k), sed_denitrif(k), denitrif(k))

       call marbl_compute_dtracer_local (autotroph_cnt, zooplankton_cnt, autotrophs, zooplankton, &
            autotroph_secondary_species(:, k), &
            zooplankton_secondary_species(:, k), &
            dissolved_organic_matter(k), &
            nitrif(k), denitrif(k), sed_denitrif(k), &
            Fe_scavenge(k) , Fe_scavenge_rate(k), &
            P_iron%remin(k), POC%remin(k), &
            P_SiO2%remin(k), P_CaCO3%remin(k), other_remin(k), &
            PON_remin(k), POP_remin(k), &
            restore_local(:, k), &
            tracer_local(o2_ind, k), &
            o2_production(k), o2_consumption(k), &
            dtracer(:, k) )

       if (ciso_on) then
          ! FIXME(bja, 2015-08) need to pull particulate share
          ! out of compute_particulate_terms!
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
       end if

       if  (k < domain%km) then
          call marbl_update_particulate_terms_from_prior_level(k+1, POC, P_CaCO3, &
               P_SiO2, dust, P_iron, QA_dust_def(:))
       endif

    end do ! k

    ! FIXME(mnl,2016-01) call store_diagnostics_interior() to call each individually
    !                    and call set_to_zero from that new routine
    call marbl_interior_diags%set_to_zero()

    call store_diagnostics_carbonate(domain, &
         carbonate, marbl_interior_diags)

    call store_diagnostics_autotrophs(domain, &
         autotroph_secondary_species, marbl_interior_diags)

    call store_diagnostics_particulates(domain, &
         POC, P_CaCO3, P_SiO2, dust,  P_iron, PON_remin, PON_sed_loss, POP_remin,  &
         POP_sed_loss, sed_denitrif, other_remin, marbl_interior_diags)

    call store_diagnostics_autotroph_sums(domain, &
         autotroph_secondary_species, marbl_interior_diags)

    call store_diagnostics_nitrification(nitrif, denitrif, marbl_interior_diags)

    call store_diagnostics_oxygen(domain, gcm_state, &
         tracer_module(o2_ind, :), o2_production, o2_consumption, marbl_interior_diags)

    call store_diagnostics_PAR(domain, PAR%col_frac(:), PAR%avg(:,:), marbl_interior_diags)

    call store_diagnostics_zooplankton(zooplankton_secondary_species, marbl_interior_diags)

    call store_diagnostics_dissolved_organic_matter(domain, &
         dissolved_organic_matter, fe_scavenge, fe_scavenge_rate, marbl_interior_diags)

    call store_diagnostics_carbon_fluxes(domain, &
         POC, P_CaCO3, dtracer, marbl_interior_diags)

    call store_diagnostics_nitrogen_fluxes(domain, &
         PON_sed_loss, denitrif, sed_denitrif, autotroph_secondary_species, dtracer, &
         marbl_interior_diags)

    call store_diagnostics_phosphorus_fluxes(domain, &
         POP_sed_loss, dtracer, marbl_interior_diags)

    call store_diagnostics_silicon_fluxes(domain, &
         P_SiO2, dtracer, marbl_interior_diags)

    call store_diagnostics_iron_fluxes(domain, &
         P_iron, dust, fesedflux, dtracer, marbl_interior_diags)

    ! store_diagnostics_restore
    do n = 1, ecosys_tracer_cnt
       marbl_restore_diags%diags(n)%field_3d(:,1) = restore_local(n,:)
    end do

    end associate

    call PAR%destruct()

  end subroutine marbl_set_interior_forcing

  !***********************************************************************

  subroutine marbl_init_particulate_terms(k, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, QA_dust_def, NET_DUST_IN)

    ! !DESCRIPTION:
    !  Set incoming fluxes (put into outgoing flux for first level usage).
    !  Set dissolution length, production fraction and mass terms.
    !
    !  The first 6 arguments are intent(inout) in
    !  order to preserve contents on other blocks.

    use marbl_share_mod, only : dust_flux        

    integer(int_kind)                  , intent(in)    :: k
    real (r8)                          , intent(in)    :: net_dust_in     ! dust flux
    type(column_sinking_particle_type) , intent(inout) :: POC             ! base units = nmol C
    type(column_sinking_particle_type) , intent(inout) :: P_CaCO3         ! base units = nmol CaCO3
    type(column_sinking_particle_type) , intent(inout) :: P_SiO2          ! base units = nmol SiO2
    type(column_sinking_particle_type) , intent(inout) :: dust            ! base units = g
    type(column_sinking_particle_type) , intent(inout) :: P_iron          ! base units = nmol Fe
    real (r8)                          , intent(inout) :: QA_dust_def(:)  ! incoming deficit in the QA(dust) POC flux (km)

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
    P_CaCO3%gamma = 0.10_r8         ! prod frac -> hard subclass
    P_CaCO3%mass  = 100.09_r8       ! molecular weight of CaCO
    P_CaCO3%rho   = 0.05_r8 * P_CaCO3%mass / POC%mass ! QA mass ratio for CaCO3

    P_SiO2%diss   = parm_SiO2_diss  ! diss. length (cm), modified by TEMP
    P_SiO2%gamma  = 0.10_r8         ! prod frac -> hard subclass
    P_SiO2%mass   = 60.08_r8        ! molecular weight of SiO2
    P_SiO2%rho    = 0.05_r8 * P_SiO2%mass / POC%mass ! QA mass ratio for SiO2

    dust%diss     = 30000.0_r8      ! diss. length (cm)
    dust%gamma    = 0.99_r8         ! prod frac -> hard subclass
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

    QA_dust_def(k) = dust%rho * (dust%sflux_out(k) + dust%hflux_out(k))

  end subroutine marbl_init_particulate_terms

  !***********************************************************************

  subroutine marbl_update_particulate_terms_from_prior_level(k, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, QA_dust_def)

    integer (int_kind)                 , intent(in)    :: k ! vertical model level
    type(column_sinking_particle_type) , intent(inout) :: POC, P_CaCO3, P_SiO2, dust, P_iron
    real(r8)                           , intent(inout) :: QA_dust_def(:) !(km)

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

    integer (int_kind), intent(in) :: k
    type(column_sinking_particle_type), intent(inout) :: sinking_particle

    ! NOTE(bja, 201504) level k influx is equal to the level k-1 outflux.
    sinking_particle%sflux_out(k) = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_out(k) = sinking_particle%hflux_out(k-1)
    sinking_particle%sflux_in(k)  = sinking_particle%sflux_out(k-1)
    sinking_particle%hflux_in(k)  = sinking_particle%hflux_out(k-1)

  end subroutine marbl_update_sinking_particle_from_prior_level

  !***********************************************************************

  subroutine marbl_compute_particulate_terms(k, domain,                       &
             marbl_particulate_share, POC, P_CaCO3, P_SiO2, dust, P_iron,     &
             PON_remin, PON_sed_loss, POP_remin, POP_sed_loss, QA_dust_def,   &
             temperature, tracer_local, carbonate, sed_denitrif, other_remin, &
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

    use marbl_parms           , only : Tref

    ! !INPUT PARAMETERS:
    integer (int_kind)                      , intent(in)    :: k                   ! vertical model level
    type(marbl_domain_type)                 , intent(in)    :: domain                              
    real (r8)                               , intent(in)    :: temperature         ! temperature for scaling functions bsi%diss
    real (r8), dimension(ecosys_tracer_cnt) , intent(in)    :: tracer_local        ! local copies of model tracer concentrations
    type(carbonate_type)                    , intent(in)    :: carbonate
    logical (log_kind)                      , intent(in)    :: lexport_shared_vars ! flag to save shared_vars or not
    real(r8)                                , intent(in)    :: fesedflux           ! sedimentary Fe input

    ! !OUTPUT PARAMETERS:
    real(r8)                                , intent(out)   :: PON_remin           ! remin of PON
    real(r8)                                , intent(out)   :: PON_sed_loss        ! loss of PON to sediments

    ! !INPUT/OUTPUT PARAMETERS:
    type(column_sinking_particle_type)      , intent(inout) :: POC                 ! base units = nmol C
    type(column_sinking_particle_type)      , intent(inout) :: P_CaCO3             ! base units = nmol CaCO3
    type(column_sinking_particle_type)      , intent(inout) :: P_SiO2              ! base units = nmol SiO2
    type(column_sinking_particle_type)      , intent(inout) :: dust                ! base units = g
    type(column_sinking_particle_type)      , intent(inout) :: P_iron              ! base units = nmol Fe
    real (r8)                               , intent(out)   :: POP_remin           ! remin of POP
    real (r8)                               , intent(out)   :: POP_sed_loss        ! loss of POP to sediments
    real (r8)                               , intent(inout) :: QA_dust_def         ! incoming deficit in the QA(dust) POC flux
    real (r8)                               , intent(out)   :: sed_denitrif        ! sedimentary denitrification (umolN/cm^2/s)
    real (r8)                               , intent(out)   :: other_remin         ! sedimentary remin not due to oxic or denitrification
    type(marbl_particulate_share_type)      , intent(inout) :: marbl_particulate_share

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
         scalelength,        & ! used to scale dissolution length scales as function of depth
         o2_scalefactor,     & ! used to scale dissolution length scales as function of o2
         flux, flux_alt,     & ! temp variables used to update sinking flux
         dz_loc, dzr_loc       ! dz, dzr at a particular i, j location

    real (r8), parameter :: &  ! o2_sf is an abbreviation for o2_scalefactor
         o2_sf_o2_range_hi = 50.0_r8, & ! apply o2_scalefactor for O2_loc less than this
         o2_sf_o2_range_lo =  5.0_r8, & ! o2_scalefactor is constant for O2_loc < this parameter
         o2_sf_val_lo_o2   =  2.5_r8    ! o2_scalefactor for O2_loc < o2_sf_o2_range_lo

    integer (int_kind) :: n     ! loop indices

    logical (log_kind) :: poc_error   ! POC error flag
    !-----------------------------------------------------------------------

    associate(                                                                         &
         column_kmt               => domain%kmt,                                       &
         delta_z                  => domain%delta_z,                                   &
         zw                       => domain%zw,                                        & 
         O2_loc                   => tracer_local(o2_ind),                             &
         NO3_loc                  => tracer_local(no3_ind),                            &
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

    DECAY_Hard     = exp(-delta_z(k) / 4.0e6_r8)
    DECAY_HardDust = exp(-delta_z(k) / 1.2e8_r8)

    !----------------------------------------------------------------------
    !   Tref = 30.0 reference temperature (deg. C)
    !-----------------------------------------------------------------------
    TfuncS = 1.5_r8**(((temperature + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

    poc_error = .false.
    dz_loc = delta_z(k)

    if (k <= column_kmt) then

       dzr_loc    = c1 / dz_loc
       poc_diss   = POC%diss
       sio2_diss  = P_SiO2%diss
       caco3_diss = P_CaCO3%diss
       dust_diss  = dust%diss

       !-----------------------------------------------------------------------
       !  increase POC diss length scale where O2 concentrations are low
       !-----------------------------------------------------------------------

       if (O2_loc < o2_sf_o2_range_hi) then
          o2_scalefactor = c1 + (o2_sf_val_lo_o2 - c1) * &
               min(c1, (o2_sf_o2_range_hi - O2_loc)/(o2_sf_o2_range_hi - o2_sf_o2_range_lo))
          poc_diss   = poc_diss   * o2_scalefactor
          sio2_diss  = sio2_diss  * o2_scalefactor
          caco3_diss = caco3_diss * o2_scalefactor
          dust_diss  = dust_diss  * o2_scalefactor
       endif

       !-----------------------------------------------------------------------
       !  apply scalelength factor to length scales
       !-----------------------------------------------------------------------

       poc_diss = scalelength * poc_diss
       sio2_diss = scalelength * sio2_diss
       caco3_diss = scalelength * caco3_diss
       dust_diss = scalelength * dust_diss

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

       PON_remin = Q * POC%remin(k)

       POP_remin = Qp_zoo_pom * POC%remin(k)

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

       ! add term for desorption of iron from sinking particles
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

       PON_remin = c0

       POP_remin = c0

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
    !  Calcite is preserved in sediments above a threshold depth, 
    !     which is based on caco3_bury_thres_opt.
    !-----------------------------------------------------------------------

    POC%sed_loss(k)     = c0
    P_SiO2%sed_loss(k)  = c0
    P_CaCO3%sed_loss(k) = c0
    P_iron%sed_loss(k)  = c0
    dust%sed_loss(k)    = c0

    PON_sed_loss        = c0

    POP_sed_loss        = c0

    if ((k == column_kmt)) then

       flux = POC%sflux_out(k) + POC%hflux_out(k)

       if (flux > c0) then
          flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day

          POC%sed_loss(k) = flux * min(0.8_r8, parm_POMbury &
               * (0.013_r8 + 0.53_r8 * flux_alt*flux_alt / (7.0_r8 + flux_alt)**2))

          PON_sed_loss = PON_bury_coeff * Q * POC%sed_loss(k)

          POP_sed_loss = POP_bury_coeff * Qp_zoo_pom * POC%sed_loss(k)

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

       if (caco3_bury_thres_iopt == caco3_bury_thres_iopt_fixed_depth) then
          if (zw(k) < caco3_bury_thres_depth) then
             P_CaCO3%sed_loss(k) = P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)
          endif
       else ! caco3_bury_thres_iopt = caco3_bury_thres_iopt_omega_calc
          if (carbonate%CO3 > carbonate%CO3_sat_calcite) then
             P_CaCO3%sed_loss(k) = P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)
          endif
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

          PON_remin = PON_remin &
               + ((Q * flux - PON_sed_loss) * dzr_loc)

          POP_remin = POP_remin &
               + ((Qp_zoo_pom * flux - POP_sed_loss) * dzr_loc)
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
      ! FIXME: use marbl_status_log
       call shr_sys_abort(subname /&
            &/ ': mass ratio of ballast ' /&
            &/ 'production exceeds POC production')
    endif

    end associate

  end subroutine marbl_compute_particulate_terms

  !***********************************************************************

  subroutine marbl_set_surface_forcing( &
       ciso_on,               &
       num_elements,          &
       marbl_forcing_input,   &
       marbl_forcing_output,  &
       marbl_forcing_share,   &
       marbl_forcing_diags)

    ! !DESCRIPTION:
    !  Compute surface forcing fluxes 

    use co2calc_column        , only : co2calc_surf
    use co2calc_column        , only : thermodynamic_coefficients_type
    use schmidt_number        , only : schmidt_co2_surf
    use marbl_oxygen          , only : schmidt_o2_surf
    use marbl_oxygen          , only : o2sat_surf
    use marbl_share_mod       , only : lflux_gas_o2
    use marbl_share_mod       , only : lflux_gas_co2
    use marbl_share_mod       , only : ndep_data_type
    use marbl_share_mod       , only : ndep_shr_stream_scale_factor
    use marbl_share_mod       , only : gas_flux_forcing_iopt_drv
    use marbl_share_mod       , only : gas_flux_forcing_iopt_file
    use marbl_share_mod       , only : gas_flux_forcing_iopt
    use marbl_share_mod       , only : nox_flux_monthly 
    use marbl_share_mod       , only : nhy_flux_monthly 
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
    use ecosys_diagnostics_mod, only : store_diagnostics_sflux

    ! !INPUT PARAMETERS:
    integer (int_kind)              , intent(in) :: num_elements
    logical (log_kind)              , intent(in) :: ciso_on ! flag to save shared_vars or not
    type(marbl_forcing_input_type)  , intent(in) :: marbl_forcing_input

    ! !INPUT/OUTPUT PARAMETERS:
    type(marbl_forcing_output_type) , intent(inout) :: marbl_forcing_output
    type(marbl_diagnostics_type),     intent(inout) :: marbl_forcing_diags
    type(marbl_forcing_share_type)  , intent(inout) :: marbl_forcing_share

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(*), parameter :: subname = 'ecosys_mod:marbl_set_surface_forcing'

    integer (int_kind) :: &
         n,               & ! loop indices
         auto_ind           ! autotroph functional group index

    real (r8), dimension(num_elements) :: &
         PHLO,         & ! lower bound for ph in solver
         PHHI,         & ! upper bound for ph in solver
         PH_NEW,       & ! computed PH from solver
         XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
         O2SAT_1atm      ! O2 saturation @ 1 atm (mmol/m^3)

    type(thermodynamic_coefficients_type), dimension(num_elements) :: co3_coeffs
    !-----------------------------------------------------------------------

    associate(                                                              &
         land_mask            => marbl_forcing_input%land_mask            , & 
         ifrac                => marbl_forcing_input%ifrac                , &
         u10_sqr              => marbl_forcing_input%u10_sqr              , &
         sst                  => marbl_forcing_input%sst                  , & ! sea surface temperature (c)
         sss                  => marbl_forcing_input%sss                  , & ! sea surface salinity (psu)
         xco2                 => marbl_forcing_input%xco2                 , & ! atmospheric co2 conc. (dry-air, 1 atm)
         xco2_alt_co2         => marbl_forcing_input%xco2_alt_co2         , & ! atmospheric alternative CO2 (dry-air, 1 atm)
         ap_used              => marbl_forcing_input%atm_press            , & ! used atm pressure (atm)
         xkw                  => marbl_forcing_input%xkw                  , &
         dust_flux_in         => marbl_forcing_input%dust_flux            , &
         iron_flux_in         => marbl_forcing_input%iron_flux            , &
         surface_vals         => marbl_forcing_input%surface_vals         , & 
         input_forcings       => marbl_forcing_input%input_forcings            , & 
         ph_prev              => marbl_forcing_input%ph_prev              , &
         ph_prev_alt_co2      => marbl_forcing_input%ph_prev_alt_co2      , &

         ph_prev_new          => marbl_forcing_output%ph_prev             , &
         ph_prev_alt_co2_new  => marbl_forcing_output%ph_prev_alt_co2     , &
         iron_flux_in_new     => marbl_forcing_output%iron_flux           , &
         flux_co2             => marbl_forcing_output%flux_co2            , &
         flux_alt_co2         => marbl_forcing_output%flux_alt_co2        , & ! (used by store_sflux)
         flux_o2              => marbl_forcing_output%flux_o2             , & ! (used by store_sflux)
         co2star              => marbl_forcing_output%co2star             , & ! (used by store_sflux)
         dco2star             => marbl_forcing_output%dco2star            , & ! (used by store_sflux)
         pco2surf             => marbl_forcing_output%pco2surf            , & ! (used by store_sflux)
         dpco2                => marbl_forcing_output%dpco2               , & ! (used by store_sflux)
         co3                  => marbl_forcing_output%co3                 , & ! (used by store_sflux)
         co2star_alt          => marbl_forcing_output%co2star_alt         , & ! (used by store_sflux)
         dco2star_alt         => marbl_forcing_output%dco2star_alt        , & ! (used by store_sflux)
         pco2surf_alt         => marbl_forcing_output%pco2surf_alt        , & ! (used by store_sflux)
         dpco2_alt            => marbl_forcing_output%dpco2_alt           , & ! (used by store_sflux)
         schmidt_co2          => marbl_forcing_output%schmidt_co2         , & ! (used by store_sflux) used schmidt number 
         schmidt_o2           => marbl_forcing_output%schmidt_o2          , & ! (used by store_sflux) used schmidt number 
         pv_o2                => marbl_forcing_output%pv_o2               , & ! (used by store_sflux) piston velocity (cm/s) 
         pv_co2               => marbl_forcing_output%pv_co2              , & ! (used by store_sflux) piston velocity (cm/s) 
         o2sat                => marbl_forcing_output%o2sat               , & ! (used by store_sflux) used O2 saturation (mmol/m^3) 
         stf_module           => marbl_forcing_output%stf_module(:,:)     , & !

         PV_SURF_fields       => marbl_forcing_share%PV_SURF_fields       , & ! IN/OUT
         DIC_SURF_fields      => marbl_forcing_share%DIC_SURF_fields      , & ! IN/OUT
         CO2STAR_SURF_fields  => marbl_forcing_share%CO2STAR_SURF_fields  , & ! IN/OUT
         DCO2STAR_SURF_fields => marbl_forcing_share%DCO2STAR_SURF_fields , & ! IN/OUT
         CO3_SURF_fields      => marbl_forcing_share%CO3_SURF_fields      , & ! IN/OUT
         dic_riv_flux_fields  => marbl_forcing_share%dic_riv_flux_fields  , & ! IN/OUT
         doc_riv_flux_fields  => marbl_forcing_share%doc_riv_flux_fields    & ! IN/OUT
         )

    !-----------------------------------------------------------------------
    !  fluxes initially set to 0
    !-----------------------------------------------------------------------

    STF_MODULE(:, :) = c0

    !-----------------------------------------------------------------------
    !  calculate gas flux quantities if necessary
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  compute CO2 flux, computing disequilibrium one row at a time
    !-----------------------------------------------------------------------
       
    if (lflux_gas_o2 .or. lflux_gas_co2) then

       !-----------------------------------------------------------------------
       !  Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
       !-----------------------------------------------------------------------

       xkw_ice(:) = (c1 - ifrac(:)) * xkw

       !-----------------------------------------------------------------------
       !  compute O2 flux
       !-----------------------------------------------------------------------

       if (lflux_gas_o2) then
          schmidt_o2(:) = schmidt_o2_surf(num_elements, sst, land_mask)

          o2sat_1atm(:) = o2sat_surf(num_elements, sst, sss, land_mask)

          where (land_mask)
             pv_o2(:) = xkw_ice(:) * sqrt(660.0_r8 / schmidt_o2(:))
             o2sat(:) = ap_used(:) * o2sat_1atm(:)
             flux_o2(:) = pv_o2(:) * (o2sat(:) - surface_vals(:, o2_ind))
             stf_module(:, o2_ind) = stf_module(:, o2_ind) + flux_o2(:)
          elsewhere
             pv_o2(:) = c0
             o2sat(:) = c0
          end where
       else
          schmidt_o2(:) = c0
          pv_o2(:)      = c0
          o2sat(:)      = c0
       endif  ! lflux_gas_o2

       !-----------------------------------------------------------------------
       !  compute CO2 flux, computing disequilibrium one row at a time
       !-----------------------------------------------------------------------

       if (lflux_gas_co2) then

          SCHMIDT_CO2 = SCHMIDT_CO2_surf(num_elements, SST, land_mask)

          where (land_mask)
             PV_CO2 = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_CO2)
          elsewhere
             PV_CO2 = c0
          end where

          !-----------------------------------------------------------------------
          !  Set FLUX_CO2
          !-----------------------------------------------------------------------

          where (PH_PREV /= c0)
             PHLO = PH_PREV - del_ph
             PHHI = PH_PREV + del_ph
          elsewhere
             PHLO = phlo_surf_init
             PHHI = phhi_surf_init
          end where

          call co2calc_surf(num_elements, land_mask, .true., co3_coeffs, SST, SSS, &
                            surface_vals(:,dic_ind), surface_vals(:,alk_ind), &
                            surface_vals(:,po4_ind), surface_vals(:,sio3_ind), &
                            PHLO, PHHI, PH_NEW, XCO2, AP_USED, &
                            CO2STAR, DCO2STAR, pCO2SURF, DpCO2, CO3)

          PH_PREV_NEW = PH_NEW

          FLUX_CO2 = PV_CO2 * DCO2STAR
 
          !-------------------------------------------------------------------
          !  The following variables need to be shared with other modules,
          !  and are now defined in marbl_share as targets.
          !-------------------------------------------------------------------

          if (ciso_on) then
             PV_SURF_fields       = PV_CO2(:)
             DIC_SURF_fields      = surface_vals(:,dic_ind)
             CO2STAR_SURF_fields  = CO2STAR(:)
             DCO2STAR_SURF_fields = DCO2STAR(:)
             CO3_SURF_fields      = CO3(:)
          endif

          !-----------------------------------------------------------------------
          !  Set FLUX_ALT_CO2
          !-----------------------------------------------------------------------

          where (PH_PREV_ALT_CO2 /= c0)
             PHLO = PH_PREV_ALT_CO2 - del_ph
             PHHI = PH_PREV_ALT_CO2 + del_ph
          elsewhere
             PHLO = phlo_surf_init
             PHHI = phhi_surf_init
          end where

          call co2calc_surf(num_elements, land_mask, .false., co3_coeffs, SST, SSS,   &
                            surface_vals(:,dic_alt_co2_ind), surface_vals(:,alk_ind), &
                            surface_vals(:,po4_ind)        , surface_vals(:,sio3_ind),&
                            PHLO, PHHI, PH_NEW, XCO2_ALT_CO2, AP_USED,          &
                            CO2STAR_ALT, DCO2STAR_ALT, pCO2SURF_ALT, DpCO2_ALT, CO3)

          PH_PREV_ALT_CO2_NEW = PH_NEW

          FLUX_ALT_CO2    = PV_CO2 * DCO2STAR_ALT

          !-----------------------------------------------------------------------
          !  set air-sea co2 gas flux named field, converting units from
          !  nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
          !-----------------------------------------------------------------------

          STF_MODULE(:, dic_ind)         = STF_MODULE(:, dic_ind)         + FLUX_CO2(:)
          STF_MODULE(:, dic_alt_co2_ind) = STF_MODULE(:, dic_alt_co2_ind) + FLUX_ALT_CO2(:)

       else
          SCHMIDT_CO2(:) = c0
          PV_CO2(:)           = c0
       endif  !  lflux_gas_co2

    endif  ! lflux_gas_o2 .or. lflux_gas_co2

    !-----------------------------------------------------------------------
    !  calculate iron and dust fluxes if necessary
    !-----------------------------------------------------------------------

    IRON_FLUX_IN_NEW(:) = IRON_FLUX_IN(:) * parm_Fe_bioavail  ! TODO: this gets moved up and out - a forcing field modify

    STF_MODULE(:, fe_ind) = STF_MODULE(:, fe_ind) + IRON_FLUX_IN_NEW(:)

    !-----------------------------------------------------------------------
    !  Add phosphate and silicate from dust after Krishnamurthy et al. (2010)
    !  factors convert from g/cm2/s to nmol/cm2/s
    !  ( P frac in dust by weight) * ( P solubility) / ( P molecular weight) * (mol->nmol)
    !  (Si frac in dust by weight) * (Si solubility) / (Si molecular weight) * (mol->nmol)
    !-----------------------------------------------------------------------

    STF_MODULE(:, po4_ind) = STF_MODULE(:, po4_ind)                           &
       + (dust_flux_in * (0.00105_r8 *  0.15_r8 / 30.974_r8 * 1.0e9_r8))
    STF_MODULE(:, sio3_ind) = STF_MODULE(:, sio3_ind)                         &
       + (dust_flux_in * (  0.308_r8 * 0.075_r8 / 28.085_r8 * 1.0e9_r8))

    !-----------------------------------------------------------------------
    !  calculate nox and nhy fluxes if necessary
    !-----------------------------------------------------------------------
       
    if (nox_flux_monthly%has_data) then
       STF_MODULE(:, no3_ind) = STF_MODULE(:, no3_ind) + input_forcings(:, ind_nox_flux)
    endif
       
    if (nhy_flux_monthly%has_data) then
       STF_MODULE(:, nh4_ind) = STF_MODULE(:, nh4_ind) + input_forcings(:, ind_nhy_flux)
    endif
       
    if (trim(ndep_data_type) == 'shr_stream') then
       where (land_mask(:))
          STF_MODULE(:, no3_ind) = STF_MODULE(:, no3_ind) &
               + ndep_shr_stream_scale_factor * input_forcings(:, ind_no3_flux)
          
          STF_MODULE(:, nh4_ind) = STF_MODULE(:, nh4_ind) &
               + ndep_shr_stream_scale_factor * input_forcings(:, ind_nh4_flux)
       endwhere
    endif

    !-----------------------------------------------------------------------
    !  calculate river bgc fluxes if necessary
    !-----------------------------------------------------------------------

    if (din_riv_flux%has_data) then
       STF_MODULE(:, no3_ind) = STF_MODULE(:, no3_ind) + input_forcings(:, ind_din_riv_flux)
    endif

    if (dip_riv_flux%has_data) then
       STF_MODULE(:, po4_ind) = STF_MODULE(:, po4_ind) + input_forcings(:, ind_dip_riv_flux)
    endif

    if (don_riv_flux%has_data) then
       STF_MODULE(:, don_ind)  = STF_MODULE(:, don_ind)  +                    &
            (input_forcings(:, ind_don_riv_flux) * (c1 - DONriv_refract))
       STF_MODULE(:, donr_ind) = STF_MODULE(:, donr_ind) +                    &
            (input_forcings(:, ind_don_riv_flux) * DONriv_refract)
    endif

    if (dop_riv_flux%has_data) then
       STF_MODULE(:, dop_ind)  = STF_MODULE(:, dop_ind)  +                    &
             (input_forcings(:, ind_dop_riv_flux) * (c1 - DOPriv_refract))
       STF_MODULE(:, dopr_ind) = STF_MODULE(:, dopr_ind) +                    &
             (input_forcings(:, ind_dop_riv_flux) * DOPriv_refract)
    endif

    if (dsi_riv_flux%has_data) then
       STF_MODULE(:, sio3_ind) = STF_MODULE(:, sio3_ind) + input_forcings(:, ind_dsi_riv_flux)
    endif

    if (dfe_riv_flux%has_data) then
       STF_MODULE(:, fe_ind) = STF_MODULE(:, fe_ind) + input_forcings(:, ind_dfe_riv_flux)
    endif

    if (dic_riv_flux%has_data) then
       STF_MODULE(:, dic_ind)         = STF_MODULE(:, dic_ind)         + input_forcings(:, ind_dic_riv_flux)
       STF_MODULE(:, dic_alt_co2_ind) = STF_MODULE(:, dic_alt_co2_ind) + input_forcings(:, ind_dic_riv_flux)
       if (ciso_on) dic_riv_flux_fields = input_forcings(:, ind_dic_riv_flux)
    endif

    if (alk_riv_flux%has_data) then
       STF_MODULE(:, alk_ind) = STF_MODULE(:, alk_ind) + input_forcings(:, ind_alk_riv_flux)
    endif

    if (doc_riv_flux%has_data) then
       STF_MODULE(:, doc_ind) = STF_MODULE(:, doc_ind) +                      &
              (input_forcings(:, ind_doc_riv_flux) * (c1 - DOCriv_refract))
       STF_MODULE(:, docr_ind) = STF_MODULE(:, docr_ind) +                    &
              (input_forcings(:, ind_doc_riv_flux) * DOCriv_refract)

       ! FIXME(ktl) sending total doc river input to ciso for now, need to separate doc and docr
       if (ciso_on) doc_riv_flux_fields = input_forcings(:, ind_doc_riv_flux)
    endif

    !-----------------------------------------------------------------------
    !  Apply NO & NH fluxes to alkalinity
    !-----------------------------------------------------------------------

    STF_MODULE(:, alk_ind) = STF_MODULE(:, alk_ind) + STF_MODULE(:, nh4_ind) - STF_MODULE(:, no3_ind)

    end associate

    call store_diagnostics_sflux(marbl_forcing_input, marbl_forcing_output,   &
         marbl_forcing_diags)

  end subroutine marbl_set_surface_forcing

  !***********************************************************************

  subroutine marbl_init_forcing_metadata()

    !-----------------------------------------------------------------------
    ! initialize surface forcing metadata
    !-----------------------------------------------------------------------

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

    call marbl_init_monthly_forcing_metadata(dust_flux)
    call marbl_init_monthly_forcing_metadata(iron_flux)
    call marbl_init_monthly_forcing_metadata(fice_file)
    call marbl_init_monthly_forcing_metadata(xkw_file)
    call marbl_init_monthly_forcing_metadata(ap_file)
    call marbl_init_monthly_forcing_metadata(nox_flux_monthly)
    call marbl_init_monthly_forcing_metadata(nhy_flux_monthly)
    call marbl_init_monthly_forcing_metadata(din_riv_flux)
    call marbl_init_monthly_forcing_metadata(dip_riv_flux)
    call marbl_init_monthly_forcing_metadata(don_riv_flux)
    call marbl_init_monthly_forcing_metadata(dop_riv_flux)
    call marbl_init_monthly_forcing_metadata(dsi_riv_flux)
    call marbl_init_monthly_forcing_metadata(dfe_riv_flux)
    call marbl_init_monthly_forcing_metadata(dic_riv_flux)
    call marbl_init_monthly_forcing_metadata(alk_riv_flux)
    call marbl_init_monthly_forcing_metadata(doc_riv_flux)

  end subroutine marbl_init_forcing_metadata

  !*****************************************************************************

  subroutine marbl_init_monthly_forcing_metadata(var)

    implicit none

    type(forcing_monthly_every_ts), intent(out) :: var

    var%interp_type = 'linear'
    var%data_type   = 'monthly-calendar'
    var%interp_freq = 'every-timestep'
    var%filename    = 'not-used-for-monthly'
    var%data_label  = 'not-used-for-monthly'

  end subroutine marbl_init_monthly_forcing_metadata

  !***********************************************************************

  subroutine marbl_init_non_autotroph_tracer_metadata(marbl_tracer_metadata, &
       non_living_biomass_ecosys_tracer_cnt)

    !-----------------------------------------------------------------------
    !  initialize non-autotroph tracer_d values and accumulate
    !  non_living_biomass_ecosys_tracer_cnt
    !-----------------------------------------------------------------------

    implicit none

    type (marbl_tracer_metadata_type) , intent(inout) :: marbl_tracer_metadata(:)             ! descriptors for each tracer
    integer (int_kind)                , intent(out)   :: non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers

    integer(int_kind) :: n

    non_living_biomass_ecosys_tracer_cnt = 0

    marbl_tracer_metadata(po4_ind)%short_name='PO4'
    marbl_tracer_metadata(po4_ind)%long_name='Dissolved Inorganic Phosphate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(no3_ind)%short_name='NO3'
    marbl_tracer_metadata(no3_ind)%long_name='Dissolved Inorganic Nitrate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(sio3_ind)%short_name='SiO3'
    marbl_tracer_metadata(sio3_ind)%long_name='Dissolved Inorganic Silicate'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(nh4_ind)%short_name='NH4'
    marbl_tracer_metadata(nh4_ind)%long_name='Dissolved Ammonia'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(fe_ind)%short_name='Fe'
    marbl_tracer_metadata(fe_ind)%long_name='Dissolved Inorganic Iron'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(o2_ind)%short_name='O2'
    marbl_tracer_metadata(o2_ind)%long_name='Dissolved Oxygen'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(dic_ind)%short_name='DIC'
    marbl_tracer_metadata(dic_ind)%long_name='Dissolved Inorganic Carbon'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(dic_alt_co2_ind)%short_name='DIC_ALT_CO2'
    marbl_tracer_metadata(dic_alt_co2_ind)%long_name='Dissolved Inorganic Carbon, Alternative CO2'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(alk_ind)%short_name='ALK'
    marbl_tracer_metadata(alk_ind)%long_name='Alkalinity'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(doc_ind)%short_name='DOC'
    marbl_tracer_metadata(doc_ind)%long_name='Dissolved Organic Carbon'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(don_ind)%short_name='DON'
    marbl_tracer_metadata(don_ind)%long_name='Dissolved Organic Nitrogen'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(dop_ind)%short_name='DOP'
    marbl_tracer_metadata(dop_ind)%long_name='Dissolved Organic Phosphorus'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(dopr_ind)%short_name='DOPr'
    marbl_tracer_metadata(dopr_ind)%long_name='Refractory DOP'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(donr_ind)%short_name='DONr'
    marbl_tracer_metadata(donr_ind)%long_name='Refractory DON'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    marbl_tracer_metadata(docr_ind)%short_name='DOCr'
    marbl_tracer_metadata(docr_ind)%long_name='Refractory DOC'
    non_living_biomass_ecosys_tracer_cnt = non_living_biomass_ecosys_tracer_cnt + 1

    do n = 1, non_living_biomass_ecosys_tracer_cnt
       if (n == alk_ind) then
          marbl_tracer_metadata(n)%units      = 'meq/m^3'
          marbl_tracer_metadata(n)%tend_units = 'meq/m^3/s'
          marbl_tracer_metadata(n)%flux_units = 'meq/m^3 cm/s'
       else
          marbl_tracer_metadata(n)%units      = 'mmol/m^3'
          marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
          marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
       endif
    end do

  end subroutine marbl_init_non_autotroph_tracer_metadata

  !***********************************************************************

  subroutine marbl_check_ecosys_tracer_count_consistency(non_living_biomass_ecosys_tracer_cnt, marbl_status_log)

    implicit none

    integer (int_kind), intent(in) :: &
         non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers
    type(marbl_log_type)           , intent(inout) :: marbl_status_log

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
       write(error_msg, "(4A)"), "ecosys_tracer_cnt = ", ecosys_tracer_cnt, &
                                 "but computed ecosys_tracer_cnt = ", n
       call marbl_status_log%log_error(error_msg, "ecosys_mod::marbl_check_ecosys_tracer_count_consistency()")
        return
    endif

  end subroutine marbl_check_ecosys_tracer_count_consistency

  !***********************************************************************

  subroutine marbl_initialize_zooplankton_tracer_metadata(marbl_tracer_metadata, &
       non_living_biomass_ecosys_tracer_cnt, n)

    !-----------------------------------------------------------------------
    !  initialize zooplankton tracer_d values and tracer indices
    !-----------------------------------------------------------------------

    type (marbl_tracer_metadata_type), dimension(:), intent(inout) :: marbl_tracer_metadata   ! descriptors for each tracer

    integer (int_kind), intent(in) :: &
         non_living_biomass_ecosys_tracer_cnt ! number of non-autotroph ecosystem tracers

    integer (int_kind), intent(inout) :: n

    integer (int_kind) :: &
         zoo_ind            ! zooplankton functional group index

    n = non_living_biomass_ecosys_tracer_cnt + 1

    do zoo_ind = 1, zooplankton_cnt
       marbl_tracer_metadata(n)%short_name = trim(zooplankton(zoo_ind)%sname) // 'C'
       marbl_tracer_metadata(n)%long_name  = trim(zooplankton(zoo_ind)%lname) // ' Carbon'
       marbl_tracer_metadata(n)%units      = 'mmol/m^3'
       marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
       marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
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

  subroutine marbl_initialize_autotroph_tracer_metadata(marbl_tracer_metadata, n)

    !-----------------------------------------------------------------------
    !  initialize autotroph tracer_d values and tracer indices
    !-----------------------------------------------------------------------

    type (marbl_tracer_metadata_type), dimension(:), intent(inout) :: marbl_tracer_metadata   ! descriptors for each tracer
    integer(int_kind), intent(inout) :: n

    integer (int_kind) :: auto_ind ! zooplankton functional group index

    do auto_ind = 1, autotroph_cnt
       marbl_tracer_metadata(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Chl'
       marbl_tracer_metadata(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Chlorophyll'
       marbl_tracer_metadata(n)%units      = 'mg/m^3'
       marbl_tracer_metadata(n)%tend_units = 'mg/m^3/s'
       marbl_tracer_metadata(n)%flux_units = 'mg/m^3 cm/s'
       autotrophs(auto_ind)%Chl_ind = n
       n = n + 1

       marbl_tracer_metadata(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'C'
       marbl_tracer_metadata(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Carbon'
       marbl_tracer_metadata(n)%units      = 'mmol/m^3'
       marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
       marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%C_ind = n
       n = n + 1

       marbl_tracer_metadata(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Fe'
       marbl_tracer_metadata(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Iron'
       marbl_tracer_metadata(n)%units      = 'mmol/m^3'
       marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
       marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
       autotrophs(auto_ind)%Fe_ind = n
       n = n + 1

       if (autotrophs(auto_ind)%kSiO3 > c0) then
          marbl_tracer_metadata(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'Si'
          marbl_tracer_metadata(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' Silicon'
          marbl_tracer_metadata(n)%units      = 'mmol/m^3'
          marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
          marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%Si_ind = n
          n = n + 1
       else
          autotrophs(auto_ind)%Si_ind = 0
       endif

       if (autotrophs(auto_ind)%imp_calcifier .or. &
            autotrophs(auto_ind)%exp_calcifier) then
          marbl_tracer_metadata(n)%short_name = trim(autotrophs(auto_ind)%sname) // 'CaCO3'
          marbl_tracer_metadata(n)%long_name  = trim(autotrophs(auto_ind)%lname) // ' CaCO3'
          marbl_tracer_metadata(n)%units      = 'mmol/m^3'
          marbl_tracer_metadata(n)%tend_units = 'mmol/m^3/s'
          marbl_tracer_metadata(n)%flux_units = 'mmol/m^3 cm/s'
          autotrophs(auto_ind)%CaCO3_ind = n
          n = n + 1
       else
          autotrophs(auto_ind)%CaCO3_ind = 0
       endif
    end do

    if (my_task == master_task) THEN
       write (stdout, *) '----- autotroph tracer indices -----'
       do auto_ind = 1, autotroph_cnt
          write (stdout, *) 'Chl_ind(', trim(autotrophs(auto_ind)%sname), ') = '   , autotrophs(auto_ind)%Chl_ind
          write (stdout, *) 'C_ind(', trim(autotrophs(auto_ind)%sname), ') = '     , autotrophs(auto_ind)%C_ind
          write (stdout, *) 'Fe_ind(', trim(autotrophs(auto_ind)%sname), ') = '    , autotrophs(auto_ind)%Fe_ind
          write (stdout, *) 'Si_ind(', trim(autotrophs(auto_ind)%sname), ') = '    , autotrophs(auto_ind)%Si_ind
          write (stdout, *) 'CaCO3_ind(', trim(autotrophs(auto_ind)%sname), ') = ' , autotrophs(auto_ind)%CaCO3_ind
       end do
       write (stdout, *) '------------------------------------'
    endif

  end subroutine marbl_initialize_autotroph_tracer_metadata

  !***********************************************************************

  subroutine marbl_setup_local_tracers(k, column_kmt, tracer_module, tracer_local)
    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------

    integer(int_kind), intent(in) :: k
    integer(int_kind), intent(in) :: column_kmt
    real (r8), dimension(ecosys_tracer_cnt), intent(in) :: tracer_module ! tracer values
    real (r8), dimension(ecosys_tracer_cnt), intent(out) :: tracer_local ! local copies of model tracer concentrations
    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: n ! tracer index
    !-----------------------------------------------------------------------

    ! FIXME(bja, 2015-06) only need to loop over non-living-biomass-ecosys-tracer-cnt. 
    ! Does it actually need to be a loop?
    do n = 1, ecosys_tracer_cnt
       tracer_local(n) = max(c0, tracer_module(n))
       if ( k > column_kmt) then
          tracer_local(n) = c0
       end if
    end do

  end subroutine marbl_setup_local_tracers

  !***********************************************************************

  subroutine marbl_setup_local_zooplankton(k, column_kmt, &
       tracer_module, zoo_cnt, zoo, zooplankton_local)

    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------

    ! FIXME(bja, 2015-07) shortening zooplankton to zoo to avoid
    ! a namespace collision with the global imported into the
    ! module. Need to fix after global is removed.

    integer (int_kind)                               , intent(in)  :: k
    integer(int_kind)                                , intent(in)  :: column_kmt
    real (r8)                                        , intent(in)  :: tracer_module(:) ! tracer values
    integer(int_kind)                                , intent(in)  :: zoo_cnt
    type(zooplankton_type), dimension(zoo_cnt)       , intent(in)  :: zoo
    type(zooplankton_local_type), dimension(zoo_cnt) , intent(out) :: zooplankton_local

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: zoo_ind, n ! tracer index
    !-----------------------------------------------------------------------

    do zoo_ind = 1, zoo_cnt
       n = zoo(zoo_ind)%C_ind
       zooplankton_local(zoo_ind)%C = max(c0, tracer_module(n))
       if (k > column_kmt) then
          zooplankton_local(zoo_ind)%C = c0
       end if
    end do

  end subroutine marbl_setup_local_zooplankton

  !***********************************************************************

  subroutine marbl_setup_local_autotrophs(k, column_kmt, &
       tracer_module, auto_cnt, auto_meta, autotroph_loc)

    !-----------------------------------------------------------------------
    !  create local copies of model tracers
    !  treat negative values as zero
    !  apply mask to local copies
    !-----------------------------------------------------------------------

    ! FIXME(bja, 2015-07) autotroph --> auto are horrible names, but
    ! can't use full name until it is removed from the global namespace!

    integer (int_kind)                              , intent(in)  :: k
    integer(int_kind)                               , intent(in)  :: column_kmt
    real (r8), dimension(:)                         , intent(in)  :: tracer_module ! tracer values
    integer(int_kind)                               , intent(in)  :: auto_cnt ! autotroph_cnt
    type(autotroph_type), dimension(auto_cnt)       , intent(in)  :: auto_meta ! autotrophs
    type(autotroph_local_type), dimension(auto_cnt) , intent(out) :: autotroph_loc

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: auto_ind, tracer_ind ! tracer index
    !-----------------------------------------------------------------------

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

       if (k > column_kmt) then
          autotroph_loc(auto_ind)%Chl = c0
          autotroph_loc(auto_ind)%C = c0
          autotroph_loc(auto_ind)%Fe = c0
          autotroph_loc(auto_ind)%Si = c0
          autotroph_loc(auto_ind)%CaCO3 = c0
       end if
    end do

  end subroutine marbl_setup_local_autotrophs

  !***********************************************************************

  subroutine marbl_autotroph_consistency_check(auto_cnt, column_kmt, auto_meta, &
       autotroph_local)

    !-----------------------------------------------------------------------
    !  If any phyto box are zero, set others to zeros.
    !-----------------------------------------------------------------------

    ! FIXME(bja, 2015-07) autotroph --> auto are horrible names, but
    ! can't use full name until it is removed from the global namespace!

    integer(int_kind), intent(in) :: auto_cnt   ! autotroph_cnt
    integer(int_kind), intent(in) :: column_kmt ! number of active model layers
    type(autotroph_type), dimension(auto_cnt), intent(in) :: auto_meta ! autotrophs
    type(autotroph_local_type), dimension(auto_cnt,column_kmt), intent(inout) :: autotroph_local

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: auto_ind, k
    logical (log_kind) :: zero_mask
    !-----------------------------------------------------------------------

    do k = 1, column_kmt
       do auto_ind = 1, auto_cnt
          zero_mask = (autotroph_local(auto_ind,k)%Chl == c0 .or. &
               autotroph_local(auto_ind,k)%C == c0 .or. &
               autotroph_local(auto_ind,k)%Fe == c0)
          if (auto_meta(auto_ind)%Si_ind > 0) then
             zero_mask = zero_mask .or. autotroph_local(auto_ind,k)%Si == c0
          end if
          if (zero_mask) then
             autotroph_local(auto_ind,k)%Chl = c0
             autotroph_local(auto_ind,k)%C = c0
             autotroph_local(auto_ind,k)%Fe = c0
          end if
          if (auto_meta(auto_ind)%Si_ind > 0) then
             if (zero_mask) then
                autotroph_local(auto_ind,k)%Si = c0
             end if
          end if
          if (auto_meta(auto_ind)%CaCO3_ind > 0) then
             if (zero_mask) then
                autotroph_local(auto_ind,k)%CaCO3 = c0
             end if
          end if
       end do
    end do

  end subroutine marbl_autotroph_consistency_check

  !***********************************************************************

  subroutine marbl_compute_autotroph_elemental_ratios(auto_cnt, auto_meta,    &
             autotroph_local, tracer_local, autotroph_secondary_species)

    use marbl_parms     , only : epsC
    use marbl_parms     , only : gQsi_0
    use marbl_parms     , only : gQsi_max
    use marbl_parms     , only : gQsi_min

    integer (int_kind)         , intent(in) :: auto_cnt
    type(autotroph_type)       , intent(in) :: auto_meta(auto_cnt)             ! autotrophs
    type(autotroph_local_type) , intent(in) :: autotroph_local(auto_cnt)
    real (r8)                  , intent(in) :: tracer_local(ecosys_tracer_cnt) ! local copies of model tracer concentrations
    type(autotroph_secondary_species_type), intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real :: cks      ! constant used in  quota modification
    real :: cksi     ! constant used in Si quota modification
    integer(int_kind) :: auto_ind
    !-----------------------------------------------------------------------

    associate(                                                 &
         Fe_loc     => tracer_local(fe_ind),                   &
         SiO3_loc   => tracer_local(sio3_ind),                 &
         auto_C     => autotroph_local(:)%C,                   &
         auto_Chl   => autotroph_local(:)%Chl,                 &
         auto_Fe    => autotroph_local(:)%Fe,                  &
         auto_Si    => autotroph_local(:)%Si,                  &
         auto_CaCO3 => autotroph_local(:)%CaCO3,               &
         thetaC     => autotroph_secondary_species(:)%thetaC , & ! local Chl/C ratio (mg Chl/mmol C)
         QCaCO3     => autotroph_secondary_species(:)%QCaCO3 , & ! CaCO3/C ratio (mmol CaCO3/mmol C)
         Qfe        => autotroph_secondary_species(:)%Qfe,     & ! init fe/C ratio (mmolFe/mmolC)
         gQfe       => autotroph_secondary_species(:)%gQfe,    & ! fe/C for growth
         Qsi        => autotroph_secondary_species(:)%Qsi,     & ! initial Si/C ratio (mmol Si/mmol C)
         gQsi       => autotroph_secondary_species(:)%gQsi     & ! diatom Si/C ratio for growth (new biomass)
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

    cks = 10._r8
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

  subroutine marbl_compute_PAR(domain, gcm_state, auto_cnt, autotroph_local, PAR)

    !-----------------------------------------------------------------------
    !  compute PAR related quantities
    !  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
    !  0.45   fraction of incoming SW -> PAR (non-dim)
    !-----------------------------------------------------------------------

    ! PAR is intent(inout) because it components, while entirely set here, are allocated elsewhere

    integer(int_kind)              , intent(in)    :: auto_cnt
    type(marbl_domain_type)        , intent(in)    :: domain  
    type(marbl_gcm_state_type)     , intent(in)    :: gcm_state
    type(autotroph_local_type)     , intent(in)    :: autotroph_local(auto_cnt, domain%km)
    type(marbl_PAR_type)           , intent(inout) :: PAR

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real (r8) :: WORK1(domain%kmt)
    integer(int_kind) :: k, subcol_ind
    !-----------------------------------------------------------------------

    associate(                                        &
         dkm              => domain%km,               &
         column_kmt       => domain%kmt,              &
         delta_z          => domain%delta_z,          &
         PAR_nsubcols     => domain%num_PAR_subcols   &
         )

    !-----------------------------------------------------------------------
    ! set depth independent quantities, sub-column fractions and PAR at surface
    ! ignore provided shortwave where col_frac == 0
    !-----------------------------------------------------------------------

    PAR%col_frac(:) = gcm_state%PAR_col_frac(:)

    where (PAR%col_frac(:) > c0)
       PAR%interface(0,:) = f_qsw_par * gcm_state%surf_shortwave(:)
    elsewhere
       PAR%interface(0,:) = c0
    endwhere

    !-----------------------------------------------------------------------
    ! avoid further computations, such as computing attenuation coefficient, if there is no light
    ! treat forcing as a single dark value, by setting col_frac(1) to 1
    !-----------------------------------------------------------------------

    if (all(PAR%interface(0,:) == c0)) then
       PAR%col_frac(:)    = c0
       PAR%col_frac(1)    = c1
       PAR%interface(:,:) = c0
       PAR%avg(:,:)       = c0
       PAR%KPARdz(:)      = c0
       return
    end if

    !-----------------------------------------------------------------------
    ! compute attenuation coefficient over column
    !-----------------------------------------------------------------------

    ! FIXME (bja, 2015-07) move calculation outside and just pass in this
    ! work array as autotroph_Chl instead of passing in all of autotroph_loc?
    WORK1(:) = max(sum(autotroph_local(:,1:column_kmt)%Chl, dim=1), 0.02_r8)

    do k = 1, column_kmt

       if (WORK1(k) < 0.13224_r8) then
          PAR%KPARdz(k) = 0.000919_r8*(WORK1(k)**0.3536_r8)
       else
          PAR%KPARdz(k) = 0.001131_r8*(WORK1(k)**0.4562_r8)
       end if

       PAR%KPARdz(k) = PAR%KPARdz(k) * delta_z(k)

    enddo

    PAR%KPARdz(column_kmt+1:dkm) = c0

    !-----------------------------------------------------------------------
    ! propagate PAR values through column, only on subcolumns with PAR>0
    ! note that if col_frac is 0, then so is PAR
    !-----------------------------------------------------------------------

    WORK1(:) = exp(-PAR%KPARdz(1:column_kmt))

    do subcol_ind = 1, PAR_nsubcols
       if (PAR%interface(0,subcol_ind) > c0) then

          ! this look will probably not vectorize
          do k = 1, column_kmt
             PAR%interface(k,subcol_ind) = PAR%interface(k-1,subcol_ind) * WORK1(k)
          enddo
          PAR%interface(column_kmt+1:dkm,subcol_ind) = c0

          do k = 1, column_kmt
             PAR%avg(k,subcol_ind) = PAR%interface(k-1,subcol_ind) * (c1 - WORK1(k)) / PAR%KPARdz(k)
          enddo
          PAR%avg(column_kmt+1:dkm,subcol_ind) = c0

       else

          PAR%interface(1:dkm,subcol_ind) = c0
          PAR%avg(1:dkm,subcol_ind) = c0

       endif
   end do

   end associate

   end subroutine marbl_compute_PAR

  !***********************************************************************

  subroutine marbl_compute_carbonate_chemistry(domain, &
       temperature, salinity, &
       tracer_local, carbonate, &
       ph_prev_3d, ph_prev_alt_co2_3d, &
       zsat_calcite, zsat_aragonite)

    use co2calc_column        , only : comp_co3terms         
    use co2calc_column        , only : comp_co3_sat_vals     
    use co2calc_column        , only : thermodynamic_coefficients_type

    type(marbl_domain_type)        , intent(in)    :: domain
    real (r8)                      , intent(in)    :: temperature(domain%km)                    ! old potential temperature (C)
    real (r8)                      , intent(in)    :: salinity(domain%km)                       ! current salinity (msu)
    real (r8)                      , intent(in)    :: tracer_local(ecosys_tracer_cnt,domain%km) ! local copies of model tracer concentrations
    type(carbonate_type)           , intent(out)   :: carbonate(domain%km)
    real(r8)                       , intent(inout) :: ph_prev_3d(domain%km)
    real(r8)                       , intent(inout) :: ph_prev_alt_co2_3d(domain%km)
    real(r8)                       , intent(inout) :: zsat_calcite(domain%km)                   ! Calcite Saturation Depth
    real(r8)                       , intent(inout) :: zsat_aragonite(domain%km)                 ! Aragonite Saturation Depth

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer :: k
    type(thermodynamic_coefficients_type), dimension(domain%km) :: co3_coeffs
    logical(log_kind) , dimension(domain%km) :: mask
    logical(log_kind) , dimension(domain%km) :: pressure_correct
    real(r8)          , dimension(domain%km) :: ph_lower_bound
    real(r8)          , dimension(domain%km) :: ph_upper_bound
    real(r8)          , dimension(domain%km) :: press_bar ! pressure at level (bars)
    real(r8)          , dimension(domain%km) :: DIC_loc
    real(r8)          , dimension(domain%km) :: DIC_ALT_CO2_loc
    real(r8)          , dimension(domain%km) :: ALK_loc
    real(r8)          , dimension(domain%km) :: PO4_loc
    real(r8)          , dimension(domain%km) :: SiO3_loc
    !-----------------------------------------------------------------------

    ! make local copies instead of using associate construct because of gnu fortran bug
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=68546

    DIC_loc(:)         = tracer_local(dic_ind,:)
    DIC_ALT_CO2_loc(:) = tracer_local(dic_alt_co2_ind,:)
    ALK_loc(:)         = tracer_local(alk_ind,:)
    PO4_loc(:)         = tracer_local(po4_ind,:)
    SiO3_loc(:)        = tracer_local(sio3_ind,:)

    associate(                                                &
         dkm               => domain%km,                      &
         column_kmt        => domain%kmt,                     &
         CO3               => carbonate(:)%CO3,               &
         HCO3              => carbonate(:)%HCO3,              &
         H2CO3             => carbonate(:)%H2CO3,             &
         pH                => carbonate(:)%pH,                &
         CO3_sat_calcite   => carbonate(:)%CO3_sat_calcite,   &
         CO3_sat_aragonite => carbonate(:)%CO3_sat_aragonite, &
         CO3_ALT_CO2       => carbonate(:)%CO3_ALT_CO2,       &
         HCO3_ALT_CO2      => carbonate(:)%HCO3_ALT_CO2,      &
         H2CO3_ALT_CO2     => carbonate(:)%H2CO3_ALT_CO2,     &
         pH_ALT_CO2        => carbonate(:)%pH_ALT_CO2         &
         )

    pressure_correct = .TRUE.
    pressure_correct(1) = .FALSE.
    do k=1,dkm

      mask(k) = (k <= column_kmt)
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
    use marbl_parms, only : c10

    real(r8), intent(in)  :: column_temperature
    real(r8), intent(out) :: Tfunc

    Tfunc = Q_10**(((column_temperature + T0_Kelvin) - (Tref + T0_Kelvin)) / c10)

  end subroutine marbl_compute_function_scaling

  !***********************************************************************

  subroutine marbl_compute_Pprime(k, domain, auto_cnt, auto_meta, &
       autotroph_local, column_temperature, autotroph_secondary_species)

    use marbl_parms           , only : thres_z1_auto
    use marbl_parms           , only : thres_z2_auto

    integer(int_kind)                      , intent(in)  :: k
    type(marbl_domain_type)                , intent(in)  :: domain
    integer(int_kind)                      , intent(in)  :: auto_cnt
    type(autotroph_type)                   , intent(in)  :: auto_meta(auto_cnt)
    type(autotroph_local_type)             , intent(in)  :: autotroph_local(auto_cnt)
    real(r8)                               , intent(in)  :: column_temperature
    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    real(r8) :: f_loss_thres
    real(r8) :: C_loss_thres
    !-----------------------------------------------------------------------

    associate(                                           &
         zt     => domain%zt(:),                         &
         Pprime => autotroph_secondary_species(:)%Pprime & ! output
         )

    !  calculate the loss threshold interpolation factor
    if (zt(k) > thres_z1_auto) then
       if (zt(k) < thres_z2_auto) then
          f_loss_thres = (thres_z2_auto - zt(k))/(thres_z2_auto - thres_z1_auto)
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

  subroutine marbl_compute_Zprime(k, domain, &
       zoo_cnt, zoo_meta, zooC, &
       Tfunc, zooplankton_secondary_species)

    use marbl_parms           , only : c1, c0
    use marbl_parms           , only : thres_z1_zoo
    use marbl_parms           , only : thres_z2_zoo

    integer(int_kind)                        , intent(in)    :: k
    type(marbl_domain_type)                  , intent(in)    :: domain
    integer(int_kind)                        , intent(in)    :: zoo_cnt
    type(zooplankton_type)                   , intent(in)    :: zoo_meta(zoo_cnt)
    real(r8)                                 , intent(in)    :: zooC(zoo_cnt)
    real(r8)                                 , intent(in)    :: Tfunc
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: zoo_ind
    real(r8) :: f_loss_thres
    real(r8) :: C_loss_thres
    !-----------------------------------------------------------------------

    associate(                                                  &
         zt       => domain%zt(:),                              & !(km)
         Zprime   => zooplankton_secondary_species(:)%Zprime,   & !(zoo_cnt)
         zoo_loss => zooplankton_secondary_species(:)%zoo_loss  & !(zoo_cnt) output
         )

    !  calculate the loss threshold interpolation factor
    if (zt(k) > thres_z1_zoo) then
       if (zt(k) < thres_z2_zoo) then
          f_loss_thres = (thres_z2_zoo - zt(k))/(thres_z2_zoo - thres_z1_zoo)
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
       tracer_local, autotroph_secondary_species)

    use marbl_parms     , only : c1

    integer(int_kind)                      , intent(in)  :: auto_cnt
    type(autotroph_type)                   , intent(in)  :: auto_meta(auto_cnt)
    real(r8)                               , intent(in)  :: tracer_local(ecosys_tracer_cnt)
    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  Get relative nutrient uptake rates for autotrophs,
    !  min. relative uptake rate modifies C fixation in the manner
    !  that the min. cell quota does in GD98.
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt

       associate(                                                             &
                 DOP_loc => tracer_local(dop_ind),                            &
                 NO3_loc => tracer_local(no3_ind),                            &
                 NH4_loc => tracer_local(nh4_ind),                            &
                 PO4_loc => tracer_local(po4_ind),                            &
                 Fe_loc   => tracer_local(fe_ind),                            &
                 SiO3_loc => tracer_local(sio3_ind),                          &
                 ! OUTPUTS
                 VNO3  => autotroph_secondary_species(auto_ind)%VNO3,         &
                 VNH4  => autotroph_secondary_species(auto_ind)%VNH4,         &
                 VNtot => autotroph_secondary_species(auto_ind)%VNtot,        &
                 VFe   => autotroph_secondary_species(auto_ind)%VFe,          &
                 f_nut => autotroph_secondary_species(auto_ind)%f_nut,        &
                 VDOP  => autotroph_secondary_species(auto_ind)%VDOP,         &
                 VPO4  => autotroph_secondary_species(auto_ind)%VPO4,         &
                 VPtot => autotroph_secondary_species(auto_ind)%VPtot,        &
                 VSiO3 => autotroph_secondary_species(auto_ind)%VSiO3,        &
                 ! AUTO_META
                 kNO3   => auto_meta(auto_ind)%kNO3,                          &
                 kNH4   => auto_meta(auto_ind)%kNH4,                          &
                 kFe    => auto_meta(auto_ind)%kFe,                           &
                 kPO4   => auto_meta(auto_ind)%kPO4,                          &
                 kDOP   => auto_meta(auto_ind)%kDOP,                          &
                 kSiO3  => auto_meta(auto_ind)%kSiO3,                         &
                 Nfixer => auto_meta(auto_ind)%Nfixer                         &
                )

       VNO3 = (NO3_loc / kNO3) / (c1 + (NO3_loc / kNO3) + (NH4_loc / kNH4))
       VNH4 = (NH4_loc / kNH4) / (c1 + (NO3_loc / kNO3) + (NH4_loc / kNH4))
       VNtot = VNO3 + VNH4
       if (Nfixer) then
          VNtot = c1
       end if

       VFe = Fe_loc / (Fe_loc + kFe)

       VPO4 = (PO4_loc / kPO4) / (c1 + (PO4_loc / kPO4) + (DOP_loc / kDOP))
       VDOP = (DOP_loc / kDOP) / (c1 + (PO4_loc / kPO4) + (DOP_loc / kDOP))
       VPtot = VPO4 + VDOP

       if (kSiO3 > c0) then
          VSiO3 = SiO3_loc / (SiO3_loc + kSiO3)
       endif

       f_nut = min(VNtot, VFe)
       f_nut = min(f_nut, VPO4)
       if (kSiO3 > c0) then
          f_nut = min(f_nut, VSiO3)
       endif

    end associate

    end do

  end subroutine marbl_compute_autotroph_uptake

  !***********************************************************************

  subroutine marbl_compute_autotroph_photosynthesis (auto_cnt, PAR_nsubcols, &
       auto_meta, autotroph_loc, temperature, Tfunc, PAR_col_frac, PAR_avg,  &
       autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !     get photosynth. rate, phyto C biomass change, photoadapt
    !-----------------------------------------------------------------------

    use marbl_parms     , only : c0, c1
    use marbl_parms     , only : epsTinv

    integer(int_kind)                      , intent(in)    :: auto_cnt
    integer(int_kind)                      , intent(in)    :: PAR_nsubcols
    type(autotroph_type)                   , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)             , intent(in)    :: autotroph_loc(auto_cnt)
    real(r8)                               , intent(in)    :: temperature
    real(r8)                               , intent(in)    :: Tfunc
    real(r8)                               , intent(in)    :: PAR_col_frac(PAR_nsubcols)
    real(r8)                               , intent(in)    :: PAR_avg(PAR_nsubcols)
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind, subcol_ind
    real(r8) :: PCmax            ! max value of PCphoto at temperature TEMP (1/sec)
    real(r8) :: light_lim_subcol ! light_lim for a sub-column
    real(r8) :: PCphoto_subcol   ! PCphoto for a sub-column
    real(r8) :: pChl_subcol      ! Chl synth. regulation term (mg Chl/mmol N)
    real(r8) :: photoacc_subcol  ! photoacc for a sub-column
    !-----------------------------------------------------------------------

    do auto_ind = 1, auto_cnt

    associate(                                                                &
         ! local Chl/C ratio (mg Chl / mmol C)
         thetaC    => autotroph_secondary_species(auto_ind)%thetaC,           &
         ! INPUTS
         f_nut     => autotroph_secondary_species(auto_ind)%f_nut,            &
         VNtot    => autotroph_secondary_species(auto_ind)%VNtot,             &
         ! OUTPUTS
         light_lim => autotroph_secondary_species(auto_ind)%light_lim,        &
         PCPhoto   => autotroph_secondary_species(auto_ind)%PCPhoto,          &
         photoC    => autotroph_secondary_species(auto_ind)%photoC,           &
         photoacc => autotroph_secondary_species(auto_ind)%photoacc,          &
         ! AUTO_META
         PCref   => auto_meta(auto_ind)%PCref,                                &
         alphaPI => auto_meta(auto_ind)%alphaPI                               &
         )

       PCmax = PCref * f_nut * Tfunc
       if (temperature < autotrophs(auto_ind)%temp_thres) then
          PCmax = c0
       end if

       if (thetaC > c0) then
          light_lim = c0
          PCphoto   = c0
          photoacc  = c0

          do subcol_ind = 1, PAR_nsubcols
             if (PAR_avg(subcol_ind) > c0) then
                light_lim_subcol = (c1 - exp((-c1 * alphaPI * thetaC * PAR_avg(subcol_ind)) / (PCmax + epsTinv)))

                PCphoto_subcol = PCmax * light_lim_subcol

                ! GD 98 Chl. synth. term
                pChl_subcol = autotrophs(auto_ind)%thetaN_max * PCphoto_subcol / &
                   (autotrophs(auto_ind)%alphaPI * thetaC * PAR_avg(subcol_ind))
                photoacc_subcol = (pChl_subcol * PCphoto_subcol * Q / thetaC) * autotroph_loc(auto_ind)%Chl

                light_lim = light_lim + PAR_col_frac(subcol_ind) * light_lim_subcol
                PCphoto   = PCphoto   + PAR_col_frac(subcol_ind) * PCphoto_subcol
                photoacc  = photoacc  + PAR_col_frac(subcol_ind) * photoacc_subcol
             end if
          end do

          photoC = PCphoto * autotroph_loc(auto_ind)%C
       else
          light_lim = c0
          PCphoto   = c0
          photoacc  = c0
          photoC    = c0
       endif

    end associate

    end do

  end subroutine marbl_compute_autotroph_photosynthesis

  !***********************************************************************

  subroutine marbl_compute_autotroph_phyto_diatoms (auto_cnt, auto_meta, &
       autotroph_loc, autotroph_secondary_species)

    !-----------------------------------------------------------------------
    !  Get nutrient uptakes by small phyto based on calculated C fixation
    !-----------------------------------------------------------------------

    use marbl_parms     , only : c0
    use marbl_parms     , only : Q

    integer(int_kind)                      , intent(in)    :: auto_cnt
    type(autotroph_type)                   , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)             , intent(in)    :: autotroph_loc(auto_cnt)
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    !-----------------------------------------------------------------------

    associate(                                               &
         gQfe     => autotroph_secondary_species(:)%gQfe,    & ! fe/C for growth
         gQsi     => autotroph_secondary_species(:)%gQsi,    & ! diatom Si/C ratio for growth (new biomass)
         VNO3     => autotroph_secondary_species(:)%VNO3,    & ! input
         VNH4     => autotroph_secondary_species(:)%VNH4,    & ! input
         VNtot    => autotroph_secondary_species(:)%VNtot,   & ! input
         VPO4     => autotroph_secondary_species(:)%VPO4,    & ! input
         VDOP     => autotroph_secondary_species(:)%VDOP,    & ! input
         VPtot    => autotroph_secondary_species(:)%VPtot,   & ! input
         photoC   => autotroph_secondary_species(:)%photoC,  & ! input
         NO3_V    => autotroph_secondary_species(:)%NO3_V,   & ! output
         NH4_V    => autotroph_secondary_species(:)%NH4_V,   & ! output
         PO4_V    => autotroph_secondary_species(:)%PO4_V,   & ! output
         DOP_V    => autotroph_secondary_species(:)%DOP_V,   & ! output
         photoFe  => autotroph_secondary_species(:)%photoFe, & ! output
         photoSi  => autotroph_secondary_species(:)%photoSi  & ! output
         )

    do auto_ind = 1, auto_cnt

       if (VNtot(auto_ind) > c0) then
          NO3_V(auto_ind) = (VNO3(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
          NH4_V(auto_ind) = (VNH4(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
       else
          NO3_V(auto_ind) = c0
          NH4_V(auto_ind) = c0
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

    use marbl_parms     , only : parm_f_prod_sp_CaCO3
    use marbl_parms     , only : CaCO3_sp_thres
    use marbl_parms     , only : CaCO3_temp_thres1
    use marbl_parms     , only : CaCO3_temp_thres2
    use marbl_parms     , only : f_photosp_CaCO3

    integer(int_kind)                      , intent(in)    :: auto_cnt
    type(autotroph_type)                   , intent(in)    :: auto_meta(auto_cnt)
    type(autotroph_local_type)             , intent(in)    :: autotroph_loc(auto_cnt)
    real(r8)                               , intent(in)    :: temperature
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    !-----------------------------------------------------------------------

    associate(                                                       &
         f_nut      => autotroph_secondary_species(:)%f_nut,     & ! input
         photoC     => autotroph_secondary_species(:)%photoC,    & ! input
         CaCO3_PROD => autotroph_secondary_species(:)%CaCO3_PROD & ! output
         )

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%imp_calcifier) then
          CaCO3_PROD(auto_ind) = parm_f_prod_sp_CaCO3 * photoC(auto_ind)
          CaCO3_PROD(auto_ind) = CaCO3_PROD(auto_ind) * f_nut(auto_ind) * f_nut(auto_ind)

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

    use marbl_parms     , only : Q
    use marbl_parms     , only : r_Nfix_photo

    integer(int_kind)                          , intent(in)  :: auto_cnt
    type(autotroph_type)                       , intent(in)  :: auto_meta(auto_cnt)
    type(autotroph_secondary_species_type) , intent(out) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    real(r8) :: work1
    !-----------------------------------------------------------------------

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

    use marbl_parms, only : dps

    integer(int_kind)                          , intent(in)    :: auto_cnt
    type(autotroph_type)                       , intent(in)    :: auto_meta(auto_cnt)
    real(r8)                                   , intent(in)    :: Tfunc
    type(autotroph_secondary_species_type) , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind
    !-----------------------------------------------------------------------

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
            auto_meta(auto_ind)%mort2 * Pprime(auto_ind)**1.75_r8)
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
    use marbl_parms     , only : grz_fnc_michaelis_menten
    use marbl_parms     , only : grz_fnc_sigmoidal
    use marbl_parms     , only : c0

    integer(int_kind)                        , intent(in)    :: auto_cnt
    integer(int_kind)                        , intent(in)    :: zoo_cnt
    integer(int_kind)                        , intent(in)    :: grazer_prey_cnt
    type(autotroph_type)                     , intent(in)    :: auto_meta(auto_cnt)
    real(r8)                                 , intent(in)    :: Tfunc
    type(zooplankton_local_type)             , intent(in)    :: zooplankton_loc(zoo_cnt)
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind, auto_ind2
    integer  :: zoo_ind, zoo_ind2
    integer  :: pred_ind
    integer  :: prey_ind
    real(r8) :: work1, work2, work3, work4
    real(r8) :: graze_rate
    !-----------------------------------------------------------------------

    associate(                                                              &
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
                     min(spc_poc_fac * (Pprime(auto_ind)+0.5_r8)**1.5_r8,    &
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

    use marbl_parms     , only : c1
    use marbl_parms     , only : Qp_zoo_pom
    use marbl_parms     , only : parm_labile_ratio

    integer(int_kind)                        , intent(in)    :: auto_cnt
    integer(int_kind)                        , intent(in)    :: zoo_cnt
    type(autotroph_type)                     , intent(in)    :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type) , intent(inout) :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(inout) :: autotroph_secondary_species(auto_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind, zoo_ind
    real(r8) :: remaining_P      ! used in routing P from autotrophs w/ Qp different from Qp_zoo_pom
    !-----------------------------------------------------------------------

    associate(                                                               &
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

  subroutine marbl_compute_dissolved_organic_matter (k, auto_cnt, zoo_cnt,    &
             PAR_nsubcols, auto_meta, zooplankton_secondary_species,          &
             autotroph_secondary_species, PAR_col_frac, PAR_in, PAR_avg,      &
             dz1, tracer_local, dissolved_organic_matter)

    use marbl_parms     , only : Qfe_zoo
    use marbl_parms     , only : Qp_zoo_pom
    use marbl_parms     , only : Q

    use marbl_parms     , only : DOC_reminR_light
    use marbl_parms     , only : DON_reminR_light
    use marbl_parms     , only : DOP_reminR_light
    use marbl_parms     , only : DOC_reminR_dark
    use marbl_parms     , only : DON_reminR_dark
    use marbl_parms     , only : DOP_reminR_dark

    use marbl_parms     , only : DOCr_reminR0
    use marbl_parms     , only : DONr_reminR0
    use marbl_parms     , only : DOPr_reminR0
    use marbl_parms     , only : DOMr_reminR_photo

    integer(int_kind)                       , intent(in)  :: k
    integer                                 , intent(in)  :: auto_cnt
    integer                                 , intent(in)  :: zoo_cnt
    integer(int_kind)                       , intent(in)  :: PAR_nsubcols
    type(autotroph_type)                    , intent(in)  :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type), intent(in)  :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)  , intent(in)  :: autotroph_secondary_species(auto_cnt)
    real(r8)                                , intent(in)  :: PAR_col_frac(PAR_nsubcols)
    real(r8)                                , intent(in)  :: PAR_in(PAR_nsubcols)
    real(r8)                                , intent(in)  :: PAR_avg(PAR_nsubcols)
    real(r8)                                , intent(in)  :: dz1
    real(r8)                                , intent(in)  :: tracer_local(ecosys_tracer_cnt)
    type(dissolved_organic_matter_type)     , intent(out) :: dissolved_organic_matter

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind, subcol_ind
    real(r8) :: work
    real(r8) :: DOC_reminR        ! remineralization rate (1/sec)
    real(r8) :: DON_reminR        ! remineralization rate (1/sec)
    real(r8) :: DOP_reminR        ! remineralization rate (1/sec)
    real(r8) :: DOCr_reminR       ! remineralization rate (1/sec)
    real(r8) :: DONr_reminR       ! remineralization rate (1/sec)
    real(r8) :: DOPr_reminR       ! remineralization rate (1/sec)
    !-----------------------------------------------------------------------

    associate(                                                               &
         DOC_loc         => tracer_local(doc_ind)                          , &
         DON_loc         => tracer_local(don_ind)                          , &
         DOP_loc         => tracer_local(dop_ind)                          , &
         DONr_loc        => tracer_local(donr_ind)                         , &
         DOPr_loc        => tracer_local(dopr_ind)                         , &
         DOCr_loc        => tracer_local(docr_ind)                         , &
         Qfe             => autotroph_secondary_species(:)%Qfe             , & ! input
         remaining_P_dop => autotroph_secondary_species(:)%remaining_P_dop , & ! input
         auto_loss_doc   => autotroph_secondary_species(:)%auto_loss_doc   , & ! input
         auto_graze_doc  => autotroph_secondary_species(:)%auto_graze_doc  , & ! input
         zoo_loss_doc    => zooplankton_secondary_species(:)%zoo_loss_doc  , & ! input
         zoo_graze_doc   => zooplankton_secondary_species(:)%zoo_graze_doc , & ! input
         DOC_prod        => dissolved_organic_matter%DOC_prod              , & ! output production of DOC (mmol C/m^3/sec)
         DOC_remin       => dissolved_organic_matter%DOC_remin             , & ! output remineralization of DOC (mmol C/m^3/sec)
         DOCr_remin      => dissolved_organic_matter%DOCr_remin            , & ! output remineralization of DOCr
         DON_prod        => dissolved_organic_matter%DON_prod              , & ! output production of DON
         DON_remin       => dissolved_organic_matter%DON_remin             , & ! output remineralization of DON
         DONr_remin      => dissolved_organic_matter%DONr_remin            , & ! output remineralization of DONr
         DOP_prod        => dissolved_organic_matter%DOP_prod              , & ! output production of DOP
         DOP_remin       => dissolved_organic_matter%DOP_remin             , & ! output remineralization of DOP
         DOPr_remin      => dissolved_organic_matter%DOPr_remin              & ! output remineralization of DOPr
         )

    !-----------------------------------------------------------------------
      !  compute terms for DOM
      !-----------------------------------------------------------------------

      DOC_prod = sum(zoo_loss_doc(:)) + sum(auto_loss_doc(:)) + sum(auto_graze_doc(:)) + sum(zoo_graze_doc(:))
      DON_prod = Q * DOC_prod * f_toDON
      DOP_prod = Qp_zoo_pom * ( sum(zoo_loss_doc(:)) + sum(zoo_graze_doc(:)) )
      do auto_ind = 1, auto_cnt
         if (auto_meta(auto_ind)%Qp == Qp_zoo_pom) then
            DOP_prod = DOP_prod + auto_meta(auto_ind)%Qp * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind))
         else
            DOP_prod = DOP_prod + remaining_P_dop(auto_ind)
         endif
      end do

      !-----------------------------------------------------------------------
      !  Different remin rates in light and dark for semi-labile pools
      !-----------------------------------------------------------------------

      DOC_reminR = c0
      DON_reminR = c0
      DOP_reminR = c0

      do subcol_ind = 1, PAR_nsubcols
         if (PAR_col_frac(subcol_ind) > c0) then
            if (PAR_avg(subcol_ind) > 1.0_r8) then
               DOC_reminR = DOC_reminR + PAR_col_frac(subcol_ind) * DOC_reminR_light
               DON_reminR = DON_reminR + PAR_col_frac(subcol_ind) * DON_reminR_light
               DOP_reminR = DOP_reminR + PAR_col_frac(subcol_ind) * DOP_reminR_light
            else
               DOC_reminR = DOC_reminR + PAR_col_frac(subcol_ind) * DOC_reminR_dark
               DON_reminR = DON_reminR + PAR_col_frac(subcol_ind) * DON_reminR_dark
               DOP_reminR = DOP_reminR + PAR_col_frac(subcol_ind) * DOP_reminR_dark
            endif
         endif
      end do

      !-----------------------------------------------------------------------
      !  Refractory remin increased in top layer from photodegradation due to UV
      !-----------------------------------------------------------------------

      DOCr_reminR = DOCr_reminR0
      DONr_reminR = DONr_reminR0
      DOPr_reminR = DOPr_reminR0

      if (k == 1) then
         do subcol_ind = 1, PAR_nsubcols
            if ((PAR_col_frac(subcol_ind) > c0) .and. (PAR_in(subcol_ind) > 1.0_r8)) then
               work = PAR_col_frac(subcol_ind) * (log(PAR_in(subcol_ind))*0.4373_r8) * (10.0e2/dz1)
               DOCr_reminR = DOCr_reminR + work * DOMr_reminR_photo
               DONr_reminR = DONr_reminR + work * DOMr_reminR_photo
               DOPr_reminR = DOPr_reminR + work * DOMr_reminR_photo
            endif
         end do
      endif

      DOC_remin  = DOC_loc  * DOC_reminR
      DON_remin  = DON_loc  * DON_reminR
      DOP_remin  = DOP_loc  * DOP_reminR
      DOCr_remin = DOCr_loc * DOCr_reminR
      DONr_remin = DONr_loc * DONr_reminR
      DOPr_remin = DOPr_loc * DOPr_reminR

    end associate

  end subroutine marbl_compute_dissolved_organic_matter

  !***********************************************************************

  subroutine marbl_compute_large_detritus(k, auto_cnt, zoo_cnt, auto_meta, &
       zooplankton_secondary_species, autotroph_secondary_species, Fe_loc, &
       POC, P_CaCO3, P_SiO2, dust, P_iron, &
       Fe_scavenge, Fe_scavenge_rate)

    use marbl_parms     , only : f_graze_CaCO3_remin
    use marbl_parms     , only : f_graze_si_remin
    use marbl_parms     , only : Qfe_zoo
    use marbl_parms     , only : parm_Fe_scavenge_rate0
    use marbl_parms     , only : dust_fescav_scale
    use marbl_parms     , only : Fe_scavenge_thres1
    use marbl_parms     , only : fe_max_scale2
    use marbl_parms     , only : yps

    ! Note (mvertens, 2016-02), all the column_sinking_partiles must be intent(inout)
    ! rather than intent(out), since if they were intent(out) they would be automatically 
    ! deallocated on entry in this routine (this is not required behavior - but is
    ! standard)

    integer                                  , intent(in)    :: k
    integer                                  , intent(in)    :: auto_cnt
    integer                                  , intent(in)    :: zoo_cnt
    type(autotroph_type)                     , intent(in)    :: auto_meta(auto_cnt)
    type(zooplankton_secondary_species_type) , intent(in)    :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(in)    :: autotroph_secondary_species(auto_cnt)
    real(r8)                                 , intent(in)    :: Fe_loc
    type(column_sinking_particle_type)       , intent(inout) :: POC
    type(column_sinking_particle_type)       , intent(inout) :: P_CaCO3
    type(column_sinking_particle_type)       , intent(inout) :: P_SiO2
    type(column_sinking_particle_type)       , intent(inout) :: dust
    type(column_sinking_particle_type)       , intent(inout) :: P_iron
    real(r8)                                 , intent(out)   :: Fe_scavenge
    real(r8)                                 , intent(out)   :: Fe_scavenge_rate

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real(r8) :: work
    integer :: auto_ind
    !-----------------------------------------------------------------------

    associate(                                                             &
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
    !  2) scale by sinking mass flux (POMx4 + Dust + bSi + CaCO3)
    !  3) increase scavenging at higher iron concentrations (>0.8nM), reduce low (< 0.3nM)
    !-----------------------------------------------------------------------

    Fe_scavenge_rate = parm_Fe_scavenge_rate0

    Fe_scavenge_rate = Fe_scavenge_rate * &
         ((POC%sflux_out(k)     + POC%hflux_out(k)    ) * 4.0_r8*12.01_r8 + &
          (P_CaCO3%sflux_out(k) + P_CaCO3%hflux_out(k)) * P_CaCO3%mass + &
          (P_SiO2%sflux_out(k)  + P_SiO2%hflux_out(k) ) * P_SiO2%mass + &
          (dust%sflux_out(k)    + dust%hflux_out(k)   ) * dust_fescav_scale)

    if (Fe_loc < 0.3e-3_r8) then
       Fe_scavenge_rate = Fe_scavenge_rate * (Fe_loc / 0.3e-3_r8)
    end if

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

  subroutine marbl_compute_nitrif(k, PAR_nsubcols, column_kmt,                &
             PAR_col_frac, PAR_in, PAR_out, KPARdz, NH4_loc, nitrif)

    !-----------------------------------------------------------------------
    !  nitrate & ammonium
    !  nitrification in low light
    !  use exponential decay of PAR across model level to compute taper factor
    !-----------------------------------------------------------------------

    use marbl_parms, only : parm_nitrif_par_lim
    use marbl_parms, only : parm_kappa_nitrif

    integer(int_kind) , intent(in)  :: k
    integer(int_kind) , intent(in)  :: PAR_nsubcols
    integer(int_kind) , intent(in)  :: column_kmt
    real(r8)          , intent(in)  :: PAR_col_frac(PAR_nsubcols)
    real(r8)          , intent(in)  :: PAR_in(PAR_nsubcols)
    real(r8)          , intent(in)  :: PAR_out(PAR_nsubcols)
    real(r8)          , intent(in)  :: kPARdz
    real(r8)          , intent(in)  :: NH4_loc
    real(r8)          , intent(out) :: nitrif

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer(int_kind) :: subcol_ind
    real(r8) :: nitrif_subcol
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! skip computations for non-active layers or NH4 is 0
    !-----------------------------------------------------------------------

    nitrif = c0

    if ((k > column_kmt) .or. (NH4_loc == c0)) return

    do subcol_ind = 1, PAR_nsubcols
       if (PAR_col_frac(subcol_ind) > c0) then
          if (PAR_out(subcol_ind) < parm_nitrif_par_lim) then
             nitrif_subcol = parm_kappa_nitrif * NH4_loc
             if (PAR_in(subcol_ind) > parm_nitrif_par_lim) then
                nitrif_subcol = nitrif_subcol * &
                   log(PAR_out(subcol_ind) / parm_nitrif_par_lim) / (-KPARdz)
             end if
             nitrif = nitrif + PAR_col_frac(subcol_ind) * nitrif_subcol
          end if
       end if
    end do

  end subroutine marbl_compute_nitrif

  !***********************************************************************

  subroutine marbl_compute_denitrif(O2_loc, NO3_loc, DOC_remin, DOCr_remin,   &
             POC_remin, other_remin, sed_denitrif, denitrif)

    !-----------------------------------------------------------------------
    !  Compute denitrification under low O2 conditions
    !-----------------------------------------------------------------------

    real(r8) , intent(in)    :: O2_loc
    real(r8) , intent(in)    :: NO3_loc
    real(r8) , intent(in)    :: DOC_remin
    real(r8) , intent(in)    :: DOCr_remin
    real(r8) , intent(in)    :: POC_remin
    real(r8) , intent(in)    :: OTHER_REMIN
    real(r8) , intent(inout) :: SED_DENITRIF
    real(r8) , intent(out)   :: denitrif

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    real(r8) :: work
    !-----------------------------------------------------------------------

    work = ((parm_o2_min + parm_o2_min_delta) - O2_loc) / parm_o2_min_delta
    work = min(max(work, c0), c1)
    denitrif = work * ((DOC_remin + DOCr_remin + POC_remin * (c1 - POCremin_refract) &
         - other_remin) / denitrif_C_N  - sed_denitrif)

    ! scale down denitrif if computed rate would consume all NO3 in 10 days
    if (NO3_loc < ((c10*spd)*(denitrif+sed_denitrif))) then
       work = NO3_loc / ((c10*spd)*(denitrif+sed_denitrif))
       denitrif = denitrif * work
       sed_denitrif = sed_denitrif * work
    endif

  end subroutine marbl_compute_denitrif

  !***********************************************************************

  subroutine marbl_compute_dtracer_local (auto_cnt, zoo_cnt, auto_meta, zoo_meta, &
       autotroph_secondary_species, &
       zooplankton_secondary_species, &
       dissolved_organic_matter, &
       nitrif, denitrif, sed_denitrif, &
       Fe_scavenge, Fe_scavenge_rate, &
       P_iron_remin, POC_remin, &
       P_SiO2_remin, P_CaCO3_remin, other_remin, PON_remin, POP_remin, &
       restore_local, &
       O2_loc, o2_production, o2_consumption, &
       dtracer)

    integer                                  , intent(in)  :: auto_cnt
    integer                                  , intent(in)  :: zoo_cnt
    type(autotroph_type)                     , intent(in)  :: auto_meta(auto_cnt)
    type(zooplankton_type)                   , intent(in)  :: zoo_meta(zoo_cnt)
    type(zooplankton_secondary_species_type) , intent(in)  :: zooplankton_secondary_species(zoo_cnt)
    type(autotroph_secondary_species_type)   , intent(in)  :: autotroph_secondary_species(auto_cnt)
    type(dissolved_organic_matter_type)      , intent(in)  :: dissolved_organic_matter
    real(r8)                                 , intent(in)  :: nitrif
    real(r8)                                 , intent(in)  :: denitrif
    real(r8)                                 , intent(in)  :: sed_denitrif
    real(r8)                                 , intent(in)  :: Fe_scavenge
    real(r8)                                 , intent(in)  :: Fe_scavenge_rate
    real(r8)                                 , intent(in)  :: P_iron_remin
    real(r8)                                 , intent(in)  :: POC_remin
    real(r8)                                 , intent(in)  :: P_SiO2_remin
    real(r8)                                 , intent(in)  :: P_CaCO3_remin
    real(r8)                                 , intent(in)  :: other_remin
    real(r8)                                 , intent(in)  :: PON_remin
    real(r8)                                 , intent(in)  :: POP_remin
    real(r8)                                 , intent(in)  :: restore_local(ecosys_tracer_cnt)
    real(r8)                                 , intent(in)  :: O2_loc
    real(r8)                                 , intent(out) :: o2_production
    real(r8)                                 , intent(out) :: o2_consumption
    real(r8)                                 , intent(out) :: dtracer(ecosys_tracer_cnt)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer  :: auto_ind, zoo_ind, n
    real(r8) :: auto_sum
    !-----------------------------------------------------------------------

    associate(                                                            &
         thetaC          => autotroph_secondary_species(:)%thetaC          , & ! local Chl/C ratio (mg Chl/mmol C)
         QCaCO3          => autotroph_secondary_species(:)%QCaCO3          , & ! CaCO3/C ratio (mmol CaCO3/mmol C)
         Qfe             => autotroph_secondary_species(:)%Qfe             , & ! init fe/C ratio (mmolFe/mmolC)
         Qsi             => autotroph_secondary_species(:)%Qsi             , & ! initial Si/C ratio (mmol Si/mmol C)
         NO3_V           => autotroph_secondary_species(:)%NO3_V           , & ! nitrate uptake (mmol NO3/m^3/sec)
         NH4_V           => autotroph_secondary_species(:)%NH4_V           , & ! ammonium uptake (mmol NH4/m^3/sec)
         PO4_V           => autotroph_secondary_species(:)%PO4_V           , & ! PO4 uptake (mmol PO4/m^3/sec)
         DOP_V           => autotroph_secondary_species(:)%DOP_V           , & ! DOP uptake (mmol DOP/m^3/sec)
         photoC          => autotroph_secondary_species(:)%photoC          , & ! C-fixation (mmol C/m^3/sec)
         photoFe         => autotroph_secondary_species(:)%photoFe         , & ! iron uptake
         photoSi         => autotroph_secondary_species(:)%photoSi         , & ! silicon uptake (mmol Si/m^3/sec)
         photoacc        => autotroph_secondary_species(:)%photoacc        , & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
         auto_loss       => autotroph_secondary_species(:)%auto_loss       , & ! autotroph non-grazing mort (mmol C/m^3/sec)
         auto_loss_dic   => autotroph_secondary_species(:)%auto_loss_dic   , & ! auto_loss routed to dic (mmol C/m^3/sec)
         auto_loss_doc   => autotroph_secondary_species(:)%auto_loss_doc   , & ! auto_loss routed to doc (mmol C/m^3/sec)
         auto_agg        => autotroph_secondary_species(:)%auto_agg        , & ! autotroph aggregation (mmol C/m^3/sec)
         auto_graze      => autotroph_secondary_species(:)%auto_graze      , & ! autotroph grazing rate (mmol C/m^3/sec)
         auto_graze_zoo  => autotroph_secondary_species(:)%auto_graze_zoo  , & ! auto_graze routed to zoo (mmol C/m^3/sec)
         auto_graze_dic  => autotroph_secondary_species(:)%auto_graze_dic  , & ! auto_graze routed to dic (mmol C/m^3/sec)
         auto_graze_doc  => autotroph_secondary_species(:)%auto_graze_doc  , & ! auto_graze routed to doc (mmol C/m^3/sec)
         CaCO3_PROD      => autotroph_secondary_species(:)%CaCO3_PROD      , & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
         Nfix            => autotroph_secondary_species(:)%Nfix            , & ! total Nitrogen fixation (mmol N/m^3/sec)
         Nexcrete        => autotroph_secondary_species(:)%Nexcrete        , & ! fixed N excretion
         remaining_P_dip => autotroph_secondary_species(:)%remaining_P_dip , & ! remaining_P from mort routed to remin

         f_zoo_detr      => zooplankton_secondary_species(:)%f_zoo_detr    , & ! frac of zoo losses into large detrital pool (non-dim)
         x_graze_zoo     => zooplankton_secondary_species(:)%x_graze_zoo   , & ! {auto, zoo}_graze routed to zoo (mmol C/m^3/sec)
         zoo_graze       => zooplankton_secondary_species(:)%zoo_graze     , & ! zooplankton losses due to grazing (mmol C/m^3/sec)
         zoo_graze_zoo   => zooplankton_secondary_species(:)%zoo_graze_zoo , & ! grazing of zooplankton routed to zoo (mmol C/m^3/sec)
         zoo_graze_dic   => zooplankton_secondary_species(:)%zoo_graze_dic , & ! grazing of zooplankton routed to dic (mmol C/m^3/sec)
         zoo_graze_doc   => zooplankton_secondary_species(:)%zoo_graze_doc , & ! grazing of zooplankton routed to doc (mmol C/m^3/sec)
         zoo_loss        => zooplankton_secondary_species(:)%zoo_loss      , & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
         zoo_loss_dic    => zooplankton_secondary_species(:)%zoo_loss_dic  , & ! zoo_loss routed to dic (mmol C/m^3/sec)
         zoo_loss_doc    => zooplankton_secondary_species(:)%zoo_loss_doc  , & ! zoo_loss routed to doc (mmol C/m^3/sec)

         DOC_prod        => dissolved_organic_matter%DOC_prod        , & ! production of DOC (mmol C/m^3/sec)
         DOC_remin       => dissolved_organic_matter%DOC_remin       , & ! remineralization of DOC (mmol C/m^3/sec)
         DOCr_remin      => dissolved_organic_matter%DOCr_remin      , & ! remineralization of DOCr
         DON_prod        => dissolved_organic_matter%DON_prod        , & ! production of DON
         DON_remin       => dissolved_organic_matter%DON_remin       , & ! remineralization of DON
         DONr_remin      => dissolved_organic_matter%DONr_remin      , & ! remineralization of DONr
         DOP_prod        => dissolved_organic_matter%DOP_prod        , & ! production of DOP
         DOP_remin       => dissolved_organic_matter%DOP_remin       , & ! remineralization of DOP
         DOPr_remin      => dissolved_organic_matter%DOPr_remin        & ! remineralization of DOPr
         )

    !-----------------------------------------------------------------------
    !  nitrate & ammonium
    !-----------------------------------------------------------------------

    dtracer(no3_ind) = restore_local(no3_ind) + nitrif - denitrif - sed_denitrif - sum(NO3_V(:))

    dtracer(nh4_ind) = -sum(NH4_V(:)) - nitrif + DON_remin + DONr_remin  &
         + Q * (sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) + sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) &
                + DOC_prod*(c1 - f_toDON)) &
         + PON_remin * (c1 - PONremin_refract)

    do auto_ind = 1, auto_cnt
       if (auto_meta(auto_ind)%Nfixer) then
          dtracer(nh4_ind) = dtracer(nh4_ind) + Nexcrete(auto_ind)
       end if
    end do

    !-----------------------------------------------------------------------
    !  dissolved iron
    !-----------------------------------------------------------------------

    dtracer(fe_ind) = P_iron_remin - sum(photofe(:)) - Fe_scavenge &
       + Qfe_zoo * ( sum(zoo_loss_dic(:)) + sum(zoo_loss_doc(:)) + sum(zoo_graze_dic(:)) + sum(zoo_graze_doc(:)) )

    do auto_ind = 1, autotroph_cnt
       dtracer(fe_ind) = dtracer(fe_ind) &
            + (Qfe(auto_ind) * (auto_loss_dic(auto_ind) + auto_graze_dic(auto_ind))) &
            + auto_graze_zoo(auto_ind) * (Qfe(auto_ind) - Qfe_zoo) &
            + (Qfe(auto_ind) * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind)))
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
         + (c1 - POPremin_refract) * POP_remin + Qp_zoo_pom * ( sum(zoo_loss_dic(:)) + sum(zoo_graze_dic(:)) )

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

    dtracer(doc_ind) = DOC_prod * (c1 - DOCprod_refract) - DOC_remin

    dtracer(docr_ind) = DOC_prod * DOCprod_refract - DOCr_remin + (POC_remin * POCremin_refract)

    dtracer(don_ind) = (DON_prod * (c1 - DONprod_refract)) - DON_remin

    dtracer(donr_ind) = (DON_prod * DONprod_refract) - DONr_remin + (PON_remin * PONremin_refract)

    dtracer(dop_ind) = (DOP_prod * (c1 - DOPprod_refract)) - DOP_remin - sum(DOP_V(:))

    dtracer(dopr_ind) = (DOP_prod * DOPprod_refract) - DOPr_remin + (POP_remin * POPremin_refract)

    !-----------------------------------------------------------------------
    !  dissolved inorganic Carbon
    !-----------------------------------------------------------------------

    dtracer(dic_ind) = &
         sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) - sum(photoC(:)) &
            + DOC_remin + POC_remin * (c1 - POCremin_refract) + sum(zoo_loss_dic(:)) &
            + sum(zoo_graze_dic(:)) + P_CaCO3_remin + DOCr_remin

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
    o2_consumption = o2_consumption * ( (POC_remin * (c1 - POCremin_refract) + DOC_remin &
         + DOCr_remin - (sed_denitrif * denitrif_C_N) - other_remin + sum(zoo_loss_dic(:)) &
         + sum(zoo_graze_dic(:)) + sum(auto_loss_dic(:)) + sum(auto_graze_dic(:)) ) &
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

    real(r8)                            , intent(in)    :: tracer_local(ecosys_tracer_cnt)
    type(carbonate_type)                , intent(in)    :: carbonate
    type(dissolved_organic_matter_type) , intent(in)    :: dissolved_organic_matter
    real(r8)                            , intent(in)    :: QA_dust_def
    type(marbl_interior_share_type)     , intent(inout) :: marbl_interior_share

    associate( &
         share => marbl_interior_share &
         )

    share%QA_dust_def    = QA_dust_def
    share%DIC_loc_fields = tracer_local(DIC_ind)
    share%DOC_loc_fields = tracer_local(DOC_ind)
    share%O2_loc_fields  = tracer_local(O2_ind)
    share%NO3_loc_fields = tracer_local(NO3_ind)


    share%CO3_fields   = carbonate%CO3
    share%HCO3_fields  = carbonate%HCO3
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

    integer(int_kind)                        , intent(in)    :: zoo_cnt
    type(zooplankton_local_type)             , intent(in)    :: zooplankton_local(zoo_cnt)
    type(zooplankton_secondary_species_type) , intent(in)    :: zooplankton_secondary_species(zoo_cnt)
    type(marbl_zooplankton_share_type)       , intent(inout) :: marbl_zooplankton_share(zoo_cnt)

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

    integer(int_kind)                      , intent(in)    :: auto_cnt
    type(autotroph_local_type)             , intent(in)    :: autotroph_local(auto_cnt)
    type(autotroph_secondary_species_type) , intent(in)    :: autotroph_secondary_species(auto_cnt)
    type(marbl_autotroph_share_type)       , intent(inout) :: marbl_autotroph_share(auto_cnt)

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

  !*****************************************************************************

  function marbl_compute_totalChl(tracer_in, nb, ne) result(compute_totalChl)

    ! use specified indices because that is what autotrophs(:)%Chl_ind uses

    integer       , intent(in) :: nb, ne
    real(kind=r8) , intent(in) :: tracer_in(nb:ne)

    real(kind=r8) :: compute_totalChl

    integer :: auto_ind, n

    compute_totalChl = c0
    do auto_ind = 1, autotroph_cnt
       n = autotrophs(auto_ind)%Chl_ind
       compute_totalChl = compute_totalChl + max(c0, tracer_in(n))
    end do

  end function marbl_compute_totalChl

end module ecosys_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
