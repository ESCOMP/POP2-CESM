module ecosys_fields

! !MODULE: ecosys_fields

!-----------------------------------------------------------------------------
!   This module contains definitions of variables that are used in both 
!   ecocys_mod.F90 and ecosys_Ciso_mod.F90. They are shared using threading 
!   with pointers, and need to be pointed to in the code.
!-----------------------------------------------------------------------------
  USE blocks, ONLY: nx_block, ny_block
  USE domain_size, ONLY: max_blocks_clinic
  USE ecosys_parms

  implicit none
  
  public
  save

   
  
   type(sinking_particle), save :: & 
      POC,            & ! base units = nmol C
      P_CaCO3           ! base units = nmol CaCO3


  real (r8), dimension(nx_block,ny_block,max_blocks_clinic), target, public :: &
      DIC_SURF_fields,       & ! surface values of DIC for solver
      CO2STAR_SURF_fields,   & ! CO2STAR from solver
      DCO2STAR_SURF_fields,  &! DCO2STAR from solver
      PV_SURF_fields,        & ! piston velocity (cm/s)
      CO3_fields,            & ! carbonate ion
      CO3_SURF_fields,       & ! Surface carbonate ion
      HCO3_fields,           & ! bicarbonate ion
      H2CO3_fields,          & ! carbonic acid
      f_zoo_detr_fields,     & ! frac of zoo losses into large detrital pool (non-dim)
      DIC_loc_fields,        & ! local copy of model DIC
      DOC_loc_fields,        & ! local copy of model DOC
!      spCaCO3_loc_fields,    & ! local copy of model spCaCO3
      zooC_loc_fields,       & ! local copy of model zooC
!      DOM_remin_fields,      & ! fraction of DOM remineralized at current TEMP
      DECAY_CaCO3_fields,    & ! scaling factor for dissolution of CaCO3
      DECAY_Hard_fields,     & ! scaling factor for dissolution of Hard Ballast
      zoo_loss_fields,       & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
      zoo_loss_doc_fields,   & ! zoo_loss routed to doc (mmol C/m^3/sec)
      zoo_loss_dic_fields,   & ! zoo_loss routed to dic (mmol C/m^3/sec)
      POC_PROD_avail_fields, & ! POC production available for excess POC flux
      PCphoto_fields           ! C-specific rate of photosynth. (1/sec)
   
    real (r8), dimension(nx_block,ny_block,autotroph_cnt,max_blocks_clinic), target, public :: &
      CaCO3_PROD_fields,        & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      QCaCO3_fields,            & ! small phyto CaCO3/C ratio (mmol CaCO3/mmol C)
      autotrophCaCO3_loc_fields,& ! local copy of model autotroph CaCO3
      autotrophChl_loc_fields,  & ! local copy of model autotroph Chl
      autotrophC_loc_fields,    & ! local copy of model autotroph C
      autotrophFe_loc_fields,   & ! local copy of model autotroph Fe
      autotrophSi_loc_fields,   & ! local copy of model autotroph Si
      auto_graze_fields,        & ! autotroph grazing rate (mmol C/m^3/sec)
      auto_graze_zoo_fields,    & ! auto_graze routed to zoo (mmol C/m^3/sec)
      auto_graze_poc_fields,    & ! auto_graze routed to poc (mmol C/m^3/sec)
      auto_graze_doc_fields,    & ! auto_graze routed to doc (mmol C/m^3/sec)
      auto_graze_dic_fields,    & ! auto_graze routed to dic (mmol C/m^3/sec)
      auto_loss_fields,         & ! autotroph non-grazing mort (mmol C/m^3/sec)
      auto_loss_poc_fields,     & ! auto_loss routed to poc (mmol C/m^3/sec)
      auto_loss_doc_fields,     & ! auto_loss routed to doc (mmol C/m^3/sec)
      auto_loss_dic_fields,     & ! auto_loss routed to dic (mmol C/m^3/sec)
      auto_agg_fields,          & ! autotroph aggregation (mmol C/m^3/sec)
      photoC_fields              ! C-fixation (mmol C/m^3/sec)
 
 
 end module ecosys_fields
