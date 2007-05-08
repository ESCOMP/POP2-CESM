#! /bin/csh -f

#    ============================================
if ( ${OCN_GRID} == gx3v5 || ${OCN_GRID} == gx3v6) then
#    ============================================

 cat >> $POP2BLDSCRIPT << EOF3

#..... hmix_del2u_nml and hmix_del2t_nml
 set am_del2_value = 3.0e9
 set ah_del2_value = 1.0e7

#..... hmix_gm_nml
 set ah_gm_value    = 0.8e7
 set ah_bolus_value = 0.8e7
 set ah_bkg_srfbl   = 0.8e7

#..... hmix_aniso_nml
 set hmix_alignment_choice =  grid
 set lvariable_hmix_aniso  =  .true.
 set lsmag_aniso           =  .false.
 set visc_para =  1.0
 set visc_perp =  1.0
 set c_para    =  0.0
 set c_perp    =  0.0
 set u_para    =  0.0
 set u_perp    =  0.0
 set vconst_1  =  1.0e7
 set vconst_2  = 24.5
 set vconst_3  =  0.2
 set vconst_4  =  1.0e-8
 set vconst_5  =  3
 set vconst_6  = 1.0e7
 set vconst_7  = 90.0

#..... forcing_sfwf_nml
 set sfwf_weak_restore = 0.092

#..... transports_nml
 set transport_reg2_names = ("'Atlantic Ocean'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay'")

EOF3

#         ====================================================================
else if ( ${OCN_GRID} == gx1v3 || ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
#         ====================================================================

 cat >> $POP2BLDSCRIPT << EOF3

#..... hmix_del2u_nml and hmix_del2t_nml
 set am_del2_value = 0.5e8
 set ah_del2_value = 0.6e7

#..... hmix_gm_nml
 set ah_gm_value    = 0.6e7
 set ah_bolus_value = 0.6e7
 set ah_bkg_srfbl   = 0.6e7

#..... hmix_aniso_nml
 set hmix_alignment_choice =  east
 set lvariable_hmix_aniso =  .true.
 set lsmag_aniso          =  .false.
 set visc_para = 50.0e7
 set visc_perp = 50.0e7
 set c_para    =  8.0
 set c_perp    =  8.0
 set u_para    =  5.0
 set u_perp    =  5.0
 set vconst_1  =  0.6e7
 set vconst_2  =  0.5
 set vconst_3  =  0.16
 set vconst_4  =  2.e-8
 set vconst_5  =  3
 set vconst_6  =  0.6e7
 set vconst_7  =  45.0

#..... forcing_sfwf_nml
 set sfwf_weak_restore = 0.0115

#..... transports_nml
 set transport_reg2_names = ("'Atlantic Ocean'","'Mediterranean Sea'","'Labrador Sea'","'GIN Sea'","'Arctic Ocean'","'Hudson Bay'")

EOF3

endif



#    ====================
if ( ${OCN_GRID} == gx3v6 ) then
#    ====================

 cat >> $POP2BLDSCRIPT << EOF3

#.... tidal_nml
 set ltidal_mixing = .true.

#.... vmix_kpp_nml
 set bckgrnd_vdc1 = 0.1
 set bckgrnd_vdc2 = 0.0
EOF3

#         ============================================
else if ( ${OCN_GRID} == gx1v4 || ${OCN_GRID} == gx1v5 ) then
#         ============================================

 cat >> $POP2BLDSCRIPT << EOF3

#.... tidal_nml
 set ltidal_mixing = .true.

#.... vmix_kpp_nml
 set bckgrnd_vdc1 = 0.1
 set bckgrnd_vdc2 = 0.0
EOF3

else

 cat >> $POP2BLDSCRIPT << EOF3

#.... tidal_nml
 set ltidal_mixing = .false.

#.... vmix_kpp_nml
 set bckgrnd_vdc1 = 0.524
 set bckgrnd_vdc2 = 0.313
EOF3

endif
