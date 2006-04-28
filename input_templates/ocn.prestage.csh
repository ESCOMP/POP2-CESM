#! /bin/csh -f
 
# This script builds a script to copy remote resolution-dependent ocean 
# files to local files.  This script is invoked from the CCSM pop2.template
# script, and the output from this script is used within 
# $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh 

cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2

echo ' -----------------------------------------'
echo ' prestage resolution-dependent ocean files'
echo ' -----------------------------------------'

    set L = ( )   # local  filenames
    set R = ( )   # remote filenames

EOF2
 
if ($OCN_GRID == gx1v3) then
#                =====
cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2
    set L = ( \$L horiz_grid )    
    set R = ( \$R \$ocndata/grid/horiz_grid_20010402.ieeer8 )

    set L = ( \$L region_mask )   
    set R = ( \$R \$ocndata/grid/region_mask_20010709.ieeei4 )

    set L = ( \$L topography )    
    set R = ( \$R \$ocndata/grid/topography_20010702.ieeei4 )

    set L = ( \$L shf )           
    set R = ( \$R \$ocndata/forcing/shf_mm_all_85-88_20010308.ieeer8 )

    set L = ( \$L sfwf )          
    set R = ( \$R \$ocndata/forcing/sfwf_mm_PHC2_salx_flxio_20011012.ieeer8 )

    set L = ( \$L chl_data )      
    set R = ( \$R \$ocndata/forcing/chl_mm_SeaWiFs97-01_20030328.ieeer8 )
EOF2

else if ($OCN_GRID == gx1v4) then
#                     =====
cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2
    set L = ( \$L horiz_grid )    
    set R = ( \$R \$ocndata/grid/horiz_grid_20010402.ieeer8 )

    set L = ( \$L region_mask )   
    set R = ( \$R \$ocndata/grid/region_mask_20060206.ieeei4 )

    set L = ( \$L topography )    
    set R = ( \$R \$ocndata/grid/topography_20060206.ieeei4 )

    set L = ( \$L shf )           
    set R = ( \$R \$ocndata/forcing/shf_mm_all_85-88_20010308.ieeer8 )

    set L = ( \$L sfwf )          
    set R = ( \$R \$ocndata/forcing/sfwf_mm_PHC2_salx_flxio_20060314.ieeer8 )

    set L = ( \$L chl_data )      
    set R = ( \$R \$ocndata/forcing/chl_mm_SeaWiFs97-01_20030328.ieeer8 )
EOF2
 
else if ($OCN_GRID == gx3v5) then
#                     =====
cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2
    set L = ( \$L horiz_grid )    
    set R = ( \$R \$ocndata/grid/horiz_grid_20030806.ieeer8 )

    set L = ( \$L region_mask )   
    set R = ( \$R \$ocndata/grid/region_mask_20040220.ieeei4 )

    set L = ( \$L topography )    
    set R = ( \$R \$ocndata/grid/topography_20040323.ieeei4 )

    set L = ( \$L shf )           
    set R = ( \$R \$ocndata/forcing/shf_20031208.ieeer8 )

    set L = ( \$L sfwf )          
    set R = ( \$R \$ocndata/forcing/sfwf_20040517.ieeer8 )

    set L = ( \$L chl_data )      
    set R = ( \$R \$ocndata/forcing/chl_mm_SeaWiFs97-01_20031205.ieeer8 )
EOF2
 
else if ($OCN_GRID == gx3v6) then
#                     =====
cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2
    set L = ( \$L horiz_grid )    
    set R = ( \$R \$ocndata/grid/horiz_grid_20030806.ieeer8 )

    set L = ( \$L region_mask )   
    set R = ( \$R \$ocndata/grid/region_mask_20041215.ieeei4 )

    set L = ( \$L topography )    
    set R = ( \$R \$ocndata/grid/topography_20041215.ieeei4 )

    set L = ( \$L shf )           
    set R = ( \$R \$ocndata/forcing/shf_20041215.ieeer8 )

    set L = ( \$L sfwf )          
    set R = ( \$R \$ocndata/forcing/sfwf_20041215.ieeer8 )

    set L = ( \$L chl_data )      
    set R = ( \$R \$ocndata/forcing/chl_mm_SeaWiFs97-01_20031205.ieeer8 )
  
    set L = ( \$L tidal_energy )      
    set R = ( \$R \$ocndata/forcing/tidal_mixing_energy_gx3v6_20050406.ieeer8 )
  
    set L = ( \$L buoyancy_freq )      
    set R = ( \$R \$ocndata/ic/buoyancy_freq_20050502.nc )
EOF2
endif
  
cat >> $CASEROOT/Buildnml_Prestage/pop2.buildnml_prestage.csh << EOF2

    @ n = 1
    while (\$n <= \$#L)
       \$UTILROOT/Tools/ccsm_getinput \$R[\$n]          \$L[\$n]           || exit 99
       \$UTILROOT/Tools/ccsm_getinput \$R[\$n]:r.readme \$L[\$n]:r.readme  || exit 99
    @ n++
    end # while

EOF2
