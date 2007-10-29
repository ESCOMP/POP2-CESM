#! /bin/csh -f
 
#=============
  set verbose
#=============

# This script copies remote, resolution-dependent ocean data files
# files to the local directory.  This script is invoked from the 
# ocean buildnml_prestage.csh script in the $CASEROOT/Buildnml_Prestage
# directory.
#
#
#

    if ($#argv < 1) then
       echo ocn.prestage.csh : command argument missing
       exit 1
    else
       set ocndata = $argv[1]
    endif

    set L = ( )   # local  filenames
    set R = ( )   # remote filenames

if ($OCN_GRID == gx1v3) then
#                =====
    set L = ( $L horiz_grid )    
    set R = ( $R $ocndata/grid/horiz_grid_20010402.ieeer8 )

    set L = ( $L region_mask )   
    set R = ( $R $ocndata/grid/region_mask_20010709.ieeei4 )

    set L = ( $L topography )    
    set R = ( $R $ocndata/grid/topography_20010702.ieeei4 )

    set L = ( $L shf )           
    set R = ( $R $ocndata/forcing/shf_mm_all_85-88_20010308.ieeer8 )

    set L = ( $L sfwf )          
    set R = ( $R $ocndata/forcing/sfwf_mm_PHC2_salx_flxio_20011012.ieeer8 )

    set L = ( $L chl_data )      
    set R = ( $R $ocndata/forcing/chl_mm_SeaWiFs97-01_20030328.ieeer8 )

else if ($OCN_GRID == gx1v4) then
#                     =====
    set L = ( $L horiz_grid )    
    set R = ( $R $ocndata/grid/horiz_grid_20010402.ieeer8 )

    set L = ( $L region_mask )   
    set R = ( $R $ocndata/grid/region_mask_20060206.ieeei4 )

    set L = ( $L topography )    
    set R = ( $R $ocndata/grid/topography_20060831.ieeei4 )

    set L = ( $L shf )           
    set R = ( $R $ocndata/forcing/shf_mm_all_85-88_20010308.ieeer8 )

    set L = ( $L sfwf )          
    set R = ( $R $ocndata/forcing/sfwf_mm_PHC2_salx_flxio_20060713.ieeer8 )

    set L = ( $L tidal_energy )      
    set R = ( $R $ocndata/forcing/tidal_mixing_energy_gx1v4_20060720.ieeer8 )
  
    set L = ( $L chl_data )      
    set R = ( $R $ocndata/forcing/chl_mm_SeaWiFs97-01_20030328.ieeer8 )
 
else if ($OCN_GRID == gx1v5) then
#                     =====
    set L = ( $L horiz_grid )    
    set R = ( $R $ocndata/grid/horiz_grid_20010402.ieeer8 )

    set L = ( $L region_mask )   
    set R = ( $R $ocndata/grid/region_mask_20061229.ieeei4 )

    set L = ( $L topography )    
    set R = ( $R $ocndata/grid/topography_20061229.ieeei4 )

    set L = ( $L shf )           
    set R = ( $R $ocndata/forcing/shf_mm_all_85-88_20010308.ieeer8 )

    set L = ( $L sfwf )          
    set R = ( $R $ocndata/forcing/sfwf_mm_PHC2_salx_flxio_20061230.ieeer8 )

    set L = ( $L tidal_energy )      
    set R = ( $R $ocndata/forcing/tidal_energy_gx1v5_20070102.ieeer8 )
  
    set L = ( $L chl_data )      
    set R = ( $R $ocndata/forcing/chl_filled_gx1v5_20061230.ieeer8 )
 
# the following files support a non-standard research option:

    set L = ( $L bathymetry  )    
    set R = ( $R $ocndata/grid/bathymetry_20070521.ieeer8 )

    set L = ( $L ts_PHC2_jan_ic_resindpt  )    
    set R = ( $R $ocndata/../res_indpt/ic/ts_PHC2_jan_ic_resindpt_20070920.nc )

else if ($OCN_GRID == gx3v5) then
#                     =====
    set L = ( $L horiz_grid )    
    set R = ( $R $ocndata/grid/horiz_grid_20030806.ieeer8 )

    set L = ( $L region_mask )   
    set R = ( $R $ocndata/grid/region_mask_20040220.ieeei4 )

    set L = ( $L topography )    
    set R = ( $R $ocndata/grid/topography_20040323.ieeei4 )

    set L = ( $L shf )           
    set R = ( $R $ocndata/forcing/shf_20031208.ieeer8 )

    set L = ( $L sfwf )          
    set R = ( $R $ocndata/forcing/sfwf_20040517.ieeer8 )

    set L = ( $L chl_data )      
    set R = ( $R $ocndata/forcing/chl_mm_SeaWiFs97-01_20031205.ieeer8 )
 
else if ($OCN_GRID == gx3v6) then
#                     =====
    set L = ( $L horiz_grid )    
    set R = ( $R $ocndata/grid/horiz_grid_20030806.ieeer8 )

    set L = ( $L region_mask )   
    set R = ( $R $ocndata/grid/region_mask_20041215.ieeei4 )

    set L = ( $L topography )    
    set R = ( $R $ocndata/grid/topography_20041215.ieeei4 )

    set L = ( $L shf )           
    set R = ( $R $ocndata/forcing/shf_20041215.ieeer8 )

    set L = ( $L sfwf )          
    set R = ( $R $ocndata/forcing/sfwf_20041215.ieeer8 )

    set L = ( $L chl_data )      
    set R = ( $R $ocndata/forcing/chl_mm_SeaWiFs97-01_20031205.ieeer8 )
  
    set L = ( $L tidal_energy )      
    set R = ( $R $ocndata/forcing/tidal_mixing_energy_gx3v6_20050406.ieeer8 )
  
    set L = ( $L buoyancy_freq )      
    set R = ( $R $ocndata/ic/buoyancy_freq_20050502.nc )
endif
  


 @ n = 1
 while ($n <= $#L)
    $UTILROOT/Tools/ccsm_getinput $R[$n]          $L[$n]           || exit 99
    $UTILROOT/Tools/ccsm_getinput $R[$n]:r.readme $L[$n]:r.readme  || exit 99
    echo ' '
 @ n++
 end # while

