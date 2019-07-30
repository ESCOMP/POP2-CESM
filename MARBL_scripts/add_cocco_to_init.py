#!/usr/bin/env python

"""
Read in previous ecosystem initial conditions, copy all the sp* tracers to cocco*
There's also some xarray settings to keep all the metadata the same
(Although I do update the global history attribute to note the cocco additions)
"""

import os
import xarray as xr

# Store old initial condition files in a dictionary
files_to_update = dict()
files_to_update['gx3v7'] = os.path.join(os.path.sep, 'glade','p','cesmdata', 'cseg', 'inputdata', 'ocn', 'pop', 'gx3v7', 'ic', 'ecosys_jan_IC_gx3v7_20180308.nc')
files_to_update['gx1v7'] = os.path.join(os.path.sep, 'glade','p','cesmdata', 'cseg', 'inputdata', 'ocn', 'pop', 'gx1v6', 'ic', 'ecosys_jan_IC_gx1v6_20180308.nc')

# Dictionary of variables to copy (key = new var, value = old)
# Note that I originally started with coccoP = scalefactor * spC,
# which is why the new variable is the key instead of the value
var_copy_dict = dict()
var_copy_dict['coccoP'] = 'spP'
var_copy_dict['coccoC'] = 'spC'
var_copy_dict['coccoChl'] = 'spChl'
var_copy_dict['coccoFe'] = 'spFe'
var_copy_dict['coccoCaCO3'] = 'spCaCO3'

# For each key in files_to_update, update the file
for grid in files_to_update:
  # open_dataset should not do any decoding
  ds = xr.open_dataset(files_to_update[grid], decode_cf=False)
  # 1. update global history attribute
  hist = ds.attrs['history']
  ds.attrs['history'] = "{}\nJuly 12 2019: mlevy adds cocco fields (copies of sp fields)".format(hist)

  # 2. Copy sp tracers to cocco
  for new_var in var_copy_dict:
    print("Updating {} in {}...".format(new_var, grid))
    ds[new_var] = ds[var_copy_dict[new_var]]

  # 3. use encoding to avoid adding _FillValue to variables that didn't have that attribute (h/t to stackoverflow)
  # https://stackoverflow.com/questions/45693688/xarray-automatically-applying-fillvalue-to-coordinates-on-netcdf-output
  for var in ds:
    if '_FillValue' not in ds[var].attrs:
      ds[var].encoding['_FillValue'] = False

  # 4. Create new netcdf file
  print("Writing {}.nc...".format(grid))
  ds.to_netcdf('{}.nc'.format(grid))
  print('')
