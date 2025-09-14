# WRFJoiner.py
#
# This script is a Python implementation of the JOINER program provided at the
# WRF Users' Page.
# <https://www2.mmm.ucar.edu/wrf/users/special_code.html>
#
# Author: Jeonghoe Kim <jeonghoekim.14@snu.ac.kr>

import numpy   as np
import netCDF4 as nc
import configparser
import datetime
import glob
import multiprocessing
import itertools
import resource
import os
import sys

# from tqdm import tqdm

#-------------------------------------------------------------------------------
# Read the namelist.ini file.

def read_namelist_file(fil_namelist):
  namelist = configparser.ConfigParser()
  namelist.read(fil_namelist)

  loc_in     = namelist["NAMELIST"]["loc_in"    ]
  loc_out    = namelist["NAMELIST"]["loc_out"   ]
  start_date = namelist["NAMELIST"]["start_date"]
  end_date   = namelist["NAMELIST"]["end_date"  ]
  interval   = namelist["NAMELIST"]["interval"  ]
  dataset    = namelist["NAMELIST"]["dataset"   ]
  domain     = namelist["NAMELIST"]["domain"    ]
  e_we       = namelist["NAMELIST"]["e_we"      ]
  e_sn       = namelist["NAMELIST"]["e_sn"      ]
  nproc      = namelist["NAMELIST"]["nproc"     ]
  verbose    = namelist["NAMELIST"]["verbose"   ]

  loc_in     = loc_in.strip()
  loc_out    = loc_out.strip()
  start_date = datetime.datetime.strptime(start_date, "%Y-%m-%d_%H:%M:%S")
  end_date   = datetime.datetime.strptime(end_date  , "%Y-%m-%d_%H:%M:%S")
  interval   = int(interval)
  dataset    = [i.strip() for i in dataset.split(",")]
  domain     = [i.strip() for i in domain.split(",")]
  e_we       = [int(i) for i in e_we.split(",")]
  e_sn       = [int(i) for i in e_sn.split(",")]
  nproc      = int(nproc)
  verbose    = bool(verbose)

  return loc_in, loc_out, start_date, end_date, interval, dataset, domain, \
    e_we, e_sn, nproc, verbose

namelist = configparser.ConfigParser()
if len(sys.argv) > 1:
  fil_namelist = sys.argv[1]
else:
  fil_namelist = "namelist.ini"

loc_in, loc_out, start_date, end_date, interval, dataset, domain, \
  e_we, e_sn, nproc, verbose = read_namelist_file(fil_namelist)

assert (len(domain) == len(e_we)) & (len(domain) == len(e_sn))
dict_e_we = { domain[i] : e_we[i] for i in range(len(domain)) } 
dict_e_sn = { domain[i] : e_sn[i] for i in range(len(domain)) } 

#-------------------------------------------------------------------------------

def read_single_patch(fil):
  proc = int(fil.split("_")[-1])
  ncin = nc.Dataset(fil, "r")

  sn_start_stag   = ncin.getncattr("SOUTH-NORTH_PATCH_START_STAG"  )
  sn_end_stag     = ncin.getncattr("SOUTH-NORTH_PATCH_END_STAG"    )
  we_start_stag   = ncin.getncattr("WEST-EAST_PATCH_START_STAG"    )
  we_end_stag     = ncin.getncattr("WEST-EAST_PATCH_END_STAG"      )
  sn_start_unstag = ncin.getncattr("SOUTH-NORTH_PATCH_START_UNSTAG")
  sn_end_unstag   = ncin.getncattr("SOUTH-NORTH_PATCH_END_UNSTAG"  )
  we_start_unstag = ncin.getncattr("WEST-EAST_PATCH_START_UNSTAG"  )
  we_end_unstag   = ncin.getncattr("WEST-EAST_PATCH_END_UNSTAG"    )

  boundaries = {
    "sn_start_stag"   : sn_start_stag  ,
    "sn_end_stag"     : sn_end_stag    ,
    "we_start_stag"   : we_start_stag  ,
    "we_end_stag"     : we_end_stag    ,
    "sn_start_unstag" : sn_start_unstag,
    "sn_end_unstag"   : sn_end_unstag  ,
    "we_start_unstag" : we_start_unstag,
    "we_end_unstag"   : we_end_unstag  ,
  }

  variables = { name : var[:] for name, var in ncin.variables.items() }

  ncin.close()

  if verbose:
    print(f"Fetched PID {proc:5d}: {sn_start_stag:4d}, {sn_end_stag:4d}, {we_start_stag:4d}, {we_end_stag:4d}")
    # tqdm.write(f"Fetched PID {proc:5d}: {sn_start_stag:4d}, {sn_end_stag:4d}, {we_start_stag:4d}, {we_end_stag:4d}")

  return (proc, boundaries, variables)

#-------------------------------------------------------------------------------

def process(dataset_i, domain_i, time_i):
  filname = f"{dataset_i}_{domain_i}_{time_i.strftime('%Y-%m-%d_%H:%M:%S')}"
  fil_out = f"{loc_out}/{filname}"
  fils = sorted(glob.glob(f"{loc_in}/{filname}_*"))
  nfil = len(fils)
  if len(fils) == 0:
    return

  print(f"Processing {filname}")

  # Create a NetCDF file which will contain the joined results, by using the
  # first MPI process' NetCDF file structure as a reference.
  e_we_i           = dict_e_we[domain_i]
  e_sn_i           = dict_e_sn[domain_i]
  west_east        = e_we_i-1
  west_east_stag   = e_we_i
  south_north      = e_sn_i-1
  south_north_stag = e_sn_i

  ncin  = nc.Dataset(fils[0], "r")
  ncout = nc.Dataset(fil_out, "w", format="NETCDF4")

  ncout.setncatts(ncin.__dict__)
  setattr(ncout, "WEST-EAST_GRID_DIMENSION"      , west_east_stag  )
  setattr(ncout, "SOUTH-NORTH_GRID_DIMENSION"    , south_north_stag)
  setattr(ncout, "WEST-EAST_PATCH_START_STAG"    , 1               )
  setattr(ncout, "WEST-EAST_PATCH_END_STAG"      , west_east_stag  )
  setattr(ncout, "SOUTH-NORTH_PATCH_START_STAG"  , 1               )
  setattr(ncout, "SOUTH-NORTH_PATCH_END_STAG"    , south_north_stag)
  setattr(ncout, "WEST-EAST_PATCH_START_UNSTAG"  , 1               )
  setattr(ncout, "WEST-EAST_PATCH_END_UNSTAG"    , west_east       )
  setattr(ncout, "SOUTH-NORTH_PATCH_START_UNSTAG", 1               )
  setattr(ncout, "SOUTH-NORTH_PATCH_END_UNSTAG"  , south_north     )

  for name, dim in ncin.dimensions.items():
    if   name == "west_east":
      ncout.createDimension(name, west_east)
    elif name == "west_east_stag":
      ncout.createDimension(name, west_east_stag)
    elif name == "south_north":
      ncout.createDimension(name, south_north)
    elif name == "south_north_stag":
      ncout.createDimension(name, south_north_stag)
    else:
      ncout.createDimension(name, None if dim.isunlimited() else len(dim))

  for name, var in ncin.variables.items():
    ncout.createVariable(name, var.datatype, var.dimensions)
    ncout.variables[name].setncatts(var.__dict__)

  ncin.close()

  # Read the patches in parallel using multiprocessing module.
  with multiprocessing.Pool(nproc) as pool:
    fetched = pool.map(read_single_patch, fils)
    # fetched = list(tqdm(pool.imap(read_single_patch, fils), total=nfil, desc="Fetching WRF output patches"))

  # Write the fetched results into a single NetCDF file.
  for proc, boundaries, variables in fetched:
  # for proc, boundaries, variables in tqdm(fetched, total=len(fetched), desc="Writing into a single file"):
    sn_start_stag   = boundaries["sn_start_stag"  ]
    sn_end_stag     = boundaries["sn_end_stag"    ]
    we_start_stag   = boundaries["we_start_stag"  ]
    we_end_stag     = boundaries["we_end_stag"    ]
    sn_start_unstag = boundaries["sn_start_unstag"]
    sn_end_unstag   = boundaries["sn_end_unstag"  ]
    we_start_unstag = boundaries["we_start_unstag"]
    we_end_unstag   = boundaries["we_end_unstag"  ]

    perc = (proc+1) / nfil * 100
    if verbose:
      print(f"Writing PID {proc:5d} / {nfil-1:5d} ({perc:5.1f}%): {sn_start_stag:4d}, {sn_end_stag:4d}, {we_start_stag:4d}, {we_end_stag:4d}")
      # tqdm.write(f"Writing PID {proc:5d} / {nfil-1:5d} ({perc:5.1f}%): {sn_start_stag:4d}, {sn_end_stag:4d}, {we_start_stag:4d}, {we_end_stag:4d}")

    for name, data in variables.items():
      var  = ncout.variables[name]
      dims = var.dimensions

      # By default, index all elements (i.e., use colons) instead of slicing.
      sliced = [slice(None)] * len(dims)

      if "west_east"        in dims:
        idx = dims.index("west_east")
        sliced[idx] = slice(we_start_unstag-1, we_end_unstag)
      if "west_east_stag"   in dims:
        idx = dims.index("west_east_stag")
        sliced[idx] = slice(we_start_stag-1, we_end_stag)
      if "south_north"      in dims:
        idx = dims.index("south_north")
        sliced[idx] = slice(sn_start_unstag-1, sn_end_unstag)
      if "south_north_stag" in dims:
        idx = dims.index("south_north_stag")
        sliced[idx] = slice(sn_start_stag-1, sn_end_stag)

      var[tuple(sliced)] = data[:]

  ncout.close()

  peak_mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2)
  if verbose:
    print(f"Peak Memory Usage: {peak_mem_usage:.3f} GB")
  
  return

#-------------------------------------------------------------------------------

if __name__ == "__main__":
  times = []
  time_now = start_date
  while (time_now <= end_date):
    times.append(time_now)
    time_now += datetime.timedelta(minutes=interval)

  for dataset_i, domain_i, time_i in itertools.product(dataset, domain, times):
    process(dataset_i, domain_i, time_i)