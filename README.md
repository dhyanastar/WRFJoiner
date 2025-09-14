# WRFJoiner

A collection of Python scripts to join WRF output patches (`io_form_history = 102`) into a single NetCDF file.

## Description

When `io_form_history` or `io_form_restart` in the WRF's `namelist.input` is set to `102`, then the WRF produces output patches assigned to each MPI process. To join these patches into a single NetCDF file, you can use the Fortran program called as "JOINER" is provided at the [WRF Users' Page](https://www2.mmm.ucar.edu/wrf/users/special_code.html).

Since the JOINER program was developed years ago, it needs a slight modification for the WRFV4 outputs. As an alternative, I made some Python scripts that do the same job. A simple verification with my WRF outputs showed that these scripts work OK. Still, please use with caution if you are going to use these scripts.

1. `WRFJoiner_Parallel.py`

This script has a feature of fetching patch files in parallel using Python's multiprocessing module. Note that writing the output NetCDF file is in serial for now.

2. `WRFJoiner_Serial.py`

This script works using only a single process. One patch file is loaded and then pushed to the output NetCDF file, so this script uses less memory than the parallel version.

3. `namelist.ini`

This is a sample configuration file for this script. You need to define the following configurations. You can set most of them same as your `namelist.input`. Here are descriptions of some entrie.s
  * `loc_in`  : Location of the patch files
  * `loc_out` : Location of the output file to write
  * `interval`: Interval of WRF outputs in minute (= `history_interval` in minute)
  * `verbose` : Whether to print more lines
  * `nproc`   : Number of processes to use (only valid for the parallel version)

## Contact
Jeonghoe Kim (jeonghoekim.14@snu.ac.kr)