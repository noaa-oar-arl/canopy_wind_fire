## Getting started

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/noaa-oar-arl/canopy-app/HEAD?labpath=python%2Fexamples.ipynb)

Try it out on Binder :point_up: or follow the steps below :point_down: to run locally.

---

Check that your active Python environment has the dependencies listed in [`environment.yml`](./environment.yml).
Otherwise, you can create a new environment using:

    conda env create -f environment.yml

And then activate it using:

    conda activate canopy-app

The primary function is `run()`, which runs canopy-app and returns an xarray Dataset.
By default, it runs the [default case](../input/namelist.canopy).
```python
from canopy_app import run

ds = run()
```

You can override default namelist options by passing a dictionary to `run()`.
```python
# Point setup
ds = run(
    config={
        "filenames": {"file_vars": "../input/point_file_20220701.sfcf000.txt"},
        "userdefs": {"infmt_opt": 1, "nlat": 1, "nlon": 1},
    },
)
```

There are also helper functions for running sets of experiments with different namelist options.
```python
# Multiple point setup cases
from canopy_app import config_cases, run_config_sens

cases = config_cases(
    file_vars="../input/point_file_20220701.sfcf000.txt",
    infmt_opt=1,
    ntime=1,
    nlat=1,
    nlon=1,
    z0ghc=[0.001, 0.01],
    lambdars=[1.0, 1.25],
    product=True,  # 4 cases if true, 2 if false
)
ds = run_config_sens(cases)
```
:point_up: We still get a single output dataset, but it has a `case` dimension.

### Generating global data

You can also generate and run global gridded canopy-app inputs using Python.  Note:  The python script here assumes that the user already has the global canopy data files (available by request or as staged on NOAA Hera or GMU Hopper) and GFS meteorological files (available by request or on NOAA Hera/HPSS) locally available (under the 'path' working directory set in the python script).  The user may generate their own gridded, meteorological NetCDF files conforming to the format provided here in the GFS examples and README. 

1. Edit python script (`global_data_process.py`)

2. Change global settings

   ```python
   '''User Options'''
   path     = '/scratch/pcampbe8/canopy-app/input'  # work directory
   ref_lev  = 10     # reference height above the canopy (m)
   frp_src  = 0      # frp data source (0: local source; 1: 12 month climatology; 2: all ones when ifcanwaf=.FALSE.)
   ```

2. Activate canopy-app Conda environment:

   ```
   conda activate canopy-app
   ```

3. Run Python script with time arguments:

   Time format: `YYYY(year)+MM(month)+DD(date)+HH(hour in UTC)+FFF(forecast hour)` Multiple time steps available. Use `,` to separate different times, no spaces.

   Example 1: 15-17 Jul 2020 12Z, forecast hour=0

   ```
   python global_data_process.py 2020071512000,2020071612000,2020071712000
   ```

   Example 2: 1 Jul 2020 12Z, forecast hour=0-4

   ```
   python global_data_process.py 2020070112000,2020070112001,2020070112002,2020070112003,2020070112004
   ```
### Running global data process and canopy-app with global data (following Example 1 above)

1. Edit namelist (`namelist.canopy`) for correct number of global lat/lon grid points

   For example:
   ```
   nlat        =  1536
   nlon        =  3072
   ```

2. Edit namelist (`namelist.canopy`) for correctly named global input files and output prefix

   For example:
   ```
   file_vars    = 'input/gfs.t12z.20220630.sfcf023.canopy.nc' 'input/gfs.t12z.20220701.sfcf000.canopy.nc' 'input/gfs.t12z.20220701.sfcf001.canopy.nc'
   file_out     = 'output/2022-07-01-11-0000_global'
   ```

3. Edit namelist (`namelist.canopy`) for correct time start/end and interval

   For example:
   ```
   time_start  = '2022-07-01-11:00:00.0000'
   time_end    = '2022-07-01-13:00:00.0000'
   ntime       =  3
   time_intvl  =  3600   !!! For hourly time steps
   ```
   `time_intvl = 86400` for daily time steps (i.e., 24*3600).

5. Run global process and canopy-app

   Running canopy-app globally requires a lot of memory, so suggest not running on head node. Slurm batch script suggestion (cpu time = ~20-30 min per model time step):
   ```
   #!/bin/bash -l
   #SBATCH --partition=bigmem             # big memory node
   #SBATCH --job-name=canopy-app          # name the job
   #SBATCH --output=canopy-%j.out         # write stdout to named file
   #SBATCH --error=canopy-%j.err          # write stdout to named file
   #SBATCH --time=0-02:00:00              # Run for max of 00 hrs, 10 mins, 00 secs
   #SBATCH --nodes=1                      # Request N nodes
   #SBATCH --ntasks=1                     # Request n tasks
   #SBATCH --mem-per-cpu=256GB            # Request nGB RAM per core
   conda activate canopy-app
   python python/global_data_process.py 2022063012023,2022070112000,2022070112001
   srun canopy
   ```
