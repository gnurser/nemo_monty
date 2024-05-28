# nemo_monty

Provides executable python script `3dfiles_M_mp.py` that can produce a variety of diagnostics from a NEMO run. In particular, it can output Montgomery function, potential density temperature, salinity etc. on surfaces of constant Boussinesq density anomaly (see Aksenov et al., 2011).

## Pre-installation

Follow the installation instructions for the `nemo_eos` package.

## Installation

First go into parent directory into which you want to download the package.

```
    cd /path/to/parent
```

 Then, assuming that the mamba python environment into which `nemo_eos` has been installed is `big`, move into this environment (if it's something else, use that instead).

```
   mamba activate big
```

Then git clone the source from github and move into new directory

```
   git clone https://github.com/gnurser/nemo_monty.git
   cd nemo_monty
```

Build and install

```
   ./install.bash
```

Check that  `3dfiles_M_mp.py` can be found:
```
   which 3dfiles_M_mp.py
```
The directory enclosing this symbolic link is added to your `$PATH` while the `big` (or otherwise named) environment is activated by `mamba`.

## Usage
Running ` 3dfiles_M_mp.py --help` will give full details of all the available options:
```
(big) ~/VC_packages/nemo_monty (main)$ 3dfiles_M_mp.py --help
usage: 3dfiles_M_mp.py [-h] [--meshdir MESHDIR] [--meshfile MESHFILE] [--infile INFILE] [--hlimits HLIMITS HLIMITS HLIMITS HLIMITS] [-t [MTRACERS ...]]
                       [--passive_s [PASSIVE_S ...]] [-x [XTRACERS ...]] [--density DENSITY] [--TS0 TS0 TS0] [--neos NEOS] [--nthreads NTHREADS] [--dims DIMS [DIMS ...]]
                       [-r RUNNAME] [--xroot [XROOTS ...]] [--refresh] [-y YEARS [YEARS ...]] [--no_bounds] [--checkmask] [-y0 YEAR00] [--pass [PASSNO ...]] [-m [MONTHS ...]]
                       [-d [DAYS ...]] [-o OUTDIR] [--restarts [RTRACERS ...]]

Output various NEMO diagnostics

options:
  -h, --help            show this help message and exit
  --meshdir MESHDIR     name of mesh directory; can be set from environment variable MESHDIR
  --meshfile MESHFILE   name of meshfile inside mesh directory
  --infile INFILE       path of data file to read
  --hlimits HLIMITS HLIMITS HLIMITS HLIMITS
                        horizontal limits
  -t [MTRACERS ...], --tracers [MTRACERS ...]
                        names of mean tracers to read
  --passive_s [PASSIVE_S ...]
                        names of output passive tracers on surfaces
  -x [XTRACERS ...], --xtracers [XTRACERS ...]
                        names of calculated tracers to output
  --density DENSITY     layer density for layer output
  --TS0 TS0 TS0         initial guess for T0 and S0 on density layer
  --neos NEOS           choose EOS: -1=> old Jackett McDougall, 0=> poly EOS-80, 2=> TEOS-10
  --nthreads NTHREADS   number of threads for EOS; 16 is good for workstations
  --dims DIMS [DIMS ...]
                        dimensions in output cdf file
  -r RUNNAME            Run name
  --xroot [XROOTS ...]  Extra root directories
  --refresh             refresh cache?
  -y YEARS [YEARS ...]  last and first years
  --no_bounds           do not output meshbounds
  --checkmask           always use mask from mask.nc
  -y0 YEAR00            first year of dataset
  --pass [PASSNO ...]   pass number
  -m [MONTHS ...], --months [MONTHS ...]
                        months to save
  -d [DAYS ...], --days [DAYS ...]
                        days to save
  -o OUTDIR             directory of output file
  --restarts [RTRACERS ...]
                        name of restart tracers to read

```

Here we focus on diagnostics of Montgomery function, potential density temperature, salinity etc. on surfaces of constant Boussinesq density anomaly (see Aksenov et al., 2011). We diagnose these fields in the Dec 2021 output file from  Alex Megann's GO8p7 run.

```
(big) ~/VC_packages/nemo_monty (main) cd /noc/msm/scratch/nemo2/agn/ARCTIC
(big) /noc/msm/scratch/nemo2/agn/ARCTIC (main) export R12=/noc/msm/working/poles_apart/GO8p7/eORCA12
(big) /noc/msm/scratch/nemo2/agn/ARCTIC (main)$ 3dfiles_M_mp.py --meshdir $R12/../../stefryn --meshfile mesh_mask_CLASS-MEDUSA \
--infile $R12/monthly/nemo_g8p7ho_1m_202112-202112_grid-T.nc \
--hlimits 2500 3000 1000 2000  --tracers ssh T S \
--xtracers z_s sigma_s sigma_med_s T_s S_s mont  \
--density 26.5 --TS0 2. 34.5 --nthreads 1 --neos 2 --dims t y x -r GO8p7
```
This outputs a file `GO8p7__z_sigma_sigma_med_T_S_mont_26.5.nc` with properties on the `r_b=26.5` isopycnal surface: `z`, `sigma`, `sigma_med` (`sigma` referenced to the mean height of the surface), `T`, `S` and Montogomery function `mont`
Note that a meshfile is required including the `dep` variables as well as `mask` and `maskutil`; the `domain_cfg.nc` file is insufficent.



