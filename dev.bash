#!/bin/bash
echo $1
function isolate_build(){
    echo doing isolated build
    #make proper build environment with mamba
    # git clone theia:~agn/Bare/nemo_eos.git
    # cd monty
    mamba install compilers numpy netcdf4 numba
    export PREFIX_OUT=$CONDA_PREFIX
    mamba env create -f monty_dev.yml
    mamba activate monty_dev
    meson setup builddir --prefix $PREFIX_OUT # /Users/agn/miniforge3/envs/$mamba_env
    meson compile -C builddir
    meson install -C builddir
    mamba remove --name monty_dev --all
}

# install_rpath needs to be specified in the definition of py.extension_module in meson.build
#   ... the line below does not work.
# meson install -C builddir -Dinstall_rpath=$CONDA_PREFIX/lib
#note the discussion on "fixing' rpaths at
# https://github.com/mesonbuild/meson-python/pull/340
# https://github.com/mesonbuild/meson/issues/6541
#

function run_build(){
    #just add build requirements to mamba run environment
    echo doing run build
    # git clone theia:~agn/Bare/nemo_eos.git
    # cd monty
    mamba install compilers numpy netcdf4 numba
    mamba install meson-python pkg-config #=>

# + pyproject-metadata    0.7.1  pyhd8ed1ab_0  conda-forge     Cached
# + colorama              0.4.6  pyhd8ed1ab_0  conda-forge     Cached
# + ninja                1.11.1  hb8565cd_0    conda-forge     Cached
# + meson                 1.4.0  pyhd8ed1ab_0  conda-forge      634kB
# + meson-python         0.15.0  pyh0c530f3_0  conda-forge     Cached
# + pkg-config           0.29.2  ha3d46e9_1008 conda-forge     Cached
    meson setup builddir --python.install-env auto
    meson compile -C builddir
    meson install -C builddir
}
if [[ $1 == isolate_build ]]
then
    isolate_build
else
    run_build
fi
