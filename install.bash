#!/bin/bash

function run_build(){
if python -c "import numba" > /dev/null 2>&1; then
    if python -c "import netcdf4" > /dev/null 2>&1; then
	:
    else
	mamba install netcdf4
    fi
else
    if python -c "import netcdf4" > /dev/null 2>&1; then
	mamba install numba
    else
	mamba install netcdf4 numba
    fi
fi

    meson setup builddir --python.install-env auto
    meson compile -C builddir
    meson install -C builddir
    echo building
}

if python -c "import nemo_eos" > /dev/null 2>&1; then
    run_build
else
    echo nemo_eos not installed in this environment
fi
