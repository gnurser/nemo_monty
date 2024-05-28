#!/bin/bash

function run_build(){
    mamba install netcdf4 numba
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
