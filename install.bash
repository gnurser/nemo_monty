#!/bin/bash

function run_build(){
    mamba install netcdf4 numba
    meson setup builddir --python.install-env auto
    meson compile -C builddir
    meson install -C builddir
    }
run_build
