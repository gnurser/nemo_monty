#!/bin/bash

function run_build(){
    meson setup builddir --python.install-env auto
    meson compile -C builddir
    meson install -C builddir
    }
