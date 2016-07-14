#!/usr/bin/env bash

rm -f output.nc
if [ ! -f hill_gauss_x_input.nc ]
then
    python hill_gauss_x_input.py
fi
