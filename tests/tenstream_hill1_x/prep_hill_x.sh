#!/usr/bin/env bash

rm -f output.nc
if [ ! -f hill1_x_input.nc ]
then
    python hill1_x_input.py
fi
