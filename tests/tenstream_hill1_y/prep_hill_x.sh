#!/usr/bin/env bash

rm -f output.nc
if [ ! -f input.nc ]
then
    python hill1_y_input.py
fi
