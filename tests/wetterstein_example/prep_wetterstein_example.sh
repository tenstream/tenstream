#!/usr/bin/env bash

rm -f output.nc
if [ ! -f wetterstein_input.nc ]
then
    python wetterstein_input.py
fi

