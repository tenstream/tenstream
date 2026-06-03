#!/usr/bin/env bash

rm -f output.nc
if [ ! -f input.nc ]
then
    python3 wetterstein_input.py
fi

