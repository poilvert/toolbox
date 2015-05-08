#!/bin/bash

SRC="resample_setup.py"
EXC=`which python`

rm -rf _resample.so
$EXC $SRC build_ext --inplace
mv -v resample/_resample.so .
rm -rvf _resample.c _resample.pyx.md5 resample/ build/
