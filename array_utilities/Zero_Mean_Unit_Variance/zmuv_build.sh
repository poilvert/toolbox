#!/bin/bash

SRC="zmuv_setup.py"
EXC=`which python`

rm -rf _zmuv.so
$EXC $SRC build_ext --inplace
mv -v zmuv/_zmuv.so .
rm -rvf _zmuv.c _zmuv.pyx.md5 zmuv/ build/
