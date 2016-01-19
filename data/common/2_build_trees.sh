#!/bin/sh

# Just a dummy wrapper around the 
# build tree script to keep everything organised

src_folder=${0%/2_build_trees.sh}
prefix=$1
use_raxml=$2

r -f $src_folder/build.trees.R --args $prefix $use_raxml
