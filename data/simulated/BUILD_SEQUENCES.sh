#!/bin/bash

Rscript ../../src/gen_guide_trees.R 1989 1989
python ../../src/fix_control.py control.unfixed.txt control.txt
~/bin/indelible
mv *.fas aligned/
