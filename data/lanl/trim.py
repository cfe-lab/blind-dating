#!/bin/python

import re, sys

filename = sys.argv[1]
trimstart = int(sys.argv[2])
trimend = int(sys.argv[3])

file_ = open(filename, 'r')
text = file_.read()
file_.close()

for name, sequence in re.findall("(>.+)\n([agtcAGTC\n-]+)", text):
	print name
	print ''.join(sequence.split('\n'))[trimstart:trimend]
