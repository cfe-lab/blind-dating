#!/usr/bin/python
import sys

# changes Ne-n to fixed floating point

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as fi:
	with open(output_file, 'w') as fo:
		for l in fi:
			new_l = l
			i = new_l.find("e-", 0)
			
			while i > -1:
				st = new_l.rfind(":", 0, i) + 1
				
				digit = '0.' + "".zfill(int(new_l[i+2:i+4]) - 1) + new_l[st] + new_l[st+2:i]
				
				print new_l[st:i+4]
				print digit
				
				new_l = new_l.replace(new_l[st:i+4], digit)
				
				i = new_l.find("e-", i + 1)
				
			fo.write(new_l)