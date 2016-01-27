#!/usr/bin/python
import sys

# changes Ne-n to fixed floating point

def find_num_start(string, st, end):
	index = end - 1

	while index >= st:
		if string[index] != '.' and  not ('0' <= string[index] <= '9'):		
			break
			
		index = index - 1
	
	return index + 1

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as fi:
	with open(output_file, 'w') as fo:
		for l in fi:
			new_l = l
			i = new_l.find("e-", 0)
			
			while i > -1:
				st = find_num_start(new_l, 0, i)
				
				digit = "0." + "".zfill(int(new_l[i+2:i+4]) - 1) + new_l[st] + new_l[st+2:i]
				
				print new_l[st:i+4] + ", " + digit
				
				new_l = new_l.replace(new_l[st:i+4], digit)
				
				i = new_l.find("e-", i + 1)
				
			fo.write(new_l)