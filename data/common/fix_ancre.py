#!/usr/bin/python

with open('out.txt') as f:
	cur_file_name = None

	for l in f:
		lstrip = l.strip()
		
		if '.' in lstrip:
			if cur_file_name != None:
				with open('../aligned_marked/' + cur_file_name, 'w') as cur_f:
					cur_f.write(cur_file)
		
			cur_file_name = lstrip
		
			with open(cur_file_name) as cur_f:
				cur_file = cur_f.read()
		else:
			if lstrip != 'REFERENCE':
				lsplit = lstrip.split(' ')
				name = lsplit[0]
				marker = '_PBMC' if len(lsplit) >= 2 and lsplit[1] == 'DNA' else '_PLASMA'
				index = cur_file.rindex('_', 0, cur_file.index('\n', cur_file.index(name)))
				cur_file = cur_file[:index] + marker + cur_file[index:]

if cur_file_name != None:
	with open('../aligned_marked/' + cur_file_name, 'w') as cur_f:
		cur_f.write(cur_file)