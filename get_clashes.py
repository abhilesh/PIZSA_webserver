import sys
import os 

cur_dir = os.getcwd()

outfile = open('clashes_train_chimera.csv', 'w')
outlist = []

for filename in os.listdir(cur_dir):
	if 'clashes' in filename and 'get' not in filename and 'chimera' not in filename:
		#print filename
		pdb_name = filename.replace('.clashes', '')
		#print pdb_name
		temp_f = open(filename, 'r')
		for line in temp_f.readlines():
			#print line
			if '3H' in line:
				print filename
			line = line.strip().split()
			#print line
			interface = line[0]
			clash = line[1]
			outlist.append([pdb_name, interface, clash])
			outfile.write(' '.join([pdb_name, interface, clash]) + '\n')

num_clash = {}

for element in outlist:
	if element[-1] not in num_clash.keys():
		num_clash.update({element[-1] : 1})
	else:
		num_clash[element[-1]] += 1


for element in num_clash:
	print element, num_clash[element]
