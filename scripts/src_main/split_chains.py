# Place this script in the folder with the multimer PDBs

# write a dictionary of the interfaces to extract
# key = PDB ID, value = chain_1 + '-' + chain_2
chain_map = {'1cdl' : 'A-E', '1dva' : 'H-X', '1dx5' : 'B-N',
			'1dx5' : 'N-J', '1ebp' : 'A-C', '1es7' : 'A-B',
			'1es7' : 'A-C', '1es7' : 'A-D', '1fak' : 'H-T',
			'1fak' : 'L-T', '1fe8' : 'A-L', '1fe8' : 'A-H',
			'1foe' : 'A-B', '1g3i' : 'A-G', '1gl4' : 'A-B',
			'1ihb' : 'A-B', '1jat' : 'A-B', '1jpp' : 'B-D',
			'1mq8' : 'A-B', '1nfi' : 'A-F', '1nfi' : 'B-F',
			'1nun' : 'A-B', '1ub4' : 'A-C', '2hhb' : 'A-B'}

def read_pdb_lines(input_pdb):

	'''open the PDB file and parse the ATOM entries'''

	pdb_lines = []

	with open(input_pdb, 'r') as fin:
		for line in fin.readlines():
			if line[:4] == "ATOM":
				line = line.strip()
				pdb_lines.append(line)

	return pdb_lines

def split_multimer(pdb_lines, chain_1, chain_2):

	'''module to subset pdb file for only the subunits specified'''

	split_file = []

	for line in pdb_lines:
		chain = line[21]
		if chain == chain_1 or chain == chain_2:
			split_file.append(line)

	return split_file

def main(chain_map):

	'''The main function for splitting the multimer'''

	for pdb in chain_map.keys():

		input_pdb = pdb + '.pdb'
		interface = chain_map[pdb]
		chain_1 = interface.split('-')[0]
		chain_2 = interface.split('-')[1]

		pdb_lines = read_pdb_lines(input_pdb)

		split_file = split_multimer(pdb_lines, chain_1, chain_2)

		with open(pdb+'_'+chain_1+'-'+chain_2+'.pdb', 'w') as fout:
			for line in split_file:
				fout.write(line + '\n')

	return

if __name__ == '__main__':
	main(chain_map)
