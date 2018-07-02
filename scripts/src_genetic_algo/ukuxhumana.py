#!/usr/bin/env python

'''Genetic algorithm to optimize protein-protein interaction matrix'''

''' Start with an initial population of 150 potential matrices, mutate 10% of the residues '''


import time
import cPickle as pickle


def addoptions():

	'''Provides the user flags to tweak the program'''

	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('input_pdb', help = 'Input PDB file for the protein complex')

	parser.add_argument('-d', '--cutoff', 
						type = float,
						choices = [4.0, 6.0, 8.0],
						default = 4.0,
						help = 'Specify the distance cut-off')

	parser.add_argument('-t', '--intertype',
						choices = ["all", "mm", "ms", "ss"],
						default = "ss", 
						help = 'Specify the interacting atom types')

	parser.add_argument('-o', '--outfile',
						help = 'Specify the output file')

	parser.add_argument('-cp', '--custom_pot',
						help = 'Specify the custom potential file in .csv format')

	parser.add_argument('-p1', '--protein_1',
						help = 'Specify the first protein of the interface')

	parser.add_argument('-p2', '--protein_2',
						help = 'Specify the second protein of the interface')

	parser.add_argument('-alscn', '--alascan',
						nargs = '*',
						default = '0',
						help = 'Specify whether Alanine Scanning is required')

	args = parser.parse_args()

	pdb_input = args.input_pdb
	cut_off = args.cutoff
	inter_type = args.intertype
	outfile = args.outfile
	protein_1 = args.protein_1
	protein_2 = args.protein_2
	custom_pot = args.custom_pot
	alscn = args.alascan

	return pdb_input, cut_off, inter_type, protein_1, protein_2, outfile, custom_pot, alscn


def get_options():

	'''Parse the user provided arguments and return a dictionary'''

	options = {}

	arguments = addoptions()

	options['input_pdb'] = arguments[0]
	options['cutoff'] = arguments[1]
	options['intertype'] = arguments[2]
	options['protein1'] = arguments[3]
	options['protein2'] = arguments[4]
	options['outfile'] = arguments[5]
	options['custom_pot'] = arguments[6]
	options['alascan'] = arguments[7]

	return options


def select_pot(cut_off, inter_type):

	'''Selects the appropriate statistical potential'''

	if cut_off == 4:
		if inter_type == 'all':
			pot = '4_all_cmpd_cifa_avg.p'
		elif inter_type == 'mm':
			pot = '4_mm_norm_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '4_ms_norm_cifa_avg.p'
		elif inter_type == 'ss':
			pot = '4_ss_norm_cifa_avg.p'
	elif cut_off == 6:
		if inter_type == 'all':
			pot = '6_all_cmpd_cifa_avg.p'
		elif inter_type == 'mm':
			pot = '6_mm_cmpd_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '6_ms_cmpd_c_ij_no_avg.p'
		elif inter_type == 'ss':
			pot = '6_ss_cmpd_cifa_no_avg.p'
	elif cut_off == 8:
		if inter_type == 'all':
			pot = '8_all_norm_cifa_avg.p'
		elif inter_type == 'mm':
			pot = '8_mm_cmpd_cifa_avg.p'
		elif inter_type == 'ms':
			pot = '8_ms_cmpd_c_ij_no_avg.p'
		elif inter_type == 'ss':
			pot = '8_ss_cmpd_cifa_avg.p'

	return pot


def parse_pdb(pdb_file):

	'''Parses the PDB file to return the number of subunits and residues in each subunit'''

	pdb_lines = []
	chain_set = set()

	pdb = open(pdb_file, 'r').readlines()

	for line in pdb:
		if line[0:4] == 'ATOM':
			chain = line[21]
			chain_set.add(chain)
			pdb_lines.append(line)

	chain_set = sorted(chain_set)
	chain_num = len(chain_set)

	return chain_num, chain_set, pdb_lines


def parse_dist_all(pdb_file):

	'''Parse dist file in case the software is run for all the interfaces of a multimeric complex'''

	filtered_pdb = altloc_check(pdb_file)

	dist_list = get_dist(filtered_pdb)

	parsed_dist = parse_dist(dist_list)


	return parsed_dist


def parse_dist(dist_file):

	'''Parses the atomic distances file for duplications, hydrogens, nucleotide bases etc.'''

	nucleotide_bases = ['A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT']

	t_parsed_dist = []
	res_seen_list = []

	parsed_dist = {}

	chain_res = {}

	for line in dist_file:
		res_1, res_2, dist = line.split()[0], line.split()[1], line.split()[2].strip()
		chain_1, chain_2 = res_1[-1], res_2[-1]

		if chain_1 not in chain_res.keys():
			chain_res.update({chain_1 : {}})

		if chain_2 not in chain_res.keys():
			chain_res.update({chain_2 : {}})

		chain_list = [chain_1, chain_2]
		chain_list.sort()

		intersig = chain_list[0] + '-' + chain_list[1]
		if intersig not in parsed_dist.keys():
			parsed_dist.update({intersig : []})

		pair_el = res_order(line)

		if pair_el not in t_parsed_dist:
			t_parsed_dist.append(pair_el)

	for el in t_parsed_dist:
		res_1, res_2, dist = el.split()[0], el.split()[1], el.split()[2]
		sorted_ch = sorted([res_1[-1], res_2[-1]])
		inter_sig = sorted_ch[0] + '-' + sorted_ch[1]
		res_pair = [res_1, res_2]

		if res_1[0] == 'H' or res_2[0] == 'H':     # Excludes Hydrogen atoms
			pass
		elif res_1[0:3] == 'OXT' or res_2[0:3] == 'OXT':    # Excludes Terminal Oxygen atoms
			pass
		elif res_pair in res_seen_list:            # Remove duplications due to alternate side chain conformations and also multiple models
			pass
		elif res_1.split(':')[1] in nucleotide_bases or res_2.split(':')[1] in nucleotide_bases:  # Remove non-standard amino acids
			pass
		elif res_1.split(':')[1] == 'UNK' or res_2.split(':')[1] == 'UNK':  # Exclude unknown amino acids
			pass
		else:
			parsed_dist[inter_sig].append(el)
			res_seen_list.append(res_pair)
			resna_1, resna_2 = res_1.split(':')[1], res_2.split(':')[1]
			ressig_1 = res_1.split(':')[2] + ':' + res_1.split(':')[3]
			ressig_2 = res_2.split(':')[2] + ':' + res_2.split(':')[3]
			chain_1, chain_2 = ressig_1[-1], ressig_2[-1]
			if ressig_1 not in chain_res[chain_1].keys():
				chain_res[chain_1].update({ressig_1 : resna_1})
			if ressig_2 not in chain_res[chain_2].keys():
				chain_res[chain_2].update({ressig_2 : resna_2})

	return parsed_dist, chain_res


def res_order(line):

	'''Function that defines the residue order in a residue pair. The residues in a pair are lexicographically sorted.'''

	res_1, res_2, dist = line.split()[0], line.split()[1], line.split()[2].strip()
	resno_1, resno_2 = res_1.split(':')[2], res_2.split(':')[2]
	chain_1, chain_2 = res_1[-1], res_2[-1]	
	ressig_1, ressig_2 = resno_1 + ':' + chain_1, resno_2 + ':' + chain_2

	if len(resno_1) == len(resno_2):
		if max(ressig_1, ressig_2) == ressig_1:
			pair_el = res_1 + '\t' + res_2 + '\t' + dist

		else:
			pair_el = res_2 + '\t' + res_1 + '\t' + dist
			
	elif len(resno_1) < len(resno_2):
		if resno_1 in resno_2:
			if resno_2.index(resno_1) == 0:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
				
			else:
				if max(ressig_1, ressig_2) == ressig_1:
					pair_el = res_1 + '\t' + res_2 + '\t' + dist
					
				else:
					pair_el = res_2 + '\t' + res_1 + '\t' + dist
					
		else:
			if max(ressig_1, ressig_2) == ressig_1:
				pair_el = res_1 + '\t' + res_2 + '\t' + dist
				
			else:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
				
	elif len(resno_2) < len(resno_1):
		if resno_2 in resno_1:
			if resno_1.index(resno_2) == 0:
				pair_el = res_1  + '\t' + res_2 + '\t' + dist
				
			else:
				if max(ressig_1, ressig_2) == ressig_1:
					pair_el = res_1 + '\t' + res_2 + '\t' + dist
					
				else:
					pair_el = res_2 + '\t' + res_1 + '\t' + dist
					
		else:
			if max(ressig_1, ressig_2) == ressig_1:
				pair_el = res_1 + '\t' + res_2 + '\t' + dist
				
			else:
				pair_el = res_2 + '\t' + res_1 + '\t' + dist
	
	return pair_el


def altloc_check(input_file):

	'''Checks if the PDB has alternative atomic coordinates.'''

	input_pdb = open(input_file, 'r')
	input_lines = input_pdb.readlines()

	alt_status = 0

	for line in input_lines:
		if line[0:4] == 'ATOM':
			line = line.strip()
			alt = line[16]
			if alt != ' ':
				alt_status = 1
				break
			else:
				continue

	if alt_status == 1:
		filtered_pdb = altloc_filter(input_lines)
		return filtered_pdb
	else:
		return input_lines


def altloc_filter(input_lines):

	'''Filters the alternative atomic coordinates in the PDB file.'''

	filtered_pdb = []
	altlocs = {}

	for line in input_lines:
		if line[0:4] == 'ATOM':
			line = line.strip()
			resna = line[16:20].strip()
			if len(resna) > 3:
				atom_name = line[12:16]
				chain = line[21:22]
				resnum = line[22:26]
				occupancy = line[55:61]
				sig = ':'.join([atom_name, resnum, chain])
				if sig in altlocs.keys():
					if occupancy > altlocs[sig][0]:
						line = line.replace(line[16:20], ' ' + line[17:20])
						altlocs[sig][0] = occupancy
						altlocs[sig][1] = line
					else:
						pass
				else:
					line = line.replace(line[16:20], ' ' + line[17:20])
					altlocs.update({sig : [occupancy, line]})
			else:
				filtered_pdb.append(line)

	for element in altlocs.keys():
		filtered_pdb.append(altlocs[element][1])

	return filtered_pdb


def readpdb(inf):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>

	'''This routine reads the pdb file given as input and returns the coordinates of the corresponding atom identities.'''
	
	#read the pdb file#
	lines=inf;
	
	#initializing the maximum variable with negative infinity and  minimum value as positive infinity
	minx=10000000.0
	maxx=-10000000.0
	miny=10000000.0
	maxy=-10000000.0
	minz=10000000.0
	maxz=-10000000.0
	
	#Initiallizing the atom_props variable
	atom_props=[]
	
	#Iterating through each line
	for l in lines:
		#if line starts with 'ATOM' identifier, then continue
		if (l[0:4]=='ATOM'):
			#Getting the coordinates from pdb files
			xc=float(l[30:38]);
			yc=float(l[38:46]);
			zc=float(l[46:54]);
			
			#Getting the corresponding atom properties. Atom name, residue name,number, chain id
			# @Abhilesh modified - Hydrogen filter 
			if l.strip()[-1] != 'H':
				atom_type=l[12:16].strip();
				res_name=l[16:20].strip();
				res_num=l[22:26].strip();
				chain_num=l[21:22].strip();
			
			#Getting the max, min of x direction
				if minx>xc:
					minx=xc;
				elif maxx<xc:
					maxx=xc;
			
			#Getting the max, min of y direction
				if miny>yc:
					miny=yc;
				elif maxy<yc:
					maxy=yc;
			#Getting the max, min of z direction
				if minz>zc:
					minz=zc;
				elif maxz<zc:
					maxz=zc;
			
			#Append the atom_coord with each atom properties
				atom_props.append([xc,yc,zc,atom_type,res_name,res_num,chain_num]);

	#for element in atom_props:
	#   print element
	
	#Return            
	return minx,maxx,miny,maxy,minz,maxz,atom_props


def create_mesh(param,loc,cutoff):

	import numpy as np

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>
	'''This routine creates the 3D mesh, such that all the atoms in the pdb file are within this mesh.'''
	
	#Geting the maximum and minimum values of each coordinates from the param argument#
	minx=param[0]
	maxx=param[1]
	miny=param[2]
	maxy=param[3]
	minz=param[4]
	maxz=param[5]
	
	#Add the cutoff to the maximum value and subtract from the minimum value. This is done to extend the mesh dimension by cutoff distance in each direction.#
	minx-=cutoff
	maxx+=cutoff
	miny-=cutoff
	maxy+=cutoff
	minz-=cutoff
	maxz+=cutoff
	
	#loc is the length of cell. Here loc is equal to the cutoff value. This can ve a variable#
	#Getting the total number of cells in each dimension#
	xn=int((maxx-minx)/loc);
	yn=int((maxy-miny)/loc);
	zn=int((maxz-minz)/loc);
	
	#Setting the origin as minimum value in each dimension#
	origin=[minx,miny,minz]
	
	#Generating the equidistant points in x, y and z direction#
	xarray = np.linspace(minx,maxx,xn)
	yarray = np.linspace(miny,maxy,yn)
	zarray = np.linspace(minz,maxz,zn)
	
	#Setting up the mesh variable#
	mesh=[xarray,yarray,zarray,xn,yn,zn,origin,loc]
	
	#Print the MEsh Statistics#
	#print "\nMesh with following no of mesh points in (x,y,z) created: ","(",xn,yn,zn,")"
	#print "Mesh Origin: ", "(",origin[0],",",origin[1],",",origin[2],")"
	#print "Length of square cell/cutoff: ", loc
	#print "Total no of cells: ",xn*yn*zn,"\n"
	
	#Return the Mesh#
	return mesh


def recBinSearch(x, nums, low, high):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>
	'''This routine implements binary search for x in the list nums'''

	#Check if low value is greater than high value
	if low > high:
		return low
	#fing the mid point
	mid = (low + high) / 2
	
	item1 = nums[mid]
	item2 = nums[mid+1]
	#if range is found
	if ((item1 <= x) & (item2 > x)):
		return mid
	#Else do the binary search on left and right parts
	elif x < item1:
		return recBinSearch(x, nums, low, mid-1)
	else: 
		return recBinSearch(x, nums, mid+1, high)


def bsearch(x, nums):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>
	'''Calls the binary search routine'''

	return recBinSearch(x, nums, 0, len(nums)-1)


def assign_atomlist_to_mesh(mesh,atom_props,cell_lst):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>
	'''This routine assigns each atom in the 3D mesh. This is done through binary search on the x,y and z arrays of mesh points.'''
	
	#get the total number of mesh points in each direction#
	xn=mesh[3];
	yn=mesh[4];
	zn=mesh[5];
	
	#storing coordinates in x-y-z major. ie. z coordinates first then y then x
	for c in atom_props:
		#Do a binary search of current x,y,z coordinates in the corresponding x,y,z array of points#
		ix=bsearch(c[0],mesh[0])
		iy=bsearch(c[1],mesh[1])
		iz=bsearch(c[2],mesh[2])
		
		#z-y-x major
		index=zn*yn*ix + iy*zn +iz;
		
		#append the mesh cell with coordinate & its properties
		cell_lst[index].append(c);
		#Return cell list
	return cell_lst


def assign_ngh(mesh):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>
	'''This routine calculates the distance between all the atoms in the cell)_list variable'''

	import math

	#get the total number of points in each direction
	xn=mesh[3];
	yn=mesh[4];
	zn=mesh[5];
	
	#initializing the neighbor dictionary
	nghdict={};
	
	#Calculating the total number of cells
	totcell=xn*yn*zn;
	
	#Iterating through total number of cells
	for lstidx in range(0,totcell):
		
		#calculate the x,y,z index using the mesh indexes. This is possible because the mesh indexes are z then y then x major. Here we are unwinding the mesh in x,y and z direction
		r=lstidx%(yn*zn);
		#get mesh x index
		mix=lstidx/(yn*zn);
		#get mesh z index
		miz=r%zn;
		#get mesh y index
		miy=r/zn;
		
		#initializing the temp neighbor list variable
		temp_nghlst=[];
		
		#All permuations/neighbors are taken by adding -1,0,-1 values in each direction. here for each cell in the mesh, its neighbors are recorded in nghdict variable
		for idsx in -1,0,1:
			for idsy in -1,0,1:
				for idsz in -1,0,1:
					#Get index in each dimension
					xid,yid,zid=(mix+idsx),(miy+idsy),(miz+idsz);
					#Check if each index is less than the maximum number of cells in that direction and if each index is greater than -1
					if((xid < xn)&(xid > -1)):
						if((yid < yn)&(yid > -1)):
							if((zid < zn)&(zid > -1)):
								#Recalculate the mesh index
								index=zn*yn*xid + yid*zn +zid;
								#Append the mesh index in the nighbor list 
								temp_nghlst.append(index)
		
		#Add the temp neighbor list in the dictionary
		nghdict[lstidx]=temp_nghlst;
		
	return nghdict


def calculate_distance_dm(cell_list,cmesh,nghdict,cutoff):

	import math
	
	#get the total number of cells in each dircetion and the total numbers
	cxn=cmesh[3]
	cyn=cmesh[4]
	czn=cmesh[5]
	cmax=cxn*cyn*czn
	
	#Calculates the square of the cutoff
	square_cutoff=cutoff*cutoff;
	outline=""
	
	#Iterate through all cells in the Mesh
	for i in range(0,cmax):
		#Get the coordinate set 1
		coord_set1=cell_list[i];
		#get the neighbors indexs for the current cell
		for idx in nghdict[i]:
			#Get the current coordinate set 2
			coord_set2=cell_list[idx];
			
			#For each atom in the coordinate sets 1,2 calculate the distance
			for coord1 in coord_set1:
				for coord2 in coord_set2:
					# @Abhilesh modified - If clause added
					if coord1[-1] != coord2[-1]:
					#Calculate the square of the distance
						square_distance=compute_distance(coord1,coord2)
					#If square of the distance is less than the square of the cutoff distnace, then calculate the distance by taking the square root.
						if (square_distance < square_cutoff):
							distance=math.sqrt(square_distance)
						
						#Format the string to right in the output file
							outline+=str(coord1[3])+":"+str(coord1[4])+":"+str(coord1[5])+":"+str(coord1[6])+"\t"+str(coord2[3])+":"+str(coord2[4])+":"+str(coord2[5])+":"+str(coord2[6])+"\t"+str(distance)+"\n";

	#return the output line
	return outline

# @Abhilesh modified - Function duplicated from above and a few lines
def calculate_distance_mm(cell_list,cmesh,nghdict,cutoff,p1,p2):

	import math
	
	#get the total number of cells in each dircetion and the total numbers
	cxn=cmesh[3]
	cyn=cmesh[4]
	czn=cmesh[5]
	cmax=cxn*cyn*czn
	
	#Calculates the square of the cutoff
	square_cutoff=cutoff*cutoff;
	outline=""
	
	#Iterate through all cells in the Mesh
	for i in range(0,cmax):
		#Get the coordinate set 1
		coord_set1=cell_list[i];
		#get the neighbors indexs for the current cell
		for idx in nghdict[i]:
			#Get the current coordinate set 2
			coord_set2=cell_list[idx];
			
			#For each atom in the coordinate sets 1,2 calculate the distance
			for coord1 in coord_set1:
				for coord2 in coord_set2:
					# @Abhilesh modified - Condition added
					if coord1[-1] != coord2[-1]:
						# @Abhilesh modified - Condition added
						if (coord1[-1] == p1 and coord2[-1] == p2) or (coord1[-1] == p2 and coord2[-1] == p1):
					#Calculate the square of the distance
							square_distance=compute_distance(coord1,coord2)
					#If square of the distance is less than the square of the cutoff distnace, then calculate the distance by taking the square root.
							if (square_distance < square_cutoff):
								distance=math.sqrt(square_distance)
						
						#Format the string to right in the output file
								outline+=str(coord1[3])+":"+str(coord1[4])+":"+str(coord1[5])+":"+str(coord1[6])+"\t"+str(coord2[3])+":"+str(coord2[4])+":"+str(coord2[5])+":"+str(coord2[6])+"\t"+str(distance)+"\n";

	#return the output line
	return outline


def compute_distance(coord1,coord2):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com
	'''This routine calculates the distance between any two coordinates'''

	return ((coord1[0]-coord2[0])*(coord1[0]-coord2[0]) + (coord1[1]-coord2[1])*(coord1[1]-coord2[1]) + (coord1[2]-coord2[2])*(coord1[2]-coord2[2]))


def get_dist(input_pdb, p1 = None, p2 = None):

	#Copyright 2013 Neelesh Soni <neelrocks4@gmail.com>; @Abhilesh modified - Flags for specifying proteins
	
	'''Main routine for cell_list implementation for calculating pairwise distances.'''

	#Open the input file for parsing.#
	inf2=input_pdb;
	#parse_pdb function returns the details of all atoms in the pdb file with maximum and minimum value of coordinates in each direction. All these values gets stored in the 'param' list#
	param=readpdb(inf2);
	#close the file#
	#inf2.close();

	#Set the distance cutoff#
	cutoff=8.0;

	#Set the coarse mesh cell length in angstrom#
	cmesh_len=cutoff;

	#function 'create_mesh' creates a 3D mesh such that all atoms should be contained within this mesh. Extra 'cutoff' length is added in each direction for boundary cells#
	cmesh=create_mesh(param,cmesh_len,cutoff);

	#getting the all atom coordinates and their identities like atom type, residue name, residue number and chain identifier. All this are stored int he 6th element of param list#
	atom_props=param[6];

	#Getting the number of total number cell in the mesh#
	totcell=cmesh[3]*cmesh[4]*cmesh[5];

	#initializing the cell_list variable. This stores the list of all atoms in each cell of the mesh#
	cell_list=[ [] for n in range(0,totcell) ]

	#Assign the atoms in each cell of the mesh and store them in cell_list#
	cell_list=assign_atomlist_to_mesh(cmesh,atom_props,cell_list);

	#get the neighbors of each cell in the mesh and store them in the neighbor dictionary "nghdict" #
	ngh_dict=assign_ngh(cmesh);

	#Calculate the dictance of all the atoms that are withion cutoff distance. This function stores the distances with atom identities in the outline string.#
	# @Abhilesh modified - If clause added
	if p1 == None and p2 == None:
		outline=calculate_distance_dm(cell_list,cmesh,ngh_dict,cutoff)
	else:
		outline=calculate_distance_mm(cell_list,cmesh,ngh_dict,cutoff,p1,p2)

	#Open the output file for writing the distances#
	dist_list = []
	#Write the outline variable to the file#
	for element in outline.split('\n')[:-1]:
		dist_list.append(element + '\n')

	return dist_list


def get_contacts(parsed_dist, interactions_filter):

	'''Takes the parsed dist and the interactions_filter and returns the interacting atoms per residue pair accordingly'''

	interface_respairs = {}
	respairs = {}
	interface_res = {}
	num_interactors = {}

	main_chain_atoms = ['C', 'CA', 'O', 'OXT', 'N']

	# For interactions involving all the atoms of both the residues
	if interactions_filter == 'all':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if res_pair not in interface_respairs[interface].keys():
					interface_respairs[interface].update({res_pair : [line]})
				else:
					interface_respairs[interface][res_pair].append(line)
				resna_1, resna_2 = res_1[0:3], res_2[0:3]
				interface_res[chain_1].update({ressig_1 : resna_1})
				interface_res[chain_2].update({ressig_2 : resna_2})
				respair = ressig_1 + '-' + ressig_2
				respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the side chains of both the residues
	elif interactions_filter == 'ss':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2
			
				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})
			
				res_pair = res_1 + '-' + res_2
				if atom_1 not in main_chain_atoms and atom_2 not in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the main chain atoms of both the residues
	elif interactions_filter == 'mm':
		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if atom_1 in main_chain_atoms and atom_2 in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})

	# For interactions involving the main chain atoms of residue 1 and the side chain atoms of residue 2
	elif interactions_filter == 'ms':

		flag_dict = {}

		for interface in parsed_dist.keys():

			if interface not in interface_respairs.keys():
				interface_respairs.update({interface : {}})

			if interface not in respairs.keys():
				respairs.update({interface : {}})

			for line in parsed_dist[interface]:
				res_1, res_2 = ':'.join(line.split()[0].split(':')[1:]), ':'.join(line.split()[1].split(':')[1:])
				atom_1, atom_2 = line.split()[0].split(':')[0], line.split()[1].split(':')[0]
				chain_1, chain_2 = res_1[-1], res_2[-1]
				ressig_1, ressig_2 = res_1.split(':')[1] + ':' + chain_1, res_2.split(':')[1] + ':' + chain_2

				if chain_1 not in interface_res.keys():
					interface_res.update({chain_1 : {}})
				if chain_2 not in interface_res.keys():
					interface_res.update({chain_2 : {}})

				res_pair = res_1 + '-' + res_2
				if atom_1 in main_chain_atoms and atom_2 not in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})
					flag_dict.update({respair : '1'})
				elif atom_1 not in main_chain_atoms and atom_2 in main_chain_atoms:
					if res_pair not in interface_respairs[interface].keys():
						interface_respairs[interface].update({res_pair : [line]})
					else:
						interface_respairs[interface][res_pair].append(line)
					resna_1, resna_2 = res_1[0:3], res_2[0:3]
					interface_res[chain_1].update({ressig_1 : resna_1})
					interface_res[chain_2].update({ressig_2 : resna_2})
					respair = ressig_1 + '-' + ressig_2
					respairs[interface].update({respair : resna_1 + '-' + resna_2})
					flag_dict.update({respair : '2'})

	for interface in respairs.keys():
		num_interactors.update({interface : len(respairs[interface])})

	if interactions_filter != 'ms':
		return interface_res, respairs, interface_respairs, num_interactors
	else:
		return interface_res, respairs, interface_respairs, num_interactors, flag_dict


def get_weights(interface_respairs, options, flag_dict = None):

	'''Calculates the weights for the interacting residue pairs, the weight being equal to the number of interacting atoms
	divided by the total number of atoms for that amino acid'''

	inter_type = options['intertype']

	aa_all_atoms = {
		'GLY' : 4,
		'PRO' : 9, 
		'ALA' : 5, 
		'VAL' : 7,
		'LEU' : 8, 
		'ILE' : 8, 
		'MET' : 8, 
		'CYS' : 6, 
		'PHE' : 11, 
		'TYR' : 12, 
		'TRP' : 14,
		'HIS' : 10, 
		'LYS' : 9, 
		'ARG' : 11, 
		'GLN' : 9, 
		'ASN' : 8, 
		'GLU' : 9, 
		'ASP' : 8, 
		'SER' : 6, 
		'THR' : 7
		}

	aa_sc_atoms = {
		'GLY' : 0,
		'PRO' : 5, 
		'ALA' : 1, 
		'VAL' : 3,
		'LEU' : 4, 
		'ILE' : 4, 
		'MET' : 4, 
		'CYS' : 2, 
		'PHE' : 7, 
		'TYR' : 8, 
		'TRP' : 10,
		'HIS' : 6, 
		'LYS' : 5, 
		'ARG' : 7, 
		'GLN' : 5, 
		'ASN' : 4, 
		'GLU' : 5, 
		'ASP' : 4, 
		'SER' : 2, 
		'THR' : 3
		}

	aa_mc_atoms = {
		'GLY' : 4,
		'PRO' : 4, 
		'ALA' : 4, 
		'VAL' : 4,
		'LEU' : 4, 
		'ILE' : 4, 
		'MET' : 4, 
		'CYS' : 4, 
		'PHE' : 4, 
		'TYR' : 4, 
		'TRP' : 4,
		'HIS' : 4, 
		'LYS' : 4, 
		'ARG' : 4, 
		'GLN' : 4, 
		'ASN' : 4, 
		'GLU' : 4, 
		'ASP' : 4, 
		'SER' : 4, 
		'THR' : 4
		}

	main_chain_atoms = ['N', 'C', 'CA', 'O', 'OXT']

	interaction_respairs = {}

	for interface in interface_respairs.keys():
		if interface not in interaction_respairs.keys():
			interaction_respairs.update({interface : {}})
		for element in interface_respairs[interface].keys():
			if element.count('-') == 1:
				res_1, res_2 = element.split('-')[0], element.split('-')[1]
			elif element.count('-') > 1:
				res_1 = ':'.join([element.split(':')[0], element.split(':')[1], element.split(':')[2][0]]) 
				res_2 = ':'.join([element.split(':')[2][-3:], element.split(':')[3], element.split(':')[4]])
			res_1_int_atoms = set()
			res_2_int_atoms = set()
			for item in interface_respairs[interface][element]:
				res_a, res_b = ':'.join(item.split()[0].split(':')[1:]), ':'.join(item.split()[1].split(':')[1:])
				atom_a, atom_b = item.split()[0].split(':')[0], item.split()[1].split(':')[0]
				if res_a == res_1 and res_b == res_2:
					res_1_int_atoms.add(atom_a), res_2_int_atoms.add(atom_b)
				elif res_a == res_2 and res_b == res_1:
					res_1_int_atoms.add(atom_b), res_2_int_atoms.add(atom_a)
				interaction_respairs[interface].update({element : [res_1, len(res_1_int_atoms), res_2, len(res_2_int_atoms)]})

	weight_dict = {}

	for interface in interaction_respairs.keys():
		if interface not in weight_dict.keys():
			weight_dict.update({interface : {}})
		for element in interaction_respairs[interface].keys():
			if element.count('-') == 1:
				res_type_a, res_type_b = element.split('-')[0].split(':')[0], element.split('-')[1].split(':')[0]
				resno_a, resno_b = element.split('-')[0].split(':')[1], element.split('-')[1].split(':')[1]
				chain_a, chain_b = element.split('-')[0][-1], element.split('-')[1][-1]
			elif element.count('-') > 1:
				res_type_a, res_type_b = element[0:3], element.split(':')[2][-3:]
				resno_a, resno_b = element.split(':')[1], element.split(':')[3]
				chain_a, chain_b = element.split(':')[2][0], element.split(':')[-1]
			ressig_a, ressig_b = resno_a + ':' + chain_a, resno_b + ':' + chain_b
			res_1_int_atoms, res_2_int_atoms = interaction_respairs[interface][element][1], interaction_respairs[interface][element][3]

			try:
				aa_all_atoms[res_type_a]
				aa_all_atoms[res_type_b]
				
				if inter_type == 'all':
					res_1_total_atoms, res_2_total_atoms = aa_all_atoms[res_type_a], aa_all_atoms[res_type_b]
				elif inter_type == 'ss':
					res_1_total_atoms, res_2_total_atoms = aa_sc_atoms[res_type_a], aa_sc_atoms[res_type_b]
				elif inter_type == 'mm':
					res_1_total_atoms, res_2_total_atoms = aa_mc_atoms[res_type_a], aa_mc_atoms[res_type_b]
				elif inter_type == 'ms':
					flag = flag_dict[ressig_a + '-' + ressig_b]
					if flag == '1':
						res_1_total_atoms, res_2_total_atoms = aa_mc_atoms[res_type_a], aa_sc_atoms[res_type_b]
					elif flag == '2':
						res_1_total_atoms, res_2_total_atoms = aa_sc_atoms[res_type_a], aa_mc_atoms[res_type_b]

				weight = min((float(res_1_int_atoms)/float(res_1_total_atoms)), float(res_2_int_atoms)/float(res_2_total_atoms))
				weight = round_val(weight)
				weight_dict[interface].update({ressig_a + '-' + ressig_b : weight})

			except KeyError:
				pass

	return weight_dict


def round_val(value):

	'''Rounds off the given value to 3 decimal places'''

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num


def get_linkages(dist_file):

	'''Get the interacting linkages for each interface'''

	linkage_pairs = []
	sub1_list = set()
	sub2_list = set()

	for interface in dist_file.keys():
		for element in dist_file[interface]:
			res_1 = ':'.join(element.split()[0].split(':')[1:])
			res_2 = ':'.join(element.split()[1].split(':')[1:])
			#res_pair_list = sorted([res_1, res_2], key = lambda t:t[-1])
			#res_pair = res_pair_list[0] + '-' + res_pair_list[1]
			res_pair = res_1 + '-' + res_2
			sub1_list.add(res_1)
			sub2_list.add(res_2)

			if res_pair not in linkage_pairs:
				linkage_pairs.append(res_pair)

	return linkage_pairs, sub1_list, sub2_list


def call_score_func(chain_res, respairs, weight_dict, pot_dict):

	'''Main module to score the native and the randomised structures'''

	score_native = score_complex(respairs, weight_dict, pot_dict)
	score_decoys = get_rand_scores(chain_res, respairs, weight_dict, pot_dict)

	return score_native, score_decoys


def score_complex(respairs, weight_dict, pot_dict):

	'''Scores the native interface'''

	native_scores = {}

	for interface in weight_dict.keys():
		if len(weight_dict[interface]) != 0:
			score_native = 0.0
			for pair_el in weight_dict[interface]:
				resna_pair = respairs[interface][pair_el]
				weight = weight_dict[interface][pair_el]
				score = pot_dict[resna_pair] * weight
				score_native += score

			score_native = round_val(score_native)
			native_scores.update({interface : score_native})
		else:
			pass

	return native_scores


def get_rand_scores(chain_res, respairs, weight_dict, pot_dict):

	'''Randomises the native interface and then scores them'''

	# This module breaks apart the interface and replaces each amino acid on the interface 
	# with another amino acid from the subunit, effectively changing the identity of the interface
	# and randomizing amino acid interactions across the interface.

	# Note: The integrity of the interface is not maintained in this algorithm

	from random import randint

	ranres = {}
	t_ranres = {}
	score_decoys = {}

	for chain in chain_res.keys():
		ranres.update({chain : {}})
		for res in chain_res[chain].keys():
			ranres[chain].update({res : chain_res[chain][res]})

	for interface in respairs.keys():
		if len(respairs[interface]) != 0:
			score_decoys.update({interface : {}})

	for counter in range(0, 1000):
		for chain in ranres.keys():
			resnalist = ranres[chain].values()
			resnolist = ranres[chain].keys()
			for i in range(len(resnalist) - 1, -1, -1):
				j = randint(0, i)
				resnalist[i], resnalist[j] = resnalist[j], resnalist[i]
			for k in range(0, len(resnolist)):
				t_ranres.update({resnolist[k] : resnalist[k]})
		r_score_el = score_rand(t_ranres, respairs, weight_dict, pot_dict)
		for el in r_score_el:
			score_decoys[el[0]].update({counter : el[1]})

	return score_decoys


def score_rand(t_ranres, respairs, weight_dict, pot_dict):

	'''Scores the randomized interfaces'''

	r_scores = []

	for interface in respairs.keys():
		if len(respairs[interface]) != 0:
			decoy_score = 0.000
			for pair in respairs[interface]:
				if pair.count('-') == 1:
					pa_res1, pa_res2 = pair.split('-')[0], pair.split('-')[1]
				else:
					split_index = pair.find('-')
					if split_index != 0:
						pa_res1, pa_res2 = pair[:split_index], pair[split_index + 1:]
					else:
						split_index = (pair.replace('-', '', 1).find('-')) + 1
						pa_res1, pa_res2 = pair[:split_index], pair[split_index + 1:]
				n_res1, n_res2 = t_ranres[pa_res1], t_ranres[pa_res2]
				n_pair = n_res1 + '-' + n_res2
				if n_res1 == 'GLY' or n_res2 == 'GLY':
					pass
				else:
					weight = weight_dict[interface][pair]
					pair_score = pot_dict[n_pair] * weight
					decoy_score += pair_score

			decoy_score = round_val(decoy_score)

			r_score_el = (interface, decoy_score)
			r_scores.append(r_score_el)

	return r_scores


def calc_zscore(score_native, score_decoys):

	'''Calculates the z-score and reports other metrics'''

	from math import sqrt
	import sys

	outvals = {}

	for interface in score_decoys.keys():
		total = 0.000
		false_pos = 0
		for score_el in score_decoys[interface].values():
			total += score_el

			if score_el <= score_native[interface]:
				false_pos += 1

			if score_el > score_native[interface]:
				try:
					min_tn
				except NameError:
					min_tn = score_el
				else:
					if score_el <= min_tn:
						min_tn = score_el

			try:
				min_score
			except NameError:
				min_score = score_el
			else:
				if score_el < min_score:
					min_score = score_el

		false_pos_rate = false_pos / len(score_decoys[interface].values())
		avg_score = total / len(score_decoys[interface].values())

		std_dev = 0
		for score_el in score_decoys[interface].values():
			std_dev += (score_el - avg_score) ** 2
	
		std_dev = std_dev / len(score_decoys[interface].values())

		if std_dev == 0:
			print "Error: stdev of 0"
			sys.exit()

		std_dev = sqrt(std_dev)

		z_bg = []
		for score_el in score_decoys[interface].values():
			z_scr = (score_el - avg_score) / std_dev
			z_bg.append(z_scr)

		z_score = (score_native[interface] - avg_score) / std_dev

		try:
			min_tn
		except NameError:
			z_min_tn = 'undef'
			z_primer = 'undef'
		else:
			z_min_tn = (min_tn - avg_score) / std_dev
			z_prime = z_score - z_min_tn

		z_min = (min_score - avg_score) / std_dev
		z_2 = z_score - z_min

		z_score = round_val(z_score)
		avg_score = round_val(avg_score)
		std_dev = round_val(std_dev)

		outvals.update({interface : (z_score, avg_score, std_dev)})

	return outvals


def call_alascan(pot_dict, weight_dict, interface_res, native_score, mode = ['1'], res_list = None):

	'''Main module to run Alanine Scanning in the mode specified by the user'''

	import sys

	if mode == ['1']:
		alscn_results = alascan(pot_dict, weight_dict, interface_res, native_score)
	elif mode == ['2']:
		alscn_results = alascan_all(pot_dict, weight_dict, interface_res, native_score)
	elif mode == ['3']:
		print "Error: Please provide Residue list for Alanine Scanning"
		sys.exit()
	elif len(mode) == 2 and mode[0] == '3':
		alscn_results = mutate_reslist(pot_dict, weight_dict, interface_res, res_list, native_score)

	return alscn_results


def alascan(pot_dict, weight_dict, interface_res, native_score):

	'''Mutates all the interface residues to Alanine and computes the mutational effect on the complex'''

	decoy_score_dict = {}
	res_mutate = {}

	for subunit in interface_res.keys():
		for res in interface_res[subunit].keys():
			for interface in weight_dict.keys():

				decoy_score = 0.0

				try:
					decoy_score_dict[interface]
				except KeyError:
					decoy_score_dict.update({interface : {}})

				for res_pair in weight_dict[interface].keys():
					if res in res_pair:
						res_index = res_pair.index(res)

						if res_pair.count('-') == 1:
							res_1 = res_pair.split('-')[0]
							res_2 = res_pair.split('-')[1]
						else:
							split_index = res_pair.find('-')
							if split_index != 0:
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]
							else:
								split_index = (res_pair.replace('-', '', 1).find('-')) + 1
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]

						resna_1 = interface_res[res_1[-1]][res_1]
						resna_2 = interface_res[res_2[-1]][res_2]

						if res_index == 0:
							sub_respair = 'ALA' + '-' + resna_2
						else:
							sub_respair = resna_1 + '-' + 'ALA'

						pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
						decoy_score += pair_score

					else:
						if res_pair.count('-') == 1:
							res_1 = res_pair.split('-')[0]
							res_2 = res_pair.split('-')[1]
						else:
							split_index = res_pair.find('-')
							if split_index != 0:
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]
							else:
								split_index = (res_pair.replace('-', '', 1).find('-')) + 1
								res_1 = res_pair[:split_index]
								res_2 = res_pair[split_index + 1:]

						resna_1 = interface_res[res_1[-1]][res_1]
						resna_2 = interface_res[res_2[-1]][res_2]

						orig_respair = resna_1 + '-' + resna_2

						pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
						decoy_score += pair_score

				decoy_score_dict[interface].update({res : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		for el_res in decoy_score_dict[interface].keys():
			r_decoy_score = decoy_score_dict[interface][el_res]
			if r_decoy_score > native_score[interface]:
				res_mutate.update({el_res : r_decoy_score - native_score[interface]})

	return res_mutate


def alascan_all(pot_dict, weight_dict, interface_res, native_score):

	'''Mutates all the interface residues to all the other amino acids and computes the mutational effect on the complex'''

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

	decoy_score_dict = {}
	res_mutate = {}

	for subunit in interface_res.keys():
		for res in interface_res[subunit].keys():
			for aa in amino_acids:
				for interface in weight_dict.keys():

					decoy_score = 0.0

					try:
						decoy_score_dict[interface]
					except KeyError:
						decoy_score_dict.update({interface : {}})

					try:
						decoy_score_dict[interface][res]
					except KeyError:
						decoy_score_dict[interface].update({res : {}})

					for res_pair in weight_dict[interface].keys():
						if res in res_pair:
							res_index = res_pair.index(res)

							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							if res_index == 0:
								if aa != resna_1:
									sub_respair = aa + '-' + resna_2
							else:
								if aa != resna_2:
									sub_respair = resna_1 + '-' + aa

							pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
							decoy_score += pair_score

						else:
							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							orig_respair = resna_1 + '-' + resna_2

							pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
							decoy_score += pair_score

					try: 
						decoy_score_dict[interface][res][aa]
					except KeyError:
						decoy_score_dict[interface][res].update({aa : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		for res in decoy_score_dict[interface].keys():
			score_diff = 0.0
			for aa in decoy_score_dict[interface][res].keys():
				decoy_score = decoy_score_dict[interface][res][aa]
				score_diff += decoy_score - native_score[interface]
			res_mutate.update({res : score_diff})


	return res_mutate


def mutate_reslist(pot_dict, weight_dict, interface_res, res_list, native_score):

	'''Mutates the residues submitted by the user to all other amino acids and computes their mutational effect on the complex'''

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

	decoy_score_dict = {}
	res_mutate = {}

	for subunit in res_list.keys():
		for res in res_list[subunit].keys():
			for aa in amino_acids:
				for interface in weight_dict.keys():

					decoy_score = 0.0

					try:
						decoy_score_dict[interface]
					except KeyError:
						decoy_score_dict.update({interface : {}})

					try:
						decoy_score_dict[interface][res]
					except KeyError:
						decoy_score_dict[interface].update({res : {}})

					for res_pair in weight_dict[interface].keys():
						if res in res_pair:
							res_index = res_pair.index(res)

							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							if res_index == 0:
								resna_1 = res_list[res_1[-1]][res_1]
								resna_2 = interface_res[res_2[-1]][res_2]
							else:
								resna_1 = interface_res[res_1[-1]][res_1]
								resna_2 = res_list[res_2[-1]][res]

							if res_index == 0:
								if aa != resna_1:
									sub_respair = aa + '-' + resna_2
							else:
								if aa != resna_2:
									sub_respair = resna_1 + '-' + aa

							pair_score = weight_dict[interface][res_pair] * pot_dict[sub_respair]
							decoy_score += pair_score

						else:
							if res_pair.count('-') == 1:
								res_1 = res_pair.split('-')[0]
								res_2 = res_pair.split('-')[1]
							else:
								split_index = res_pair.find('-')
								if split_index != 0:
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]
								else:
									split_index = (res_pair.replace('-', '', 1).find('-')) + 1
									res_1 = res_pair[:split_index]
									res_2 = res_pair[split_index + 1:]

							resna_1 = interface_res[res_1[-1]][res_1]
							resna_2 = interface_res[res_2[-1]][res_2]

							orig_respair = resna_1 + '-' + resna_2

							pair_score = weight_dict[interface][res_pair] * pot_dict[orig_respair]
							decoy_score += pair_score

					try: 
						decoy_score_dict[interface][res][aa]
					except KeyError:
						decoy_score_dict[interface][res].update({aa : round_val(decoy_score)})

	for interface in decoy_score_dict.keys():
		for res in decoy_score_dict[interface].keys():
			score_diff = 0.0
			for aa in decoy_score_dict[interface][res].keys():
				decoy_score = decoy_score_dict[interface][res][aa]
				score_diff += decoy_score - native_score[interface]
			res_mutate.update({res : score_diff})


	return res_mutate


def get_reslist(filename):

	'''This routine parses the residue file for alanine scanning provided by user'''

	import os

	fn = open(os.getcwd() + '/input/' + filename, 'r')

	reslist = {}

	for line in fn.readlines():
		line = line.strip()
		split_res = line.split(':')
		resno = split_res[0]
		resna = split_res[1]
		chain = split_res[2]

		try:
			reslist[chain]
		except KeyError:
			reslist.update({chain : {}})

		reslist[chain].update({resno + ':' + chain : resna})

	return reslist


def write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, respairs, alscn = None):

	'''Writes out the output files'''

	import sys

	pdb_name = options['input_pdb'].split('/')[-1].replace('.pdb', '')

	out_score_file = output_file.replace('.out', '_scores.out')
	out_respairs_file = output_file.replace('.out', '_respairs.out')
	out_alascan_file = output_file.replace('.out', '_alascan.out')


	# Write out the main values (z-score, raw score etc.) in the score file
	sys.stdout = open(out_score_file, 'w')

	for interface in out_val.keys():
		
		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		print "{0:<31s} {1:1s} {2:<10s}".format('Input_pdb', ':', pdb_name + '.pdb')
		print "{0:<31s} {1:1s} {2:<15s}".format('Interface', ':', pdb_name + '_' + p1 + '-' + pdb_name + '_' + p2)
		print "{0:<31s} {1:1s} {2:<5d}".format('Number of Interacting Pairs', ':', num_interactors[interface])
		print "{0:<31s} {1:1s} {2:<2.2f}".format('Distance Threshold', ':', options['cutoff'])
		print "{0:<31s} {1:1s} {2:<2s}".format('Potential Type', ':', options['intertype'])
		print "{0:<31s} {1:1s} {2:<6.3f}".format('Raw Score', ':', native_score[interface])
		print "{0:<31s} {1:1s} {2:<3.3f}".format('Z-score', ':', out_val[interface][0])
		print "{0:<31s} {1:1s} {2:<6.3f}".format('Average Background Score', ':', out_val[interface][1])
		print "{0:<31s} {1:1s} {2:<6.3f}".format('Standard Deviation', ':', out_val[interface][2])
		print "{0:<31s} {1:1s} {2:<3.3f}".format('Z-score threshold', ':', z_threshold)
		print "\n"
		if out_val[interface][0] < z_threshold:
			print "Congratulations! It BINDS!"
		else:
			print "Oops, Seems like these two proteins do not like each other!"
		print "\n"


	# Write out the interacting residue pairs for each interface in the respairs file
	sys.stdout = open(out_respairs_file, 'w')

	for interface in respairs.keys():

		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		print "{0:<10s} {1:1s} {2:<15s}".format('Interface', ':', pdb_name + '_' + p1 + '-' + pdb_name + '_' + p2)
		print "\n"

		for res_pair in respairs[interface].keys():
			if res_pair.count('-') == 1:
				res_1 = res_pair.split('-')[0]
				res_2 = res_pair.split('-')[1]
			else:
				split_index = res_pair.find('-')
				if split_index != 0:
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
				else:
					split_index = (res_pair.replace('-', '', 1).find('-')) + 1
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]

			resno_1, resno_2 = res_1.split(':')[0], res_2.split(':')[0]
			chain_1, chain_2 = res_1.split(':')[1], res_2.split(':')[1]
			resna = respairs[interface][res_pair].split('-')

			resna_1, resna_2 = resna[0], resna[1]

			out_res_1 = resno_1 + ':' + resna_1 + ':' + chain_1
			out_res_2 = resno_2 + ':' + resna_2 + ':' + chain_2

			print "{0:<15s} {1:4s} {2:<15s}".format(out_res_1, "    ", out_res_2)

		print "\n"


	# Write out the alanine scanning predictions in the alascan file
	if alscn != None:

		sys.stdout = open(out_alascan_file, 'w')

		print "{0:<5s} {1:3s} {2:<11s} {3:3s} {4:<6s} {5:3s} {6:<10s}".format('Rank', '   ', 'Residue no.', '   ', 'Chain', '   ', 'Score Diff')

		alscn_list = sorted(alscn, key = alscn.get, reverse = True)
		for element in alscn_list:
			rank = alscn_list.index(element) + 1
			res_no = element.split(':')[0]
			chain = element.split(':')[1]
			diff = alscn[element]
			print "{0:3d} {1:3s} {2:>11s} {3:3s} {4:>6s} {5:3s} {6:>10.3f}".format(rank, '   ', res_no, '   ', chain, '   ', diff)


	# Write out the respair files for our webserver developer
	for interface in respairs.keys():
	
		sys.stdout = open(output_file.replace('.out', '_' + interface + '_respairs.out'), 'w')

		interface_res = set()

		p1 = interface.split('-')[0]
		p2 = interface.split('-')[1]

		for res_pair in respairs[interface].keys():
			if res_pair.count('-') == 1:
				res_1 = res_pair.split('-')[0]
				res_2 = res_pair.split('-')[1]
			else:
				split_index = res_pair.find('-')
				if split_index != 0:
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]
				else:
					split_index = (res_pair.replace('-', '', 1).find('-')) + 1
					res_1 = res_pair[:split_index]
					res_2 = res_pair[split_index + 1:]

			resno_1, resno_2 = res_1.split(':')[0], res_2.split(':')[0]
			chain_1, chain_2 = res_1.split(':')[1], res_2.split(':')[1]
			resna = respairs[interface][res_pair].split('-')

			resna_1, resna_2 = resna[0], resna[1]

			out_res_1 = resno_1 + ':' + resna_1 + ':' + chain_1
			out_res_2 = resno_2 + ':' + resna_2 + ':' + chain_2

			interface_res.add(out_res_1)
			interface_res.add(out_res_2)

		for element in interface_res:
			print "{0:<15s}".format(element)


		sys.stdout = open("/dev/stdout", "w")

 
	sys.stdout = open("/dev/stdout", "w")

	return


def predict_binding_interface(input_file, pot_dict, z_threshold, output_file, options, p1 = None, p2 = None, f_alascan = None):

	'''Main module to run the software for dimers or for specific interfaces in a multimeric complex'''

	import sys
	import cPickle as pickle

	pdb_name = output_file.split('/')[-1].replace('.out', '')

	filtered_pdb = altloc_check(input_file)

	dist_list = get_dist(filtered_pdb, p1, p2)

	parsed_dist = parse_dist(dist_list)

	interface_residue_pairs = get_contacts(parsed_dist[0], options['intertype'])

	if options['intertype'] == 'ms':
		weight_dict = get_weights(interface_residue_pairs[2], options, flag_dict = interface_residue_pairs[4])
	else:
		weight_dict = get_weights(interface_residue_pairs[2], options)

	num_interactors = interface_residue_pairs[3]

	scores = call_score_func(parsed_dist[1], interface_residue_pairs[1], weight_dict, pot_dict)

	native_score = scores[0]
	decoy_scores = scores[1]

	out_val = calc_zscore(native_score, decoy_scores)

	linkage_pairs = get_linkages(parsed_dist[0])

	if f_alascan == None:
		pass
	elif len(f_alascan) == 2:
		reslist = get_reslist(f_alascan[1])
		alascan = call_alascan(pot_dict, weight_dict, interface_residue_pairs[0], native_score, mode = f_alascan, res_list = reslist)
	else:
		alascan = call_alascan(pot_dict, weight_dict, interface_residue_pairs[0], native_score, mode = f_alascan)

	if f_alascan == None:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1])
	else:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], alscn = alascan)

	return


def predict_binding_all(input_file, pot_dict, z_threshold, output_file, options, f_alascan = None):

	'''Main module to run the software on all the interfaces of a complex'''

	import sys

	pdb_name = output_file.split('/')[-1].replace('.out', '')

	parsed_dist = parse_dist_all(input_file)

	interface_residue_pairs = get_contacts(parsed_dist[0], options['intertype'])

	weight_dict = get_weights(interface_residue_pairs[2], options)

	num_interactors = interface_residue_pairs[3]

	scores = call_score_func(parsed_dist[1], interface_residue_pairs[1], weight_dict, pot_dict)

	native_score = scores[0]
	decoy_scores = scores[1]
	
	out_val = calc_zscore(native_score, decoy_scores)

	linkage_pairs = get_linkages(parsed_dist[0])

	if f_alascan == None:
		pass
	elif len(f_alascan) == 2:
		reslist = get_reslist(f_alascan[1])
		alascan = call_alascan(pot_dict, weight_dict, interface_residue_pairs[0], native_score, mode = f_alascan, res_list = reslist)
	else:
		alascan = call_alascan(pot_dict, weight_dict, interface_residue_pairs[0], native_score, mode = f_alascan)

	if f_alascan == None:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1])
	else:
		output = write_out(output_file, out_val, options, num_interactors, native_score, z_threshold, interface_residue_pairs[1], alscn = alascan)

	return


def main():

	'''Main routine that decides the modules to run based on user input and/or complex type'''

	import os 
	import sys 
	import time
	import cPickle as cPickle

	parent_dir = os.getcwd()
	data_dir = parent_dir + '/data/'
	#input_dir = parent_dir + '/input/'
	#output_dir = parent_dir + '/output/'

	start_time = time.time()

	options = get_options()

	input_file = options['input_pdb']
	if options['outfile'] == None:
		output_file = options['input_pdb'].replace('.pdb', '.out')
	else:
		output_file = (options['outfile'] + options['input_pdb'].split('/')[-1]).replace('.pdb', '.out')

	parsed_pdb = parse_pdb(input_file)
	chain_set = parsed_pdb[1]
	chain_num = parsed_pdb[0]

	if chain_num < 2:
		sys.exit("Status: Aborting Execution \nReason: Input file does not contain a Protein Complex")

	operating_points = pickle.load(open(os.path.join(data_dir, 'operating_points.p'), 'rb'))


	if options['custom_pot'] == None:
		pot = select_pot(options['cutoff'], options['intertype'])
		pot_dict = pickle.load(open(os.path.join(data_dir, pot), 'rb'))

		z_threshold = operating_points[pot.replace('.p', '').replace('_', '.')]
	else:
		pot_file = open(options['custom_pot'], 'r').readlines()
		pot_dict = convert_pot(pot_file)

		z_threshold = - 0.7

	if options['alascan'] == '0':
		if len(chain_set) > 2:
			if options['protein1'] == None and options['protein2'] == None:
				predict_binding_all(input_file, pot_dict, z_threshold, output_file, options)
			else:
				predict_binding_interface(input_file, pot_dict, z_threshold, output_file, options, p1 = options['protein1'], p2 = options['protein2'])
		else:
			predict_binding_interface(input_file, pot_dict, z_threshold, output_file, options, p1 = chain_set[0], p2 = chain_set[1])
	else:
		if len(chain_set) > 2:
			if options['protein1'] == None and options['protein2'] == None:
				predict_binding_all(input_file, pot_dict, z_threshold, output_file, options, f_alascan = options['alascan'])
			else:
				predict_binding_interface(input_file, pot_dict, z_threshold, output_file, options, \
				 p1 = options['protein1'], p2 = options['protein2'], f_alascan = options['alascan'])
		else:
			predict_binding_interface(input_file, pot_dict, z_threshold, output_file, options, \
			 p1 = chain_set[0], p2 = chain_set[1], f_alascan = options['alascan'])

	end_time = time.time()

	return


# Uncomment this block for running the main program

# Run the actual program
#start_time = time.time()
#main()
#end_time = time.time()

#print "Run Completed! Check the 'output' folder for results!"
#print "Time elapsed: ", round(end_time - start_time, 3), 's'



