#!/usr/bin/env python

''' Genetic algorithm to optimize the statistical potential used to discriminate between native and non-native structures.

'''

from __future__ import division
from ukuxhumana import get_options, parse_pdb, parse_dist_all, get_contacts, get_weights, round_val
import collections
import time

def init_pop():

	'''Get the initial population of statistical potentials.'''

	import os
	import random
	import cPickle as pickle

	parent_dir = os.getcwd()
	data_dir = parent_dir + '/data/'

	pot_dict = pickle.load(open(os.path.join(data_dir, '4_ss_norm_cifa_avg.p'), 'rb'))

	amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
	'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

	aa_pairs = [a + '-' + b for a in amino_acids for b in amino_acids]

	pot_list = []
	pot_dict_db = []

	for pair in aa_pairs:
		pot_list.append(pot_dict[pair])

	initial_population = []
	initial_population.append(pot_list)

	while len(initial_population) < 150:

		pot_list = []

		for pair in aa_pairs:
			pot_list.append(pot_dict[pair])

		new_indv = pot_list 

		for num in random.sample(range(400), 40):
			alt = random.choice(range(-6,6,1))
			new_indv[num] = round(new_indv[num] + alt, 3)

 		initial_population.append(new_indv)

 	for pot in initial_population:
 		temp_pot_dict = {}
 		for pair in aa_pairs:
 			temp_weight = pot[aa_pairs.index(pair)]
 			temp_pot_dict.update({pair : temp_weight})

 		pot_dict_db.append(temp_pot_dict)

 	pickle.dump(pot_dict_db, open('/work/meiklejohn/abhilesh/init_pop.p', 'wb'))

 	return pot_dict_db 


def aa_pairs():

	amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
	'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

	aa_pairs = [a + '-' + b for a in amino_acids for b in amino_acids]

	return aa_pairs


def score_ga(respairs, weight_dict, pot_dict):

	native_scores = []

	for interface in weight_dict.keys():
		num_pairs = len(weight_dict[interface])
		if num_pairs != 0:
			score_native = 0.0
			for pair_el in weight_dict[interface]:
				resna_pair = respairs[interface][pair_el]
				weight = weight_dict[interface][pair_el]
				score = pot_dict[resna_pair] * weight
				score_native += score

			score_native = round_val(score_native / num_pairs) 

			#print score_native, num_pairs, round_val(score_native / num_pairs)

		else:
			pass

	return score_native


def scores_per_structure(pot_dict_db, input_file):

	'''Uses multiple potentials to score a list of structures'''

	import sys

	options = get_options()

	#input_file = options['input_pdb']

	filename = input_file.split('/')[-1].replace('.pdb', '')
	
	if options['outfile'] == None:
		output_file = options['input_pdb'].replace('.pdb', '.out')
	else:
		output_file = (options['outfile'] + options['input_pdb'].split('/')[-1]).replace('.pdb', '.out')

	parsed_pdb = parse_pdb(input_file)
	chain_set = parsed_pdb[1]
	chain_num = parsed_pdb[0]

	if chain_num < 2:
		sys.exit("Status: Aborting Execution \nReason: Input file does not contain a Protein Complex")

	pot_list = init_pop()

	parsed_dist = parse_dist_all(input_file)

	interface_residue_pairs = get_contacts(parsed_dist[0], options['intertype'])

	weight_dict = get_weights(interface_residue_pairs[2], options)

	num_interactors = interface_residue_pairs[3]

	scores_list = []

	for pot in pot_dict_db:
		score = score_ga(interface_residue_pairs[1], weight_dict, pot)
		scores_list.append(score)

	return filename, scores_list


def score_all_structures(files_list, work_dir, pot_dict_db):

	# Might wanna change this section in accordance with file naming conventions
	# Rename the files as XXXX_og/dc_X.pdb for this module to work correctly 

	og_scores_db = {}
	dc_scores_db = {}

	for filename in files_list:
		if filename[0] != '.':
			print filename
			sc_per_pdb = scores_per_structure(pot_dict_db, work_dir + filename)
			pdb_id = filename.split('_')[0]
			case_id = '_'.join(filename.split('_')[1:]).replace('.pdb', '')
			db_key = (pdb_id, case_id)
			print db_key 
			if 'og' in case_id:
				og_scores_db.update({pdb_id : sc_per_pdb[1]})
			elif 'dc' in case_id:
				if pdb_id not in dc_scores_db:
					dc_scores_db.update({pdb_id : [sc_per_pdb[1]]})
				else:
					dc_scores_db[pdb_id].append(sc_per_pdb[1])
			#if pdb_id not in scores_db.keys():
			#	scores_db.update({pdb_id : {}})
			#scores_db[pdb_id].update({case_id : sc_per_pdb[1]})
			#scores_db.update({db_key : sc_per_pdb[1]})

	return og_scores_db, dc_scores_db


def calculate_fitness(og_scores_db, dc_scores_db):

	fitness = []
	avg_Fitness = []

	for pdb in og_scores_db:
		len_potList = len(og_scores_db[pdb])
		for i in range(len_potList):
			dc_scores = []
			og_score = og_scores_db[pdb][i]
			for decoy in dc_scores_db[pdb]:
				dc_scores.append(decoy[i])
			rank = sorted(dc_scores + [og_score]).index(og_score)
			fit = round(1 - (rank / (len(dc_scores_db[pdb]))), 1)
			if len(fitness) < len_potList:
				fitness.insert(i, [fit])
			else:
				fitness[i].append(fit)

	for element in fitness:
		avg_fit = round(sum(element)/len(element), 2)
		avg_Fitness.append(avg_fit)

	return fitness, avg_Fitness


def next_gen_pop(avg_Fitness, pot_list, len_new_pop, aa_pairs):

	'''Select the parents for the next generation, cross them over, mutate a little and voila we have new potentials'''

	import random

	def RWS(avg_Fitness, len_new_pop):

		'''Roulette Wheel Selection to select the parents for the next generation'''

		max_fitness = max(avg_Fitness)
		index_List = range(len(avg_Fitness))
		sel_parents = []

		while len(sel_parents) < len_new_pop:
			rnd = random.random()
			rnd_choice = random.choice(index_List)
			if rnd < (avg_Fitness[rnd_choice] / max_fitness) and rnd_choice not in sel_parents:
				sel_parents.append(rnd_choice)
	
		return sel_parents

	new_parents = RWS(avg_Fitness, len_new_pop)
	parent_pots = []

	for i in new_parents:
		parent_pots.append(pot_list[i])

	#print parent_pots

	def crossover(parent_pots, aa_pairs):

		gen_F1 = []

		for i in range(int(len(parent_pots) / 2)):
			rnd_pot_1 = random.choice(parent_pots)
			rnd_pot_2 = random.choice(parent_pots)

			crsovr_pt = random.choice(range(len(rnd_pot_1)))

			#aa_sorted = sorted(rnd_pot_1.keys())

			val_1 = []
			val_2 = []
			new_pot_1 = {}
			new_pot_2 = {}

			for aa in aa_pairs:
				val_1.append(rnd_pot_1[aa])
				val_2.append(rnd_pot_2[aa])

			new_list_1 = val_1[:crsovr_pt] + val_2[crsovr_pt:]
			new_list_2 = val_2[:crsovr_pt] + val_1[crsovr_pt:]

			for i in range(len(aa_pairs)): 
				new_pot_1.update({aa_pairs[i] : val_1[i]})
				new_pot_2.update({aa_pairs[i] : val_2[i]})

			gen_F1.append(new_pot_1)
			gen_F1.append(new_pot_2)

		return gen_F1

	def mutate(gen_F1, aa_pairs):

		mut_gen_F1 = []

		for pot in gen_F1:
			for num in random.sample(range(400), 40):
				pair_aa = aa_pairs[num]
				val = pot[pair_aa]
				alt = random.choice(range(-6,6,1))
				pot.update({pair_aa : round(val + alt, 3)})

			mut_gen_F1.append(pot)

		return mut_gen_F1

	new_potentials = mutate(crossover(parent_pots, aa_pairs), aa_pairs)

	return new_potentials


def main_ga():

	import os
	import sys
	import cPickle as pickle

	initial_populations = init_pop()
	pot_dict_db = initial_populations
	amino_acid_pairs = aa_pairs()

	###files_list = os.listdir(sys.argv[1])
	files_list = []
	files_list.append(raw_input()) #file name in string format bash input
    
	work_dir = sys.argv[1]

	running_scores = score_all_structures(files_list, work_dir, pot_dict_db)
	og_scores_db = running_scores[0]
	dc_scores_db = running_scores[1]

	fitness_param = calculate_fitness(og_scores_db, dc_scores_db)
	fitness = fitness_param[0]
	avg_Fitness = fitness_param[1]

    ###technically avg_fitness == fitness since it is calculated for a single pdb
	pickle.dump(fitness, open('/work/meiklejohn/abhilesh/fit_{0}.p'.format(files_list[0]), 'wb'))
    pickle.dump(running_scores, open('/work/meiklejohn/abhilesh/rscores_{0}.p'.format(files_list[0]), 'wb'))
	

	#next_gen = next_gen_pop(avg_Fitness, pot_dict_db, 30, amino_acid_pairs)
    #run next_gen_pop separately after calculating avg_fitness 
	return


start_time = time.time()
main_ga()
end_time = time.time()

print 'Time elapsed: ', end_time - start_time, 's'

