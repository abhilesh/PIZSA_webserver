import cPickle as pickle

def aa_pairs():

	amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
	'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

	aa_pairs = [a + '-' + b for a in amino_acids for b in amino_acids]

	return aa_pairs


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


init_pop = pickle.load(open('init_pop.p', 'rb'))
avg_Fitness = pickle.load(open('avg_fitness.p', 'rb'))
aar_pairs = aa_pairs()
len_new_pop = 100

new_potentials = next_gen_pop(avg_Fitness, init_pop, len_new_pop, aar_pairs)
pickle.dump(new_potentials, open('next_gen_pop.p', 'wb'))