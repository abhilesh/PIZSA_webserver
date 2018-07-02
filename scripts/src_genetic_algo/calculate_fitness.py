from __future__ import division
import os
import cPickle as pickle
import time

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
		avg_fit = round(sum(element)/len(element), 1)
		avg_Fitness.append(avg_fit)

	return fitness, avg_Fitness

start_time = time.time()
cur_dir = os.getcwd()

og_scores_db = {}
dc_scores_db = {}

for filename in os.listdir(cur_dir):
	if 'og' in filename:
		og_file = pickle.load(open(filename, 'rb'))
		og_scores_db.update({ og_file.keys()[0][-4:] : og_file.values()[0]})
	if 'dc' in filename:
		dc_file = pickle.load(open(filename, 'rb'))
		dc_scores_db.update({ dc_file.keys()[0][-4:] : dc_file.values()[0]})

print 'calc_fitness did run'

Fit_param = calculate_fitness(og_scores_db, dc_scores_db)
pickle.dump(Fit_param[1], open('avg_fitness.p', 'wb'))

print 'calc_fitness_output done'
stop_time = time.time()
print 'Time elappsed:', stop_time - start_time
