import os
import sys
import shutil
import subprocess
import time
import cPickle as pickle

work_dir = "/work/meiklejohn/abhilesh/"

num_iters = 100 

iteration = 2

while iteration < num_iters:
	working_calc_fit = True
	working_nxt_gen_pop = True
	print "Iteration: ", iteration
	iter_dir_name = 'iter_' + str(iteration)
	if os.path.exists(iter_dir_name):
		pass
	else:
		os.makedirs(iter_dir_name)
	running_data_files = len(os.listdir('running_data'))
	while running_data_files < 282:
		running_data_files = len(os.listdir('running_data'))
		print running_data_files
		time.sleep(600)
	for filename in os.listdir('running_data'):
		os.rename(work_dir+"running_data/"+filename, work_dir+iter_dir_name+'/'+filename)
	shutil.copy(work_dir+'calculate_fitness.py', work_dir+iter_dir_name+'/'+'calculate_fitness.py')
	time.sleep(3)
	calculate_fitness_cmd = "/work/meiklejohn/abhilesh/"+iter_dir_name+"/calculate_fitness.py"
	subprocess.call(calculate_fitness_cmd)
	while working_calc_fit:
		if os.path.exists(work_dir+iter_dir_name+"/avg_fitness.p"):
			working_calc_fit = False
		else:
			time.sleep(15)
	shutil.copy(work_dir+iter_dir_name+'/avg_fitness.p', work_dir+'avg_fitness.p')
	next_gen_pop_cmd = "python /work/meiklejohn/abhilesh/next_gen_pop.py"
	os.system(next_gen_pop_cmd)
	while working_nxt_gen_pop:
		if os.path.exists(work_dir+"next_gen_pop.py"):
			working_nxt_gen_pop = False
		else:	
			time.sleep(10)
	shutil.copy(work_dir+'next_gen_pop.p', work_dir+iter_dir_name+'/'+'next_gen_pop.p')
	os.rename(work_dir+'next_gen_pop.p', work_dir+'init_pop.p')
	os.system("sbatch brun_ukuxhumana_mc.sh")
	iteration += 1

