"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine calculates the z-score and other metrics from the score of the native structure and the scores of the decoy structures.

INPUT - Score of the native structure & Scores of the randomized score_decoys
OUTPUT - Z-score, Avg. Background score, 
'''

"""
MODIFICATIONS - added calc_zscore_norm that calculates the zscore from normalized raw_scores
"""


# Calculates the Z-score and other reports other metrics
def calc_zscore(score_native, score_decoys, output_file):

	from math import sqrt
	from err_warn_report import err_report
	import sys

	outvals = {}

	for interface in score_decoys.keys():
		total = 0.000
		false_pos = 0
		num_decoys = len(score_decoys[interface].values())
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

		false_pos_rate = false_pos / num_decoys
		avg_score = total / num_decoys

		std_dev = 0
		for score_el in score_decoys[interface].values():
			std_dev += (score_el - avg_score) ** 2
	
		std_dev = std_dev / num_decoys

		if std_dev == 0:
			error_code = 1
			err_report(error_code, output_file)
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
			z_prime = 'undef'
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


# Calculates the Z-score from the normalized raw_scores
def calc_zscore_norm(score_native_norm, score_decoys_norm, output_file):

	from math import sqrt
	from err_warn_report import err_report
	import sys

	outvals_norm = {}

	for interface in score_decoys_norm.keys():
		total_norm = 0.000
		false_pos_norm = 0
		num_decoys_norm = len(score_decoys_norm[interface].values())
		for score_el in score_decoys_norm[interface].values():
			total_norm += score_el

			if score_el <= score_native_norm[interface]:
				false_pos_norm += 1

			if score_el > score_native_norm[interface]:
				try:
					min_tn_norm
				except NameError:
					min_tn_norm = score_el
				else:
					if score_el <= min_tn_norm:
						min_tn_norm = score_el

			try:
				min_score_norm
			except NameError:
				min_score_norm = score_el
			else:
				if score_el < min_score_norm:
					min_score_norm = score_el

		false_pos_rate_norm = false_pos_norm / num_decoys_norm
		avg_score_norm = total_norm / num_decoys_norm

		std_dev_norm = 0 
		for score_el in score_decoys_norm[interface].values():
			std_dev_norm += (score_el - avg_score_norm) ** 2

		std_dev_norm = std_dev_norm / num_decoys_norm

		if std_dev_norm == 0:
			error_code = 1
			err_report(error_code, output_file)
			sys.exit()

		std_dev_norm = sqrt(std_dev_norm)

		z_bg_norm = []
		for score_el in score_decoys_norm[interface].values():
			z_scr_norm = (score_el - avg_score_norm) / std_dev_norm
			z_bg_norm.append(z_scr_norm)

		z_score_norm = (score_native_norm[interface] - avg_score_norm) / std_dev_norm

		try:
			min_tn_norm
		except NameError:
			z_min_tn_norm = 'undef'
			z_primer_norm = 'undef'
		else:
			z_min_tn_norm = (min_tn_norm - avg_score_norm) / std_dev_norm
			z_prime_norm = z_score_norm - z_min_tn_norm

		z_min_norm = (min_score_norm - avg_score_norm) / std_dev_norm
		z_2 = z_score_norm - z_min_norm

		z_score_norm = round_val(z_score_norm)
		avg_score_norm = round_val(avg_score_norm)
		std_dev_norm = round_val(std_dev_norm)

		outvals_norm.update({interface : (z_score_norm, avg_score_norm, std_dev_norm)})	

	return outvals_norm

# Rounds off the given value to 3 decimal digits.
def round_val(value):

	from math import ceil

	num = value
	num = ceil(num * 1000) / 1000

	return num
