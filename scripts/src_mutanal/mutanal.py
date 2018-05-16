"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""


def mutanal_test(pot_dict, weight_dict, respairs_scores, interface_res, res_list, native_score):

	from collections import OrderedDict

	amino_acids = ['GLY', 'PRO', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'CYS', 'PHE', 'TYR', \
 				   'TRP', 'HIS', 'LYS', 'ARG', 'GLN', 'ASN', 'GLU', 'ASP', 'SER', 'THR']

 	score_diff_dict = OrderedDict()

 	for mutares in res_list.keys():

 		mutares_split = mutares.split(':')

 		mutares_id = mutares_split[0] + ':' + mutares_split[2]

 		mutant_list = res_list[mutares]

 		affctd_interfaces = []

 		for interface in respairs_scores:
 			if mutares_id[-1] in interface:
 				affctd_interfaces.append(interface)

 		for mutant_el in mutant_list:

 			t_score_diff = {}

 			score_diff = 0.0

 			mutate_param = mutate_parameters(mutares_id, affctd_interfaces, respairs_scores, interface_res)
 			decoy_score_dict = mutate_param[0]
 			mutant_dict = mutate_param[1]

 			f_dc_score_dict = mutate_scores(decoy_score_dict, mutant_dict, mutant_el, mutares_id, pot_dict, weight_dict)

 			for interface in f_dc_score_dict:
 				score_diff = f_dc_score_dict[interface] - native_score[interface]
 				t_score_diff.update({interface : score_diff})

 			cum_score_diff = sum(t_score_diff.values())
 			score_diff_dict.update({(mutares, mutant_el) : [cum_score_diff, t_score_diff]})

 				
 		#for element in score_diff_dict:
 		#	if element[0].split(':')[1] != element[1]:
 		#		print element[0], element[1], score_diff_dict[element][0]

	return score_diff_dict


def mutate_scores(decoy_score_dict, mutant_dict, mutate_to, mutares_id, pot_dict, weight_dict):

	for interface in mutant_dict:
		for element in mutant_dict[interface]:
			pair_el = element[0]
			partner_resna = element[1]
			m_weight = weight_dict[interface][pair_el]
			if pair_el.index(mutares_id) == 0:
				m_pair = mutate_to + '-' + partner_resna
			else:
				m_pair = partner_resna + '-' + mutate_to

			if 'GLY' not in m_pair:
				m_score = pot_dict[m_pair] * m_weight
				decoy_score_dict[interface] += m_score
			else:
				pass

	return decoy_score_dict


def mutate_parameters(mutares_id, affctd_interfaces, respairs_scores, interface_res):

	decoy_score_dict = {}
	mutant_dict = {}

	for interface in affctd_interfaces:
		decoy_score_dict.update({interface : 0.0 })
		for pair_el in respairs_scores[interface]:
			split_p_el = pair_el.split('-')
			chain_sort = sorted([split_p_el[0][-1], split_p_el[1][-1]])
			sp_interface = chain_sort[0] + '-' + chain_sort[1]
			if mutares_id not in pair_el:
				decoy_score_dict[sp_interface] += respairs_scores[sp_interface][pair_el]
			elif pair_el.split('-')[0] == mutares_id or pair_el.split('-')[1] == mutares_id:
				mutant_dict.update({sp_interface : []})
				partner_res = pair_el.replace(mutares_id, '').replace('-', '')
				partner_resna = interface_res[partner_res[-1]][partner_res]
				mutant_dict[sp_interface].append([pair_el, partner_resna])

	return decoy_score_dict, mutant_dict
