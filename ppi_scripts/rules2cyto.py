#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
from pathlib import Path
import re
import sys

sys.path.insert(0,
                r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/ppi_scripts/go_tools')
from go_tools import obo_tools

try:
    input_file = Path(sys.argv[1])
except IndexError:
    print('No input file was defined.')
try:
    output_file = Path(sys.argv[2])
except IndexError:
    print('No output file was defined.')

go_dict = obo_tools.importOBO(
    r'/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/host-pathogen-ppi-fim/go_data/go.obo')

# with input_file.open() as f:
#     rule_dict = defaultdict(int)
#     for line in f:
#         rule = line.split(', ')[0].strip('(\'')
#         print(rule)
#         # v@GO0003677,v@GO0005515,v@GO0016032>h@GO0043226,h@GO0044237,h@GO0044238,h@GO0044464
#         split_rule = rule.split('>')
#         # print(split_rule)
#         # consequent = '>' + split_rule[1]
#         consequents = ['>' + i for i in split_rule[1].split(',')]
#         print(consequents)
#         antecedents = split_rule[0].split(',')
#         print(antecedents)
#         # pairwise_rules = [i + consequent for i in antecedents]
#         pairwise_rules = [i + j for i in antecedents for j in consequents]
#         print(pairwise_rules)
#         for i in pairwise_rules:
#             rule_dict[i] += 1

with input_file.open() as f:
    rule_dict = defaultdict(int)
    freq_dict = defaultdict(int)
    conf_dict = defaultdict(int)
    for line in f:
        # ('v@GO0006355,v@GO0042025>h@GO0005488,h@GO0043226', 0, (0.9119718309859155, 1036))
        rule = line.split(', ')
        terms = rule[0].split('>')
        antecedents = terms[0].strip('(\'').split(',')
        consequents = terms[1].strip('\'').split(',')
        pairwise_rules = [i + '>' + j for i in antecedents for j in consequents]
        # print(pairwise_rules)
        for i in pairwise_rules:
            rule_dict[i] += 1
            conf_dict[i] += float(rule[2].strip('('))
            freq_dict[i] += int(rule[3].strip('))\n'))

    for key, item in conf_dict.items():
        conf_dict[key] = item / float(rule_dict[key])

with output_file.open('w') as o:
    o.write('antecedent\tconsequent\tcount\tmean_confidence\tfreq_sum\n')
    for key, item in rule_dict.items():
        ant, con = key.split('>')
        o.write(ant + '\t' + con + '\t' + str(item) + '\t' + str(round(conf_dict.get(key),4)) + '\t' + str(freq_dict.get(key)) + '\n')
