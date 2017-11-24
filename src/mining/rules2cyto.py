#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to generate pairwise subrules from association rules.
"""

import sys
from collections import defaultdict
from pathlib import Path

try:
    input_file = Path(sys.argv[1])
except IndexError:
    print('No input file was defined.')
try:
    output_file = Path(sys.argv[2])
except IndexError:
    print('No output file was defined.')

with input_file.open() as f:
    rule_dict = defaultdict(int)
    # freq_dict = defaultdict(float)
    conf_dict = defaultdict(float)
    max_conf_dict = defaultdict(float)
    lift_dict = defaultdict(float)
    max_lift_dict = defaultdict(float)
    for line in f:
        # h@GO0000502,h@GO0043230>v@GO0050794;0.5159235668789809;0.17136733359396894
        rule = line.split(';')
        terms = rule[0].split('>')
        antecedents = terms[0].split(',')
        consequents = terms[1].split(',')
        pairwise_rules = [i + '>' + j for i in antecedents for j in consequents]
        for i in pairwise_rules:
            rule_dict[i] += 1
            confidence = float(rule[1])
            conf_dict[i] += confidence
            max_conf_dict[i] = max(confidence, max_conf_dict[i])
            lift = float(rule[2])
            lift_dict[i] += lift
            max_lift_dict[i] = max(lift, max_lift_dict[i])
            # freq_dict[i] += int(rule[3].strip('))\n'))

    # Average confidence and lift
    for key, item in conf_dict.items():
        conf_dict[key] = item / float(rule_dict[key])
    for key, item in lift_dict.items():
        lift_dict[key] = item / float(rule_dict[key])

with output_file.open('w') as o:
    o.write('antecedent\tconsequent\tcount\tmean_confidence\tmax_confidence\tmean_lift\tmax_lift\n')
    for key, item in rule_dict.items():
        ant, con = key.split('>')
        o.write(ant + '\t' + con + '\t' + str(item) + '\t' + str(round(conf_dict.get(key), 4)) + '\t' +
                str(round(max_conf_dict.get(key), 4)) + '\t' + str(round(lift_dict.get(key), 4)) + '\t' +
                str(round(max_lift_dict.get(key), 4)) + '\n')
