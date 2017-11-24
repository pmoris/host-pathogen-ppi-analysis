#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

"""

Tool to format lines of the format:
('v@GO0003677,v@GO0005515,v@GO0030683>h@GO0044238', 157, (0.7804054054054054, 231))
h@GO1902578,h@IPR016024>v@GO0033647;0.5474452554744526;0.26993124967785675


to 

h@GO name, h@IPR016024 > v@GO0name; 0.5474452554744526; 0.26993124967785675

('v@"GO name",v@"GO name",v@"GO name">h@"GO name"', 157, (0.7804054054054054, 231))

NOTE: paths to obo-tools are hardcoded.
"""


import os
import re
import sys
sys.path.append(os.path.abspath('..'))

from pathlib import Path
from go_tools import obo_tools

try:
    input_file = Path(sys.argv[1])
except IndexError:
    print('No input file was defined.')
try:
    output_file = Path(sys.argv[2])
except IndexError:
    print('No output file was defined.')

go_dict = obo_tools.importOBO(r'../../data/raw/go_data/go.obo')

def create_interpro_dict(file):
    interpro_dict = {}
    regex = re.compile('(IPR)(\d*)', re.IGNORECASE)
    with Path(file).open('r') as f:
        for line in f:
            l = line.split('\t')
            id = l[1]
            id = regex.sub(r'\1' + ':' + r'\2', id)
            explanation = l[2]
            interpro_dict[id] = explanation
    return interpro_dict

interpro_dict = create_interpro_dict(r'../../data/interim/interpro_data/protein2ipr_filtered.txt')

def grab_name(match):
    # l = element.split('@')
    # print(l)
    # go_name = go_dict.get([1]).name
    # return l[0]+'-'+go_name
    id = match.group(1) + ':' + match.group(2)
    if id in go_dict:
        return 'GO-' + go_dict[id].name
    elif id in interpro_dict:
        return 'IP-' + interpro_dict[id]
    else:
        print('{} was not found in GO dictionary'.format(id))
        return id

with input_file.open() as f:
    with output_file.open('w') as o:
        regex = re.compile('(GO|IPR)(\d*)', re.IGNORECASE)
        regex_hv = re.compile('(h|v)@', re.IGNORECASE)
        for line in f:
            substituted = regex.sub(grab_name, line)
            substituted = regex_hv.sub(lambda x: x.group(1).upper() + '-', substituted)
            o.write(substituted)
