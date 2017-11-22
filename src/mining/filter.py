#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import re
import sys

try:
    input_file = Path(sys.argv[1])
except IndexError:
    print('No input file was defined.')
try:
    output_file = Path(sys.argv[2])
except IndexError:
    print('No output file was defined.')

# exclusion = {'GO:0005623': 'cell', 'GO:0005488': 'binding', 'GO:0043226': 'organelle',
#                  'GO:0044422': 'organelle part', 'GO:0044464': 'cell part', 'GO:0033643': 'host cell part',
#                  'GO:0033646': 'host intracellular part', 'GO:0043656': 'intracellular region of host',
#                  'GO:0043657': 'host cell', 'GO:0018995': 'host', 'GO:0044424': 'intracellular part',
#                  'GO:0016032': 'viral process', 'GO:0044215': 'other organism',
#                  'GO:0050789': 'regulation of biological process', 'GO:0005515': 'protein binding',
#                  'GO:0019012': 'virion', 'GO:0044423': 'virion part'}

# reduced filter
# exclusion = {'GO:0005623': 'cell', 'GO:0005488': 'binding', 
#              'GO:0044464': 'cell part', 'GO:0033643': 'host cell part',
#              'GO:0043657': 'host cell', 'GO:0018995': 'host', 
#              'GO:0016032': 'viral process', 'GO:0044215': 'other organism',
#              'GO:0050789': 'regulation of biological process', 'GO:0005515': 'protein binding',
#              'GO:0019012': 'virion', 'GO:0044423': 'virion part'}

# biological process
exclusion = {'GO:0005623': 'cell', 'GO:0005488': 'binding', 
             'GO:0044464': 'cell part', 'GO:0033643': 'host cell part',
             'GO:0043657': 'host cell', 'GO:0018995': 'host', 
             'GO:0016032': 'viral process', 'GO:0044215': 'other organism',
             'GO:0005515': 'protein binding',
             'GO:0019012': 'virion', 'GO:0044423': 'virion part'}

regex = re.compile(r'')

with input_file.open() as f:
    with output_file.open('w') as o:
        for line in f:
            print(line)

            s = line.split(', ')
            terms_s = s[0].split('>')
            ant_s = terms_s[0].strip('(\'').split(',')
            con_s = terms_s[1].strip('\'').split(',')

            print(ant_s)
            print(con_s)

            for i in ant_s[:]:
                if 'GO' in i:
                    term = 'GO:' + i.split('@GO')[1]
                    if term in exclusion:
                        print('removing',i)
                        ant_s.remove(i)

            if not ant_s:
                print('omitting full rule')
                continue

            for i in con_s[:]:
                if 'GO' in i:
                    term = 'GO:' + i.split('@GO')[1]
                    if term in exclusion:
                        print('removing',i)
                        con_s.remove(i)

            if not con_s:
                print('omitting full rule')
                continue

            new_line = "('" + ','.join(ant_s) + '>' + ','.join(con_s) + '\', ' + ', '.join(s[1:])

            o.write(new_line)



