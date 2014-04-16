#!/usr/bin/python
"""
The purpose of this script is to produce diffs between the four
different codebases used in the toric-code-decoding project:
 + perfect toric
 + perfect concatenated
 + noisy toric
 + noisy concatenated
I'm expressly leaving out dijkstra.cc and dijkstra.h, because I hope
they'll be rendered unnecessary by a more direct lattice distance 
calculation.
"""
import os
import subprocess as sbp

files_to_diff = ['constants.h', 'lattice.cc', 'lattice.h',
                    'Makefile', 'operator.cc', 'operator.h', 
                    'qecsim.cc', 'qecsim_test', 
                    'README.txt', 'simulation.cc', 'simulation.h']

dirs = ['diffs/pc_pt/','diffs/nc_nt/','diffs/pc_nc/','diffs/pt_nt/']

prefixes = [['perfect_concatenated', 'perfect_toric'],
['noisy_concatenated', 'noisy_toric'],
['perfect_concatenated', 'noisy_concatenated'],
['perfect_toric', 'noisy_toric']]

for idx, diff_dir in enumerate(dirs):
    f_dir_1, f_dir_2 = prefixes[idx]
    f_dir_1 = 'qecsim_' + f_dir_1
    f_dir_2 = 'qecsim_' + f_dir_2
    for filename in files_to_diff:
        diff_filename = diff_dir + '_'.join(filename.split('.')) + '.diff'
        with open(diff_filename, 'w') as diff_file:
            sbp.call(['diff', '-y', '/'.join([f_dir_1, filename]), 
                '/'.join([f_dir_2, filename])], stdout=diff_file)