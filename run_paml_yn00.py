#!/usr/bin/python
import glob
import os
import sys
from Bio.Phylo.PAML import yn00

input_path = sys.argv[1]
output_file = sys.argv[2]

file_counter = 0
for input_file in glob.glob(os.path.join(input_path,'*.txt')):
        with open(input_file, 'rU'):
		runyn = yn00.Yn00()
		runyn.alignment = input_file
		file_counter += 1
		output_file_name = 'output_file'+'%d' % file_counter
		runyn.out_file = output_file_name
        	
print "%d files processed" % file_counter