#!/usr/bin/python
import glob
import os
import sys
from Bio.Phylo.PAML import yn00

# input_path = sys.argv[1]
# output_file = sys.argv[2]

def run_yn00(input_phy,yn00_binary):
    binary = yn00_binary
    yn = yn00.Yn00()
    yn.alignment = input_phy
    yn.out_file = input_phy+".out"
    yn.working_dir = "./"
    yn.set_options(commonf3x4 = 1)
    print "Analyzing "+input_phy
    run_result = yn.run(command=binary,verbose=False)
    return run_result