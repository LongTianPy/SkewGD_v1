#!/usr/bin/python
import glob
import os
import sys
from Bio.Phylo.PAML import yn00

# input_path = sys.argv[1]
# output_file = sys.argv[2]

def run_yn00(input_phy):
    binary = "/Users/longtian/Desktop/paml4.8/bin/yn00"
    yn = yn00.Yn00()
    yn.alignment = input_phy
    yn.out_file = input_phy+".out"
    yn.working_dir = "./"
    yn.set_options(commonf3x4 = 0)
    yn.print_options()
    run_result = yn.run(command=binary)
    return run_result