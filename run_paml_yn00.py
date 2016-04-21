#!/usr/bin/python
import glob
import os
import sys
from Bio.Phylo.PAML import yn00

# input_path = sys.argv[1]
# output_file = sys.argv[2]

def run_yn00(input_phy):
    binary = "/Users/longtian/Desktop/paml4.8/bin/yn00"
    os.mkdir("yn00_tmp")
    yn = yn00.Yn00()
    yn.alignment = input_phy
    yn.set_options(seqtype = 1)
    yn.set_options(CodonFreq = 2)
    yn.set_options(model = 0)
    yn.out_file = input_phy+".out"
    yn.working_dir("yn00_tmp")
    run = yn.run(command=binary)
    return run