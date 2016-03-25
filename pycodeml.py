#!/usr/bin/python
"""
This script is written for getting the results from codeml.
"""

### IMPORT
from Bio.Phylo.PAML import yn00
import sys

def perform_codeml():
    yn00 = yn00.Yn00()
    yn00.alignment = "TF105351.Eut.3.phy"
    yn00.tree = "TF105351.Eut.3.53876.tree"
    yn00.out_file = "yn_result.out"
    yn00.working_dir = "./"
    yn00.read_ctl_file = "TF105351.Eut.3.53876.ctl"
    run = yn00.run(command="/Users/longtian/Desktop/paml4.8/bin/yn00")
    return run

if __name__=='__main__':
    perform_codeml(ctl_file)