#!/usr/bin/python
"""
This script is written for getting the results from codeml.
"""

### IMPORT
from Bio.Phylo.PAML import yn00
import sys

def perform_codeml():
    yn = yn00.Yn00()
    yn.alignment = "TF105351.Eut.3.phy"
    yn.tree = "TF105351.Eut.3.53876.tree"
    yn.out_file = "yn_result.out"
    yn.working_dir = "./"
    yn.read_ctl_file = "TF105351.Eut.3.53876.ctl"
    run = yn.run(command="/Users/longtian/Desktop/paml4.8/bin/yn00")
    return run

if __name__=='__main__':
    perform_codeml(ctl_file)