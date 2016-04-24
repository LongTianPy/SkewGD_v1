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
    yn.out_file = "yn_result.out"
    yn.working_dir = "./"
    yn.print_options()
    run = yn.run(command="/Users/longtian/Desktop/paml4.8/bin/yn00")
    return run

if __name__=='__main__':
    perform_codeml(ctl_file)