#!/usr/bin/python

# This script implements everything into one and multiprocess tasks that can be paralleled.

import optparse
import prot_to_cds
import run_yn00
import os

def Andrew_wrapper(prot_cluster_file, nucleotide_file, prot_to_cds_out):
    """
    This function contains Andrew's part
    :param prot_cluster_file:
    :param nucleotide_file:
    :param prot_to_cds_out:
    :return:
    """
    cmd = "python prot_to_cds.py fasta {0} {1} {2}".format(prot_cluster_file, nucleotide_file, prot_to_cds_out)
    os.system(cmd)


