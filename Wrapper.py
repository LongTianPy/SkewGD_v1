#!/usr/bin/python

# This script implements everything into one and multiprocess tasks that can be paralleled.

# IMPORT
# BUILT-IN MODULES
import optparse
# SCRIPTS
import prot_to_cds
import run_paml_yn00
import ks_correction

def Andrew_wrapper(prot_cluster_file, nucleotide_file):
    """
    This function contains Andrew's part
    :param prot_cluster_file: Sequences of clusters
    :param nucleotide_file: cds file
    :return:
    """
    prot_to_cds_out = prot_cluster_file+".phy"
    prot_to_cds.write_align(prot_align_file=prot_cluster_file, nuc_fasta_file=nucleotide_file, nuc_align_file=prot_to_cds_out)
    run = run_paml_yn00.run_yn00(prot_to_cds_out)
    return run

def Long_wrapper(yn00_obj):
    ks_df = ks_correction.correct_ks(yn00=yn00_obj)
    ks_correction.draw_histo(ks_df=ks_df)







