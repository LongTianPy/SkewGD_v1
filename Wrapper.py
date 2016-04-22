#!/usr/bin/python

# This script implements everything into one and multiprocess tasks that can be paralleled.

# IMPORT
# BUILT-IN MODULES
import argparse
import os
from os import listdir
from os.path import isfile, join
import multiprocessing as mp
# SCRIPTS
import prot_to_cds
import run_paml_yn00
import ks_correction
import convert1
import process_blast
import process_cluster_all
import run_muscle

# FUNCTIONS

# Arguments
def get_parsed_args():
    """
    Parse the command line arguments
    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.
    :return:
        args: An argparse.Namespace object containing all parsed the arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate kS distrbution histogram to detect Whole Genome Duplication (WGD) events. "+
                    "Taking the full coding sequences of an organism as input.")
    parser.add_argument("-i", dest=nucleotide_cds, help="Full coding sequences of the organism of interest.")
    parser.add_argument("-o", dest=outout_pref, help="Prefix for the MCL clustered files.")
    parser.add_argument("-d", dest=working_dir, default="./", help="Working directory to store intermediate files of each step. Default: ./ .")
    parser.add_argument("--identity", dest="identity", type=int, default=50, help="Threshold of percentage identity in BLAST result. Default: 50 .")
    parser.add_argument("--coverage", dest="coverage", type=int, default=30, help="Threshold of percentage alignment coverage in BLAST result. Default: 30 .")

    args = parser.parse_args()
    return args

# Individual wrappers
def Hong_wrapper(nucleotide_cds,output_prefix,identity,coverage,working_dir):
    protein_cds = nucleotide_cds+".protein"
    convert1.convert(nucleotide_cds,protein_cds)
    process_blast.run_blast(protein_cds=protein_cds,identiy=identity,coverage=coverage)
    mcl_out = protein_cds+".mcl_out"
    process_cluster_all.process_cluster(mcl_out=mcl_out, protein_cds=protein_cds, output_prefix=output_prefix,working_dir=working_dir)
    cluster_file_list = [join(working_dir,) for f in listdir(working_dir) if isfile(join(working_dir,f)) and f.endswith(".txt")]
    pool_size = 8
    pool = mp.Pool(processes=pool_size)
    pool.map(run_muscle.muscle, cluster_file_list)
    return cluster_file_list

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

# MAIN
def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir
    nucleotide_cds = args.nucleotide_cds
    out_prefix = args.output_pref
    identity = args.identity
    coverage = args.coverage
    cluster_file_list = Hong_wrapper(nucleotide_cds=nucleotide_cds, identity=identity, coverage=coverage, output_prefix=out_prefix,working_dir=working_dir)






#if __name__ == "__main__":







