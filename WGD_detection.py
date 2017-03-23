#!/usr/bin/python

# This script implements everything into one and multiprocess tasks that can be paralleled.

# IMPORT
# BUILT-IN MODULES
import argparse
import os
from os import listdir
from os.path import isfile, join
import multiprocessing as mp
import pandas as pd
from datetime import datetime
from distutils.spawn import find_executable
import sys
from functools import partial
import shutil
# SCRIPTS
import prot_to_cds
import run_paml_yn00
import ks_correction
import convert1
import process_blast
import process_cluster_all
import run_muscle
import config

# FUNCTIONS

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
    out_prefix = args.output_pref
    if find_executable(config.YN00_DEFAULT) is not None:
        yn00_binary = config.YN00_DEFAULT
    else:
        if find_executable(args.yn00_path) is not None:
            yn00_binary = args.yn00_path
        else:
            print("Cannot find yn00 executable, please check if it has been installed.")
    if find_executable(config.BLASTP_DEFAULT) is not None:
        blastp_exe = config.BLASTP_DEFAULT
    else:
        if find_executable(args.blastp) is not None:
            blastp_exe = args.blastp
        else:
            print("Cannot find blastp executable, please check if it has been installed.")
    if find_executable(config.MAKEBLASTDB_DEFAULT) is not None:
        makeblastdb_exe = config.MAKEBLASTDB_DEFAULT
    else:
        if find_executable(args.makeblastdb) is not None:
            makeblastdb_exe = args.makeblastdb
        else:
            print("Cannot find makeblastdb executable, please check if it has been installed.")
    if find_executable(config.MUSCLE_DEFAULT) is not None:
        muscle_exe = config.MUSCLE_DEFAULT
    else:
        if find_executable(args.muscle) is not None:
            muscle_exe = args.muscle
        else:
            print("Cannot find MUSCLE executable, please check if it has been installed.")
    if find_executable(config.MCL_DEFAULT) is not None:
        mcl_exe = config.MCL_DEFAULT
    else:
        if find_executable(args.mcl) is not None:
            mcl_exe = args.mcl
        else:
            print("Cannot find MCL executable, please check if it has been installed.")
    identity = args.identity
    coverage = args.coverage
    blastp_threads = args.blastp_threads
    mcl_threads = args.mcl_threads
    mcl_inflation = args.mcl_inflation
    cluster_aln_threads=args.cluster_aln_threads
    print "=================================================="
    print "Welcome to SkewGD pipeline"
    print "Current time:", datetime.now()
    print "==================================================\n"
    if args.nucleotide_cds and not args.cds_folder:
        nucleotide_cds = args.nucleotide_cds
        pipeline_single_cds(nucleotide_cds=nucleotide_cds,output_prefix=out_prefix,identity=identity,coverage=coverage,
                            working_dir=working_dir, blastp_threads=blastp_threads,mcl_threads=mcl_threads,
                            mcl_inflation=mcl_inflation,cluster_aln_threads=cluster_aln_threads,yn00_path=yn00_binary,
                            blastp_exe=blastp_exe,makeblastdb_exe=makeblastdb_exe,muscle_exe=muscle_exe,mcl_exe=mcl_exe)
        clean(working_dir=working_dir,prefix=out_prefix)
    elif args.cds_folder and not args.nucleotide_cds:
        cds_folder = args.cds_folder
        nucleotide_cds_list = [join(cds_folder, f) for f in listdir(cds_folder) if isfile(join(cds_folder,f))]
        clean(working_dir=working_dir, prefix=out_prefix)
    elif args.nucleotide_cds and args.cds_folder:
        print "-i and -I cannot be used at the same time.\nExiting..."
        sys.exit()
    elif not args.cds_folder and not args.nucleotide_cds:
        print "Either one CDS file or a directory with several CDS files should be provided with -i or -I, respectively."
        print "Exiting..."
        sys.exit()



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
    parser.add_argument("-i", dest='nucleotide_cds', help="Full coding sequences of the organism of interest.")
    parser.add_argument("-I", dest="cds_folder",help="A directory with CDS files of different organisms only. NOTE: This"
                                                     " option cannot be used with -i at the same time. Options for threads"
                                                     "need to be set to reasonable number since a maximum of 2 files can"
                                                     "be running at the same time.")
    parser.add_argument("-o", dest='output_pref', help="Prefix for the MCL clustered files.")
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ .")
    parser.add_argument("--blastp", dest="blastp", help="File path to blastp executable. "
                                                                                   "Default: /usr/bin/blastp .")
    parser.add_argument("--makeblastdb", dest="makeblastdb",help="File path to makeblastdb executable. "
                                                                 "Default: /usr/bin/makeblastdb .")
    parser.add_argument("--muscle",dest="muscle",help="File path to MUSCLE executable.")
    parser.add_argument("--mcl",dest="mcl",help="File path to MCL executable. Default: /usr/bin/mcl .")
    parser.add_argument("--yn00", dest="yn00_path", help="File path to yn00 executable. Default: /usr/bin/yn00 .")
    parser.add_argument("--identity", dest="identity", type=int, default=50, help="Threshold of percentage identity in "
                                                                                  "BLAST result. Default: 50 .")
    parser.add_argument("--coverage", dest="coverage", type=int, default=30, help="Threshold of percentage alignment "
                                                                                  "coverage in BLAST result. Default: 30 .")
    parser.add_argument("--blastp_threads", dest="blastp_threads", type=int, default=8, help="Number of threads for "
                                                                                             "running BLASTp. Default: 8 .")
    parser.add_argument("--mcl_threads", dest="mcl_threads", type=int, default=1,help="Number of threads for running "
                                                                                      "MCL. Default: 1 .")
    parser.add_argument("--mcl_inflation", dest="mcl_inflation", type=float, default=2.0, help="Tune the granularity of "
                                                                                               "clustering. Usually "
                                                                                               "choose from "
                                                                                               "the range of [1.2, 5.0]."
                                                                                               " 5.0 makes it finely "
                                                                                               "grained and 1.2 "
                                                                                               "makes clustering "
                                                                                               "coarsed. Default: 2.0 .")
    parser.add_argument("--cluster_aln_threads", dest="cluster_aln_threads",type=int,default=8,help="Number of threads "
                                                                                                    "for parallelling "
                                                                                                    "the alignment of "
                                                                                                    "clusters. Default: "
                                                                                                    "8 .")

    args = parser.parse_args()
    return args

# Clean working directory
def clean(working_dir,prefix):
    os.mkdir(join(working_dir+"output"))
    files = [join(working_dir,file) for file in listdir(working_dir) if isfile(join(working_dir,file))]
    for file in files:
        if file.startswith(prefix):
            shutil.move(file,join(working_dir,"output/"))



# Individual wrappers
def Hong_wrapper(nucleotide_cds,output_prefix,identity,coverage,working_dir,blastp_threads,mcl_threads,mcl_inflation,
                 cluster_aln_threads,blastp_exe,makeblastdb_exe,muscle_exe,mcl_exe):
    protein_cds = nucleotide_cds+".protein"
    print "Step 1 of 9: Translating CDS to protein sequences...", datetime.now()
    convert1.convert(nucleotide_cds)
    process_blast.run_blast(protein_cds=protein_cds,blastp_threads=blastp_threads,blastp_exe=blastp_exe,makeblastdb_exe=makeblastdb_exe)
    mcl_out = protein_cds+".mcl_out"
    process_blast.process_blast_out(protein_cds=protein_cds,identity=identity,coverage=coverage,mcl_threads=mcl_threads,
                                    mcl_inflation=mcl_inflation,mcl_exe=mcl_exe)
    process_cluster_all.process_cluster(mcl_out=mcl_out, protein_cds=protein_cds, output_prefix=output_prefix,working_dir=working_dir)
    cluster_file_list = [join(working_dir,f) for f in listdir(working_dir) if isfile(join(working_dir,f)) and f.endswith(".txt")]
    print "Step 6 of 10: Aligning protein sequences within each cluster...", datetime.now()
    pool_size = cluster_aln_threads
    pool = mp.Pool(processes=pool_size)
    partial_muscle = partial(run_muscle.muscle,muscle_exe=muscle_exe)
    pool.map(partial_muscle, cluster_file_list)
    afa_file_list = [cluster_file+'.afa' for cluster_file in cluster_file_list]
    return afa_file_list

def Andrew_wrapper(prot_cluster_file, nucleotide_file, yn00_binary):
    """
    This function contains Andrew's part
    :param prot_cluster_file: Sequences of clusters
    :param nucleotide_file: cds file
    :return:
    """
    prot_to_cds_out = prot_cluster_file+".phy"
    prot_to_cds.write_align(prot_align_file=prot_cluster_file, nuc_fasta_file=nucleotide_file, nuc_align_file=prot_to_cds_out)
    prot_to_cds_out_sub = prot_to_cds_out+"_sub" # Subtitute dot to 2 spaces
    pattern = 's/\./  /g'
    cmd = "sed '{0}' {1} > {2}".format(pattern, prot_to_cds_out, prot_to_cds_out_sub)
    os.system(cmd)
    run = run_paml_yn00.run_yn00(prot_to_cds_out_sub,yn00_binary)
    return run

def pipeline_single_cds(nucleotide_cds,output_prefix,identity,coverage,working_dir,blastp_threads,mcl_threads,mcl_inflation,
                 cluster_aln_threads,yn00_path,blastp_exe,makeblastdb_exe,muscle_exe,mcl_exe):
    afa_file_list = Hong_wrapper(nucleotide_cds=nucleotide_cds, identity=identity, coverage=coverage,
                                 output_prefix=output_prefix, working_dir=working_dir, blastp_threads=blastp_threads,
                                 mcl_threads=mcl_threads, mcl_inflation=mcl_inflation,cluster_aln_threads=cluster_aln_threads,
                                 blastp_exe=blastp_exe,makeblastdb_exe=makeblastdb_exe,muscle_exe=muscle_exe,mcl_exe=mcl_exe)
    kS_df_total = pd.DataFrame()
    nucleotide_cds_trunc = nucleotide_cds + '_trunc'
    print "Step: 7, 8 and 9 of 10: Reverse translating proteins based on provided CDS, running YN00 to calculate kS " \
          "and correcting kS...", datetime.now()
    # cluster_file_list = [join(working_dir,f) for f in listdir(working_dir) if isfile(join(working_dir,f)) and f.endswith(".txt")]
    # afa_file_list = [cluster_file+'.afa' for cluster_file in cluster_file_list]
    for afa_file in afa_file_list:
        try:
            run = Andrew_wrapper(afa_file, nucleotide_cds_trunc, yn00_path)
            ks_df = ks_correction.correct_ks(run)
            kS_df_total = kS_df_total.append(ks_df)
        except:
            print "Error in " + afa_file + ", skipped."
            continue
    # Draw histogram
    print "Step 10 out of 10: Generating results and plotting on canvas...", datetime.now()
    kS_df_total.to_csv(working_dir + "ks.csv")
    ks_correction.draw_histo(kS_df_total, working_dir)
    print "=================================================="
    print "Thank you for using SkewGD, your analysis has been completed."
    print "Current time:", datetime.now()
    print "==================================================\n"
    return kS_df_total

# def pipeline_folder_cds(cds_folder,output_prefix,identity,coverage,working_dir,blastp_threads,mcl_threads,mcl_inflation,
#                  cluster_aln_threads,yn00_path):
#     pass



if __name__ == "__main__":
    main()