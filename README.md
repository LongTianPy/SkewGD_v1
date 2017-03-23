# SkewGD
Summation of Ks pairs for Exploration of Whole Genome Duplications -- Whole Genome Duplication detection pipeline

INTRODUCTION
------------

WGD_detection is a script implementing several bioinformatics software and calculation to visualize gene duplication events by processing the full set of coding sequences (CDS) of an organism.

INPUT: CDS file in FASTA format;

OUTPUT: kS distribution data in csv format and histogram in user-indicated working directory (set by -d).


WORKFLOW
--------

0. CDS translation to protein sequences;
1. Pairwise self-BLASTP by BLASTP;
2. Extraction of pairs of comparison by identity (default: 50%) and coverage (default: 30%);
3. Markov chain clustering by MCL;
4. Sequence Alignment of each cluster by MUSCLE;
5. For each alignment, reverse translation of protein sequences back to nucleotide sequences according to the input CDS;
6. Maximum likelihood phylogenetic analysis on each nucleotide sequence alignment by yn00 from PAML;
7. kS correction and gene duplication event clustering;
8. Data visualization


DEPENDENCIES AND REQUIREMENTS
-----------------------------

WGD_detection is developed in Python 2.x with modules and external software, and is Python 3 compatible.

While running this pipeline, a dependency check is at first performed to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below.

* [**Python 2.x**](https://www.python.org/)
*   Modules can be installed using [pip](https://pip.pypa.io/en/stable/installing/) `pip install [module_name]`
*   **Pandas** v0.16.2
*   **BioPython** v1.64
*   **Seaborn** v0.7.0
* **BLAST for LINUX** ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+ v2.2.26 -- [Installation guide](http://www.ncbi.nlm.nih.gov/books/NBK52640/) (Ubuntu users can directly install by ```sudo apt-get install ncbi-blast+```
* [**MCL**](http://micans.org/mcl/) v14-137
* [**MUSCLE**](http://www.drive5.com/muscle/) v3.8.31
* [**YN00**](http://abacus.gene.ucl.ac.uk/software/paml.html#download) v4.8 Ubuntu users can directly install by ```sudo apt-get install paml```


USAGE
-----

```
usage: WGD_detection.py [-h] [-i NUCLEOTIDE_CDS] [-I CDS_FOLDER]
                        [-o OUTPUT_PREF] [-d WORKING_DIR] 
                        [--blastp BLASTP] [--makeblastdb MAKEBLASTDB]
                        [--muscle MUSCLE] [--mcl MCL]
                        [-yn00 YN00_PATH]
                        [--identity IDENTITY] [--coverage COVERAGE]
                        [--blastp_threads BLASTP_THREADS]
                        [--mcl_threads MCL_THREADS]
                        [--mcl_inflation MCL_INFLATION]
                        [--cluster_aln_threads CLUSTER_ALN_THREADS]

Generate kS distrbution histogram to detect Whole Genome Duplication (WGD)
events. Taking the full coding sequences of an organism as input.

optional arguments:
  -h, --help            show this help message and exit
  -i NUCLEOTIDE_CDS     Full coding sequences of the organism of interest.
  -I CDS_FOLDER         A directory with CDS files of different organisms
                        only. NOTE: This option cannot be used with -i at the
                        same time. Options for threadsneed to be set to
                        reasonable number since a maximum of 2 files canbe
                        running at the same time.
  -o OUTPUT_PREF        Prefix for the MCL clustered files.
  -d WORKING_DIR        Working directory to store intermediate files of each
                        step. Default: ./ .
  --blastp BLASTP       File path to blastp executable. Default:
                        /usr/bin/blastp .
  --makeblastdb MAKEBLASTDB
                        File path to makeblastdb executable. Default:
                        /usr/bin/makeblastdb .
  --muscle MUSCLE       File path to MUSCLE executable.
  --mcl MCL             File path to MCL executable. Default: /usr/bin/mcl .
  --yn00 YN00_PATH      File path to yn00 executable. Default: /usr/bin/yn00 .
  --identity IDENTITY   Threshold of percentage identity in BLAST result.
                        Default: 50 .
  --coverage COVERAGE   Threshold of percentage alignment coverage in BLAST
                        result. Default: 30 .
  --blastp_threads BLASTP_THREADS
                        Number of threads for running BLASTp. Default: 8 .
  --mcl_threads MCL_THREADS
                        Number of threads for running MCL. Default: 1 .
  --mcl_inflation MCL_INFLATION
                        Tune the granularity of clustering. Usually choose
                        from the range of [1.2, 5.0]. 5.0 makes it finely
                        grained and 1.2 makes clustering coarsed. Default: 2.0
                        .
  --cluster_aln_threads CLUSTER_ALN_THREADS
                        Number of threads for parallelling the alignment of
                        clusters. Default: 8 .
```
