# WGD
Whole Genome Duplication detection pipeline

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

WGD_detection is developed in Python 2.x with modules and external software. Except Python modules, make sure all external software can be executed by directly using software name, i.e. `python`, `blastp`, `mcl` and `muscle`.
For YN00 from PAML, file path of the YN00 binary is required to be set from command line.

* [**Python 2.x**](https://www.python.org/)
*   Modules can be installed using [pip](https://pip.pypa.io/en/stable/installing/) `pip install [module_name]`
*   **Pandas**
*   **BioPython**
*   **Seaborn**
* [**BLAST for LINUX**] ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0 -- [Installation guide](http://www.ncbi.nlm.nih.gov/books/NBK52640/)
* [**MCL**](http://micans.org/mcl/)
* [**MUSCLE**](http://www.drive5.com/muscle/)
* [**YN00**](http://abacus.gene.ucl.ac.uk/software/paml.html#download)


USAGE
-----

```
usage: WGD_detection.py [-h] [-i NUCLEOTIDE_CDS] [-o OUTPUT_PREF]
                        [-d WORKING_DIR] [-y YN00_PATH] [--identity IDENTITY]
                        [--coverage COVERAGE]

Generate kS distrbution histogram to detect Whole Genome Duplication (WGD)
events. Taking the full coding sequences of an organism as input.

  -h, --help           show this help message and exit
  -i NUCLEOTIDE_CDS    Full coding sequences of the organism of interest.
  -o OUTPUT_PREF       Prefix for the MCL clustered files.
  -d WORKING_DIR       Working directory to store intermediate files of each
                       step. Default: ./ .
  -y YN00_PATH         File path to yn00 executable.
  --identity IDENTITY  Threshold of percentage identity in BLAST result.
                       Default: 50 .
  --coverage COVERAGE  Threshold of percentage alignment coverage in BLAST
                       result. Default: 30 .
```
