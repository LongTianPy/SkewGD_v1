#!/usr/bin/env python

import sys
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Data.CodonTable import ambiguous_generic_by_id

if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.5"
    sys.exit(0)

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

def check_trans(identifier, nuc, prot, table):
    """Returns nucleotide sequence if works (can remove trailing stop)"""
    if len(nuc) % 3:
        print("Nucleotide sequence for %s is length %i (not a multiple of three)"
                 % (identifier, len(nuc)))

    p = str(prot).upper().replace("*", "X")
    t = str(nuc.translate(table)).upper().replace("*", "X")
    if len(t) == len(p) + 1:
        if str(nuc)[-3:].upper() in ambiguous_generic_by_id[table].stop_codons:

            t = t[:-1]
            nuc  = nuc[:-3]
    if len(t) != len(p):
        err = ("Inconsistent lengths for %s, ungapped protein %i, "
               "tripled %i vs ungapped nucleotide %i." %
               (identifier, len(p), len(p) * 3, len(nuc)))
        if t.endswith(p):
            err += "\nThere are %i extra nucleotides at the start." % (len(t) - len(p))
        elif t.startswith(p):
            err += "\nThere are %i extra nucleotides at the end." % (len(t) - len(p))
        elif p in t:

            err += "\nHowever, protein sequence found within translated nucleotides."
        elif p[1:] in t:
            err += "\nHowever, ignoring first amino acid, protein sequence found within translated nucleotides."
        print(err)


    if t == p:
        return nuc
    elif p.startswith("M") and "M" + t[1:] == p:
        
        if str(nuc[0:3]).upper() in ambiguous_generic_by_id[table].start_codons:
            return nuc
        else:
            print("Translation check failed for %s\n"
                     "Would match if %s was a start codon (check correct table used)\n"
                     % (identifier, nuc[0:3].upper()))
    else:
        
        m = "".join("." if x==y else "!" for (x,y) in zip(p,t))
        if len(prot) < 70:
            sys.stderr.write("Protein:     %s\n" % p)
            sys.stderr.write("             %s\n" % m)
            sys.stderr.write("Translation: %s\n" % t)
        else:
            for offset in range(0, len(p), 60):
                sys.stderr.write("Protein:     %s\n" % p[offset:offset+60])
                sys.stderr.write("             %s\n" % m[offset:offset+60])
                sys.stderr.write("Translation: %s\n\n" % t[offset:offset+60])
        print("Translation check failed for %s\n" % identifier)

def sequence_back_translate(aligned_protein_record, unaligned_nucleotide_record, gap, table=0):
    
    if not gap or len(gap) != 1:
        raise ValueError("Please supply a single gap character")

    alpha = unaligned_nucleotide_record.seq.alphabet
    if hasattr(alpha, "gap_char"):
        gap_codon = alpha.gap_char * 3
        assert len(gap_codon) == 3
    else:
        from Bio.Alphabet import Gapped
        alpha = Gapped(alpha, gap)
        gap_codon = gap*3

    ungapped_protein = aligned_protein_record.seq.ungap(gap)
    ungapped_nucleotide = unaligned_nucleotide_record.seq
    if table:
        ungapped_nucleotide = check_trans(aligned_protein_record.id, ungapped_nucleotide, ungapped_protein, table)
    elif len(ungapped_protein) * 3 != len(ungapped_nucleotide):
        print("Inconsistent lengths for %s, ungapped protein %i, "
                 "tripled %i vs ungapped nucleotide %i" %
                 (aligned_protein_record.id,
                  len(ungapped_protein),
                  len(ungapped_protein) * 3,
                  len(ungapped_nucleotide)))

    seq = []
    nuc = str(ungapped_nucleotide)
    for amino_acid in aligned_protein_record.seq:
        if amino_acid == gap:
            seq.append(gap_codon)
        else:
            seq.append(nuc[:3])
            nuc = nuc[3:]
    assert not nuc, "Nucleotide sequence for %r longer than protein %r" \
        % (unaligned_nucleotide_record.id, aligned_protein_record.id)

    aligned_nuc = unaligned_nucleotide_record[:] 
    aligned_nuc.letter_annotation = {} 
    aligned_nuc.seq = Seq("".join(seq), generic_dna)
    aligned_nuc.id = aligned_protein_record.id
    assert len(aligned_protein_record.seq) * 3 == len(aligned_nuc)
    return aligned_nuc

def back_translate(protein_alignment, nucleotide_dict):
    back_translated_records = []
    for protein in protein_alignment:
        print str(protein.seq)
        print str(nucleotide_dict[protein.id].seq)
        triplets = re.findall('...', str(nucleotide_dict[protein.id].seq))
        print len(triplets)
        aa_pattern = re.compile(r"[A-Z]")
        aa_index = [m.start(0) for m in re.finditer(aa_pattern, str(protein.seq))]
        print aa_index
        print len(aa_index)
        #nucleotide_seq = [triplets[aa_idx] for aa_idx in aa_index]
        #nucleotide_seq = ''.join(nucleotide_seq)
        #back_translated_records.append(SeqIO.SeqRecord(Seq(nucleotide_seq,alphabet=generic_dna),id=protein.id))
    return back_translated_records



def alignment_back_translate(protein_alignment, nucleotide_records, key_function=None, gap=None, table=0):
    """Thread nucleotide sequences onto a protein alignment."""
    
    if key_function is None:
        key_function = lambda x: x
    if gap is None:
        gap = "-"

    aligned = []
    for protein in protein_alignment:
        try:
            nucleotide = nucleotide_records[key_function(protein.id)]
        except KeyError:
            raise ValueError("Could not find nucleotide sequence for protein %r" \
                             % protein.id)
        aligned.append(sequence_back_translate(protein, nucleotide, gap, table))
    return MultipleSeqAlignment(aligned)


# if len(sys.argv) == 4:
#     align_format, prot_align_file, nuc_fasta_file = sys.argv[1:]
#     nuc_align_file = sys.stdout
#     table = 0
# elif len(sys.argv) == 5:
#     align_format, prot_align_file, nuc_fasta_file, nuc_align_file = sys.argv[1:]
#     table = 0
# elif len(sys.argv) == 6:
#     align_format, prot_align_file, nuc_fasta_file, nuc_align_file, table = sys.argv[1:]
# else:
#     sys_exit("""
#
# Example usage:
#
# $ python prot_to_cds.py fasta prot_cluster_file.fasta nucleotide_file.fasta output_file.fasta
#
# """)
#
# try:
#     table = int(table)
# except:
#     sys_exit("Bad table argument %r" % table)
#
# prot_align = AlignIO.read(prot_align_file, align_format, alphabet=generic_protein)
# nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
# nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-", table=table)
# AlignIO.write(nuc_align, nuc_align_file, align_format)

def write_align(prot_align_file,nuc_fasta_file,nuc_align_file):
    """
    Put the above functions together
    :param prot_align_file: clusters
    :param nuc_fasta_file: whole cds file
    :param nuc_align_file: output
    :return: physical output file
    """
    prot_align = AlignIO.read(prot_align_file,"fasta",alphabet=generic_protein)
    nuc_dict = SeqIO.index(nuc_fasta_file, "fasta")
    #nuc_align = [nuc_dict[i.id] for i in prot_align]
    nuc_align = alignment_back_translate(prot_align, nuc_dict, gap="-", table=0)
    # nuc_align = back_translate(prot_align, nuc_dict)
    # f = open(nuc_align_file,"w")
    AlignIO.write(nuc_align,f,"phylip")
    # f.write("\t{0}\t{1}\n".format(len(nuc_align),len(nuc_align[0].seq)))
    # for i in nuc_align:
    #     f.write(i.id+"  "+i.seq+"\n")
    # f.close()
    nuc_dict.close()

