from Bio import SeqIO
from Bio.Alphabet import generic_dna
#Step 1: Convert cds to protein sequences

def convert(nucleotide_cds):
    protein_out = nucleotide_cds+".protein"
    myfile=open(nucleotide_cds,'r')
    # doc={}
    # for line in myfile:
    #     #line=line.rstrip('\n')
    #     if '>' in line:
    #         header = line
    #         doc[header]={'sequence':''}
    #     else:
    #         doc[header]['sequence']+=line.strip()
        #protein=translator.translate_dna_single(line)
        #slist=slist.append(protein)
        #line=myfile.readline()
        #allseq=''.join(slist)
        #output.write(header+'\n')
        #output.write(allseq+'\n')
    records_trunc = list(SeqIO.parse(myfile,"fasta"))
    myfile.close()
    myfile=open(nucleotide_cds,'r')
    records_prot = list(SeqIO.parse(myfile,"fasta"))
    myfile.close()
    count = 1
    for record in records_trunc:
        record.id = "gene"+str(count)+"."
        count += 1
        record.seq = record.seq[:-3]
    nucleotide_cds_trunc = nucleotide_cds+'_trunc'
    with open(nucleotide_cds_trunc,"w") as f:
        SeqIO.write(records_trunc,f,"fasta")
    # for i in doc:
    #     protein=translator.translate_dna_single(doc[i]['sequence'])
    #     if(protein.startswith('M') and protein.endswith('_')):
    #         output.write(i)
    #         output.write(protein[:-1]+'\n')
    # myfile.close()
    # output.close()
    #for i in doc:
        #print i, doc[i]['sequence']
    count = 1
    for i in records_prot:
        i.id = "gene"+str(count)+"."
        i.seq = i.seq.translate()
        i.seq = i.seq[:-1]
        count += 1
    with open(protein_out,'w') as outfile:
        SeqIO.write(records_prot,outfile,'fasta')
    print "{0} CDS have been tranlated to proteins".format(len(records_trunc))





