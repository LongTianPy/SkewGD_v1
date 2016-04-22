import translator
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
    records = list(SeqIO.parse(myfile,"fasta"))
    myfile.close()
    # for i in doc:
    #     protein=translator.translate_dna_single(doc[i]['sequence'])
    #     if(protein.startswith('M') and protein.endswith('_')):
    #         output.write(i)
    #         output.write(protein[:-1]+'\n')
    # myfile.close()
    # output.close()
    #for i in doc:
        #print i, doc[i]['sequence']
    for i in records:
        i.seq = i.seq.translate()
    with open(protein_out,'w') as outfile:
        SeqIO.write(records,outfile,'fasta')
    outfile.close()




