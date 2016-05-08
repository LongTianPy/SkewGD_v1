import numpy
import re
import os
from datetime import datetime

#Run blast
def run_blast(protein_cds):
    print "Step 2 of 10: Building BLAST database...", datetime.now()
    cmd1 = "makeblastdb -in {0} -dbtype prot".format(protein_cds)
    os.system(cmd1)
    print "Step 3 of 10: Performing self-BLAST...", datetime.now()
    blast_out = protein_cds+'.blast_out'
    cmd2 = "blastp -db {0} -query {1} -num_threads 8 -outfmt '6 qseqid sseqid pident qcovs' > {2}".format(protein_cds,protein_cds,blast_out)
    os.system(cmd2)

#Process blast output

def process_blast_out(protein_cds,identity,coverage):
    blast_out = protein_cds+'.blast_out'
    blast_processed = protein_cds+'.blast_processed'
    mcl_out = protein_cds+'.mcl_out'
    myfile=open(blast_out,'r')
    doc={}
    print "Step 4 of 10: Filtering BLAST result, using identity threshold of {0} and coverage threshold of {1}".format(identity,coverage), datetime.now()
    for line in myfile:
        listsq=line.split("\t")
        if listsq[0]!=listsq[1]:
            key=listsq[0]+'\t'+listsq[1]
            rkey=listsq[1]+'\t'+listsq[0]
            if doc.has_key(key)==0 and doc.has_key(rkey)==0:
                doc[key]=listsq[2]+'\t'+listsq[3].rstrip('\n')
            elif doc.has_key(key)==1:
                doc[key]=doc[key]+'\n'+listsq[2]+'\t'+listsq[3].rstrip('\n')  #rstrip is to remove last space line
    for i in doc:
        value=doc[i].split("\n")
        if len(value)>1: #compare identity pick the highest on
            ident_list=[]
            for j in value:
                ident=j.split('\t')[0]
                ident_list.append(ident)
                index=ident_list.index(max(ident_list))
                new_value=value[index]
                doc[i]=new_value
    iden_thre=identity
    cov_thre=coverage
    output=open(blast_processed,'w')
    for i in doc:
        doc_split=doc[i].split('\t')
        ident=float(doc_split[0])
        cov=int(doc_split[1])
        if ident>=iden_thre and cov>cov_thre:
            output.write(i+'\t'+doc[i]+'\n')
    output.close()
    myfile.close()
    #Clustering using mcl
    #Step 4: Clustering
    print "Step 5 of 10: Clustering proteins with Markov clustering algorithm...", datetime.now()
    os.system("mcl {0} --abc -o {1}".format(blast_processed, mcl_out))
