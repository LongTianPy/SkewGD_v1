import numpy
import re
import os

#Run blast
def run_blast(protein_cds):
    cmd1 = "makeblastdb -in {0} -dbtype prot".format(protein_cds)
    os.system(cmd1)
    blast_out = protein_cds+'.blast_out'
    cmd2 = "blastp -db {0} -query {1} -num_threads 8 -outfmt '6 qseqid sseqid pident qcovs' > {2}".format(protein_cds,protein_cds,blast_out)
    os.system(cmd2)

#Process blast output

def process_blast(protein_cds,identity,coverage):
    blast_out = protein_cds+'.blast_out'
    blast_processed = protein_cds+'.blast_processed'
    mcl_out = protein_cds+'.mcl_out'
    myfile=open(blast_out,'r')
    doc={}
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
    iden_thre=identiy
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
    os.system("mcl {0} --abc -o {1}".format(blast_processed, mcl_out))
