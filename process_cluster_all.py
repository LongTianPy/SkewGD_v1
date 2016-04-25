from Bio import SeqIO
import os

def process_cluster(mcl_out,protein_cds,output_prefix,working_dir):
    myfile=open(mcl_out,'r') #all.TAIR10protein.cluster.txt, 7496 clusters
    cluster={}
    line_count=1;
    for line in myfile:
        line=line.rstrip('\n')
        listsq=line.split('\t')
        if len(listsq)>1:
            cluster[line_count]=line
            line_count=line_count+1
    print line_count
    myfile.close()
    proteins_dict = SeqIO.index(protein_cds,"fasta")
    # protein={}
    # for line in myfile:
    #     if '>' in line:
    #         header=line.replace('>','')
    #         header=header.rstrip('\n')
    #         protein[header]=''
    #     else:
    #         protein[header]+=line.strip()
    #for p in protein:
        #print p, protein[p]
    #From each cluster and ids in cluster, print out files
    #num=random.sample(range(1,7496),750)
    #with open('random_list', 'wb') as f:
    #    pickle.dump(num, f)
    current_dir = os.getcwd()
    try:
        os.chdir(working_dir)
    except:
        os.mkdir(working_dir)
        os.chdir(working_dir)
    for i in cluster:
        # output=open(output_prefix+str(i)+'.txt','w') #f = open("file_"+str(i)+".dat","w")
        ids=cluster[i].split("\t")
        cluster_records = [proteins_dict[j] for j in ids]
        with open(output_prefix+str(i)+'.txt','w') as output:
            SeqIO.write(cluster_records,output,"fasta")
    proteins_dict.close()
    os.chdir(current_dir)



