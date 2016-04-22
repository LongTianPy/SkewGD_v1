import os
import pickle

def muscle(protein_cluster):
    aligned_out = protein_cluster+'.afa'
    os.system('muscle -in {0} -out {1}'.format(protein_cluster,aligned_out))


