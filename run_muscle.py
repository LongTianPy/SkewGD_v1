import os

def muscle(protein_cluster,muscle_exe):
    aligned_out = protein_cluster+'.afa'
    os.system('{2} -in {0} -out {1} -quiet'.format(protein_cluster,aligned_out,muscle_exe))


