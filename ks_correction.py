#!/usr/bin/python
"""
This script take the result of CODEML as input
"""

### IMPORT
import pycodeml
import pandas as pd

### FUNCTIONS
def loadcluster(df,seq_names):
    """
    Load clusters with correponding pairs of samples and ks of them
    :param df: the dataframe with all ks, filtered (0<ks<5)
    :param seq_names: row names (index)
    :return: a dictionary with cluster IDs and their component pairs as well as their ks
    """
    # There are n sequences/samples, so we are expecting n-1 events
    # cluster_ID = ['Cluster'+str(i+1) for i in range(len(seq_names))]
    # clusters = {}
    pairs = []
    # for i in cluster_ID:
    #     clusters[i]=[]
    for i in range(len(seq_names)):
        tmp = df[seq_names[i]].dropna()
        tmp = tmp[tmp.index>i]
        tmp_idx = tmp.index
        for j in tmp_idx:
            pairs.append([seq_names[i],seq_names[j],tmp.loc[j]])
    return pairs


def correct_ks():
    yn00 = pycodeml.perform_codeml()
    seq_names = yn00.keys()
    seq_names.sort()
    dS_df = pd.DataFrame()
    for i in range(len(seq_names)):
        targets = yn00[seq_names[i]].keys()
        targets.sort()
        dS_targets = [yn00[seq_names[i]][j]['NG86']['dS'] for j in targets]
        dS_targets.insert(i,0)
        dS_df[seq_names[i]] = dS_targets
    # dS_df.index=seq_names
    dS_df_filtered = dS_df[dS_df<5][dS_df>0]
    clusters = loadcluster(dS_df_filtered,seq_names)

    k=1 # Count the number of clusters to be added









    average_ks_per_cluster = []
    for i in range(len(clusters.keys())):
        tmp = [i[2] for i in clusters[clusters.keys()[i]]]
        average_ks_per_cluster.append(sum(tmp)/len(tmp))
    # Sort clusters by clade sizes
    cluster_size = [len(clusters[i])+1 for i in clusters.keys()]
    cluster_size_idx = [[j,i] for i,j in enumerate(cluster_size)]
    cluster_size_idx.sort()
    kS={}
    for i in range(len(cluster_size_idx)):
        idx = cluster_size_idx[i][1]
        avg_ks = average_ks_per_cluster[cluster_size_idx[i][1]]
        components = [j[0] for j in clusters[clusters.keys()[idx]]]+[clusters.keys()[idx]]
        kS['cluster'+str(i+1)] = [avg_ks,components]
    return kS




