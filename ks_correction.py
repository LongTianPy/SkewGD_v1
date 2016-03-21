#!/usr/bin/python
"""
This script take the result of CODEML as input
"""

### IMPORT
import pycodeml
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

### FUNCTIONS
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
    dS_df.index=seq_names
    dS_df_filtered = dS_df[dS_df<5][dS_df>0]
    # There are n sequences/samples, so we are expecting n-1 events
    cluster_ID = ['Cluster'+str(i+1) for i in range(len(seq_names))]
    clusters = {}
    for i in cluster_ID:
        clusters[i]=[]
    for i in seq_names:
        tmp = dS_df_filtered[i].dropna()
        tmp_idx = tmp.index()









    filtered_ks_w_idx = [i for i in total_ks_w_idx if i[1]<5]
    return filtered_ks_w_idx

