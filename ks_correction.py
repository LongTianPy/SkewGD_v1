#!/usr/bin/python
"""
This script take the result of CODEML as input
"""

### IMPORT
import pandas as pd
from operator import itemgetter
import sys
try:
    import seaborn
except:
    print("Please install Python module seaborn.")
    sys.exit()
import matplotlib.pyplot as plt
from scipy import stats
from numpy import log2

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
    pairs = sorted(pairs,key=itemgetter(2))
    index = []
    pairs_dict = {}
    for i in pairs:
        sort_key = i[:2]
        sort_key.sort()
        sort_key='-'.join(sort_key)
        index.append(sort_key)
        pairs_dict[sort_key]=i[2]
    return index, pairs_dict

def cluster_finder(index, pairs_dict):
    tax_cluster = {}
    cluster_tax = {}
    cluster_ks = {}

    def get_biggest_cluster(node):
        cluster_list = tax_cluster[node]
        cluster_size = [len(cluster_tax[each]) for each in cluster_list]
        cluster_size_w_idx = [[size,idx] for idx,size in enumerate(cluster_size)]
        cluster_size_w_idx.sort()
        idx_biggest_cluster = cluster_size_w_idx[-1][1]
        biggest_cluster = cluster_list[idx_biggest_cluster]
        components_of_biggest_cluster = cluster_tax[biggest_cluster]
        return components_of_biggest_cluster

    for i in index:
        left = i.split('-')[0]
        right = i.split('-')[1]
        if left not in tax_cluster and right not in tax_cluster:
            NewCluster = 'Cluster'+str(len(cluster_tax.keys())+1) #Create a new cluster label
            tax_cluster[left] = [NewCluster] # Create a key for this taxa with the new cluster
            tax_cluster[right] = [NewCluster]
            cluster_tax[NewCluster] = [left,right] # Add components to this new cluster
            cluster_ks[NewCluster] = [pairs_dict[i]] # And the kS value of this pair
        elif left not in tax_cluster and right in tax_cluster: # If one of the node is not in the record, starts with left missing
            components_of_biggest_cluster = get_biggest_cluster(right)
            NewCluster = 'Cluster'+str(len(cluster_tax.keys())+1) # Create a new cluster label for later
            tax_cluster[left] = [NewCluster] # Add a new record for the missing left node
            cluster_tax[NewCluster] = components_of_biggest_cluster+[left] # For the new cluster, add all nodes including the new left node
            cluster_ks[NewCluster]=[] # Create a new key for kS values of this cluster
            for j in components_of_biggest_cluster: # For each components
                tax_cluster[j].append(NewCluster) # First add the newly created cluster to their record
                tmp_pair = [j,left] # Combine the two nodes
                tmp_pair.sort() # Sort them
                tmp_pair = '-'.join(tmp_pair) # Make them the format as can be searched from paris_dict
                try:
                    cluster_ks[NewCluster].append(pairs_dict[tmp_pair]) # Add the kS value
                except:
                    continue
        elif left in tax_cluster and right not in tax_cluster:
            # Similar case as above condition, only now the right node is not in record
            components_of_biggest_cluster = get_biggest_cluster(left)
            NewCluster = 'Cluster'+str(len(cluster_tax.keys())+1)
            tax_cluster[right] = [NewCluster]
            cluster_tax[NewCluster] = components_of_biggest_cluster+[right]
            cluster_ks[NewCluster] = []
            for j in components_of_biggest_cluster:
                tax_cluster[j].append(NewCluster)
                tmp_pair = [j,right]
                tmp_pair.sort()
                tmp_pair = '-'.join(tmp_pair)
                try:
                    cluster_ks[NewCluster].append(pairs_dict[tmp_pair])
                except:
                    continue
        else: ### left in tax_cluster and right in tax_cluster:
            # For this condition, both nodes are already recorded
            # But there are still two scenarios
            # One is that the relationship of the two biggest clusters corresponding to them are the same cluster, do nothing
            components_of_biggest_cluster_left = get_biggest_cluster(left)
            components_of_biggest_cluster_right = get_biggest_cluster(right)
            if components_of_biggest_cluster_left == components_of_biggest_cluster_right:
                cluster_list = tax_cluster[left]
                cluster_size = [len(cluster_tax[each]) for each in cluster_list]
                cluster_size_w_idx = [[size,idx] for idx,size in enumerate(cluster_size)]
                cluster_size_w_idx.sort()
                idx_biggest_cluster = cluster_size_w_idx[-1][1]
                biggest_cluster = cluster_list[idx_biggest_cluster] # To be updated
                cluster_ks[biggest_cluster].append(pairs_dict[i])
                continue
            else:
                NewCluster = 'Cluster'+str(len(cluster_tax.keys())+1)
                tax_cluster[left].append(NewCluster)
                tax_cluster[right].append(NewCluster)
                cluster_tax[NewCluster] = components_of_biggest_cluster_left + components_of_biggest_cluster_right
                cluster_ks[NewCluster] = []
                for each_right in components_of_biggest_cluster_right:
                    tax_cluster[each_right].append(NewCluster)
                for each_left in components_of_biggest_cluster_left:
                    tax_cluster[each_left].append(NewCluster)
                    for each_right in components_of_biggest_cluster_right:
                        tmp_pair = [each_left,each_right]
                        tmp_pair.sort()
                        tmp_pair = '-'.join(tmp_pair)
                        try:
                            cluster_ks[NewCluster].append(pairs_dict[tmp_pair])
                        except:
                            continue
    return tax_cluster, cluster_tax, cluster_ks

def correct_ks(yn00):
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
    index, pairs_dict = loadcluster(dS_df_filtered,seq_names)
    tax_cluster, cluster_tax, cluster_ks = cluster_finder(index,pairs_dict)
    cluster_keys = cluster_ks.keys()
    ks_per_cluster = {}
    for i in cluster_keys:
        ks_per_cluster[i] = sum(cluster_ks[i])/len(cluster_ks[i])
    cluster_summary = {}
    for i in cluster_keys:
        cluster_summary[i] = [cluster_tax[i],ks_per_cluster[i]]
    ks_df = pd.DataFrame()
    ks_df['clusters'] = cluster_tax.values()
    ks_df['kS_values'] = ks_per_cluster.values()
    return ks_df

def draw_histo(ks_df, working_dir, output_prefix):
    # print ks_df
    if "DISPLAY" not in os.environ:
        print("The DISPLAY environment variable is not set. Plotting is aborted.\n"
              "The Ks distribution has been stored in your working directory in CSV format.\n"
              "Once you have your DISPLAY variable set, you can come back and re-plot the distribution\n"
              "using your method of interest, or use the plotting function of this pipeline by: \n"
              "python {Directory to WGD}/ks_correction.py {INPUT FILE PATH}")
    else:
        size = len(ks_df)
        ax = seaborn.distplot(a=ks_df["kS_values"], kde=False,kde_kws={"color":"red"}, axlabel=False,bins=100)
        ax.set(xlabel="kS", ylabel= "Frequency", title="kS distribution",xlim=(0, 5))
        plt.savefig(working_dir+"{0}_ks_distribution.pdf".format(output_prefix),format="pdf")
        plt.clf()

if __name__ == "__main__":
    ks_file = sys.argv[1]
    if "/" in ks_file:
        working_dir = "/".join(ks_file.split("/")[:-1])
    else:
        working_dir = "./"
    output_prefix = ".".join(ks_file.split(".")[:-1])
    ks_df = pd.read_csv(ks_file,header=0,index_col=0)
    draw_histo(ks_df,working_dir,output_prefix)
    print("The Ks distribution histogram has been plotted and saved in the same directory of your input file.")




