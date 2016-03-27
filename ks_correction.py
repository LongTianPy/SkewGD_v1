#!/usr/bin/python
"""
This script take the result of CODEML as input
"""

### IMPORT
import pycodeml
import pandas as pd
from operator import itemgetter

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
    return pairs

def cluster_finder(pairs):
    tip_left = [i[0] for i in pairs]
    tip_right = [i[1] for i in pairs]
    ks_value = [i[2] for i in pairs]
    connection = {} # To store connections that have been iterated
    closed_cluster = []
    n_pairs = range(len(tip_right))
    # for i in n_pairs:
    #     print i
    #     if tip_left[i] not in connection and tip_right[i] not in connection:
    #         print "Both not in record."
    #         connection[tip_left[i]] = [tip_right[i]]
    #         connection[tip_right[i]] = [tip_left[i]]
    #         closed_cluster.append([tip_left[i],tip_right[i],ks_value[i]])
    #         print "Cluster ", tip_left[i], " ", tip_right[i], " recorded"
    #     elif tip_left[i] not in connection and tip_right[i] in connection:
    #         print "Left tip not recorded, right tip in record."
    #         connection_to_find = connection[tip_right[i]]
    #         connection[tip_left[i]] = [tip_right[i]]
    #         connection[tip_right[i]].append(tip_left[i])
    #         tmp = [tip_left[i],tip_right[i],ks_value]
    #         n_connection_to_find = len(connection_to_find)-1
    #         # Now we are looking for pairs of tip_left[i]-each_one_of_connection_to_find or each_one_of_connection_
    #         # to_find-tip_left[i]
    #         for each_connection in connection_to_find:
    #             for j in n_pairs[i+1:]:
    #                 if tip_left[j] == tip_left[i] and tip_right[j] == each_connection:
    #                     tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                     connection[tip_left[j]].append(tip_right[j])
    #                     connection[tip_right[j]].append(tip_left[j])
    #                     tip_left.pop(j)
    #                     tip_right.pop(j)
    #                     ks_value.pop(j)
    #                     n_connection_to_find = n_connection_to_find - 1
    #                     break
    #                 elif tip_right[j] == tip_left[i] and tip_left[j] == each_connection:
    #                     tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                     connection[tip_left[j]].append(tip_right[j])
    #                     connection[tip_right[j]].append(tip_left[j])
    #                     n_connection_to_find = n_connection_to_find - 1
    #                     tip_left.pop(j)
    #                     tip_right.pop(j)
    #                     ks_value.pop(j)
    #                     break
    #                 else:
    #                     continue
    #         closed_cluster.append(tmp)
    #         print tmp
    #     elif tip_left[i] in connection and tip_right[i] not in connection:
    #         connection[tip_right[i]]=[tip_left[i]]
    #         connection_to_find = connection[tip_left[i]]
    #         connection[tip_left[i]].append(tip_right[i])
    #         tmp = [tip_left[i],tip_right[i],ks_value[i]]
    #         print tmp
    #         n_connection_to_find = len(connection_to_find)-1
    #         print n_connection_to_find
    #         for each_connection in connection_to_find:
    #             for j in n_pairs[i+1:]:
    #                 if tip_left[j] == tip_left[i] and tip_right[j] == each_connection:
    #                     tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                     connection[tip_left[j]].append(tip_right[j])
    #                     connection[tip_right[j]].append(tip_left[j])
    #                     n_connection_to_find = n_connection_to_find - 1
    #                     tip_left.pop(j)
    #                     tip_right.pop(j)
    #                     ks_value.pop(j)
    #                     break
    #                 elif tip_right[j] == tip_left[i] and tip_left[j] == each_connection:
    #                     tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                     connection[tip_left[j]].append(tip_right[j])
    #                     connection[tip_right[j]].append(tip_left[j])
    #                     n_connection_to_find = n_connection_to_find - 1
    #                     tip_left.pop(j)
    #                     tip_right.pop(j)
    #                     ks_value.pop(j)
    #                     break
    #                 else:
    #                     continue
    #         closed_cluster.append(tmp)
    #         print tmp
    #     elif tip_left[i] in connection and tip_right[i] in connection:
    #         print "Both in record"
    #         connection_to_left_to_find = connection[tip_left[i]]
    #         connection_to_right_to_find = connection[tip_right[i]]
    #         connection[tip_left[i]].append(tip_right[i])
    #         connection[tip_right[i]].append(tip_left[i])
    #         tmp = [tip_left[i],tip_right[i],ks_value[i]]
    #         n_connection_to_find = len(connection_to_left_to_find)*len(connection_to_right_to_find)-1
    #         for each_connection_to_left in connection_to_left_to_find:
    #             for each_connection_to_right in connection_to_right_to_find:
    #                 for j in n_pairs[i+1:]:
    #                     if tip_left[j] == each_connection_to_left and tip_right[j] == each_connection_to_right:
    #                         tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                         connection[tip_left[j]].append(tip_right[j])
    #                         connection[tip_right[j]].append(tip_left[j])
    #                         n_connection_to_find = n_connection_to_find - 1
    #                         tip_left.pop(j)
    #                         tip_right.pop(j)
    #                         ks_value.pop(j)
    #                         break
    #                     elif tip_right[j] == each_connection_to_right and tip_left[j] == each_connection_to_left:
    #                         tmp.append([tip_left[j],tip_right[j],ks_value[j]])
    #                         connection[tip_left[j]].append(tip_right[j])
    #                         connection[tip_right[j]].append(tip_left[j])
    #                         n_connection_to_find = n_connection_to_find - 1
    #                         tip_left.pop(j)
    #                         tip_right.pop(j)
    #                         ks_value.pop(j)
    #                         break
    #                     else:
    #                         continue
    #         closed_cluster.append(tmp)
    #         print tmp
    for i in n_pairs:
        print i
        print n_pairs
        if tip_left[i] not in connection and tip_right[i] not in connection:
            connection[tip_left[i]] = [tip_right[i]]
            connection[tip_right[i]] = [tip_left[i]]
            closed_cluster.append([tip_left[i],tip_right[i],ks_value[i]])
        elif tip_left[i] in connection and tip_right[i] not in connection:
            connection_to_find = connection[tip_left[i]]
            connection[tip_right[i]]=[tip_left[i]]
            connection[tip_left[i]].append(tip_right[i])
            tmp_tips = []
            tmp_ks = []
            tmp_tips.append(tip_left[i])
            tmp_tips.append(tip_right[i])
            tmp_ks.append(ks_value[i])
            try:
                for each_connection in connection_to_find:
                    idx = n_pairs.index(i)
                    for j in n_pairs[idx+1:]:
                        if tip_left[j] == each_connection and tip_right[j] == tip_right[idx]:
                            print "If",j
                            tmp_tips.append(tip_left[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            print "Here"
                            break
                        elif tip_left[j]==tip_right[i] and tip_right[j]==each_connection:
                            print "Elif",j
                            tmp_tips.append(tip_right[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            break
                        elif i == n_pairs[-1]:
                            break
                        else:
                            print "Else",j
                            continue
                closed_cluster.append([tmp_tips,tmp_ks])
            except:
                closed_cluster.append([tmp_tips,tmp_ks])
        elif tip_left[i] not in connection and tip_right[i] in connection:
            connection_to_find = connection[tip_right[i]]
            connection[tip_left[i]] = [tip_right[i]]
            connection[tip_right[i]].append(tip_left[i])
            tmp_tips = []
            tmp_ks = []
            tmp_tips.append(tip_left[i])
            tmp_tips.append(tip_right[i])
            tmp_ks.append(ks_value[i])
            try:
                for each_connection in connection_to_find:
                    idx = n_pairs.index(i)
                    for j in n_pairs[idx+1:]:
                        if tip_left[j] == each_connection and tip_right[j] == tip_right[i]:
                            tmp_tips.append(tip_left[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            break
                        elif tip_left[j]==tip_right[i] and tip_right[j]==each_connection:
                            tmp_tips.append(tip_right[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            break
                        else:
                            continue
                closed_cluster.append([tmp_tips,tmp_ks])
            except:
                closed_cluster.append([tmp_tips,tmp_ks])
        else: # When both tips are in record, chances are that this could be the last one in the list
            connection_to_find_right = connection[tip_right[i]]
            connection_to_find_left = connection[tip_left[i]]
            connection[tip_left[i]] = [tip_right[i]]
            connection[tip_right[i]].append(tip_left[i])
            tmp_tips = []
            tmp_ks = []
            tmp_tips.append(tip_left[i])
            tmp_tips.append(tip_right[i])
            tmp_ks.append(ks_value[i])
            # If there are records behind
            try:
                for each_connection in connection_to_find:
                    idx = n_pairs.index(i)
                    for j in n_pairs[idx+1:]:
                        if tip_left[j] == each_connection and tip_right[j] == tip_right[i]:
                            tmp_tips.append(tip_left[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            break
                        elif tip_left[j]==tip_right[i] and tip_right[j]==each_connection:
                            tmp_tips.append(tip_right[j])
                            tmp_ks.append(ks_value[j])
                            connection[tip_left[j]].append(tip_right[j])
                            connection[tip_right[j]].append(tip_left[j])
                            n_pairs.pop(j)
                            # tip_left.pop(j)
                            # tip_right.pop(j)
                            # ks_value.pop(j)
                            break
                        else:
                            continue
                closed_cluster.append([tmp_tips,tmp_ks])
            except:
                closed_cluster.append([tmp_tips,tmp_ks])
    return connection, closed_cluster






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




