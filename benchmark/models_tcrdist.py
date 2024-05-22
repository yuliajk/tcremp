import numpy as np
import pandas as pd
import time, os, sys


from tcrdist.repertoire import TCRrep
from tcrdist.rep_funcs import  compute_pw_sparse_out_of_memory2
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN

import tcremb.ml_utils as ml_utils

tcr_columns_paired = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ']}

def tcrdist3_dist(data_df, chain, n_clusters=None, chunk=True, label='antigen.epitope'):
    df = data_df.copy()
    cdr3a = 'cdr3_a_aa'
    va = 'v_a_gene'
    ja = 'j_a_gene'
    cdr3b = 'cdr3_b_aa'
    vb = 'v_b_gene'
    jb = 'j_b_gene'

    if not 'count' in df.columns:
        df['count']=[1]*len(df)

    if chain == 'TRA':
        df = df.rename({tcr_columns_paired[chain][0]:cdr3a,tcr_columns_paired[chain][1]:va,tcr_columns_paired[chain][2]:ja},axis=1)
        #df = df[df['chain']==chain].reset_index(drop=True)
        #df = df.rename(columns={'cdr3aa': cdr3a,
        #                                        'v': va,
        #                                        'j':ja})
        df = df[[cdr3a, va, ja,label, 'count','data_id']]

    elif chain == 'TRB':
        df = df.rename({tcr_columns_paired[chain][0]:cdr3b,tcr_columns_paired[chain][1]:vb,tcr_columns_paired[chain][2]:jb},axis=1)
        #df = df[df['chain']==chain].reset_index(drop=True)
        #df= df.rename(columns={'cdr3aa': cdr3b,
        #                                        'v': vb,
        #                                        'j':jb})
        df = df[[cdr3b, vb, jb,label, 'count','data_id']]

    elif chain == 'TRA_TRB':
        df = df.rename({tcr_columns_paired['TRA'][0]:cdr3a,tcr_columns_paired['TRA'][1]:va,tcr_columns_paired['TRA'][2]:ja,tcr_columns_paired['TRB'][0]:cdr3b,tcr_columns_paired['TRB'][1]:vb,tcr_columns_paired['TRB'][2]:jb},axis=1)
        df = df[[cdr3a,va,ja,cdr3b,vb,jb,label, 'count','data_id']]
    else:
        pass


    seqs = df.drop(columns=[label,'data_id'], axis=1).reset_index(drop=True)

    if chain =='TRA':
        chain1 = ['alpha']
    elif chain =='TRB':
        chain1 = ['beta']
    else:
        chain1=['alpha','beta']

    # Run tcrdist3

    print('\n*** Tcrdist3 clustering %s %s chains ' % (len(seqs), chain))

    t0 = time.time()

    # Create tcrdist object
    tr = TCRrep(cell_df=seqs,   # Input data
                organism='human',   # Organism
                chains=chain1,       # Chain selection
                infer_all_genes=True,   # Infer genes
                infer_cdrs=True,        # Infer CDRs
                infer_index_cols=True,  # Infer column names
                store_all_cdr=True,     # Store CDRs
                deduplicate=False,      # Deduplicate
                compute_distances=True) 
                #compute_distances=False) 
    # Compute distances
    return tr, df
    
    # Compute tcrdist distances using sparse rect distances
def tcrdist3_compute(tr, chain, cpus, radius=50):
    if chain =='TRA':
        chain1 = ['alpha']
        name = 'alpha'
    elif chain=='TRB':
        chain1 = ['beta']
        name = 'beta'
    else:
        chain1 = ['alpha','beta']
        name = 'alpha_beta'


    S, _ = compute_pw_sparse_out_of_memory2(tr = tr,
        row_size      = 50,
        pm_processes  = cpus,
        pm_pbar       = True,
        max_distance  = radius,
        reassemble    = True,
        cleanup       = True,
        assign        = True)
    print(S.keys())
    S=S[name]

    return S   


def silhouette_clusters(X_data):    
    data_len = len(X_data)
    
    range_n_clusters = [
        #round(data_len*0.005),
        round(data_len*0.01), round(data_len*0.05) , round(data_len*0.1), round(data_len*0.15)
        ,round(data_len*0.2), round(data_len*0.25), round(data_len*0.3)
        , round(data_len*0.4), round(data_len*0.5), round(data_len*0.6), round(data_len*0.7), round(data_len*0.8), round(data_len*0.9)]
    
    silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)
    
    n_to_try = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]
    
    range_n_clusters = [round(n_to_try*0.7), round(n_to_try*0.8), round(n_to_try*0.9), n_to_try, round(n_to_try*1.1) , round(n_to_try*1.2) , round(n_to_try*1.3)]
    silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)

    silhouette_n_clusters = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]

# Cluster tcrdist matrix
def tcrdist3_cluster_kmeans(S, chain, n_clusters):
    if not n_clusters:
        n_clusters=500
    kmeans = KMeans(init='random',
                    n_init=10,
                    n_clusters=int(n_clusters)).fit(S)
    labels = kmeans.labels_
    
    pd.DataFrame(labels).to_csv()

    return labels

def tcrdist3_cluster_dbscan(S, chain, eps):
    dbscan = DBSCAN(eps=eps, min_samples=2).fit(S)
    labels = dbscan.labels_
    
    pd.DataFrame(labels).to_csv()

    return labels

def run_tcrdist3(vdjdb_data_tcrdist3, chain, output_suf, label = 'antigen.epitope', eps = None):
    #n_clusters
    
    t0 = time.time()
    tcrdist3_tr, tcrdist_data = tcrdist3_dist(vdjdb_data_tcrdist3,chain)
    tcrdist3_s = tcrdist3_compute(tcrdist3_tr,chain,2)
    if eps:
        tcrdist3_labels = tcrdist3_cluster_dbscan(tcrdist3_s, chain, eps)
    else:
        silhouette_n_clusters = silhouette_avg_scores_kmeans(tcrdist3_s)
        tcrdist3_labels = tcrdist3_cluster_kmeans(tcrdist3_s, chain, silhouette_n_clusters)
    t1 = time.time()
    t = t1 - t0
    
    tcrdist_data['cluster']= tcrdist3_labels
    binom_res = ml_utils.binominal_test(tcrdist_data, 'cluster', label)
    binom_res = binom_res.rename({label:'label_cluster'},axis=1)
    tcrdist_data = tcrdist_data.merge(binom_res)
    tcrdist_data.to_csv(f'benchmark/outputs/tcrdist_res_{chain}_{output_suf}.txt',sep='\t')
    return tcrdist_data, t

def run_tcrdist3_set(vdjdb_data_tcrdist3, chain, output_suf, label = 'antigen.epitope', eps = [1,2]):
    tcrdist_data_dict = {}
    
    tcrdist3_tr, tcrdist_data = tcrdist3_dist(vdjdb_data_tcrdist3,chain)
    #tcrdist3_s = tcrdist3_compute(tcrdist3_tr,chain,2)
    if chain =='TRA':
        tcrdist3_s = pd.DataFrame(tcrdist3_tr.pw_alpha)
    if chain =='TRB':
        tcrdist3_s = pd.DataFrame(tcrdist3_tr.pw_beta)
    if chain =='TRA_TRB':
        tcrdist3_s = pd.concat([pd.DataFrame(tcrdist3_tr.pw_alpha),pd.DataFrame(tcrdist3_tr.pw_beta)],axis=1)
    for e in eps:
        tcrdist_data_e = tcrdist_data.copy()
        tcrdist3_labels = tcrdist3_cluster_dbscan(tcrdist3_s, chain, e)
        tcrdist_data_e['cluster']= tcrdist3_labels
        binom_res = ml_utils.binominal_test(tcrdist_data_e, 'cluster', label)
        binom_res = binom_res.rename({label:'label_cluster'},axis=1)
        tcrdist_data_e = tcrdist_data_e.merge(binom_res)
        tcrdist_data_dict[e]=tcrdist_data_e
    return tcrdist_data_dict

