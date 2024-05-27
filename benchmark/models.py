import numpy as np
import pandas as pd
import time, os, sys
import matplotlib.pyplot as plt
import seaborn as sns

import tcremb.ml_utils as ml_utils

tcr_columns_paired = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ']}

def run_GIANA(data_df, chain, output_suf,cpus=2, label = 'antigen.epitope'):
    '''Run GIANA clustering algorithm
'''
    df = data_df[~data_df[tcr_columns_paired[chain][0]].isna()][[tcr_columns_paired[chain][0],tcr_columns_paired[chain][1],label,'data_id']].reset_index(drop=True)
    df = df.rename({tcr_columns_paired[chain][0]:'CDR3',tcr_columns_paired[chain][1]:'V'},axis=1)
    # Reformat input for GIANA
    #seqs = df[['CDR3','V']]
    seqs = df[['CDR3','V']].drop_duplicates()

    #save data for GIANA
    cdir = os.getcwd()
    giana_path = os.path.join(cdir, 'benchmark/GIANA/')
    os.chdir(giana_path)
    seqs.to_csv('data.txt', index=False, header=False, sep='\t')
    print('GIANA clustering of {} sequences.'.format(len(df)))

    # Run GIANA
    t0 = time.time()
    os.system('python GIANA4.1.py -f data.txt -O data_clustered.txt -v True -N {}'.format(cpus))
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(giana_path, 'data_clustered.txt'), 'r') as f:
        results = f.read().splitlines()[3:]
        results = pd.DataFrame([x.split('\t') for x in results], columns=['CDR3',
                                                                            'cluster',
                                                                            'V',
                                                                         ])   
    os.chdir(cdir)
    
    giana_data = pd.merge(df, results.drop_duplicates(),on=['CDR3','V'])
    
    binom_res = ml_utils.binominal_test(giana_data, 'cluster', label)
    binom_res = binom_res.rename({label:'label_cluster'},axis=1)
    giana_data = giana_data.merge(binom_res)
    giana_data = df.merge(giana_data, how = 'left') ##new
    giana_data['is_cluster'] = giana_data['is_cluster'].fillna(0) ## new
    giana_data.to_csv(f'benchmark/outputs/giana_res_{chain}_{output_suf}.txt',sep='\t', index=False)

    return giana_data, t


def run_ismart(data_df, chain, output_suf, cpus=2, label = 'antigen.epitope'):
    
    df = data_df[~data_df[tcr_columns_paired[chain][0]].isna()][[tcr_columns_paired[chain][0],tcr_columns_paired[chain][1],label,'data_id']].reset_index(drop=True)
    df = df.rename({tcr_columns_paired[chain][0]:'CDR3',tcr_columns_paired[chain][1]:'V'},axis=1)
    
    # Reformat input for iSMART
    #seqs = df[['CDR3','V']]
    seqs = df[['CDR3','V']].drop_duplicates()

    #save data for iSMART
    cdir = os.getcwd()
    ismart_path = os.path.join(cdir, 'benchmark/iSMART/')
    os.chdir(ismart_path)
    seqs.to_csv('data.txt', index=False, header=False, sep='\t')
    print('Clustering {} sequences with iSMART.'.format(len(df)))
    
    
    t0 = time.time()
    #os.system('python iSMARTf3.py -f input.txt -v True -N {}'.format(cpus))
    os.system('python iSMARTf3.py -f data.txt')
    t1 = time.time()
    t = t1 - t0
    
    
    print('Elapsed time: {} seconds.'.format(t))

    with open(os.path.join(ismart_path, 'data_clustered_v3.txt'), 'r') as f:
        results = f.read().splitlines()[3:]
        results = pd.DataFrame([x.split('\t') for x in results], columns=['CDR3',
                                                                            'V',
                                                                            'cluster',
                                                                         ])   
    os.chdir(cdir)
    
    ismart_data = pd.merge(df, results.drop_duplicates(),on=['CDR3','V'])
    
    binom_res = ml_utils.binominal_test(ismart_data, 'cluster', label)
    binom_res = binom_res.rename({label:'label_cluster'},axis=1)
    ismart_data = ismart_data.merge(binom_res)
    ismart_data = df.merge(ismart_data, how = 'left') ##new
    ismart_data['is_cluster'] = ismart_data['is_cluster'].fillna(0) ## new
    ismart_data.to_csv(f'benchmark/outputs/ismart_res_{chain}_{output_suf}.txt',sep='\t')

    return ismart_data, t


def run_tcremb(data_path, chain, output_suf, skip_scores=False, label = 'antigen.epitope', model='kmeans' ):
    run_name = f'compare_{output_suf}'
    label_cl = label
    
    #cdir = os.getcwd()
    print('TCRemb clustering')

    # Run TCRemb
    t0 = time.time()
    if skip_scores == True:
            command = f'python tcremb_run.py --input {data_path} --runname {run_name} --chain {chain} --label {label_cl} --clstr_model {model} --data_id data_id --skip_scores True'
    else:
        command = f'python tcremb_run.py --input {data_path} --runname {run_name} --chain {chain} --label {label_cl} --clstr_model {model} --data_id data_id'
    print(command)
    os.system(command)
    t1 = time.time()
    t = t1 - t0

    print('Elapsed time: {} seconds.'.format(t))

    tcremb_data = pd.read_csv(f'tcremb_outputs/{run_name}/tcremb_clstr_res_{chain}.txt')
    #tcremb_data.to_csv(f'benchmark/outputs/tcremb_res_{chain}_{output_suf}.txt',sep='\t', index=False)

    return tcremb_data, t