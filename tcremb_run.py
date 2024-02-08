import argparse
import math
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats
import statistics

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


sys.path.append("mirpy/")
sys.path.append("../")
import tcremb.TCRemb as TCRemb
import tcremb.ml_utils as ml_utils
import tcremb.data_proc as data_proc



tcr_columns = ['cdr3aa','v','j','chain']
clonotype_id_column = 'cloneId'
data_id= 'data_id'
annotation_id = 'annotId'
#pairing_id = 'barcode'
annotation_tcr_id_columns_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}


def clustering(args, tcremb,outputs_path):
    kmeans = TCRemb.TCRemb_clustering('KMeans')
    kmeans.clstr(args.chain, tcremb, args.label)
    kmeans.clstr_labels[args.chain].merge(tcremb.annot[args.chain][[annotation_id,args.label]]).to_csv(f'{outputs_path}tcremb_clstr_res_{args.chain}.txt', sep='\t', index=False)
    
def classsification(args, tcremb,outputs_path):
    clf = TCRemb.TCRemb_clf('RF')
    clf.clf(args.chain,tcremb, args.label)
    
def clusteting_pred():
    pass
    
def clussification_pred():
    pass
    


def main():
    parser = argparse.ArgumentParser(description='TCRemb is clustering and classification method for TCRs')

    parser.add_argument('--input', type=str,required=True,
                        help='Input file of TCRs')
    
    parser.add_argument('--runname', type=str,required=True,
                        help='Run name. Outputs will be stored in corresponding folder in tcremb_outputs/runname')
    
    parser.add_argument('--label', type=str,
                        help='label of clsasses of data for clustering and classification')
    
    parser.add_argument('--chain', type=str,required=True,
                        help='TRA or TRB if single, TRA_TRB if paired')
    
    parser.add_argument('--mode', type=str, default='clstr',
                        help='clstr - clusteting, clsf - classification, clstr_clsf - beoth, clstr_pred - clustering with train and pred datasets, clsf_pred - classification with traning and prediction of pred dataset, clstr_clsf_pred, scores - only scores count')
    
    parser.add_argument('--skip_scores', type=bool, default=False
                        , help='If score are already calculated pass skip_scores True')
    
    
    args = parser.parse_args()
    
    outputs_path= "tcremb_outputs/" + args.runname + '/'
    
    
    data_preped = pd.read_csv(args.input,sep='\t')
    
    print(args.mode)
    
    tcremb = TCRemb.TCRemb(args.runname, data_preped)
    tcremb.tcremb_clonotypes(args.chain)
    
    if not args.skip_scores:
        print(f'calculating dist scores for {args.chain} chain {args.input}')
        tcremb.tcremb_dists_count(args.chain)
    else:
        print('skip_scores was passed as parametr. continue withot scores calculation')

    if args.mode!='scores':
        print(f'calculating pca of dists scores for {args.chain} chain {args.input}')
        tcremb.tcremb_dists(args.chain)
        tcremb.tcremb_pca(args.chain)
        ##tcremb.tcremb_tsne(args.chain)
    else:
        print('Finished. only scores were calculated')
        
    if (args.mode=='clstr') or (args.mode=='clstr_clsf'):
        if args.label is None:
            print('Cant continue with clustering. Please provide label')
        else:
            clustering(args, tcremb, outputs_path)
            
    if (args.mode=='clsf') or (args.mode=='clstr_clsf'):
        if args.label is None:
            print('Cant continue with classification. Please provide label')
        else:
            clustering(args, tcremb, outputs_path)
    
    if (args.mode=='clstr_pred') or (args.mode=='clstr_clsf_pred'):
        if args.label is None:
            print('Cant continue with clustering. Please provide label')
        else:
            print('TODO')
            
    if (args.mode=='clsf') or (args.mode=='clstr_clsf_pred'):
        if args.label is None:
            print('Cant continue with traning of classificator. Please provide label')
        else:
            print('TODO')
    
    print(f'Results are in {outputs_path}')


if __name__ == '__main__':
    main()