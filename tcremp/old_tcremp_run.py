import argparse
from pathlib import Path, PurePath
import numpy as np
import pandas as pd
from time import strftime, gmtime
import os

import sys
sys.path.append("../")
sys.path.append("../mirpy/mirpy/")
import tcremp.TCRemP as TCRemP
import tcremp.tcremp_cluster as tcremp_cluster
import tcremp.ml_utils as ml_utils
import tcremp.data_proc as data_proc

from datetime import datetime
import logging



tcr_columns = ['cdr3aa','v','j','chain']
#tcr_columns_chain = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ'],'TRA_TRB':['a_cdr3aa','TRAV','TRAJ', 'b_cdr3aa','TRBV','TRBJ']}
tcr_columns_chain = {'TRA':['a_cdr3aa','a_v','a_j'],'TRB':['b_cdr3aa','b_v','b_j'],'TRA_TRB':['a_cdr3aa','a_v','a_j', 'b_cdr3aa','b_v','b_j']}
label_cluster = 'label_cluster'

annotation_tcr_id_columns_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
    
def clustering(args, tcremp, outputs_path, output_columns):
    #output_columns = [tcremp.annotation_id,tcremp.clonotype_id] + tcr_columns_chain[args.chain]
    model = tcremp_cluster.TcrempClustering(model_name = args.clstr_model)
    model.build_clusters(chain= args.chain, data= tcremp, label_cl=args.label, model = args.clstr_model)

    #df = kmeans.clstr_labels[args.chain].merge(tcremp.annot[args.chain][[tcremp.clonotype_index, tcremp.annotation_id,tcremp.clonotype_id,args.label, 'clone_size']])
    #df = model.clstr_labels[args.chain].merge(tcremp.annot[args.chain][output_columns])
    #df = tcremp.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', label_cluster,tcremp.annotation_id]])
    if args.label:
        df = tcremp.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', label_cluster,tcremp.annotation_id]])
        model.clstr_metrics_calc(args.chain, tcremp)
        ##print(f"purity:{model.clstr_metrics[args.chain]['purity']}")
        print(f"retention:{model.clstr_metrics[args.chain]['retention']}")
        print(f"f1-score:{model.clstr_metrics[args.chain]['f1-score']}")
        print(f"total pairs TCR-epitope:{model.clstr_metrics[args.chain]['total pairs TCR-epitope']}")
        print(f"total unique epitopes:{model.clstr_metrics[args.chain]['total unique epitopes']}")
        logging.info(f"purity:{model.clstr_metrics[args.chain]['purity']}")
        logging.info(f"retention:{model.clstr_metrics[args.chain]['retention']}")
        logging.info(f"f1-score:{model.clstr_metrics[args.chain]['f1-score']}")
        logging.info(f"total pairs TCR-epitope:{model.clstr_metrics[args.chain]['total pairs TCR-epitope']}")
        logging.info(f"total unique epitopes:{model.clstr_metrics[args.chain]['total unique epitopes']}")
    else:
        df = tcremp.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', tcremp.annotation_id]])
    
    #df = df.merge(tcremp.tsne[args.chain])
    #df = df.merge(tcremp.tsne_clones[args.chain].rename({'DM1':'DM1_clones','DM2':'DM2_clones'},axis=1))
    df.to_csv(f'{outputs_path}tcremp_clstr_res_{args.chain}.txt', sep='\t', index=False)
    
    

def main():
    parser = argparse.ArgumentParser(description='tcremp dists')

    parser.add_argument('-i','--input', type=str,required=True,
                        help='Input file of TCRs')
    
    parser.add_argument('-o','--output', type=str,#required=True,
                        help='Output directory path. Outputs will be stored in corresponding directory. If None -  tcremp_outputs/input_filename')
        
    parser.add_argument('--clonotype_index', type=str,
                        help='column with id of ypur input data. if you would like this id to be added to output of clustering or classifications')
    
    parser.add_argument('-c','--chain', type=str,required=True,
                        help='TRA or TRB if single, TRA_TRB if paired')

    parser.add_argument('-n','--n', type=int,
                        help='number pf prototypes randomly selected for embedding')
    
    parser.add_argument('-s','--species', type=str, default='HomoSapiens',
                        help='subset of prototypes to be used. HomoSapiens, MacacaMulatta')

    parser.add_argument('--unique_clonotypes', type=bool, default=False
                        , help='Left only unique clonotypes in annotation dataset')
    
    parser.add_argument('-r','--random', type=int,
                        help='random seed for randomlly selecting n prototypes. If not provided, first n prototypes are selected.')    
    parser.add_argument('-p','--prototypes_path', type=str,
                        help='path to input prototypes, if user would like to use custom prototypes')

    parser.add_argument('-m','--clstr_model', type=str, default='dbscan',
                        help='name of clustering algorithm, dbscan or kmeans. to skip clustering use none')
    parser.add_argument('-l','--label', type=str,
                        help='label of clsasses of data for clustering and classification. If provided metrics will be calculated')
    
    
    args = parser.parse_args()
    
    if args.output:
        outputs_path= args.output
    else:
        output = 'tcremp_' + PurePath(args.input).name.replace('.','')
        outputs_path= "tcremp_outputs/" + output
    
    outputs_path = os.path.join(outputs_path, '')
    Path(outputs_path).mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(filename=f'{outputs_path}tcremp_log.log', level=logging.DEBUG)
    logging.info(f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} start of tcremp run {outputs_path}')

    #print(f'calculating dists scores, pca and tsne for {args.chain} chain {args.input}')      
    print(f'results and temp files will be in {outputs_path}')
    
    data_preped = pd.read_csv(args.input,sep='\t')
    
    print(f'Running tcremp with: input data {args.input},output directory {outputs_path}, chain {args.chain}')
    tcremp = TCRemP.TCRemP(run_name = outputs_path, input_data = data_preped, clonotype_index=args.clonotype_index, prototypes_path= args.prototypes_path, n = args.n, species = args.species, prototypes_chain = args.chain, random_seed=args.random)
    
    print('Stage: Data cleaning and clonotypes extraction')
    tcremp.tcremp_clonotypes(args.chain, args.unique_clonotypes)
    
    ## output columns
    output_columns = [tcremp.annotation_id,tcremp.clonotype_id] + tcr_columns_chain[args.chain]
    if tcremp.clonotype_index:
        output_columns.append(tcremp.clonotype_index)
    if args.label:
        output_columns.append(args.label)
    
    ## count and save dists
    print('Stage: Distance scores calculation')
    tcremp.tcremp_dists_count(args.chain)
    tcremp.tcremp_dists(args.chain)        
    tcremp.annot[args.chain][output_columns].merge(tcremp.annot_dists[args.chain]).to_csv(f'{outputs_path}tcremp_dists_{args.chain}.txt', sep='\t', index=False)
    #dist_df = tcremp.annot[args.chain][output_columns].merge(tcremp.annot_dists[args.chain])
    #dist_df.to_csv(f'{outputs_path}tcremp_dists_{args.chain}.txt', sep='\t', index=False)
    
    ## pca
    print('Stage: PCA calculation')
    tcremp.tcremp_pca(args.chain)
    tcremp.annot[args.chain][output_columns].merge(tcremp.pca[args.chain]).to_csv(f'{outputs_path}tcremp_pca_{args.chain}.txt', sep='\t', index=False)
    
    ## tsne
    print('Stage: TSNE calculation')
    tcremp.tcremp_tsne(args.chain)
    tcremp.annot[args.chain][output_columns].merge(tcremp.tsne[args.chain]).to_csv(f'{outputs_path}tcremp_tsne_{args.chain}.txt', sep='\t', index=False)
    
    
    if args.clstr_model!='none':
        logging.info(f'Clustering with model {args.clstr_model}')
        print(f'Stage: Clustering with model {args.clstr_model}')
        clustering(args, tcremp, outputs_path, output_columns)
    else:
        print('Finished without clustering')
    
    print(f'Results are in {outputs_path}')


if __name__ == '__main__':
    main()