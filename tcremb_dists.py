import argparse
from pathlib import Path, PurePath
import numpy as np
import pandas as pd
from time import strftime, gmtime

import sys
sys.path.append("../")
sys.path.append("../mirpy/mirpy/")
import tcremb.TCRemb as TCRemb
import tcremb.TCRemb_clstr as TCRemb_clstr
import tcremb.ml_utils as ml_utils
import tcremb.data_proc as data_proc



tcr_columns = ['cdr3aa','v','j','chain']
tcr_columns_chain = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ'],'TRA_TRB':['a_cdr3aa','TRAV','TRAJ', 'b_cdr3aa','TRBV','TRBJ']}
label_cluster = 'label_cluster'

annotation_tcr_id_columns_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
    
def clustering(args, tcremb, outputs_path, output_columns):
    #output_columns = [tcremb.annotation_id,tcremb.clonotype_id] + tcr_columns_chain[args.chain]
    model = TCRemb_clstr.TCRemb_clustering(model_name = args.clstr_model)
    model.clstr(chain= args.chain, data= tcremb, label_cl=args.label, model = args.clstr_model)

    #df = kmeans.clstr_labels[args.chain].merge(tcremb.annot[args.chain][[tcremb.data_id, tcremb.annotation_id,tcremb.clonotype_id,args.label, 'clone_size']])
    #df = model.clstr_labels[args.chain].merge(tcremb.annot[args.chain][output_columns])
    df = tcremb.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', label_cluster,tcremb.annotation_id]])
    if args.label:
        model.clstr_metrics_calc(args.chain, tcremb)
        #print(f"purity:{model.clstr_metrics[args.chain]['purity']}")
        print(f"retention:{model.clstr_metrics[args.chain]['retention']}")
        print(f"f1-score:{model.clstr_metrics[args.chain]['f1-score']}")
        print(f"total pairs TCR-epitope:{model.clstr_metrics[args.chain]['total pairs TCR-epitope']}")
        print(f"total unique epitopes:{model.clstr_metrics[args.chain]['total unique epitopes']}")
    
    #df = df.merge(tcremb.tsne[args.chain])
    #df = df.merge(tcremb.tsne_clones[args.chain].rename({'DM1':'DM1_clones','DM2':'DM2_clones'},axis=1))
    df.to_csv(f'{outputs_path}tcremb_clstr_res_{args.chain}.txt', sep='\t', index=False)
    
    

def main():
    parser = argparse.ArgumentParser(description='TCRemb dists')

    parser.add_argument('-i','--input', type=str,required=True,
                        help='Input file of TCRs')
    
    parser.add_argument('-o','--output', type=str,#required=True,
                        help='Output directory path. Outputs will be stored in corresponding directory. If None -  tcremb_outputs/input_filename')
        
    parser.add_argument('--data_id', type=str,
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
        outputs_path= args.output + '/'
    else:
        output = 'tcremb_' + PurePath(args.input).name.replace('.','')
        outputs_path= "tcremb_outputs/" + output + '/'

    print(f'calculating dists scores, pca and tsne for {args.chain} chain {args.input}')      
    print(f'results and temp files will be in {outputs_path}')
    
    data_preped = pd.read_csv(args.input,sep='\t')
    
    
    tcremb = TCRemb.TCRemb(run_name = outputs_path, input_data = data_preped, data_id=args.data_id, prototypes_path= args.prototypes_path, n = args.n, species = args.species, prototypes_chain = args.chain, random_seed=args.random)
    
    tcremb.tcremb_clonotypes(args.chain, args.unique_clonotypes)
    
    ## output columns
    output_columns = [tcremb.annotation_id,tcremb.clonotype_id] + tcr_columns_chain[args.chain]
    if tcremb.data_id:
        output_columns.append(tcremb.data_id)
    if args.label:
        output_columns.append(args.label)
    
    ## count and save dists
    tcremb.tcremb_dists_count(args.chain)
    tcremb.tcremb_dists(args.chain)        
    tcremb.annot[args.chain][output_columns].merge(tcremb.annot_dists[args.chain]).to_csv(f'{outputs_path}tcremb_dists_{args.chain}.txt', sep='\t', index=False)
    #dist_df = tcremb.annot[args.chain][output_columns].merge(tcremb.annot_dists[args.chain])
    #dist_df.to_csv(f'{outputs_path}tcremb_dists_{args.chain}.txt', sep='\t', index=False)
    
    ## pca
    tcremb.tcremb_pca(args.chain)
    tcremb.annot[args.chain][output_columns].merge(tcremb.pca[args.chain]).to_csv(f'{outputs_path}tcremb_pca_{args.chain}.txt', sep='\t', index=False)
    
    ## tsne
    tcremb.tcremb_tsne(args.chain)
    tcremb.annot[args.chain][output_columns].merge(tcremb.tsne[args.chain]).to_csv(f'{outputs_path}tcremb_tsne_{args.chain}.txt', sep='\t', index=False)
    
    if args.clstr_model!='none':
        print(f'Clustering with model {args.clstr_model}')
        clustering(args, tcremb, outputs_path, output_columns)
    else:
        print('Finished without clustering')
    
    print(f'Results are in {outputs_path}')


if __name__ == '__main__':
    main()