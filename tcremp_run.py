import argparse, os
from pathlib import Path, PurePath
import numpy as np
import pandas as pd
from time import strftime, gmtime

import sys
sys.path.append("../")
sys.path.append("../mirpy/mirpy/")
import tcremp.TCRemP as TCRemP
import tcremp.TCRemP_clstr as TCRemP_clstr
import tcremp.ml_utils as ml_utils
import tcremp.data_proc as data_proc

import logging


tcr_columns = {'TRA':['a_cdr3aa','a_v','a_j'],
                    'TRB':['b_cdr3aa','b_v','b_j'],
                    'TRA_TRB':['a_cdr3aa','a_v','a_j', 'b_cdr3aa','b_v','b_j']}
tcr_columns_flat = ['cdr3aa','v','j','chain']
clone_label = 'unlabeled'
clone_index_columns = {'TRA': 'cloneId', 'TRB': 'cloneId', 'TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
species_glossary = {'homosapiens': 'HomoSapiens', 'human': 'HomoSapiens'}
    
def clustering(args, tcremp, output_path, output_columns):
    #output_columns = [tcremp.annotation_id,tcremp.clonotype_id] + tcr_columns_chain[args.chain]
    model = TCRemP_clstr.TCRemP_clustering(model_name = args.clstr_model)
    model.clstr(chain= args.chain, data= tcremp, label_cl=args.label, model = args.clstr_model)

    #df = kmeans.clstr_labels[args.chain].merge(tcremp.annot[args.chain][[tcremp.clonotype_index, tcremp.annotation_id,tcremp.clonotype_id,args.label, 'clone_size']])
    #df = model.clstr_labels[args.chain].merge(tcremp.annot[args.chain][output_columns])
    #df = tcremp.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', label_cluster,tcremp.annotation_id]])
    if args.label:
        df = tcremp.annot[args.chain][output_columns].merge(model.clstr_labels[args.chain][['cluster', clone_label, tcremp.annotation_id]])
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
    df.to_csv(f'{output_path}tcremp_clstr_res_{args.chain}.txt', sep='\t', index=False)
    
    

def main():
    parser = argparse.ArgumentParser(description='Basic TCRemP pipeline implementation')

    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to input file containing a clonotype (clone) table')
    
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output folder')
    
    parser.add_argument('-e', '--prefix', type=str,
                        help='Output prefix. Defaults to input clonotype table filename')
        
    parser.add_argument('-x', '--index-col', type=str,
                        help='(optional) Name of a column in the input table containing user-specified IDs that will be transfered to output tables')
    
    parser.add_argument('-c', '--chain', type=str, required=True,
                        help='"TRA" or "TRB" for single-chain input data (clonotypes), for paired-chain input (clones) use "TRA_TRB". Used in default prototype set selection')

    parser.add_argument('-p', '--prototypes-path', type=str,
                        help='Path to user-specified folder. If not set, will use default pre-built prototype tables, "$tcremp_path/data/data_prebuilt"')
    
    parser.add_argument('-n', '--n-prototypes', type=int,
                        help='Number of prototypes to select for clonotype "triangulation" during embedding. The total number of co-ordinates will be (number of chains) * (3 for V, J and CDR3 distances) * (n). Will use all available prototypes if not set')
    
    parser.add_argument('-s', '--species', type=str, default='HomoSapiens',
                        help='Prototype set species specification. Currently only "HomoSapiens" is supported')

    parser.add_argument('-u', '--unique',
                        help='Speed-up the analysis by running for unique clonotypes (clones) in the input table')
    
    parser.add_argument('-r', '--random-seed', type=int, default=42,
                        help='Random seed for prototype sampling and other rng-based procedures')    

    parser.add_argument('-a', '--cluster-algo', type=str, default='DBSCAN',
                        help='Embedding clustering algorithm: "DBSCAN", "K-means" or "None" to skip the step')
    
    parser.add_argument('-l', '--labels-col', type=str,
                        help='(optional) Name of a column in the input table containing clonotype labels. If provided, labels will be transferred to the output and various statistics will be calculated by comparing user-provided labels with inferred cluster labels')
    
    
    args = parser.parse_args()

    # IO setup
    output_path = Path(args.o)
    output_path.mkdir(parents=True, exist_ok=True)
    input_path = Path(args.i)
    output_prefix = args.e
    if not output_prefix:
        output_prefix = input_path.stem

    # Logging
    formatter = logging.Formatter('[%(asctime)s\t%(name)s\t%(levelname)s] %(message)s')
    logging.basicConfig(filename=f'{output_path}/{output_prefix}.log', 
                        format=formatter,
                        level=logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)    
    logging.getLogger().addHandler(handler)
    logging.info(f'Running TCRemP for i="{input_path.resolve()}", writing to o="{output_path.resolve()}/" under prefix="{output_prefix}"')

    # Load input
    if not input_path.is_file:
        logging.error(f'Missing input file "{input_path}"')
        sys.exit('Parameter error')    
    logging.info("Loading data and initializing the pipeline")
    data = pd.read_csv(input_path, sep='\t')

    # Check remaining parameters
    species = species_glossary.get(str(args.species).lower())
    if not species:
        logging.error(f'Bad species "{args.species}"')
        sys.exit('Parameter error')
    n_prototypes = args.n
    if n_prototypes:
        logging.debug(f'Will use {n_prototypes} prototypes')
        if n_prototypes> 3000:
            logging.warn('More than 3000 prototypes selected, may run very slowly')
    else:
        logging.debug('Will use all available prototypes')
    index_col = args.x
    if index_col:
        if index_col in data:
            logging.debug(f'Using {index_col} as clonotype index')
        else:
            logging.error(f'Index column "{index_col}" is missing in input data')
            sys.exit('Bad input')
    label_col = args.l
    if label_col:
        if label_col in data:
            logging.debug(f'Using {label_col} as clonotype labels')
        else:
            logging.error(f'Label column "{label_col}" is missing in input data')
            sys.exit('Bad input')

    # Setup pipeline
    tcremp = TCRemP.TCRemP(run_name = output_path, input_data = data, 
                           clonotype_index = args.clonotype_index,
                           prototypes_path = args.prototypes_path,
                           n = args.n, 
                           species = args.species,
                           prototypes_chain = args.chain, 
                           random_seed=args.random)
    
    logging.info("Checking input and extracting clonotypes")
    tcremp.tcremp_clonotypes(args.chain, args.unique_clonotypes)
    
    ## output columns
    output_columns = [tcremp.annotation_id, tcremp.clonotype_id] + tcr_columns[args.chain]
    if tcremp.clonotype_index:
        output_columns.append(tcremp.clonotype_index)
    if args.label:
        output_columns.append(args.label)
    
    ## count and save dists
    print('Stage: Distance scores calculation')
    tcremp.tcremp_dists_count(args.chain)
    tcremp.tcremp_dists(args.chain)        
    tcremp.annot[args.chain][output_columns].merge(tcremp.annot_dists[args.chain]).to_csv(f'{output_path}tcremp_dists_{args.chain}.txt', sep='\t', index=False)
    #dist_df = tcremp.annot[args.chain][output_columns].merge(tcremp.annot_dists[args.chain])
    #dist_df.to_csv(f'{output_path}tcremp_dists_{args.chain}.txt', sep='\t', index=False)
    
    ## pca
    print('Stage: PCA calculation')
    tcremp.tcremp_pca(args.chain)
    tcremp.annot[args.chain][output_columns].merge(tcremp.pca[args.chain]).to_csv(f'{output_path}tcremp_pca_{args.chain}.txt', sep='\t', index=False)
    
    ## tsne
    print('Stage: TSNE calculation')
    tcremp.tcremp_tsne(args.chain)
    tcremp.annot[args.chain][output_columns].merge(tcremp.tsne[args.chain]).to_csv(f'{output_path}tcremp_tsne_{args.chain}.txt', sep='\t', index=False)
    
    
    if args.clstr_model!='none':
        logging.info(f'Clustering with model {args.clstr_model}')
        print(f'Stage: Clustering with model {args.clstr_model}')
        clustering(args, tcremp, output_path, output_columns)
    else:
        print('Finished without clustering')
    
    print(f'Results are in {output_path}')


if __name__ == '__main__':
    main()