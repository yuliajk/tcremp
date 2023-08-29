import pandas as pd
import numpy as np
import os


def aminoacid_to_nt(cdr3):
    gene_code = {"A":"GCT", "C":"TGT", "T":"ACT", "S":"TCT", "F": "TTT", "L":"TTA", "I":"ATT", 
                 "M":"ATG", "V":"GTT", "P":"CCT", "Y":"TAT", "H":"CAT","Q": "CAA", "N":"AAT", 
                 "K":"AAA", "D":"GAT", "E":"GAA", "W":"TGG", "R":"CGT", "G":"GGT"}
    res = ''.join([gene_code.get(i) for i in list(cdr3)])
    return res

def filter_table(data,chain):
    if chain== "TRA":
        data = data[data['chain'] == "TRA"]
    
    if chain== "TRB":
        data = data[data['chain'] == "TRB"]
    #cols_input = list(data.columns)
    data = data[-data['v_gene'].isna()]
    data = data[-data['j_gene'].isna()]
    data = data[-(data['v_gene']== "None")]
    data = data[-(data['j_gene']== "None")]
    data = data[-(data['j_gene']== "")]
    data = data[-(data['v_gene']== "")]
    data['v_gene'] = data['v_gene'].str.split(',', 1, expand=True)[0]
    data['j_gene'] = data['j_gene'].str.split(',', 1, expand=True)[0]
    data['v_gene'] = data['v_gene'].str.split('*', 1, expand=True)[0]
    data['j_gene'] = data['j_gene'].str.split('*', 1, expand=True)[0]
    data = data[-data['v_gene'].str.contains(',')]
    data = data[-data['j_gene'].str.contains(',')]

    data = data[-(data['cdr3']== "")]
    data = data[-(data['cdr3']== "None")]
    data = data[-data['cdr3'].isna()]
    data = data[-data['cdr3'].str.contains('\.')]
    data = data[-data['cdr3'].str.contains('\*')]
    
    #data = data[data['high_confidence']==True]
    data = data.reset_index(drop=True)
    return data

def filter_table_vdjdb(data,chain):
    if chain== "TRA":
        data = data[data['chain'] == "TRA"]
    
    if chain== "TRB":
        data = data[data['chain'] == "TRB"]
    data = data[-data['v.segm'].isna()]
    data = data[-data['j.segm'].isna()]
    data = data[-(data['v.segm']== "None")]
    data = data[-(data['j.segm']== "None")]
    data = data[-(data['j.segm']== "")]
    data = data[-(data['v.segm']== "")]
    data['v.segm'] = data['v.segm'].str.split(',', 1, expand=True)[0]
    data['j.segm'] = data['j.segm'].str.split(',', 1, expand=True)[0]
    data['v.segm'] = data['v.segm'].str.split('*', 1, expand=True)[0]
    data['j.segm'] = data['j.segm'].str.split('*', 1, expand=True)[0]
    data = data[-data['v.segm'].str.contains(',')]
    data = data[-data['j.segm'].str.contains(',')]

    data = data[-(data['cdr3']== "")]
    data = data[-(data['cdr3']== "None")]
    data = data[-data['cdr3'].isna()]
    data = data[-data['cdr3'].str.contains('\.')]
    data = data[-data['cdr3'].str.contains('\*')]
    
    data = data[-data['antigen.epitope'].isna()]
    
    #data = data[data['high_confidence']==True]
    data = data.reset_index(drop=True)
    return data


def columns_prep(data):

    data = data.assign(d_gene = ".")
    data = data.assign(count = 1)
    data['DStart'] =  -1
    data['DEnd'] =  -1
    data['VEnd'] =  -1
    data['JStart'] =  -1
    data['freq'] =  -1
    data['cdr3nt'] = data['cdr3'].apply(lambda x: aminoacid_to_nt(x))
    
    data = data.rename(columns={"cdr3": "cdr3aa"})
    data = data.rename(columns={"v_gene": "v"})
    data = data.rename(columns={"j_gene": "j"})
    data = data.assign(subset = ".")
    #data = data.rename(columns={"value": "subset"})
    #data = data.rename(columns={"cdr3_nt": "cdr3nt"})

    cols = ["count", "freq", "cdr3aa", "cdr3nt", "v", "d", "j", "VEnd", "DStart", "DEnd","JStart","contig_id","reads","umis","length","subset"]
    data.drop(set(data.columns) - set(cols), axis=1, inplace=True)
    return data

def columns_prep_vdjdb(data):

    data = data.assign(d_gene = ".")
    data = data.assign(count = 1)
    data['DStart'] =  -1
    data['DEnd'] =  -1
    data['VEnd'] =  -1
    data['JStart'] =  -1
    data['freq'] =  -1
    data['cdr3nt'] = data['cdr3'].apply(lambda x: aminoacid_to_nt(x))
    
    data = data.rename(columns={"cdr3": "cdr3aa"})
    data = data.rename(columns={"v.segm": "v"})
    data = data.rename(columns={"j.segm": "j"})
    data = data.assign(subset = ".")
    #data = data.rename(columns={"value": "subset"})
    #data = data.rename(columns={"cdr3_nt": "cdr3nt"})

    cols = ["count", "freq", "cdr3aa", "cdr3nt", "v", "d", "j", "VEnd", "DStart", "DEnd","JStart","contig_id","reads","umis","length","subset"]
    data.drop(set(data.columns) - set(cols), axis=1, inplace=True)
    return data


def prototypes_prep(data):
    data.columns = ['cdr3nt','cdr3aa','v','j']
    cols = ["count", "freq", "cdr3aa", "cdr3nt", "v", "d", "j", "VEnd", "DStart", "DEnd","JStart","contig_id","reads","umis","length","subset"]
    data = data.assign(d = ".")
    data = data.assign(count = 1)
    data['DStart'] =  -1
    data['DEnd'] =  -1
    data['VEnd'] =  -1
    data['JStart'] =  -1
    data['freq'] =  -1
    data['reads'] =  1
    data['umis'] =  1
    data['length'] =  1
    data['subset'] =  '.'
    return data

#def data_filter_for_prototypes(data, prototypes_path):
#    prototypes=pd.read_csv(prototypes_path, sep = '\t')
#    data = data[data['v_gene'].isin(list(prototypes['v']))].reset_index(drop=True)
#    data = data[data['j_gene'].isin(list(prototypes['j']))].reset_index(drop=True)
#    return data
    
    
def data_filter_for_prototypes(data, prototypes_path, v_gene='v_gene',j_gene='j_gene'):
    prototypes=pd.read_csv(prototypes_path, sep = '\t')
    data = data[data[v_gene].isin(list(prototypes['v']))].reset_index(drop=True)
    data = data[data[j_gene].isin(list(prototypes['j']))].reset_index(drop=True)
    return data
    

def mir_clac(data, file_path_prefix, file_date, chain, prototypes_path, need_columns_prep=True):
    if need_columns_prep:
        data = columns_prep(data)
    file_path = 'data_scripts/' + file_path_prefix + file_date +'.txt'
    data.to_csv(file_path, sep='\t', index = False)
    
    mir_path = "mir-1.0-SNAPSHOT.jar"
    species = "Human"
    input_data_path = file_path
    output_path = 'data_scripts/' +file_path_prefix+ file_date
    
    command = "java -Xmx100G -cp " + mir_path + " com.antigenomics.mir.scripts.Examples cdr3aavj-pairwise-dist " + "-S " + species + " -G " + chain + " -F VDJtools " + "-I " + prototypes_path + " " + input_data_path + " -O " + output_path
    os.system(command)

#def mir_dists_format(data_dists_raw, data):
#    if len(data_dists_raw['id1'].drop_duplicates()) < len(data_dists_raw['id2'].drop_duplicates()):
#        first_index = 'id1'
#        second_index = 'id2'
#    else:
#        first_index = 'id2'
#        second_index = 'id1'            
#    data_dists_raw['cdr3_idx'] = 'cdr3_' + data_dists_raw[first_index].astype(str)
#    data_dists_raw['v_idx'] = 'v_' + data_dists_raw[first_index].astype(str)
#    data_dists_raw['j_idx'] = 'j_' + data_dists_raw[first_index].astype(str)
#    data_dists_raw = pd.concat([data_dists_raw.pivot(index=second_index,columns='cdr3_idx',values='cdr3.score').reset_index(),
#                     data_dists_raw.pivot(index=second_index,columns='v_idx',values='v.score').reset_index(),
#                     data_dists_raw.pivot(index=second_index,columns='j_idx',values='j.score').reset_index()], axis=1)
#    data_dists_raw = data_dists_raw.drop(second_index,axis=1)
#    data_dists_raw = data_dists_raw.set_index(data['barcode'],drop = True)
#    return data_dists_raw

def mir_dists_format(data_dists_raw, data, index_col,smaller = False):
    if ((len(data_dists_raw['id1'].drop_duplicates()) < len(data_dists_raw['id2'].drop_duplicates())) and (not smaller)) or ((len(data_dists_raw['id1'].drop_duplicates()) > len(data_dists_raw['id2'].drop_duplicates())) and smaller):
        first_index = 'id1'
        second_index = 'id2'
    else:
        first_index = 'id2'
        second_index = 'id1'            
    data_dists_raw['cdr3_idx'] = 'cdr3_' + data_dists_raw[first_index].astype(str)
    data_dists_raw['v_idx'] = 'v_' + data_dists_raw[first_index].astype(str)
    data_dists_raw['j_idx'] = 'j_' + data_dists_raw[first_index].astype(str)
    data_dists_raw = pd.concat([data_dists_raw.pivot(index=second_index,columns='cdr3_idx',values='cdr3.score').reset_index(),
                     data_dists_raw.pivot(index=second_index,columns='v_idx',values='v.score').reset_index(),
                     data_dists_raw.pivot(index=second_index,columns='j_idx',values='j.score').reset_index()], axis=1)
    data_dists_raw = data_dists_raw.drop(second_index,axis=1)
    data_dists_raw = data_dists_raw.set_index(data[index_col],drop = True)
    return data_dists_raw
   