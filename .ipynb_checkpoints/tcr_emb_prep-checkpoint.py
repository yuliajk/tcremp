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
    data = data[-data['v_gene'].str.contains(',')]
    data = data[-data['j_gene'].str.contains(',')]

    data = data[-(data['cdr3']== "")]
    data = data[-(data['cdr3']== "None")]
    data = data[-data['cdr3'].isna()]
    data = data[-data['cdr3'].str.contains('\.')]
    data = data[-data['cdr3'].str.contains('\*')]
    
    data = data[data['high_confidence']==True]
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
    
   
    