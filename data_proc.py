import numpy as np
import pandas as pd

# Data processing 
##work on filter_clones_data
def annot_id(data, annotation_id_str):
    df = data.copy()
    df[annotation_id_str]=df.index
    return df

def remove_asterisk(data, tcr_columns):
    df = data.copy()
    df[tcr_columns[1]] = df[tcr_columns[1]].str.split('*',n=1,expand=True)[0]
    df[tcr_columns[2]] = df[tcr_columns[2]].str.split('*',n=1,expand=True)[0]
    return df

def remove_backslash(data, tcr_columns):
    df = data.copy()
    df[tcr_columns[1]] = df[tcr_columns[1]].str.replace('/','')
    df[tcr_columns[2]] = df[tcr_columns[2]].str.replace('/','')
    return df

def filter_clones_data(df_clones, tcr_columns):
    print(df_clones.shape)
    df_clones = df_clones[-df_clones[tcr_columns[0]].isna()]
    df_clones = df_clones[-df_clones[tcr_columns[1]].isna()]
    df_clones = df_clones[-df_clones[tcr_columns[2]].isna()]
    df_clones = df_clones[-df_clones[tcr_columns[3]].isna()]
    df_clones = df_clones[-df_clones[tcr_columns[0]].str.contains(',')]
    df_clones = df_clones[-df_clones[tcr_columns[1]].str.contains(',')]
    df_clones = df_clones[-df_clones[tcr_columns[2]].str.contains(',')]
    df_clones = df_clones[-df_clones[tcr_columns[3]].str.contains(',')]
    df_clones = df_clones[-df_clones[tcr_columns[0]].str.contains('\.')]
    df_clones = df_clones[-df_clones[tcr_columns[1]].str.contains('\.')]
    df_clones = df_clones[-df_clones[tcr_columns[2]].str.contains('\.')]
    df_clones = df_clones[-df_clones[tcr_columns[3]].str.contains('\.')]
    df_clones = df_clones[-df_clones[tcr_columns[0]].str.contains('\_')]
    df_clones = df_clones[-df_clones[tcr_columns[0]].str.contains('\*')]
    df_clones = df_clones.reset_index(drop=True)
    print(df_clones.shape)
    return df_clones


## 10x proc
def read_barcodes(barcodes_file):
    barcodes = pd.read_csv(barcodes_file, sep = '\t',header = None)
    barcodes.columns = ['barcode']
    barcodes.index += 1
    barcodes['barcode_id'] = barcodes.index
    return barcodes

def read_features(features_file):
    features = pd.read_csv(features_file, sep = '\t',header = None)
    features.index += 1
    features.columns = ['feature_code','value','type']
    features['feature_id'] = features.index
    features['feature_id'] = pd.to_numeric(features['feature_id'])
    return features


def read_matrix(matrix_file):
    matrix = pd.read_csv(matrix_file,sep = '\t')
    matrix = matrix.drop([0])
    matrix[['feature_id','barcode_id', 'count']] = matrix['%%MatrixMarket matrix coordinate integer general'].str.split(expand=True)
    matrix = matrix.drop(['%%MatrixMarket matrix coordinate integer general'], axis = 1)
    matrix = matrix.apply(pd.to_numeric)
    matrix['count']=matrix['count'].astype(int)
    return matrix

def get_value_matrix(matrix):
    v_t = matrix['value'].str.split('_',n=1,expand=True)
    matrix['value']=v_t[0]
    matrix['value_type']=v_t[1]
    return matrix


def get_barcode_top_tetramer(matrix):
    matrix = matrix.sort_values(by=['count'],ascending=False)
    tetramers = matrix.drop_duplicates('barcode')
    tetramers['top_tetramer'] = tetramers['value']
    tetramers =  tetramers[['barcode','top_tetramer','count']]
    return tetramers

def norm_logp(data, count_col):
    data[count_col]= data[count_col].apply(lambda x : math.log1p(x))
    return data

def pivot_data(data):
    data = data[['count','barcode','tetramer']]
    data = data.pivot_table('count','barcode','tetramer')
    data = data.fillna(0)
    return data