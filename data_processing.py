import pandas as pd
import numpy as np
import os
from scipy import stats
from statistics import mean
from sklearn.metrics.pairwise import euclidean_distances
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn import neighbors
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import label_binarize
from sklearn import preprocessing
from sklearn.feature_selection import f_classif

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

def get_features_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()
        features = pd.DataFrame(feature_ref)
        features['id'] = features['id'].str.decode('utf-8')
        features['name'] = features['name'].str.decode('utf-8')
        features['feature_type'] = features['feature_type'].str.decode('utf-8')
        features['genome'] = features['genome'].str.decode('utf-8')
        
        return features


def read_matrix(matrix_file):
    matrix = pd.read_csv(matrix_file,sep = '\t')
    matrix = matrix.drop([0])
    matrix[['feature_id','barcode_id', 'count']] = matrix['%%MatrixMarket matrix coordinate integer general'].str.split(expand=True)
    matrix = matrix.drop(['%%MatrixMarket matrix coordinate integer general'], axis = 1)
    matrix = matrix.apply(pd.to_numeric)
    matrix['count']=matrix['count'].astype(int)
    return matrix

def merge_matrix(matrix,barcodes,features):
    matrix = pd.merge(matrix, barcodes, on="barcode_id")
    matrix = pd.merge(matrix, features, on="feature_id")
    #matrix = matrix.drop(['barcode_id','feature_id'],axis=1)
    return matrix


def get_value_matrix(matrix):
    v_t = matrix['value'].str.split('_',n=1,expand=True)
    matrix['value']=v_t[0]
    matrix['value_type']=v_t[1]
    return matrix


def get_tetramer_matrix(matrix):
    matrix = matrix.sort_values(by=['count'],ascending=False)
    tetramers = matrix.drop_duplicates('barcode')
    tetramers = tetramers[tetramers['count']>5]
    tetramers['tetramer'] = tetramers['value']
    tetramers =  tetramers[['barcode_id','tetramer']]
    matrix = pd.merge(matrix, tetramers, on="barcode_id")
    return matrix[['barcode','value','tetramer','count']]

def get_barcode_tetramer(matrix):
    matrix = matrix.sort_values(by=['count'],ascending=False)
    tetramers = matrix.drop_duplicates('barcode')
    tetramers['tetramer'] = tetramers['value']
    tetramers =  tetramers[['barcode','tetramer','count']]
    return tetramers

def prep_annot_data(data,chain='None'):
    data_c = data[data['is_cell']==True].reset_index(drop=True)
    if chain=='TRB':
        data_c = data_c[data_c['chain'] == "TRB"]
    if chain=='TRA':
        data_c = data_c[data_c['chain'] == "TRA"]
    
    data_c = data_c[-data_c['v_gene'].isna()]
    data_c = data_c[-data_c['j_gene'].isna()]
    
    data_c = data_c[-data_c['cdr3'].isna()]
    #data_no_10x = data_no_10x[-data_no_10x['antigen.epitope'].isna()]
    data_c = data_c[-data_c['v_gene'].str.contains(',')]
    data_c = data_c[-data_c['j_gene'].str.contains(',')]
    data_c = data_c[-data_c['cdr3'].str.contains('\*')]
    data_c = data_c[data_c['cdr3']!='None']
    data_c = data_c[data_c['v_gene']!='None']
    data_c = data_c[data_c['j_gene']!='None']
    return data_c

def merge_anot_matrix(data_c,matrix):
    data_c = pd.merge(data_c,matrix, on='barcode')
    #data_c = data_c[['cdr3','count','tetramer','value','barcode']]
    #data_c = data_c[data_c['cdr3']!='None']
    return data_c


## Normalization
def sum_counts(data, cell_col_name, count_col_name, cell):
    return sum(data[data[cell_col_name]==cell][count_col_name])

def norm_logp_with_total_average(data, cell_id_col, count_col):
    # total count for cell
    cell_all_counts = pd.DataFrame(data[cell_id_col].drop_duplicates())
    cell_all_counts['counts']= cell_all_counts[cell_id_col].apply(lambda x: sum_counts(data_c_all,cell_id_col,count_col,x))
    
    cell_all_counts_dict = pd.Series(cell_all_counts['counts'].values,index=cell_all_counts[cell_id_col]).to_dict()
    ## average total count for cell
    average_total_count = cell_all_counts['counts'].mean()
    ## normalize log(1 + count / (total count for cell) * (average total count for cell)) 
    data['count']= data.apply(lambda row : math.log(1 + row[count_col]/cell_all_counts_dict[row[cell_id_col]]*average_total_count),axis=1)
    return data
    
    
def norm_logp(data, count_col):
    data[count_col]= data[count_col].apply(lambda x : math.log1p(x))
    return data    

def pivot_data(data):
    data = data[['count','barcode','value']]
    data = data.pivot_table('count','barcode','value')
    data = data.fillna(0)
    return data

## analysis
def pca(data, n):
    
    if type(n) == int:
        x = StandardScaler().fit_transform(data)
        pca = PCA(n_components = n )
        principalComponents = pca.fit_transform(x)
        df_pca = pd.DataFrame(data = principalComponents)
        return df_pca, pca
    

def tsne(data,n,random_s,p):
    if type(n) == int:
        X_embedded = TSNE(n_components=n,init='pca',
                         random_state=random_s, perplexity=p).fit_transform(data)
        return X_embedded
    
def knn(n_neighbors,X_train,y_train,X_test):
    knn = neighbors.KNeighborsClassifier(n_neighbors)
    knn.fit(X_train, y_train)
    
    pred = knn.predict(X_test)
    return pred

def count_roc_auc(y_test_curv,y_pred_curv):
    n_classes = y_test_curv.shape[1]
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_curv[:,i],y_pred_curv[:,i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    return roc_auc

def cat_lable(labels):
    le = preprocessing.LabelEncoder()
    le.fit(labels)
    return le.transform(labels)

def pc_anova(data,pc_n, group):
    pc_anova = pd.DataFrame(columns=['pc','F','pvalue'])

    for n in range(0,pc_n):
        F,pval = f_classif(np.array(data[n]).reshape(-1,1),data[group])
        pc_anova = pc_anova.append({'pc': n ,'F': float(F),'pvalue':float(pval)},ignore_index=True)
    
    return pc_anova

def binominal_test(binom_df, cluster, group):
    binom_df['total_cluster'] = binom_df.groupby(cluster)[cluster].transform('count')
    binom_df['total_group'] = binom_df.groupby(group)[group].transform('count')
    binom_df['count_mached'] = binom_df.groupby([group,cluster])[group].transform('count')
    binom_df['fraction_mached'] = binom_df['count_mached']/binom_df['total_cluster']
    binom_df['fraction_mached_exp'] = binom_df['total_group']/len(binom_df.index)
    binom_df['p_value'] = binom_df.apply(lambda row: stats.binom_test(row['count_mached'], n=row['total_cluster'], p=row['fraction_mached_exp'], alternative='greater'),axis=1)
    binom_df = binom_df[binom_df['fraction_mached']>binom_df['fraction_mached_exp']]
    binom_df_cluster = binom_df[[group, cluster,'total_cluster','total_group','count_mached','fraction_mached','fraction_mached_exp','p_value']].drop_duplicates().sort_values('p_value')
    binom_df_cluster = binom_df_cluster.sort_values(['fraction_mached'],ascending=False)
    binom_df_cluster = binom_df_cluster.drop_duplicates('cluster',keep='first')
    return binom_df_cluster