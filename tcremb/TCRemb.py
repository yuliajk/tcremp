from pathlib import Path
import numpy as np
import pandas as pd
import math

import sys
sys.path.append("../")
sys.path.append("mirpy/")
from mir.common import parser, Repertoire, SegmentLibrary
from mir.distances import ClonotypeAligner, GermlineAligner
from mir.comparative import DenseMatcher

from collections import Counter
import time

import statistics
from scipy.spatial.distance import pdist, squareform, cdist
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize

import tcremb.data_proc as data_proc
import tcremb.ml_utils as ml_utils

import tcremb.motif_logo as motif_logo
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

class TCRemb:
    clonotype_id = 'cloneId'
    annotation_id = 'annotId'
    random_state = 7
    
    def __init__(self,run_name, input_data, data_id = None):
        self.clonotypes={}
        #self.clonotype_label_pairs = {}
        self.annot={}
        self.dists={}
        self.pca_clones={}
        self.pca={}
        self.pca_clone_label={}
        self.tsne={}
        self.tsne_clones={}
        self.clstr_labels ={}
        self.clsf_labels = {}
        
        self.tcr_columns = ['cdr3aa','v','j','chain']
        self.tcr_columns_paired = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ']}
        self.__rename_tcr_columns_paired = {'TRA':{'a_cdr3aa':'cdr3aa','TRAV':'v','TRAJ':'j','cloneId_TRA':'cloneId'},'TRB':{'b_cdr3aa':'cdr3aa','TRBV':'v','TRBJ':'j','cloneId_TRB':'cloneId'}}
        self.clonotype_id = 'cloneId'
        self.clonotyoe_label_id = 'pairId'
        self.input_id= 'inputId'
        self.annotation_id = 'annotId'
        self.clonotype_id_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
        #self.__prototypes_path = { 'TRA' :'mirpy/notebooks/assets/olga_humanTRA_3000.txt', 'TRB' : 'mirpy/notebooks/assets/olga_humanTRB_3000.txt'}
        self.__prototypes_path = { 'TRA' :'data/data_preped/olga_humanTRA.txt', 'TRB' : 'data/data_preped/olga_humanTRB.txt'}
        
        self.__n_components = 50
        self.__tsne_init = 'pca'
        self.__tsne_perplexity = 15
        self.__random_state = 7
        
        
        self.outputs_path = "tcremb_outputs/" + run_name + '/'
        Path(self.outputs_path).mkdir(parents=True, exist_ok=True)
        
        self.clonotypes_path = { 'TRA' : self.outputs_path + 'clonotypes_TRA.txt', 'TRB' : self.outputs_path + 'clonotypes_TRB.txt'}
        self.dists_res_path = {'TRA' : self.outputs_path + 'res_TRA.txt', 'TRB': self.outputs_path + 'res_TRB.txt'}
        
        self.data_id = data_id
        self.input_data = input_data.copy()
        self.input_data = self.__annot_id(self.input_data, self.input_id)

    def __annot_id(self, data, annotation_id_str):
        df = data.copy()
        df[annotation_id_str]=df.index
        return df
    
#    def __filter_segments(self, chain, df_clones,segments_path='mir/resources/segments.txt'):
#        segs = pd.read_csv(segments_path,sep='\t')
#        segs = segs[segs['organism']=='HomoSapiens']
#        segs_ids = list(segs['id'].drop_duplicates())
#        if (chain=='TRA') or (chain=='TRB'):
#            df_clones = df_clones[df_clones['v'].isin(segs_ids)]
#            df_clones = df_clones[df_clones['j'].isin(segs_ids)]
#        if chain=='TRA_TRB':
#            df_clones = df_clones[df_clones['TRAV'].isin(segs_ids)]
#            df_clones = df_clones[df_clones['TRAJ'].isin(segs_ids)]
#            df_clones = df_clones[df_clones['TRBV'].isin(segs_ids)]
#            df_clones = df_clones[df_clones['TRBV'].isin(segs_ids)]
#        df_clones = df_clones.reset_index(drop=True)
#        return df_clones

    def __assign_clones_ids(self, data):
        df = data.copy()
        df[self.clonotype_id]=df.groupby(self.tcr_columns,dropna=False).ngroup()
        return df
    
    def __assign_clones_ids_paired(self, data, chain):
        df = data.copy()
        if (chain=='TRA') or (chain=='TRB'):
            df[self.clonotype_id_dict['TRA_TRB'][chain]]=df.groupby(self.tcr_columns_paired[chain],dropna=False).ngroup()
        if chain=='TRA_TRB':
            df[self.clonotype_id]=df.groupby(self.tcr_columns_paired['TRA']+self.tcr_columns_paired['TRB'],dropna=False).ngroup()
        return df
    
#    def __assign_clone_label_pairs_ids(self, chain, label):
#        if chain=='TRA_TRB':
#            self.annot[chain][self.clonotyoe_label_id] = self.annot[chain].groupby([[label] + list(self.clonotype_id_dict[chain].values())],dropna=False).ngroup()
#            
#        else:
#            self.annot[chain][self.clonotyoe_label_id] = self.annot[chain].groupby([label,self.clonotype_id],dropna=False).ngroup()  
    
    
    #def __clonotypes_prep(self, clones_df,clonotypes_path, chain, tcr_columns, clonotype_id_str):
    def __clonotypes_prep(self, clones_df, chain, tcr_columns, clonotype_id_str):
        clonotypes = clones_df[clones_df['chain']==chain]
        clonotypes = data_proc.remove_asterisk(clonotypes, tcr_columns)
        clonotypes = data_proc.remove_backslash(clonotypes, tcr_columns)
        clonotypes = data_proc.filter_clones_data(clonotypes, tcr_columns)
        clonotypes = data_proc.filter_segments(clonotypes)
        
        clonotypes = clonotypes[tcr_columns + [clonotype_id_str]].drop_duplicates().reset_index(drop=True)
        #clonotypes.to_csv(clonotypes_path, sep='\t')
        return clonotypes

    
    def tcremb_clonotypes(self,chain):
        #data_tt = self.__filter_segments(chain, self.input_data)
        data_tt = self.input_data.copy()
        if (chain=='TRA') or (chain=='TRB'):
            data_tt = self.__assign_clones_ids(data_tt)
            data_tt['clone_size'] = data_tt.groupby(self.clonotype_id)[self.input_id].transform('count')
            #self.clonotypes[chain] = self.__clonotypes_prep(data_tt, self.clonotypes_path[chain], chain, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain] = self.__clonotypes_prep(data_tt,  chain, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain].to_csv(self.clonotypes_path[chain], sep='\t')
            
            data_tt = data_tt[data_tt[self.clonotype_id].isin(self.clonotypes[chain][self.clonotype_id])]
            self.annot[chain] = self.__annot_id(data_tt[data_tt['chain']==chain].reset_index(drop=True), self.annotation_id)
            
        elif chain=='TRA_TRB':
            chain_1 = 'TRA'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_chain_1 = data_tt.copy()
            data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            data_chain_1['chain']=chain_1
            #self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, chain_1, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain_1].to_csv(self.clonotypes_path[chain_1], sep='\t')
            self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)                  
            
            chain_1 = 'TRB'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_chain_1 = data_tt.copy()
            data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            data_chain_1['chain']=chain_1
            #self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, chain_1, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain_1].to_csv(self.clonotypes_path[chain_1], sep='\t')
            self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)
            
            data_tt = self.__assign_clones_ids_paired(data_tt, chain)
            data_tt['clone_size'] = data_tt.groupby(self.clonotype_id)[self.input_id].transform('count')
            
            chain_1 = 'TRA'
            data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain_1][self.clonotype_id])]
            chain_1 = 'TRB'
            data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain_1][self.clonotype_id])]
            
            self.annot[chain] = self.__annot_id(data_tt.reset_index(drop=True), self.annotation_id)

        else:
            print('Error. Chain is incorrect. Must be TRA, TRB or TRA_TRB')
     
   
    def __data_parse_mirpy(self, chain, olga_human_path, clonotypes_path):
        lib = SegmentLibrary.load_default(genes=chain)
        db = Repertoire.load(parser=parser.OlgaParser(), path=olga_human_path)
    
        pars = parser.ClonotypeTableParser(lib=lib,
                                  )
        data_parse = pars.parse(source=clonotypes_path)
        data_parse = [x for x in data_parse if len(x.cdr3aa) in range(7, 23)]
        print(data_parse[0:10])
        return lib, db, data_parse

    def __mir_launch(self, chain, lib, db, data_parse, nproc, chunk_sz):
        aligner = ClonotypeAligner.from_library(lib=lib)
        matcher = DenseMatcher(db, aligner)
        
        start = time.time()
        res = matcher.match_to_df(data_parse, nproc=nproc)
        end = time.time()
        print(np.shape(res))
        print(end - start)
        return res

    def tcremb_dists_count(self, chain, nproc= None, chunk_sz=100):
        lib, db, data_parse = self.__data_parse_mirpy(chain, self.__prototypes_path[chain],self.clonotypes_path[chain])
        res = self.__mir_launch(chain, lib, db, data_parse, nproc, chunk_sz)
        res.to_csv(self.dists_res_path[chain], sep='\t', index = False)
    
    def __mir_results_proc(self, chain, res_path_chain, clonotypes_path_chain, clonotype_id_str):
        res_df = pd.read_csv(res_path_chain,sep='\t')
        clonotypes = pd.read_csv(clonotypes_path_chain, sep='\t')
        clonotypes['id']=clonotypes.index
        res_df = res_df.merge(clonotypes[['id',clonotype_id_str]], on='id').drop('id',axis=1)
        return res_df
            
    def tcremb_palette(labels_list):
        self.palette = ml_utils.make_custom_palette(labels_list)
    
    def tcremb_dists(self, chain):
        self.dists[chain] = self.__mir_results_proc(chain, self.dists_res_path[chain], self.clonotypes_path[chain], self.clonotype_id)
        self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.dists[chain][self.clonotype_id]))].reset_index(drop=True)
        #if len(self.clonotype_label_pairs.values()) != 0:
        #    self.clonotype_label_pairs[chain] = self.clonotype_label_pairs[chain][self.clonotype_label_pairs[chain][self.clonotype_id].isin(list(self.dists[chain][self.clonotype_id]))].reset_index(drop=True)

    def tcremb_pca(self, chain):
        if (chain == 'TRA') or (chain == 'TRB'):
            self.pca_clones[chain] = ml_utils.pca_proc(self.dists[chain], self.clonotype_id, self.__n_components)
            self.pca[chain] = self.pca_clones[chain].merge(self.annot[chain][[self.clonotype_id,self.annotation_id]]).drop(self.clonotype_id, axis=1, errors =
                                                                                                                               'ignore').sort_values(self.annotation_id).reset_index(drop=True)
            self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.pca_clones[chain][self.clonotype_id]))].reset_index(drop=True)
        
        elif chain=='TRA_TRB':
            dists_data = self.annot[chain][[self.annotation_id, self.clonotype_id] +list(self.clonotype_id_dict[chain].values())]
            annot_clones = dists_data[[self.annotation_id, self.clonotype_id]]
            
            
            
            dists_data = dists_data.drop(self.annotation_id,axis=1).drop_duplicates().reset_index(drop=True)
            chain_1 = 'TRA'
            dists_data_a = dists_data.merge(self.dists[chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            dists_data_a = dists_data_a.drop(list(self.clonotype_id_dict[chain].values()),axis=1)
            chain_1 = 'TRB'
            dists_data_b = dists_data.merge(self.dists[chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            dists_data_b = dists_data_b.drop(list(self.clonotype_id_dict[chain].values()),axis=1)
            
            dists_data = dists_data_a.merge(dists_data_b, on = self.clonotype_id)
            
            self.pca_clones[chain] = ml_utils.pca_proc(dists_data, self.clonotype_id, self.__n_components)
            self.pca[chain] = self.pca_clones[chain].merge(annot_clones).drop(self.clonotype_id,axis=1)
            
            self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.pca_clones[chain][self.clonotype_id]))].reset_index(drop=True)
            
            #dists_data = self.annot[chain][[self.annotation_id] + list(self.clonotype_id_dict[chain].values())]
            #chain_1 = 'TRA'
            #dists_data = dists_data.merge(self.dists[chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            #chain_1 = 'TRB'
            #dists_data = dists_data.merge(self.dists[chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            #dists_data = dists_data.drop(self.clonotype_id_dict[chain].values(), axis=1, errors ='ignore')
            
            #self.pca[chain] = ml_utils.pca_proc(dists_data, self.annotation_id, self.__n_components).sort_values(self.annotation_id).reset_index(drop=True)
            #self.pca_clones[chain] = self.pca[chain].merge(self.annot[chain][[self.annotation_id] + list(self.clonotype_id_dict[chain].values())]).drop(self.annotation_id,axis=1).drop_duplicates()
            ##self.annot[chain] = self.annot[chain][self.annot[chain][self.annotation_id].isin(list(dists_data[self.annotation_id]))].reset_index(drop=True)
            
            
    def tcremb_tsne(self,chain):
        self.tsne[chain] = ml_utils.tsne_proc(self.pca[chain] , self.annotation_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        self.tsne_clones[chain] = ml_utils.tsne_proc(self.pca_clones[chain] , self.clonotype_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        
class TCRemb_vdjdb(TCRemb):
    def __init__(self,run_name, input_data, data_id = None):
        TCRemb.__init__(self,run_name, input_data, data_id)
        self.clonotypes_pred={}
        self.vdjdb_path = 'data/VDJdb_pred/vdjdb_data_with_cloneId.txt'
        self.vdjdb_clonotypes_path = {'TRA': 'data/VDJdb_pred/clonotypes_TRA_VDJdb.txt', 'TRB' : 'data/VDJdb_pred/clonotypes_TRB_VDJdb.txt'}
        self.vdjdb_dists_res_path = {'TRA': 'data/VDJdb_pred/res_TRA_VDJdb.txt', 'TRB' : 'data/VDJdb_pred/res_TRB_VDJdb.txt'}
        self.vdjdb = pd.read_csv(self.vdjdb_path,sep='\t')
        self.vdjdb_clonotype_id = 'vdjdb_cloneId'
        self.data_type = 'data_type'
        self.label = 'antigen.epitope_freq'
        #self.label = 'antigen.epitope'
        
        self.dists_dubs = {}
        
    def __assign_clone_ids_with_vdjdb(self, df):
        df = df.merge(df.merge(self.vdjdb[self.tcr_columns + [self.vdjdb_clonotype_id]].drop_duplicates())[[self.input_id,self.vdjdb_clonotype_id]], how='left').reset_index(drop=True)
        t = df[df[self.vdjdb_clonotype_id].isna()]
        t = self._TCRemb__assign_clones_ids(t)
        t[self.clonotype_id] = t[self.clonotype_id] + (max(self.vdjdb[self.vdjdb_clonotype_id] +1))
        
        df = df.merge(t[[self.input_id,self.clonotype_id]],how='left')
        df[self.clonotype_id] = df.apply(lambda x: x[self.vdjdb_clonotype_id] if math.isnan(x[self.clonotype_id])  else x[self.clonotype_id],axis=1)
        return df
    
    
    def tcremb_clonotypes(self,chain):
        df = self.input_data.copy()
        if (chain=='TRA') or (chain=='TRB'):
            df = self.__assign_clone_ids_with_vdjdb(df)
            
            df['clone_size'] = df.groupby(self.clonotype_id)[self.input_id].transform('count')
            
            self.clonotypes_pred[chain] = self._TCRemb__clonotypes_prep(df, chain, self.tcr_columns, self.clonotype_id)
            self.clonotypes_pred[chain].to_csv(self.clonotypes_path[chain], sep='\t')
            
            self.vdjdb[self.clonotype_id] = self.vdjdb[self.vdjdb_clonotype_id]
            self.vdjdb[self.data_type]='train'
            df[self.label]='no'
            df[self.data_type]='pred'
            df = pd.concat([self.vdjdb, df])
            
            self.clonotypes[chain] = self._TCRemb__clonotypes_prep(df, chain, self.tcr_columns, self.clonotype_id)
            
            df = df[df[self.clonotype_id].isin(self.clonotypes[chain][self.clonotype_id])]
            self.annot[chain] = self._TCRemb__annot_id(df[df['chain']==chain].reset_index(drop=True), self.annotation_id)
            
        
        ## TO DO!!
        
        
        #elif chain=='TRA_TRB':
        #    chain_1 = 'TRA'
        #    data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
        #    data_chain_1 = data_tt.copy()
        #    data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
        #    data_chain_1['chain']=chain_1
        #    self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
        #    self.clonotypes[chain_1].to_csv(self.clonotypes_path[chain_1], sep='\t')
        #    self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)                  
        #    
        #    chain_1 = 'TRB'
        #    data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
        #    data_chain_1 = data_tt.copy()
        #    data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
        #    data_chain_1['chain']=chain_1
        #    self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
        #    self.clonotypes[chain_1].to_csv(self.clonotypes_path[chain_1], sep='\t')
        #    self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)
        #    
        #    data_tt = self.__assign_clones_ids_paired(data_tt, chain)
        #    data_tt['clone_size'] = data_tt.groupby(self.clonotype_id)[self.input_id].transform('count')
        #    
        #    chain_1 = 'TRA'
        #    data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain_1][self.clonotype_id])]
        #    chain_1 = 'TRB'
        #    data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain_1][self.clonotype_id])]
        #    
        #    self.annot[chain] = self.__annot_id(data_tt.reset_index(drop=True), self.annotation_id)

        else:
            print('Error. Chain is incorrect. Must be TRA, TRB or TRA_TRB')
            
    
    def tcremb_dists(self, chain):
        
        dists_vdjdb = self._TCRemb__mir_results_proc(chain, self.vdjdb_dists_res_path[chain], self.vdjdb_clonotypes_path[chain], self.clonotype_id)
        
        dists_pred = self._TCRemb__mir_results_proc(chain, self.dists_res_path[chain], self.clonotypes_path[chain], self.clonotype_id)
        self.dists_dubs[chain] = pd.concat([dists_vdjdb,dists_pred])
        self.dists[chain] = self.dists_dubs[chain].drop_duplicates(self.clonotype_id).reset_index(drop=True)
        #self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.dists[chain][self.clonotype_id]))].reset_index(drop=True)

    
        
class TCRemb_clustering():
    def __init__(self, model_name, threshold=0.7):
        self.annotation_id = 'annotId'
        self.model_name = model_name
        self.threshold = threshold
        self.clstr_labels = {}
        self.binom_res = {}
        self.clstr_metrics = {}
        self.model = {}
        #self.n_clusters = {}
        self.silhouette_n_clusters = {}
    
    def clstr(self, chain, data, label_cl=None, model='kmeans'):
        
        annot_clones = data.annot[chain][[data.clonotype_id, data.annotation_id]]
        X_data = data.pca_clones[chain].copy()
        
        check_between = False
        model_name = None
        if (model == 'kmeans'):
            self.silhouette_clusters(X_data.drop(data.clonotype_id,axis=1), chain)
            model_name = model
            model =  KMeans(n_clusters=self.silhouette_n_clusters[chain], random_state=7)
            check_between = True
            
        if model=='dbscan':
            model_name = model
            model = DBSCAN(eps=3, min_samples=2)        
        
        clstr_labels, self.model[chain] = ml_utils.clstr_model(model, X_data , data.clonotype_id)
        self.clstr_labels[chain] = clstr_labels.merge(annot_clones).drop(data.clonotype_id, axis=1)
        
        if label_cl is not None:
            self.binom_res[chain] = ml_utils.binominal_test(pd.merge(self.clstr_labels[chain],data.annot[chain]), 'cluster', label_cl, self.threshold)
            
            if model_name=='dbscan':
                self.binom_res[chain]['is_cluster'] = self.binom_res[chain].apply(lambda x: x.is_cluster if x.cluster != -1 else 0,axis=1)
        
            #self.binom_res[chain]['is_cluster']= self.binom_res[chain]['total_cluster'].apply(lambda x: 1 if x>1 else 0)
            #self.binom_res[chain]['enriched_clstr'] =self.binom_res[chain].apply(lambda x:1 
            #                                                                     if (x.fraction_matched>=self.threshold)
            #                                                                     and (x.is_cluster==1) else 0,axis=1)
            self.clstr_metrics[chain] = ml_utils.clstr_metrics(data.annot[chain][label_cl],self.clstr_labels[chain]['cluster'])
        
            self.clstr_labels[chain] = pd.merge(self.clstr_labels[chain], self.binom_res[chain].rename({label_cl:'label_cluster'},axis=1)
                                            , on='cluster',how='left').sort_values(self.annotation_id).reset_index(drop=True)
            
            if check_between:
                self.__is_cluster_by_between_metric(data, chain)
            
            self.purity = ml_utils.count_clstr_purity(self.binom_res[chain])
            #self.mean_fraction_matched = statistics.mean(self.binom_res[chain][self.binom_res[chain]['is_cluster']==1]['fraction_matched'])
            #self.median_fraction_matched = statistics.median(self.binom_res[chain][self.binom_res[chain]['is_cluster']==1]['fraction_matched'])
            #print(f'mean fraction_matched only clusters: {self.mean_fraction_matched}')
            #print(f'median fraction_matched only clusters: {self.median_fraction_matched}')
            print(f'purity:{self.purity}')
        
        if chain == 'TRA_TRB':
            self.clstr_labels[chain] = pd.merge(self.clstr_labels[chain], data.annot[chain][[self.annotation_id] 
                                                                                   + list(data.tcr_columns_paired['TRA'].values()) 
                                                                                  + list(data.tcr_columns_paired['TRB'].values())])
        else: 
            self.clstr_labels[chain] = pd.merge(self.clstr_labels[chain],data.annot[chain][[self.annotation_id] + data.tcr_columns])
        
        
        
    def silhouette_clusters(self, X_data, chain):
        X_data = X_data.drop(self.annotation_id, axis=1, errors = 'ignore')
    
        data_len = len(X_data)
    
        range_n_clusters = [round(data_len*0.005),round(data_len*0.01), round(data_len*0.05) , round(data_len*0.1), round(data_len*0.15)
                        ,round(data_len*0.2), round(data_len*0.25), round(data_len*0.3)
                        , round(data_len*0.4), round(data_len*0.5), round(data_len*0.6), round(data_len*0.7), round(data_len*0.8), round(data_len*0.9)]
    
        silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)
    
        n_to_try = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]
    
        range_n_clusters = [round(n_to_try*0.7), round(n_to_try*0.8), round(n_to_try*0.9), n_to_try, round(n_to_try*1.1) , round(n_to_try*1.2) , round(n_to_try*1.3)]
        silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)

        self.silhouette_n_clusters[chain] = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]
        
        
    def __is_cluster_by_between_metric(self, data,chain):
        df = self.clstr_labels[chain].copy()
        df = df[(df['total_cluster']>1)]
        all_cl = list(df['cluster'].drop_duplicates())
    
        df_pca_cl = df[['cluster','annotId']].merge(data.pca[chain])
        df_centrs = pd.DataFrame(self.model[chain].cluster_centers_)
    
        min_between_cl = min(pdist(df_centrs))
        mean_between_cl = statistics.mean(df_centrs)

        all_dists = []
        for i in all_cl:
            df_cl = df_pca_cl[df_pca_cl['cluster']==i].drop(['cluster','annotId'],axis=1)
            centr = df_centrs.iloc[i]
            cl_dists = []
            for _,j in df_cl.iterrows():
                cl_dists.append(cdist(np.array([centr]),np.array([j]),'euclidean')[0][0])
            all_dists.append(statistics.mean(cl_dists))

        cl_mean_dist = statistics.mean(all_dists)
    
        cl_dists = pd.DataFrame([all_cl,all_dists]).T.rename({0:'cluster',1:'dist'},axis=1)

        br = self.binom_res[chain].copy()
        br = br.merge(cl_dists.merge(self.binom_res[chain]),how='left')

        try_dist = statistics.mean([cl_mean_dist,min_between_cl])
        #try_dist = min_between_cl
        #try_dist = statistics.mean([min_between_cl, mean_between_cl])
        br['pred_enriched_min_between']=br.apply(lambda r: 1 if (r.is_cluster == 1) & (r.dist<=try_dist) else 0, axis = 1)

        pred_enriched_between_cls = list(br[br['pred_enriched_min_between']==1]['cluster'])

        self.clstr_labels[chain]['is_cluster_between'] = self.clstr_labels[chain]['cluster'].apply(lambda x: 1 if x in pred_enriched_between_cls else 0)

        
    def __plot_logo(self, clstr_data, chain, c, list_ax):
        
        cluster_df = clstr_data[clstr_data['cluster']==c]
        lengs = cluster_df['cdr3aa_len'].drop_duplicates()
        fr_matched = cluster_df['fraction_matched'].drop_duplicates().reset_index(drop=True)[0]
        epi = cluster_df['label_cluster'].drop_duplicates().reset_index(drop=True)[0]
        total_cl = cluster_df['total_cluster'].drop_duplicates().reset_index(drop=True)[0]
        alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']
        for l in lengs:
            seqs = cluster_df[cluster_df['cdr3aa_len']==l]['cdr3aa'].reset_index(drop=True)
            if len(seqs) > 4:
                motif_logo.plot_amino_logo(seqs, 'title',ax = list_ax[0])
                list_ax[0].set_title(f"{chain}. Cluster: {c} {epi}\nFraction matched:{round(fr_matched,2)}\nCount of cdr3aa: {len(seqs)}")
        plot_v_j = clstr_data[clstr_data['cluster']==c]

        
        plot_v_j = pd.DataFrame(plot_v_j.groupby('v')['cdr3aa'].count().reset_index())
        list_ax[1].pie(plot_v_j['cdr3aa'],labels=plot_v_j['v'])


        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j = pd.DataFrame(plot_v_j.groupby('j')['cdr3aa'].count().reset_index())
        list_ax[2].pie(plot_v_j['cdr3aa'],labels=plot_v_j['j'])       
    
    def __plot_logo_paired(self, clstr_data,chain, c, list_ax):
    
        cluster_df = clstr_data[clstr_data['cluster']==c]
        fr_matched = cluster_df['fraction_matched'].drop_duplicates()['fraction_matched'].reset_index(drop=True)[0]
        epi = cluster_df['label_cluster'].drop_duplicates().reset_index(drop=True)[0]
        total_cl = cluster_df['total_cluster'].drop_duplicates().reset_index(drop=True)[0]
        lengs = cluster_df['a_cdr3aa_len'].drop_duplicates()
        alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']
        for l in lengs:
            seqs = cluster_df[cluster_df['a_cdr3aa_len']==l]['a_cdr3aa']
            if len(seqs) > 4:
                freq = np.zeros((len(alphabet), l))
                for pos in range(l):
                    for s in seqs:
                        freq[alphabet.index(s[pos]), pos] +=1
                freq_res = freq/freq.sum(axis=0, keepdims=True)
                motif = pd.DataFrame(freq_res,index=alphabet)
                motif_logo.plot_amino_logo(motif, 'title',ax = list_ax[0])
                list_ax[0].set_title(f"{chain}. Cluster: {c} {epi}\nFraction matched:{round(fr_matched,2)}\nCount of cdr3aa: {len(seqs)}")
    
        cluster_df = clstr_data[clstr_data['cluster']==c]
        lengs = cluster_df['b_cdr3aa_len'].drop_duplicates()
        alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']
        for l in lengs:
            seqs = cluster_df[cluster_df['b_cdr3aa_len']==l]['b_cdr3aa']
            if len(seqs) > 4:
                freq = np.zeros((len(alphabet), l))
                #print("cdr3aa length: {} \nCount of cdr3aa: {}".format(l,len(seqs)))
                for pos in range(l):
                    for s in seqs:
                        freq[alphabet.index(s[pos]), pos] +=1
                freq_res = freq/freq.sum(axis=0, keepdims=True)
                motif = pd.DataFrame(freq_res,index=alphabet)
                motif_logo.plot_amino_logo(motif, 'title',ax = list_ax[1])
                list_ax[1].set_title(f"TRB. cdr3aa length: {l}")
    
        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('TRAV').transform('size')
        plot_v_j = plot_v_j.sort_values('TRAV')
        sns.histplot(plot_v_j[['count','TRAV']],y='TRAV',ax=list_ax[2])
        #list_ax[2].set_title('V genes count in cluster, cluster ' + str(c))

        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('TRAJ').transform('size')
        plot_v_j = plot_v_j.sort_values('TRAJ')
        sns.histplot(plot_v_j[['count','TRAJ']],y='TRAJ',ax=list_ax[3])
        #list_ax[3].set_title('J genes count in cluster, cluster ' + str(c))

        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('TRBV').transform('size')
        plot_v_j = plot_v_j.sort_values('TRBV')
        sns.histplot(plot_v_j[['count','TRBV']],y='TRBV',ax=list_ax[4])

        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('TRBJ').transform('size')
        plot_v_j = plot_v_j.sort_values('TRBJ')
        sns.histplot(plot_v_j[['count','TRBJ']],y='TRBJ',ax=list_ax[5])
        
        
    def clstrs_motif(self, data, chain, n_head_clstrs, sfig=None):
        if (chain=='TRA') or (chain=='TRB'):
            plt_clusters = list(self.binom_res[chain].sort_values('p_value').head(n_head_clstrs)['cluster'])
            clstr_data = pd.merge(self.clstr_labels[chain],data.annot[chain])
        
            clstr_data['cdr3aa_len'] = clstr_data['cdr3aa'].apply(len)
        
            n_rows = math.ceil(len(plt_clusters)/4)
            
            if sfig is None:
                sfig = plt.figure(figsize=(28,6*n_rows))
                outer_grid = gridspec.GridSpec(n_rows, 4,figure=sfig)        
            else:
                outer_grid = gridspec.GridSpec(n_rows, 4,figure=sfig)


            gs = []
            ax_list = []
            for i in range(len(plt_clusters)):
                cl_row = math.ceil((i+1)/4)-1
                cl_col = i%4
                gs.append(outer_grid[cl_row,cl_col].subgridspec(6, 2))
                ax_list.append([sfig.add_subplot(gs[i][:3, :2])
                                ,sfig.add_subplot(gs[i][4:6, 0])
                                ,sfig.add_subplot(gs[i][4:6, 1])])
    
                self.__plot_logo(clstr_data, chain, plt_clusters[i], ax_list[i])

        elif chain=='TRA_TRB':
            plt_clusters = list(self.binom_res[chain].sort_values('p_value').head(n_head_clstrs)['cluster'])
            clstr_data = pd.merge(self.clstr_labels[chain],data.annot[chain])
    
            clstr_data['a_cdr3aa_len'] = clstr_data['a_cdr3aa'].apply(len)
            clstr_data['b_cdr3aa_len'] = clstr_data['b_cdr3aa'].apply(len)
    
            n_rows = math.ceil(len(plt_clusters)/4)

            if sfig is None:
                sfig = plt.figure(figsize=(28,8*n_rows))
                outer_grid = gridspec.GridSpec(n_rows, 4,figure=sfig)        
            else:
                outer_grid = gridspec.GridSpec(n_rows, 4,figure=sfig)
        
            gs = []
            ax_list = []
            for i in range(len(plt_clusters)):
                cl_row = math.ceil((i+1)/4)-1
                cl_col = i%4
                gs.append(outer_grid[cl_row,cl_col].subgridspec(14, 2))
                ax_list.append([sfig.add_subplot(gs[i][1:4, :2])
                                ,sfig.add_subplot(gs[i][5:8, :2])
                                ,sfig.add_subplot(gs[i][9:11, 0])
                                ,sfig.add_subplot(gs[i][9:11, 1])
                                ,sfig.add_subplot(gs[i][12:14, 0])
                                ,sfig.add_subplot(gs[i][12:14, 1])])
    
                self.__plot_logo_paired(clstr_data, chain, plt_clusters[i], ax_list[i])
        else:
            print('chain is incorrect')

class TCRemb_clustering_pred(TCRemb_clustering):
    def __init__(self, model_name, threshold=0.7):
        TCRemb_clustering.__init__(self,model_name)
        self.annotation_id = 'annotId'
        self.data_type = 'data_type'
        #self.n_clusters = {}
        self.binom_res_train = {}
        self.clstr_metrics_train = {}
    
    def clstr_pred(self, chain, data, label_cl, model=None):
        
        self.clstr(chain, data, label_cl, model)
        
        y_data_train = data.annot[chain][data.annot[chain][self.data_type]=='train'][label_cl]
                
        clstrs_labels_data_annot_train = pd.merge(self.clstr_labels[chain], data.annot[chain][data.annot[chain][self.data_type]=='train'])
        
        self.binom_res_train[chain] = ml_utils.binominal_test(clstrs_labels_data_annot_train, 'cluster', label_cl, self.threshold)
        
        #self.binom_res_train[chain]['is_cluster_train']= self.binom_res_train[chain]['total_cluster'].apply(lambda x: 1 if x>1 else 0)
        #self.binom_res_train[chain]['enriched_clstr_train'] =self.binom_res_train[chain].apply(lambda x:1 
        #                                                                     if (x.fraction_matched>=self.threshold)
        #                                                                     and (x.is_cluster==1) else 0,axis=1)
        self.clstr_metrics[chain] = ml_utils.clstr_metrics(clstrs_labels_data_annot_train[label_cl], clstrs_labels_data_annot_train['cluster'])
        
        self.clstr_labels[chain] = pd.merge(self.clstr_labels[chain], self.binom_res_train[chain].rename({label_cl:'label_cluster_train'
                                                                                                          ,'fraction_matched':'fraction_matched_train'
                                                                                                          ,'p_value':'p_value_train'
                                                                                                          ,'is_cluster':'is_cluster_train'
                                                                                                         ,'enriched_clstr':'enriched_clstr_train'},axis=1)[['cluster',
                                                                                                                                              'label_cluster_train',
                                                                                                                                              'fraction_matched_train',
                                                                                                                                              'p_value_train',
                                                                                                                                              'is_cluster_train',
                                                                                                                                              'enriched_clstr_train']]
                                            , on='cluster',how='left').sort_values(self.annotation_id).reset_index(drop=True)
        
        if chain == 'TRA_TRB':
            self.clstr_labels[chain] = self.clstr_labels[chain].merge(data.annot[chain][[self.annotation_id, self.data_type] 
                                                                                   + list(self.tcr_columns_paired['TRA'].values()) 
                                                                                  + list(self.tcr_columns_paired['TRB'].values())])
        else: 
            self.clstr_labels[chain] = self.clstr_labels[chain].merge(data.annot[chain][[self.annotation_id, self.data_type] + data.tcr_columns])
        
        self.train_purity = ml_utils.count_clstr_purity(self.binom_res_train[chain].rename({'is_cluster_train':'is_cluster'},axis=1))
        #self.train_mean_fraction_matched = statistics.mean(self.binom_res_train[chain][self.binom_res_train[chain]['is_cluster_train']==1]['fraction_matched'])
        #self.train_median_fraction_matched = statistics.median(self.binom_res_train[chain][self.binom_res_train[chain]['is_cluster_train']==1]['fraction_matched'])
        #print(f'train mean fraction_matched only clusters: {self.train_mean_fraction_matched}')
        #print(f'train median fraction_matched only clusters: {self.train_median_fraction_matched}')
        print(f'train purity:{self.train_purity}')    


class TCRemb_clf():
    def __init__(self, model_name):
        self.model_name = model_name
        self.annotation_id = 'annotId' 
        self.label = {}
        self.__max_depth = 20
        self.y_train = {}
        self.y_test = {}
        self.X_train = {}
        self.X_test = {}
        self.model = {}
        self.test_pred = {}
        self.test_pred_proba = {}
        self.clsf_metrics = {}
        self.roc_auc_proba_df = {}
        self.roc_auc_proba = {}
    
    def clf(self, chain, data, label_cl, model=None , test_size = 0.3):
        if model is None:
            model =  RandomForestClassifier(max_depth=self.__max_depth, random_state=7)
        self.model[chain] = model
        self.label[chain] = label_cl
        y_data = data.annot[chain][label_cl]
        X_data = data.pca[chain].drop(self.annotation_id, axis=1, errors = 'ignore')
        self.__clf_model(chain, y_data, X_data, test_size)
        
    def __clf_model(self, chain, y_data, X_data, test_size):
        
        self.y_train[chain], self.y_test[chain] , self.X_train[chain], self.X_test[chain]  = train_test_split(y_data, X_data ,test_size=test_size)
        
        self.model[chain].fit(self.X_train[chain], self.y_train[chain])
        
        self.test_pred[chain] = self.model[chain].predict(self.X_test[chain])
        self.test_pred_proba[chain] = self.model[chain].predict_proba(self.X_test[chain])
        
        self.clsf_metrics[chain]={}
        self.clsf_metrics[chain]['train_acc'] = self.model[chain].score(self.X_train[chain],self.y_train[chain])
        self.clsf_metrics[chain]['test_acc'] = self.model[chain].score(self.X_test[chain],self.y_test[chain])
        self.clsf_metrics[chain]['f1_macro'] = ml_utils.clsf_metrics(self.test_pred[chain],self.y_test[chain],'macro')
        self.clsf_metrics[chain]['f1_weighted'] = ml_utils.clsf_metrics(self.test_pred[chain],self.y_test[chain],'weighted')
        #self.clsf_metrics[chain]['cross_val_score'] = cross_val_score(self.model[chain], X_data, y_data, cv=5, scoring="accuracy")
        #self.clsf_metrics[chain]['cross_val_f1_weighted'] = cross_val_score(self.model[chain], X_data, y_data, cv=5, scoring="f1_weighted")
        
        print(f"Train accuracy: {self.clsf_metrics[chain]['train_acc']}")
        print(f"Test accuracy: {self.clsf_metrics[chain]['test_acc']}")
        print(f"macro: {self.clsf_metrics[chain]['f1_macro']}")
        print(f"weighted: {self.clsf_metrics[chain]['f1_weighted']}")
        
        #print(f"cross_val_score: {self.clsf_metrics[chain]['cross_val_score']}")
        #print(f"cross_val_score f1_weighted: {self.clsf_metrics[chain]['cross_val_f1_weighted']}")
        
        
    def roc_auc(self, chain, ax=None, show_legend=True, custom_palette=None):
        classes_list = list(self.model[chain].classes_)
        y_test_curv = label_binarize(self.y_test[chain], classes=classes_list)
        y_pred_curv = label_binarize(self.test_pred[chain], classes=classes_list)
        
        roc_auc_proba = ml_utils.roc_auc_count(y_test_curv,self.test_pred_proba[chain])
        self.roc_auc_proba[chain]=roc_auc_proba
        self.roc_auc_proba_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc_proba.values()})
        
        ml_utils.plot_roccurve_multi(classes_list, y_test_curv, self.test_pred_proba[chain],f'ROC curves, {chain}', ax, custom_palette=custom_palette, test_acc=self.clsf_metrics[chain]['test_acc'], f1_weighted=self.clsf_metrics[chain]['f1_weighted'],show_legend=show_legend)
            

class TCRemb_clf_bind(TCRemb_clf):  
    def __init__(self, model_name):
        TCRemb_clf.__init__(self,model_name)
        self.annotation_id = 'annotId'
        self.class_col = 'bind'
        self.pn_pairs = {}
        
    def __positive_negative_samle(self, chain, data, label_cl):
        positive_pairs = data.clonotype_label_pairs[chain].copy()
        n = len(positive_pairs)*2
        negative_pairs = ml_utils.generate_negative_pairs(positive_pairs, n , data.clonotype_id, label_cl)
        
        index_start = max(data.clonotype_label_pairs[chain][data.clonotyoe_label_id]) + 1
        negative_pairs = negative_pairs.reset_index().rename({'index':data.clonotyoe_label_id},axis=1)
        
        positive_pairs[self.class_col] = 1
        negative_pairs[self.class_col] = 0
        
        return pd.concat([positive_pairs, negative_pairs]).reset_index(drop=True)
        
    
    def clf(self, chain, data, label_cl, model=None , labels_list=None, test_size = 0.3):
        data_tt = self.__positive_negative_samle(chain, data, label_cl)
        
        if labels_list is not None:
            data_tt = data_tt[data_tt[label_cl].isin(labels_list)].reset_index(drop=True)
        
        if model is None:
            model =  RandomForestClassifier(max_depth=self.__max_depth, random_state=7)
        self.model[chain] = model
        self.label[chain] = label_cl
        
        y_data = data_tt[self.class_col]
        
        if (chain=='TRA') or (chain=='TRB'):
            X_data = pd.merge(data_tt[data.clonotype_id], data.pca_clones[chain], how='left').drop(data.clonotype_id, axis=1, errors = 'ignore')
        
        #if chain == 'TRA_TRB':
        #    X_data = 
        
        self.y_data = y_data
        self.pca_clones_labels = X_data
        self.pn_pairs[chain] = data_tt
        
        self._TCRemb_clf__clf_model(chain, y_data, X_data, test_size)        

        
    def roc_auc(self, chain, labels_list, ax=None, show_legend=True):
        test_and_pred = pd.DataFrame(self.y_test[chain])
        test_and_pred['pred'] = self.test_pred[chain]
        
        pred_pairs = pd.concat([self.pn_pairs[chain].drop('bind',axis=1),test_and_pred],axis=1)
        pred_pairs = pred_pairs[-pred_pairs['pred'].isna()]
        roc_auc_list = []
        
        for l in labels_list:
            pred_pairs_l = pred_pairs[pred_pairs[label]==l]
            y_test = pred_pairs_l[clf.class_col]
            y_test_pred = pred_pairs_l['pred']
            roc_auc_list.append({'label': l, 'auc': ml_utils.roc_auc_count_binary(y_test, y_test_pred)})
        
        self.roc_auc_proba_df[chain] = pd.DataFrame(roc_auc_list)
    
    

class TCRemb_clf_pred(TCRemb_clf):
    def __init__(self, model_name):
        TCRemb_clf.__init__(self,model_name)
        self.annotation_id = 'annotId'
        self.X_pred = {}
        self.pred = {}
        self.y_pred_real = {}
        self.pred_proba = {}
        self.clsf_metrics_pred = {}
        self.roc_auc_proba_pred_df = {}
        
    def clf(self, chain, data, label_cl, model = None, test_size = 0.3):
        if model is None:
            model =  RandomForestClassifier(max_depth=self._TCRemb_clf__max_depth, random_state=7)
        self.model[chain] = model
        y_data = data.annot[chain][data.annot[chain]['data_type']=='train'][label_cl]
        X_data = data.pca[chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        X_data = X_data[X_data['data_type']=='train'].drop([self.annotation_id, 'data_type'], axis=1, errors = 'ignore')
        self._TCRemb_clf__clf_model(chain, y_data, X_data, test_size)
        
    def clf_pred(self, chain, data, label_cl = None):
        X_pred = data.pca[chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        self.X_pred[chain] = X_pred[X_pred['data_type']=='pred'].drop([self.annotation_id, 'data_type'], axis=1, errors = 'ignore') 
        self.pred[chain] = self.model[chain].predict(self.X_pred[chain])
        self.pred_proba[chain] = self.model[chain].predict_proba(self.X_pred[chain])
        
        self.clsf_metrics_pred[chain]={}
        
        if label_cl is None:
            self.y_pred_real[chain]=None
        else: 
            self.y_pred_real[chain] = data.annot[chain][data.annot[chain]['data_type']=='pred'][label_cl]
            self.clsf_metrics_pred[chain]['pred_acc'] = self.model[chain].score(self.X_pred[chain],self.y_pred_real[chain])
            self.clsf_metrics_pred[chain]['f1_macro'] = ml_utils.clsf_metrics(self.pred[chain],self.y_pred_real[chain],'macro')
            self.clsf_metrics_pred[chain]['f1_weighted'] = ml_utils.clsf_metrics(self.pred[chain],self.y_pred_real[chain],'weighted')
            
            print(f"Prediction accuracy: {self.clsf_metrics_pred[chain]['pred_acc']}")
            print(f"macro: {self.clsf_metrics[chain]['f1_macro']}")
            print(f"weighted: {self.clsf_metrics[chain]['f1_weighted']}")
        
        
    def roc_auc_pred(self, chain, ax=None, show_legend=True, custom_palette=None):
        classes_list = list(self.model[chain].classes_)
        y_real_curv = label_binarize(self.y_pred_real[chain], classes=classes_list)
        y_pred_curv = label_binarize(self.pred[chain], classes=classes_list)
        
        roc_auc_proba = ml_utils.roc_auc_count(y_real_curv,self.pred_proba[chain])
        self.roc_auc_proba_pred_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc_proba.values()})
        
        ml_utils.plot_roccurve_multi(classes_list, y_real_curv, self.pred_proba[chain],f'ROC curves, {chain}', ax, custom_palette=custom_palette, show_legend=show_legend)
        
        #roc_auc = ml_utils.roc_auc_count(y_real_curv,y_pred_curv)
        #roc_auc_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc.values()})