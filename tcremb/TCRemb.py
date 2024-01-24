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
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize

import tcremb.ml_utils as ml_utils

import tcremb.motif_logo
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

class TCRemb:
    clonotype_id = 'cloneId'
    annotation_id = 'annotId'
    random_state = 7
    
    def __init__(self,run_name):
        print(sys.path)
        self.clonotypes={}
        self.clonotype_label_pairs = {}
        self.annot={}
        self.dists={}
        self.pca_clones={}
        self.pca={}
        self.tsne={}
        self.tsne_clones={}
        self.clstr_labels ={}
        self.clsf_labels = {}
        
        self.__tcr_columns = ['cdr3aa','v','j','chain']
        self.__tcr_columns_paired = {'TRA':['a_cdr3aa','TRAV','TRAJ'],'TRB':['b_cdr3aa','TRBV','TRBJ']}
        self.__rename_tcr_columns_paired = {'TRA':{'a_cdr3aa':'cdr3aa','TRAV':'v','TRAJ':'j','cloneId_TRA':'cloneId'},'TRB':{'b_cdr3aa':'cdr3aa','TRBV':'v','TRBJ':'j','cloneId_TRB':'cloneId'}}
        self.__clonotype_id = 'cloneId'
        self.__data_id= 'data_id'
        self.__annotation_id = 'annotId'
        self.__clonotype_id_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
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

    def __annot_id(self, data, annotation_id_str):
        df = data.copy()
        df[annotation_id_str]=df.index
        return df
    
    def __filter_segments(self, chain, df_clones,segments_path='mir/resources/segments.txt'):
        segs = pd.read_csv(segments_path,sep='\t')
        segs = segs[segs['organism']=='HomoSapiens']
        segs_ids = list(segs['id'].drop_duplicates())
        if (chain=='TRA') or (chain=='TRB'):
            df_clones = df_clones[df_clones['v'].isin(segs_ids)]
            df_clones = df_clones[df_clones['j'].isin(segs_ids)]
        if chain=='TRA_TRB':
            df_clones = df_clones[df_clones['TRAV'].isin(segs_ids)]
            df_clones = df_clones[df_clones['TRAJ'].isin(segs_ids)]
            df_clones = df_clones[df_clones['TRBV'].isin(segs_ids)]
            df_clones = df_clones[df_clones['TRBV'].isin(segs_ids)]
        df_clones = df_clones.reset_index(drop=True)
        return df_clones

    def __assign_clones_ids(self, data):
        df = data.copy()
        df[self.__clonotype_id]=df.groupby(self.__tcr_columns,dropna=False).ngroup()
        return df
    
    def __assign_clones_ids_paired(self, data, chain):
        df = data.copy()
        df[self.__clonotype_id_dict['TRA_TRB'][chain]]=df.groupby(self.__tcr_columns_paired[chain],dropna=False).ngroup()
        return df
    

    def __clonotypes_prep(self, clones_df,clonotypes_path, chain, tcr_columns, clonotype_id_str):
        clonotypes = clones_df[clones_df['chain']==chain]
        clonotypes = clonotypes[tcr_columns + [clonotype_id_str]].drop_duplicates().reset_index(drop=True)
        clonotypes.to_csv(clonotypes_path, sep='\t')
        return clonotypes

    
    def tcremb_clonotypes(self,chain, data_preped):
        data_tt = self.__filter_segments(chain, data_preped)
        if (chain=='TRA') or (chain=='TRB'):
            data_tt = self.__assign_clones_ids(data_tt)
            self.clonotypes[chain] = self.__clonotypes_prep(data_tt, self.clonotypes_path[chain], chain, self.__tcr_columns, self.__clonotype_id)
            self.annot[chain] = self.__annot_id(data_tt[data_tt['chain']==chain].reset_index(drop=True), self.__annotation_id)
            
        elif chain=='TRA_TRB':
            chain_1 = 'TRA'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_chain_1 = data_tt.copy()
            data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            data_chain_1['chain']=chain_1
            self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.__tcr_columns, self.__clonotype_id)
            self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.__annotation_id)
            
            chain_1 = 'TRB'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_chain_1 = data_tt.copy()
            data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            data_chain_1['chain']=chain_1
            self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.__tcr_columns, self.__clonotype_id)
            self.annot[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.__annotation_id)
            
            self.annot[chain] = self.__annot_id(data_tt.reset_index(drop=True), self.__annotation_id)

        else:
            print('Error. Chain is incorrect. Must be TRA, TRB or TRA_TRB')
     
    
    def tcremb_clonotype_label_pairs(self, chain, label):
        if (chain=='TRA') or (chain=='TRB'):
            self.clonotype_label_pairs[chain] = self.annot[chain][[self.__clonotype_id,label]].drop_duplicates().rest_index(drop=True)
        elif chain=='TRA_TRB':
            self.clonotype_label_pairs[chain] = self.annot[chain][list(self.__clonotype_id_dict[chain].values()) + [label]].drop_duplicates().rest_index(drop=True)
        
#    def __data_parse_mirpy(self, chain, olga_human_path, clonotypes_path):
#        lib = SegmentLibrary.load_default(genes={chain})
#        db = Repertoire(parser.parse_olga_aa(olga_human_path, lib=lib))#, n=3000))
#        data_proc = parser.parse_from_file_simple(clonotypes_path, lib=lib, gene=chain,
#                               warn=False)
#        data_proc = [x for x in data_proc if len(x.cdr3aa) in range(7, 23)]
#        print(data_proc[0:10])
#        return lib, db, data_proc
   
    def __data_parse_mirpy(self, chain, olga_human_path, clonotypes_path):
        lib = SegmentLibrary.load_default(genes=chain)
        db = Repertoire.load(parser=parser.OlgaParser(), path=olga_human_path)
    
        pars = parser.ClonotypeTableParser(lib=lib,
                                  )
        data_proc = pars.parse(source=clonotypes_path)
        data_proc = [x for x in data_proc if len(x.cdr3aa) in range(7, 23)]
        print(data_proc[0:10])
        return lib, db, data_proc
    
#    def __mir_launch(self, chain, lib, db, data_proc):
#        valign = AlignGermline.from_seqs(lib.get_seqaas(gene=chain, stype='V'))
#        jalign = AlignGermline.from_seqs(lib.get_seqaas(gene=chain, stype='J'))
#        aligner = ClonotypeAligner(v_aligner=valign, j_aligner=jalign)
#        matcher = DenseMatch(db, aligner)
#    
#        start = time.time()
#        res = matcher.match_to_df(data_proc, 64, 384)
#        end = time.time()
#        print(np.shape(res))
#        print(end - start)
#        return res

    def __mir_launch(self, chain, lib, db, data_proc):
        aligner = ClonotypeAligner.from_library(lib=lib)
        matcher = DenseMatcher(db, aligner)
        
        start = time.time()
        res = matcher.match_to_df(data_proc)
        end = time.time()
        print(np.shape(res))
        print(end - start)
        return res

    def tcremb_dists_count(self, chain):
        lib, db, data_proc = self.__data_parse_mirpy(chain, self.__prototypes_path[chain],self.clonotypes_path[chain])
        res = self.__mir_launch(chain, lib, db, data_proc)
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
        self.dists[chain] = self.__mir_results_proc(chain, self.dists_res_path[chain], self.clonotypes_path[chain], self.__clonotype_id)
        self.annot[chain] = self.annot[chain][self.annot[chain][self.__clonotype_id].isin(list(self.dists[chain][self.__clonotype_id]))].reset_index(drop=True)

    def tcremb_pca(self, chain):
        if (chain == 'TRA') or (chain == 'TRB'):
            self.pca_clones[chain] = ml_utils.pca_proc(self.dists[chain], self.__clonotype_id, self.__n_components)
            self.pca[chain] = self.pca_clones[chain].merge(self.annot[chain][[self.__clonotype_id,self.__annotation_id]]).drop(self.__clonotype_id, axis=1, errors =
                                                                                                                               'ignore').sort_values(self.__annotation_id).reset_index(drop=True)
        elif chain=='TRA_TRB':
            dists_data = self.annot[chain][[self.__annotation_id] + list(self.__clonotype_id_dict[chain].values())]
            chain_1 = 'TRA'
            dists_data = dists_data.merge(self.dists[chain_1].rename({self.__clonotype_id_dict[chain_1]:self.__clonotype_id_dict[chain][chain_1]},axis=1))
            chain_1 = 'TRB'
            dists_data = dists_data.merge(self.dists[chain_1].rename({self.__clonotype_id_dict[chain_1]:self.__clonotype_id_dict[chain][chain_1]},axis=1))
            dists_data = dists_data.drop(self.__clonotype_id_dict[chain].values(), axis=1, errors ='ignore')
            
            self.pca[chain] = ml_utils.pca_proc(dists_data, self.__annotation_id, self.__n_components).sort_values(self.__annotation_id).reset_index(drop=True)
            self.annot[chain] = self.annot[chain][self.annot[chain][self.__annotation_id].isin(list(dists_data[self.__annotation_id]))].reset_index(drop=True)
            
            
    def tcremb_tsne(self,chain):
        self.tsne[chain] = ml_utils.tsne_proc(self.pca[chain] , self.__annotation_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        #self.tsne_clones[chain] = ml_utils.tsne_proc(self.pca_clones[chain] , self.__clonotype_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        
class TCRemb_clustering():
    def __init__(self, model_name):
        self.__annotation_id = 'annotId'
        self.model_name = model_name
        self.threshold = 0.6
        self.clstr_labels = {}
        self.binom_res = {}
        self.clstr_metrics = {}
        self.model = {}
        #self.n_clusters = {}
        self.silhouette_n_clusters = {}
    
    def clstr(self, chain, data, label_cl, model=None):
        y_data = data.annot[chain][label_cl]
        X_data = data.pca[chain]#.drop(self.__annotation_id, axis=1, errors = 'ignore')
        
        if model is None:
            self.silhouette_clusters(y_data, X_data, chain,label_cl)
            model =  KMeans(n_clusters=self.silhouette_n_clusters[chain], random_state=7)
        
        
        self.clstr_labels[chain], self.model[chain] = ml_utils.clstr_model(model, X_data , self.__annotation_id)
        #self.__clstr_metrics(chain, data, label_cl)
        
    #def __clstr_metrics(self, chain, data, label_cl):
        self.binom_res[chain] = ml_utils.binominal_test(pd.merge(self.clstr_labels[chain],data.annot[chain]), 'cluster', label_cl)
        
        self.binom_res[chain]['is_cluster']= self.binom_res[chain]['total_cluster'].apply(lambda x: 1 if x>1 else 0)
        self.binom_res[chain]['enriched_clstr'] =self.binom_res[chain].apply(lambda x:1 
                                                                             if (x.fraction_matched>=self.threshold)
                                                                             and (x.is_cluster==1) else 0,axis=1)
        self.clstr_metrics[chain] = ml_utils.clstr_metrics(data.annot[chain][label_cl],self.clstr_labels[chain]['cluster'])
        
        self.clstr_labels[chain] = pd.merge(self.clstr_labels[chain], self.binom_res[chain].rename({label_cl:'label_cluster'},axis=1)
                                            , on='cluster',how='left').sort_values(self.__annotation_id).reset_index(drop=True)
        
        self.purity = ml_utils.count_clstr_purity(self.binom_res[chain])
        self.mean_fraction_matched = statistics.mean(self.binom_res[chain][self.binom_res[chain]['is_cluster']==1]['fraction_matched'])
        self.median_fraction_matched = statistics.median(self.binom_res[chain][self.binom_res[chain]['is_cluster']==1]['fraction_matched'])
        print(f'mean fraction_matched only clusters: {self.mean_fraction_matched}')
        print(f'median fraction_matched only clusters: {self.median_fraction_matched}')
        print(f'purity:{self.purity}')
        
    def silhouette_clusters(self, y_data, X_data, chain, label):
        X_data = X_data.drop(self.__annotation_id, axis=1, errors = 'ignore')
    
        data_len = len(y_data)
    
        range_n_clusters = [round(data_len*0.005),round(data_len*0.01), round(data_len*0.05) , round(data_len*0.1), round(data_len*0.15)
                        ,round(data_len*0.2), round(data_len*0.25), round(data_len*0.3)
                        , round(data_len*0.4), round(data_len*0.5), round(data_len*0.6), round(data_len*0.7), round(data_len*0.8), round(data_len*0.9)]
    
        silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)
    
        n_to_try = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]
    
        range_n_clusters = [round(n_to_try*0.7), round(n_to_try*0.8), round(n_to_try*0.9), n_to_try, round(n_to_try*1.1) , round(n_to_try*1.2) , round(n_to_try*1.3)]
        silhouette_avg_scores = ml_utils.silhouette_avg_scores_kmeans(X_data,range_n_clusters)

        self.silhouette_n_clusters[chain] = list(silhouette_avg_scores.keys())[list(silhouette_avg_scores.values()).index(max(silhouette_avg_scores.values()))]
        
    def __plot_logo(self, clstr_data, chain, c, list_ax):
        
        cluster_df = clstr_data[clstr_data['cluster']==c]
        lengs = cluster_df['cdr3aa_len'].drop_duplicates()
        fr_matched = cluster_df['fraction_matched'].drop_duplicates().reset_index(drop=True)[0]
        epi = cluster_df['label_cluster'].drop_duplicates().reset_index(drop=True)[0]
        total_cl = cluster_df['total_cluster'].drop_duplicates().reset_index(drop=True)[0]
        alphabet = [aa for aa in 'ARNDCQEGHILKMFPSTWYVBZX-']
        for l in lengs:
            seqs = cluster_df[cluster_df['cdr3aa_len']==l]['cdr3aa']
            if len(seqs) > 4:
                freq = np.zeros((len(alphabet), l))
                for pos in range(l):
                    for s in seqs:
                        freq[alphabet.index(s[pos]), pos] +=1
                freq_res = freq/freq.sum(axis=0, keepdims=True)
                motif = pd.DataFrame(freq_res,index=alphabet)
                motif_logo.plot_amino_logo(motif, 'title',ax = list_ax[0])
                list_ax[0].set_title(f"{chain}. Cluster: {c} {epi}\nFraction matched:{round(fr_matched,2)}\nCount of cdr3aa: {len(seqs)}")
        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('v').transform('size')
        plot_v_j = plot_v_j.sort_values('v')
        sns.histplot(plot_v_j[['count','v']],y='v',ax=list_ax[1])
        #list_ax[2].set_title('V genes count in cluster, cluster ' + str(c))

        plot_v_j = clstr_data[clstr_data['cluster']==c]
        plot_v_j['count']=plot_v_j.groupby('j').transform('size')
        plot_v_j = plot_v_j.sort_values('j')
        sns.histplot(plot_v_j[['count','j']],y='j',ax=list_ax[2])
        #list_ax[3].set_title('J genes count in cluster, cluster ' + str(c))        
    
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
        
        
    def clstrs_motif(self, data, chain, n_head_clstrs):
        if (chain=='TRA') or (chain=='TRB'):
            plt_clusters = list(self.binom_res[chain].sort_values('p_value').head(n_head_clstrs)['cluster'])
            clstr_data = pd.merge(self.clstr_labels[chain],data.annot[chain])
        
            clstr_data['cdr3aa_len'] = clstr_data['cdr3aa'].apply(len)
        
            n_rows = math.ceil(len(plt_clusters)/4) +1
            fig = plt.figure(figsize=(8*n_rows,28))
            outer_grid = gridspec.GridSpec(n_rows, 4,figure=fig)
            gs = []
            ax_list = []
            for i in range(len(plt_clusters)):
                cl_row = math.ceil((i+1)/4)-1
                cl_col = i%4
                gs.append(outer_grid[cl_row,cl_col].subgridspec(7, 2))
                ax_list.append([fig.add_subplot(gs[i][1:4, :2])
                                ,fig.add_subplot(gs[i][5:7, 0])
                                ,fig.add_subplot(gs[i][5:7, 1])])
    
                self.__plot_logo(clstr_data, chain, plt_clusters[i], ax_list[i])

        elif chain=='TRA_TRB':
            plt_clusters = list(self.binom_res[chain].sort_values('p_value').head(n_head_clstrs)['cluster'])
            clstr_data = pd.merge(self.clstr_labels[chain],data.annot[chain])
    
            clstr_data['a_cdr3aa_len'] = clstr_data['a_cdr3aa'].apply(len)
            clstr_data['b_cdr3aa_len'] = clstr_data['b_cdr3aa'].apply(len)
    
            fig = plt.figure(figsize=(35,28))
            n_rows = math.ceil(len(plt_clusters)/4) +1
            outer_grid = gridspec.GridSpec(n_rows, 4,figure=fig)
            gs = []
            ax_list = []
            for i in range(len(plt_clusters)):
                cl_row = math.ceil((i+1)/4)-1
                cl_col = i%4
                gs.append(outer_grid[cl_row,cl_col].subgridspec(14, 2))
                ax_list.append([fig.add_subplot(gs[i][1:4, :2])
                                ,fig.add_subplot(gs[i][5:8, :2])
                                ,fig.add_subplot(gs[i][9:11, 0])
                                ,fig.add_subplot(gs[i][9:11, 1])
                                ,fig.add_subplot(gs[i][12:14, 0])
                                ,fig.add_subplot(gs[i][12:14, 1])])
    
                self.__plot_logo_paired(clstr_data, chain, plt_clusters[i], ax_list[i])
        else:
            print('chain is incorrect')

    

class TCRemb_clf():
    def __init__(self, model_name):
        self.model_name = model_name
        self.__annotation_id = 'annotId'
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
    
    def clf(self, chain, data, label_cl, model=None , test_size = 0.3):
        if model is None:
            model =  RandomForestClassifier(max_depth=self.__max_depth, random_state=7)
        self.model[chain] = model
        y_data = data.annot[chain][label_cl]
        X_data = data.pca[chain].drop(self.__annotation_id, axis=1, errors = 'ignore')
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
        
        
    def roc_auc(self, chain, ax=None, show_legend=True):
        classes_list = list(self.model[chain].classes_)
        y_test_curv = label_binarize(self.y_test[chain], classes=classes_list)
        y_pred_curv = label_binarize(self.test_pred[chain], classes=classes_list)
        
        roc_auc_proba = ml_utils.roc_auc_count(y_test_curv,self.test_pred_proba[chain])
        self.roc_auc_proba_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc_proba.values()})
        
        ml_utils.plot_roccurve_multi(classes_list, y_test_curv, self.test_pred_proba[chain],f'ROC curves, {chain}', ax, self.clsf_metrics[chain]['test_acc'], self.clsf_metrics[chain]['f1_weighted'],show_legend)
            
            
        
class TCRemb_clf_pred(TCRemb_clf):
    def __init__(self, model_name):
        TCRemb_clf.__init__(self,model_name)
        self.__annotation_id = 'annotId'
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
        X_data = data.pca[chain].merge(data.annot[chain][[self.__annotation_id, 'data_type']])
        X_data = X_data[X_data['data_type']=='train'].drop([self.__annotation_id, 'data_type'], axis=1, errors = 'ignore')
        self._TCRemb_clf__clf_model(chain, y_data, X_data, test_size)
        
    def clf_pred(self, chain, data, label_cl = None):
        X_pred = data.pca[chain].merge(data.annot[chain][[self.__annotation_id, 'data_type']])
        self.X_pred[chain] = X_pred[X_pred['data_type']=='pred'].drop([self.__annotation_id, 'data_type'], axis=1, errors = 'ignore') 
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
        
        
    def roc_auc_pred(self, chain, ax=None, show_legend=True):
        classes_list = list(self.model[chain].classes_)
        y_real_curv = label_binarize(self.y_pred_real[chain], classes=classes_list)
        y_pred_curv = label_binarize(self.pred[chain], classes=classes_list)
        
        roc_auc_proba = ml_utils.roc_auc_count(y_real_curv,self.pred_proba[chain])
        self.roc_auc_proba_pred_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc_proba.values()})
        
        ml_utils.plot_roccurve_multi(classes_list, y_real_curv, self.pred_proba[chain],f'ROC curves, {chain}', ax,show_legend)
        
        #roc_auc = ml_utils.roc_auc_count(y_real_curv,y_pred_curv)
        #roc_auc_df[chain] = pd.DataFrame({'class':classes_list ,'roc_auc':roc_auc.values()})