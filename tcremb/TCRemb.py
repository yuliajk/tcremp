from pathlib import Path
import numpy as np
import pandas as pd
import math

import sys
#sys.path.append("../")

from collections import Counter
import time

from time import gmtime, strftime

import statistics
#from scipy.spatial.distance import pdist, squareform, cdist
from sklearn.model_selection import train_test_split

import tcremb.data_proc as data_proc
import tcremb.ml_utils as ml_utils

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

#sys.path.append("../mirpy/")
sys.path.append("../mirpy/mirpy/")
#from mir.common import parser, Repertoire, SegmentLibrary
from mir.common.repertoire import Repertoire
from mir.common.segments import SegmentLibrary
from mir.common import parser
from mir.distances import ClonotypeAligner, GermlineAligner
from mir.comparative import DenseMatcher

class TCRemb:
    clonotype_id = 'cloneId'
    annotation_id = 'annotId'
    random_state = 7
    
    def __init__(self,run_name, input_data, data_id = None, prototypes_path=None, n=None, species='HomoSapiens', prototypes_chain='TRA_TRB'):
        self.__prototypes_path_subsets = {'HomoSapiens': { 'TRA' :'data/data_preped/olga_humanTRA.txt', 'TRB' : 'data/data_preped/olga_humanTRB.txt'}}
        self.run_name = run_name
        self.species = species
        self.clonotypes={} ## extracted clonotypes
        #self.clonotype_label_pairs = {}
        self.annot_input={} ## raw input
        self.annot={} ## processed input table (cleaned clonotypes, added annotation id and clonotype id)
        self.dists={} ## dists data for clonotypes with clonotype id
        self.annot_dists = {} ## annotation id with dists data
        self.pca_clones={}
        self.pca={}
        self.pca_clone_label={}
        self.pca_ad={'clones':{},'all':{}}
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
        #self.annotation_id = 'annotId'
        self.annotation_id = 'id' ## index
        self.clonotype_id_dict = {'TRA': 'cloneId','TRB': 'cloneId','TRA_TRB': {'TRA':'cloneId_TRA', 'TRB':'cloneId_TRB'}}
        
        self.prototypes_path = self.__prototypes_path_subsets[species]
        #self.__prototypes_path = prototypes_path
        
        self.__n_components = 50
        
        self.__tsne_init = 'pca'
        self.__tsne_perplexity = 15
        self.__random_state = 7
        self.time_dict = {}
        
        
        self.outputs_path = "tcremb_outputs/" + run_name + '/'
        Path(self.outputs_path).mkdir(parents=True, exist_ok=True)
        
        self.clonotypes_path = { 'TRA' : self.outputs_path + 'clonotypes_TRA.txt', 'TRB' : self.outputs_path + 'clonotypes_TRB.txt',
                               'TRA_TRB': {'TRA' : self.outputs_path + 'clonotypes_paired_TRA.txt', 'TRB' : self.outputs_path + 'clonotypes_paired_TRB.txt'}}
        self.dists_res_path = {'TRA' : self.outputs_path + 'res_TRA.txt', 'TRB': self.outputs_path + 'res_TRB.txt',
                              'TRA_TRB':{'TRA' : self.outputs_path + 'res_paired_TRA.txt', 'TRB': self.outputs_path + 'res_paired_TRB.txt'}}
        
        self.data_id = data_id
        self.raw_input_data = input_data.copy()
        #050624self.input_data = input_data.copy()
        self.check_proc_input_data() #050624
        self.input_data = self.__annot_id(self.input_data, self.input_id)
        
        if prototypes_path:
            self.prototypes_path= { 'TRA' : self.outputs_path + 'prptotypesTRA.txt', 'TRB' : self.outputs_path + 'prptotypesTRB.txt'}
            self.prototypes_prep(prototypes_path)
        
        if n:
            new_prototypes_path = { 'TRA' : self.outputs_path + f'prptotypesTRA_{n}.txt', 'TRB' : self.outputs_path + f'prptotypesTRB_{n}.txt'}
            try:
                if prototypes_chain=='TRA_TRB':
                    self.prototypes_n(n, self.prototypes_path['TRA'], new_prototypes_path['TRA'])
                    self.prototypes_n(n, self.prototypes_path['TRB'], new_prototypes_path['TRB'])
                else:
                    self.prototypes_n(n, self.prototypes_path[prototypes_chain], new_prototypes_path[prototypes_chain])
            except ValueError:
                print('n is greater than number of clonotypes in prototypes file')
            self.prototypes_path= new_prototypes_path
            if n*2<50:
                self.__n_components =n*2
        else:
            n=3000
        
        self.dist_cols_dist={'TRA':[f'a_{xs}_{x}' for xs in range(n) for x in ['v','j','cdr3']],
               'TRB':[f'b_{xs}_{x}' for xs in range(n) for x in ['v','j','cdr3']]}

    def check_proc_input_data(self):
        data_proc.check_columns(self.raw_input_data, self.tcr_columns_paired['TRA'] + self.tcr_columns_paired['TRB'])
        self.input_data = data_proc.clean_at_least_cdr3a_or_cdr3b(self.raw_input_data, self.tcr_columns_paired['TRA'][0], self.tcr_columns_paired['TRB'][0], self.outputs_path)
        self.input_data =  data_proc.remove_asterisk(self.input_data, self.tcr_columns_paired['TRA']) 
        self.input_data =  data_proc.remove_asterisk(self.input_data, self.tcr_columns_paired['TRB'])
        self.input_data = data_proc.remove_backslash(self.input_data, self.tcr_columns_paired['TRA'])
        self.input_data = data_proc.remove_backslash(self.input_data, self.tcr_columns_paired['TRB'])
        self.input_data = data_proc.add_allele(self.input_data,self.tcr_columns_paired['TRA'])
        self.input_data = data_proc.add_allele(self.input_data,self.tcr_columns_paired['TRB']) 
        
    
    def prototypes_prep(self, input_file):
        pass
    
    def prototypes_n(self, n, old_path, new_path):
        #pd.read_csv(self.prototypes_path['TRA'],sep='\t',header=None).sample_n(n).reset_index(drop=True).to_csv('')
        pd.read_csv(old_path,sep='\t',index_col=0).sample(n).reset_index(drop=True).to_csv(new_path, sep='\t')
    
    def __annot_id(self, data, annotation_id_str):
        df = data.copy()
        df[annotation_id_str]=df.index
        return df

    def __assign_clones_ids(self, data, chain):
        df = data.copy()
        df[self.clonotype_id]=df.groupby(self.tcr_columns_paired[chain],dropna=False).ngroup()
        return df
    
    def __assign_clones_ids_paired(self, data, chain):
        df = data.copy()
        if (chain=='TRA') or (chain=='TRB'):
            df[self.clonotype_id_dict['TRA_TRB'][chain]]=df.groupby(self.tcr_columns_paired[chain],dropna=False).ngroup()
        if chain=='TRA_TRB':
            df[self.clonotype_id]=df.groupby(self.tcr_columns_paired['TRA']+self.tcr_columns_paired['TRB'],dropna=False).ngroup()
        return df
   
    
    
    #def __clonotypes_prep(self, clones_df,clonotypes_path, chain, tcr_columns, clonotype_id_str):
    def __clonotypes_prep_old(self, clones_df, chain, tcr_columns, clonotype_id_str):
        clonotypes = clones_df[clones_df['chain']==chain]
        clonotypes = clones_df.copy()
        clonotypes = data_proc.remove_asterisk(clonotypes, tcr_columns)
        clonotypes = data_proc.remove_backslash(clonotypes, tcr_columns)
        clonotypes = data_proc.filter_clones_data(clonotypes, tcr_columns)
        clonotypes = data_proc.filter_segments(clonotypes, segments_path='../mirpy/mirpy/mir/resources/segments.txt', organism=self.species)
        
        clonotypes = clonotypes[tcr_columns + [clonotype_id_str]].drop_duplicates().reset_index(drop=True)
        #clonotypes.to_csv(clonotypes_path, sep='\t')
        return clonotypes

    def __clonotypes_prep(self, clones_df, chain):
        clonotypes = clones_df.copy()
        clonotypes = clonotypes.rename(self.__rename_tcr_columns_paired[chain],axis=1)
        clonotypes['chain']=chain       
        clonotypes = clonotypes[self.tcr_columns + [self.clonotype_id]].drop_duplicates().reset_index(drop=True)
        clonotypes['d']='.' ## 070624
        return clonotypes
    
    def __clonotypes_data_clean(self, data, chain):
        df = data.copy()
        df = data_proc.filter_clones_data(df, self.tcr_columns_paired[chain], file_dir=self.outputs_path)
        #segments_path = mir.__file__.replace('__init__.py','resources/segments.txt')
        segments_path = '../mirpy/mirpy/mir/resources/segments.txt'
        df = data_proc.filter_segments(df, segments_path=segments_path, v = self.tcr_columns_paired[chain][1], j = self.tcr_columns_paired[chain][2], organism=self.species, file_dir=self.outputs_path)
        return df
    
    def tcremb_clonotypes(self,chain, unique_clonotypes=False):
        #self.time_dict[chain] = {}
        print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        
        #data_tt = self.__filter_segments(chain, self.input_data)
        data_tt = self.input_data.copy()
        if (chain=='TRA') or (chain=='TRB'):
            data_tt = data_tt[~data_tt[self.tcr_columns_paired[chain][0]].isna()].reset_index(drop=True)
            data_tt = self.__clonotypes_data_clean(data_tt, chain) #050624
            
            data_tt = self.__assign_clones_ids(data_tt, chain)
            if unique_clonotypes:
                data_tt = data_tt.drop_duplicates(self.clonotype_id).reset_index(drop=True)
            data_tt['clone_size'] = data_tt.groupby(self.clonotype_id)[self.input_id].transform('count')
            #050624data_chain_1 = data_tt.rename(self.__rename_tcr_columns_paired[chain],axis=1)
            #050624data_chain_1['chain']=chain
            #self.clonotypes[chain] = self.__clonotypes_prep(data_tt, self.clonotypes_path[chain], chain, self.tcr_columns, self.clonotype_id)
            #050624self.clonotypes[chain] = self.__clonotypes_prep(data_chain_1,  chain, self.tcr_columns, self.clonotype_id)
            self.clonotypes[chain] = self.__clonotypes_prep(data_tt, chain) #050624
            self.clonotypes[chain].to_csv(self.clonotypes_path[chain], sep='\t')
            
            data_tt = data_tt[data_tt[self.clonotype_id].isin(self.clonotypes[chain][self.clonotype_id])]
            self.annot_input[chain] = self.__annot_id(data_tt, self.annotation_id)
            
        elif chain=='TRA_TRB':
            data_tt = data_tt[~data_tt[self.tcr_columns_paired['TRA'][0]].isna()].reset_index(drop=True)
            data_tt = data_tt[~data_tt[self.tcr_columns_paired['TRB'][0]].isna()].reset_index(drop=True)
            self.clonotypes[chain] = {}
            
            chain_1 = 'TRA'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_tt = self.__clonotypes_data_clean(data_tt, chain_1) #050624
            #050624data_chain_1 = data_tt.copy()
            #050624data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            #050624data_chain_1['chain']=chain_1
            #self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
            #050624self.clonotypes[chain][chain_1] = self.__clonotypes_prep(data_chain_1, chain_1, self.tcr_columns, self.clonotype_id)
            
            self.clonotypes[chain][chain_1] = self.__clonotypes_prep(data_tt, chain_1) #050624
            self.clonotypes[chain][chain_1].to_csv(self.clonotypes_path[chain][chain_1], sep='\t')
            #self.annot_input[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)                  
            
            chain_1 = 'TRB'
            data_tt = self.__assign_clones_ids_paired(data_tt, chain_1)
            data_tt = self.__clonotypes_data_clean(data_tt, chain_1) #050624
            #050624data_chain_1 = data_tt.copy()
            #050624data_chain_1 = data_chain_1.rename(self.__rename_tcr_columns_paired[chain_1],axis=1)
            #050624data_chain_1['chain']=chain_1
            #self.clonotypes[chain_1] = self.__clonotypes_prep(data_chain_1, self.clonotypes_path[chain_1], chain_1, self.tcr_columns, self.clonotype_id)
            #050624self.clonotypes[chain][chain_1] = self.__clonotypes_prep(data_chain_1, chain_1, self.tcr_columns, self.clonotype_id)
            
            self.clonotypes[chain][chain_1] = self.__clonotypes_prep(data_tt, chain_1) #050624
            self.clonotypes[chain][chain_1].to_csv(self.clonotypes_path[chain][chain_1], sep='\t')
            #self.annot_input[chain_1] = self.__annot_id(data_chain_1.reset_index(drop=True), self.annotation_id)
            
            data_tt = self.__assign_clones_ids_paired(data_tt, chain)
            if unique_clonotypes:
                data_tt = data_tt.drop_duplicates(self.clonotype_id).reset_index(drop=True)
            data_tt['clone_size'] = data_tt.groupby(self.clonotype_id)[self.input_id].transform('count')
            
            chain_1 = 'TRA'
            data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain][chain_1][self.clonotype_id])]
            chain_1 = 'TRB'
            data_tt = data_tt[data_tt[self.clonotype_id_dict[chain_1]].isin(self.clonotypes[chain][chain_1][self.clonotype_id])]
            
            self.annot_input[chain] = self.__annot_id(data_tt.reset_index(drop=True), self.annotation_id)

        else:
            print('Error. Chain is incorrect. Must be TRA, TRB or TRA_TRB')
        
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        #self.time_dict[chain]['clonotypes'] = {end - start}
        print(f'Clonotypes extraction time: {end - start}')
     
   
    def __data_parse_mirpy(self, chain, olga_human_path, clonotypes_path):
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        lib = SegmentLibrary.load_default(genes=chain)
        db = Repertoire.load(parser=parser.OlgaParser(), path=olga_human_path)
    
        pars = parser.ClonotypeTableParser(lib=lib,
                                  )
        data_parse = pars.parse(source=clonotypes_path)
        data_parse = [x for x in data_parse if len(x.cdr3aa) in range(7, 23)]
        print(data_parse[0:10])
        
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        #self.time_dict[chain]['mir_parse'] = {end - start}
        print(f'parse data for mir: {end - start}')
        return lib, db, data_parse

    def __mir_launch(self, chain, lib, db, data_parse, nproc, chunk_sz):
        aligner = ClonotypeAligner.from_library(lib=lib)
        matcher = DenseMatcher(db, aligner)
        
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        res = matcher.match_to_df(data_parse, nproc=nproc)
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        print(np.shape(res))
        #self.time_dict[chain]['mir_launch'] = {end - start}
        print(f'Mir launch time: {end - start}')
        return res

    def tcremb_dists_count(self, chain, nproc= None, chunk_sz=100):
        if (chain=='TRA') or (chain=='TRB'):
            lib, db, data_parse = self.__data_parse_mirpy(chain, self.prototypes_path[chain],self.clonotypes_path[chain])
            res = self.__mir_launch(chain, lib, db, data_parse, nproc, chunk_sz)
            res.to_csv(self.dists_res_path[chain], sep='\t', index = False)
        elif chain=='TRA_TRB':
            chain_1 = 'TRA'
            lib, db, data_parse = self.__data_parse_mirpy(chain_1, self.prototypes_path[chain_1],self.clonotypes_path[chain][chain_1])
            res = self.__mir_launch(chain, lib, db, data_parse, nproc, chunk_sz)
            res.to_csv(self.dists_res_path[chain][chain_1], sep='\t', index = False)
            chain_1 = 'TRB'
            lib, db, data_parse = self.__data_parse_mirpy(chain_1, self.prototypes_path[chain_1],self.clonotypes_path[chain][chain_1])
            res = self.__mir_launch(chain, lib, db, data_parse, nproc, chunk_sz)
            res.to_csv(self.dists_res_path[chain][chain_1], sep='\t', index = False)
    
    def __mir_results_proc(self, chain, res_path_chain, clonotypes_path_chain, clonotype_id_str):
        res_df = pd.read_csv(res_path_chain,sep='\t')
        print(res_df.columns)
        res_df = res_df.set_axis(['id']+self.dist_cols_dist[chain],axis=1)
        print(res_df.columns)
        clonotypes = pd.read_csv(clonotypes_path_chain, sep='\t')
        clonotypes['id']=clonotypes.index
        res_df = res_df.merge(clonotypes[['id',clonotype_id_str]], on='id').drop('id',axis=1)
        return res_df
            
    def tcremb_palette(labels_list):
        self.palette = ml_utils.make_custom_palette(labels_list)
    
    def tcremb_dists(self, chain):
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        if (chain=='TRA') or (chain=='TRB'):
            self.dists[chain] = self.__mir_results_proc(chain, self.dists_res_path[chain], self.clonotypes_path[chain], self.clonotype_id)
            self.annot[chain] = self.annot_input[chain][self.annot_input[chain][self.clonotype_id].isin(list(self.dists[chain][self.clonotype_id]))].reset_index(drop=True)
            
            self.annot_dists[chain] = self.dists[chain].merge(self.annot[chain][[self.clonotype_id,self.annotation_id]]).drop(self.clonotype_id, axis=1, errors = 'ignore').sort_values(self.annotation_id).reset_index(drop=True) ##230524
            
            #if len(self.clonotype_label_pairs.values()) != 0:
            #    self.clonotype_label_pairs[chain] = self.clonotype_label_pairs[chain][self.clonotype_label_pairs[chain][self.clonotype_id].isin(list(self.dists[chain][self.clonotype_id]))].reset_index(drop=True)
        elif chain=='TRA_TRB':
            ## add clonotype_id to dists and filter not processed clones from annot
            self.dists[chain] = {}
            chain_1 = 'TRA'
            self.dists[chain][chain_1] = self.__mir_results_proc(chain_1, self.dists_res_path[chain][chain_1], self.clonotypes_path[chain][chain_1], self.clonotype_id)
            self.annot[chain] = self.annot_input[chain][self.annot_input[chain][self.clonotype_id_dict[chain][chain_1]].isin(list(self.dists[chain][chain_1][self.clonotype_id]))].reset_index(drop=True)
            
            chain_1 = 'TRB'
            self.dists[chain][chain_1] = self.__mir_results_proc(chain_1, self.dists_res_path[chain][chain_1], self.clonotypes_path[chain][chain_1], self.clonotype_id)
            self.annot[chain] = self.annot_input[chain][self.annot_input[chain][self.clonotype_id_dict[chain][chain_1]].isin(list(self.dists[chain][chain_1][self.clonotype_id]))].reset_index(drop=True)
            
            ## add annotation id to dists
            dists_data = self.annot[chain][[self.annotation_id, self.clonotype_id] +list(self.clonotype_id_dict[chain].values())] ##230524
            annot_clones = dists_data[[self.annotation_id, self.clonotype_id]] ##230524            
            
            dists_data = dists_data.drop(self.annotation_id,axis=1).drop_duplicates().reset_index(drop=True) ##230524
            
            chain_1='TRA' ##230524
            dists_data_a = dists_data.merge(self.dists[chain][chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1)) ##230524
            dists_data_a = dists_data_a.drop(list(self.clonotype_id_dict[chain].values()),axis=1) ##230524
            
            chain_1='TRB'
            dists_data_b = dists_data.merge(self.dists[chain][chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1)) ##230524
            dists_data_b = dists_data_b.drop(list(self.clonotype_id_dict[chain].values()),axis=1) ##230524
            
            self.dists[chain]['joined'] = dists_data_a.merge(dists_data_b, on = self.clonotype_id) ##230524
            self.annot_dists[chain] = self.dists[chain]['joined'].merge(annot_clones).drop(self.clonotype_id,axis=1) ##230524
            
            
        else: 
            print('Error. Chain is incorrect. Must be TRA, TRB or TRA_TRB')
        
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        #self.time_dict[chain]['dist_proc'] = {end - start}
        print(f'dist_proc: {end - start}')

    def tcremb_pca(self, chain, n_components = None):
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        if n_components is None:
            n_components = self.__n_components
        if (chain == 'TRA') or (chain == 'TRB'):
            self.pca_clones[chain] = ml_utils.pca_proc(self.dists[chain], self.clonotype_id, n_components)
            self.pca[chain] = self.pca_clones[chain].merge(self.annot[chain][[self.clonotype_id,self.annotation_id]]).drop(self.clonotype_id, axis=1, errors = 'ignore').sort_values(self.annotation_id).reset_index(drop=True)
            self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.pca_clones[chain][self.clonotype_id]))].reset_index(drop=True)
        
        elif chain=='TRA_TRB':
            dists_data = self.annot[chain][[self.annotation_id, self.clonotype_id] +list(self.clonotype_id_dict[chain].values())]
            annot_clones = dists_data[[self.annotation_id, self.clonotype_id]]
            
            
            
            dists_data = dists_data.drop(self.annotation_id,axis=1).drop_duplicates().reset_index(drop=True)
            chain_1 = 'TRA'
            self.pca_ad['clones'][chain_1] = ml_utils.pca_proc(self.dists[chain][chain_1], self.clonotype_id, round(n_components)).rename({self.clonotype_id:self.clonotype_id_dict[chain][chain_1]},axis=1)
            self.pca_ad['all'][chain_1] = self.pca_ad['clones'][chain_1].merge(self.annot[chain][[self.clonotype_id_dict[chain][chain_1],self.annotation_id]]).drop(self.clonotype_id_dict[chain][chain_1], axis=1, errors = 'ignore').sort_values(self.annotation_id).reset_index(drop=True)
            
            
            ##230524dists_data_a = dists_data.merge(self.dists[chain][chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            ##230524dists_data_a = dists_data_a.drop(list(self.clonotype_id_dict[chain].values()),axis=1)
            chain_1 = 'TRB'
            self.pca_ad['clones'][chain_1] = ml_utils.pca_proc(self.dists[chain][chain_1], self.clonotype_id, round(n_components)).rename({self.clonotype_id:self.clonotype_id_dict[chain][chain_1]},axis=1)
            self.pca_ad['all'][chain_1] = self.pca_ad['clones'][chain_1].merge(self.annot[chain][[self.clonotype_id_dict[chain][chain_1],self.annotation_id]]).drop(self.clonotype_id_dict[chain][chain_1], axis=1, errors = 'ignore').sort_values(self.annotation_id).reset_index(drop=True)
            
            ##230524dists_data_b = dists_data.merge(self.dists[chain][chain_1].rename({self.clonotype_id_dict[chain_1]:self.clonotype_id_dict[chain][chain_1]},axis=1))
            ##230524dists_data_b = dists_data_b.drop(list(self.clonotype_id_dict[chain].values()),axis=1)
            
            self.pca_ad['clones'][chain] = pd.merge(dists_data.merge(self.pca_ad['clones']['TRA']).drop(list(self.clonotype_id_dict[chain].values()),axis=1)
                                                    ,dists_data.merge(self.pca_ad['clones']['TRB']).drop(list(self.clonotype_id_dict[chain].values()),axis=1),on = self.clonotype_id)
            self.pca_ad['all'][chain] = self.pca_ad['clones'][chain].merge(annot_clones).drop(self.clonotype_id,axis=1)
            
            ##230524dists_data = dists_data_a.merge(dists_data_b, on = self.clonotype_id)
            
            ##230524self.dists[chain]['joined'] = dists_data.copy()
            
            self.pca_clones[chain] = ml_utils.pca_proc(self.dists[chain]['joined'], self.clonotype_id, self.__n_components)
            self.pca[chain] = self.pca_clones[chain].merge(annot_clones).drop(self.clonotype_id,axis=1)
            
            self.annot[chain] = self.annot[chain][self.annot[chain][self.clonotype_id].isin(list(self.pca_clones[chain][self.clonotype_id]))].reset_index(drop=True)
            
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        #self.time_dict[chain]['pca'] = {end - start}
        print(f'pca: {end - start}')    
            
    def tcremb_tsne(self,chain):
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        start = time.time()
        self.tsne[chain] = ml_utils.tsne_proc(self.pca[chain] , self.annotation_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        self.tsne_clones[chain] = ml_utils.tsne_proc(self.pca_clones[chain] , self.clonotype_id, self.__tsne_init, self.__random_state, self.__tsne_perplexity)
        
        #print(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        end = time.time()
        #self.time_dict[chain]['tsne'] = {end - start}
        print(f'tsne: {end - start}')    
