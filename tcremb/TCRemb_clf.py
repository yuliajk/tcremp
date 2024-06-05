from pathlib import Path
import numpy as np
import pandas as pd
import math

import sys
sys.path.append("../")
sys.path.append("mirpy/")

import time

import statistics
from scipy.spatial.distance import pdist, squareform, cdist
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.neighbors import NearestNeighbors

import tcremb.data_proc as data_proc
import tcremb.ml_utils as ml_utils


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
    
    def clf(self, chain, data, label_cl, model=None , test_size = 0.3, pca_sep_chains=False):        
        if model is None:
            model =  RandomForestClassifier(max_depth=self.__max_depth, random_state=7)
        self.model[chain] = model
        self.label[chain] = label_cl
        y_data = data.annot[chain][label_cl]
        if pca_sep_chains:
            X_data = data.pca_ad['all'][chain].drop(self.annotation_id, axis=1, errors = 'ignore')
        else:
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
        
    def clf(self, chain, data, label_cl, model = None, test_size = 0.3, pca_sep_chains=False):
        self.pca_sep_chains=pca_sep_chains
        if model is None:
            model =  RandomForestClassifier(max_depth=self._TCRemb_clf__max_depth, random_state=7)
        self.model[chain] = model
            
        y_data = data.annot[chain][data.annot[chain]['data_type']=='train'][label_cl]
        if self.pca_sep_chains:
            X_data = data.pca_ad['all'][chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        else:
            X_data = data.pca[chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
                
        X_data = X_data[X_data['data_type']=='train'].drop([self.annotation_id, 'data_type'], axis=1, errors = 'ignore')
        self._TCRemb_clf__clf_model(chain, y_data, X_data, test_size)
        
    def clf_pred(self, chain, data, label_cl = None):
        
        if self.pca_sep_chains:
            X_pred = data.pca_ad['all'][chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        else:
            X_pred = data.pca[chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        
        #X_pred = data.pca[chain].merge(data.annot[chain][[self.annotation_id, 'data_type']])
        
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