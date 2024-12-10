#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 20:37:47 2020

@author: yisupeng
"""

import json
import csv
import numpy as np
import matplotlib.pyplot as plt
import os

#%%

#data_source = 'PRIDE'
#data_source = 'HeLa'
data_source = 'NIST'

if data_source == 'PRIDE':
    species_list = [
            'A.thaliana',
#            'C.elegans',
            'D.melanogaster',
            'E.coli',
            'H.sapiens2',
            'H.sapiens3',
            'M.musculus',
            'M.musculus2',
            'M.musculus3',
#            'S.cerevisiae',
            'S.cerevisiae2',
            'S.cerevisiae3',
        ]
    result_dir = 'test_search/est_results_pride/'
    data_source = 'PRIDE'
    figure_dir = 'figures/'

elif data_source == 'HeLa':
    species_list = [
            'HeLa01ng',
            'HeLa01ng.2',
            'HeLa01ng.3',
            'HeLa1ng',
            'HeLa1ng.2',
            'HeLa1ng.3',
            'HeLa10ng',
            'HeLa10ng.2',
            'HeLa10ng.3',
            'HeLa50ng',
            'HeLa50ng.2',
            'HeLa50ng.3',
            'HeLa100ng',
            'HeLa100ng.2',
            'HeLa100ng.3',
        ]
    psm_dir = 'test_search/pride/'
    result_dir = 'test_search/est_results_hela/'
    data_source = ''

elif data_source == 'NIST':
    species_list = [
            'c_elegans_tda',
            'h_sapiens_tda',
            'm_musculus_tda',
            's_cerevisiae_tda',
#            'drosophila',
#            'e_coli',
#            'human',
#            'mouse',
#            'yeast',
        ]
    psm_dir = 'data/nist/tsv_result_tda/'
    result_dir = 'data/nist/fdr_result_tda/'
    data_source = 'NIST'

method_list = [
        '_1s2ca',
        '_1s2c',
        '_2s3ci',
    ]

#psm_dir = 'test_search/pride/'
#result_dir = 'test_search/est_results_pride/'

#psm_dir = 'test_search/nist/'
#result_dir = 'test_search/est_results_nist/'

#%%
def get_tda_fdr(species):
#    species_dir = 'test_search/est_results/'+species#+'/'
    
    matches = open(psm_dir+'%s_d.tsv'%species)
    matches_csv = csv.DictReader(matches, delimiter='\t')
    
    evalue_field = 'SpecEValue'

    
    s1_tda = []
    fdr_tda = []
    flag = False
    thres = None
    prevspec = None
    for row in matches_csv:
        
        if row['Protein'][0:3] == 'REV':
            continue
        
        ind = int(row['SpecID'].split('=')[1])
        if prevspec == ind:
            continue
        prevspec = ind
        
        qval = float(row['QValue'])
        
        s = -np.log(float(row[evalue_field]))

        if not flag and qval >= 0.01:
            thres = s
            flag = True
        
        if qval == 1:
            continue
            
        s1_tda.append(s)
        fdr_tda.append(qval)
    
    print(len(fdr_tda))
    obj = {'ds': species, 'algo': 'TDA', 's1': s1_tda, 'fdr': fdr_tda, 't1p': thres}
    return obj

#%%
tda_fdrs = []
for species in species_list:
    obj = get_tda_fdr(species)
    tda_fdrs.append(obj)
#    plt.plot(obj['s1'], obj['fdr'])
    plt.plot(obj['fdr'])
#    print(obj['t1p'])
#    break

#%%
json_dir = result_dir+'json/'
if not os.path.exists(json_dir):
    os.mkdir(json_dir)
json.dump(tda_fdrs, open(result_dir+'json/tda_fdr.json', 'w'))
