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

#%%

#species_list = [
#        'A.thaliana',
##        'C.elegans',
#        'D.melanogaster',
#        'E.coli',
#        'H.sapiens2',
#        'H.sapiens3',
#        'M.musculus',
#        'M.musculus2',
#        'M.musculus3',
##        'S.cerevisiae',
#        'S.cerevisiae2',
#        'S.cerevisiae3',
#    ]

species_list = [
        'HeLa01ng',
        'HeLa1ng',
        'HeLa10ng',
        'HeLa50ng',
        'HeLa100ng',
        'HeLa01ng.2',
        'HeLa1ng.2',
        'HeLa10ng.2',
        'HeLa50ng.2',
        'HeLa100ng.2',
        'HeLa01ng.3',
        'HeLa1ng.3',
        'HeLa10ng.3',
        'HeLa50ng.3',
        'HeLa100ng.3',
    ]

#species_list = [
#        'c_elegans',
#        'drosophila',
#        'e_coli',
#        'human',
#        'mouse'
#    ]

psm_dir = 'test_search/pride/'
#psm_dir = 'nist/'

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
json.dump(tda_fdrs, open('test_search/est_results/json/tda_fdr.json', 'w'))
#json.dump(tda_fdrs, open('test_search/est_results_nist/json/tda_fdr.json', 'w'))
