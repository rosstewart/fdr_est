#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:02:44 2020

@author: yisupeng
"""
import os
import tarfile

import csv
import numpy as np
import matplotlib.pyplot as plt

import json
from collections import defaultdict


#%%

#data_source = 'PRIDE'
data_source = 'HeLa'
#data_source = 'NIST'

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
    res_dir = 'test_search/est_results_pride/'
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
    res_dir = 'test_search/est_results_hela/'
    data_source = ''

elif data_source == 'NIST':
    species_list = [
            'c_elegans',
            'drosophila',
            'e_coli',
            'human',
            'mouse',
            'yeast',
        ]
    res_dir = 'test_search/est_results_nist/'
    data_source = 'NIST'

method_list = [
        '_1s2ca',
        '_1s2c',
        '_2s3ci',
    ]

#%%
#%%

method_map = {
        '1S2D gamma & gaussian': '_1s2ca',
        '1S2D skew normal': '_1s2c',
        '2S3D skew normal': '_2s3ci',
        'SNMax1': 'SNMax1',
        'SNMax2': 'SNMax2',
        'TDA': 'tda',
    }
algo_map = {v:k for k,v in method_map.items()}

algo_name_map = {
        'tda': 'TDA',
        '1S2D gamma & gaussian': 'GG',
        '1S2D skew normal': '1SMix',
        '2S3D skew normal': '2SMix',
        'SNMax1': 'SNMax1',
        'SNMax2': 'SNMax2',
    }

#%%

tdajson = json.load(open(res_dir+'json/tda_fdr.json'))
tdajson = {res['ds']:res for res in tdajson}


#%%
def plot_fdr_matches(fdr, nmatches):
    print(species, method)
#    print(estres.keys())
    
    fdr = np.array(fdr)
    nmatches = np.array(nmatches)
    flags = fdr < 0.1
    
    plt.plot(fdr[flags], nmatches[flags])

    
def plot_fdr_matches_for_species(species):
    species_dir = res_dir + species + '/'
    estjson = json.load(open(res_dir+'json/'+species+'.json'))
    estjson = {method_map[res['algo']]:res for res in estjson}
    legends = []
    for method in method_list:
        estres = estjson[method]
        fdr = estres['fdr']
        nmatches = list(reversed(range(1,len(estres['s1'])+1)))
        plot_fdr_matches(fdr, nmatches)
        legends.append(algo_name_map[estres['algo']])
    
    tdares = tdajson[species]
    nmatches = list(range(1,len(tdares['s1'])+1))
    plot_fdr_matches(tdares['fdr'], nmatches)
    legends.append('TDA')
    
    plt.legend(legends)
    plt.xlabel('FDR')
    plt.ylabel('# Matches')
    plt.title(data_source+' '+species)

    fdrcurv_dir = species_dir + 'fdrcurv/'
    
    plt.savefig(fdrcurv_dir + 'fdrmatches.png', dpi=320)
    
    plt.show()
    

    
#%%
for species in species_list:
    plot_fdr_matches_for_species(species)
#    break
        
    

