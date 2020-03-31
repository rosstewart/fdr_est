#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 23:37:43 2020

@author: yisupeng
"""

import csv
import matplotlib.pyplot as plt
import numpy as np

#%%
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
    result_dir = 'test_search/est_results/'
    data_source = 'PRIDE'

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
    result_dir = 'test_search/est_results/'
    data_source = ''
    psm_dir = 'test_search/pride/'

elif data_source == 'NIST':
    species_list = [
            'c_elegans',
            'drosophila',
            'e_coli',
            'human',
            'mouse'
        ]
    result_dir = 'test_search/est_results_nist/'
    data_source = 'NIST'


#%%
legends = []
def plot_repmis(species):
    species_dir = result_dir + species + '/'
    
    repmis_file = species_dir + 'repmis.csv'
    repmis_file = open(repmis_file)
    repmis_csv = csv.reader(repmis_file)
    
    pl = []
    c = 0
    n = 0
    for m, s in repmis_csv:
        m = int(m)
        s = float(s)
        if m == 0:
            c += 1
#        if m != -1:
        if m == 0 or m == 1:
            n += 1
        p = (n - c) / n
#        print(p)
        pl.append((s, p))
    
    pl = np.array(pl)

    plt.plot(pl[:,0], pl[:,1])
    legends.append(species.replace('_', '.').capitalize())

#%%
fig = plt.figure(figsize=[7,5], dpi=200)
for species in species_list:
    plot_repmis(species)
#    break
plt.xlabel('-log(EValue)')
plt.ylabel('Error rate')
plt.legend(legends)










