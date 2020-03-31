#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:48:57 2020

@author: yisupeng
"""


import json
import csv
import numpy as np
import matplotlib.pyplot as plt

#%%

species_list = [
        'A.thaliana',
#        'C.elegans',
        'D.melanogaster',
        'E.coli',
        'H.sapiens2',
        'H.sapiens3',
        'M.musculus',
        'M.musculus2',
        'M.musculus3',
#        'S.cerevisiae',
        'S.cerevisiae2',
        'S.cerevisiae3',
        
        'HeLa01ng_2',
        'HeLa1ng',
        'HeLa10ng',
        'HeLa50ng',
        'HeLa100ng'
    ]

#species_list = [
#        'E.coli',
#    ]

psm_dir = 'test_search/pride/'

res_dir = 'test_search/est_results_full/'
#res_dir = 'test_search/est_results_8inits/'
#res_dir = 'test_search/est_results_part/'

#%%
def load_data(species):
#    species_dir = 'test_search/est_results/'+species#+'/'
    
    matches = open(psm_dir+'%s_d.tsv'%species)
    matches_csv = csv.DictReader(matches, delimiter='\t')
    
    evalue_field = 'SpecEValue'

    decoy_flags = []
    
    s1_tda = []
    fdr_tda = []
    flag = False
    thres = None
    
    data_mat = []
    for row in matches_csv:
        if row['Protein'][0:3] == 'REV':
            decoy_flag = 1
        else:
            decoy_flag = 0
        
        s = -np.log(float(row[evalue_field]))
        
        qval = float(row['QValue'])

        data_mat.append([s, decoy_flag, qval])
    
    return np.array(data_mat)
#%%
ind_s = 0
ind_rev = 1
ind_qval = 2
def calc_qvalue(data):
    n = len(data)
    decoyn = 0
    qstarti = 0
    qval = 0
#    data = np.hstack([data, np.zeros([n,1])])
    for i in range(n):
        row = data[i]
        if row[ind_rev]:
            qval = max(decoyn / (i+1 - decoyn), qval)
#            qval = decoyn / (i+1 - decoyn)
            decoyn += 1
            data[qstarti:i, ind_qval] = qval
            qstarti = i
    data[qstarti:n, ind_qval] = qval
    return data

def get_thres_fdr(data, fdr_thres):
    return 1
#%%
method = '_tda'
def bootstrapping(species, nboots, boot_ratio):
    data = load_data(species)
    species_dir = res_dir + species + '/'
    fdr_dir = species_dir + 'fdr/' + method + '/bootstrap/'
    n = len(data)
    print(n)
    boot_size = int(boot_ratio * n)
    thres_list = []
    for boot_num in range(nboots):
        print(species, boot_num)
        ind = sorted(np.random.choice(range(n), boot_size))
#        print(ind)
        boot_data = data[ind, :]
#        plt.plot(boot_data[:,0], boot_data[:,2])
        
        boot_data = calc_qvalue(boot_data)
        
        for i in range(boot_size):
            fdr = boot_data[i, ind_qval]
            thres = boot_data[i, ind_s]
            if fdr >= 0.01:
                print(fdr)
                break
        
        print(thres)
        thres_list.append(thres)
        
#        plt.plot(boot_data[:,0], boot_data[:,2])
                
#        fdr_file = fdr_dir + str(boot_num) + '.csv'
#        savetxt(fdr_file, data, delimiter=',')
    return thres_list

#%%
res = []
for species in species_list:
    thres_list = bootstrapping(species, 100, 0.1)
    obj = {'ds': species, 'algo': 'TDA', 't1p': thres_list}
#    thres_mat.append(thres_list)
#    plt.legend(species_list)
    res.append(obj)
json.dump(res, open(res_dir+'json/bootstrap/tda.json', 'w'))

#%%
thres_mat = [obj['t1p'] for obj in res]
plt.boxplot(thres_mat)
#%%
#bootstrapping('H.sapiens3', 100, 0.1)

#%%
#species = 'H.sapiens3'
#data = load_data(species)
##plt.plot(data[:,2])
#data2 = calc_qvalue(data)
#
#qdiff = data2[:,2] - data[:,2]
#cmp_mat = np.hstack([qdiff.reshape([-1,1]), data, data2])
