#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:07:09 2020

@author: yisupeng
"""

import json

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

method_list = [
#        'SNMax1',
#        'SNMax2',
    ]

result_dir = 'test_search/est_results/'
#result_dir = 'test_search/est_results_nist/'

thres_output = open(result_dir + 'thresholds2.txt', 'w')

tdajson = json.load(open(result_dir+'json/tda_fdr.json'))
tdajson = {res['ds']:res for res in tdajson}
print(tdajson.keys())

#species = 'HeLa01ng'

#method = 'SNMax1'

estjson = None

#%%
def calc_thres_charge(species, charge, method):
    tdares = tdajson[species]
    species = species + '.c' + str(charge)
    estres = estjson[method]
    alpha = estres['pi_C']
    curv = zip(tdares['fdr'], tdares['s1'])
    print(species, method)
    #print(list(curv))
    for fdr,s in curv:
        if fdr > 0.01:
            break
        thres = s
    thres_cor = 60;
    for fdr,s in curv:
        if (1-alpha)*fdr > 0.01:
            break
        thres_cor = s
    
    thres2 = estres['t1p']
    
    thres_output.write(species+'\t'+method+'\t')
    thres_output.write(str(thres)+'\t')
    thres_output.write(str(thres_cor)+'\t')
    thres_output.write(str(thres2)+'\n')

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


#%%
for species in species_list:
    for charge in [2,3]:
        #estjson = json.load(open('test_search/shantanu/json/'+species+'.json'))
        estjson = json.load(open(result_dir+'json/'+species+'.c'+str(charge)+'.json'))
        estjson = {method_map[res['algo']]:res for res in estjson}
        #for method in method_list:
        for method, res in estjson.items():
            calc_thres_charge(species, charge, method)

#%%
thres_output.close()

    
