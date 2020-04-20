#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 21:07:09 2020

@author: yisupeng
"""

import json

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
    result_dir = 'test_search/est_results_hela/'
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
    result_dir = 'test_search/est_results_nist/'
    data_source = 'NIST'

method_list = [
        '_1s2ca',
        '_1s2c',
        '_2s3ci',
    ]

thres_output = open(result_dir + 'thresholds2.txt', 'w')
thres_output.write('species\tmethod\tthres\tthres_cor\tthres2\n')

tdajson = json.load(open(result_dir+'json/tda_fdr.json'))
tdajson = {res['ds']:res for res in tdajson}
print(tdajson.keys())

#species = 'HeLa01ng'

#method = 'SNMax1'

estjson = None

#%%
def calc_thres(species, method):
    print(species, method)
    species_dir = result_dir + species + '/'
    
    tdares = tdajson[species]
    estres = estjson[method]
    alpha = estres['pi_C']
    tdacurv = zip(tdares['fdr'], tdares['s1'])
    
    curv = np.array(list(reversed(list(zip(estres['fdr'], estres['s1'])))))
    
    if data_source == 'NIST':
        fdr_correction_file = species_dir + 'fdr_correction.csv'
        fdr_correction = np.genfromtxt(fdr_correction_file)
        print(fdr_correction.shape)
    
        curv[:,0] += fdr_correction
#    print(tdacurv)
#    plt.plot(tdacurv[:,0], tdacurv[:,1])
#    return
    for fdr,s in tdacurv:
        thres = s
        if fdr > 0.01:
            break
    
#    thres_cor = thres
    for fdr,s in tdacurv:
        thres_cor = s
        if (1-alpha)*fdr > 0.01:
            break

    for fdr,s in curv:
        thres2 = s
        if fdr > 0.01:
            break
        
    print(thres2, estres['t1p'])
#    thres2 = estres['t1p']
    
    thres_output.write(species+'\t'+method+'\t')
    thres_output.write(str(thres)+'\t')
    thres_output.write(str(thres_cor)+'\t')
    thres_output.write(str(thres2)+'\n')


method_map = {
        '1S2D gamma & gaussian': '_1s2ca',
        '1S2D skew normal': '_1s2c',
        '2S3D skew normal': '_2s3ci',
        'SNMax1': 'SNMax1',
        'SNMax2': 'SNMax2',
        'TDA': 'tda',
    }
algo_map = {v:k for k,v in method_map.items()}


for species in species_list:
    #estjson = json.load(open('test_search/shantanu/json/'+species+'.json'))
    estjson = json.load(open(result_dir+'json/'+species+'.json'))
    estjson = {method_map[res['algo']]:res for res in estjson}
    #for method in method_list:
    for method, res in estjson.items():
        calc_thres(species, method)

#%%
thres_output.close()

    
