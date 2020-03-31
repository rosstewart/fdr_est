#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:25:55 2020

@author: yisupeng
"""

import json

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
    ]

algo_list = [
        'GG',
        'SNMixOne',
        'SNMixTwo',
        ]

algo_short_map = {
        '1S2D gamma & gaussian': 'GG',
        '1S2D skew normal': 'SNMixOne',
        '2S3D skew normal': 'SNMixTwo'
        }
algo_map = {
        'GG': '1S2D gamma & gaussian',#: 'GG',
        'SNMixOne': '1S2D skew normal',# : 'SNMixOne'
        'SNMixTwo': '2S3D skew normal',# : 'SNMixTwo'
        }

#%%
terms = [
        'deltaCdf',
        't01p',
        't1p',
        't10p',
        'll',
        
        ]
#%%
#res_arrs = {}
#for term in terms:
#    res_arrs[term] = {}
#    for algo in algo_list:
#        res = []
#        res_arrs[term][algo] = res
#        for species in species_list:
#            objarr = json.load(open('test_search/est_results/json/'+species+'.json'))
#            
#            for obj in objarr:
#                algo = algo_short_map[obj['algo']]
#                if obj['algo'] == algo_map[algo]:
#                    print('%.4f' % (obj[term]))
#                    res.append(obj[term])
#%%
from collections import defaultdict
res_arrs = defaultdict(lambda : defaultdict(list))
for species in species_list:
    objarr = json.load(open('test_search/est_results/json/'+species+'.json'))
    
    for obj in objarr:
        algo = algo_short_map[obj['algo']]
        for term in terms:
            res = res_arrs[term][algo]
            print('%.4f' % (obj[term]))
            res.append('%.4f'%(obj[term]))
        
#%%
for term in terms:
    print(term)
    for algo in algo_list:
        print('&\\' + '&'.join([algo]+res_arrs[term][algo]) + '\\\\')
    print('')





