#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:55:20 2020

@author: yisupeng
"""



import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from matplotlib.patches import Polygon
from statistics import mean, stdev
import json
import os

#%%
method_list = [
        '_1s2ca',
        '_1s2c',
#        '_2s2ci',
        '_2s3ci',
#        '_2s3ct',
#        '_3s4ci'
    ]

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
        
#        'HeLa01ng',
#        'HeLa1ng',
#        'HeLa10ng',
#        'HeLa50ng',
#        'HeLa100ng',
    ]

#species_list = [
#        'c_elegans',
#        'drosophila',
#        'e_coli',
#        'human',
#        'mouse',
#    ]
#
#species_list = [
#        'E.coli',
#    ]

results_dir = 'test_search/est_results_full/'
#results_dir = 'test_search/est_results_8inits/'
#results_dir = 'test_search/est_results_part/'
#results_dir = 'test_search/est_results_nist/'


#%%

nboots = 200
def boxplot_thres(species, method):
    thres_list = []
    for bootstrap_num in range(nboots):
        species_dir = results_dir + '/'+species+'/'
        fdr_dir = species_dir + 'fdr/'
        bootstrap_dir = fdr_dir + method+'/bootstrap/'
        fdr_file = bootstrap_dir + str(bootstrap_num+1)+'.csv'
#        print(fdr_file)
        
        fdr_est = np.genfromtxt(fdr_file, delimiter=',')
        if len(fdr_est.shape) == 1:
            thres_list.append(max(thres_list))
            continue
        
#        plt.plot(fdr_est[:,0])
#        print(fdr_est.shape)
        
        eval_fdr_map = {}
        
#        print(bootstrap_num, fdr_file, fdr_est.shape)
        thres2 = 0
        for fdr,n,s in fdr_est:
            eval_fdr_map[s] = fdr
        for fdr,n,s in fdr_est:
            nmatches2 = n
            thres2 = s
#            print(s)
            if fdr > 0.01:
                break
        thres_list.append(thres2)
    plt.plot(thres_list)
    return thres_list


def boxplot_dcdf(species, method):
    thres_list = []
    for bootstrap_num in range(nboots):
        species_dir = results_dir + '/'+species+'/'
        sdcdf_dir = species_dir + 'sdcdf/'
        bootstrap_dir = sdcdf_dir + method+'/bootstrap/'
        sdcdf_file = bootstrap_dir + str(bootstrap_num+1)+'.mat'
        
        sdcdf_est = sio.loadmat(sdcdf_file)['sdcdf'][0,0]
        
        thres_list.append(sdcdf_est)
    return thres_list

#thres_list = boxplot_thres('S.cerevisiae3', '_1s2c')

#%%

ticks = []
pos = 1
thres_mat = []
dcdf_mat = []
for species in species_list:
    for method in method_list:
        if method == '_1s2c':
            ticks.append(species+'_1S2D skew normal')
        elif method == '_1s2ca':
            ticks.append(species+'_1S2D gamma & gaussion')
        elif method == '_2s3ci':
            ticks.append(species+'_2S3D skew normal')
        thres_list = boxplot_thres(species, method)
        dcdf_list = boxplot_dcdf(species, method)
#        ax.boxplot(thres_list, positions=[pos])
        thres_mat.append(thres_list)
        dcdf_mat.append(dcdf_list)
        print(species, method)
        pos += 1

#%%
i = 0
for method in method_list:
    if method == '_1s2c':
        algo = ('1S2D skew normal')
    elif method == '_1s2ca':
        algo = ('1S2D gamma & gaussion')
    elif method == '_2s3ci':
        algo = ('2S3D skew normal')
    print(',\t'+algo, end='')
print('')

algo_map = {'_1s2c': '1S2D', '_1s2ca': '1S2D GG', '_2s3ci': '2S3D'}
bootstrap_dir = results_dir + '/json/bootstrap/'
if not os.path.isdir(bootstrap_dir):
    os.mkdir(bootstrap_dir)
for species in species_list:
    print(species, end=',\t')
    dataset_arr = []
    for method in method_list:
        obj = {'ds': species, 'algo': algo_map[method], 't1p': thres_mat[i]}
        dataset_arr.append(obj)
        mu = mean(thres_mat[i])
        sigma = stdev(thres_mat[i])
        print("%.2f"%mu, "±", "%.2f"%sigma, end=',\t')
        
        i += 1
    print('')
    json.dump(dataset_arr, open(bootstrap_dir+species+'.json', 'w'))
        

