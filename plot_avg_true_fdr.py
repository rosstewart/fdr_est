#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 22:49:38 2020

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

#species = 'c_elegans'
#species = 'drosophila'
#species = 'e_coli'
#species = 'human'
#species = 'mouse'
#species = 'yeast'


species_list = [
        'c_elegans',
#        'drosophila',
#        'e_coli',
        'human',
        'mouse',
        'yeast',
    ]

#%%
    
res_dir = 'test_search/est_results_nist/'

#%%

tdajson = json.load(open(res_dir+'json/tda_fdr.json'))
tdajson = {res['ds']:res for res in tdajson}

#%%
def split_comments(comment):
    ret = ''
    quote = False
    for c in comment:
        if c == ' ' and not quote:
            yield ret
            ret = ''
        elif c == '"':
            quote = not quote
        else:
            ret += c
            

#%%

def parse_msp(msplist):
    peaks = []
    item = {}
    for l in msplist:
        if not l:
            item['Peaks'] = peaks
            
            props = {k:v for k,v in map(lambda x: x.split('=', 1), split_comments(item['Comment']))}
            item['Props'] = props
#            print(len(peaks))
            yield item
            item = {}
            peaks = []
            continue
            
        if ': ' not in l:
            peak = l.split('\t')
    #        print(peak)
            peaks.append(peak)
            continue
    #        peaks.append()
        
        [k, v] = l.split(': ', 1)
    #    print(k,v)
        if not item:
            assert(k == 'Name')
        item[k] = v
    yield item
    

#%%

#%%


def get_first_psms(ss):
    prevspec = None
        
    for row in ss:
    #    print(row)
        ind = int(row['SpecID'].split('=')[1])
        if prevspec == ind:
            continue
        prevspec = ind
        yield (ind, row['Peptide'], row['SpecEValue'])
#    break

def match_peptide(p_nist, p_msgf):
    i=0
    j=0
#    print(p_nist, p_msgf)
    if p_msgf[1] == '.':
        while p_msgf[j] != '.':
            j += 1
        j += 1
    while True:
#        print(p_nist[i])
        if p_nist[i] == '/':
            return True
        if p_msgf[j] in '-+.0123456789':
#            print(p_msgf[j], j, len(p_msgf))
            j += 1
            continue
        
        if p_nist[i] != p_msgf[j]:
            print(False, p_nist, p_msgf)
            return False
        
        i += 1
        j += 1
    
    return True



#%%

method_map = {
        '1S2D gamma & gaussian': '_1s2ca',
        '1S2D skew normal': '_1s2c',
        '2S3D skew normal': '_2s3c',
        'SNMax1': 'SNMax1',
        'SNMax2': 'SNMax2',
    }

#log_scale = False
log_scale = True

skip_nomatch = True
#skip_nomatch = False

if log_scale:
    xgrid = np.power(10, np.arange(-5, 0, 0.01))
else:
    xgrid = np.arange(0, 0.25, 0.001)

def get_fdrcurvs(species):
    curvs = {}
    
    species_dir = res_dir + species + '/'
    
    fdr_dir = species_dir + 'fdr/'
    truefdr_file = fdr_dir + 'true.csv'
    
    true_fdr = np.genfromtxt(truefdr_file)
    

    fdr_dir = species_dir + 'fdr/'
    truefdr_tda_file = fdr_dir + 'true_tda.csv'
    true_fdr_tda = np.genfromtxt(truefdr_tda_file)
    print(true_fdr_tda.shape)
    
    fdr_correction_file = species_dir + 'fdr_correction.csv'
    fdr_correction = np.genfromtxt(fdr_correction_file)
    
    if skip_nomatch:
        nomatch = np.genfromtxt(species_dir + 'nomatch.csv', dtype='int32')
    
    def get_fdrcurv():
        estres = json.load(open(res_dir+'json/'+species+'.json'))
            
        max_fdr = max(true_fdr)
        
#        fig = plt.figure(figsize=[7,5], dpi=200)
#        fig = plt.figure(figsize=[7,5])
#        ax = plt.axes()
#        ax = fig.add_axes()
        
        for res in estres:
            method = method_map[res['algo']]
            est_fdr = list(reversed(res['fdr']))
            est_fdr = np.array(est_fdr)
            if skip_nomatch:
                flags = np.ones(est_fdr.shape[0],dtype=bool)
                flags[nomatch] = 0
                est_fdr = est_fdr[flags]
                
            est_fdr += fdr_correction
            max_fdr = np.maximum(max_fdr, max(est_fdr))
            
            curvs[res['algo']] = np.interp(xgrid, est_fdr, true_fdr)
        
#        shares = json.load(open('test_search/shantanu/'+'json/'+species+'.json'))
#        
#        legends = []
#        for res in shares:
#            method = method_map[res['algo']]
#            est_fdr = list(reversed(res['fdr']))
#            max_fdr = np.maximum(max_fdr, max(est_fdr))
#            ax.plot(est_fdr, true_fdr, linewidth=1)
#            legends.append(res['algo'])
        
        tdares = tdajson[species]
        est_fdr = tdares['fdr']
        print(len(est_fdr))
        
        curvs['tda'] = np.interp(xgrid, est_fdr, true_fdr_tda)
    
    get_fdrcurv()
    
    return curvs
        
#%%
    
def plot_avg_true_fdr():
    curvs = defaultdict(list)
    for species in species_list:
        sc = get_fdrcurvs(species)
        for method, curv in sc.items():
            curvs[method].append(curv)
    
    #%%
    fig = plt.figure(figsize=[7,5], dpi=200)
    ax = plt.axes()
    legends = []
    lines = []
    for method, mat in curvs.items():
        print(method)
        mat = np.array(mat)
        print(mat.shape)
        
        mcurv = mat.mean(0)
        curvstd = mat.std(0)
        
        ucurv = mcurv + curvstd
        dcurv = mcurv - curvstd
    
        mc = plt.plot(xgrid, mcurv)
        lines.append(mc[0])
        plt.plot(xgrid, ucurv, color=mc[0].get_c(), linewidth=1, alpha=.4)
        plt.plot(xgrid, dcurv, color=mc[0].get_c(), linewidth=1, alpha=.4)
        plt.fill_between(xgrid, dcurv, ucurv, alpha=.2)
        
        legends.append(method)
    
    mc = plt.plot([0,0.2], [0,0.2], linewidth=1, color='darkred')
    lines.append(mc[0])
    legends.append('ground truth')
    
    ax.set_xlim(1e-3, 0.2)
    ax.set_ylim(1e-3, 0.2)
    
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')
    
    ax.set_aspect('equal')
    ax.set_xlabel('Estimated FDR')
    ax.set_ylabel('Percent of mismatches v.s. NIST')
    plt.legend(lines, legends)
    
    if log_scale:
        plt.savefig('test_search/est_results_nist/fdrcmp_log.png')
    else:
        plt.savefig('test_search/est_results_nist/fdrcmp.png')


#%%
log_scale = False
plot_avg_true_fdr()

log_scale = True
plot_avg_true_fdr()

