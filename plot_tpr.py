#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 00:19:17 2020

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
        'drosophila',
        'e_coli',
        'human',
        'mouse',
    ]

methods_list = [
        '_1s2c',
        '_1s2ca',
        '_2s3c',
        'tda',
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
        
        if row['Protein'][0:3] == 'REV':
            continue
        
        ind = int(row['SpecID'].split('=')[1])
        if prevspec == ind:
            continue
        prevspec = ind
        
#        qval = float(row['QValue'])
#        if qval == 1:
#            continue
        
        s = -np.log(float(row['SpecEValue']))
        yield (ind, row['Peptide'], s)
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
        if j >= len(p_msgf):
            return False
        if p_msgf[j] in '+.0123456789':
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
        'TDA': 'tda',
    }
algo_map = {v:k for k,v in method_map.items()}

xgrid = np.arange(0, 0.25, 0.001)

def get_tprcurvs(species):
    curvs = {}
    
    numspec = open('nist/'+species+'_cnt.txt')
    numspec = int(numspec.readline())
    print(numspec)
    
    species_dir = res_dir + species + '/'
    
    def get_peps(species):
        #mspfile = open('nist/human_consensus_final_true_lib.msp')
        mspfile = tarfile.open('nist/'+species+'_consensus_final_true_lib.tar.gz', 'r|gz')
        tarinfo = mspfile.next()
        mspfile = mspfile.extractfile(tarinfo)
        
        msplist = (l.decode("utf-8").rstrip() for l in mspfile)
        
        peps = []
        for item in parse_msp(msplist):
            if not item:
                continue
            peps.append(item['Name'])
        return peps
    
    peps = get_peps(species)
    
    def get_curv(species, method):
        if method == 'tda':
            f = open('nist/'+species+'_d.tsv')
        else:
            f = open('nist/'+species+'_nod.tsv')
        psmcsv = csv.DictReader(f, delimiter='\t')
        
        psms = get_first_psms(psmcsv)
        matches = []
        for spec, pep, score in psms:
        #    print(peps[spec], pep)
            matches.append((match_peptide(peps[spec], pep), score))
        #    print(match_peptide(peps[spec], pep))
            
        n = 0
        fc = 0
        curv = []
        ncorr = 0
        nwrong = 0
        for m, score in matches:
        #    print(row)
            n += 1
            if not m:
                fc += 1
                nwrong += 1
            else:
                ncorr += 1
            tpr = ncorr / numspec
            curv.append((score, tpr))
        
        print('correct vs wrong', ncorr, nwrong)
        
        return np.array(curv)
    
    for method in methods_list:
        curvs[algo_map[method]] = get_curv(species, method)


#    fdr_dir = species_dir + 'fdr/'
#    truefdr_file = fdr_dir + 'true.csv'
#    
#    np.savetxt(truefdr_file, true_fdr)
    
    return curvs
        
#%%
#def plot_tda_tpr(species):
#    tdares = tdajson[species]
#    
#    
#    s1 = tdares['s1']
#    n = len(s1)
#    tpr = []
#    for i in range(0,n):
#        tpr.append(i / numspec)
#    plt.plot(s1, tpr)

#%%
legends = []
for species in species_list:
    curvs = get_tprcurvs(species)

    fig = plt.figure(figsize=[7,5], dpi=200)

    for method, curv in curvs.items():
        plt.plot(curv[:,0], curv[:,1], linewidth=.5)
        legends.append(method)
#        plt.show()
    plt.legend(legends)
    
    species_dir = res_dir + species + '/'
    
    plt.savefig(species_dir+'tprcurv.png')

#plt.ylim(-0.05,1.05)