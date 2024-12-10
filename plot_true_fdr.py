#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 22:10:26 2020

@author: yisupeng
"""

import os
import tarfile

import csv
import numpy as np
import matplotlib.pyplot as plt

import json


#%%

#species = 'c_elegans'
#species = 'drosophila'
#species = 'e_coli'
#species = 'human'
#species = 'mouse'
#species = 'yeast'


species_list = [
        'c_elegans',
        'h_sapiens',
        'm_musculus',
        's_cerevisiae',
#        'drosophila',
#        'e_coli',
#        'human',
#        'mouse',
#        'yeast',
    ]

#%%
    
#res_dir = 'test_search/est_results_nist/'
res_dir = 'data/nist/true_fdr_25ppm/'

#%%

#tdajson = json.load(open(res_dir+'json/tda_fdr.json'))
#tdajson = {res['ds']:res for res in tdajson}

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

log_scale = False
#log_scale = True

#skip_nomatch = True
skip_nomatch = False

def plot_fdrcurv(species):
    species_dir = res_dir + species + '/'
    
    fdr_dir = species_dir + 'fdr/'
    truefdr_file = fdr_dir + 'true.csv'
    
    true_fdr = np.genfromtxt(truefdr_file)
    

    fdr_dir = species_dir + 'fdr/'
    #truefdr_tda_file = fdr_dir + 'true_tda.csv'
    #true_fdr_tda = np.genfromtxt(truefdr_tda_file)
    #print(true_fdr_tda.shape)
    
    fdr_correction_file = species_dir + 'fdr_correction.csv'
    fdr_correction = np.genfromtxt(fdr_correction_file)
    
    if skip_nomatch:
        nomatch = np.genfromtxt(species_dir + 'nomatch.csv', dtype='int32')
    
    def draw_fdrcurv():
        estres = json.load(open(res_dir+'json/'+species+'.json'))
        #estres = json.load(open('data/nist/matdata_25ppm/'+species+'_rep.json'))
            
        max_fdr = max(true_fdr)
        
        fig = plt.figure(figsize=[7,5], dpi=200)
#        fig = plt.figure(figsize=[7,5])
        ax = plt.axes()
#        ax = fig.add_axes()
        
        legends = []
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
            ax.plot(est_fdr, true_fdr, linewidth=1)
            legends.append(res['algo'])
        
#        shares = json.load(open('test_search/shantanu/'+'json/'+species+'.json'))
#        
#        legends = []
#        for res in shares:
#            method = method_map[res['algo']]
#            est_fdr = list(reversed(res['fdr']))
#            max_fdr = np.maximum(max_fdr, max(est_fdr))
#            ax.plot(est_fdr, true_fdr, linewidth=1)
#            legends.append(res['algo'])
        
        #tdares = tdajson[species]
        #est_fdr = tdares['fdr']
        #print(len(est_fdr))
        
        #ax.plot(est_fdr, true_fdr_tda, linewidth=1)
        legends.append('tda')
        
        ax.plot([0,max_fdr], [0,max_fdr], linewidth=1, color='darkred',linestyle='dashed')

        ax.set_aspect('equal')
        ax.set_xlabel('Estimated FDR')
        ax.set_ylabel('Percent of mismatches vs. NIST')
        
        ax.set_xlim(1e-3, 0.1)
        ax.set_ylim(1e-3, 0.1)
        
        if log_scale:
            plt.xscale('log')
            plt.yscale('log')
        
#        legends.append('ground truth')
        ax.legend(legends)
        
        #curv_dir = species_dir + 'fdrcmp/'
        #
        #if not os.path.exists(curv_dir):
        #    os.mkdir(curv_dir)
        
        fdrcurv_dir = 'data/nist/fdr_result_25ppm/' + species + '/fdrcurv/'
        
        if log_scale:
            plt.savefig(fdrcurv_dir + 'fdrcmp_log.png')
        else:
            plt.savefig(fdrcurv_dir + 'fdrcmp.png')
        
        plt.show()
    
    draw_fdrcurv()
        
#%%

log_scale = False
for species in species_list:
    plot_fdrcurv(species)
log_scale = True
for species in species_list:
    plot_fdrcurv(species)

