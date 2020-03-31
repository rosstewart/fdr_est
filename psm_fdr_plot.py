#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:24:01 2019

@author: yisupeng
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sio

#method = '_2mix'
#method = '_3mix'
method = '_1s2c'
#method = ''

method_list = [
        '_1s2c',
        '_1s2ca',
#        '_2s2ci',
        '_2s3ci',
#        '_2s3ct',
#        '_3s4ci'
    ]

species = 'A.thaliana'
#species = 'C.elegans'
#species = 'D.melanogaster'
#species = 'E.coli'
#species = 'H.sapiens2'
#species = 'H.sapiens3'
#species = 'M.musculus'
#species = 'M.musculus2'
#species = 'M.musculus3'
#species = 'S.cerevisiae'
#species = 'S.cerevisiae2'
#species = 'S.cerevisiae3'

#species_list = [
#        'A.thaliana',
##        'C.elegans',
#        'D.melanogaster',
#        'E.coli',
#        'H.sapiens2',
#        'H.sapiens3',
#        'M.musculus',
#        'M.musculus2',
#        'M.musculus3',
##        'S.cerevisiae',
#        'S.cerevisiae2',
#        'S.cerevisiae3',
#    ]

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
psm_dir = 'test_search/pride/'
#psm_dir = 'test_search/pride/fragger/'
#psm_dir = 'nist/'

#result_dir = 'test_search/est_results_nist/'
result_dir = 'test_search/est_results/'

thres_output = open(result_dir + 'thresholds.txt', 'w')
thres_output.write('species\tmethod\tthres\tthres_cor\tthres2\n')

def plot_fdrcurv(species, method):

    #matches = open('test_search/Adult_Adrenalgland_Gel_Elite_49_f01_d.tsv')
    #matches = open('test_search/isb02_t1.tsv')
    #matches = open('test_search/pride/H.sapiens_d.tsv')
    matches = open(psm_dir+'%s_d.tsv'%species)
    matches_csv = csv.DictReader(matches, delimiter='\t')
    
    #evalue_field = 'EValue'
    evalue_field = 'SpecEValue'
    
    n = 0
    
    curv = []
    
    spec_pep_map = {}
    
    for row in matches_csv:
    #    print(row)
        if row['Protein'][0:3] == 'REV':
            continue
        qval = float(row['QValue'])
        if qval <= 0.1:
        #    s = float(row['MSGFScore'])
            s = -np.log(float(row[evalue_field]))
            curv.append((qval, n, s))
        if qval <= 0.01:
    #        spec = row['ScanNum']
            spec = row['SpecID']
            pep = row['Peptide']
            spec_pep_map[spec] = pep
        n += 1
    curv = np.array(curv)
    
    #%% pep curve
    #matches = open('test_search/isb02_t1.tsv')
    #matches = open('test_search/pride/H.sapiens_d.tsv')
    matches = open(psm_dir+'%s_d.tsv'%species)
    matches_csv = csv.DictReader(matches, delimiter='\t')
    
    pepcurv = []
    
    peps = set()
    
    peps_above = set()
    
    npep = 0
    nrep = 0
    ndecoy = 0
    qval = 0 
    for row in matches_csv:
    #    spec = row['ScanNum']
        spec = row['SpecID']
        pep = row['Peptide']
        if row['Protein'][0:3] == 'REV':
            qval = ndecoy / npep
            if pep not in peps:
                ndecoy += 1
                peps.add(pep)
                npep += 1
            continue
        if pep not in peps:
            npep += 1
            peps.add(pep)
        else:
            nrep += 1
            continue
        
        if qval <= 0.1:
            s = -np.log(float(row[evalue_field]))
            pepcurv.append((qval, npep, s))
        if qval <= 0.01:
    #        print(pep, qval)
            peps_above.add(pep)
        n += 1
    
    
    pepcurv = np.array(pepcurv)
    
    #%%
    species_dir = result_dir+species#+'/'
    twomix = open(species_dir+'/fdr/'+method+'.csv')
    twomix = np.genfromtxt(species_dir+'/fdr/'+method+'.csv', delimiter=',')
    #%%
    fig = plt.figure()
    plt.plot(curv[:,0], curv[:,1])
    plt.plot(twomix[:,0], twomix[:,1])
    
    plt.legend(['TDA approach', 'two mixture approach'])
    plt.show()
    fig.savefig(species_dir+'/fdrcurv/%s.png'%(method))
    
    #%%
    paramfile = result_dir+species+'/params/'+method+'.mat'
    theta = sio.loadmat(paramfile)['theta']
    alpha = theta['alpha'][0,0][0,0]
    print(alpha)
    
    #%% 1% fdr
    nmatches = 0
    for fdr,n,s in curv:
        if fdr > 0.01:
            break
        nmatches = n
        thres = s
    
    for fdr,n,s in curv:
        if (1-alpha)*fdr > 0.01:
            break
        nmatches = n
        thres_cor = s
    
    nmatches2 = 0
    eval_fdr_map = {}
    
    thres2 = 60
    for fdr,n,s in twomix:
        eval_fdr_map[s] = fdr
    for fdr,n,s in twomix:
        if fdr > 0.01:
            break
        nmatches2 = n
        thres2 = s
    
    
    thres_output.write(species+'\t'+method+'\t')
    thres_output.write(str(thres)+'\t')
    thres_output.write(str(thres_cor)+'\t')
    thres_output.write(str(thres2)+'\n')
    
    #%% Decoy free matches
    
    #matches = open('test_search/Adult_Adrenalgland_Gel_Elite_49_f01_d.tsv')
    #matches = open('test_search/isb02_t3_nodecoy.tsv')
    #matches = open('test_search/pride/H.sapiens_nod.tsv')
    matches = open(psm_dir+'%s_nod.tsv'%species)
    matches_csv = csv.DictReader(matches, delimiter='\t')
    
    n2 = 0
    
    npep2 = 0
    peps2 = set()
    specs2 = set()
    peps_above2 = set()
    
    spec_pep_map2 = {}
    pepcurv2 = []
    
    for row in matches_csv:
    #    print(row)
        if row['Protein'][0:3] == 'REV':
            print('unexpected decoy')
            continue
    #    spec = row['ScanNum']
        spec = row['SpecID']
        pep = row['Peptide']
        
        s = -np.log(float(row[evalue_field]))
        
        if s not in eval_fdr_map:
            continue
        
        if spec in specs2:
            continue
        
        specs2.add(spec)
        
        fdr = eval_fdr_map[s]
    #    print(spec, s, fdr)
        if fdr <= 0.1:
            pepcurv2.append((fdr, npep2, s))
        if s > thres2:
            peps_above2.add(pep)
            if (spec not in spec_pep_map2):
                spec_pep_map2[spec] = pep        
        n2 += 1
        if pep not in peps2:
            peps2.add(pep)
            npep2 += 1
    
    pepcurv2 = np.array(pepcurv2)
    #%% pep Venn diagram
    fig = plt.figure()
    plt.plot(pepcurv[:,0], pepcurv[:,1])
    plt.plot(pepcurv2[:,0], pepcurv2[:,1])
    
    plt.legend(['TDA approach', 'two mixture approach'])
    plt.show()
    fig.savefig(species_dir+'/fdrcurv/%s_pep.png'%(method))
    #%% Comp matches
    missed_matches2 = []
    same_matches = []
    diff_matches = []
    for spec,pep2 in spec_pep_map2.items():
        if spec not in spec_pep_map:
            missed_matches2.append((spec,pep2))
            continue
        pep = spec_pep_map[spec]
        if pep == pep2:
            same_matches.append((spec, pep))
        else:
            diff_matches.append((spec,(pep,pep2)))
    
    missed_matches = []
    for spec,pep in spec_pep_map.items():
        if spec not in spec_pep_map2:
            missed_matches.append((spec,pep))
            continue
    
    n_pepmatches = 0
    for pep in peps_above:
        if pep in peps_above2:
            n_pepmatches += 1
    
    #%%
    print('matches in nodecoy found in decoy %d' % (len(same_matches)+len(diff_matches)))
    print('same matches %d' % (len(same_matches)))
    print('diff matches %d' % (len(diff_matches)))
    print('matches in nodecoy not in decoy %d' % len(missed_matches2))
    print('matches in decoy not in nodecoy %d' % len(missed_matches))


#%%
for species in species_list:
    for method in method_list:
        plot_fdrcurv(species, method)
        

#%%
thres_output.close()