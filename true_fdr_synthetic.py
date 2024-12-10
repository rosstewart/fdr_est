#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 22:10:26 2020

@author: yisupeng-->rossstewart
"""

import os
import tarfile

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import json


#%%

#species = 'c_elegans'
#species = 'drosophila'
#species = 'e_coli'
#species = 'human'
#species = 'mouse'
#species = 'yeast'


species_list = [f'synthetic_{i}' for i in range(1,11)]
#         'c_elegans',
# #        'drosophila',
# #        'e_coli',
# #        'human',
# #        'mouse',
#         's_cerevisiae',#'yeast',
#         'h_sapiens',#'human_hcd',
#         'm_musculus',#'mouse_hcd',
    # ]

# libmap = {
#         # 'c_elegans': 'data/nist/true_fdr_25ppm/2011_05_24_c_elegans_consensus_final_true_lib.tar.gz',#c_elegans_consensus_final_true_lib.tar.gz',
#         # # 'drosophila': 'drosophila_consensus_final_true_lib.tar.gz',
#         # # 'e_coli': 'e_coli_consensus_final_true_lib.tar.gz',
#         # # 'human': 'human_consensus_final_true_lib.tar.gz',
#         # # 'mouse': 'mouse_consensus_final_true_lib.tar.gz',
#         # 's_cerevisiae': 'data/nist/true_fdr_25ppm/2012_04_06_yeast_consensus_final_true_lib.tar.gz',#yeast_consensus_final_true_lib.tar.gz',
#         # 'h_sapiens': 'data/nist/true_fdr_25ppm/2014_05_29_human_consensus_final_true_lib.tar.gz',#human_hcd_selected.msp.tar.gz',
#         # 'm_musculus': 'data/nist/true_fdr_25ppm/2013_05_20_mouse_consensus_final_true_lib.tar.gz',#cptac2_mouse_hcd_selected.msp.tar.gz',
#         #'m_musculus': 'cptac2_mouse_hcd_selected.msp.tar.gz',
#     }

#%%
    
#res_dir = 'test_search/est_results_nist/'
res_dir = 'synthetic/fdr_result/'

data_dir = 'synthetic/matdata/'

#%%
def split_comments(comment):
    ret = ''
    quote = False
    for c in comment:
        if c == ' ' and not quote:
            if '=' in ret:
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
    yielded_item = {}
    target_sequence_range_count = 0
    other_sequence_range_count = 0
    for l in msplist:
        if not l:
            item['Peaks'] = peaks
            
            props = {k:v for k,v in map(lambda x: x.split('=', 1), split_comments(item['Comment']))}
            item['Props'] = props
#            print(len(peaks))
            # only return sequence lengths of 7-16
            if len(item['Name'].split('/')[0]) >= 7 and len(item['Name'].split('/')[0]) <= 16:
                target_sequence_range_count += 1
                yielded_item[item['Name']] = True
                yield item
            else:
                other_sequence_range_count += 1
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
            # print('k:',k,'v:',v)
            assert(k == 'Name')
        item[k] = v
    # if 'Name' not in item.keys():
    #     print(item)
    # if len(item['Name'].split('/')[0]) >= 7 and len(item['Name'].split('/')[0]) <= 16 and not yielded_item[item['Name']]:
    # this item is an empty dict, will be skipped. not sure why yisu had this
    print(f'{target_sequence_range_count} with sequence length 7-16, {other_sequence_range_count} skipped (not in range)')
    yield item
    

#%%

#%%


# def get_first_psms(ss):
#     prevspec = None
        
#     for row in ss:
#     #    print(row)
#         ind = int(row['SpecID'].split('=')[1])
#         if prevspec == ind:
#             continue
#         prevspec = ind
        
#         qval = float(row['QValue'])
#         if qval == 1:
#             continue
#         yield (ind, row['Peptide'], row['SpecEValue'])
#    break

def get_first_psms(ss):
    # psm_file = psm_dir + species + '.txt'
    # with open(psm_file,'r') as file:
    prevspec = None
    i = 0
    while i < len(ss):
        line = ss[i]
        i += 1
        # if not line:
        #     break
        if line[0] == '>': # ex: >> 0 34671 YYYDHSK/2_0 (SQS 0.772)
            ind = int(line.split()[2]) # assume score files concatenated correctly
            line = ss[i] # ex: #Index	RnkScr	PnvScr	N-Gap	C-Gap	[M+H]	Charge	Sequence
            i += 1
            if line[0:6] != '#Index':
                continue
            line = ss[i] # ex: 0	6.254	74.297	0.000	0.000	975.515	2	YYYDHSK
            i += 1
            vals = line.strip().split()
            if not vals: # likely did not have 10 candidate peptides in spectrum graph
                print(line,end='')
                print(ss[i],end='')
                i += 1
                exit()
            peptide = vals[-1]
            score = vals[1]
            if float(score) < 0:
                print(float(score))
            yield (ind, peptide, max(0.1,float(score))) # shift everything to be positive, max is used for rare -99 scores (i'm lazy)

# def get_psm_lists(ss):
#     prevspec = None
#     psmlist = []
#     for row in ss:
#     #    print(row)
#         ind = int(row['SpecID'].split('=')[1])
#         if prevspec != ind:
#             if psmlist:
#                 yield (prevspec, psmlist)
#                 psmlist = []
#         prevspec = ind
        
# #        qval = float(row['QValue'])
# #        if qval == 1:
# #            continue
#         s = -np.log(float(row['SpecEValue']))
#         psmlist.append((row['Peptide'], s))
#     if psmlist:
#         yield (prevspec, psmlist)

def get_psm_lists(ss):
    # psm_file = psm_dir + species + '.txt'
    # with open(psm_file,'r') as file:
    prevspec = None
    psmlist = []
    i = 0
    while i < len(ss):
        line = ss[i]
        i += 1
        # if not line:
        #     break
        if line[0] == '>': # ex: >> 0 34671 YYYDHSK/2_0 (SQS 0.772)
            ind = int(line.split()[2]) # assume score files concatenated correctly
            line = ss[i] # ex: #Index	RnkScr	PnvScr	N-Gap	C-Gap	[M+H]	Charge	Sequence
            i += 1
            if line[0:6] != '#Index':
                continue
            for j in range(10):
                line = ss[i] # ex: 0	6.254	74.297	0.000	0.000	975.515	2	YYYDHSK
                i += 1
                vals = line.strip().split()
                if not vals: # likely did not have 10 candidate peptides in spectrum graph
                    if j == 0 or j == 1:
                        print(line,end='')
                        print(ss[i],end='')
                        i += 1
                        exit()
                    # print('only',i,'candidate peptides')
                    break
                peptide = vals[-1]
                score = vals[1]
                if float(score) < 0:
                    print(float(score))
                psmlist.append((peptide, max(0.1,float(score)))) # shift everything to be positive, max is used for rare -99 scores (i'm lazy)
            yield (ind, psmlist)
            psmlist = []
        
def match_aminoacid(a, b):
    if a == b:
        return True
    # if a == 'I' and b == 'L':
    #     return True
    # if a == 'L' and b == 'I':
    #     return True
    # if a == 'K' and b == 'Q':
    #     return True
    # if a == 'Q' and b == 'K':
    #     return True
    
    return False

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
            #print(j, p_msgf, p_nist)
        if p_msgf[j] in '-+.0123456789':
#            print(p_msgf[j], j, len(p_msgf))
            j += 1
            continue
        
#        if p_nist[i] != p_msgf[j]:
#            if p_nist[i] not in ['I', 'L'] or p_msgf[j] not in ['I', 'L']:
##                print(False, p_nist, p_msgf)
##                print(p_nist[i], p_msgf[j])
#                return False
##            print(p_nist[i], p_msgf[j])
        if not match_aminoacid(p_nist[i], p_msgf[j]):
            return False
        
        i += 1
        j += 1
    
    return True

def match_peptide_list(p_nist, psmlist):
#    print(p_nist)
#    print(psmlist)
    for i, (pep, score) in enumerate(psmlist):
        if p_nist == pep:
            return i
        # if match_peptide(p_nist, pep):
        #     return i
    return -1


#%%

method_map = {
        '1S2D gamma & gaussian': '_1s2ca',
        '1S2D skew normal': '_1s2c',
        '2S3D skew normal': '_2s3c',
    }

#skip_nomatch = True
skip_nomatch = False

def get_truefdr(species):
    species_dir = res_dir + species + '/'
    
#    repmatch = np.genfromtxt('test_search/matdata_nist/'+species+'_rep.csv')
    repmatch = json.load(open(data_dir+species+'_rep.json'))
    
#     def get_peps(species):
#         #mspfile = open('nist/human_consensus_final_true_lib.msp')
# #        mspfile = tarfile.open('test_search/nist/'+species+'_consensus_final_true_lib.tar.gz', 'r|gz')
#         mspfile = libmap[species]
#         mspfile = tarfile.open(mspfile, 'r|gz')
#         tarinfo = mspfile.next()
#         mspfile = mspfile.extractfile(tarinfo)
        
#         msplist = (l.decode("utf-8").rstrip() for l in mspfile)
        
#         peps = []
#         for item in parse_msp(msplist):
#             if not item:
#                 continue
#             peps.append(item['Name'])
#         return peps
    
#     peps = get_peps(species)
    # print('peps',species,':',peps)
    #peps = pd.read_csv(data_dir+species+'_peps.csv', delimiter='\t')
    #peps = peps['Pep'].values
    
    reps = []
    mismatches = []
    mis_rep = []
    p_mis = []
    nomatch = []
    def get_truefdr(species):
        f = open('synthetic/scores/'+species+'.txt')
        # psmcsv = csv.DictReader(f, delimiter='\t')

        # if skip_nomatch:
        #     pepfound = pd.read_csv(species_dir+'matches.csv', names=['pepfound', 'protfound', 'protid'])
        #     pepfound = pepfound['pepfound'].values
        
#        psms = get_first_psms(psmcsv)
        psmlists = get_psm_lists(f.readlines())
        f.close()
        
        matches = []
        samescores = []
        pepindb = 0
        # false_streak = 0
#        i = 0
        for i, (spec, psmlist) in enumerate(psmlists):
            print(spec,psmlist)
            pep, score = psmlist[0]
            rep = repmatch[str(spec)]
            # print(spec,len(peps))
            # print('spec:',spec,', ground truth',peps[spec])
            # print('psmlist:',psmlist)
            # m = match_peptide_list(peps[spec], psmlist)
            m = match_peptide_list('ABCDEFG', psmlist)
            
            nrep = len(rep)
            if skip_nomatch and m == -1:
                if not pepfound[spec]:
                    nomatch.append([i])
                    continue
                else:
                    pepindb += 1
            if m == 0:
                matches.append((True, score))
                # print('1',end='')
                # false_streak = 0
            else:
                matches.append((False, score))
                # print('0',end='')
                # false_streak += 1
                # if false_streak >= 1000:
                #     print(f'\nspec: {spec}, pep: {peps[spec]}, pred: {psmlist[0]}, score: {score}')
                mismatches.append([m, 'ABCDEFG',] + list(psmlist[0]))
                for j in range(1, len(psmlist)):
                    mismatches.append(['', '',] + list(psmlist[j]))
            if len(rep) > 1:
                reps.append([m, (score), 'ABCDEFG', '\n'.join([pep for pep,s in rep])])
                mis_rep.append([m, (score)])
            p_mis.append(1-1/nrep)
    
            #    print(match_peptide(peps[spec], pep))
#            i += 1
        
        print('\npep in db but not in match:', pepindb)
            
        n = 0
        fc = 0
        qstart = 0
        qvals = np.zeros(len(matches))
        curv = []
        ncorr = 0
        nwrong = 0
        # sort of descending score
        matches = np.array(matches)
        matches = matches[matches[:, 1].argsort()[::-1]]
        for m, score in matches:
        #    print(row)
            n += 1
            if not m:
                qval = fc / n
                curv.append((qval, n, score))
                # print(m, qval, n, score)
                qvals[qstart:n-1] = qval
                fc += 1
                qstart = n-1
                nwrong += 1
            else:
                # print(m,score)
                ncorr += 1
        
        print('correct vs wrong', ncorr, nwrong)
        
        qval = fc / n
        qvals[qstart:n] = qval
        curv.append((qval, n, score))
        curv = np.array(curv, dtype=float)
#        curv[:,2] = -np.log(curv[:,2])
        true_fdr = np.array(qvals)
        
#        return true_fdr
    
#    true_fdr = get_truefdr(species)
    
#    def get_fdr_correction():
        fdr_correction = []
        n = 0
        nwrong = 0
        pwrong = 0
        for p in p_mis:
            nwrong += p
            n += 1
#            if nwrong >= 5:
            # print(pwrong,nwrong,n)
            pwrong = max(pwrong, nwrong / n)
            fdr_correction.append(pwrong)
        fdr_correction = np.array(fdr_correction)
        np.savetxt(species_dir + 'fdr_correction.csv', fdr_correction)
#        return fdr_correction
#    fdr_correction = get_fdr_correction()
    
#    def get_true_thres():
        true_thres = 0
        for i, (fdr, n, s) in enumerate(curv):
            true_thres = s
            if fdr + fdr_correction[i] >= 0.01:
                break
        #    nmatches = n
        
        ttfile = open(species_dir + 'true_thres.txt', 'w')
        ttfile.write("%f\n"%true_thres)
        ttfile.close()
        
        return true_fdr
#    get_true_thres()
    true_fdr = get_truefdr(species)
    
    repcsv = csv.writer(open(species_dir + 'rep.csv', 'w'))
    repcsv.writerows(reps)
    mismatch_csv = csv.writer(open(species_dir + 'mismatch.csv', 'w'))
    mismatch_csv.writerows(mismatches)
    
    repmis_csv = csv.writer(open(species_dir + 'repmis.csv', 'w'))
    repmis_csv.writerows(mis_rep)
    
    nomatch_csv = csv.writer(open(species_dir + 'nomatch.csv', 'w'))
    nomatch_csv.writerows(nomatch)
    
    
    fdr_dir = species_dir + 'fdr/'
    truefdr_file = fdr_dir + 'true.csv'
    
    np.savetxt(truefdr_file, true_fdr)
    
    return true_fdr
        
#%%
for species in species_list:
    print(species)
    get_truefdr(species)




