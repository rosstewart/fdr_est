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
import matplotlib.pyplot as plt

import json
from collections import defaultdict

import scipy.io as sio


#%%
method_list = [
        '_1s2c',
        '_1s2ca',
#        '_2s2ci',
        '_2s3ci',
#        '_2s3ct',
#        '_3s4ci'
    ]

species_list = ['c_elegans','h_sapiens','m_musculus','s_cerevisiae']

psm_dir = 'pepnovo/scores_trypsin_correct_pm_7-16/'
data_dir = 'pepnovo/nist/matdata/'

if not os.path.exists(psm_dir):
    os.mkdir(psm_dir)
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

#%%
    
class Entity(object):
    def __init__(self, rawdata):
        self.__dict__['_raw'] = rawdata
    def __getattr__(self, key):
        if key == '_raw':
            return
        if key in self._raw:
            return self._raw[key]
        else:
            return super().__getattribute__(key)
    def __setattr__(self, key, val):
        self._raw[key] = val
    def __getitem__(self, key):
        return self._raw[key]
    def __repr__(self):
        return str(self._raw)

class PSM(Entity):
    def __repr__(self):
        return "%s: %s %s" % (self.scannum, self.expectscore, self.peptide)


#%%
# evalue_field = 'SpecEValue'

# def get_psms(species):
#     psm_file = psm_dir + species + '_nod.tsv'
#     ss = csv.DictReader(open(psm_file), delimiter='\t')
    
#     psmlist = []
#     previd = None
#     for row in ss:
#         specid = row['SpecID']
#         if specid != previd:
#             if psmlist:
#                 yield psmlist
#             previd = specid
#             psmlist = []
#         psm = PSM(row)
#         psmlist.append(psm)
#     if psmlist:
#         yield psmlist

def get_psm_lists(species):
    psm_file = psm_dir + species + '.txt'
    with open(psm_file,'r') as file:
        prevspec = None
        psmlist = []
        spectra_seen = 0

        while True:
            line = file.readline()
            if not line:
                break
            if line[0] == '>': # ex: >> 0 34671 YYYDHSK/2_0 (SQS 0.772)
                ind = int(line.split()[2]) # assume score files concatenated correctly
                line = file.readline() # ex: #Index	RnkScr	PnvScr	N-Gap	C-Gap	[M+H]	Charge	Sequence
                if line[0:6] != '#Index':
                    continue
                for i in range(10):
                    line = file.readline() # ex: 0	6.254	74.297	0.000	0.000	975.515	2	YYYDHSK
                    vals = line.strip().split()
                    if not vals: # likely did not have 10 candidate peptides in spectrum graph
                        if i == 0 or i == 1:
                            print(line,end='')
                            print(file.readline(),end='')
                            exit()
                        # print('only',i,'candidate peptides')
                        break
                    peptide = vals[-1]
                    score = vals[1]
                    if float(score) + 20 < 0:
                        print(float(score))
                    psmlist.append((peptide, max(0.1,float(score) + 20.0))) # shift everything to be positive, max is used for rare -99 scores (i'm lazy)
                yield (ind, psmlist)
                psmlist = []

#%%
def extract_mat(species):
    psm_lists = get_psm_lists(species)
        
    i = 0
#    mat = np.zeros([len(psms), 10])
#    smat = np.zeros([len(psms), 10])
    slen = []
    mat = []
    smat = []
        
    for spec, psmlist in psm_lists:
    #    print(scannum, len(psmlist))
        
        slist = [0]*10
        sslist = [0]*10
        
        slen.append(len(psmlist))
        j = 0
#        prevs = 0
        for psm in psmlist:
            peptide, score = psm
    #    for j, psm in enumerate(psmlist):
            if j >= 10:
                break
            if j > 0 and slist[j-1] == score:
                continue
            slist[j] = score
            sslist[j] = np.e**(-score)
            j += 1
        
        mat.append(slist)
        smat.append(sslist)
        i += 1

    print(species, i)
#    return

    slen = np.array(slen)
    
    mat = np.array(mat)
    smat = np.array(smat)
    omat = mat
    print(omat)
    mat2 = omat[slen>=2,0:2]
    mat3 = omat[slen>=3,0:3]
    matobj = {'mat': mat, 'omat': omat, 'mat2': mat2, 'mat3': mat3, 'smat': smat,
              'species': species}
#        sio.savemat('test_search/matdata_hela/'+species+'_c'+c+'_data.mat', matobj)
    sio.savemat(data_dir+species+'_data.mat', matobj)


def extract_rep_mat(species):
    psm_lists = get_psm_lists(species)
        
    i = 0
#    mat = np.zeros([len(psms), 10])
#    smat = np.zeros([len(psms), 10])
    slen = []
    matches = {}
    
    for spec, psmlist in psm_lists:
    #    print(scannum, len(psmlist))
        
#        slist = [0]*10
#        sslist = [0]*10
        slist = []
        
        slen.append(len(psmlist))
        j = 0
#        prevs = 0
        prevs = 0
        for pep, s in psmlist:
            
            if prevs and prevs != s:
                break
            
            prevs = s
            slist.append((pep, s))
            
            j += 1
        
        matches[spec] = slist
        i += 1
    print(species, i)
    #return

    
#    mat = np.array(mat)
    
#    np.savetxt(data_dir+species+'_rep.csv', mat)
    json.dump(matches, open(data_dir+species+'_rep.json', 'w'))



#%%
for species in species_list:
    print(species)
    extract_mat(species)
    extract_rep_mat(species)

