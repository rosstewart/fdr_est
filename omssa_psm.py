#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:10:18 2019

@author: yisupeng
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict

import scipy.io as sio

#method = '_2mix'
#method = '_3mix'
method = '_3s4c'
#method = ''

psm_dir = 'test_search/pride/omssa/'

#species = 'M.musculus'
#species = 'M.musculus2'
#species = 'M.musculus3'
species = 'H.sapiens2'
#species = 'H.sapiens3'
#species = 'C.elegans'
#species = 'D.melanogaster'
#species = 'S.cerevisiae'
#species = 'S.cerevisiae2'
#species = 'S.cerevisiae3'
#species = 'E.coli'
#species = 'A.thaliana'

matches = open(psm_dir+'%s_nod.tsv'%species)
matches_csv = csv.DictReader(matches)
print(matches_csv.fieldnames)

#%%
evalue_field = 'P-value'


psm_list = []

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
        return "%s: %s %s" % (self['Filename/id'], self['E-value'], self['Peptide'])

for row in matches_csv:
    psm = PSM(row)
    psm_list.append(psm)
#    print(row)
#curv = np.array(curv)

psm_list = sorted(psm_list, key=lambda psm: psm[evalue_field])


#%%
n = 0
ndecoy = 0
qval = 0

curv = []

spec_pep_map = {}


for psm in psm_list:
    if psm.protein[0:3] == 'REV':
        ndecoy += 1
        qval = ndecoy / n
        psm.qval = qval
        continue
#    qval = float(psm.QValue)
#    print(psm)
    psm.qval = qval
    if qval <= 0.1:
    #    s = float(psm.MSGFScore)
        s = -np.log(float(psm[evalue_field]))
        curv.append((qval, n, s))
    if qval <= 0.01:
#        spec = psm.ScanNum
        spec = psm.scannum
        if spec in spec_pep_map:
            print('spec', spec, 'already recorded, and here is a new psm', psm)
        pep = psm.peptide
        spec_pep_map[spec] = psm
    n += 1
curv = np.array(curv)
    
print(spec_pep_map)

fig = plt.figure()
plt.plot(curv[:,0], curv[:,1])

#%%
matches = open(psm_dir+'%s_nod.tsv'%species)
matches_csv = csv.DictReader(matches, delimiter=',')
print(matches_csv.fieldnames)

psm_list = []
for row in matches_csv:
    psm = PSM(row)
    psm_list.append(psm)

psms = defaultdict(list)
for psm in psm_list:
#    psms[psm.scannum].append(-np.log(float(psm[evalue_field])))
    psms[psm['Filename/id']].append(psm)
    

#%%
i = 0
mat = np.zeros([len(psms), 10])
smat = np.zeros([len(psms), 10])
slen = []
#mat2 = np.zeros(len(psms), 2)
#mat3 = np.zeros(len(psms), 3)
#omat = np.zeros(len(psms), 10)
for psmid, psmlist in psms.items():
#    print(scannum, len(psmlist))
    slen.append(len(psmlist))
    j = 0
    prevs = 0
    for psm in psmlist:
#    for j, psm in enumerate(psmlist):
        if j >= 10:
            break
        if j > 0 and mat[i, j-1] == -np.log(float(psm[evalue_field])):
            continue
        mat[i, j] = -np.log(float(psm[evalue_field]))
        smat[i, j] = psm[evalue_field]
        j += 1
    i += 1
    
slen = np.array(slen)

omat = mat
mat2 = omat[slen>=2,0:2]
mat3 = omat[slen>=3,0:3]
matobj = {'mat': mat, 'omat': omat, 'mat2': mat2, 'mat3': mat3, 'smat': smat, 'species': species}
sio.savemat('test_search/matdata/omssa/'+species+'_data.mat', matobj)

    
    
    