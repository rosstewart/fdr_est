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

shantanu_m_l = ['SNMax1',]
#
#method_list += shantanu_m_l

#data_source = 'HeLa'
data_source = 'NIST'

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
    result_dir = 'test_search/est_results/'
    data_source = 'PRIDE'

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
    result_dir = 'test_search/est_results/'
    data_source = ''
    psm_dir = 'test_search/pride/'
    data_dir = 'test_search/matdata_hela/'

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
    psm_dir = 'test_search/nist/'
    data_dir = 'test_search/matdata_nist/'

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
evalue_field = 'SpecEValue'

def get_psms(species):
    psm_file = psm_dir + species + '_nod.tsv'
    ss = csv.DictReader(open(psm_file), delimiter='\t')
    
    psmlist = []
    previd = None
    for row in ss:
        specid = row['SpecID']
        if specid != previd:
            if psmlist:
                yield psmlist
            previd = specid
            psmlist = []
        psm = PSM(row)
        psmlist.append(psm)
    if psmlist:
        yield psmlist

def get_psm_lists(species):
    psm_file = psm_dir + species + '_nod.tsv'
    ss = csv.DictReader(open(psm_file), delimiter='\t')
    
    prevspec = None
    psmlist = []
    for row in ss:
    #    print(row)
        ind = int(row['SpecID'].split('=')[1])
        if prevspec != ind:
            if psmlist:
                yield (prevspec, psmlist)
                psmlist = []
        prevspec = ind
        
        psmlist.append((row['Peptide'], float(row['SpecEValue'])))
    if psmlist:
        yield (prevspec, psmlist)

#%%
def extract_mat(species):
    psms = get_psms(species)
        
    i = 0
#    mat = np.zeros([len(psms), 10])
#    smat = np.zeros([len(psms), 10])
    slen = []
    mat = []
    smat = []
        
    for psmlist in psms:
    #    print(scannum, len(psmlist))
        
        slist = [0]*10
        sslist = [0]*10
        
        slen.append(len(psmlist))
        j = 0
#        prevs = 0
        for psm in psmlist:
    #    for j, psm in enumerate(psmlist):
            if j >= 10:
                break
            if j > 0 and slist[j-1] == -np.log(float(psm[evalue_field])):
                continue
            slist[j] = -np.log(float(psm[evalue_field]))
            sslist[j] = float(psm[evalue_field])
            j += 1
        
        mat.append(slist)
        smat.append(sslist)
        i += 1

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

    
#    mat = np.array(mat)
    
#    np.savetxt(data_dir+species+'_rep.csv', mat)
    json.dump(matches, open(data_dir+species+'_rep.json', 'w'))



#%%
for species in species_list:
    print(species)
    extract_mat(species)
    extract_rep_mat(species)

