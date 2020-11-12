#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 22:10:26 2020

@author: yisupeng
"""

import os
import tarfile
from glob import glob

import csv
import numpy as np
import matplotlib.pyplot as plt

import json
from collections import defaultdict

import scipy.io as sio


#%%
method_list = [
        '_1s3c_xl',
    ]

shantanu_m_l = ['SNMax1',]
#
#method_list += shantanu_m_l

data_source = 'XL'

if data_source == 'XL':
    species_list = [
            'ecoli',
        ]
    result_dir = 'test_search/est_results_xl/'
    data_source = 'XL'
    psm_dir = 'test_search/xlms/'
    data_dir = 'test_search/matdata_xl/'

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
evalue_field = 'E-value'

#def get_psms(species):
#    psm_file = psm_dir + species + '_nod.tsv'
def get_psms(psm_file):
    ss = csv.DictReader(open(psm_file), delimiter='\t')
    
    psmlist = []
    previd = None
    for row in ss:
        if row['Obs Mass'] == '0': # skip 0 mass spec
            continue
        if row['Protein #2'] == '-': # skip no 2nd pep match
            continue
        if row['Peptide #2 Score'] == '0': # skip no 2nd pep match
            continue
        if 'REV_' in row['Protein #1'] or 'REV_' in row['Protein #2']: # skip decoy match
            continue
        specid = row['Scan Number']
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
        
    i = 0
#    mat = np.zeros([len(psms), 10])
#    smat = np.zeros([len(psms), 10])
    slen = []
    mat = []
    smat = []
    
    files = glob(psm_dir + species + "/*.csv")
    for psmfile in files:
        print(psmfile)
        psms = get_psms(psmfile)
        for psmlist in psms:
        #    print(scannum, len(psmlist))
            
            slist = [0]*10
            sslist = [0]*10
            
            j = 0
    #        prevs = 0
            for psm in psmlist:
        #    for j, psm in enumerate(psmlist):
                if j >= 10:
                    break
                if j > 0 and slist[j-1] == -np.log(float(psm[evalue_field])):
                    # same score at top two psms
                    continue
                slist[j] = -np.log(float(psm[evalue_field]))
                sslist[j] = float(psm[evalue_field])
                j += 1
            
            slen.append(len(slist))
            mat.append(slist)
            smat.append(sslist)
            i += 1

    print(species, i)
    # return

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
    files = glob.glob(psm_dir + species + '/*.csv')
    print(files)
    exit(0)
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
    # return

    
#    mat = np.array(mat)
    
#    np.savetxt(data_dir+species+'_rep.csv', mat)
    json.dump(matches, open(data_dir+species+'_rep.json', 'w'))



#%%
for species in species_list:
    print(species)
    extract_mat(species)
    #extract_rep_mat(species)

