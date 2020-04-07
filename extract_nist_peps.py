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
        'drosophila',
        'e_coli',
        'human',
        'mouse',
        'yeast',
    ]

#%%
nist_dir = 'test_search/nist/'
data_dir = 'test_search/matdata_nist/'

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
    if item:
        yield item


#%%
def extract_peps(species):
    mspfile = tarfile.open(nist_dir+species+'_consensus_final_true_lib.tar.gz', 'r|gz')
    tarinfo = mspfile.next()
    mspfile = mspfile.extractfile(tarinfo)
    
    msplist = (l.decode("utf-8").rstrip() for l in mspfile)

    pepfile = open(data_dir+species+'_peps.csv', 'w')
    pepfile.write('Pep')
    pepfile.write('\t')
    pepfile.write('Prot')
    pepfile.write('\n')

    #writer = csv.DictWriter(pepfile, fieldnames=['Name', 'Protein'])
    n = 0
    for item in parse_msp(msplist):
        if 'Name' not in item:
            print(item)
            break
        pep = item['Name'].split('/')[0]
        prot = item['Props']['Protein'].split(' ')[0]
        pepfile.write(pep)
        pepfile.write('\t')
        pepfile.write(prot)
        pepfile.write('\n')
        n += 1

    print('number of peps:', n)
    pepfile.close()

#%%
for species in species_list:
    print(species)
    extract_peps(species)



