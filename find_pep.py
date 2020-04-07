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
from Bio import SeqIO


#%%

species_list = [
        #'c_elegans',
        #'drosophila',
        #'e_coli',
        #'human',
        'mouse',
        'yeast',
    ]

#%%
res_dir = 'test_search/est_results_nist/'

data_dir = 'test_search/matdata_nist/'

#db_dir = '../ms/db/nist/'
db_dir = '../ms/db/'

seqfile_map = {
        'c_elegans': 'uniprot-taxonomy_6239.fasta',
        'drosophila': 'uniprot-taxonomy_7227.fasta',
        'e_coli': 'uniprot-proteome_UP000000625.fasta',
        'human': 'uniprot-taxonomy_9606.fasta',
        'mouse': 'uniprot-taxonomy_10090.fasta',
        'yeast': 'uniprot-taxonomy_559292.fasta',
    }

#%%
def get_prots(species):
    seqfile = db_dir + seqfile_map[species]
    recs = SeqIO.parse(seqfile, "fasta")
    # for i, rec in enumerate(recs):
    #     if i >= 5:
    #         break
    #     print(rec.id, rec.seq)
    prots = {rec.id: rec.seq for rec in recs}
    return prots

#%%
def find_pep_in_seq(seq, pep):
    return pep in seq

def find_pep_in_prots(prots, pep):
    for seq in prots.values():
        if pep in seq:
            return True
    return False

#%%
def find_peptides_in_db(species):
    species_dir = res_dir + species + '/'

    matches = []
    prots = get_prots(species)

    peps = open(data_dir+species+'_peps.csv')
    peps_csv = csv.DictReader(peps, delimiter='\t')

    match_csv = csv.writer(open(species_dir + 'matches.csv', 'w'))
    for row in peps_csv:
        pep = row['Pep']
        prot = row['Prot']

        if prot in prots:
            m = find_pep_in_seq(prots[prot], pep)
            rec = [m, True, prot]
            #print('finding pep:', pep, 'in prot:', prot, prots[prot], m)
        else:
            print(prot)
            m = find_pep_in_prots(prots, pep)
            rec = [m, False, prot]
            #print('finding pep:', pep, 'from all', m)
        matches.append(rec)

        match_csv.writerow(rec)

    return matches

#%%
for species in species_list:
    matches = find_peptides_in_db(species)



