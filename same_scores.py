#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:16:43 2019

@author: yisupeng
"""


import csv
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


#species = 'M.musculus'
#species = 'M.musculus2'
#species = 'M.musculus3'
#species = 'H.sapiens2'
#species = 'H.sapiens3'
#species = 'C.elegans'
#species = 'D.melanogaster'
#species = 'S.cerevisiae'
species = 'S.cerevisiae2'
#species = 'S.cerevisiae3'
#species = 'E.coli'
#species = 'A.thaliana'

#%%
matches = open('test_search/pride/%s_nod.tsv'%species)
matches_csv = csv.DictReader(matches, delimiter='\t')


class Entity:
    def __init__(self, raw):
        self.raw = raw

    def __getattr__(self, item):
        return self.raw[item]

    def __repr__(self):
        return str(self.raw)

spec_scores = defaultdict(list)
peps = defaultdict(list)
proteins = defaultdict(list)
for row in matches_csv:
    psm = Entity(row)
    spec_scores[psm.SpecID].append(psm.SpecEValue)
    peps[psm.SpecID].append(psm.Peptide)
    proteins[psm.SpecID].append(psm.Protein)

rawmat = spec_scores.values()
slen = [len(sl) for sl in rawmat]
slen = np.array(slen)
mat = np.array(list(rawmat))
mat = mat[slen>=3]
mat3 = [row[0:3] for row in mat]
mat3 = np.array(mat3)

s1 = mat3[:,0]
s2 = mat3[:,1]
s3 = mat3[:,2]

#%%
flags = s1 == s2

#for i, spec in enumerate(spec_scores.keys()):
#    if i < len(flags) and flags[i]:
#        print(flags[i])
ls = ''
for spec, scores in spec_scores.items():
    if len(scores) >= 3 and scores[0] == scores[1]:
#        print("%s,%s,%s,%s,%s"%tuple(spec_scores[spec][0:1] + peps[spec][0:2] + proteins[spec][0:2]))
#        ls.append("%s,%s,%s,%s,%s"%tuple(spec_scores[spec][0:1] + peps[spec][0:2] + proteins[spec][0:2]))
        ls += ("%s\t%s\t%s\t%s\t%s\n"%tuple(spec_scores[spec][0:1] + peps[spec][0:2] + proteins[spec][0:2]))

        

