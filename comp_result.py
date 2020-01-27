#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 12:52:41 2019
@author: yisupeng
"""

import csv
import numpy as np

m1 = open('test_search/isb02_t1.tsv')
m1 = csv.DictReader(m1, delimiter='\t')

n = 0
scores1 = []

msgfscore_eval_curv = []

for row in m1:
#    print(row)
    if row['Protein'][0:3] == 'REV':
        continue
    qval = float(row['EValue'])
#    if qval > 0.1:
#        break
    s = float(row['EValue'])
    scores1.append((n, s))
    ms = float(row['MSGFScore'])
    msgfscore_eval_curv.append((s,ms))
    n += 1

#%% plot
import matplotlib.pyplot as plt
msgfscore_eval_curv = np.array(msgfscore_eval_curv)
plt.plot(np.log(msgfscore_eval_curv[:,0]),msgfscore_eval_curv[:,1])
#%%
m2 = open('test_search/isb02_t3_nodecoy.tsv')
m2 = csv.DictReader(m2, delimiter='\t')

prevscan = 0
scores2 = []
n2 = 0
for row in m2:
#    print(row)
    if row['Protein'][0:3] == 'REV':
        continue
    scan = row['ScanNum']
    if scan == prevscan:
        continue
    prevscan = scan
#    print(scan)
#    qval = float(row['QValue'])
#    if qval > 0.1:
#        break
    s = float(row['EValue'])
    scores2.append((n2, s))
    n2 += 1