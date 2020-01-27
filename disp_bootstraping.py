#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 22:24:44 2020

@author: yisupeng
"""


import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from matplotlib.patches import Polygon

#%%
method_list = [
        '_1s2c',
        '_1s2ca',
#        '_2s2ci',
#        '_2s3ci',
#        '_2s3ct',
#        '_3s4ci'
    ]

species_list = [
        'A.thaliana',
#        'C.elegans',
        'D.melanogaster',
        'E.coli',
        'H.sapiens2',
        'H.sapiens3',
        'M.musculus',
        'M.musculus2',
        'M.musculus3',
#        'S.cerevisiae',
        'S.cerevisiae2',
#        'S.cerevisiae3',
    ]

#%%

nboots = 200
def boxplot_thres(species, method):
    thres_list = []
    for bootstrap_num in range(nboots):
        species_dir = 'test_search/est_results/'+species+'/'
        fdr_dir = species_dir + 'fdr/'
        bootstrap_dir = fdr_dir + method+'/bootstrap/'
        fdr_file = bootstrap_dir + str(bootstrap_num+1)+'.csv'
        
        fdr_est = np.genfromtxt(fdr_file, delimiter=',')
        
        eval_fdr_map = {}
            
        for fdr,n,s in fdr_est:
            eval_fdr_map[s] = fdr
        for fdr,n,s in fdr_est:
            if fdr > 0.01:
                break
            nmatches2 = n
            thres2 = s
        thres_list.append(thres2)
    return thres_list

#thres_list = boxplot_thres('S.cerevisiae3', '_1s2c')

ax = plt.axes()

ticks = []
pos = 1
thres_mat = []
for species in species_list:
    for method in method_list:
        if method == '_1s2c':
            ticks.append(species+'_1S2D skew normal')
        elif method == '_1s2ca':
            ticks.append(species+'_1S2D gamma & gaussion')
        thres_list = boxplot_thres(species, method)
#        ax.boxplot(thres_list, positions=[pos])
        thres_mat.append(thres_list)
        print(species, method)
        pos += 1


#%%
bp = ax.boxplot(thres_mat)
ax.set_xticklabels(ticks, rotation=45, ha='right')
ax1 = ax

box_colors = ['darkkhaki', 'royalblue']
num_boxes = len(thres_mat)
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    # Alternate between Dark Khaki and Royal Blue
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))
    # Now draw the median lines back over what we just filled in
    med = bp['medians'][i]
    medianX = []
    medianY = []
    for j in range(2):
        medianX.append(med.get_xdata()[j])
        medianY.append(med.get_ydata()[j])
        ax1.plot(medianX, medianY, 'k')
    medians[i] = medianY[0]
    # Finally, overplot the sample averages, with horizontal alignment
    # in the center of each box
    ax1.plot(np.average(med.get_xdata()), np.average(thres_mat[i]),
             color='w', marker='*', markeredgecolor='k')
    
plt.show()