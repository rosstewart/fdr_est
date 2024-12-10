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
from matplotlib.ticker import FormatStrFormatter
from statistics import mean, stdev
import json
import os

#%%
method_list = [
        '_1s2ca',
        '_1s2c',
#        '_2s2ci',
        '_2s3ci',
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
        'S.cerevisiae3',
    ]
species_list = [f'synthetic_{i}' for i in range(1,31)]

#species_list = [
#        'HeLa01ng',
#        'HeLa1ng',
#        'HeLa10ng',
#        'HeLa50ng',
#        'HeLa100ng'
#    ]

#species_list = [
#        'E.coli',
#    ]

#results_dir = 'test_search/est_results_8inits/'
#results_dir = 'test_search/est_results_full/'
#results_dir = 'test_search/est_results_part/'
results_dir = 'synthetic/est_results_full/'

#%%

#ticks = []
#pos = 1
#thres_mat = []
#dcdf_mat = []
#for species in species_list:
#    for method in method_list:
#        if method == '_1s2c':
#            ticks.append(species+'_1S2D skew normal')
#        elif method == '_1s2ca':
#            ticks.append(species+'_1S2D gamma & gaussion')
#        elif method == '_2s3ci':
#            ticks.append(species+'_2S3D skew normal')
#        thres_list = boxplot_thres(species, method)
#        dcdf_list = boxplot_dcdf(species, method)
##        ax.boxplot(thres_list, positions=[pos])
#        thres_mat.append(thres_list)
#        dcdf_mat.append(dcdf_list)
#        print(species, method)
#        pos += 1

#%%
i = 0
for method in method_list:
    if method == '_1s2c':
        algo = ('1S2D skew normal')
    elif method == '_1s2ca':
        algo = ('1S2D gamma & gaussion')
    elif method == '_2s3ci':
        algo = ('2S3D skew normal')
    print(',\t'+algo, end='')
print('')



#%%
bootstrap_dir = results_dir + '/json/bootstrap/'

#tda_thres = json.load(open(bootstrap_dir+'tda.json'))
#tda_thres = {obj['ds']:obj['t1p'] for obj in tda_thres}

def boxplot_species(species):
    thres_1 = json.load(open(bootstrap_dir+species+'.json'))
    thres_1 = np.array([obj['t1p'] for obj in thres_1])
    print(thres_1.shape)
   # thres_2 = np.array(tda_thres[species])
    
    thres_mat = np.vstack([thres_1]).T
    ticks = ['TDA', 'GG', '1SMix', '2SMix']
    
    fig = plt.figure(figsize=[4,3])
#    plt.subplots_adjust(left=0.2, bottom=0.25)
    ax = plt.axes()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    bp = ax.boxplot(thres_mat)
    ax.set_xticklabels(ticks, rotation=45, ha='right')
#    plt.ylabel('Threshold')
#    plt.xlabel('Methods')
    ax1 = ax
    
    box_colors = [[0, 0.4470, 0.7410],
                  [0.8500, 0.3250, 0.0980],
                  [0.9290, 0.6940, 0.1250],
                  [0.4940, 0.1840, 0.5560],
                  [0.4660, 0.6740, 0.1880],]
    
    num_boxes = thres_mat.shape[1]
    medians = np.empty(num_boxes)
    for i in range(num_boxes):
        print(num_boxes, i)
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        # Alternate between Dark Khaki and Royal Blue
        ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % len(box_colors)]))
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
#        ax1.plot(np.average(med.get_xdata()), np.average(thres_mat[i]),
#                 color='w', marker='*', markeredgecolor='k')
    
    plt.tight_layout()
    
    plt.savefig(results_dir + '/bootstrap/'+species+'.png', dpi=320)
    plt.savefig(results_dir + '/bootstrap/'+species+'.eps', dpi=320)
    plt.show()

#%%
for species in species_list:
    boxplot_species(species)
    

#%% boxplot delta_CDF
#ax = plt.axes()
#bp = ax.boxplot(dcdf_mat)
#ax.set_xticklabels(ticks, rotation=45, ha='right')
#ax1 = ax
#
#box_colors = ['darkkhaki', 'royalblue', 'darkred']
#num_boxes = len(dcdf_mat)
#medians = np.empty(num_boxes)
#for i in range(num_boxes):
#    box = bp['boxes'][i]
#    boxX = []
#    boxY = []
#    for j in range(5):
#        boxX.append(box.get_xdata()[j])
#        boxY.append(box.get_ydata()[j])
#    box_coords = np.column_stack([boxX, boxY])
#    # Alternate between Dark Khaki and Royal Blue
#    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 3]))
#    # Now draw the median lines back over what we just filled in
#    med = bp['medians'][i]
#    medianX = []
#    medianY = []
#    for j in range(2):
#        medianX.append(med.get_xdata()[j])
#        medianY.append(med.get_ydata()[j])
#        ax1.plot(medianX, medianY, 'k')
#    medians[i] = medianY[0]
#    # Finally, overplot the sample averages, with horizontal alignment
#    # in the center of each box
#    ax1.plot(np.average(med.get_xdata()), np.average(dcdf_mat[i]),
#             color='w', marker='*', markeredgecolor='k')
#    
#plt.show()

