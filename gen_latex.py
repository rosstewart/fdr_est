#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:16:38 2019

@author: yisupeng
"""

import scipy.io as sio

latex_output = open('latex_snippet.txt', 'w')

#%%
method_list = [
        '_1s2c',
        '_1s2ca',
#        '_2s2ci',
        '_2s3ci',
        '_2s3ct',
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

#%%
subfigure_width = 0.5
template = open('subfigure_template').read()

#%%
import csv
thres_map = csv.DictReader(open('test_search/est_results/thresholds.txt'), delimiter='\t')
thres_map = { row['species']+row['method'] : (float(row['thres']), float(row['thres_cor']), float(row['thres2'])) for row in thres_map}
#%%
def list_params(params):
    return ", ".join(["$%s=%.2f$"%(k,v) for k,v in params.items()])
def gen_latex(species, method):
    global num_S, num_D, model, thres, thres_cor, thres2, params, ll, mean, sdcdf
    num_S = method[1]
    num_D = method[3]
    type_D = 's'
    if len(method) > 5 and method[5] != 'i':
        type_D = method[5]
    print(method, num_S, num_D, type_D)
    if type_D == 's':
        model = 'skew normal'
    elif type_D == 't':
        model = 'skew normal, $\mu$ truncated'
    elif type_D == 'a':
        model = 'gamma \& gaussian'
    
    thres, thres_cor, thres2 = thres_map[species+method]
    thres = "%.2f"%thres
    thres_cor = "%.2f"%thres_cor
    thres2 = "%.2f"%thres2
    
    content = template
    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        content = content.replace('!'+k, v)
    
    replace_var('subfigure_width')
    replace_var('species')
    replace_var('method')
    replace_var('num_S')
    replace_var('num_D')
    replace_var('model')
    
    params = {}
    paramfile = 'test_search/est_results/'+species+'/params/'+method+'.mat'
    theta = sio.loadmat(paramfile)['theta']
#    if 'theta_i' in theta:
    
    params['\\alpha'] = theta['alpha'][0,0][0,0]
    if int(method[3]) in [3,]:
        params['\\beta'] = theta['beta'][0,0][0,0]
    
    def put_param(p, s):
        v = theta['theta_'+s][0,0][p][0,0][0,0]
        if p == 'u':
            p = '\mu'
        params['%s_{%s}'%(p, s.upper())] = v
    print(method, species)
    put_param('u', 'c')
    if type_D in ['s', 't']:
        put_param('u', 'i')
    elif type_D == 'a':
        put_param('a', 'i')
    
    if int(num_D) > 2:
        put_param('u', 'i2')
    if int(num_D) > 3:
        put_param('u', 'i3')
    params = list_params(params)
    print(params)
    
    llfile = 'test_search/est_results/'+species+'/ll/'+method+'.mat'
    sdcdffile = 'test_search/est_results/'+species+'/sdcdf/'+method+'.mat'
    meanfile = 'test_search/est_results/'+species+'/mean/'+method+'.mat'

    lls = sio.loadmat(llfile)
    ll = []
    for k,v in lls.items():
        if 'll' in k:
            k = k.split('_')
            k = '%s_{%s}'%tuple(k)
            ll.append('$%s=%.4f$'%(k,v[0,0]))
    ll = ', '.join(ll)
    
    sdcdf = sio.loadmat(sdcdffile)
    sdcdf = sdcdf['sdcdf'][0,0]
    sdcdf = '$\\delta_\\text{cdf} = %.2f$' % sdcdf 
    replace_var('sdcdf')
    
    meanf = sio.loadmat(meanfile)
    mean = []
    for k,v in meanf.items():
        if 'E' in k:
            k = k.split('_')
            k = '%s_{%s}'%(k[0], k[1].upper())
            mean.append('$%s=%.4f$'%(k,v[0,0]))
    mean = ', '.join(mean)
    replace_var('thres2')
    replace_var('thres_cor')
    replace_var('thres')
    replace_var('params')
    replace_var('ll')
    replace_var('mean')
#    print(content)
    latex_output.write(content)
    latex_output.write("\n")
    
#gen_latex('H.sapiens2', '_2s3ci')

#%%
for species in species_list:
    for method in method_list:
        gen_latex(species, method)
    latex_output.write("\n%%===================================================================\n")
latex_output.close()