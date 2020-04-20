#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:16:38 2019

@author: yisupeng
"""

import scipy.io as sio

latex_output = open('latex_snippet.txt', 'w')

#%%

#data_source = 'PRIDE'
data_source = 'HeLa'
#data_source = 'NIST'

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
    result_dir = 'test_search/est_results_pride/'
    data_source = 'PRIDE'
    figure_dir = 'figures/'

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
    result_dir = 'test_search/est_results_hela/'
    data_source = ''
    figure_dir = 'figures/'

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
    figure_dir = 'figures/nist/'

method_list = [
        '_1s2ca',
        '_1s2c',
        '_2s3ci',
    ]

#%%
subfigure_width = 0.5
template = open('template_subfigure.txt').read()

#%%
import csv
thres_map = csv.DictReader(open(result_dir+'thresholds.txt'), delimiter='\t')
thres_map = { row['species']+row['method'] : (float(row['thres']), float(row['thres_cor']), float(row['thres2'])) for row in thres_map}
#%%
def list_params(params):
    return ", ".join(["$%s=%.2f$"%(k,v) for k,v in params.items()])
def gen_latex(species, method):
    global num_S, num_D, model, thres, thres_cor, thres2, params, ll, mean, sdcdf
    global thres_true
    num_S = method[1]
    num_D = method[3]
    type_D = 's'
    if len(method) > 5 and method[5] != 'i':
        type_D = method[5]
#    print(method, num_S, num_D, type_D)
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
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)
    
    replace_var('figure_dir')
    replace_var('subfigure_width')
#    replace_var('data_source')
    replace_var('species')
    replace_var('method')
    replace_var('num_S')
    replace_var('num_D')
    replace_var('model')
    
    params = {}
    paramfile = result_dir+species+'/params/'+method+'.mat'
    
    if method in shantanu_m_l:
        content.replace('!params, ', '')
    else:
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
#        print(method, species)
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
#        print(params)
        
        replace_var('params')
    
    llfile = result_dir+species+'/ll/'+method+'.mat'
    sdcdffile = result_dir+species+'/sdcdf/'+method+'.mat'
    meanfile = result_dir+species+'/mean/'+method+'.mat'

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
    
    def get_val(s):
        return s.strip('$').split('=')[1]
    if method == '_1s2ca':
        print(thres, end='\t')
    print(get_val(sdcdf), get_val(ll.split(', ')[0]), thres2, sep='\t', end='\t')
    if method == '_2s3ci':
        print(get_val(ll.split(', ')[1]), end='')
    
    
    replace_var('thres2')
    replace_var('thres_cor')
    if data_source == 'NIST':
        thres_true = open(result_dir+species+'/true_thres.txt')
        thres_true = float(thres_true.read().rstrip())
        thres_true = "%.2f"%thres_true
        replace_var('thres_true')
    else:
        content = content.replace('$\\tau_{\\text{TRUE}}=!thres_true$, ', '')
    replace_var('thres')
    replace_var('ll')
    replace_var('mean')
#    print(content)
    latex_output.write(content)
    latex_output.write("\n")
    
#gen_latex('H.sapiens2', '_2s3ci')

#%%
for species in species_list:
    print(species, end='\t')
    for method in method_list:
        gen_latex(species, method)
    print('')
    latex_output.write("\n%%===================================================================\n")
latex_output.close()