#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:16:38 2019

@author: yisupeng
"""

import scipy.io as sio

latex_output = open('latex_slides_snippet.txt', 'w')

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

elif data_source == 'NIST':
    species_list = [
            'c_elegans',
            'drosophila',
            'e_coli',
            'human',
            'mouse'
        ]
    result_dir = 'test_search/est_results_nist/'
    data_source = 'NIST'


#%%
subfigure_width = 0.5
if data_source == 'NIST':
    template = open('template_slide_subflot_truethres.txt').read()
else:
    template = open('template_slide_subflot.txt').read()

#%%
import csv
thres_map = csv.DictReader(open(result_dir+'thresholds.txt'), delimiter='\t')
thres_map = { row['species']+row['method'] : (float(row['thres']), float(row['thres_cor']), float(row['thres2'])) for row in thres_map}
print(thres_map.keys())

#%%


#%%
def list_params(params):
    return ", ".join(["$%s=%.2f$"%(k,v) for k,v in params.items()])

def gen_latex(species, method):
    global title
    global num_S, num_D, model, dep, thres, thres_true, thres_cor, thres2, params, ll, mean, sdcdf
    num_S = method[1]
    num_D = method[3]
    if method == 'SNMax1':
        num_S = str(1)
        num_D = str(2)
        dep = 'Dependent'
    else:
        dep = 'Marginal'
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
    
    if method in shantanu_m_l:
        model = 'skew normal, Max'
    
    thres, thres_cor, thres2 = thres_map[species+method]
    thres = "%.2f"%thres
    thres_cor = "%.2f"%thres_cor
    thres2 = "%.2f"%thres2
    
    title = data_source + ' \emph{' + species.replace('_', '.~').capitalize() + '}'
    
    content = ''
    content += '\\begin{frame}{!title}\n'
    content += template
    content += '\\end{frame}\n'
    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)
    
    replace_var('title')
    replace_var('subfigure_width')
    replace_var('data_source')
    replace_var('species')
    replace_var('method')
    replace_var('dep')
    replace_var('num_S')
    replace_var('num_D')
    replace_var('model')

    params = {}
    def put_param(p, s):
        v = theta['theta_'+s][0,0][p][0,0][0,0]
        if p == 'u':
            p = '\mu'
        params['%s_{%s}'%(p, s.upper())] = v
        
    if method in shantanu_m_l:
#        content.replace('!params, ', '')
#        print(content.find(', !mean'))
        content = content.replace(', !mean', '')
        print(content)
        
        estjson = json.load(open('test_search/shantanu/json/'+species+'.json'))
        estjson = {res['algo']:res for res in estjson}
        estres = estjson[method]
        alpha = estres['pi_C']

        params['\\alpha'] = alpha
    else:
        paramfile = result_dir+species+'/params/'+method+'.mat'
        theta = sio.loadmat(paramfile)['theta']
    #    if 'theta_i' in theta:
        
        params['\\alpha'] = theta['alpha'][0,0][0,0]
        if int(method[3]) in [3,]:
            params['\\beta'] = theta['beta'][0,0][0,0]
        
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
        
        meanfile = result_dir+species+'/mean/'+method+'.mat'
        meanf = sio.loadmat(meanfile)
        mean = []
        for k,v in meanf.items():
            if 'E' in k:
                k = k.split('_')
                k = '%s_{%s}'%(k[0], k[1].upper())
                mean.append('$%s=%.4f$'%(k,v[0,0]))
        mean = ', '.join(mean)
            
        replace_var('mean')
    

    params = list_params(params)
    print(params)
    
    replace_var('params')
    
    
    llfile = result_dir+species+'/ll/'+method+'.mat'
    sdcdffile = result_dir+species+'/sdcdf/'+method+'.mat'

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
    
    
    
    replace_var('thres2')
    replace_var('thres_cor')
    if data_source == 'NIST':
        thres_true = open(result_dir+species+'/true_thres.txt')
        thres_true = float(thres_true.read().rstrip())
        thres_true = "%.2f"%thres_true
        replace_var('thres_true')
    else:
        content = content.replace(', !thres_true', '')
    replace_var('thres')
    replace_var('ll')
#    print(content)
    latex_output.write(content)
    latex_output.write("\n")

#%%

def gen_latex_for_charge(species, charge, method):
    global title
    global num_S, num_D, model, dep, thres, thres_true, thres_cor, thres2, params, ll, mean, sdcdf
    num_S = method[1]
    num_D = method[3]
    if method == 'SNMax1':
        num_S = str(1)
        num_D = str(2)
        dep = 'Dependent'
    else:
        dep = 'Marginal'
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
    
    if method in shantanu_m_l:
        model = 'skew normal, Max'
    
    title = data_source + ' \emph{' + species.replace('_', '.~').capitalize() + '}' \
        + ' Charge ' + str(charge) + '+'
    species = species + '.c' + str(charge)
    
    thres, thres_cor, thres2 = thres_map[species+method]
    thres = "%.2f"%thres
    thres_cor = "%.2f"%thres_cor
    thres2 = "%.2f"%thres2
    
    content = ''
    content += '\\begin{frame}{!title}\n'
    content += template
    content += '\\end{frame}\n'
    def replace_var(k):
        nonlocal content, species
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)
    
    replace_var('title')
    replace_var('subfigure_width')
    replace_var('data_source')
    replace_var('species')
    replace_var('method')
    replace_var('dep')
    replace_var('num_S')
    replace_var('num_D')
    replace_var('model')

    params = {}
    def put_param(p, s):
        v = theta['theta_'+s][0,0][p][0,0][0,0]
        if p == 'u':
            p = '\mu'
        params['%s_{%s}'%(p, s.upper())] = v
        
    if method in shantanu_m_l:
#        content.replace('!params, ', '')
#        print(content.find(', !mean'))
        content = content.replace(', !mean', '')
        print(content)
        
        estjson = json.load(open('test_search/shantanu/json/'+species+'.json'))
        estjson = {res['algo']:res for res in estjson}
        estres = estjson[method]
        alpha = estres['pi_C']

        params['\\alpha'] = alpha
    else:
        paramfile = result_dir+species+'/params/'+method+'.mat'
        theta = sio.loadmat(paramfile)['theta']
    #    if 'theta_i' in theta:
        
        params['\\alpha'] = theta['alpha'][0,0][0,0]
        if int(method[3]) in [3,]:
            params['\\beta'] = theta['beta'][0,0][0,0]
        
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
        
        meanfile = result_dir+species+'/mean/'+method+'.mat'
        meanf = sio.loadmat(meanfile)
        mean = []
        for k,v in meanf.items():
            if 'E' in k:
                k = k.split('_')
                k = '%s_{%s}'%(k[0], k[1].upper())
                mean.append('$%s=%.4f$'%(k,v[0,0]))
        mean = ', '.join(mean)
            
        replace_var('mean')
    

    params = list_params(params)
    print(params)
    
    replace_var('params')
    
    
    llfile = result_dir+species+'/ll/'+method+'.mat'
    sdcdffile = result_dir+species+'/sdcdf/'+method+'.mat'

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
    
    
    
    replace_var('thres2')
    replace_var('thres_cor')
    if data_source == 'NIST':
        thres_true = open(result_dir+species+'/true_thres.txt')
        thres_true = float(thres_true.read().rstrip())
        thres_true = "%.2f"%thres_true
        replace_var('thres_true')
    else:
        content = content.replace(', !thres_true', '')
    replace_var('thres')
    replace_var('ll')
#    print(content)
    latex_output.write(content)
    latex_output.write("\n")
    
#gen_latex('H.sapiens2', '_2s3ci')
#%%
def gen_latex_fdrcurv(species):
#    content = ''
#    content += '\\begin{frame}{!title}\n'
#    content += '\\begin{figure}\n'
#    content += '\\subfloat[linear scale]{'
#    content += '  \\includegraphics[width=.5\\textwidth]{figures/nist/{!species}/fdrcurv/fdrcmp.png}\n'
#    content += '}'
#    content += '\\subfloat[log scale]{'
#    content += '  \\includegraphics[width=.5\\textwidth]{figures/nist/{!species}/fdrcurv/fdrcmp_log.png}\n'
#    content += '}'
#    content += '  \\caption{Estimated FDR vs. True FDR}\n'
#    content += '\\end{figure}\n'
#    content += '\\end{frame}\n'

    content = ''
    content += '\\begin{frame}{!title}\n'
    content += '\\begin{figure}\n'
    content += '  \\includegraphics[width=.8\\textwidth]{figures/nist/{!species}/fdrcurv/fdrcmp.png}\n'
    content += '  \\caption{Estimated FDR vs. True FDR}\n'
    content += '\\end{figure}\n'
    content += '\\end{frame}\n'
    
    content += '\\begin{frame}{!title}\n'
    content += '\\begin{figure}\n'
    content += '  \\includegraphics[width=.8\\textwidth]{figures/nist/{!species}/fdrcurv/fdrcmp_log.png}\n'
    content += '  \\caption{Estimated FDR vs. True FDR in log scale}\n'
    content += '\\end{figure}\n'
    content += '\\end{frame}\n'
    
    title = data_source + ' \emph{' + species.replace('_', '.~').capitalize() + '}'

    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)

    replace_var('title')
    replace_var('data_source')
    replace_var('species')
    
    
    latex_output.write(content)
    latex_output.write("\n")

#%%
def gen_latex_tprcurv(species):
    content = ''
    content += '\\begin{frame}{!title}\n'
    content += '\\begin{figure}\n'
    content += '  \\includegraphics[width=\\textwidth]{figures/nist/{!species}/tprcurv.png}\n'
    content += '  \\caption{True Positive Rate at different scores}\n'
    content += '\\end{figure}\n'
    content += '\\end{frame}\n'
    
    title = data_source + ' \emph{' + species.replace('_', '.~').capitalize() + '}'

    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)

    replace_var('title')
    replace_var('data_source')
    replace_var('species')
    
    
    latex_output.write(content)
    latex_output.write("\n")

#%%
#bootstrap_result_dir = 'test_search/est_results_full/'

def gen_latex_bootstrap(species):
    content = ''
    content += '\\begin{frame}{!title}\n'
    content += '\\begin{figure}\n'
    content += '  \\includegraphics[width=\\textwidth]{figures/bootstrap/{!species}.png}\n'
    content += '  \\caption{1\\% FDR Threshold in Bootstraping}\n'
    content += '\\end{figure}\n'
    content += '\\end{frame}\n'
        
    title = data_source + ' \emph{' + species.replace('_', '.~').capitalize() + '}'

    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)
    
    replace_var('title')
    replace_var('data_source')
    replace_var('species')
    
    
    latex_output.write(content)
    latex_output.write("\n")

#%%
def gen_latex_charges(species):
    global title
    content = ''
    content += '\\begin{frame}{!title}\n'
    content += '\\begin{figure}\n'
    content += '  \\includegraphics[width=0.6\\textwidth]{figures/{!species}/charges_hist.png}\n'
    content += '  \\caption{Histogram for PSM of different charges}\n'
    content += '\\end{figure}\n'
    content += '\\end{frame}\n'
    
    title = species + ' Charges +2/+3'

    def replace_var(k):
        nonlocal content
        v = str(eval(k))
        if k == 'species':
            v = v.replace('_', '.')
        content = content.replace('!'+k, v)
    
    replace_var('title')
    replace_var('data_source')
    replace_var('species')
    
    
    latex_output.write(content)
    latex_output.write("\n")

#%%
bootstraps = {
        'HeLa01ng',
        'HeLa1ng',
        'HeLa10ng',
        'HeLa50ng',
        'HeLa100ng',
    }
for species in species_list:
    if data_source == '':
        gen_latex_charges(species)
    for method in method_list:
        gen_latex(species, method)
        if data_source == '':
            for charge in [2, 3]:
                gen_latex_for_charge(species, charge, method)
        
    if data_source == 'NIST':
        gen_latex_fdrcurv(species)
        gen_latex_tprcurv(species)
#        gen_latex_tprcurv_log(species)
    else:
        if species in bootstraps:
            gen_latex_bootstrap(species)
    
    latex_output.write("\n%%===================================================================\n")
latex_output.close()
