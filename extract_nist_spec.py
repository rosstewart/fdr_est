#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 14:22:05 2020

@author: yisupeng
"""

import tarfile

#%%

#species = 'human'
#species = 'mouse'
#species = 'drosophila'
#species = 'e_coli'
#species = 'c_elegans'
#species = 'yeast'

species_list = [
        'c_elegans',
        'h_sapiens',
        'm_musculus',
        's_cerevisiae',
]

#libmap = {
#        'human_hcd': 'human_hcd_selected.msp.tar.gz',
#        'mouse_hcd': 'cptac2_mouse_hcd_selected.msp.tar.gz',
#    }

#data_dir = 'test_search/nist_hcd/'
data_dir = 'data/nist/'

#mspfile = open('nist/human_consensus_final_true_lib.msp')
#mspfile = tarfile.open('nist/'+species+'_consensus_final_true_lib.tar.gz', 'r|gz')

#%%
def split_comments(comment):
    ret = ''
    quote = False
    for c in comment:
        if c == ' ' and not quote:
            if '=' in ret:
                yield ret
            else:
                print(ret)
            ret = ''
        elif c == '"':
            quote = not quote
        else:
            ret += c
            

#%%

items = []

def parse_msp(msplist):
    peaks = []
    item = {}
    for l in msplist:
        if not l:
            item['Peaks'] = peaks
            
            props = {k:v for k,v in map(lambda x: x.split('=', 1), split_comments(item['Comment']))}
            item['Props'] = props
#            print(len(peaks))
#            items.append(item)
            yield item
    #        if len(items) >= 5:
    #            break
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
#    yield item
    


#%%
        
def write_mgf(f, item):
    props = item['Props']
    f.write('BEGIN IONS\n')
    f.write('PEPMASS=%s\n'%props['Parent'])
    f.write('CHARGE=%s+\n'%props['Charge'])
    f.write('TITLE=%s\n'%item['Name'])
    for peak in item['Peaks']:
        f.write('%s %s\n'%(peak[0], peak[1]))
    
    f.write('END IONS\n\n')

for species in species_list:
    mspfile_str = f'{data_dir}{species}/{species}.msp'
    #mspfile = tarfile.open(mspfile, 'r|gz')
    #tarinfo = mspfile.next()
    #mspfile = mspfile.extractfile(tarinfo)
    
    mspfile = open(mspfile_str, 'r')
    msplist = (l.rstrip() for l in mspfile)
    
    #for i in range(10):
    #    print(next(msplist))

    mgffile = open(f'{data_dir}{species}/{species}.mgf', 'w')
    n = 0
    for item in parse_msp(msplist):
        if not item:
            continue
        write_mgf(mgffile, item)
        n += 1
        print(n)

    mspfile.close()
    mgffile.close()    

