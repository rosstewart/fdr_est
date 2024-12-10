'''
author: Ross Stewart
created: Oct 27 2024
'''

import subprocess
import os
import sys

assert len(sys.argv) == 2, f"Usage: python {sys.argv[0]} <dataset_name>"

datasets = ["pride","nist"]
assert sys.argv[1] in datasets, f"Error: {sys.argv[1]} not in {datasets}"
dataset = sys.argv[1]

if dataset == "pride":
    species_list = 'a_thaliana	d_melanogaster	e_coli		h_sapiens	h_sapiens2	m_musculus	m_musculus2	m_musculus3	s_cerevisiae	s_cerevisiae2'.split()
    instruments = [0,1,0,1,1,3,3,3,3,3]
    fragmentation_methods = [0,0,0,0,0,3,3,3,0,0]
elif dataset == "nist":
    species_list = ['c_elegans','h_sapiens','m_musculus','s_cerevisiae','h_sapiens1','h_sapiens2']
    instruments = [0,0,0,0,0,0]
    fragmentation_methods = [0,0,0,0,0,0]

path_to_script = os.path.realpath(__file__).replace(sys.argv[0].split("/")[-1],"")
ppm = 25

for i,species in enumerate(species_list):
    if not os.path.exists(f'{path_to_script}data/{dataset}/{species}/{species}.mzid'):
        print(f'starting {species}',flush=True)
        subprocess.Popen(f'java -Xmx8G -jar {path_to_script}msgf/MSGFPlus.jar -s {path_to_script}data/{dataset}/{species}/{species}.mgf -d {path_to_script}data/reference_dbs/{species.replace("1","").replace("2","").replace("3","")}.fasta -n 10 -t {ppm}ppm -inst {instruments[i]} -m {fragmentation_methods[i]}'.split())
