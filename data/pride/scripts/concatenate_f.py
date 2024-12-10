import sys
import glob


assert len(sys.argv) == 2, f"Usage: python {sys.argv[0]} <species_dir>"
   
fs = glob.glob(f'{sys.argv[1]}/*.mgf')
#print(fs)

with open(f'{sys.argv[1]}/{sys.argv[1].replace("/","")}.mgf','w') as f_w:
    for f in fs:
        with open(f,'r') as f_r:
            f_w.writelines(f_r.readlines())
