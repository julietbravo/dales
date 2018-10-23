import glob
import os

all_files = glob.glob('*')
exclude   = ['namoptions.001','prof.inp.001','clean.py','lscale.inp.001']

for f in all_files:
    if f not in exclude:
        os.remove(f)
