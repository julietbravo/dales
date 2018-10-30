import glob
import os

all_files = glob.glob('*')
exclude   = ['namoptions.001','prof.inp.001','lscale.inp.001','run.PBS']

for f in all_files:
    if f not in exclude and '.py' not in f:
        os.remove(f)
