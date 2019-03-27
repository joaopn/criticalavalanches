import h5py
import numpy as np
from time import gmtime, strftime
from itertools import product
import sys
import os
import matplotlib.pyplot as plt
import mrestimator as mre
import numpy as np

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))


l_neur = [160000]
l_m = [0.9, 0.98, 1.0]
l_h = np.logspace(np.log10(1.e-7), np.log10(1e-5), 9)
l_T = [1e6]
l_rep = range(0,3)[0]
for m in l_m:
    for h in l_h:
        for r in l_rep:
            path = f"./dat/coalcomp/m{m:.5f}_h{h:.3e}_r{r:d}.hdf5"
            print(path)
            try:
                f = h5py.File(path, "r")

                dat = f["/data/activity"][5000:]
                act = np.mean(dat)/f["meta/num_neur"][0]/2.*1000.
                m = f["meta/m_micro"][0]
                h = f["meta/h_prob"][0]
                print(f"\tact: {act:.4} Hz m: {m:5} h: {h:5e}")

                mre.full_analysis(dat, 2, 2000, dtunit='ms', fitfuncs=['e', 'eo'],
                    title="act: {act:.4} Hz m: {m:5} h: {h:5e}")
            except:
                pass
