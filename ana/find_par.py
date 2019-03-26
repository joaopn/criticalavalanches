import h5py
import numpy as np
from time import gmtime, strftime
import sys
import os
import matplotlib.pyplot as plt
import mrestimator as mre
import numpy as np


# https://stackoverflow.com/questions/41231678/obtaining-a-exclusive-lock-when-writing-to-an-hdf5-file

for name in ['non', 'rev', 'crit1', 'e5']:
    f = h5py.File(f"/Users/paul/owncloud/mpi/simulation/critical_avalanches/gh_test_src/prelim/test_{name}.hdf5", "r")

    dat = f["/data/activity"][5000:]
    plt.plot(dat)
    act = np.mean(dat)/f["meta/num_neur"][0]/2.*1000.
    m = f["meta/m_micro"][0]
    h = f["meta/h_prob"][0]
    print(f"{name} act: {act:.4} Hz m: {m:5} h: {h:5e}")



    mre.full_analysis(dat, 2, 2000, dtunit='ms', fitfuncs=['e', 'eo'], title=name)
