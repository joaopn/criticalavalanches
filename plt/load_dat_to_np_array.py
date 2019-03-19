import os
import numpy as np
from itertools import product
import sys
import glob
import mrestimator as mre
import matplotlib.pyplot as plt
from time import gmtime, strftime

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

num_elec = 100
num_neur_list = [144000,72000,36000,18000]
h_list = [1.1e-4]
m_list = [0.9]

num_neur=num_neur_list[-1]
m = m_list[-1]
current_dir = glob.glob(f"../dat/N{num_neur:06d}_m{m:.3f}*/*")[0]
current_dir = os.path.dirname(os.path.abspath(current_dir))
electrodes = []
for idx in range(0, num_elec):
    electrodes.append(np.loadtxt(f"{current_dir}/elec{idx:03d}.gz"))
elec_realisation = np.array(electrodes)

