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

num_neur=10000
for file in glob.glob(f"../dat/N{num_neur:06d}*/*"):
    print(file)
