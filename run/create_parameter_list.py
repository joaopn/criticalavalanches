import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# for calls eg. with qsub -t 1-12500
start_id = 1

# l_neur = np.logspace(1,3,7)*256
l_neur = [256000]
l_m = [0.9]
l_h = [1.1e-4]
l_T = [5e5]
l_rep = range(0,100)
arg_list = product(l_neur, l_m, l_h, l_T, l_rep)

count = 0
with open('./parameter_list.tsv', 'w') as f:

    fmt = "-s %06d -N %06d -m %.3e -h %.3e -T %.3e -r %02d -o %s\n"
    f.write("#" + fmt)

    for i in arg_list:
        i = (start_id,) + i + (f"./dat/{start_id:06d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
