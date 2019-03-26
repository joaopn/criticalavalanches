import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# for calls eg. with qsub -t 1-12500
start_id = 1

# l_neur = np.logspace(1,3,7)*256
l_neur = [40000]
# l_m = [0.9]
l_m = np.logspace(np.log10(1.006), np.log10(1.06), 10)
# l_h = [1.1e-4]
l_h = np.linspace(1e-7, 2.8e-5, 6)
l_T = [5e6]
l_rep = range(0,3)
arg_list = product(l_neur, l_m, l_h, l_T, l_rep)

count = 0
with open('./find_par.tsv', 'w') as f:

    fmt = "-s %06d -N %06d -m %.3e -h %.3e -T %.3e -r %02d -o %s -c 1e4 -dn 50 -de 4\n"
    f.write("#" + fmt)

    for i in arg_list:
        i = (start_id,) + i + (f"./dat/crit/m{i[1]:.3f}_{i[2]:.3e}_r{i[4]:d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
