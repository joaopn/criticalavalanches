import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# for calls eg. with qsub -t 1-12500
start_id = 1

# l_neur = np.logspace(1,3,7)*256
l_neur = [160000]
# l_m = [0.9]
l_m = [0.92, 1.006, 1.03]
l_h = [1.8e-4, 3.0e-5, 2.2e-6]
# l_h = np.linspace(1e-7, 2.8e-5, 6)
l_T = [1e5]
l_rep = range(0,3)
arg_list = product(l_neur, l_m, l_h, l_T, l_rep)

count = 0
with open('./prelim.tsv', 'w') as f:

    fmt = "-s %06d -N %06d -m %.3e -h %.3e -T %.3e -r %02d -o %s -c 1e5 -dn 50 -de 4 -t 1e3\n"
    f.write("#" + fmt)

    for i in arg_list:
        if (i[1] == .92   and i[2] != 1.8e-4) \
        or (i[1] == 1.006 and i[2] != 3.0e-5) \
        or (i[1] == 1.03  and i[2] != 2.2e-6):
            continue
        i = (start_id,) + i + (f"./dat/prelim/m{i[1]:.5f}_h{i[2]:.3e}_r{i[4]:d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
