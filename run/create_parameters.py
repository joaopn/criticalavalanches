import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# seeds will  not match if you rerun this, data has been copied.
start_id = 0

# settings for reverberating regime without cc
l_neur = [160000]
l_m = [   0.0,    0.9,   0.98,  0.996,  0.999]
l_h = [2.0e-3, 2.0e-4, 4.0e-5, 8.0e-6, 2.0e-6]
l_T = [1e6]
l_de = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
l_rep = range(0,50)

arg_list = product(l_neur, l_m, l_h, l_de, l_T, l_rep)

count = 0
with open('./parameters.tsv', 'w') as f:

    fmt = "-s %06d -N %06d -m %.3e -h %.3e -de %02d -dn 50 -T %.3e -r %02d -o %s -g 6.0 -c 1e5 -t 1e4\n"
    f.write("#" + fmt)

    for i in arg_list:
        if (l_m.index(i[1]) != l_h.index(i[2])): continue
        i = (start_id,) + i + (f"./dat/m{i[1]:.5f}_h{i[2]:.3e}_d{i[3]:02d}_r{i[-1]:02d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
