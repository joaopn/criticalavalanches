import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# seeds will not match if you rerun this, data has been copied.
start_id = 25000

# settings for reverberating regime without cc
l_exe = ["./exe/cc_random"]
l_folder = ["random_local", "random_nonlocal"]
l_neur = [160000]
l_m = [   0.9,   0.98]
l_h = [2.0e-4, 4.0e-5]
l_T = [1e6]
l_de = [4]
l_ga = [-1.0]
l_g  = [6.0, 0.0]
l_rep = range(0,50)

arg_list = product(l_exe, l_neur, l_m, l_h, l_de, l_T, l_ga, l_g, l_rep)

count = 0
with open('./parameters.tsv', 'w') as f:

    fmt = "%s -s %06d -N %06d -m %.3e -h %.3e -de %02d -dn 50 -T %.3e -ga %.2f -g %.2f -r %02d -o %s -c 1e5 -t 1e4\n"
    f.write("#" + fmt)

    for i in arg_list:
        if (l_m.index(i[2]) != l_h.index(i[3])): continue
        folder = l_folder[l_g.index(i[7])]
        i = (i[0],) + (start_id,) + i[1:] + (f"./dat/{folder}/m{i[2]:.5f}_h{i[3]:.3e}_de{i[4]:02d}_ga{i[6]:.2f}_r{i[-1]:02d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
