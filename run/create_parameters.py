import os
import numpy as np
from itertools import product

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

# seeds will not match if you rerun this, data has been copied.
start_id = 0

# settings for reverberating regime without cc
l_exe = ["./exe/cc", "./exe/cc_orlandi"]
l_folder = ["gauss", "orlandi"]
l_neur = [160000]
l_m = [   0.9,   0.98]
l_h = [2.0e-4, 4.0e-5]
l_T = [1e6]
l_de = [2, 4]
l_ga = [-1.0, -1.25, -1.5, -1.75, -2.0]
l_rep = range(0,50)

arg_list = product(l_exe, l_neur, l_m, l_h, l_de, l_T, l_ga, l_rep)

count = 0
with open('./parameters.tsv', 'w') as f:

    fmt = "%s -s %06d -N %06d -m %.3e -h %.3e -de %02d -dn 50 -T %.3e -ga %.2f -r %02d -o %s -g 6.0 -c 1e5 -t 1e4\n"
    f.write("#" + fmt)

    for i in arg_list:
        if (l_m.index(i[2]) != l_h.index(i[3])): continue
        folder = l_folder[l_exe.index(i[0])]
        i = (i[0],) + (start_id,) + i[1:] + (f"./dat/{folder}/m{i[2]:.5f}_h{i[3]:.3e}_de{i[4]:02d}_ga{i[6]:.2f}_r{i[-1]:02d}.hdf5",)
        start_id += 1
        count += 1
        f.write(fmt % i)


print(f"number of argument combinations: {count}")
