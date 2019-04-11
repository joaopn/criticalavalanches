import os
import numpy as np
from itertools import product
import multiprocessing as mp

def run_cmd(m, h, de):

    if (l_m.index(m) != l_h.index(h)): return

    mask = f"'m{m:.5f}_h{h:.3e}_d{de:02d}'"

    # cmd = f"python ./run_analysis.py --mode threshold --bw_filter True --reps 50 --datafolder ../dat/ --datamask {mask}"
    # os.system(cmd)
    # print(f"{cmd} done")

    # cmd = f"python ./run_analysis.py --mode save_ps -b '1,2,4,8,16,32' --bw_filter True --datafolder ../dat/ --datamask {mask}"
    # os.system(cmd)
    # print(f"{cmd} done")

    cmd = f"python ./run_analysis.py --mode save_ps_alpha -b '1,2,4,8,16,32' --bw_filter True --datafolder ../dat/ --datamask {mask}"
    os.system(cmd)
    print(f"{cmd} done")    

    return mask

# set directory to the location of this script file to use relative paths
os.chdir(os.path.dirname(__file__))

l_m = [   0.0,    0.9,   0.98,  0.996,  0.999]
l_h = [2.0e-3, 2.0e-4, 4.0e-5, 8.0e-6, 2.0e-6]
l_de = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

arg_list = product(l_m, l_h, l_de)
print(arg_list)

num_proc = 4
if mp.cpu_count() is not None:
    num_proc = mp.cpu_count()

with mp.Pool(processes=num_proc) as pool:
    results = pool.starmap(run_cmd, arg_list)
