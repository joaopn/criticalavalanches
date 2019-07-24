import os
import h5py
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

import networkx as nx

from circles import *

file = h5py.File("/Users/paul/Desktop/random.hdf5",'r')
try:
    n_x   = file['/neurons/pos_x'][:]
    n_y   = file['/neurons/pos_y'][:]
    n_R_s = 7.5
    n_R_d = file['/neurons/radius_dendritic_tree'][:]
    # n_l = file['/neurons/axon_length'][:]
except: pass

try:
    axon_segments_x = file['/axons/segments_x'][:]
    axon_segments_y = file['/axons/segments_y'][:]
    # no nans in hdf5, 0 is the default padding
    axon_segments_x = np.where(axon_segments_x==0, np.nan, axon_segments_x)
    axon_segments_y = np.where(axon_segments_y==0, np.nan, axon_segments_y)
except: pass

try:
    e_x   = file['/electrodes/pos_x'][:]
    e_y   = file['/electrodes/pos_y'][:]
    e_d_zone = file['/meta/elec_dead_zone'][:]
except: pass

try:
    mtx   = file['/connectivity_matrix'][:]
except: pass
file.close()

try:
    fig, ax = plt.subplots(figsize=[6.4, 6.4])
    # for i in range(1000):
    for i in range(len(axon_segments_x)):
        ax.plot(axon_segments_x[i], axon_segments_y[i], color='black', lw=0.25, zorder=0, alpha=0.3)
except:pass


# pos = {0: (40, 20), 1: (20, 30), 2: (40, 30), 3: (30, 10)}
# pos = {}
# G=nx.Graph()
# G.add_nodes_from(np.arange(len(n_x)))
# for i in range(len(n_x)):
#     pos[i] = (n_x[i], n_y[i])
    # G.node[i]['pos'] = (n_x[i], n_y[i])
# nx.draw(G, pos=pos, ax=ax, node_size=15, node_color='black')

try:
    circles(n_x, n_y, n_R_s, ax=ax, fc='white', ec='none', alpha=1,
        lw=0.5, zorder=4)
    circles(n_x, n_y, n_R_s, ax=ax, fc='none',  ec='black', alpha=0.3,
        lw=0.5, zorder=5)
    # circles(n_x, n_y, n_R_d, ax=ax, fc='gray',  ec='none', alpha=0.001,
    #     zorder=-5)
except: pass

# plot electrodes
try:
    circles(e_x, e_y, e_d_zone, ax=ax, fc='orange',  ec='none', alpha=0.9, zorder=+5)
except:
    pass



# plot one guy especiall
color = -1
def highlight_single(n=0):
    try:
        global color
        color = (color+1)%10
        ax.plot(axon_segments_x[n], axon_segments_y[n],
            color=f'C{color}', lw=0.7, zorder=10)
        circles(n_x[n], n_y[n], n_R_s,    ax=ax,
            ls='-',  lw=0.7, alpha=1.0, zorder=15, fc='white', ec=f'C{color}')
        circles(n_x[n], n_y[n], n_R_d[n], ax=ax,
            ls='--', lw=0.7, alpha=0.9, zorder=7,  fc='none', ec=f'C{color}')
        circles(n_x[n], n_y[n], n_R_d[n], ax=ax,
            ls='--', lw=0.0, alpha=0.1, zorder=6,  fc=f'C{color}',   ec='none')
    except: pass

def highlight_some():
    try:
        highlight_single(1)
        highlight_single(3)
        highlight_single(12)
        highlight_single(13)
        highlight_single(123)
        highlight_single(138)
    except: pass

def highlight_connected(n=0):
    try:
        connected = np.where(mtx[n]==1)
        rejected  = np.where(mtx[n]==2)
        color=0

        ax.plot(axon_segments_x[n], axon_segments_y[n],
            color=f'C{color}', lw=0.7, zorder=10)
        circles(n_x[n], n_y[n], n_R_s,    ax=ax,
            ls='-',  lw=0.7, alpha=1.0, zorder=15, fc='white', ec=f'C{color}')

        for i in connected:
            circles(n_x[i], n_y[i], n_R_s,    ax=ax,
                ls='-',  lw=0.4, alpha=1.0, zorder=15, fc='white', ec=f'C{color}')
            # circles(n_x[i], n_y[i], n_R_d[i], ax=ax,
            #     ls='--', lw=0.7, alpha=0.9, zorder=7,  fc='none', ec=f'C{color}')
            circles(n_x[i], n_y[i], n_R_d[i], ax=ax,
                ls='--', lw=0.0, alpha=0.1, zorder=6,  fc=f'C{color}',   ec='none')

        color=1
        for i in rejected:
            circles(n_x[i], n_y[i], n_R_s,    ax=ax,
                ls='-',  lw=0.4, alpha=1.0, zorder=15, fc='white', ec=f'C{color}')
            # circles(n_x[i], n_y[i], n_R_d[i], ax=ax,
            #     ls='--', lw=0.7, alpha=0.9, zorder=7,  fc='none', ec=f'C{color}')
            circles(n_x[i], n_y[i], n_R_d[i], ax=ax,
                ls='--', lw=0.0, alpha=0.01, zorder=6,  fc=f'C{color}',   ec='none')
    except: pass


# highlight_connected(99)
highlight_some()




