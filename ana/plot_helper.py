from analysis import avalanche, plot
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
import glob


def plot_dist_mean(ax, data, **kwargs):
    kwargs.setdefault("lw", 0.5)
    ax.plot(data[0], data[1], **kwargs)


def plot_dist_err(ax, data, **kwargs):
    kwargs.setdefault("alpha", 0.1)
    kwargs.setdefault("lw", 0.0)
    ax.fill_between(data[0], data[1]-data[2]/2., data[1]+data[2]/2., **kwargs)

if __name__ == "__main__":

    fig = plt.figure(figsize=(3.5, 1.6), dpi=300)
    # ax = fig.add_axes([0.1, 0.1, 1, 1])
    ax = fig.add_subplot(1,1,1)

    colorid = 0
    for b in [1, 2, 4, 8, 16]:
        foo = np.loadtxt(f"/Users/paul/owncloud/mpi/simulation/critical_avalanches/nocc_01/dat/analyzed_filtered/m0.99400_h3.000e-05_d05_th3.0_rep50/pS_coarse_b{b:02d}.tsv")

        plot_dist_mean(ax, foo, zorder=0, color=f"C{colorid}", label=f"$b = {b:02d}$")
        plot_dist_err(ax, foo, zorder=-5, color=f"C{colorid}")
        colorid += 1


    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend()
