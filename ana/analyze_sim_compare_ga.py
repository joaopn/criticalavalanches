# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-07-19 23:08:49
# @Last Modified by:   joaopn
# @Last Modified time: 2019-07-28 18:16:59

"""
Plots multiple --ga parameters.

Runs the avalanche analysis for all datasets in a folder, averaging over unique filenames ending with '_rXX.hdf5'. 
- Skips averaging for -reps 0.
- Doesn't depend on specific filename syntax.
"""

from analysis import parser, dataset, plot
import argparse, os, pickle
import matplotlib.pyplot as plt

def color_picker(color_str):

	if color_str in ['red', 'r']:
		color_rgb = [255,51,51]
	elif color_str in ['blue', 'b']:
		color_rgb = [51,102,153]
	elif color_str in ['green', 'g']:
		color_rgb = [0,153,51]
	elif color_str in ['gray', 'grey']:
		color_rgb = [128,128,128]

	return np.array(color_rgb)/255	

def parametersDefault():

	#default Parameters
	stateDefault = None
	binsizeDefault = 2
	gaDefault = "1.0,1.5,2.0"
	deDefault = 2
	thresholdDefault = 3
	repsDefault = 1
	datafolderDefault = 'dat/'
	bw_filterDefault = True

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("--state",
		type=str, nargs='?', const=1, default=stateDefault)
	parser.add_argument("-b", "--binsize",
		type=int, nargs='?', const=1, default=binsizeDefault)
	parser.add_argument("--de",
		type=int, nargs='?', const=1, default=deDefault)
	parser.add_argument("-t", "--threshold",
		type=float, nargs='?', const=1, default=thresholdDefault)
	parser.add_argument("--reps",
		type=int,   nargs='?', const=1, default=repsDefault)
	parser.add_argument("--datafolder",
		type=str,   nargs='?', const=1, default=datafolderDefault)
	parser.add_argument("--datamask",
		type=str,   nargs='?', const=1, default=None)
	parser.add_argument("--bw_filter",
		type=bool,  nargs='?', const=1, default=bw_filterDefault)

	#Adds string of binsizes
	parser.add_argument("--ga",type=str,nargs='?',const=1,default=gaDefault)
	args = parser.parse_args()
	args.ga = [float(item) for item in args.ga.split(',')]

	return args

if __name__ == "__main__":

	#Parse input
	args       = parametersDefault()
	state 	= args.state
	binsize    = args.binsize
	de = args.de
	threshold  = args.threshold
	reps       = args.reps
	datafolder = args.datafolder
	bw_filter  = args.bw_filter
	ga = args.ga

	#Parameters
	fig_pS_size = [6.5,6.5] #in cm
	color_sub = [1,0.647,0] #orange
	gradient_brightness = 0.6 #how bright the lightest color is
	

	#Figure save path
	figure_folder = 'fig/' + datafolder
	if not os.path.exists(figure_folder):
		os.makedirs(figure_folder)

	#Hard-coded state values (for practicity)
	states = {"subcritical": {"h": 2e-4, "m": 0.900, "color":'blue'}, "critical": {"h": 2e-06, "m": 0.999, "color":'red'},
	"poisson": {"h": 1e-3, "m": 0.000, "color":'gray'}, "reverberating": {"h": 4e-5, "m": 0.980, "color":'green'}}

	#Checks if the state is one of the given ones
	if state not in states.keys():
		raise ValueError('Invalid --state. Valid options: critical, reverberating, subcritical, poisson')

	#Creates plot figure
	fig = plt.figure(figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))

	#Builds filelist path
	filepath_list = parser.sim_build_filename(states[state]['m'], states[state]['h'],de,ga, prefix=datafolder)
	filepath_noga = parser.sim_build_filename(states[state]['m'], states[state]['h'],de,None)
	filepath_all = []

	#Sets colors for plots
	plt_color = plot._color_gradient_rgba(states[state]["color"],len(ga),brightness_lim=gradient_brightness)
	plt_color = plot._convert_rgba_rgb(plt_color)

	#Plot zorder
	zorder = list(range(len(ga)))
	zorder.reverse()
	zorder = [zi+2 for zi in zorder] #reserves zorder=1 for sub

	#Plots the coarse-sampled data for each ga, and all reps
	for i in range(len(ga)):
		print('[coarse]: ' + filepath_list[i])

		#Builds rep filepaths
		filepath_reps = parser.sim_add_reps(filepath_list[i],reps)
		filepath_all.append(filepath_reps)
		
		#Plots coarse-sampled data
		str_plot = r'$\gamma$ = {:0.1f}'.format(ga[i])
		dataset.sim_plot_pS(filepath_reps,binsize,'coarse',str_plot,threshold,bw_filter,save_fig=None, color=plt_color[i,:], show_error=False, zorder=zorder[i])

	#Plots the sub-sampled data, averaged over all ga and rep
	filepath_all_flat = [item for sublist in filepath_all for item in sublist]
	print('[sub]: ' + filepath_noga[0])
	dataset.sim_plot_pS(filepath_all_flat,binsize,'sub','spikes',threshold,bw_filter,save_fig=None, color=color_sub, lineType='-', show_error=False,zorder=-1)

	#Touches up figure and saves it
	plt.legend(loc=3,fontsize=8)
	filepath_figure = figure_folder + filepath_noga[0] + '_b{:02d}'.format(binsize)
	fig.savefig(filepath_figure+'.pdf',bbox_inches="tight")
	fig.savefig(filepath_figure+'.png',bbox_inches="tight")
	with open(filepath_figure + '.pkl','wb') as file:
		pickle.dump(fig, file)	