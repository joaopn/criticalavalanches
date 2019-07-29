# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-07-19 23:08:49
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-29 12:41:23

"""
Plots multiple dynamic states

Runs the avalanche analysis for all datasets in a folder, averaging over unique filenames ending with '_rXX.hdf5'. 
- Skips averaging for -reps 0.
- Doesn't depend on specific filename syntax.
"""

from analysis import parser, dataset, plot
import argparse, os, pickle
import matplotlib.pyplot as plt
import numpy as np

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
	statesDefault = "subcritical,reverberating,critical"
	legDefault = None
	binsizeDefault = 2
	gaDefault = 1
	deDefault = 2
	thresholdDefault = 3
	repsDefault = 1
	datafolderDefault = 'dat/'
	bw_filterDefault = True

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("-b", "--binsize",
		type=int, nargs='?', const=1, default=binsizeDefault)
	parser.add_argument("--ga",
		type=float, nargs='?', const=1, default=gaDefault)
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
	parser.add_argument("--leg",
		type=str,   nargs='?', const=1, default=legDefault)

	#Adds string of binsizes
	parser.add_argument("--states",type=str,nargs='?',const=1,default=statesDefault)
	args = parser.parse_args()
	args.states = [item for item in args.states.split(',')]

	return args

if __name__ == "__main__":

	#Parse input
	args       = parametersDefault()
	states 	= args.states
	binsize    = args.binsize
	de = args.de
	threshold  = args.threshold
	reps       = args.reps
	datafolder = args.datafolder
	bw_filter  = args.bw_filter
	ga = args.ga
	leg = args.leg

	#Parameters
	fig_pS_size = [6.5,6.5] #in cm	

	#Figure save path
	figure_folder = 'fig/' + datafolder
	if not os.path.exists(figure_folder):
		os.makedirs(figure_folder)

	#Creates plot figures
	fig_coarse = plt.figure("coarse",figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))	
	fig_sub = plt.figure("sub",figsize=(fig_pS_size[0]/2.54,fig_pS_size[1]/2.54))		

	#Hard-coded state values (for practicity)
	states_params = {"subcritical": {"h": 2e-4, "m": 0.900, "color":'blue', "name":"Subcritical"}, "critical": {"h": 2e-06, "m": 0.999, "color":'red', "name":"Critical"},
	"poisson": {"h": 1e-3, "m": 0.000, "color":'gray',"name":"Poisson"}, "reverberating": {"h": 4e-5, "m": 0.980, "color":'green', "name":"Reveberating"}}


	for state in states:
		#Checks if the state is one of the given ones
		if state not in states_params.keys():
			raise ValueError('Invalid --state. Valid options: critical, reverberating, subcritical, poisson')

		#Builds filelist path
		filepath_norep = parser.sim_build_filename(states_params[state]['m'], states_params[state]['h'],de,ga, prefix=datafolder)[0]
		filepath_reps = parser.sim_add_reps(filepath_norep,reps)

		print('Analyzing ' + state + ', de = {:02d}, b = {:02d}, {:02d} reps: '.format(de,binsize,reps) + filepath_norep)

		#Selects color and legend
		plt_color = color_picker(states_params[state]["color"])
		str_plot = states_params[state]['name']

		#Plots coarse-sampled data
		plt.figure('coarse')
		dataset.sim_plot_pS(filepath_reps,binsize,'coarse',str_plot,threshold,bw_filter,save_fig=None, color=plt_color, show_error=False)

		#Plots sub-sampled data
		plt.figure('sub')
		dataset.sim_plot_pS(filepath_reps,binsize,'sub',str_plot,threshold,bw_filter,save_fig=None, color=plt_color, show_error=False)		

	#Touches up figures
	plt.figure('coarse')
	#plt.title('coarse-sampled')
	ax = plt.gca()
	ax.get_legend().set_visible(False)
	plt.figure('sub')
	ax = plt.gca()
	ax.get_legend().set_visible(False)
	#plt.title('sub-sampled')

	#Adds legend if required
	if leg=='coarse':
		plt.figure('coarse')
		plt.legend(loc=1,fontsize=8)
	elif leg=='sub':
		plt.figure('sub')
		plt.legend(loc=1,fontsize=8)

	#Saves figures
	figure_name = 'comparison_de{:02d}_ga{:0.2f}_b{:02d}'.format(de,ga,binsize)
	filepath_fig = figure_folder + figure_name

	fig_coarse.savefig(filepath_fig+'_coarse.pdf',bbox_inches="tight")
	fig_coarse.savefig(filepath_fig+'_coarse.png',bbox_inches="tight")
	with open(filepath_fig + '_coarse.pkl','wb') as file:
		pickle.dump(fig_coarse, file)	

	fig_sub.savefig(filepath_fig+'_sub.pdf',bbox_inches="tight")
	fig_sub.savefig(filepath_fig+'_sub.png',bbox_inches="tight")
	with open(filepath_fig + '_sub.pkl','wb') as file:
		pickle.dump(fig_sub, file)	

