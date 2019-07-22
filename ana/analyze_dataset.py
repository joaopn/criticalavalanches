# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-07-19 23:08:49
# @Last Modified by:   joaopn
# @Last Modified time: 2019-07-22 20:10:23

"""
Runs the avalanche analysis for all datasets in a folder, averaging over unique filenames ending with '_rXX.hdf5'. 
- Skips averaging for -reps 0.
- Doesn't depend on specific filename syntax.
"""

from analysis import parser, dataset
import argparse

def parametersDefault():

	#default Parameters
	binsizeDefault="1,2,4,8"
	thresholdDefault = 3
	repsDefault = 1
	datatypeDefault = 'coarse'
	datafolderDefault = 'dat/'
	bw_filterDefault = True
	fit_limDefault = 50

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("-t", "--threshold",
		type=float, nargs='?', const=1, default=thresholdDefault)
	parser.add_argument("--reps",
		type=int,   nargs='?', const=1, default=repsDefault)
	parser.add_argument("--datatype",
		type=str,   nargs='?', const=1, default=datatypeDefault)
	parser.add_argument("--datafolder",
		type=str,   nargs='?', const=1, default=datafolderDefault)
	parser.add_argument("--datamask",
		type=str,   nargs='?', const=1, default=None)
	parser.add_argument("--bw_filter",
		type=bool,  nargs='?', const=1, default=bw_filterDefault)
	parser.add_argument("--fit_lim",
		type=int,   nargs='?', const=1, default=fit_limDefault)


	#Adds string of binsizes
	parser.add_argument("-b","--binsize",type=str,nargs='?',const=1,default=binsizeDefault)
	args = parser.parse_args()
	args.binsize = [int(item) for item in args.binsize.split(',')]

	return args

if __name__ == "__main__":

	#Parse input
	args       = parametersDefault()
	binsize    = args.binsize
	threshold  = args.threshold
	reps       = args.reps
	datatype   = args.datatype
	datafolder = args.datafolder
	datamask   = args.datamask
	bw_filter  = args.bw_filter
	fit_lim = args.fit_lim

	#Gets the unique datasets in the folder
	datasets_unique = parser.sim_find_unique(datafolder)

	#Runs it for each unique dataset
	for dataset in datasets_unique:

		print('Running: '+dataset)

		#Rebuilds dataset list
		filepath = parser.sim_add_reps(datafolder+dataset,reps)

		#Runs analysis
		fig_loc = datafolder + dataset
		dataset.sim_plot_deltaT(filepath=filepath,deltaT=binsize, datatype=datatype,threshold=threshold,S_fit_max=fit_lim,bw_filter=bw_filter,save_fig=fig_loc)