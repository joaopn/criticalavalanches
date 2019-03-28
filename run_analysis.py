# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   Joao PN
# @Last Modified time: 2019-03-28 15:14:02

from analysis import *
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
import glob

def parse_sim_data(datadir):

	#Lists all files
	homedir = os.getcwd()
	os.chdir(datadir)
	datafiles = [f.partition('.hdf5')[0] for f in glob.glob('*.hdf5')]
	os.chdir(homedir)

	#Removes last part ('_r%i')
	data_removed = [datafiles[i][:-3] for i in range(len(datafiles))]
	data_unique = list(set(data_removed))

	return data_unique

def run_save_plot(data_dir,filename,binsize,threshold,datatype,reps):

	#Parameters
	fig_dir = data_dir + '/plot/'

	#Plots results
	plt.figure()
	plt.title(filename)
	plt.xlabel('S')
	plt.ylabel('p(S)')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(1,1e3)
	plt.ylim(1e-5,1)

	#Runs analysis for all binsizes
	for bs in binsize:
		print(filename + ' ,running for b = {}'.format(bs))

		S_list = []
		for rep in range(reps):

			#Creates filepath
			filepath = data_dir + filename + '_r{:d}.hdf5'.format(rep)

			#Analyzes rep
			data_binned = avalanche.run_analysis(
				filepath=filepath,
				binsize=bs,
				threshold=threshold,
				datatype=datatype,
				channels=64
				)

			S_list.append(avalanche.get_S(data_binned))

		plot.pS_mean(S_list,label='bs = {}'.format(bs))

	#Creates save dir
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	#Saves fig
	str_fig = fig_dir + datatype + '_pS_' + filename + '_th{}'.format(threshold) + '.png'
	plt.savefig(str_fig)

def parametersDefault():

	#default Parameters
	binsizeDefault="1,2,4"
	thresholdDefault = 3
	repsDefault = 1
	datatypeDefault = 'coarse'
	datafolderDefault = 'dat/'

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("-t","--threshold",type=float,nargs='?',const=1,default=thresholdDefault)
	parser.add_argument("--reps",type=int,nargs='?',const=1,default=repsDefault)
	parser.add_argument("--datatype",type=str,nargs='?',const=1,default=datatypeDefault)
	parser.add_argument("--datafolder",type=str,nargs='?',const=1,default=datafolderDefault)

	#Adds string of binsizes
	parser.add_argument("-b","--binsize",type=str,nargs='?',const=1,default=binsizeDefault)
	args = parser.parse_args()
	args.binsize = [int(item) for item in args.binsize.split(',')]

	return args

if __name__ == "__main__":

	#Parse input
	args = parametersDefault()
	binsize = args.binsize
	threshold = args.threshold
	reps = args.reps
	datatype = args.datatype
	datafolder = args.datafolder

	#Selects all datasets (without the reps) in a folder to analyze
	dataset_list = parse_sim_data(datafolder)

	#Analyzes and saves plot for all datasets
	for dataset in dataset_list:
			run_save_plot(data_dir=datafolder,
				filename=dataset,
				binsize=binsize,
				threshold=threshold,
				datatype=datatype,
				reps=reps)
