# -*- coding: utf-8 -*-
# @Author: Joao PN
# @Date:   2019-03-25 16:45:25
# @Last Modified by:   joaopn
# @Last Modified time: 2019-04-02 00:27:40

from analysis import avalanche, plot
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
import glob
import h5py

def parse_sim_data(datadir, datamask = None):

	#Lists all files
	homedir = os.getcwd()
	os.chdir(datadir)
	datafiles = [f.partition('.hdf5')[0] for f in glob.glob('*.hdf5')]
	os.chdir(homedir)

	# only keep files for which path contains the mask
	if datamask is not None:
		datafiles = [item for item in datafiles if datamask in item]
		# print(datafiles)

	#Removes last part ('_r00')
	data_removed = [datafiles[i][:-4] for i in range(len(datafiles))]
	data_unique = list(set(data_removed))

	return data_unique

def parse_sim_thresholded(datadir, bw_filter=False, datamask = None):

	if bw_filter:
		dir_thresholded = 'thresholded_filtered/'
	else:
		dir_thresholded = 'thresholded_unfiltered/'

	#Lists all files
	homedir = os.getcwd()
	os.chdir(datadir + dir_thresholded)
	datafiles = [f.partition('.hdf5')[0] for f in glob.glob('*.hdf5')]
	os.chdir(homedir)

    # only keep files for which path contains the mask
    if datamask is not None:
        datafiles = [item for item in datafiles if datamask in item]
        # print(datafiles)

	return datafiles

def parametersDefault():

	#default Parameters
	binsizeDefault="1,2,4"
	thresholdDefault = 3
	repsDefault = 2
	datatypeDefault = 'coarse'
	datafolderDefault = 'dat/'
	modeDefault = 'save_plot'
	bw_filterDefault = False

	#Parse input
	parser = argparse.ArgumentParser()

	#Add single parameters
	parser.add_argument("--mode",
		type=str,   nargs='?', const=1, default=modeDefault)
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


	#Adds string of binsizes
	parser.add_argument("-b","--binsize",type=str,nargs='?',const=1,default=binsizeDefault)
	args = parser.parse_args()
	args.binsize = [int(item) for item in args.binsize.split(',')]

	return args

def save_plot(data_dir,filename,threshold,datatype,reps,binsize,bw_filter):

	#Save location
	if bw_filter:
		fig_dir = data_dir + 'plot_filtered/'
	else:
		fig_dir = data_dir + 'plot_unfiltered/'

	#Plots results
	plt.figure()
	plt.title(filename)
	plt.xlabel('S')
	plt.ylabel('p(S)')
	plt.yscale('log')
	plt.xscale('log')
	plt.xlim(1,1e3)
	plt.ylim(1e-5,1)
	X = np.arange(5,64)
	Y = np.power(X,-1.5)
	Y = 5*Y/np.sum(Y)
	plt.plot(X,Y,linestyle='--',color='black',label=r'$p(S)~S^{-1.5}$')

	#Runs analysis for all binsizes
	for bs in binsize:
		print(filename + ' ,running for b = {}'.format(bs))

		S_list = []
		for rep in range(reps):

			#Creates filepath
			filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)

			#Analyzes rep
			data_thresholded = avalanche.analyze_sim_raw(
				filepath=filepath,
				threshold=threshold,
				datatype=datatype,
				bw_filter=bw_filter
				)
			data_binned = avalanche.bin_data(data=data_thresholded,binsize=bs)

			S_list.append(avalanche.get_S(data_binned))

		plot.pS_mean(S_list,label='bs = {}'.format(bs))

	#Creates save dir
	if not os.path.exists(fig_dir):
		os.makedirs(fig_dir)

	#Saves fig
	str_fig = fig_dir + datatype + '_pS_' + filename + '_th{}'.format(threshold) + '.png'
	plt.savefig(str_fig)
	plt.close()

def save_thresholded(data_dir,filename,threshold,reps,bw_filter,timesteps=None):

	#Parameters
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
	data_save = data_dir + dir_threshold
	str_savefile = data_save + filename + '_th{:0.1f}.hdf5'.format(threshold)
	datatypes = ['coarse', 'sub']

	#Gets timesteps
	if timesteps is None:
		filepath_0 = data_dir + filename + '_r00.hdf5'
		file_raw = h5py.File(filepath_0,'r')
		timesteps = file_raw['data/activity'].shape[0]
		file_raw.close()

	#Creates path to .h5 save file
	if not os.path.exists(data_save):
		os.makedirs(data_save)
	file = h5py.File(str_savefile,'w') #Overwrites old file


	#Saves thresholded data ('coarse' and 'sub')
	for datatype in datatypes:

		#Creates .h5 dataset
		file.create_dataset(datatype,shape=(reps,timesteps),dtype=int,chunks=(1,timesteps),compression=1,maxshape=(None,timesteps))
		file[datatype].attrs['threshold'] = threshold

		for rep in range(reps):

			#Creates filepath
			filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)

			#Analyzes rep
			data_thresholded = avalanche.analyze_sim_raw(
				filepath=filepath,
				threshold=threshold,
				datatype=datatype,
				bw_filter=bw_filter
				)
			file[datatype][rep,:] = np.int32(data_thresholded)

	#Copies population activity data
	file.create_dataset('activity',shape=(reps,timesteps),dtype=float,chunks=(1,timesteps),compression=1,maxshape=(None,timesteps))
	for rep in range(reps):
		filepath = data_dir + filename + '_r{:02d}.hdf5'.format(rep)
		file_raw = h5py.File(filepath,'r')
		file['activity'][rep,:] = file_raw['data/activity'][:]
		file_raw.close()

	#Flushes and closes file
	file.close()

if __name__ == "__main__":

	#Parse input
	args       = parametersDefault()
	mode       = args.mode
	binsize    = args.binsize
	threshold  = args.threshold
	reps       = args.reps
	datatype   = args.datatype
	datafolder = args.datafolder
	datamask   = args.datamask
	bw_filter  = args.bw_filter

	#Does the requested operation
	if mode == 'save_plot':
		dataset_list = parse_sim_data(datafolder, datamask)
		print('Analyzing and plotting')
		for i in range(len(dataset_list)):
			save_plot(data_dir=datafolder,
				filename=dataset_list[i],
				binsize=binsize,
				threshold=threshold,
				datatype=datatype,
				reps=reps,
				bw_filter=bw_filter)

	elif mode == 'threshold':
		dataset_list = parse_sim_data(datafolder, datamask)
		for i in range(len(dataset_list)):
			print('Thresholding and saving: ' + dataset_list[i])
			save_thresholded(data_dir=datafolder,
				filename=dataset_list[i],
				threshold=threshold,
				reps=reps,
				bw_filter=bw_filter)

	elif mode == 'save_ps':
		dataset_list = parse_sim_thresholded(datafolder, bw_filter, datamask)
		for i in range(len(dataset_list)):
			print('Analysing thresholded data: ' + dataset_list[i])
			for j in range(len(binsize)):
				avalanche.save_sim_pS(data_dir=datafolder,
					dataset=dataset_list[i],
					binsize=binsize[j],
					bw_filter=bw_filter
					)
