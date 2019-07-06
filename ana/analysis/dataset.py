# -*- coding: utf-8 -*-
# @Author: Joao
# @Date:   2019-07-05 17:56:44
# @Last Modified by:   Joao
# @Last Modified time: 2019-07-06 13:59:53

from analysis import avalanche, fitting, plot

def plot_pS(filepath,deltaT,datatype,str_leg='Data', threshold=3,bw_filter=True,timesteps=None,channels=None):

	#Loads and thresholds data
	data_th = avalanche.analyze_sim_raw(filepath, threshold, datatype, bw_filter, timesteps, channels)

	#Bins data
	data_binned = avalanche.bin_data(data=data_th,binsize=deltaT)

	#Gets S and plots it
	S = avalanche.get_S(data_binned)
	plot.pS(S,label=str_leg)


def plot_deltaT(filepath,deltaT,datatype,timesteps=None,channels=None):

	#Parameters
	fs = 500
	bw_freqs = [0.1, 200]

	#File location
	h5_dir = 'data/'
	h5_full_dir = 'data/activity'

	file = h5py.File(filepath,'r')

	#Gets timesteps and channels from dataset
	if timesteps is None:
		timesteps = file[h5_full_dir].shape[0]
	if channels is None:
		channels = file[h5_dir+datatype].shape[0]

	#Analyses channel by channel and stacks results
	data_th = np.zeros(timesteps)
	for ch in range(channels):
		
		#Loads data
		data_ch = file[h5_dir+datatype][ch,:]

		#Analyzes sub and coarse channel data accordingly
		if datatype == 'coarse':
			if bw_filter:
				data_ch = avalanche.filter_bw_ch(data_ch,bw_freqs,fs)		
			data_th_ch = avalanche.threshold_ch(data_ch,threshold)			
		elif datatype == 'sub':
			data_th_ch = avalanche.convert_timestamps(data_ch,timesteps)

		data_th = data_th + data_th_ch

	#Thresholds data for each binsize
	for bs in deltaT:
		data_binned = avalanche.bin_data(data=data_th,binsize=bs)