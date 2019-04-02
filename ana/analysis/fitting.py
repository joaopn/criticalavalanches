# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-04-01 01:44:18
# @Last Modified by:   joaopn
# @Last Modified time: 2019-04-02 11:12:55

import numpy as np
import h5py

def tau_linear(data,deltaT = 2):

	#Defines X and Y
	X = data[:-1] #A(t)
	Y = data[1:] #A(t+1)

	fit = np.polyfit(X,Y,1)
	tau = -deltaT/np.log(fit[0])

	return tau

def tau_sim_dataset(m,h,d,threshold,data_dir,bw_filter):
	
	#Sets up filepaths
	if bw_filter:
		dir_threshold = 'thresholded_filtered/'
	else:
		dir_threshold = 'thresholded_unfiltered/'
	filename = 'm{:0.5f}_h{:0.3e}_d{:02d}_th{:0.1f}.hdf5'.format(m,h,d,threshold)
	str_load = data_dir + dir_threshold + filename

	#Loads file
	file = h5py.File(str_load,'r')
	reps = int(file['activity'].shape[0])

	tau_rep = np.zeros(reps)

	for i in range(reps):
		tau_rep[i] = tau_linear(file['activity'][i,:])

	return tau_rep

def m_avalanche(data)

	#Finds crossings of the signal to positive
	id_cross=np.where(np.sign(data[:-1]) != np.sign(data[1:]))[0] + 1
	id_cross_plus = id_cross[data[id_cross]>0]

	