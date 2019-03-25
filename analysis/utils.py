"""
Module for utilities.

"""
# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-03-22 12:58:11
# @Last Modified by:   joaopn
# @Last Modified time: 2019-03-25 10:08:17

import h5py
import numpy as np

def load_hdf5(filepath,datafolder, ch_id):
	"""Loads HDF5 data.
	
	Arguments:
		filepath {str} -- path to file

		datafolder {str} -- folder structure to dataset. Ex: 'data/data01'

		ch_id {int} -- if of channels to load.
	
	Returns:
		data -- data array
	"""

	file = h5py.File(filepath,'r')
	data = file[datafolder][ch_id,:]

	return data

def bin_data(data, binsize):
	"""Bins 1D data
	
	Arguments:
		data {float} -- numpy 1D array with timesteps.
		binsize {float} -- size of each bin, in units of timesteps.
	
	Returns:
		data_binned -- data binned in bins of size binsize.
	"""


	bin_array = np.arange(0,data.max(),binsize)
	data_hist = np.histogram(data,bin=bin_array)
	data_binned = np.asmatrix(data_hist[0])

	return data_binned

