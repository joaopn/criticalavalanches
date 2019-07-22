# -*- coding: utf-8 -*-
# @Author: joaopn
# @Date:   2019-07-19 12:33:04
# @Last Modified by:   joaopn
# @Last Modified time: 2019-07-22 19:58:17

"""
Parses file strings
"""

import os, glob

def sim_build_filename(m,h,de=None,ga=None,version='new'):
	"""Builds a filename list using the syntax from the simulation (new version)
	
	Args:
		m (float): Branching parameter
		h (float): drive per neuron
		de (int, optional): Inter-electrode distance
		ga (float, optional): electrode FOV exponent (Signal~1/R^ga)
		version (str, optional): Syntax version [new/old]
	
	Returns:
		list: List of filenames created from the combination of parameters
	"""

	#Make sure everything is a list
	if type(m) is not list:
		m = [m]
	if type(h) is not list:
		h = [h]
	if type(de) is not list:
		de = [de]
	if type(ga) is not list:
		ga = [ga]

	#Parse input
	if len(m) != len(h):
		raise ValueError('m and h have different lengths')

	#Builds paths
	filenames = []
	for i in range(len(m)):
			for de_i in de:
				for ga_i in ga:
					str_file = 'm{:.5f}_h{:.3e}_de{:02d}_ga-{:0.2f}'.format(m[i],h[i],de_i,ga_i)
					filenames.append(str_file)

	return filenames


def sim_add_reps(filepath_base,reps):
	"""Builds a list of filepaths in the format '[filepath]_rXX.hdf5'
	
	Args:
	    filepath_base (str): Base filepath
	    reps (int): Number of files
	"""	
	if reps == 0:
		return filepath_base

	filepath = []
	for i in range(reps):
		filepath.append(filepath_base+'_r{:02d}.hdf5'.format(i))

	return filepath

def sim_find_unique(datadir, datamask = None):
	"""Parses simulation data from a folder, returning the unique names without "_rXX.hdf5"

	
	Args:
		datadir (str): Folder to search
		datamask (str, optional): Adds a mask, excluding some files
	
	Returns:
		str: List of unique dataset names in the folder.
	"""

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

def sim_find_unique_no_d(datadir):
	#Returns the unique filenames with no e.g. "_d00_th3.0.hdf5". 
	#ALL datasets must have same values for "dxx"

	#Parameters
	d_check_max = 30

	#Lists all files
	homedir = os.getcwd()
	os.chdir(datadir)
	datafiles = [f.partition('.hdf5')[0] for f in glob.glob('*.hdf5')]
	os.chdir(homedir)	

	#Checks all strings and adds present dxx to list
	d_check_list = ["d{:02d}".format(d_i) for d_i in range(1,d_check_max+1)]
	d_list = []
	for d_i in d_check_list:
		if any(d_i in data_i for data_i in datafiles):
			d_list.append(int(d_i[1:]))

	#Removes dxx_thyy from filename
	for d in d_list:
		datafiles = [f.partition("d{:02d}".format(d))[0] for f in datafiles]

	#Gets unique strings
	datafiles = list(set(datafiles))

	return datafiles,d_list

def sim_find_thresholded(datadir, bw_filter=False, datamask = None):
	"""Parses thresholded simulation data from a folder
	
	Args:
		datadir (str): Folder to search
		bw_filter (bool, optional): Whether the thresholded data is filtered
		datamask (None, optional): Adds a mask, excluding some files
	
	Returns:
		str: List of thresholded datasets in the folder
	"""
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

def sim_find_thresholded_no_d(datadir, bw_filter=False, datamask = None):
	"""Parses thresholded simulation data from a folder, without termination
	
	Args:
		datadir (TYPE): Description
		bw_filter (bool, optional): Description
		datamask (None, optional): Description
	
	Returns:
		TYPE: Description
	"""
	#Returns the unique filenames with no e.g. "_d00_th3.0.hdf5". 
	#ALL datasets must have same values for "dxx"


	d_check_max = 30

	#Gets datafiles
	datafiles = parser.sim_find_thresholded(datadir,bw_filter,datamask)

	#Checks all strings and adds present dxx to list
	d_check_list = ["d{:02d}".format(d_i) for d_i in range(1,d_check_max+1)]
	d_list = []
	for d_i in d_check_list:
		if any(d_i in data_i for data_i in datafiles):
			d_list.append(int(d_i[1:]))

	#Removes dxx_thyy from filename
	for d in d_list:
		datafiles = [f.partition("d{:02d}".format(d))[0] for f in datafiles]

	#Gets unique strings
	datafiles = list(set(datafiles))

	return datafiles,d_list