"""
@author: Jérémie Gaffarel & Rudi Smits - TU Delft - Faculty of Aerospace Engineering - BSc 2 (2018-2019)

This module is meant to open MAV flight measurements
"""

import numpy as np

def get_data(flight=1):
	""" Get Data:
	Input (optional):
		* flight (int between 1 and 4): fligh id, default is flight 1 (with 2 MAVs)
	Outputs:
		* gt/onb: list containing arrays, containing the desired data extracted from the files
			Example:	gt[0] =		["x", "y", "vx", "vy", "psi", "z"] of MAV 1 (gt[1] for MAV 2)
						onb[0] =	["id", "id_other", "rssi", "x_rel", "y_rel", "vx", "vy", "vx_other", "vy_other", "psi", "psi_other", "z_rel"] of MAV 1 (onb[1] for MAV 2)
		* onb_avg: list containing the average of the ["rssi", "x_rel", "y_rel", z_rel] measurements of the two MAVs
		* gt_rel: ground true relative ["x", "y", "z"] positions of the two MAVs
	"""

	flight = 1 if flight < 0 or flight > 4 else flight	# If the specified flight number is not in (1, 4), return flight 1 data
	txt = []
	gt = []
	onb = []
	# Filenames are encoded in a list, makes it easier to edit; {} will be replaced by the flight number
	filenames = ["n2_f{}_data_gt_mav1", "n2_f{}_data_gt_mav2", "n2_f{}_data_onboard_mav1", "n2_f{}_data_onboard_mav2"]
	# Specify if flights are on-board (False) or ground-truth (True) measurements
	types = [True, True, False, False]
	# Specify what each column in the onb (on-board) and gt(ground-truth) files is
	onb_col = ["time", "id", "id_other", "rssi", "x_rel", "y_rel", "vx", "vy", "vx_other", "vy_other", "psi", "psi_other", "z_rel"]
	gt_col = ["time", "x", "y", "vx", "vy", "psi", "z"]

	# Get the content of each data/measurment file
	for filename in filenames:
		with open("data/" + filename.format(flight) + ".txt", "r", encoding="utf8") as code:
			txt.append(code.read())

	i = 0
	# Extract all of the measurements from the text files to put them in lists
	for type in types:
		if type:	# If data is from ground-truth
			gt.append(extract_from_txt(txt[i], gt_col))		# gt.append [time, vx, vy, psi, z] for instance
		else:		# Data is from on-board
			onb.append(extract_from_txt(txt[i], onb_col))	# onb.append [time, vx, vy, rssi] for instance
		i += 1

	# Define the commom time = to the ground truth time of MAV 1
	time = gt[0]["time"]
	# Linearly interpolare all of the flights to make them use the common time
	for i in range(len(gt)):
		gt[i] = same_size(time, gt[i])
	for i in range(len(onb)):
		onb[i] = same_size(time, onb[i])
	
	# Compute the average of the on-board measurements
	onb_avg = avg(onb)
	# Get the relative position of the 2 MAVs in x/y/z-directions, from their absolute positions
	gt_rel = abs_to_rel(gt, onb_avg)

	return gt, onb, onb_avg, gt_rel, time

def abs_to_rel(d, measurements):
	"""	Absolute to Relative:
	Inputs:
		* d (list): data we want to compute the relative positions from
		* measurements (list): measurements of absolute x/y/z relative positions
	Outputs:
		* gt_rel (dictionary): x/y/z relative positions
	"""
	d_rel = []
	coord = ["x", "y", "z"]
	for i_d in coord:							# for each coordinate...
		mes = measurements[i_d + "_rel"]		# get the measurement list for that coordinate
		d_rel_t01 = d[0][i_d] - d[1][i_d]		# try to get the relative measurements by doing mes_1 - mes_2...
		d_rel_t10 = d[1][i_d] - d[0][i_d]		# or mes_2 - mes_1
		# Compute error for 1 - 2 and 2 - 1
		e01 = sum(np.fabs(mes - d_rel_t01))	
		e10 = sum(np.fabs(mes - d_rel_t10))
		# Define the relative measurements as the ones giving the lowest error
		if e01 < e10 and i_d != "y" or e01 > e10:	
			d_rel.append(d_rel_t10)
		else:
			d_rel.append(d_rel_t01)
	return dict(zip(coord, d_rel))

def same_size(time, d):
	""" Same size
	Inputs:
		- time [np.array]: time array, on which the measurements will me aligned
		- np [np.ndarray]: array of measurment out of line: [x, y, z]
	Output: [np.darray]: array of measurment, based on the same time as "time" input
	"""
	rtrn = dict()
	for key, x in d.items():
		if key != "time":
			# Linearly interpolate the data to sync it on the defined time list
			sync = np.interp(time, d["time"], x) 
			rtrn[key] = sync
	return rtrn

def avg(d):
	""" Average
	Input: data [dict]: dictionary of measurements
	Outputs: averaged measurements [dict]
	"""
	rtrn = dict()
	for key, x in d[0].items():
		# for the relative measurements, the relative average is done by doing (mes - mes) / 2
		# this is because the relative measurement of the other MAV is negative (coordinate system the other way around)
		if key in ["x_rel", "y_rel", "z_rel"]:
			rtrn[key] = (d[0][key] - d[1][key]) / 2
		elif key in ["rssi", "vx", "vy"]:
			rtrn[key] = (d[0][key] + d[1][key]) / 2
	return rtrn

def extract_from_txt(txt, col):
	""" Extract from text:
	Inputs:
		* txt (string): the opened measurement file
		* pos ([[int], [int]]): position of the interesting measurements, and of the velocities
	Outputs:
		* 2D array containing the extracted data
	"""
	v_pos = [col.index("vx"), col.index("vy")]								# Indexes of the velocities
	rslt = []																# PLaceholder for the return list
	last_const = False
	this_const = False
	last_d = []
	i = 0
	lines = txt.split("\n")
	for line in lines:														# For every line in the text...
		d_t = line.split("\t")												# Split the columns
		if "previous" not in locals():
			previous = np.zeros(len(d_t))
		if len(d_t) != 1:													# Ignore if the line was empty
			d = [float(i) for i in d_t]										# Convert to floats
			this_const, previous = remove_stationary(d, v_pos, previous)	# Check if the drone is in motion
			if i > 0 and not(this_const) and not(last_const):				# If the drone is currently in a moving phase...
				rslt.append(last_d)											# append the measurements to the return list
			i += 1
			last_d = d
			last_const = this_const
	rtrn_dic = dict(zip(col, np.array(rslt).T))
	return rtrn_dic

def remove_stationary(d, pos, previous):
	"""
	Remove measurements if they are equal to the previous one, that is when the MAV is static (before/after the actual test started)
	"""
	return np.array_equal(np.array(d)[pos], previous), np.array(d)[pos]
	return False, np.array(d)[pos]
