#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:44:27 2020

@author: pichugin
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from copy import deepcopy

def GroupByMask(Data, Mask):
	df0 = []
	df1 = []
	for i in np.arange(len(Mask)):
		if Mask[i] == 0:
			df0.append(Data[i])
		elif Mask[i] == 1:
			df1.append(Data[i])
	Result0 = 0
	if len(df0)>0:
		Result0 = pd.concat(df0, sort=False)
	Result1 = 0
	if len(df1)>0:
		Result1 = pd.concat(df1, sort=False)
	return [Result0, Result1]



def MakeHistFigure(Data_LD, Data_SD, Data_LC, Data_SC, BinsList, FileName):
	""" make figure """
	fig, ax = plt.subplots(1,1)
	HTYPE = 'step'
	CLR_LD = '#7fc97f'#(1.0, 0.0, 0.0, 1.0)
	CLR_SD = '#beaed4'#(1.0, 0.0, 1.0, 1.0)
	CLR_LC = '#fdc086'#(0.0, 0.0, 1.0, 1.0)
	CLR_SC = '#386cb0'#(0.0, 0.0, 0.0, 1.0)
	plt.hist(Data_LD['Deviation'], BinsList, histtype = HTYPE, color = CLR_LD, linewidth = 2)
	plt.hist(Data_SD['Deviation'], BinsList, histtype = HTYPE, color = CLR_SD, linewidth = 2)
	plt.hist(Data_LC['Deviation'], BinsList, histtype = HTYPE, color = CLR_LC, linewidth = 2)
	plt.hist(Data_SC['Deviation'], BinsList, histtype = HTYPE, color = CLR_SC, linewidth = 2)
	plt.legend(['Linear death', 'Step death', 'Linear cleavage', 'Step cleavage'], loc = 'upper left')
	
	LabelsFontSize = 15
	plt.xlabel('Deviation value', fontsize = LabelsFontSize)
	plt.ylabel('Observations', fontsize = LabelsFontSize)
	if flag_save_a_figure == 1: 
		plt.savefig(FileName, dpi=300)
	plt.show()
	
	return 0

""" Flags """
flag_save_a_figure = 1
ReloadData = 1


""" Read data """
if ReloadData == 1:
	Folder_name = '../02 Simulation results/'
	
	FileNamesList = []
	
	FileNamesList.append('Cleavage_constant_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_proportional_model_last_entry_nofilter.txt')
	
	FileNamesList.append('Cleavage_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_step_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_threshold_model_last_entry_nofilter.txt')
	
	FileNamesList.append('Cleavage_bottom_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_general_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_michaelis_menten_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_quadratic_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_saturating_exponent_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_sigmoid_model_last_entry_nofilter.txt')
	
	
	""" Masks to identify certain classes of models """
	MaskDeath = np.asarray( 12*[0] )# + 6*[1])
	MaskCleavage = np.asarray( 12*[1] )# + 6*[0])
	Mask1Par = np.asarray( [0, 1, 0, 1, 0, 0])
	Mask2Par = np.asarray( [1, 0, 1, 0, 1, 1])
	Mask1ParC = np.asarray( [0, 1, 0, 1, 0, 0] )# + 6*[0])
	Mask2ParC = np.asarray( [1, 0, 1, 0, 1, 1] )# + 6*[0])
	MaskAll = np.asarray( 6*[1] )
	
	""" Labels to name the models """
	Labels = [
			'Constant',
			'Proportional',
			'Top capped',
			'Breaking point',
			'Step',
			'Fracture',
			'Bottom capped',
			'Linear',
			'Michaelis Menten',
			'Quadratic',
			'Saturating exponent',
			'Sigmoid'
			]
	
	DataFromModels = len(FileNamesList)*[0]
	for i in np.arange(len(DataFromModels)):
		DataFromModels[i] = pd.read_csv(Folder_name + FileNamesList[i], sep = ',') 


""" 
NoGrowthDeviation is what would be the the deviation value in the limit case of stationary population
"""
NoGrowthDeviation = 16174.949307178855

DataD, DataC = GroupByMask(DataFromModels, MaskCleavage)	
DevC = DataC['Deviation']/NoGrowthDeviation


#
""" Deviation cdfs for _all_ cleavage models """
MinErr = min(DevC)
MaxErr = 0.35
BinsList = np.linspace(MinErr, MaxErr, num = 30)
FileName = '06_cleavage_models_all_cdf.png'
""" make figure """
fig, ax = plt.subplots(1,1)
HTYPE = 'step'
CLR = ['#fdbf6f', '#b2df8a', '#1f78b4', '#a6cee3', '#ff7f00', '#b15928', '#e31a1c', '#fb9a99', '#6a3d9a', '#ffff99', '#cab2d6', '#33a02c']



counter = 0
LegendText = []
for i in np.arange(len(DataFromModels)):
	if MaskCleavage[i] == 1:
		Dev = DataFromModels[i]['Deviation']/NoGrowthDeviation
		X_set = [0]+list(np.sort(Dev))
		Y_set = [0]+list(np.linspace(0, 1, len(Dev), endpoint=False))
		plt.plot(X_set, Y_set, color = CLR[counter], linewidth = 3)
		counter += 1
		LegendText.append(Labels[i])
		print('Minimal dev in ', Labels[i], ' is ', min(Dev))


plt.legend(LegendText, loc = 'upper left', bbox_to_anchor=(1.05, 1))

counter = 0
LineHeight = [1.0, 1.0, 0.8, 1.0, 1.0, 1.0, 0.6, 0.5, 1.0, 1.0, 0.9, 0.9]
for i in np.arange(len(DataFromModels)):
	if MaskCleavage[i] == 1:
		Dev = DataFromModels[i]['Deviation']/NoGrowthDeviation
		plt.plot([min(Dev), min(Dev)], [0, LineHeight[counter]], '--', color = CLR[counter], linewidth = 3)
		counter += 1

plt.xlim([0.14, 0.36])


ax.set_aspect(1.0/ax.get_data_ratio())
LabelsFontSize = 15
plt.xlabel('Absolute regression error', fontsize = LabelsFontSize)
plt.ylabel('Cumulative \n fraction of outcomes', fontsize = LabelsFontSize)
if flag_save_a_figure == 1: 
	plt.savefig(FileName, dpi=300, bbox_inches='tight')
plt.show()








""" Deviation cdfs for _selected_ cleavage models """

MaskSelected = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1]


MinErr = min(DevC)
MaxErr = 0.35
BinsList = np.linspace(MinErr, MaxErr, num = 30)
FileName = '06_cleavage_models_selected_cdf_v1_1.png'
""" make figure """
fig, ax = plt.subplots(1,1)
HTYPE = 'step'
CLR = ['#fdbf6f', '#b2df8a', '#1f78b4', '#a6cee3', '#ff7f00', '#b15928', '#e31a1c', '#fb9a99', '#6a3d9a', '#ffff99', '#cab2d6', '#33a02c']



counter = 0
LegendText = []
for i in np.arange(len(DataFromModels)):
	if MaskCleavage[i] == 1:
		Dev = DataFromModels[i]['Deviation']/NoGrowthDeviation
		X_set = [0]+list(np.sort(Dev))
		Y_set = [0]+list(np.linspace(0, 1, len(Dev), endpoint=False))
		if MaskSelected[i]==1:
			plt.plot(X_set, Y_set, color = CLR[counter], linewidth = 3)
			LegendText.append(Labels[i])
		counter += 1



counter = 0
LineHeight = [0.8, 1.0, 0.8, 1.0, 0.9, 1.0, 0.6, 0.5, 0.45, 1.0, 0.9, 0.6]
for i in np.arange(len(DataFromModels)):
	if MaskCleavage[i] == 1:
		Dev = DataFromModels[i]['Deviation']/NoGrowthDeviation
		if MaskSelected[i]==1:
			plt.plot([min(Dev), min(Dev)], [0, LineHeight[counter]], '--', color = CLR[counter], linewidth = 3)
		counter += 1
plt.xlim([0.14, 0.36])

ax.set_aspect(1.0/ax.get_data_ratio())
LabelsFontSize = 15
plt.xlabel('Absolute regression error', fontsize = LabelsFontSize)
plt.ylabel('Cumulative \n fraction of outcomes', fontsize = LabelsFontSize)
if flag_save_a_figure == 1: 
	plt.savefig(FileName, dpi=300, bbox_inches='tight')
plt.show()

