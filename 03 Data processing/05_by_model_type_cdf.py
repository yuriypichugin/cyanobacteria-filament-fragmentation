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
	Folder_name_consumption = '../02 Simulation results/'
	
	FileNamesList = []
	
	FileNamesList.append('Cleavage_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_constant_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_proportional_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_step_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_threshold_model_last_entry_nofilter.txt')
	
	FileNamesList.append('Cleavage_bottom_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_general_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_michaelis_menten_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_quadratic_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_saturating_exponent_model_last_entry_nofilter.txt')
	FileNamesList.append('Cleavage_sigmoid_model_last_entry_nofilter.txt')
	
	FileNamesList.append('Death_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_constant_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_proportional_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_step_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_threshold_model_last_entry_nofilter.txt')

	FileNamesList.append('Death_bottom_capped_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_general_linear_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_michaelis_menten_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_quadratic_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_saturating_exponent_model_last_entry_nofilter.txt')
	FileNamesList.append('Death_sigmoid_model_last_entry_nofilter.txt')
	
	FileNamesList_consumption = []
	
	FileNamesList_consumption.append('Consumption_constant_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_exponent_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_inverse_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_linear_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_quadratic_concave_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_quadratic_convex_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_sigmoid_model_last_entry_nofilter.txt')
	FileNamesList_consumption.append('Consumption_step_model_last_entry_nofilter.txt')
	
	
	""" Masks to identify certain classes of models """
	MaskDeath = np.asarray( 12*[0] + 12*[1] + 8*[0])
	MaskCleavage = np.asarray( 12*[1] + 12*[0] + 8*[0])
	MaskConsumption = np.asarray( 12*[0] + 12*[0] + 8*[1])

	
	""" Labels to name the models """
	Labels = [
			'Cleavage capped',
			'Cleavage constant',
			'Cleavage linear',
			'Cleavage proportional',
			'Cleavage step',
			'Cleavage threshold',
			'Cleavage bottom capped',
			'Cleavage general linear',
			'Cleavage Michaelis Menten',
			'Cleavage quadratic',
			'Cleavage saturating exponent',
			'Cleavage sigmoid',
			'Death capped',
			'Death constant',
			'Death linear',
			'Death proportional',
			'Death step',
			'Death threshold',
			'Death bottom capped',
			'Death general linear',
			'Death Michaelis Menten',
			'Death quadratic',
			'Death saturating exponent',
			'Death sigmoid',
			'Consumption constant',
			'Consumption exponent',
			'Consumption inverse',
			'Consumption linear',
			'Consumption q. concave',
			'Consumption q. convex',
			'Consumption sigmoid',
			'Consumption step'
			]
	
	DataFromModels = (len(FileNamesList)+len(FileNamesList_consumption))*[0]
	for i in np.arange(len(FileNamesList)):
		DataFromModels[i] = pd.read_csv(Folder_name + FileNamesList[i], sep = ',')
	for i in np.arange(len(FileNamesList_consumption)):
		DataFromModels[i+len(FileNamesList)] = pd.read_csv(Folder_name_consumption + FileNamesList_consumption[i], sep = ',') 


""" 
What would be the the deviation value in the limit case of stationary population
see Focused Models/Control/Control_no_growth 
"""
NoGrowthDeviation = 16174.949307178855

zzz, DataC = GroupByMask(DataFromModels, MaskCleavage)	
zzz, DataD = GroupByMask(DataFromModels, MaskDeath)	
zzz, DataCNS = GroupByMask(DataFromModels, MaskConsumption)	
DevC = DataC['Deviation']/NoGrowthDeviation
DevD = DataD['Deviation']/NoGrowthDeviation
DevCNS = DataCNS['Deviation']/NoGrowthDeviation





""" Deviation CDFs for cleavage and death models """

DevC_sorted = np.sort(DevC)
DevD_sorted = np.sort(DevD)
DevCNS_sorted = np.sort(DevCNS)


MinErr = min(min(DevC), min(DevD), min(DevCNS))
MaxErr = 1.0
print(min(DevC))
print(min(DevD))
print(min(DevCNS))
BinsList = np.linspace(MinErr, MaxErr, num = 300)
FileName = '05_by_model_type_cdf_v1_1.png'
""" make figure """
fig, ax = plt.subplots(1,1)
HTYPE = 'step'
CLR_C = '#fc8d62'
CLR_D = '#8da0cb'
CLR_CNS = '#66c2a5'


X_C = [0]+list(np.sort(DevC))
Y_C = [0]+list(np.linspace(0, 1, len(DevC), endpoint=False))

X_D = [0]+list(np.sort(DevD))
Y_D = [0]+list(np.linspace(0, 1, len(DevD), endpoint=False))

X_CNS = [0]+list(np.sort(DevCNS))
Y_CNS = [0]+list(np.linspace(0, 1, len(DevCNS), endpoint=False))

plt.plot(X_C, Y_C, color = CLR_C, linewidth = 3)
plt.plot(X_D, Y_D, color = CLR_D, linewidth = 3)
plt.plot(X_CNS, Y_CNS, color = CLR_CNS, linewidth = 3)

plt.plot([min(DevC), min(DevC)], [0, 0.8], '--',color = CLR_C, linewidth = 2)
plt.plot([min(DevD), min(DevD)], [0, 0.8], '--',color = CLR_D, linewidth = 2)
plt.plot([min(DevCNS), min(DevCNS)], [0, 0.8], '--',color = CLR_CNS, linewidth = 2)

plt.xlim([0,1.05])

ax.set_aspect(1.0/ax.get_data_ratio())
LabelsFontSize = 15
plt.xlabel('Absolute regression error', fontsize = LabelsFontSize)
plt.ylabel('Cumulative \n fraction of outcomes', fontsize = LabelsFontSize)
if flag_save_a_figure == 1: 
	plt.savefig(FileName, dpi=300, bbox_inches='tight')
plt.show()
