#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 14:31:09 2020

@author: pichugin
"""


import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from copy import deepcopy
import os

def ExtractData(FolderName, DataFileMask, File2WriteName):
	MaxNum = 250

#	FileList = []
	DataPath = FolderName+'Results/' + DataFileMask
#	ErrorPath = FolderName+'errors/error.LC_'

	Index = 0
	for i in np.arange(MaxNum):
#		ErrorFile = ErrorPath+str(i)
		DataFile = DataPath+str(i)+'.txt'
		if (os.stat(DataFile).st_size > 0):
			
			FileData = pd.read_csv(DataFile, delimiter = ', ', engine='python')
			if Index == 0:
				LastEntryData = FileData.iloc[-1]
			else:
				LastEntryData = pd.concat([LastEntryData, FileData.iloc[-1]], ignore_index=True, axis = 1)
			Index+=1
	""" transpose to format the dataframe in the right way """
	LastEntryData = LastEntryData.T
	""" get rid of nan-s """
	LastEntryData = LastEntryData[ ~( LastEntryData['Deviation'].isna() ) ]
	""" drop empty columns """
#	LastEntryData = LastEntryData.drop(columns = ['empty_1', 'empty_2', 'empty_3', 'empty_4'])
	Column_names = LastEntryData.columns
	Empty_column_indices = [i for i, item in enumerate(Column_names) if "empty" in item]
	Empty_column_names = []
	for i in Empty_column_indices:
		Empty_column_names.append(Column_names[i])
	LastEntryData = LastEntryData.drop(columns = Empty_column_names)
	
	LastEntryData.to_csv(File2WriteName, sep = ',')


def ExtractInitialData(FolderName, DataFileMask, File2WriteName):
	MaxNum = 250

#	FileList = []
	DataPath = FolderName+'Results/' + DataFileMask
#	ErrorPath = FolderName+'errors/error.LC_'

	Index = 0
	for i in np.arange(MaxNum):
#		ErrorFile = ErrorPath+str(i)
		DataFile = DataPath+str(i)+'.txt'
		if (os.stat(DataFile).st_size > 0):
			
			FileData = pd.read_csv(DataFile, delimiter = ', ', engine='python')
			if Index == 0:
				LastEntryData = FileData.iloc[0]
			else:
				LastEntryData = pd.concat([LastEntryData, FileData.iloc[0]], ignore_index=True, axis = 1)
			Index+=1
	""" transpose to format the dataframe in the right way """
	LastEntryData = LastEntryData.T
	""" get rid of nan-s """
	LastEntryData = LastEntryData[ ~( LastEntryData['Deviation'].isna() ) ]
	""" drop empty columns """
#	LastEntryData = LastEntryData.drop(columns = ['empty_1', 'empty_2', 'empty_3', 'empty_4'])
	Column_names = LastEntryData.columns
	Empty_column_indices = [i for i, item in enumerate(Column_names) if "empty" in item]
	Empty_column_names = []
	for i in Empty_column_indices:
		Empty_column_names.append(Column_names[i])
	LastEntryData = LastEntryData.drop(columns = Empty_column_names)
	
	
	LastEntryData.to_csv(File2WriteName, sep = ',')




""" --- main --- """

FolderName = 'Cleavage capped model/'
File2WriteName = 'Cleavage_capped_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Capped_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_capped_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)


FolderName = 'Cleavage constant model/'
File2WriteName = 'Cleavage_constant_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_constant_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage linear model/'
File2WriteName = 'Cleavage_linear_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Linear_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_linear_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage proportional model/'
File2WriteName = 'Cleavage_proportional_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Proportional_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_proportional_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage step model/'
File2WriteName = 'Cleavage_step_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Step_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_step_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage threshold model/'
File2WriteName = 'Cleavage_threshold_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_threshold_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage general linear model/'
File2WriteName = 'Cleavage_general_linear_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_General_Linear_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_general_linear_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage bottom capped model/'
File2WriteName = 'Cleavage_bottom_capped_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Bottom_Capped_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_bottom_capped_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage sigmoid model/'
File2WriteName = 'Cleavage_sigmoid_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Sigmoid_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_sigmoid_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage Michaels Menten model/'
File2WriteName = 'Cleavage_michaelis_menten_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Michaelis_Menten_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_michaelis_menten_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage quadratic model/'
File2WriteName = 'Cleavage_quadratic_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Quadratic_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_quadratic_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Cleavage saturating exponent model/'
File2WriteName = 'Cleavage_saturating_exponent_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Saturating_Exponent_Cleavage_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Cleavage_saturating_exponent_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)










FolderName = 'Death capped model/'
File2WriteName = 'Death_capped_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Capped_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Capped_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_capped_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death constant model/'
File2WriteName = 'Death_constant_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Constant_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_constant_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death linear model/'
File2WriteName = 'Death_linear_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Linear_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Linear_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_linear_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death proportional model/'
File2WriteName = 'Death_proportional_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Proportional_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_proportional_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death step model/'
File2WriteName = 'Death_step_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Step_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Step_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_step_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death threshold model/'
File2WriteName = 'Death_threshold_model_last_entry_nofilter.txt'
#DataFileMask = 'CyanoFrag_Fitting_Threshold_Cleavage_model_'
DataFileMask = 'CyanoFrag_Fitting_Threshold_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_threshold_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death general linear model/'
File2WriteName = 'Death_general_linear_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_General_Linear_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_general_linear_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death bottom capped model/'
File2WriteName = 'Death_bottom_capped_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Bottom_Capped_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_bottom_capped_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death sigmoid model/'
File2WriteName = 'Death_sigmoid_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Sigmoid_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_sigmoid_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death Michaels Menten model/'
File2WriteName = 'Death_michaelis_menten_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Michaelis_Menten_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_michaelis_menten_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death quadratic model/'
File2WriteName = 'Death_quadratic_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Quadratic_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_quadratic_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)

FolderName = 'Death saturating exponent model/'
File2WriteName = 'Death_saturating_exponent_model_last_entry_nofilter.txt'
DataFileMask = 'CyanoFrag_Fitting_Sigmoid_Death_model_'
ExtractData(FolderName, DataFileMask, File2WriteName)
File2WriteFirstName = 'Death_saturating_exponent_model_first_entry.txt'
ExtractInitialData(FolderName, DataFileMask, File2WriteFirstName)