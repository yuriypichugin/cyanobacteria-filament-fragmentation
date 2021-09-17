#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 08:56:39 2020

@author: pichugin
"""

import numpy as np
#import csv
#import matplotlib.pyplot as plt
#import pandas as pd
#import seaborn as sns
from copy import deepcopy
import Cyano_model_support_v3_1 as CMS
import Lib_ExperimetalDeviation_v2_4 as LDE
import sys
import io

#import sys

def Compute_dev(Parameters, NMAX, DataFolder, Model):
	Delta3 = LDE.Dev_fig3(Parameters, NMAX, DataFolder, Model)
#	Delta4 = LDE.Dev_fig4(Parameters, NMAX, DataFolder, Model)
	Delta5 = LDE.Dev_fig5(Parameters, NMAX, DataFolder, Model)
	Delta8 = LDE.Dev_fig8(Parameters, NMAX, DataFolder, Model)
	Delta12 = LDE.Dev_fig12(Parameters, NMAX, DataFolder, Model)
	
	Total_dev = np.sqrt(Delta3+Delta5+Delta8+Delta12)

	return Total_dev

def RecordResults(File2Write, Parameters, Total_dev):
	
	String2Write = ""
	for elem in Parameters:
		String2Write += str(elem) + ', '
	String2Write += str(Total_dev) + '\n'
	File2Write.write(String2Write)
	File2Write.flush()
	return 0

""" parameters and constants """
flag_write2file = 0
flag_DebugMode = 1

""" how many data points to test """
Repeats = 1
#Repeats = 1
""" maximal group size """
NMAX = 32


DataFolder = 'Data/'

JobID = int(sys.argv[1]) 
#JobID = 0
FileName = "Results/CyanoFrag_Fitting_Linear_Death_model_"+str(JobID)+".txt"
File2Write = open(FileName, "w")
String2Write = "W, K, Conc_treshold, Slope, empty_1, empty_2, empty_3, empty_4, P_comp, D_comp, Deviation"+'\n'
File2Write.write(String2Write)

""" initialization """

TheModel = CMS.LinearDeathModel

W = np.power(10, np.random.uniform(-3.0, 0.0))
K = np.power(10, np.random.uniform(3, 4))
P_comp = 1.0
D_comp = np.power(10, np.random.uniform(-2.0, 0.5))
Conc_treshold = np.power(10, np.random.uniform(0, 4))
Slope = np.power(10, np.random.uniform(-6, -2))
Parameters = [W, K, Conc_treshold, Slope, 0, 0, 0, 0, P_comp, D_comp]



Deviation = Compute_dev(Parameters, NMAX, DataFolder, TheModel)
	
RecordResults(File2Write, Parameters, Deviation)

ParametersBest = deepcopy(Parameters)
DeviationBest = deepcopy(Deviation)

PertScale = 2e0
while PertScale > 1e-3:
	BetterFound = 0
	
	CurrentBestParams = deepcopy(ParametersBest)
	CurrentBestDeviation = deepcopy(DeviationBest)
	for i in [0,1,2,3,9]:
		Parameters_Plus = deepcopy(ParametersBest)
#		Parameters_Plus[i] = Parameters_Plus[i] * (1 + PertScale)
		Parameters_Plus[i] = Parameters_Plus[i] * np.exp( PertScale)
		TestDeviation = Compute_dev(Parameters_Plus, NMAX, DataFolder, TheModel)
		if (1 - TestDeviation/CurrentBestDeviation) > 1e-3:
			CurrentBestParams = deepcopy(Parameters_Plus)
			CurrentBestDeviation = deepcopy(TestDeviation)
			BetterFound = 1
		
		Parameters_Minus = deepcopy(ParametersBest)
#		Parameters_Minus[i] = Parameters_Minus[i] * (1 - PertScale)
		Parameters_Minus[i] = Parameters_Minus[i] * np.exp( - PertScale)
		TestDeviation = Compute_dev(Parameters_Minus, NMAX, DataFolder, TheModel)
		if (1 - TestDeviation/CurrentBestDeviation) > 1e-3:
			CurrentBestParams = deepcopy(Parameters_Minus)
			CurrentBestDeviation = deepcopy(TestDeviation)
			BetterFound = 1
	
	if BetterFound == 1:
		ParametersBest = deepcopy(CurrentBestParams)
		DeviationBest = deepcopy(CurrentBestDeviation)
		RecordResults(File2Write, ParametersBest, DeviationBest)
		PertScale = min(2.0, PertScale * 1.1)
	else:
		PertScale = PertScale / 2.0


File2Write.close()