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

JobID = 0


""" initialization """

TheModel = CMS.LinearCleavageProjectionMatrix

W = 0.0
K = np.power(10, np.random.uniform(3, 4))
P_comp = 0.0
D_comp = np.power(10, np.random.uniform(-2.0, 0.5))
c_intercept = np.power(10, np.random.uniform(-6, -3))
c_slope = np.power(10, np.random.uniform(-6, -3))
Parameters = [W, K, c_intercept, c_slope, 0, 0, 0, 0, P_comp, D_comp]

Deviation = Compute_dev(Parameters, NMAX, DataFolder, TheModel)
print(Deviation)
