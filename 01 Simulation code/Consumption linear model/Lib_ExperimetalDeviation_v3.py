#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 08:47:06 2020

@author: pichugin

v3 added consumption models
"""
import numpy as np
import Cyano_model_support_v4 as CMS


def Discretize(Demography):
	"""
	Computes the disctretization of the detailed demography
	Discrete[0] - 1 and 2 cells
	Discrete[1] - 4 cells
	Discrete[2] - 8 cells
	Discrete[3] - 16 cells
	Discrete[4] - 16+ cells
	Discrete[5] - 3,5,6,7,9-15 cells (non-power-of-2)
	"""
	Discrete = np.zeros(6)
	Discrete[0] = Demography[0]+2*Demography[1]
	Discrete[1] = 4*Demography[3]
	Discrete[2] = 8*Demography[7]
	Discrete[3] = 16*Demography[15]
	Discrete[4] = sum(np.multiply(np.arange(17,33), Demography[16:32]))
	Discrete[5] = 3*Demography[2] + sum( np.multiply(np.arange(5,8), Demography[4:7]) )
	Discrete[5] += sum( np.multiply(np.arange(9,16), Demography[8:15]) )
	return(Discrete)

def Dev_fig3(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
	
	
	
	""" Fig.3 processing """
	Data_3 = np.loadtxt(DataFolder+'Fig3data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig3 = 0
	for i in np.arange(len(Data_3)):
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_3[i][1]
		Comp_init = 0
		SampleTimes = np.asarray([48])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_3[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig3 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig3


def Dev_fig4(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'
	
	""" Fig.4 processing """
	Data_4 = 	np.loadtxt(DataFolder+'Fig4data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig4 = 0
	Pop_init = np.zeros(NMAX)
	Pop_init[0] = Data_4[0][2]
	Comp_init = 0
	SampleTimes = Data_4[:,1]
	SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
	for i in np.arange(1, len(SimResults[0])):
		DiscResult = Discretize(SimResults[0][i])
		ExpData = np.append(Data_4[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig4 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig4

def Dev_fig5(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'

	""" Fig.5 processing """	
	Deviation_fig5 = 0
	for replica in np.arange(3):
		Data_5 = np.loadtxt(DataFolder+'Fig5_replica_'+str(replica)+'.csv', skiprows = 1, delimiter = ',' )
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_5[0][3]
		Pop_init[3] = Data_5[0][4]/4.0
		Pop_init[7] = Data_5[0][5]/8.0
		Pop_init[15] = Data_5[0][6]/16.0
		Comp_init = 0
		SampleTimes = np.asarray([48])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_5[1][3:8],[0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig5 +=  LocalDeviation
	
	return Deviation_fig5


def Dev_fig8(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'

	""" Fig.8 processing """
	Deviation_fig8 = 0
	for replica in np.arange(18):
		Data_8 = np.loadtxt(DataFolder+'Fig8_replica_'+str(replica)+'.csv', skiprows = 1, delimiter = ',' )
		""" prepare initial compound """
		if np.isclose(Data_8[0][4], 0):
			""" fresh media """
			Comp_init = 0
		elif np.isclose(Data_8[0][4], 1):
			""" filamentation inhibitor """
			""" prepare population to harvest compound """
			H_Pop_init = np.zeros(NMAX)
			H_Pop_init[0] = 500
			H_Comp_init = 0
			H_SampleTimes = np.asarray([24])
			""" harvest """
			Comp_init = CMS.GrowthSimulation(H_Pop_init, H_Comp_init, Parameters, H_SampleTimes, Model)[1][0]
		elif np.isclose(Data_8[0][4], 2):
			""" fragmentation inducer """
			""" prepare population to harvest compound """
			H_Pop_init = np.zeros(NMAX)
			H_Pop_init[0] = 92
			H_Comp_init = 0
			H_SampleTimes = np.asarray([72])
			""" harvest """
			Comp_init = CMS.GrowthSimulation(H_Pop_init, H_Comp_init, Parameters, H_SampleTimes, Model)[1][0]
		else:
			Comp_init = float('nan')
		
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_8[0][5]
		Pop_init[3] = Data_8[0][6]/4.0
		Pop_init[7] = Data_8[0][7]/8.0
		Pop_init[15] = Data_8[0][8]/16.0
		
		SampleTimes = np.asarray([24])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)	
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_8[1][5:9],[0, 0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig8 +=  LocalDeviation
		
	return Deviation_fig8


def Dev_fig12(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'
	
	""" Fig.12 processing """
	Data_12 = 	np.loadtxt(DataFolder+'Fig12data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig12 = 0
	Pop_init = np.zeros(NMAX)
	Pop_init[0] = Data_12[0][2]
	Comp_init = 0
	SampleTimes = Data_12[:,1]
	SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
	for i in np.arange(1, len(SimResults[0])):
		DiscResult = Discretize(SimResults[0][i])
		ExpData = np.append(Data_12[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig12 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig12



""" --- consumption models --- """

def Dev_consumption_fig3(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
	""" Fig.3 processing """
	Data_3 = np.loadtxt(DataFolder+'Fig3data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig3 = 0
	for i in np.arange(len(Data_3)):
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_3[i][1]
		Comp_init = Parameters[8]
		SampleTimes = np.asarray([48])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_3[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig3 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig3


def Dev_consumption_fig4(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'
	
	""" Fig.4 processing """
	Data_4 = 	np.loadtxt(DataFolder+'Fig4data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig4 = 0
	Pop_init = np.zeros(NMAX)
	Pop_init[0] = Data_4[0][2]
	Comp_init = Parameters[8]
	SampleTimes = Data_4[:,1]
	SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
	for i in np.arange(1, len(SimResults[0])):
		DiscResult = Discretize(SimResults[0][i])
		ExpData = np.append(Data_4[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig4 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig4

def Dev_consumption_fig5(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'

	""" Fig.5 processing """	
	Deviation_fig5 = 0
	for replica in np.arange(3):
		Data_5 = np.loadtxt(DataFolder+'Fig5_replica_'+str(replica)+'.csv', skiprows = 1, delimiter = ',' )
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_5[0][3]
		Pop_init[3] = Data_5[0][4]/4.0
		Pop_init[7] = Data_5[0][5]/8.0
		Pop_init[15] = Data_5[0][6]/16.0
		Comp_init = Parameters[8]
		SampleTimes = np.asarray([48])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_5[1][3:8],[0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig5 +=  LocalDeviation
	
	return Deviation_fig5


def Dev_consumption_fig8(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'

	""" Fig.8 processing """
	Deviation_fig8 = 0
	for replica in np.arange(18):
		Data_8 = np.loadtxt(DataFolder+'Fig8_replica_'+str(replica)+'.csv', skiprows = 1, delimiter = ',' )
		""" prepare initial compound """
		if np.isclose(Data_8[0][4], 0):
			""" fresh media """
			Comp_init = Parameters[8]
		elif np.isclose(Data_8[0][4], 1):
			""" filamentation inhibitor """
			""" prepare population to harvest compound """
			H_Pop_init = np.zeros(NMAX)
			H_Pop_init[0] = 500
			H_Comp_init = Parameters[8]
			H_SampleTimes = np.asarray([24])
			""" harvest """
			Comp_init = CMS.GrowthSimulation(H_Pop_init, H_Comp_init, Parameters, H_SampleTimes, Model)[1][0]
		elif np.isclose(Data_8[0][4], 2):
			""" fragmentation inducer """
			""" prepare population to harvest compound """
			H_Pop_init = np.zeros(NMAX)
			H_Pop_init[0] = 92
			H_Comp_init = Parameters[8]
			H_SampleTimes = np.asarray([72])
			""" harvest """
			Comp_init = CMS.GrowthSimulation(H_Pop_init, H_Comp_init, Parameters, H_SampleTimes, Model)[1][0]
		else:
			Comp_init = float('nan')
		
		Pop_init = np.zeros(NMAX)
		Pop_init[0] = Data_8[0][5]
		Pop_init[3] = Data_8[0][6]/4.0
		Pop_init[7] = Data_8[0][7]/8.0
		Pop_init[15] = Data_8[0][8]/16.0
		
		SampleTimes = np.asarray([24])
		SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)	
		DiscResult = Discretize(SimResults[0][0])
		ExpData = np.append(Data_8[1][5:9],[0, 0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig8 +=  LocalDeviation
		
	return Deviation_fig8


def Dev_consumption_fig12(Parameters, NMAX, DataFolder, Model):
	"""
	Compute the discrepancy between the simulated and observed dynamics
	"""
	
#	DataFolder = '../../../Data/ExperimentalData/DataAnalysis/'
	
	""" Fig.12 processing """
	Data_12 = 	np.loadtxt(DataFolder+'Fig12data.csv', skiprows = 1, delimiter = ',' )
	Deviation_fig12 = 0
	Pop_init = np.zeros(NMAX)
	Pop_init[0] = Data_12[0][2]
	Comp_init = Parameters[8]
	SampleTimes = Data_12[:,1]
	SimResults = CMS.GrowthSimulation(Pop_init, Comp_init, Parameters, SampleTimes, Model)
	for i in np.arange(1, len(SimResults[0])):
		DiscResult = Discretize(SimResults[0][i])
		ExpData = np.append(Data_12[i][2:6],[0,0])
		LocalDeviation = sum(np.power(DiscResult-ExpData, 2))
		Deviation_fig12 +=  LocalDeviation
	
	""" there were three replicas for these data so the outcome has weight 3 """
	return 3*Deviation_fig12
