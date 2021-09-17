#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 08:56:39 2020

@author: pichugin

v3.2 added new models
v3.3 added new models
v4   added consumption models
"""

import numpy as np
#import csv
#import matplotlib.pyplot as plt
#import pandas as pd
#import seaborn as sns
from copy import deepcopy
from scipy.integrate import solve_ivp


def ExponentialSigmoid(x):
	"""
	raw sigmoid
	"""
	x = max(x, -100)
	return 1.0/(1.0+np.exp(-x))

def Sigmoid(x, H, x_0, sigma):
	""" 
	Computes the exponential sigmoid S(x) with boundary conditions:
	S(0) = 0
	S(+\infty) = H
	
	Arguments:
		(float) x - value of the argument
		(float) H - limit value at the plus infinity
		(float) x_0 - the position of the inflection point
		(float) sigma - the smoothness of the sigmoid
	"""
	if x_0 / sigma > -20:
		Numerator = ExponentialSigmoid((x-x_0)/sigma) - ExponentialSigmoid(-x_0 / sigma)
		Denominator = 1 - ExponentialSigmoid(-x_0 / sigma)
	else:
		Numerator = 1 - np.exp(-x/sigma)
		Denominator = 1
	
	return H*Numerator/Denominator

""" --- Models --- """

def CoreProjectionMatrix(Pop, Parameters):
	""" generate projection matrix of a population without any fragmentation effects 
		this part is the common first step of matrix generation for all models
	"""
	W = Parameters[0]
	K = Parameters[1]
	
	NMAX = len(Pop)
	A = np.zeros((NMAX, NMAX))
	
	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
	EffectiveW = W * (1 - TotalCellCount / K)

	for i in np.arange(NMAX):
		A[i,i] -= EffectiveW 
	""" maximal size stop growing """
	A[NMAX-1, NMAX-1] += EffectiveW
		
	for i in np.arange(int((NMAX+1)/2+1)):
		if 2*(i+1)-1 < NMAX:
			A[2*(i+1)-1, i] += EffectiveW	
	return A

def GenericCleavageModel(Pop, A, CleavageRate):
	"""
	Fills the cleavage model matrix with fragmentation terms for a given CleavageRate
	This is a shared step for all cleavage models
	"""
	NMAX = len(Pop)
	for i in np.arange(NMAX):
		A[i,i] -= (i * CleavageRate )
		
	for i in np.arange(NMAX):
		for j in np.arange(i+1, NMAX):
			A[i, j] += 2*(CleavageRate)	
	return A

def GenericDeathModel(Pop, A, ToxicityRate):
	"""
	Fills the death model matrix with fragmentation terms for a given ToxicityRate
	This is a shared step for all death models
	"""
	NMAX = len(Pop)
	for i in np.arange(NMAX):
		A[i,i] -= ((i+1) * ToxicityRate )
		
	for i in np.arange(NMAX):
		for j in np.arange(i+1, NMAX):
			A[i, j] += 2*(ToxicityRate)	
	return A

""" Cleavage set of models """

def LinearCleavageModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Slope = Parameters[3]
	CleavageRate = (CompConc > Conc_treshold)*(Slope * (CompConc - Conc_treshold))

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A


def StepCleavageModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = (CompConc > Conc_treshold)*Max_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def CappedCleavageModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = Max_effect * np.min([ CompConc/Conc_treshold, 1 ])

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def ThresholdCleavageModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Min_effect = Parameters[3]
	CleavageRate = (CompConc > Conc_treshold) * Min_effect * CompConc/Conc_treshold

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def ProportionalCleavageModel(Pop, CompConc, Parameters):
	
	Slope = Parameters[3]
#	Min_effect = Parameters[3]
	CleavageRate = Slope * CompConc

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def ConstantCleavageModel(Pop, CompConc, Parameters):
	
	Effect = Parameters[3]
#	Min_effect = Parameters[3]
	CleavageRate = Effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def BottomCappedCleavageModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = Max_effect * np.max([ CompConc/Conc_treshold, 1 ])

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
	
	return A

def GeneralLinearCleavageModel(Pop, CompConc, Parameters):
	
	Slope = Parameters[2]
	Min_effect = Parameters[3]
	CleavageRate = Slope * CompConc + Min_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def SigmoidCleavageModel(Pop, CompConc, Parameters):
	
	Max_effect = Parameters[2]
	TransitionWidth = Parameters[3]
	TransitionConc = Parameters[4]
#	Slope = Parameters[3]
#	Min_effect = Parameters[4]
	CleavageRate = Max_effect/(1 + np.exp(-(CompConc - TransitionConc)/TransitionWidth))

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def MichaelisMentenCleavageModel(Pop, CompConc, Parameters):
	
	Kin_const = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = Max_effect*CompConc/(CompConc + Kin_const)

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def QuadraticCleavageModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	CleavageRate = Min_effect + (CompConc/Scale)**2

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def SaturatedExponentCleavageModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	Max_effect = Parameters[4]
	CleavageRate = Max_effect - (Max_effect - Min_effect)*np.exp(- CompConc/Scale)

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A





""" Death set of models """



def LinearDeathModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Slope = Parameters[3]
	ToxicityRate = (CompConc > Conc_treshold)*(Slope * (CompConc - Conc_treshold))

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A


def StepDeathModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	ToxicityRate = (CompConc > Conc_treshold)*Max_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def CappedDeathModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	ToxicityRate = Max_effect * np.min([ CompConc/Conc_treshold, 1 ])

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def ThresholdDeathModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Min_effect = Parameters[3]
	ToxicityRate = (CompConc > Conc_treshold) * Min_effect * CompConc/Conc_treshold

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def ProportionalDeathModel(Pop, CompConc, Parameters):
	
	Slope = Parameters[3]
#	Min_effect = Parameters[3]
	ToxicityRate = Slope * CompConc

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def ConstantDeathModel(Pop, CompConc, Parameters):
	
	Effect = Parameters[3]
#	Min_effect = Parameters[3]
	ToxicityRate = Effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def BottomCappedDeathModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	ToxicityRate = Max_effect * np.max([ CompConc/Conc_treshold, 1 ])

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def GeneralLinearDeathModel(Pop, CompConc, Parameters):
	
	Slope = Parameters[2]
	Min_effect = Parameters[3]
	ToxicityRate = Slope * CompConc + Min_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def SigmoidDeathModel(Pop, CompConc, Parameters):
	
	Max_effect = Parameters[2]
	TransitionWidth = Parameters[3]
	TransitionConc = Parameters[4]
	
	ToxicityRate = Max_effect/(1 + np.exp(-(CompConc - TransitionConc)/TransitionWidth))

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def MichaelisMentenDeathModel(Pop, CompConc, Parameters):
	
	Kin_const = Parameters[2]
	Max_effect = Parameters[3]
	ToxicityRate = Max_effect*CompConc/(CompConc + Kin_const)
	
	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def QuadraticDeathModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	ToxicityRate = Min_effect + (CompConc / Scale)**2
	
	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A

def SaturatedExponentDeathModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	Max_effect = Parameters[4]
	ToxicityRate = Max_effect - (Max_effect - Min_effect)*np.exp(- CompConc/Scale)

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericDeathModel(Pop, A, ToxicityRate)
		
	return A









""" Consumption set of models """

def GeneralLinearConsumptionModel(Pop, CompConc, Parameters):
	
	Slope = Parameters[2]
	Min_effect = Parameters[3]
	Max_concentration = Parameters[8]
	CleavageRate = Slope * (Max_concentration - CompConc) + Min_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def StepConsumptionModel(Pop, CompConc, Parameters):
	
	Conc_treshold = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = (CompConc < Conc_treshold)*Max_effect

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A


def InverseConsumptionModel(Pop, CompConc, Parameters):
	
	Conc_level = Parameters[2]
	Scale = Parameters[3]
	CleavageRate = Scale/(1.0 + CompConc/Conc_level)

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def QuadraticConvexConsumptionModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	CleavageRate = Min_effect + Scale*(CompConc -1.0)**2

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def QuadraticConcaveConsumptionModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Max_effect = Parameters[3]
	CleavageRate = np.max([0, Max_effect - Scale*CompConc**2])

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def SigmoidConsumptionModel(Pop, CompConc, Parameters):
	
	Max_effect = Parameters[2]
	TransitionWidth = Parameters[3]
	TransitionConc = Parameters[4]
#	Slope = Parameters[3]
#	Min_effect = Parameters[4]
	CleavageRate = Max_effect/(1 + np.exp((CompConc - TransitionConc)/TransitionWidth))

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A

def DecayingExponentConsumptionModel(Pop, CompConc, Parameters):
	
	Scale = Parameters[2]
	Min_effect = Parameters[3]
	Max_effect = Parameters[4]
	CleavageRate = Min_effect + (Max_effect - Min_effect)*np.exp(- CompConc/Scale)

	A = CoreProjectionMatrix(Pop, Parameters)
	A = GenericCleavageModel(Pop, A, CleavageRate)
		
	return A


""" --- old models to delete --- """
#def LinearDeathProjectionMatrix(Pop, CompConc, Parameters):
#	
#	W = Parameters[0]
#	K = Parameters[1]
#	d_intercept = Parameters[2]
#	d_slope = Parameters[3]
#	
##	H_tox = Parameters[2]
##	x0_tox = Parameters[3]
##	sigma_tox = Parameters[4]
##	H_clv = Parameters[5]
##	x0_clv = Parameters[6]
##	sigma_clv = Parameters[7]
#
#	ToxicityRate = d_intercept + d_slope * CompConc
##	ToxicityRate = Sigmoid(CompConc, H_tox, x0_tox, sigma_tox)
##	CleavageRate = Sigmoid(CompConc, H_clv, x0_clv, sigma_clv)
#	
#	NMAX = len(Pop)
#	A = np.zeros((NMAX, NMAX))
#	
#	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
#	EffectiveW = W * (1 - TotalCellCount / K)
##	EffectiveW = (EffectiveW + abs(EffectiveW))/2.0
#	
#	
#	for i in np.arange(NMAX):
#		A[i,i] -= (EffectiveW + (i+1) * ToxicityRate )
#	""" maximal size stop growing """
#	A[NMAX-1, NMAX-1] += EffectiveW
#		
#	for i in np.arange(int((NMAX+1)/2+1)):
#		if 2*(i+1)-1 < NMAX:
#			A[2*(i+1)-1, i] += EffectiveW
#		
#	for i in np.arange(NMAX):
#		for j in np.arange(i+1, NMAX):
#			A[i, j] += 2*(ToxicityRate)
#		
#	return A	
#
#def StepDeathProjectionMatrix(Pop, CompConc, Parameters):
#	
#	W = Parameters[0]
#	K = Parameters[1]
#	Conc_treshold = Parameters[2]
#	Max_effect = Parameters[3]
#	
##	H_tox = Parameters[2]
##	x0_tox = Parameters[3]
##	sigma_tox = Parameters[4]
##	H_clv = Parameters[5]
##	x0_clv = Parameters[6]
##	sigma_clv = Parameters[7]
#
#	ToxicityRate = (CompConc > Conc_treshold) * Max_effect
##	ToxicityRate = Sigmoid(CompConc, H_tox, x0_tox, sigma_tox)
##	CleavageRate = Sigmoid(CompConc, H_clv, x0_clv, sigma_clv)
#	
#	NMAX = len(Pop)
#	A = np.zeros((NMAX, NMAX))
#	
#	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
#	EffectiveW = W * (1 - TotalCellCount / K)
##	EffectiveW = (EffectiveW + abs(EffectiveW))/2.0
#	
#	
#	for i in np.arange(NMAX):
#		A[i,i] -= (EffectiveW + (i+1) * ToxicityRate )
#	""" maximal size stop growing """
#	A[NMAX-1, NMAX-1] += EffectiveW
#		
#	for i in np.arange(int((NMAX+1)/2+1)):
#		if 2*(i+1)-1 < NMAX:
#			A[2*(i+1)-1, i] += EffectiveW
#		
#	for i in np.arange(NMAX):
#		for j in np.arange(i+1, NMAX):
#			A[i, j] += 2*(ToxicityRate)
#		
#	return A
#
#def LinearCleavageProjectionMatrix(Pop, CompConc, Parameters):
#	
#	W = Parameters[0]
#	K = Parameters[1]
#	c_intercept = Parameters[2]
#	c_slope = Parameters[3]
#
#	CleavageRate = c_intercept + c_slope * CompConc
##	ToxicityRate = Sigmoid(CompConc, H_tox, x0_tox, sigma_tox)
##	CleavageRate = Sigmoid(CompConc, H_clv, x0_clv, sigma_clv)
#	
#	NMAX = len(Pop)
#	A = np.zeros((NMAX, NMAX))
#	
#	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
#	EffectiveW = W * (1 - TotalCellCount / K)
##	EffectiveW = (EffectiveW + abs(EffectiveW))/2.0
#	
#	
#	for i in np.arange(NMAX):
#		A[i,i] -= (EffectiveW + i * CleavageRate )
#	""" maximal size stop growing """
#	A[NMAX-1, NMAX-1] += EffectiveW
#		
#	for i in np.arange(int((NMAX+1)/2+1)):
#		if 2*(i+1)-1 < NMAX:
#			A[2*(i+1)-1, i] += EffectiveW
#		
#	for i in np.arange(NMAX):
#		for j in np.arange(i+1, NMAX):
#			A[i, j] += 2*(CleavageRate)
#		
#	return A
#
#
#def StepCleavageProjectionMatrix(Pop, CompConc, Parameters):
#	
#	W = Parameters[0]
#	K = Parameters[1]
#	Conc_treshold = Parameters[2]
#	Max_effect = Parameters[3]
#
#	CleavageRate = (CompConc > Conc_treshold) * Max_effect
##	ToxicityRate = Sigmoid(CompConc, H_tox, x0_tox, sigma_tox)
##	CleavageRate = Sigmoid(CompConc, H_clv, x0_clv, sigma_clv)
#	
#	NMAX = len(Pop)
#	A = np.zeros((NMAX, NMAX))
#	
#	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
#	EffectiveW = W * (1 - TotalCellCount / K)
##	EffectiveW = (EffectiveW + abs(EffectiveW))/2.0
#	
#	
#	for i in np.arange(NMAX):
#		A[i,i] -= (EffectiveW + i * CleavageRate )
#	""" maximal size stop growing """
#	A[NMAX-1, NMAX-1] += EffectiveW
#		
#	for i in np.arange(int((NMAX+1)/2+1)):
#		if 2*(i+1)-1 < NMAX:
#			A[2*(i+1)-1, i] += EffectiveW
#		
#	for i in np.arange(NMAX):
#		for j in np.arange(i+1, NMAX):
#			A[i, j] += 2*(CleavageRate)
#		
#	return A
#
#
#def ComputeProjectionMatrix(Pop, CompConc, Parameters):
#	
#	W = Parameters[0]
#	K = Parameters[1]
#	H_tox = Parameters[2]
#	x0_tox = Parameters[3]
#	sigma_tox = Parameters[4]
#	H_clv = Parameters[5]
#	x0_clv = Parameters[6]
#	sigma_clv = Parameters[7]
#
#	
#	ToxicityRate = Sigmoid(CompConc, H_tox, x0_tox, sigma_tox)
#	CleavageRate = Sigmoid(CompConc, H_clv, x0_clv, sigma_clv)
#	
#	NMAX = len(Pop)
#	A = np.zeros((NMAX, NMAX))
#	
#	TotalCellCount = sum( np.arange(1, NMAX+1) * Pop )
#	EffectiveW = W * (1 - TotalCellCount / K)
##	EffectiveW = (EffectiveW + abs(EffectiveW))/2.0
#	
#	
#	for i in np.arange(NMAX):
#		A[i,i] -= (EffectiveW + (i+1) * ToxicityRate + i * CleavageRate )
#	""" maximal size stop growing """
#	A[NMAX-1, NMAX-1] += EffectiveW
#		
#	for i in np.arange(int((NMAX+1)/2+1)):
#		if 2*(i+1)-1 < NMAX:
#			A[2*(i+1)-1, i] += EffectiveW
#		
#	for i in np.arange(NMAX):
#		for j in np.arange(i+1, NMAX):
#			A[i, j] += 2*(ToxicityRate + CleavageRate)
#		
#	return A
""" ---Simulations--- """


def Diff_Eq_RHS(t, y, ProjectionMatrixModel, Parameters):
	"""
	Computes the right hand side of differential equation and has a format required by ODE solver
	Arguments:
		(float) t - time (not used but needed by the format)
		(list of floats) y - dynamic variables: 
							y[0:-1] - Pop
							y[-1] - Comp_conc
		(callable) ProjectionMatrixModel - which model to use for projection matrix computation
		(list of floats) Params - parameters of the dynamics
	"""
	NMAX = len(y) - 1
	P_comp = Parameters[8]
	D_comp = Parameters[9]
	
	A = ProjectionMatrixModel(y[0:-1], y[-1], Parameters)
	dPop = np.dot(A, y[0:-1])
	dComp_conc = ( P_comp * sum(np.arange(1,NMAX+1) * y[0:-1]) - D_comp*y[-1] )
	
	RHS = np.append(dPop, dComp_conc)
	
	return RHS

def Diff_Eq_RHS_consumption(t, y, ProjectionMatrixModel, Parameters):
	"""
	Computes the right hand side of differential equation and has a format required by ODE solver
	Arguments:
		(float) t - time (not used but needed by the format)
		(list of floats) y - dynamic variables: 
							y[0:-1] - Pop
							y[-1] - Comp_conc
		(callable) ProjectionMatrixModel - which model to use for projection matrix computation
		(list of floats) Params - parameters of the dynamics
	"""
	NMAX = len(y) - 1
#	P_comp = Parameters[8]
	D_comp = Parameters[9]
	
	A = ProjectionMatrixModel(y[0:-1], y[-1], Parameters)
	dPop = np.dot(A, y[0:-1])
	dComp_conc = - D_comp * np.max([0, y[-1]]) * sum(np.arange(1,NMAX+1) * y[0:-1])
	
	RHS = np.append(dPop, dComp_conc)
	
	return RHS


def GrowthSimulation(Pop_init, Comp_init, Parameters, RecordTimes, ProjectionMatrixModel):
	"""
	Simulate a single round of growth
	Parameters: 
		(numpy array) Pop_init - initial state of population
		(float) Comp_init - initial concentration of compound
		(list of floats) Params - parameters of the dynamics
		(numpy array) RecordTimes - times at which the record should be taken
		(function name) ProjectionMatrixModel - which model to use for projection matrix computation
	"""
	dt = 0.1
	TMAX = max(RecordTimes) + 2*dt	
	tspan = [0, TMAX]
	
	y0 = np.append(Pop_init, Comp_init)
	
	if Parameters[7]==0:
		Differential_equation = Diff_Eq_RHS
	elif Parameters[7]==1:
		Differential_equation = Diff_Eq_RHS_consumption
	else:
		print('The type of model is not specified')
		Differential_equation = 'None'
	
	solution = solve_ivp(Differential_equation, tspan, y0, t_eval = RecordTimes, args = (ProjectionMatrixModel, Parameters))
	
	""" manage output """
	NMAX = len(Pop_init)
	PopRecord = np.zeros((len(RecordTimes), NMAX))
	CompRecord = np.zeros(len(RecordTimes))
	
	PopRecord = np.transpose(solution.y[0:-1])
	CompRecord = solution.y[-1]

	return [PopRecord, CompRecord]
