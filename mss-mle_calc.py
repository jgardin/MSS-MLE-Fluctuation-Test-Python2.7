#!/usr/bin/python2.7

import numpy as np
import sympy as sp
import argparse
from scipy.optimize import minimize
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='MSS Maximum Likelihood Calculatior for Fluctuation Tests',formatter_class=RawTextHelpFormatter)
parser.add_argument('-f', type=str, nargs=1, help='Please provide the list of observed mutants as a text file with one observed number of mutants per line. An example is the following:\n1\n0\n0\n23\n1234\n.\n.\n.\n13\n\n')
parser.add_argument('-n', type=int, nargs=1, help='Please provide the total number of cells plated for each experiment as an integer. For the MSS-MLE method to be valid, the total number of cells plated must be the same for every experiment.\n\n')
parser.add_argument('-p', type=float, nargs=1, help='Please provide the total percentage of cells plated from the grown culture as a decimal place number. For the MSS-MLE method to be valid, the total number of cells plated must be the same for every experiment.\n\n')

args=parser.parse_args()


def leeCoulson(nparray):
	median=np.median(nparray)
	x=sp.Symbol('x')
	M_est=sp.solve(sp.Eq(-x*sp.log(x) - 1.24*x + median,0),x)
	return float(M_est[0])

def ctArray(nparray,max):
	list=[0] * int(max+1)
	for i in range(int(max)+1):
		list[i]=nparray.count(i)
	return list





data=np.genfromtxt(str(*args.f),delimiter=',')
mVal=int(max(data))
ctArray_=ctArray(np.ndarray.tolist(data),mVal)

def mssCalc(estM,max=mVal,count=ctArray_):
	def rec(pi,r):
		pr=(estM/r)*sum([(pi[i]/(r-i+1)) for i in range(0,r)])
		return pr
	prod=1
	pi=[0]*(max+1)
	pi[0]=np.exp(-1.0*estM)
	for r in range(1,max+1):
		pi[r]=rec(pi,r)
		prod=prod*(pi[r]**count[r])
	
	return -1*prod

finalM=minimize	(mssCalc,leeCoulson(data),method='nelder-mead',options={'xtol':1e-3,'disp':True})
adjM=np.true_divide(float(*finalM.x),float(*args.p))
mutR=np.true_divide(float(np.log(2)*adjM),float(*args.n))
stdDev=(1.225*(np.power(adjM,-0.315)))/np.sqrt(len(data))
up95M = np.exp( np.log(adjM) + 1.96 * stdDev * ( np.power( np.exp( 1.96 * stdDev ) , -0.315 ) ) )
low95M =np.exp( np.log(adjM) - 1.96 * stdDev * ( np.power( np.exp( 1.96 * stdDev ) , 0.315 ) ) )
upMutR=np.true_divide(float(np.log(2)*up95M),float(*args.n))
lowMutR=np.true_divide(float(np.log(2)*low95M),float(*args.n))
print	"Adjusted M = ",adjM,"\n", \
	"Mutation Rate = ",mutR,"\n" \
	"95% Confidence Interval = [",upMutR," , ",lowMutR," ]\n"

