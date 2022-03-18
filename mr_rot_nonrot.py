import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *

import numpy as np 
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil

#This code applies the eqs 1-3 of Suleimanov et al 2020.


#M=1.4
#R_0=12.0
#R=R_0

nup = 600.0
#M_bar = M/1.4
#nu_cr = 1278*sqrt(M_bar)*(10/R)**1.5
#nu_bar = nup/nu_cr  
#R_e = R_0*(0.9766 + 0.025/(1.07- nu_bar)+0.07*M_bar**1.5*nu_bar**2)

#a1 = 0.001*M_bar**1.5
#a0 = 1.0 - a1/1.1
#a2 = 10*a1
#M_prime = M*(a0 + a1/(1.1-nu_bar)+ a2*nu_bar**2)
##print(M_prime/M)
##print(R_e/R_0)
#R_g=M_prime*2.95

#Solve the equation to the other direction:

def R_eq(M,R):
	M_bar = M/1.4
	nu_cr = 1278*sqrt(M_bar)*(10/R)**1.5
	nu_bar = nup/nu_cr  
	#print(nup,nu_cr,nu_bar)
	R_e = R*(0.9766 + 0.025/(1.07- nu_bar)+0.07*M_bar**1.5*nu_bar**2)
	if(nu_bar > 0.99):
		return 1e-5
	return R_e

def M_obl(M,R):
	M_bar = M/1.4
	nu_cr = 1278*sqrt(M_bar)*(10/R)**1.5
	nu_bar = nup/nu_cr 
	a1 = 0.001*M_bar**1.5
	a0 = 1.0 - a1/1.1
	a2 = 10*a1
	M_prime = M*(a0 + a1/(1.1-nu_bar)+ a2*nu_bar**2)
	if(nu_bar > 0.99):
		return 1e10
	return M_prime


#from sympy import symbols, Eq, solve
#x, y = symbols('M R')
#eq1 = Eq(R_eq(M,R)-12.0,0)
#eq2 = Eq(M_obl(M,R)-1.4,0)
#res = solve((eq1,eq2), (M, R))
#print(res)


#def func(x):
#	return x + 2*cos(x)

#def func2(x):
#	out = [x[0]*cos(x[1]) - 4]
#	out.append(x[1]*x[0] - x[1] - 5)
#	return out

def func_2eq(x):
	out = [R_eq(x[0],x[1]) - 12.0]
	out.append(M_obl(x[0],x[1]) - 1.4)
	return out


from scipy.optimize import fsolve
#x0 = fsolve(func, 0.3)
#print(x0)
#x02 = fsolve(func2, [1, 1])
#print(x02)

#mr = fsolve(func_2eq, [1.4, 12.0])
#m_nrot, r_nrot = mr[0], mr[1]
#print(mr)

#checking that the solution is correct
#print(M_obl(m_nrot,r_nrot),R_eq(m_nrot,r_nrot))







