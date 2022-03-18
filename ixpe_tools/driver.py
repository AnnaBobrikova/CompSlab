###############################
#This program can be used to calculate and save polarized pulse profiles from accreting millisecond pulsars
###############################
import time
from polpulse import compf

mass = 1.4 #NS mass
rad = 12.0 #NS radius
incl = 60.0 #observer inclination
theta = 20.0 #spot co-latitude
rho = 1.0 #spot angular size
pol = 0.0 #0.1171 #0.0  #USE NOW ONLY 0.0, this parameter is outdated


from numpy import logspace, zeros, fromfile, linspace
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
NEnergy = 50 
NPhase = 150
x_l, x_u = -3.7 , -1.2 # lower and upper bounds of the log_10 energy span 
evere=.5109989e6 # electron volts in elecron rest energy 
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
x,x_weight=IntEnergy #energies
phase =linspace(0,1,num=NPhase,endpoint=True,retstep=False) #input phase points
energy_keV = x*evere/1e3 # input energies in keV

Flux = compf(mass,rad,incl,theta,rho,pol,energy_keV,phase,spherical=False,antipodal=False,spath="pulses/pulse_testX",savePulse=False,atmos_path="atmos_thom/") 

import numpy as np
print("Pulse profile computation for one set of paramters finished with the following result:")
print(np.shape(Flux))
#print(Flux[0,:,0]) #Spectrum at phase 0 for Stokes I
#print(Flux[0,:,1]) #Specrtrum at phase 0 for Stokes Q
#print(Flux[0,:,2]) #Specrtrum at phase 0 for Stokes U  

print("Pulse profiles at lowest energy for I, Q, and U")
print(Flux[:,0,0]) #Pulse profile at lowest energy bin for Stokes I #sum over 0th axis to get bolometric profile
print(Flux[:,0,1]) #Pulse profile at lowest energy bin for Stokes Q 
print(Flux[:,0,2]) #Pulse profile at lowest energy bin for Stokes U 
print("phases:",phase)
print("Observed energies:",energy_keV)

