
import sys
import numpy as np
from numpy import logspace, zeros, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from scipy.interpolate import interp1d


def shift_phase(phi,shift):
	return (phi + shift) % 1 

def compute_logl(phi,PA,PA_obs,phaseshift):

	sigma_tot2 = 225.0#4.0#225.0#abs(PA[t])#1.0#PA+insigma**2+(0.005*PA)**2 #(error expected/guessed in PA)**2 = 15**2 = 225, or 2**2 = 4
	norm = 0.0#0.5*log(sigma_tot2)
	phi_new = shift_phase(phi,phaseshift)
	PA_interp = interp1d(phi_new,PA,fill_value = 'extrapolate')
	PA_new = PA_interp(phi)
	loglik = 0.0
	insigma = 0.0
	NPhase = len(phi)
	#print(PA_new)
	for t in range(NPhase):
		#ignore points where 180 shift might happen, and therefore large errors:
		#print(abs(PA_new[t]-PA_obs[t]))
		if(abs(PA_new[t]-PA_obs[t]) < 30.0):
			loglik = loglik - (PA_new[t]-PA_obs[t])**2/(2.0*sigma_tot2)-norm
	return loglik


def find_best_phshift(phi,PA,phi_obs,PA_obs):
	#use bisection:

	#interpolate PA_obs to have the same resolution than PA and phi:
	PA_obs_interp = interp1d(phi_obs,PA_obs,fill_value='extrapolate')# workd with newer scipy
	PA_obs = PA_obs_interp(phi) 

	#print(PA,phi)
	#print(PA_obs,phi_obs)
	#quit()

	phgrid = 100#50
	phaseshift = 0.0
	ph_min = 0.9#0.0#0.1 #Hacks helping to find correct point
	ph_mid = 0.95#0.5#0.35
	ph_max = 1.0#0.5
	gf_min = compute_logl(phi,PA,PA_obs,ph_min)
	gf_mid = compute_logl(phi,PA,PA_obs,ph_mid)
	gf_max = compute_logl(phi,PA,PA_obs,ph_max)

	ph_grid = np.zeros((phgrid+3))
	gf_grid = np.zeros((phgrid+3))
	#update grids
	ph_grid[0] = ph_min
	gf_grid[0] = gf_min
	ph_grid[1] = ph_mid
	gf_grid[1] = gf_mid
	ph_grid[2] = ph_max
	gf_grid[2] = gf_max

	for ishift in range(0,phgrid):
		if (ishift%2 == 0):
			phaseshift = 0.5*(ph_min + ph_mid)
			gf =  compute_logl(phi,PA,PA_obs,phaseshift)
			ph_grid[ishift+3] = phaseshift
			gf_grid[ishift+3] = gf

			if (gf < gf_mid):
				ph_min = phaseshift
				gf_min = gf
			else:
				ph_max = ph_mid
				gf_max = gf_mid
				ph_mid = phaseshift
				gf_mid = gf
		else:
			phaseshift = 0.5*(ph_mid + ph_max)
			gf = compute_logl(phi,PA,PA_obs,phaseshift)                  
			ph_grid[ishift+3] = phaseshift
			gf_grid[ishift+3] = gf

			if (gf < gf_mid):
				ph_max = phaseshift
				gf_max = gf
			else:
				ph_min = ph_mid
				gf_min = gf_mid
				ph_mid = phaseshift
				gf_mid = gf
		#print(phaseshift,gf)

	best_phaseshift = ph_mid

	return best_phaseshift, gf
	#return 0.1
