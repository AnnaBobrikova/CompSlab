#Running emcee for polarization angles

import sys
import numpy as np
from numpy import logspace, zeros, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from scipy.interpolate import interp1d
#print (np.__file__)
#import numpy.random as random
#import scipy
try:
	import matplotlib
	matplotlib.use('Agg')
	##import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	import corner
	import kombine
	from pyf_read_data import pyf_read_data
except ImportError:
	pass
import emcee
#from emcee import PTSampler
#from emcee.utils import MPIPool
import math, os
#import pymultinest
import time
from cs_ts2_func import compf
#from multiprocessing import Pool


#spath = "res/B/"#str(sys.argv[1])
print ("Outputs will be saved in (should be e.g., res/B/test_***): ") 
print (str(sys.argv[1])) 

spath = str(sys.argv[1])

NEnergy = 281
NPhase = 150
NPhase_extp = 10
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
evere=.5109989e6 # electron volts in elecron rest energy 
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
x,x_weight=IntEnergy #energies
ene_index = 118
energy_keV = x*evere/1e3 #energies in keV
#print(energy_keV[110:160])
print("The energy chosen = ", energy_keV[ene_index], " keV")
param_names = ["mass","rad","incl","theta","rho"]
#params_true = [1.4,12.0,40.0,60.0,1.0]
params_true = [1.4,12.0,40.0,60.0,10.0]
#params_true =  [ 1.0444313,  10.29601718, 40.52374309, 59.42754502,  5.28796299] #crashed once with these params
#params_true = [2.0,12.0,40.0,60.0]
#params_true = [1.6,12.0,50.0,50.0]
low_limit = [1.0, 4.0, 0.0, 0.0,1.0]
high_limit = [2.0, 18.0, 90.0, 90.0,40.0]
stlow_limit = [1.3, 11.0, 30.0, 50.0, 5.0]
sthigh_limit = [1.5, 13.0, 50.0, 70.0, 15.0]
#stlow_limit = np.copy(low_limit)
#sthigh_limit = np.copy(high_limit)
restart = True #False
#restart_file = "res/B/oblsph_emcee.dat"
restart_file = "emcee_res/oblobl_le_burst2_emcee.dat"#"res/B/oblobl_sp2_pb10_rphf_emcee.dat" #"res/B/oblobl_emcee.dat"


def shift_phase(phi,shift):
	return (phi + shift) % 1 


def readdata():
	#PulsName='res/B/B0Prho10' #used in the old results
	PulsName='pOS_pulses/lbb_rho10_sp2_f600_obl_burst2_dt'#'res/B/lbb_rho10_sp2_f600_obl_burst2_dt2_nrmr'
	#PulsName='res/B/B0P2'
	#PulsName='res/B/B0P1'
	inFlux = open(PulsName+'FF.bin')
	inphi = open(PulsName+'ff.bin')
	Flux1 = fromfile(inFlux)
	phi = fromfile(inphi)
	fluxlcurve_Iene = Flux1[0+ene_index*3:len(Flux1):3*NEnergy]
	fluxlcurve_Qene = Flux1[1+ene_index*3:len(Flux1):3*NEnergy]
	fluxlcurve_Uene = Flux1[2+ene_index*3:len(Flux1):3*NEnergy]
	Flux2 = np.array([fluxlcurve_Iene, fluxlcurve_Qene, fluxlcurve_Uene])
	return phi, Flux2


def compute_logl(phi,PA,PA_obs,phaseshift):

	sigma_tot2 = 4.0#225.0#abs(PA[t])#1.0#PA+insigma**2+(0.005*PA)**2 #(error expected/guessed in PA)**2 = 15**2 = 225, or 2**2 = 4
	norm = 0.0#0.5*log(sigma_tot2)
	phi_new = shift_phase(phi,phaseshift)
	PA_interp = interp1d(phi_new,PA,fill_value = 'extrapolate')
	PA_new = PA_interp(phi)
	loglik = 0.0
	insigma = 0.0
	phase_factor = int(NPhase/NPhase_extp)
	for t in range(NPhase_extp):
		loglik = loglik - (PA_new[t*phase_factor]-PA_obs[t*phase_factor])**2/(2.0*sigma_tot2)-norm
	return loglik


def find_best_phshift(phi,PA,PA_obs):
	#use bisection:

	phgrid = 50
	phaseshift = 0.0
	ph_min = 0.1
	ph_mid = 0.35
	ph_max = 0.5
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


#likelyhoods for emcee:
def lnprob(modelpar, low_limit, high_limit):

	#check boundaries:
	for ii in range(0,len(modelpar)):
		if(modelpar[ii]<low_limit[ii] or modelpar[ii] > high_limit[ii]):
			#print(modelpar[ii],low_limit[ii],high_limit[ii])
			#quit()
			return -np.inf

	mass = modelpar[0]
	rad = modelpar[1]
	incl = modelpar[2]
	theta = modelpar[3]
	rho = modelpar[4]

	sph = False

	#Uncomment the following if want to use the non-rotating M&R as fitting parameters for oblate star
	#if not sph:
	#	#values for rotating star:
	#	from mr_rot_nonrot import R_eq, M_obl #these assuming now that freq=600!
	#	mass0, rad0 = mass, rad
	#	rad = R_eq(mass0,rad0)
	#	mass = M_obl(mass0,rad0)

	if(mass/rad > 0.96*1.0/(2.95*1.52)): #checking causality                
		return -np.inf


	Flux = compf(mass,rad,incl,theta,rho,spherical=sph)
	#print(Flux)
	phi,Flux_obs = readdata()

	#remove the last element (that appears twice in the list)
	Flux = Flux[:,0:NPhase-1]
	Flux_obs = Flux_obs[:,0:NPhase-1]
	phi = phi[0:len(phi)-1]

	ene = ene_index #the chosen energy
	I=zeros(NPhase-1)
	Q=zeros(NPhase-1)
	U=zeros(NPhase-1)
	for t in range(NPhase-1):
		I[t],Q[t],U[t]=Flux[t,ene]
	#p=sqrt(Q**2+U**2)/I*100
	PA=arctan2(-U,-Q)*90/pi+90

	for t in range(NPhase-1):
		I[t],Q[t],U[t]=Flux_obs[0,t], Flux_obs[1,t], Flux_obs[2,t]
	PA_obs=arctan2(-U,-Q)*90/pi+90

	#apply a phase shift
	#################################################
	phi_obs_new = shift_phase(phi,0.2) #artificial phaseshift set by hand
	PA_obs_interp = interp1d(phi_obs_new,PA_obs,fill_value = 'extrapolate')
	PA_obs_new = PA_obs_interp(phi)

	#find the best phaseshift and corresponding log-likelyhood:
	phshift, loglikk = find_best_phshift(phi,PA,PA_obs_new)
	if(phshift < 0.11 or phshift > 0.49):#print when phshift is clearly wrong
		print("Phshift: ",phshift)
		print("With these parameters: ",modelpar)
	##############################################
	#loglik = compute_logl(phi,PA,PA_obs,0.0) #last argument is phase shift
	#print(loglik, loglikk)
	#quit()
	return loglikk


#Flux = compf(1.0,1.0,1.0,1.0)
print("Running emcee")

nwalkers = 20#14#50
nparams = len(param_names)
ndim = nparams
p0 = np.zeros((nwalkers,ndim))


#test call to lnprob:
#start = time.time()
##Originally a crash with these parameters:
#params_true = [ 1.43880055, 17.33520256, 38.19160194, 64.39232485, 12.23075588]
#cloglik = lnprob(params_true,low_limit,high_limit)
#print(cloglik)
##print("chi^2/d.o.f.=",-1.0*cloglik/(NPhase-nparams))
#print("chi^2/d.o.f.=",-1.0*cloglik/(NPhase_extp-nparams))
#end = time.time()
#print ("Time spend for one fit: ")
#print(end - start)
#quit()

jj = 0
# Choose an initial set of positions for the walkers.
#p0 = [np.random.rand(ndim) for i in xrange(nwalkers)]
#print(p0[:][0])

def readstart(fname,nwalk):
	#read output from emcmc/kombine run:
	npars = nparams

	try:
		full_chain= [[] for x in xrange(npars+1)]
	except NameError:
		full_chain= [[] for x in range(npars+1)]

	datafile = fname
	Nchain_size = sum(1 for line in open(datafile))
	input = open(datafile, 'r')
	lines = input.readlines()
	input.close()
	#read only the last nwalk lines:
	c_lines = Nchain_size-nwalk

	for j in range(0,len(full_chain)):
		for i in range(c_lines,Nchain_size): #not reading comment lines
			parts = lines[i].split()
			full_chain[j].append(float(parts[j]))
		parts = lines[c_lines].split()

	full_chain = np.array(full_chain)
	#remove the walker numbers from array:
	full_chain = full_chain[1:,:]
	#return the last positions of walkers:
	return full_chain[:,:].T

if(restart):
	p0 = readstart(restart_file,nwalkers)
	print ("Starting position read from file: ", restart_file) 
	#print (p0)
else:
	while(jj < nwalkers):
		#print jj
		causal = False
		for ii in range(0,ndim):
			p0[jj][ii] = stlow_limit[ii]+np.random.rand()*(sthigh_limit[ii]-stlow_limit[ii])
		#checking additional constraints: (if not fullfiled, guess new starting point to all params)
		if(p0[jj][0]/p0[jj][1] > 0.96*1.0/(2.95*1.52)): #checking causality
			continue
		jj = jj+1


print("walkers initialized succesfully")
#print(p0)

## Initialize the sampler with the chosen specs.
start = time.time()
print("Initializing the sampler")

for abc in range(0,1):
#with Pool() as pool:

	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[low_limit,high_limit],threads=2)#pool=pool)#, args=[means, icov])
	#note: threads ignored when pool used for parallelization
	#set pool to None if want to use multiprocessing parallelization
	end = time.time()
	print ("Time spend: ")
	print(end - start)

	f = open(spath + "emcee.dat", "w")
	f.close()
	start = time.time()
	## Run some steps as a burn-in.
	#pos, prob, state = sampler.run_mcmc(p0, 50)
	nsteps = 50#2#3#2#50

	#result = sampler.sample(p0, iterations=nsteps)

	#sampler.run_mcmc(p0, nsteps)

	#print("hello??")
	#quit()

	print("Start burn-in computation:")
	for i, result in enumerate(sampler.sample(p0, iterations=nsteps)): #or use sampler.run()?:
		print(i)	
		if (i+1) % 10 == 0:
			print("{0:5.1%}".format(float(i) / nsteps))
			print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
		pos, prob, state = result
		position = pos #result[0]
		#print(position)
		f = open(spath + "emcee.dat", "a")
		for k in range(position.shape[0]):
			f.write("{0:4d} {1:s}".format(k, " ".join(map(str,position[k]))))
			f.write(" ")
			f.write("{:f}\n".format(prob[k]))
		f.close()
		if(i == 0):
			end = time.time()
			print ("Time spend on the first step: ")
			print(end - start)


	## Reset the chain to remove the burn-in samples.
	print("Reset chain to remove the burn-in samples")
	sampler.reset()


	## Starting from the final position in the burn-in chain, sample for 100
	## steps.
	#sampler.run_mcmc(pos, 1000, rstate0=state) #rstate0 = The state of the random number generator.

	nsteps = 4000#2#20000
	for i, result in enumerate(sampler.sample(pos, iterations=nsteps, rstate0=state)): #or use sampler.run()?:
		if (i+1) % 100 == 0:
			print("{0:5.1%}".format(float(i) / nsteps))
			print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
		pos, prob, state = result
		position = pos #result[0]
		#print(position)
		f = open(spath + "emcee.dat", "a")
		for k in range(position.shape[0]):
			f.write("{0:4d} {1:s}".format(k, " ".join(map(str,position[k]))))
			f.write(" ")
			f.write("{:f}\n".format(prob[k]))
		f.close()


##samples = sampler.chain[:, 50:, :].reshape((-1,ndim))
samples = sampler.chain[:, 0:, :].reshape((-1,ndim))

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

limits =  zip(low_limit,high_limit)
fig = corner.corner(samples,labels=param_names,truths=params_true,range=limits)
fig.savefig(spath+"emcmc_triangle.pdf")

#print lnprob()









