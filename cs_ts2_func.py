#!/usr/bin/python

# coding: utf-8

# In[17]:



# timer:
#from time import time
#time0=time()
#def pr(s,gn,message = ''): 
#      #print(gn,' :'+message+': ',time()-s)
#      #return time(),gn+1
#       return time0,gn+1
#s,gn=pr(time0,0)

Spectrum={
      1:'Compton',
      2:'Burst',
      3:'Thomson',
      4:'Bbody',
      0:'FromFile'
#}[0]
#}[4]
}[2]

oblateness='AlGendy'

AtmName='res/B/B0_new' # the prefix for all result files related to the set of parameters
#PulsName=AtmName+'P2'
#PulsName='res/B/B0P1'
PulsName='res/B/lbb_rho10x_sp2_f600_obl_burst2_dt2_test'#B0Prho10sphsp2_test'#lbb_rho10_sp2_f600_sph_test_dmr2' #B0Ptest'
#PulsName='res/B/B0Prho10'
computePulse= True
plotAtm=False#True
plotPulse=True


# In[18]:exit()


#import:
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
# m=fromfile(AtmName+'m.bin')
# print(m)
# exit()
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil, argmin
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d#,CubicSpline 
from scipy.interpolate import CubicSpline 
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect

#import numpy as np
#import matplotlib.pyplot as plt
#s,gn=pr(s,gn, 'importing')

def Planck(x,T):
      """   Planck function for Intensity of black body radiation 
      The only argument x is the energy of a photon in units of electron rest energy ( h \\nu / m_e c^2 ) 
      The photon temperature is given by T also in units of electron rest mass
      Planck returns the intensity of  BB radiation
      """
      e=x/T
      C=1.  # some dimension constant.
      R=C*e*e # Rayleigh Jeans law'
      I=.0 if e>5e2 else R*e/(exp(e)-1.) if e > 1e-5 else R
      return I


def compf(mass,eqrad,incl_deg,theta_deg,rho_deg,spherical=False):




	# In[19]:




	colors=['xkcd:brownish red',
		'xkcd:red',
		'xkcd:orange',
		'xkcd:dark yellow',
		'xkcd:dark yellow green',
		'xkcd:deep green',
		'xkcd:dark cyan',
		'xkcd:blue',
		'xkcd:purple'   
	] # 'rygcbm' # Rainbow 
	NColors=len(colors)

	#physical constants:
	evere=.5109989e6 # electron volts in elecron rest energy 
	G=13275412528e1 # G*M_sol in km^3/s^2 
	c=299792458e-3 # speed of light in km/s

	# Atmosphere parameters: 
	tau_T= .6 # Thomson optical depth of thermalization 
	x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
	Theta = 0.1  # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
	T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2 #~  1.0219978 keV

	#precomputations :
	ScatterNum = 23 # total number of scatterings
	NGamma= 7# number of Lorenz factor points (\gamma)
	NAzimuth= 10 # 12 # numbers of azimuth angles (\phi) [0,pi]
	NEnergy = 281 # 50# 101 # number of energy points (x)
	NDepth = 100# 101  # number of optical depth levels (\tau)
	NMu = 22 # 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
	NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]

	IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
	IntAzimuth = leggauss(NAzimuth*2) # sample points and weights for computing azimuth-averaged matrix
	IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
	IntDepth = linspace(0,tau_T,num=NDepth,retstep=True)  # sample points and weights for integrations over the optical depth computing intencity 
	IntZenith = leggauss(NZenith) #  sample points and weights for integrations over zenith angle in positive and negative directions together

	# IntZenith = linspace(1e-6-1,1-1e-6,num=NZenith), array([2.]+[2.]*(NZenith-2)+[2.])/(NZenith)
	# IntZenith = linspace(1/NZenith - 1,1 - 1/NZenith,num=NZenith), array([2.]*(NZenith))/(NZenith)
	# mu,mu_weight=IntZenith
	# mu=array( ([1-1e-7]+list(mu)+[1-1e-7]) )
	# mu_weight=array([1e-5]+list(mu_weight)+[1e-5])
	# IntZenith=(mu, mu_weight)
	# NZenith=NZenith+2
	# NMu=NMu+1


	K2Y = kn(2,1./Theta) # second modified Bessel function of reversed dimensionless temperature       

	mu,mu_weight=IntZenith
	x,x_weight=IntEnergy
	tau,tau_weight=IntDepth

	outParams = open(AtmName+'.dat','w')
	outParams.write(Spectrum+'\n')
	outParams.write(str(tau_T)+' = total depth tau_T\n')
	outParams.write(str(x_l)+' '+str(x_u)+' = log_10 energy span from to\n')
	outParams.write(str(Theta)+' = dimensionless electron gas temperature in m_ec^2\n')
	outParams.write(str(T)+' = dimensionless initial bb photons temperature in m_ec^2\n')
	outParams.write(str(ScatterNum)+' = Number of scatterings computed\n')
	outParams.write(str(NMu)+' = NMu, NZenith=2 NMu, number of points in mu grid\n')
	outParams.write(str(NEnergy)+' = NEnergy, number of points in x grid\n')
	outParams.write(str(NGamma)+' '+str(NAzimuth)+' '+str(NDepth)+' = NGamma NAzimuth NDepth parameters\n')



	#s,gn=pr(s,gn,'precomps')

	# In[ ]:
	if Spectrum=='Bbody' : # Initializing Stokes vectors arrays, computing zeroth scattering 
		Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
		for e in range(NEnergy):
			#print("energy, T =", x[e]*evere/1e3,T*evere/1e3)
			E=x[e]/T
			for d in range(NZenith):
				constbb= 5.039617e22 
				ex=exp(-x[e]/T)
				Intensity[e][d][0]=constbb*x[e]**3*ex/(1.0-ex)
				Intensity[e][d][1]=0.0
				#s,gn=pr(s,gn,'I0')




	if Spectrum=='Burst' : # Initializing Stokes vectors arrays, computing zeroth scattering 
		Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
		for e in range(NEnergy):
			# print(e,log(Planck(x[e-1])/Planck(x[e]))/log(x[e]/x[e-1]),'           ',log(x[e])/log(10),'  ')
			for d in range(NZenith):
				#Intensity[e,d,0]=Planck(x[e],T)#*(1 + 2.06*mu[d])#TS TESTING BLACKBODY energy spectrum with burst-beaming only in polarization
                                #burst beaming approx normalized so that int over mus from 0 to 1 gives same as it will for istropic case (1/2):
				#Intensity[e,d,0]=Planck(x[e],T)*(3.0/7.0)*(1.0+2.0*mu[d])
				Intensity[e,d,0]=Planck(x[e],T)*(0.421+0.868*mu[d]) #more accurate version
				Intensity[e,d,1]=Intensity[e,d,0]*0.1171*(mu[d] - 1.)/(1. + 3.582*mu[d])#abs(mu[d]))
				#Intensity[e,d,1]=Intensity[e,d,0]#*0.1171*(mu[d] - 1.)/(1. + 3.82*abs(mu[d]))
		#s,gn=pr(s,gn,'I0')
		# exit()



	# In[ ]:




	if Spectrum=='FromFile' :
		inI = open(AtmName+'I.bin')
		inx = open(AtmName+'x.bin')
		inm = open(AtmName+'m.bin')
		x=fromfile(inx)
		mu=fromfile(inm)
		NEnergy=len(x)
		NZenith=len(mu)
		NMu=NZenith//2
		Intensity=fromfile(inI).reshape((NEnergy,NZenith,2))
		#s,gn=pr(s,gn,'I is read')
	else: #Make sure this save wont overwrite your old tables!
		pass
		#print("Not saving the spectrum")
		#outI = open(AtmName+'I.bin','w')
		#outx = open(AtmName+'x.bin','w')
		#outm = open(AtmName+'m.bin','w')
		#Intensity.tofile(outI,format="%e")
		#x.tofile(outx,format="%e")
		#mu.tofile(outm,format="%e")



	# In[28]:
	if computePulse:


		NPhi = 120 #500 #120 # Number of equidistant phase points
		NPhase = 150 #500# 150 # Number of observation phases
		NBend= 20 # Number of knots in light bending integrations
		NAlpha= 200#1000 # 10000 # Number of psi/aplha grid points 
		IntBend = leggauss(NBend)
		NZenithBig=100
		#NZenithBig = NZenith

		phi=linspace(0,2*pi,num=NPhi,endpoint=False,retstep=False)
		phase =linspace(0,1,num=NPhase,endpoint=True,retstep=False)
		phase_obs=zeros(NPhi)
		nu=350.0 #600.0#1.0#100 #600 # star rotation frequency in Hz
		#M=1.4 # star mass in solar masses
		M=mass #input param
		R_g=M*2.95325 # gravitational Schwarzschild radius #TS: Made this more accurate
		#R_e=12.0 # equatorial radius of the star in kilometers
		R_e=eqrad #input param
      
		NRho=20#40#20#4#2#8
		NVarphi=20#40#20#6#4

		# IntVarphi = linspace(0,2*pi,num=NVarphi,endpoint=False,retstep=True)
		IntVarphi = leggauss(NVarphi)
		IntRho = leggauss(NRho)
	      
		if oblateness=='AlGendy': # from AlGendy et. al. (2014)
			Omega_bar=2*pi*nu*sqrt(2*R_e**3/R_g)/c
			#print('_O_^2',Omega_bar**2,'_O_',Omega_bar)
			flattening=(0.788-0.515*R_g/R_e)*Omega_bar**2 
			#print(R_e*(1-flattening))
		elif oblateness=='Sphere':
			flattening=0.0
		else:
			flattening=oblateness
		#exit()
		#print(flattening)
		if(spherical):#TS: from input param, 0.0 not working in this code
			#print("spherical star")
			flattening = 1e-8


		outParams = open(PulsName+'.dat','w')
		outParams.write(AtmName+'.dat is the name of file with some corresponding atmosphere model\n')
		outParams.write(str(R_e)+' = equatorial radius R_e\n')
		outParams.write(str(M)+' = star mass M in solar masses\n')
		outParams.write(str(nu)+' = star rotation frequency nu in Hz\n')
		outParams.write(str(NPhi)+' '+str(NBend)+' '+str(NAlpha)+' = NPhi NBend NAlpha parameters\n')



		def Beloborodov(cos_psi):
		    """Beloborodov's approximation for cos_alpha(cos_psi) light bending function
		    takes the cos psi 
		    returns the cos alpha and its derivative
		    """
		    return 1. + (cos_psi - 1.)/redshift**2 ,1./redshift**2

		def Schwarzschild(R,alpha):
		    """Schwarzschild exact relation between the \psi and \\alpha angles, where
		    \\alpha is the angle between radius vector of the spot and the direction of the outgoing photon near the surface
		    and \psi is the angle between normal and light propagation at the limit of infinite distance.
		    For given distance from the mass center and the emission angle \\alpha 
		    this function returns two numbers: 
		          the corresponding angle \psi 
		          and the time lag over against the fotons emited with zero impact parameter at the radius.
		    """
		    kx,wx=IntBend
		    eps=(1+kx[0])/4e2
		    u=R_g/R 
		    b=sin(alpha)/sqrt(1-u)*R # impact parameter
		    if 2*alpha>pi+eps:
		          cos_3eta=sqrt(27)*R_g/2/b
		          if cos_3eta > 1:
		                return pi+2*eps,0 # the timelag 
		          closest_approach=-2*b/sqrt(3)*cos(arccos(cos_3eta)/3 + 2*pi/3)
		          psi_max, lag_max= Schwarzschild(closest_approach,pi/2.)
		          psi_min, lag_min= Schwarzschild(R,pi-alpha)
		          psi=2*psi_max - psi_min    
		          lag=2*lag_max - lag_min # + 2*(R - closest_approach + R_g*log((R - R_g)/(closest_approach - R_g)))/c 
		          if psi>pi:
		                return pi+eps,lag
		    else:
		          psi=0
		          R_ref=R_e*(1 - flattening*(cos(theta_deg*pi/180.0))**2)
		          #print("R_ref, alpha, R=",R_ref, alpha, R)
		          #exit()
		          lag=(R_ref - R + R_g*log( (R_ref - R_g)/(R - R_g) ) )/c
		          for i in range(NBend):
		                ex=(kx[i]+1)/2
		                q=(2. - ex*ex - u*(1 - ex*ex)**2/(1 - u))*sin(alpha)**2
		                sr=sqrt(cos(alpha)**2+ex*ex*q)
		                if  2*alpha>pi-eps:
		                      dpsi=b/R/sqrt(q)*wx[i] #*2/2
		                else:
		                      dpsi=ex*b/R/sr*wx[i] #*2/2
		                dlag=dpsi*b/c/(1+sr) #*2/2
		                psi+= dpsi
		                lag+= dlag
		    return psi,lag
	      # flattening=0

		NRadius=2 + int(flattening*R_e/1e-1)
		#print(NRadius)
		NRadius=4+ int(flattening*R_e/1e-1)#4#
		r, dr = linspace(R_e*(1 - flattening),R_e,num=NRadius,retstep=True)
		print(r)
		alpha, dalpha = linspace(0,arccos(-1/sqrt(2*r[0]/R_g/3)),NAlpha,retstep=True)
		psi=zeros((NRadius,NAlpha))
		dt=zeros((NRadius,NAlpha))
		for d in range(NRadius):
			#print(d)
			for a in range(NAlpha):
				psi[d,a],dt[d,a]=Schwarzschild(r[d],alpha[a])

		#s,gn=pr(s,gn,'psidt')      


		Flux=zeros((NPhase,NEnergy,3))
		Flux_obs=zeros((NPhi,NEnergy,3))

		i=pi/180.0*incl_deg#pi*7/18    # line of sight colatitude

		intgspot_version = "new"
		#Even newer version for integration over spot from Vlad:
		if(intgspot_version=="new"):
			#i=pi*7/18    # line of sight colatitude

			#rho_total=0.15 # pi/5 #pi*5/180 # radius of the spot
			rho_total=rho_deg*pi/180.0#pi/5 #pi*5/180 # radius of the spot
			theta_center=pi/180.0*theta_deg#pi/4.1

			antipodal = True

			NSpots=0
			varphi,dvarphi=IntVarphi[0]*pi,IntVarphi[1] *pi
			rho,drho=(IntRho[0]+1)*rho_total/2,(IntRho[1])*rho_total/2
			#       print(IntVarphi)
			#       print(rho,drho)

			l=[]
			theta=[]
			dS=[]

			for v in range(NVarphi):
				for rh in range(NRho):
					NSpots+=1   
					cos_theta=cos(theta_center)*cos(rho[rh])+sin(theta_center)*sin(rho[rh])*cos(varphi[v])
					sin_l=sin(rho[rh])*sin(varphi[v])/sqrt(1- cos_theta**2)
					cos_l=sqrt(1- sin_l**2)
					if cos_theta*cos(theta_center)> cos(rho[rh]) : 
						cos_l=-cos_l      
					l.append(arctan2(-sin_l,-cos_l) + pi)
					theta.append(arccos(cos_theta))  
					dS.append(drho[rh]*dvarphi[v]*sin(rho[rh]))
					#print(v,rh,cos_theta,cos_l,sin_l,l[-1],theta[-1],dS[-1])
					if antipodal :
						NSpots+=1
						l.append(arctan2(sin_l,cos_l) + pi)
						theta.append(pi- theta[-1])
						dS.append(dS[-1])
					#print(v,rh,cos_theta,cos_l,sin_l,l[-1],theta[-1],dS[-1])

			#print(NSpots)


		#New version for integration over spot from Vlad:
		if(intgspot_version=="old"):
			rho=rho_deg*pi/180.0#pi/5 #pi*5/180 # radius of the spot
			Nrho = 10#12#4#15#4#10#4
			drho=rho/Nrho
			antipodal = True#False#True
			l=[0]
			theta=[pi/180.0*theta_deg]#[pi/4.1]
			NSpots=1
			if antipodal :
				l.append(pi)
				theta.append(pi- theta[0])
				NSpots+=1
			for tr in range(Nrho):
				for tph in range(8*tr):
					varphi=(2*tph+1)*pi/8/tr
					cos_theta=cos(theta[0])*cos(drho*tr)+sin(theta[0])*sin(drho*tr)*cos(varphi)
					sin_l=sin(drho*tr)*sin(varphi)/sqrt(1- cos_theta**2)
					cos_l=sqrt(1- sin_l**2)
					if cos_theta*cos(theta[0])> cos(drho*tr) : 
						cos_l=-cos_l
						# print(tr,(2*tph+1)*pi/8/tr,l[-1])
					l.append(arctan2(-sin_l,-cos_l) + pi)
					theta.append(arccos(cos_theta))  
					NSpots+=1
					if antipodal :
						l.append(arctan2(sin_l,cos_l) + pi)
						theta.append(pi- theta[-1])
						NSpots+=1
			###############################

			#NSpots= 2*4 # * somewhat
			## theta = [pi/3,2*pi/3] # spot colatitude
			## l=[0,pi] # spot longitude
			## dS=[1,1] # some arbitrary units
			#little=10.0*pi/180.0#1e-8#1e-2 #TS: made this smaller to match better infitesimal spot
			dS=[1]*NSpots
			#l=[0,0,little,little,pi,pi,pi+little,pi+little]
			## pi*=2/3
			#theta=[pi/3,pi/3+little,pi/3,pi/3+little,2*pi/3,2*pi/3+little,2*pi/3,2*pi/3+little]
			## pi*=3/2

			###############################
			#i=pi/180.0*incl_deg#40.0
			#thettta = pi/180.0*theta_deg#60.0
			#theta=[thettta,thettta+little,thettta,thettta+little,2*thettta,2*thettta+little,2*thettta,2*thettta+little]

			#In case of 1 spot:
			####i=pi*4/18#1.5    # line of sight colatitude #in radians
			###i=pi/180.0*40.0
			#NSpots= 1 #2 # * somewhat
			###theta = [pi/3,pi-pi/3] # spot colatitude
			#theta = [thettta,pi-thettta] # spot colatitude
			#l=[0,pi] # spot longitude
			#dS=[1,1] # some arbitrary units
			###############################


		outParams.write(str(i)+' = sight colatitude i\n')
		outParams.write(str(theta)+' = spot colatitudes theta\n')
		outParams.write(str(l)+' = spot longitudes l\n')

		sin_i=sin(i)
		cos_i=cos(i)

		BoloFlux=zeros((NPhase,3))
		z=cos(linspace(-pi/2,pi/2,num=NZenithBig))
		logIntensity=zeros((NEnergy,NZenithBig,3))


		#Setting the phase delay 0 for photons emitted at phi=0 for 1 small spot):
		cos_psi_ph0=cos_i*cos(theta_center) + sin_i*sin(theta_center)
		psi_ph0 = arccos(cos_psi_ph0)

		R=R_e*(1 - flattening*(cos(theta_center)**2)) 
		r1=bisect(r[1:-1],R) 
		r2=r1 + 1
		dr1=(R - r[r1])/dr
		dr2=(r[r2] - R)/dr

		a1=bisect(psi[r1], psi_ph0)
		a2=bisect(psi[r2], psi_ph0)
		if(a1 >= len(psi[r1])):
			a1=len(psi[r1])-1
		if(a2 >= len(psi[r2])):
			a2=len(psi[r2])-1
		psi1=psi[r1,a1]					
		psi2=psi[r2, a2]
		dpsi1=psi1 - psi[r1, a1 - 1]
		dpsi2=psi2 - psi[r2, a2 - 1]
		dpsi=dpsi1*dr2 + dpsi2*dr1
		dalpha1 = dalpha*(psi1 - psi_ph0)/dpsi1
		dalpha2 = dalpha*(psi2 - psi_ph0)/dpsi2
		alpha1=alpha[a1] - dalpha1
		alpha2=alpha[a2] - dalpha2
		dt1=dt[r1,a1 - 1]*dalpha1/dalpha + dt[r1,a1]*(1. - dalpha1/dalpha)
		dt2=dt[r2,a2 - 1]*dalpha2/dalpha + dt[r2,a2]*(1. - dalpha2/dalpha)
		dphase_ph0=(dt1*dr2 + dt2*dr1)*nu


		for e in range(NEnergy):
			IntInt=CubicSpline(mu,Intensity[e,:,0],extrapolate=True) # interpolate intensity
			#IntInt = interp1d(mu,Intensity[e,:,0])#,fill_value="extrapolate")
			IQ=CubicSpline(mu,Intensity[e,:,1],extrapolate=True) #
			#IQ=interp1d(mu,Intensity[e,:,1])#,fill_value="extrapolate") #


   

			for d in range(NZenithBig):

				#TS: Prevent the code from computing log(0):
				#log0 = -100.0
				#if(IntInt(z[d])< 1e-8 or absolute(IQ(z[d])) < 1e-8):
				#	if(absolute(IQ(z[d])) < 1e-8):
				#		if(IntInt(z[d]) < 1e-8):
				#			logIntensity[e,d] = log0,log0,sign(IQ(z[d]))
				#		else: 
				#			logIntensity[e,d] = log(max(0,IntInt(z[d]))),log0,sign(IQ(z[d]))
				#	else:
				#		if(IntInt(z[d]) < 1e-8):
				#			logIntensity[e,d] = log0,log(absolute(IQ(z[d]))),sign(IQ(z[d]))
				#else:
				#	logIntensity[e,d] = log(max(0,IntInt(z[d]))),log(absolute(IQ(z[d]))),sign(IQ(z[d]))
				logIntensity[e,d] = log(max(0,IntInt(z[d]))),log(absolute(IQ(z[d]))),sign(IQ(z[d]))

		mu=z.copy()   

   
		#s,gn=pr(s,gn,'second precomp')
		for p in range(NSpots):
			sin_theta=sin(theta[p])
			cos_theta=cos(theta[p])

			R=R_e*(1 - flattening*cos_theta**2) 
			dR=2*R_e*flattening*cos_theta*sin_theta # dR / d\theta

			r1=bisect(r[1:-1],R) 
			r2=r1 + 1
			dr1=(R - r[r1])/dr
			dr2=(r[r2] - R)/dr

			redshift=1.0/sqrt(1.0 - R_g/R) # 1/sqrt(1-R_g/R) = 1+ z = redshift
			f=redshift/R*dR
			sin_gamma=f/sqrt(1 + f**2) # angle gamma is positive towards the north pole 
			cos_gamma=1.0/sqrt(1 + f**2)
			beta=2*pi*nu*R*redshift*sin_theta/c
			Gamma=1.0/sqrt(1.0 - beta**2)
			Gamma1= (1.0-sqrt(1.0 - beta**2) )/ beta
		          
			for t in range(NPhi):
				if True: # find mu
					phi0=phi[t]+l[p]
					sin_phi=sin(phi0)
					cos_phi=cos(phi0)
					cos_psi=cos_i*cos_theta + sin_i*sin_theta*cos_phi
					sin_psi=sqrt(1. - cos_psi**2)


					psi0=arccos(cos_psi) 
					a1=bisect(psi[r1], psi0)
					a2=bisect(psi[r2], psi0)
					#print(psi0)
					#print(len(psi[r1]),a1,a2)
					#TS: bug correction ....#####
					if(a1 >= len(psi[r1])):
						a1=len(psi[r1])-1
					if(a2 >= len(psi[r2])):
						a2=len(psi[r2])-1
					########
					psi1=psi[r1,a1]					
					psi2=psi[r2, a2]
					dpsi1=psi1 - psi[r1, a1 - 1]
					dpsi2=psi2 - psi[r2, a2 - 1]
					dpsi=dpsi1*dr2 + dpsi2*dr1
					dalpha1 = dalpha*(psi1 - psi0)/dpsi1
					dalpha2 = dalpha*(psi2 - psi0)/dpsi2
					alpha1=alpha[a1] - dalpha1
					alpha2=alpha[a2] - dalpha2

					cos_alpha = cos(alpha2*dr1 + alpha1*dr2) # linear interpolation of alpha(psi)
					sin_alpha = sqrt(1. - cos_alpha**2)
					sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift
					dcos_alpha=sin_alpha_over_sin_psi *dalpha/dpsi # d cos\alpha \over d \cos \psi

					# cos_alpha, dcos_alpha=Beloborodov(cos_psi) # insert exact formula here
					# sin_alpha = sqrt(1. - cos_alpha**2)
					# sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift

					dt1=dt[r1,a1 - 1]*dalpha1/dalpha + dt[r1,a1]*(1. - dalpha1/dalpha)
					dt2=dt[r2,a2 - 1]*dalpha2/dalpha + dt[r2,a2]*(1. - dalpha2/dalpha)

					dphase=(dt1*dr2 + dt2*dr1)*nu # \delta\phi = \phi_{obs} - \phi
					#print(dphase_ph0,dphase,t,theta[p]*180.0/pi)#dphase = 0
					#exit()
					phase_obs[t]=( phi[t]/2/pi+dphase-dphase_ph0)%1  #( phi[t]/2/pi+dphase)%1.
					#print(dphase-dphase_ph0,t,theta[p]*180.0/pi)
                                        
					cos_xi = - sin_alpha_over_sin_psi*sin_i*sin_phi
					delta = 1./Gamma/(1.-beta*cos_xi)
					#print("delta=",delta)
					#exit()
					cos_sigma = cos_gamma*cos_alpha + sin_alpha_over_sin_psi*sin_gamma*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)

					sin_sigma = sqrt(1. - cos_sigma**2)
					mu0=delta*cos_sigma # cos(sigma')
					Omega=dS[p]*mu0*redshift**2*dcos_alpha*Gamma*R*R/cos_gamma   
					# Omegaarray[t]=max(Omega,0)
					#print(t,' : \t',mu0,' \t ',dcos_alpha,'\t',cos_alpha,cos_psi,Omega)
					if mu0<0: # this only for speeding up. the backwards intensity is usually zero
						Flux_obs[t]=0
						continue 


					if True: # find chi
						sin_chi_0= - sin_theta*sin_phi # times sin psi
						cos_chi_0=sin_i*cos_theta - sin_theta*cos_i*cos_phi # times sin psi 
						chi_0=arctan2(sin_chi_0,cos_chi_0)

						sin_chi_1=sin_gamma*sin_i*sin_phi*sin_alpha_over_sin_psi #times sin alpha sin sigma 
						cos_chi_1=cos_gamma - cos_alpha*cos_sigma  #times sin alpha sin sigma 
						chi_1=arctan2(sin_chi_1,cos_chi_1)

						sin_lambda=sin_theta*cos_gamma - sin_gamma*cos_theta
						cos_lambda=cos_theta*cos_gamma + sin_theta*sin_gamma
						cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos_phi + cos_psi*sin_gamma) - cos_alpha*sin_gamma
						# this line is the longest one
						# alt_cos_eps=(cos_sigma*cos_gamma - cos_alpha)/sin_gamma # legit! thanks God I checked it!
						# sin_chi_prime=cos_eps*mu0*Gamma*beta # times something
						#sin_chi_prime=cos_eps*mu0*delta*Gamma*beta*(1-Gamma1*cos_xi)# times something
						# cos_chi_prime=1. - cos_sigma**2 /(1. - beta*cos_xi) # times the samething
						#cos_chi_prime=sin_sigma**2 - Gamma*mu0**2*beta*cos_xi*(1 - Gamma1*cos_xi)  # times the samething
							
						#OR THE ORIGINAL VERSION:
						sin_chi_prime=cos_eps*mu0*Gamma*beta#*the_thing#*delta**3#*delta
						cos_chi_prime=(1. - cos_sigma**2 /(1. - beta*cos_xi))#*the_thing

						chi_prime=arctan2(sin_chi_prime,cos_chi_prime)
                                                
						#cos_eps_sph = sin_alpha_over_sin_psi*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)
						#sin_chi_prime_sph=cos_eps_sph*mu0*delta*Gamma*beta*(1-Gamma1*cos_xi)# times something
						#cos_chi_prime_sph=sin_alpha**2 - Gamma*mu0**2*beta*cos_xi*(1 - Gamma1*cos_xi)  # times the samething
						#chi_prime_sph= arctan2(sin_chi_prime_sph,cos_chi_prime_sph)


						chi=chi_0+chi_1+chi_prime #oblate formula
						#chi=chi_0+chi_prime_sph #spherical formula

						#  print(chi,'\t',chi_0/chi,'\t',chi_1/chi ,'\t',  chi_prime/chi )

						sin_2chi=sin(2*chi)
						cos_2chi=cos(2*chi)
						# print(chi_prime,' \t',cos_chi_prime )


				d2=bisect(mu[:-1],mu0)
				d1=d2-1
				# print(mu0,mu[d2],mu[d1],d2,d1,'       ')
				mu1,mu2=mu[d1],mu[d2]
				dmu, dmu1, dmu2 = mu2 - mu1, mu0 - mu1, mu2 - mu0
				shift=delta/redshift


				for e in (118,119):# compute only for one observed energy (4.94 keV):
				#for e in range(NEnergy): 
					x0=x[e]/shift
					e1=bisect(x[1:-1],x0) # not the fastest way? anybody cares? ## seems, that light bending is more time consuming anyways
					e2=e1+1

					x1, x2 = x[e1], x[e2]
					dx, dx1, dx2 = x2 - x1, x0 - x1, x2 - x0
					logIQ = (
					      dx2*dmu2*logIntensity[e1, d1] + 
					      dx2*dmu1*logIntensity[e1, d2] +
					      dx1*dmu2*logIntensity[e2, d1] +
					      dx1*dmu1*logIntensity[e2, d2] # bilinear interpolation of the Stokes parameters
					)/dx/dmu 
					I,Q=exp(logIQ[:2])* shift**3 * Omega
					Q*=logIQ[2]

					if I<0: ############
						print('never')
					Flux_obs[t,e]=[I, Q*cos_2chi, Q*sin_2chi]
					#Flux_obs[t,e]=[I, -I*0.1171*cos_2chi*(1.0 - mu0)/(1.0 + 3.82*abs(mu0)),-0.1171*I*sin_2chi*(1.0 - mu0)/(1.0 + 3.82*abs(mu0))]
					#Flux_obs[t,e]=[I, (1.0 - cos_alpha),(1.0 - cos_alpha)]
				#print(t,mu0,-0.1171*(1.0 - mu0)/(1.0 + 3.82*abs(mu0)))
					#Flux_obs[t,e]=[I, cos_2chi, sin_2chi]

			for t in range(NPhase):
				phase0=phase[t]
				for t2 in range(NPhi):
					t1=t2-1
					phase2=phase_obs[t2]
					phase1=phase_obs[t1]
					if phase0>phase1 and phase0<phase2:
						break
				else :
					#phase0=0
					if(phase0>max(phase_obs)):
						phase0=phase0-1
					t2=argmin(phase_obs)
					t1=t2-1
					phase2=phase_obs[t2]
					phase1=phase_obs[t1]-1
					#print("else happened: ",t1,t2,phase1,phase2,phase0)
					#print(phase_obs)

				dphase1=phase0-phase1
				dphase2=phase2-phase0
				dphase=phase2-phase1
				Flux[t-1]+=(Flux_obs[t2]*dphase1+Flux_obs[t1]*dphase2)/dphase #TS: changed here t to t-1 to match the phaseshift to my fortran results
				#print(Flux[t-1,0], t-1, NPhi)

		#s,gn=pr(s,gn,'curves done ')
	    
	# In[31]:
	savePulse = False#True
	if savePulse:
		outF = open(PulsName + 'FF.bin','w')
		outf = open(PulsName + 'ff.bin','w')
		Flux.tofile(outF,format="%e")
		phase.tofile(outf,format="%e")


	#for e in range(0,NEnergy,20): # too many pictures
	#	I=zeros(NPhase)
	#	Q=zeros(NPhase)
	#	U=zeros(NPhase)
	#	for t in range(NPhase):
	#		I[t],Q[t],U[t]=Flux[t,e]
	#	p=sqrt(Q**2+U**2)/I*100
	#	PA=arctan2(-U,-Q)*90/pi+90
	#	print(Q[0], U[0], Flux[0,0], x[e]*evere/1e3)

	#print('end')
	return Flux

#compf(1.0,1.0,1.0,1.0)


