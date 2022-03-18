##To plot polarization position angle profile comparisons between our Obl+Schw and arcmancer

#import:
import matplotlib
matplotlib.use('Agg')
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect

def shift_phase(phi,shift):
	return (phi + shift) % 1 


#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

NPhase = 150#500#150#128*2#120#128#120 # Number of equidistant phase points
NBend= 20 # Number of knots in light bending integrations
NBendPhase= 1000 # Number of psi/aplha grid points
IntBend = leggauss(NBend)
NEnergy = 281 # 50# 101 # number of energy points (x)

phi,phi_weight=linspace(0,2*pi,num=NPhase,endpoint=False,retstep=True) #Size of spacing between samples = phi_weight
nu=600#1#100#600#700#600 # star rotation frequency in Hz
M=1.4#1.6#1.4 # star mass in solar masses
R_g=M*2.95 # gravitational Schwarzschild radius
R_e=12.0 # equatorial radius of the star in kilometers      
i=pi*4/18#50.0*pi/180.0#pi*4/18#1.5    # line of sight colatitude #in radians
#i=pi/180.0*80.0
theta = [pi/3,pi-pi/3]#[50.0*pi/180.0]#[pi/3,pi-pi/3] # spot colatitude
#theta = [pi/180.0*15.0,pi-pi/180.0*15.0] # spot colatitude
rho = 10#1#10

#Try switching:
#theta_new = [i,pi-i] # spot colatitude
#i = theta[0] # spot colatitude
#theta = theta_new


tau_T= .6 # Thomson optical depth of thermalization 
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
Theta = 0.1  # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
x,x_weight=IntEnergy


#shapes = ["Sphere","AlGendy"]

#colors = ["yellow","black","red"]
colors = ["yellow","blue","green"]
#colors = ["blue"]
shapes = np.copy(colors)

figA=figure(2,figsize=(10,14))
figA.clear()

plot_only_I = False#True
plot_QU = False#True

if(plot_only_I):
	#THIS OPTION SHOULD NOT BE USED IN THIS VERSION
	figA=figure(1,figsize=(14,14))
	figA.clear()
else:	
	plotAc=figA.add_subplot(2,1,1)      #
	plotAd=figA.add_subplot(2,1,2)      #

for ish in range(1,3):#len(shapes)):

	#oblateness='AlGendy'#'Sphere'#'AlGendy'
	oblateness=shapes[ish]
	print(ish)
	#AtmName='res/B/C1obl' # the prefix for all result files related to the set of parameters
	#AtmName='res/B/B0P2' # the prefix for all result files related to the set of parameters
	if(ish == 0):
		PulsName='res/B/lbb_rhoinf_chi-1'
	if(ish == 1):
		#PulsName='res/B/lbb_rhoinf_sp1_f001_p100_ires'#THIS one was the first shown to work with arcmancer
		PulsName='res/B/lbb_rho10_sp1_f600_obl'
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lbb_rhoinf_chi0'
	if(ish == 2):
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lbb_rhoinf_sp1_f001_p100_ires_swit'
		PulsName='res/B/lbb_rho10_sp1_f600_sph'
	#PulsName=AtmName+'P1'
	computePulse= True
	plotAtm=not True
	plotPulse=True
	mod=True

	#outF = open(PulsName + 'F.bin','w')
	#outf = open(PulsName + 'f.bin','w')
	#Flux.tofile(outF,format="%e")
	#phi.tofile(outf,format="%e")
	print(PulsName)

	inFlux = open(PulsName+'FF.bin')
	inphi = open(PulsName+'ff.bin')
	Flux1 = fromfile(inFlux)
	phi = fromfile(inphi)
	#print(phi, Flux1)	
	fluxlcurve0 = Flux1[0:len(Flux1):3*NEnergy] #light curve with lowest E
	fluxspec0 = Flux1[0:3*NEnergy:3] #spectrum at phase=0
	#print(fluxlcurve0)
	#print(" ")
	#print(fluxspec0)

	ene = 118#166#140#166#140 #The chosen energy index
	print("The chosen energy (keV): ", x[ene]*evere/1e3)
	fluxlcurve_Iene = Flux1[0+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Qene = Flux1[1+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Uene = Flux1[2+ene*3:len(Flux1):3*NEnergy]
            	            
	#Flux=zeros((NPhase,NEnergy,3))
	#print(fluxlcurve_Iene)
	#print(fluxlcurve_Qene)
	#print(fluxlcurve_Uene)
	
	labelsize=20
	fontsize=25


	#phase=list(phi/2/pi)+[1.]
	phase=list(phi)+[1.]
	I=zeros(NPhase+1)
	Q=zeros(NPhase+1)
	U=zeros(NPhase+1)
	for t in range(NPhase+1):
		#I[t],Q[t],U[t]=Flux[t-1,e]*x[e] 
		I[t],Q[t],U[t]=fluxlcurve_Iene[t-1]*x[ene] ,fluxlcurve_Qene[t-1]*x[ene] ,fluxlcurve_Uene[t-1]*x[ene] 

	p=sqrt(Q**2+U**2)/I*100
	#PA=arctan2(-U,-Q)*90/pi+90
	PA=arctan2(U,Q)*90/pi+90
       	

	figA.suptitle(r'$\nu={:5.0f}Hz$'.format(nu)+
	              r'$,\,R_e={:5.1f}km$'.format(R_e)+
	              #r'$,\,R_e=11,12,14km$'.format(R_e)+
	              #r'$,\,M=1.0, 1.5, 2.0$'+r'$M_{\odot}$'+',\n'+
	              r'$,\,M=$'+str(M)+r'$M_{\odot}$'+',\n'+
	              r'$\,\theta={:5.1f}\degree$'.format(theta[0]*180/pi)+#r'$,{:5.1f}\degree$'.format(40.0)+
	              r'$,\,i={:5.1f}\degree$'.format(i*180/pi)+',\n'#r'$,{:5.1f}\degree$'.format(60.0)+',\n'
	              r'$\rho={:5.1f}\degree$'.format(rho)+',\n'+
	              r'$\,E={:6.2f}keV$'.format(x[ene]*evere/1e3),fontsize=fontsize)  


	if not(plot_only_I):
		if(plot_QU):
			print("plot_QU option not valid in this version!")
			quit()	
		else:

			plotAc.set_xlim(0,1)
			plotAc.set_ylim(0,180)
			#plotAc.set_ylim(40,140)
			#plotAc.set_yticks([0,30,60,90,120,150,180])
			#plotAc.set_ylim(-180,180)
			#plotAc.set_yticks([0,30,60,90,120,150,180])
			plotAc.tick_params(axis='both', which='major', labelsize=labelsize)
			plotAc.set_ylabel(r'$\chi\,[\degree]$',fontsize=fontsize)
			#plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

	#col=colors[(e*NColors)//NEnergy]
	col = colors[ish]

	PA_VP04 = PA
	phase_VP04 = phase
	if(ish == 1): #arcmancer is phase shifted relative to ish=2, so here shift phase of ish=1 to fit best arcmancer (stupid way to do this)
		phshift = 0.0#0.001#0.008#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(np.array(phase),phshift)
		for ipha in range(0,len(phase_new)-1):
			if(phase_new[ipha+1] > phase_new[ipha]):
				plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-o",color=col,markersize="1.0")
		PA0_VP04 = PA
		phase0_VP04 = phase
		phshift0 = phshift
	else:
		#plotAc.plot(phase,PA,color=col,marker="o",markersize=1.0)
		print("...")


compare_to_arcmancer = True
if(compare_to_arcmancer): 
	#colors = ["green","blue","black"]
	colors = ["black","black"]
	for ic in range(0,1):
		if(ic == 0):
			#datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d60_i40_x10_agm.csv"# (copy).csv"
			datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d60_i40_x10_obl.csv"# (copy).csv"
		if(ic == 1):
			datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d40_i60_x01_sph.csv"
			#datafile = "../arcmancer/out3/polar_f700_bb_r12_m1.6_d50_i50_x05.csv"
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x01.csv"# (copy).csv"
		#input = file(datafile, 'r')
		input = open(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = sum(1 for line in open(datafile))
		#c_lines = 21#28
		c_lines = 1
		#egrid = 3
		egrid = 6#5#3
		#full_chain= [[] for x in xrange(egrid+1)]
		full_chain= [[] for x in range(egrid+1)]

		#input = file(datafile, 'r')
		input = open(datafile, 'r')
		lines = input.readlines()
		input.close()

		for j in range(0,len(full_chain)):
			for i in range(c_lines,Nchain_size): #not reading comment lines
				parts = lines[i].split(",")
				#print parts
				full_chain[j].append(float(parts[j]))
			parts = lines[c_lines].split(",")

		full_chain = np.array(full_chain)

		#print full_chain[0,:]
		#print full_chain[30,:]
		#energy_keV = [2.0,6.0,12.0]
		energy_keV = [4.94]
		phase = full_chain[0,:]
		norm_obsF = np.zeros((len(phase), egrid))
		for i in range(1,egrid+1):
			#norm_obsF[:,i-1] = full_chain[i,:]*energy_keV[i-1]/np.max(full_chain[i,:]*energy_keV[i-1])
			norm_obsF[:,i-1] = full_chain[i,:]/np.max(full_chain[i,:])


		#for i in range(0,egrid):
		#print energy_keV
		ene = 0
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.2517#0.25#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(phase,phshift)

		#phase_new = shift_phase(phase,-phshift)
		col = colors[ic]
		if not(plot_only_I):
			if(plot_QU):
				print("plot_QU option not valid in this version!")
				quit()	
			else:
				#p = full_chain[3,:]*100.0
				#PA = full_chain[2,:]*180.0/pi+90.0
				p = full_chain[5,:]*100.0
				PA = full_chain[4,:]*180.0/pi+90.0
				use_PA_avg = False#True
				if use_PA_avg:
					PA = full_chain[6,:]*180/pi+90.0#full_chain[4,:]*180.0/pi+90.0
					#PA = np.array(list(reversed(full_chain[6,:])))*180/pi+90.0#full_chain[4,:]*180.0/pi+90.0
					for qwe in range(0,len(PA)):
						if(PA[qwe] < 0.0):
							PA[qwe] = PA[qwe]+180.0
						if(PA[qwe] > 180.0):
							PA[qwe] = PA[qwe]-180.0
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-o",color=col,markersize="1.0")
				if(ic == 0):
					PA_acm0 = PA#norm_obsF[:,ene]#PA
					F_acm0 = norm_obsF[:,ene]
					phase_acm0 = phase_new

				if(ic == 1):
					PA_acm = PA#norm_obsF[:,ene]#PA
					phase_acm = phase_new


plot_PA_residuals = True
if(plot_QU):
	plot_PA_residuals = False
if(plot_PA_residuals): 
	col = "black"
	plotAd.set_xlim(0,1)
	#plotAd.set_ylim(0.97,1.06)
	#plotAd.set_ylim(0.0,0.11)
	#plotAd.set_ylim(-4.0,3.0)
	#plotAd.set_ylim(-0.02,0.02)
	plotAd.set_ylim(-0.004,0.004)
	#plotAd.set_ylim(-0.01,0.01)
	#plotAd.set_ylim(-0.1,0.1)
	#plotAd.set_yticks([0,30,60,90,120,150,180])
	#plotAd.tick_params(axis='both', which='major', labelsize=labelsize)
	plotAd.tick_params(axis='y', which='major', labelsize=8)
	plotAd.tick_params(axis='x', which='major', labelsize=labelsize)
	#plotAd.set_ylabel(r'$\chi_{\mathrm{acm}}/\chi_{\mathrm{vp}}$',fontsize=fontsize)
	#plotAd.set_ylabel(r'$\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}} [\degree]$',fontsize=fontsize)
	plotAd.set_ylabel(r'$\frac{\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}}}{\chi_{\mathrm{max}}-\chi_{\mathrm{min}}}$',fontsize=fontsize)
	plotAd.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

	#print(len(phase_VP04),len(phase_acm))
	PA_VP04_interp = interp1d(phase_VP04,PA_VP04)#,"cubic")
	#PA0_VP04_interp = interp1d(phase0_VP04,PA0_VP04)#,"cubic")
	print(shift_phase(np.array(phase0_VP04),phshift0),PA0_VP04)
	#quit()
	#PA0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift0),PA0_VP04)#,"cubic")
	PA0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift0),PA0_VP04,fill_value='extrapolate')# workd with newer scipy

	norm = abs(np.max(PA_acm0)-np.min(PA_acm0))
	for ipha in range(0,len(phase_new)-1):
		if(phase_new[ipha+1] > phase_new[ipha]):
			plotAd.plot(phase_acm0[ipha:ipha+2],(PA_acm0[ipha:ipha+2]-PA0_VP04_interp(phase_acm0[ipha:ipha+2]))/norm,color="blue")
			#plotAd.plot(phase_acm0[ipha:ipha+2],(PA_acm0[ipha:ipha+2]-PA_VP04_interp(phase_acm0[ipha:ipha+2]))/norm,color="green")


#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))


#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))
figA.savefig('res/B/plot.pdf')#.format(e))
figA.clf()



