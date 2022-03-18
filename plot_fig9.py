#import:
import matplotlib
matplotlib.use('Agg')
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil, ma
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
#i=pi*4/18#50.0*pi/180.0#pi*4/18#1.5    # line of sight colatitude #in radians
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
incl = 40.0*pi/180.0

#shapes = ["Sphere","AlGendy"]

#colors = ["yellow","black","red"]
colors = ["darkorange","blue","green"]
#colors = ["blue"]
shapes = np.copy(colors)

#figA=figure(2,figsize=(10,14))
#figA.clear()


#rc("text", usetex=True)
figA = figure(figsize=(14,18), dpi=300) #8,6                                                                                                  
#figA = figure(figsize=(15,16), dpi=300) #8,6
#rc("font", family="serif")
#rc("font",serif="Times")
matplotlib.pyplot.figure(1)
lbfontsz = 30#35#30#25 
lwidth= 2.5#3.0#2.5#2.0#1.5 
lpad = 10 #15 #12 #10 #12 

rc("xtick", labelsize=lbfontsz)
rc("ytick", labelsize=lbfontsz)
rc("axes", linewidth=lwidth)
#figA.clear()
matplotlib.pyplot.rcParams.update({'axes.titlesize': lbfontsz})
matplotlib.pyplot.rcParams.update({'font.size': lbfontsz})
matplotlib.pyplot.rcParams.update({'lines.linewidth': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'xtick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'xtick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'font.family': 'serif'})
#matplotlib.pyplot.rcParams.update({'font.serif': 'Times'})

plot_only_I = False#True
plot_QU = False#True

def maskchi(c):
	cn=np.zeros((len(c)))
	for i in range(len(c)):
		if abs(c[i]-c[i-1])>90:
			cn[i] = ma.masked
		else:
			cn[i] = c[i]
	return cn



if(plot_only_I):
	figA=figure(1,figsize=(14,14))
	figA.clear()
	plotAF=figA.add_subplot(1,1,1,yscale='linear') 
else:	
	matplotlib.pyplot.subplots_adjust(wspace=0, hspace=0)
	plotAF=figA.add_subplot(3,1,1,yscale='linear') 
	plotAp=figA.add_subplot(3,1,2)      #
	plotAc=figA.add_subplot(3,1,3)      #
	#plotAd=figA.add_subplot(4,1,4)      #

for ish in range(1,3):#len(shapes)):

	#oblateness='AlGendy'#'Sphere'#'AlGendy'
	oblateness=shapes[ish]
	print(ish)
	#AtmName='res/B/C1obl' # the prefix for all result files related to the set of parameters
	#AtmName='res/B/B0P2' # the prefix for all result files related to the set of parameters
	if(ish == 0):
		#AtmName='res/C2/C1sph'
		#PulsName='res/C2/sph2test'
		#PulsName='res/B/lcpol_sph_r2'
		#PulsName='res/B/lcpol_sph_sp2_11km'
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lcpol_obl_sp1_10msun'
		#PulsName='res/B/lbb_rhoinf'
		PulsName='res/B/lbb_rho10_sp2_f600_obl_burst2_dt2_m1'#'res/B/lbb_rhoinf_chi-1'
	if(ish == 1):
		#AtmName='res/C2/C1sph'
		#PulsName='res/C2/sph2test'
		#PulsName='res/B/lcpol_sph_sp1'
		#PulsName='res/B/lcpol_obl_sp1_15msun'
		#PulsName='res/B/lbb_rhoinf_sp1'
		#PulsName='res/B/lbb_rhoinf_sp1_f001'
		#PulsName='res/B/lbb_rhoinf_sp1_f100_p100_ires'
		#PulsName='res/B/lbb_rhoinf_sp1_f001_p100_ires'#THIS one was shown to work with arcmancer
		#PulsName='res/B/lbb_rho10_sp1_f600_obl'
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/B0Prho10'
		#PulsName='res/B/B0Prho10oblsp2' #in the old version
		PulsName='pOS_pulses/lbb_rho10_sp2_f600_obl_burst2_dt2'
	if(ish == 2):
		#AtmName='res/C2/C1obl'  
		#PulsName='res/C2/obl2test'
		#PulsName='res/B/lcpol_sph_r3'
		#PulsName='res/B/lcpol_obl_sp1_20msun'
		#PulsName='res/B/lbb3_rho10'
		#PulsName='res/B/lbb3_polburst_sp1_rho10_f100_'
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lbb_rhoinf_sp1_f001_p100_ires_swit'
		#PulsName='res/B/lbb_rho10_sp1_f600_sph'
		#PulsName='res/B/B0Prho10sphsp2'
		#PulsName='res/B/B0Prho10sphsp2_dmr' #this used in the old version
		PulsName='pOS_pulses/lbb_rho10_sp2_f600_sph_burst2_dt2'#dmr2_burst'
		#PulsName='res/B/lbb_rho10_sp2_f600_obl_burst2_dt2_m2'
                
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
	
	labelsize=lbfontsz#20
	fontsize=lbfontsz#25

	if(ish == 2):
		NPhase = 150#500
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

	#if(ish > 0):
	#	plotAF.plot(phase,I/I.max())#,color=col)
	#	plotAp.plot(phase,p)#,color=col)
	#	plotAc.plot(phase,PA)#,color=col)	 
	#	break	           	

	#figA.suptitle(r'$\nu={:5.0f}Hz$'.format(nu)+
	#              r'$,\,R_e={:5.1f}km$'.format(R_e)+
	#              #r'$,\,R_e=11,12,14km$'.format(R_e)+
	#              #r'$,\,M=1.0, 1.5, 2.0$'+r'$M_{\odot}$'+',\n'+
	#              r'$,\,M=$'+str(M)+r'$M_{\odot}$'+',\n'+
	#              r'$\,\theta={:5.1f}\degree$'.format(theta[0]*180/pi)+#r'$,{:5.1f}\degree$'.format(40.0)+
	#              r'$,\,i={:5.1f}\degree$'.format(incl*180/pi)+',\n'#r'$,{:5.1f}\degree$'.format(60.0)+',\n'
	#              r'$\rho={:5.1f}\degree$'.format(rho)+', '+
	#              r'$\,E={:6.2f}keV$'.format(x[ene]*evere/1e3),fontsize=fontsize)   

	# plotAF.set_xlim([x[0],x[-1]])
	plotAF.set_xlim(0,1)

	# plotAF.locator_params(axis='y', nbins=10)
	plotAF.set_ylabel(r"$F_{\mathrm{I}}/F_{\mathrm{I}}^{\mathrm{max}}$",fontsize=fontsize)
	plotAF.tick_params(axis='both', which='major', labelsize=labelsize,direction='in',pad=lpad,top=True,right = True)
	# plotAF.plot(x,xIinx,'k-.')

	# plotAF.set_xlim([x[0],x[-1]])
	if not(plot_only_I):
		if(plot_QU):
			plotAp.set_xlim(0,1)
			plotAp.tick_params(axis='both', which='major', labelsize=labelsize,direction='in',pad=lpad,top=True,right = True)
			plotAp.set_ylabel(r'$Q/Q_{max}$',fontsize=fontsize)
			plotAc.tick_params(axis='both', which='major', labelsize=labelsize,direction='in',pad=lpad,top=True,right = True)
			plotAc.set_ylabel(r'$U/U_{max}$',fontsize=fontsize)
			plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)		
		else:
			plotAp.set_xlim(0,1)
			plotAp.tick_params(axis='both', which='major', labelsize=labelsize,direction='in',pad=lpad,top=True,right = True)
			plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
			#plotAp.set_ylim(99.8,100.1)
			plotAc.set_xlim(0,1)
			plotAc.set_ylim(0,180)
			#plotAc.set_ylim(40,140)
			#plotAc.set_yticks([0,30,60,90,120,150,180])
			#plotAc.set_ylim(-180,180)
			#plotAc.set_yticks([0,30,60,90,120,150,180])
			plotAc.tick_params(axis='both', which='major', labelsize=labelsize,direction='in',pad=lpad,top=True,right = True)
			plotAc.set_ylabel(r'$\chi\,[\mathrm{deg}]$',fontsize=fontsize)
			plotAc.set_xlabel(r'$\varphi/(2\pi)$',fontsize=fontsize)

	#col=colors[(e*NColors)//NEnergy]
	col = colors[ish]

	plotAF.plot(phase,I/I.max(),color=col,linewidth=lwidth)
	#plotAF.plot(phase,I,color=col)
	PA_VP04 = PA
	phase_VP04 = phase
	if(ish == 1):
		PA0_VP04 = PA#I/I.max()#PA
		F0_VP04 = I/I.max()
		phase0_VP04 = phase

	if not(plot_only_I):
		if(plot_QU):
			plotAp.plot(phase,Q/np.max(abs(Q)),color=col)
			plotAc.plot(phase,U/np.max(abs(U)),color=col)
		else:
			#plotAp.plot(phase,p,color=col,marker="o",markersize=1.0)
			#plotAc.plot(phase,PA,color=col,marker="o",markersize=1.0)
			plotAp.plot(phase,p,color=col,linewidth=lwidth)
			plotAc.plot(phase,maskchi(PA),color=col,linewidth=lwidth)

compare_to_tslc = False#True
if(compare_to_tslc):
	colors = ["cyan", "orange", "red"]#"red"]
	#colors = ["red"]
	for ic in range(2,3):#3):

		datafile=""
		if(ic == 0):
			#datafile = "../../MCMC/results_sph_test.txt"
			#datafile = "../../MCMC/results_sph_11kmsp2.txt"
			datafile = "../../MCMC/results_obl_rhoinf.txt"
			#datafile = "../../MCMC/results_obl_rho10_sp1.txt"
		if(ic == 1):
			#datafile = "../../MCMC/results_obl2_test.txt"
			#datafile = "../../MCMC/results_obl2_simpl.txt"
			#datafile = "../../MCMC/results_sph_12kmsp2.txt"
			datafile = "../../MCMC/results_obl_rho1.txt"
		if(ic == 2):
			#datafile = "../../MCMC/results_sph_14kmsp2.txt"
			datafile = "../../MCMC/results_obl_rho10.txt"
			#datafile = "../../MCMC/results_obl_rhoinf.txt"
		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = 256#128#16
		#c_lines = 21#28
		c_lines = 29
		egrid = 40
		full_chain= [[] for x in xrange(egrid+1)]

		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		for j in range(0,len(full_chain)):
			for i in range(c_lines,c_lines+Nchain_size): #not reading comment lines
				parts = lines[i].split()
				full_chain[j].append(float(parts[j]))
			parts = lines[c_lines].split()

		full_chain = np.array(full_chain)

		#print full_chain[0,:]
		#print full_chain[30,:]
		energy_keV = np.arange(2.0,25.00001, (25.0-2.0)/int(egrid-1))
		phase = full_chain[0,:]
		norm_obsF = np.zeros((len(phase), egrid))
		for i in range(1,egrid+1):
			norm_obsF[:,i-1] = full_chain[i,:]*energy_keV[i-1]/np.max(full_chain[i,:]*energy_keV[i-1])


		#for i in range(0,egrid):
		#print energy_keV
		ene = 5#0#5#14#38##14
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.0#-0.005#0.025
		for i in range(ene,ene+1):
			plotAF.plot(phase-phshift,norm_obsF[:,i],markersize=2,label=str("%.2f" % energy_keV[i])+' keV',color=colors[ic])
		#or plot bolometric energy flux:
		#bolomflux = full_chain[len(full_chain)-1,:]/np.max(full_chain[len(full_chain)-1,:])
		##print bolomflux
		#plotAF.plot(phase-phshift,bolomflux,markersize=2,label=str("bolomflux"),color=colors[ic])

compare_to_fortran_vlad = False#True
if(compare_to_fortran_vlad): 
	#colors = ["green","blue","black"]
	colors = ["orange"]
	for ic in range(0,1):
		#datafile = "../../MCMC/results_obl_vlad_rhotest.txt"
		datafile = "../../MCMC/results_obl_vlad_rho10f600.txt"
		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = sum(1 for line in open(datafile))
		c_lines = 29
		egrid = 5#3
		full_chain= [[] for x in xrange(egrid+1)]

		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		for j in range(0,len(full_chain)):
			for i in range(c_lines,Nchain_size): #not reading comment lines
				parts = lines[i].split()
				full_chain[j].append(float(parts[j]))
			parts = lines[c_lines].split()

		full_chain = np.array(full_chain)

		energy_keV = [4.94]
		phase = full_chain[0,:]
		norm_obsF = np.zeros((len(phase), egrid))
		for i in range(1,egrid+1):
			#norm_obsF[:,i-1] = full_chain[i,:]*energy_keV[i-1]/np.max(full_chain[i,:]*energy_keV[i-1])
			norm_obsF[:,i-1] = full_chain[i,:]/np.max(full_chain[i,:])


		ene = 0
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.0#-0.01#-0.02#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(phase,phshift)
		#print phase_new
		for i in range(ene,ene+1):
		#for i in range(0,egrid):
			for ipha in range(0,len(phase_new)-1):
				if(phase_new[ipha+1] > phase_new[ipha]):
					plotAF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],"-",markersize=5,label=str("%.2f" % energy_keV[i])+' keV',color=colors[i])

		#phase_new = shift_phase(phase,-phshift)
		col = colors[ic]
		if not(plot_only_I):
			if(plot_QU):
				Q = -full_chain[2,:]
				U = -full_chain[3,:]
				#phase_new = shift_phase(phase,0.372)
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/np.max(abs(Q)),color=col)
						plotAc.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/np.max(abs(U)),color=col)
			else:
				#p = full_chain[3,:]*100.0
				#PA = full_chain[2,:]*180.0/pi+90.0
				p = full_chain[5,:]
				PA = arctan2(full_chain[3,:],full_chain[2,:])*90/pi+90#full_chain[4,:]*180.0/pi+90.0#arctan2(full_chain[3,:],full_chain[2,:])*90/pi+90#full_chain[4,:]*180.0/pi+90.0
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAp.plot(phase_new[ipha:ipha+2],p[ipha:ipha+2],color=col)
						plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],color=col)



compare_to_arcmancer = False#True
if(compare_to_arcmancer): 
	#colors = ["green","blue","black"]
	colors = ["blue","blue"]
	for ic in range(0,1):
		if(ic == 0):
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x01_sph.csv"# (copy).csv"
			datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d60_i40_x10_obl.csv"# (copy).csv"
			#datafile = "../arcmancer/out3/pol_angle_test.csv"# (copy).csv"
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x01.csv"# (copy).csv"
		if(ic == 1):
			datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d40_i60_x01_sph.csv"
			#datafile = "../arcmancer/out3/polar_f700_bb_r12_m1.6_d50_i50_x05.csv"
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x01.csv"# (copy).csv"
		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = sum(1 for line in open(datafile))
		#c_lines = 21#28
		c_lines = 1
		#egrid = 3
		egrid = 6#5#3
		full_chain= [[] for x in xrange(egrid+1)]

		input = file(datafile, 'r')
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
		#print phase_new
		for i in range(ene,ene+1):
		#for i in range(0,egrid):
			for ipha in range(0,len(phase_new)-1):
				if(phase_new[ipha+1] > phase_new[ipha]):
					plotAF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],label=str("%.2f" % energy_keV[i])+' keV',color=colors[i],marker="o",markersize="1.0")

		#phase_new = shift_phase(phase,-phshift)
		col = colors[ic]
		if not(plot_only_I):
			if(plot_QU):
				Q = full_chain[2,:]
				U = full_chain[3,:]
				#phase_new = shift_phase(phase,0.372)
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/np.max(abs(Q)),color=col)
						plotAc.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/np.max(abs(U)),color=col)
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
						plotAp.plot(phase_new[ipha:ipha+2],p[ipha:ipha+2],"-o",color=col,markersize="1.0")
						plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-o",color=col,markersize="1.0")
				if(ic == 0):
					PA_acm0 = PA#norm_obsF[:,ene]#PA
					F_acm0 = norm_obsF[:,ene]
					phase_acm0 = phase_new

				if(ic == 1):
					PA_acm = PA#norm_obsF[:,ene]#PA
					phase_acm = phase_new


plot_PA_residuals = False#True
if(plot_QU):
	plot_PA_residuals = False
if(plot_PA_residuals): 
	col = "black"
	plotAd.set_xlim(0,1)
	#plotAd.set_ylim(0.97,1.06)
	#plotAd.set_ylim(0.0,0.11)
	#plotAd.set_ylim(-4.0,3.0)
	#plotAd.set_ylim(-0.02,0.02)
	#plotAd.set_ylim(-0.004,0.004)
	#plotAd.set_ylim(-0.01,0.01)
	#plotAd.set_ylim(-0.1,0.1)
	plotAd.set_ylim(0.0,0.5)
	#plotAd.set_yticks([0,30,60,90,120,150,180])
	#plotAd.tick_params(axis='both', which='major', labelsize=labelsize)
	plotAd.tick_params(axis='y', which='major', labelsize=labelsize,direction='in',top=True,right = True)
	plotAd.tick_params(axis='x', which='major', labelsize=labelsize,direction='in',top=True,right = True)
	#plotAd.set_ylabel(r'$\chi_{\mathrm{acm}}/\chi_{\mathrm{vp}}$',fontsize=fontsize)
	#plotAd.set_ylabel(r'$\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}} [\degree]$',fontsize=fontsize)
	#plotAd.set_ylabel(r'$\frac{\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}}}{\chi_{\mathrm{max}}-\chi_{\mathrm{min}}}$',fontsize=fontsize)
	plotAd.set_ylabel(r'$|\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}}|\,[\degree]$',fontsize=fontsize)
	plotAd.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

	#print(len(phase_VP04),len(phase_acm))
	PA_VP04_interp = interp1d(phase_VP04,PA_VP04)#,"cubic")
	PA0_VP04_interp = interp1d(phase0_VP04,PA0_VP04)#,"cubic")
	F0_VP04_interp = interp1d(phase0_VP04,F0_VP04)#,"cubic")
	#PAacm = interp1d(phase_acm,PA_acm)
	#plotAd.plot(phase_VP04,PA_VP04,color="red")
	#norm1 = abs(np.max(PA_acm)-np.min(PA_acm))
	norm2 = abs(np.max(PA_acm0)-np.min(PA_acm0))
	norm3 = abs(np.max(F_acm0)-np.min(F_acm0))
	for ipha in range(0,len(phase_new)-1):
		if(phase_new[ipha+1] > phase_new[ipha]):
			#if(abs(PA_VP04_interp(phase_acm[ipha])-PA_VP04_interp(phase_acm[ipha+1]))/PA_VP04_interp(phase_acm[ipha+1]) < 0.2):
			#plotAd.plot(phase_acm[ipha:ipha+2],abs(PA_acm[ipha:ipha+2]-PA_VP04_interp(phase_acm[ipha:ipha+2]))/PA_VP04_interp(phase_acm[ipha:ipha+2]),color=col)
			#plotAd.plot(phase_acm[ipha:ipha+2],PA_acm[ipha:ipha+2]/PA_VP04_interp(phase_acm[ipha:ipha+2]),color=col)
			#plotAd.plot(phase_acm[ipha:ipha+2],PA_VP04_interp(phase_acm[ipha:ipha+2])/PA_acm[ipha:ipha+2],color=col)
			#plotAd.plot(phase_acm[ipha:ipha+2],(PA_acm[ipha:ipha+2]-PA_VP04_interp(phase_acm[ipha:ipha+2]))/norm1,color="orange")
			#plotAd.plot(phase_acm0[ipha:ipha+2],(PA_acm0[ipha:ipha+2]-PA0_VP04_interp(phase_acm0[ipha:ipha+2]))/norm2,color="blue")
			plotAd.plot(phase_acm0[ipha:ipha+2],abs(PA_acm0[ipha:ipha+2]-PA0_VP04_interp(phase_acm0[ipha:ipha+2])),color="blue")


			#plotAd.plot(phase_acm0[ipha:ipha+2],(F_acm0[ipha:ipha+2]-F0_VP04_interp(phase_acm0[ipha:ipha+2]))/norm3,color="cyan") #error in normal flux

#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))



plot_analytic_estimates = False#True
plot_analytic_estimates2 = False#True
if(plot_analytic_estimates):
	i = 40.0*pi/180.0
	theta = 60.0*pi/180.0
	#print sin(90*pi/180.0)
	phi = (np.array(phase_VP04))*2.0*pi
	tanchi0 = - (sin(theta)*sin(phi))/(sin(i)*cos(theta)-cos(i)*sin(theta)*cos(phi))
	#print(tanchi0)
	chi0 = np.arctan(tanchi0)*180/pi+90
	#print(chi0)
	plotAc.plot(phase_VP04,chi0,color="yellow")
if(plot_analytic_estimates2):
	i = 60.0*pi/180.0
	theta = 40.0*pi/180.0
	#print sin(90*pi/180.0)
	phi = (np.array(phase_VP04))*2.0*pi
	tanchi0 = - (sin(theta)*sin(phi))/(sin(i)*cos(theta)-cos(i)*sin(theta)*cos(phi))
	#print(tanchi0)
	chi0 = np.arctan(tanchi0)*180/pi+90
	#print(chi0)
	plotAc.plot(phase_VP04,chi0,color="yellow")


plotAF.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAp.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
#plotAc.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
#plotAd.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())


#align manually ylabels:
labelx = -0.1 #-0.16
plotAF.yaxis.set_label_coords(labelx, 0.5)
plotAc.yaxis.set_label_coords(labelx, 0.5)
plotAp.yaxis.set_label_coords(labelx, 0.5)

plotAF.margins(x =0,y=0)
plotAc.margins(x =0)
plotAp.margins(x =0,y=0)

#adjust limits and ticks:
plotAF.set_ylim(0.2,1.1)
#plotAF.set_yticks([0.2,0.4,0.6,0.8,1.0])
#plotAF.set_yticklabels(["0.2","0.4","0.6","0.8","1.0"],fontstyle="normal")

plotAp.set_ylim(0.0,4.0)
plotAp.set_yticks([0,1,2,3])
plotAp.set_yticklabels(["0","1","2","3"],fontstyle="normal")


plotAc.set_ylim(10.0,180.0)

#figA.tight_layout()
figA.subplots_adjust(wspace=0, hspace=0)
figA.subplots_adjust(left=0.15)

#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))
figA.savefig('res/B/plot.pdf',bbox_inches='tight')#.format(e))
figA.clf()



