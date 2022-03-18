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
import find_best_phshift

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
incl=pi*4/18#50.0*pi/180.0#pi*4/18#1.5    # line of sight colatitude #in radians
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



##figA=figure(2,figsize=(10,14))
#figA=figure(2,figsize=(14,18))
#figA.clear()


#rc("text", usetex=True)
figA = figure(figsize=(14,30), dpi=300)
#figA = figure(figsize=(14,18), dpi=300) #8,6
#rc("font", family="serif")
#rc("font",serif="Times")
matplotlib.pyplot.figure(1)
lbfontsz = 35 #25 #35 #25 #30#25 
lwidth= 3.0 #2.5 #3.5 #3.0 #2.5 #3.0#2.5#2.0#1.5 
lpad = 10 #12 

labelsize=lbfontsz#40#30#20
fontsize=lbfontsz#35#25
ticksize=lbfontsz#25

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

#matplotlib.pyplot.tight_layout()

plot_only_I = False
plot_all = True
plot_QU = False#True
two_spots = False#True
plot_resid = True

matplotlib.pyplot.subplots_adjust(wspace=0, hspace=0)


if(plot_all):
	if not (plot_resid):
	        plotAF=figA.add_subplot(4,1,1,yscale='linear') 
	        plotAp=figA.add_subplot(4,1,2)      #
	        #plotAc=figA.add_subplot(4,1,3)      #
	        #plotAd=figA.add_subplot(4,1,4)      #
	        plotAd=figA.add_subplot(4,1,3)      #
	        plotAc=figA.add_subplot(4,1,4)      #
	else:
		plotAF=figA.add_subplot(5,1,1,yscale='linear')
		plotAp=figA.add_subplot(5,1,2)
		plotAd=figA.add_subplot(5,1,3)
		plotAc=figA.add_subplot(5,1,4)
		plotAr=figA.add_subplot(5,1,5)
else:	
	plotAc=figA.add_subplot(2,1,1)      #
	plotAd=figA.add_subplot(2,1,2)      #



#shapes = ["Sphere","AlGendy"]

#colors = ["yellow","black","red"]
colors = ["yellow","blue","green"]
#colors = ["blue"]
shapes = np.copy(colors)


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
		#PulsName='res/B/lbb_rho10_sp1_f600_obl'
		#PulsName='res/B/lbb_rho10_sp1_f600_obl_eqnew'
		#PulsName='res/B/lbb_rho10_sp1_f001_obl_accspot'
		#PulsName='res/B/lbb_rho10_sp1_f600_obl'#obl_accspot used in old results
		PulsName='pOS_pulses/lbb_rho10_sp1_f600_obl_burst2_dt'#test'
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lbb_rhoinf_chi0'
	if(ish == 2):
		#PulsName='res/B/B0Ptest'
		#PulsName='res/B/lbb_rhoinf_sp1_f001_p100_ires_swit'
		#PulsName='res/B/lbb_rho10_sp1_f600_sph_dmr2'#_accspot' #used in the previous results
		PulsName='pOS_pulses/lbb_rho10_sp1_f600_sph_burst2_dt'#_dmr2_test'
		#PulsName='res/B/B0Ptest'
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
	#print(x)
	#print(evere)
	print("The chosen energy (keV): ", x[ene]*evere/1e3)
	fluxlcurve_Iene = Flux1[0+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Qene = Flux1[1+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Uene = Flux1[2+ene*3:len(Flux1):3*NEnergy]
            	            
	#Flux=zeros((NPhase,NEnergy,3))
	#print(fluxlcurve_Iene)
	#print(fluxlcurve_Qene)
	#print(fluxlcurve_Uene)


	#phase=list(phi/2/pi)+[1.]
	#phase=list(phi)#+[1.]
	phil = list(phi)
	phase=zeros(NPhase)
	I=zeros(NPhase)#+1)
	Q=zeros(NPhase)#+1)
	U=zeros(NPhase)#+1)
	for t in range(NPhase):#+1):
		#I[t],Q[t],U[t]=Flux[t-1,e]*x[e] 
		I[t],Q[t],U[t]=fluxlcurve_Iene[t-1]*x[ene] ,fluxlcurve_Qene[t-1]*x[ene] ,fluxlcurve_Uene[t-1]*x[ene]
		phase[t] = phil[t]
		#I[t],Q[t],U[t]=fluxlcurve_Iene[t]*x[ene] ,fluxlcurve_Qene[t]*x[ene] ,fluxlcurve_Uene[t]*x[ene]

	p=sqrt(Q**2+U**2)/I*100
	#PA=arctan2(-U,-Q)*90/pi+90
	PA=arctan2(U,Q)*90/pi+90
       	
	#print("ish=",ish)
	#for t in range(NPhase):
	#	print(phase[t], I[t])
        
	#figA.suptitle(r'$\nu={:5.0f}Hz$'.format(nu)+
	#              r'$,\,R_e={:5.1f}km$'.format(R_e)+
	#              #r'$,\,R_e=11,12,14km$'.format(R_e)+
	#              #r'$,\,M=1.0, 1.5, 2.0$'+r'$M_{\odot}$'+',\n'+
	#              r'$,\,M=$'+str(M)+r'$M_{\odot}$'+',\n'+
	#              r'$\,\theta={:5.1f}\degree$'.format(theta[0]*180/pi)+#r'$,{:5.1f}\degree$'.format(40.0)+
	#              r'$,\,i={:5.1f}\degree$'.format(incl*180/pi)+',\n'#r'$,{:5.1f}\degree$'.format(60.0)+',\n'
	#              r'$\rho={:5.1f}\degree$'.format(rho)+', '+
	#              r'$\,E={:6.2f}keV$'.format(x[ene]*evere/1e3),fontsize=fontsize)  


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
			plotAc.tick_params(axis='both', which='major', labelsize=fontsize,direction='in',pad=lpad,top=True,right = True)
			plotAc.set_ylabel(r'$\chi\,[\mathrm{deg}]$',fontsize=fontsize)
			#plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

			if(plot_all):
				plotAp.set_xlim(0,1)
				plotAp.tick_params(axis='both', which='major', labelsize=fontsize,direction='in',pad=lpad,top=True,right = True)
				#plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
				#plotAp.set_ylabel(r'$|\frac{F_{\mathrm{vp}}-F_{\mathrm{acm}}}{F_{\mathrm{vp}}}|$',fontsize=fontsize)
				#plotAp.set_ylabel(r'$F_{\mathrm{Q}}(\varphi)/F_{\mathrm{Q}}^{\mathrm{max}}$',fontsize=fontsize)
				plotAF.set_xlim(0,1)
				# plotAF.locator_params(axis='y', nbins=10)
				#plotAF.set_ylabel(r"$F_{x}(\varphi)/F_{x}^{\mathrm{max}}$",fontsize=fontsize)
				plotAF.set_ylabel(r'$F_{\mathrm{I}}/F_{\mathrm{I}}^{\mathrm{max}}$',fontsize=fontsize)
				#plotAd.set_ylabel(r'$F_{\mathrm{U}}(\varphi)/F_{\mathrm{U}}^{\mathrm{max}}$',fontsize=fontsize)
				plotAp.set_ylabel(r'$F_{\mathrm{Q}}/F_{\mathrm{I}}$',fontsize=fontsize)
				plotAd.set_ylabel(r'$F_{\mathrm{U}}/F_{\mathrm{I}}$',fontsize=fontsize) 
				plotAd.tick_params(axis='both', which='major', labelsize=fontsize,direction='in',pad=lpad,top=True,right = True)
				plotAd.set_xlim(0,1)
				plotAF.tick_params(axis='both', which='major', labelsize=fontsize,direction='in',pad=lpad,top=True,right = True)

				#plotAp.yaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
				#plotAp.set_yticks([-0.9,0.0,0.9])
				#plotAp.set_yticklabels(["-0.9","0.0","0.9"],fontstyle="normal")

				plotAF.set_ylim(0.0,1.1)#(0.1,1.1)
				plotAF.yaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
				#plotAF.set_yticks([0.25,0.5,0.75,1.0])
				#plotAF.set_yticklabels(["0.25","0.50","0.75","1.00"],fontstyle="normal")
				plotAF.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])
				plotAF.set_yticklabels(["0","0.2","0.4","0.6","0.8","1.0"],fontstyle="normal")


                                
	#col=colors[(e*NColors)//NEnergy]
	col = colors[ish]

	PA_VP04 = PA
	phase_VP04 = phase
	if(ish == 1): 
		#find best-phasehift:
		#phshift1, gf1 = find_best_phshift.find_best_phshift(np.array(phase),PA,phase_acm0,PA_acm0)
		#phshift1 = phshift1 -1.0
		#print("Phaseshift for oblate star: ",phshift1)
		#quit()
		#in the end, setting the shift by hand seems still to produce better results		
		#phshift1 = 0.0 #0.019#0.0#0.019#0.195#-0.07#-0.2517#0.0#0.001#0.008#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phshift1 = -0.018#695#  -0.016830883531243623#-0.025#-0.048315
		phase_new = shift_phase(np.array(phase),phshift1)
		for ipha in range(0,len(phase_new)-1):
			if(phase_new[ipha+1] > phase_new[ipha]):
				if(PA[ipha] > 50.0*PA[ipha+1]):
					plotAc.plot(phase_new[ipha:ipha+2],[PA[ipha],180.0],"-",color="blue",linewidth=lwidth)
				else:
					plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-",color="blue",linewidth=lwidth)
				if(plot_all):
					plotAF.plot(phase_new[ipha:ipha+2],I[ipha:ipha+2]/I.max(),color=col)
					##plotAF.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/Q.max(),color="darkblue")
					##plotAF.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/U.max(),color="lightblue")
					#plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/Q.max(),color=col,linewidth=lwidth)
					#plotAd.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/U.max(),color=col,linewidth=lwidth)

					plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/I[ipha:ipha+2],color=col,linewidth=lwidth)
					plotAd.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/I[ipha:ipha+2],color=col,linewidth=lwidth)
		PA0_VP04 = PA
		F0_VP04 = I/I.max()
		Q0_VP04 = Q/Q.max()
		U0_VP04 = U/U.max()
		phase0_VP04 = phase

	if(ish == 2): 
		#phshift2 = 0.0
		#phase_new = shift_phase(np.array(phase),phshift)
		##find best-phasehift:
		#phshift2, gf2 = find_best_phshift.find_best_phshift(np.array(phase),PA,phase_acm0,PA_acm0)
		phshift2 = -0.022371326512027778 #phshift2 -1.0
		print("phaseshift for spherical star: ",phshift2)
		#phshift2 = 0.0 #0.1 #-0.2517#0.0#0.001#0.008#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(np.array(phase),phshift2)
		print(len(phase_new))
		for ipha in range(0,len(phase_new)-1):
			if(phase_new[ipha+1] > phase_new[ipha]):
				if(PA[ipha] > 30.0*PA[ipha+1]):
					plotAc.plot(phase_new[ipha:ipha+2],[PA[ipha],180.0],"-",color="blue",linewidth=lwidth)
				else:                                
					plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-",color=col,linewidth=lwidth)
				if(plot_all):
					plotAF.plot(phase_new[ipha:ipha+2],I[ipha:ipha+2]/I.max(),"--",color=col,linewidth=lwidth)
					###plotAF.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/Q.max(),color="darkgreen")
					###plotAF.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/U.max(),color="lightgreen")
					##plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/Q.max(),color=col,linewidth=lwidth)
					##plotAd.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/U.max(),color=col,linewidth=lwidth)#,linestyle="dashed")

					plotAp.plot(phase_new[ipha:ipha+2],Q[ipha:ipha+2]/I[ipha:ipha+2],"-",color=col,linewidth=lwidth)
					plotAd.plot(phase_new[ipha:ipha+2],U[ipha:ipha+2]/I[ipha:ipha+2],"-",color=col,linewidth=lwidth)
                                        
                                        
		PA0_VP04_2 = PA
		F0_VP04_2 = I/I.max()
		Q0_VP04_2 = Q/Q.max()
		U0_VP04_2 = U/U.max()
		phase0_VP04_2 = phase
	#else:
	#	#plotAc.plot(phase,PA,color=col,marker="o",markersize=1.0)
	#	print("...")



compare_to_arcmancer = True
if(compare_to_arcmancer): 
	#colors = ["green","blue","black"]
	colors = ["black","black","black"]
	for ic in range(0,1):
		if(ic == 0):
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x10_agm.csv"# (copy).csv"
			#datafile = "arcman_res/polar_f600_bb_r12_m1.4_d60_i40_x10_agm.csv"# used in old results
			datafile = "../arcmancer/out3/polar_acc_f600_burst_r12_m1.4_d60_i40_x10_agm.csv"
			#datafile = "../arcmancer/out3/polar_f001_bb_r12_m1.4_d60_i40_x10_obl_img1000.csv"# (copy).csv"
			#datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d60_i40_x10_obl.csv"# (copy).csv"
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
		full_chain= [[] for ix in range(egrid+1)]

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
		for i in range(2,egrid+1): #(1,egrid+1):
			#norm_obsF[:,i-1] = full_chain[i,:]*energy_keV[i-1]/np.max(full_chain[i,:]*energy_keV[i-1])
			#norm_obsF[:,i-1] = full_chain[i,:]/np.max(full_chain[i,:])
			norm_obsF[:,i-1] = full_chain[i,:]/full_chain[1,:]
		norm_obsF[:,0] = full_chain[1,:]/np.max(full_chain[1,:])


		#for i in range(0,egrid):
		#print energy_keV
		ene = 0
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.0#0.2517#0.25#0.2421#0.2517#0.2535#0.069#0.0#-0.195#-0.18#-0.172#0.0
		phase_new = shift_phase(phase,phshift)
		fluxsp1 = full_chain[1,:]


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

				if(two_spots):
					U1 = full_chain[3,:]
					Q1 = full_chain[2,:]

					datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d120_i40_x10_obl2.csv"# (copy).csv"
					input = open(datafile, 'r')
					lines = input.readlines()
					input.close()
					Nchain_size = sum(1 for line in open(datafile))
					print(Nchain_size)
					egrid = 6#5#3
					full_chain= [[] for ix in range(egrid+1)]
					for j in range(0,len(full_chain)):
						for i in range(c_lines,Nchain_size): #not reading comment lines
							parts = lines[i].split(",")
							#print parts
							full_chain[j].append(float(parts[j]))
						parts = lines[c_lines].split(",")
					full_chain = np.array(full_chain)

					U2 = full_chain[3,:]
					Q2 = full_chain[2,:]

					#make 0.5 phase shift
					U2_1 = U2[0:int(len(U2)/2)]
					U2_2 = U2[int(len(U2)/2):len(U2)]
					Q2_1 = Q2[0:int(len(Q2)/2)]
					Q2_2 = Q2[int(len(Q2)/2):len(Q2)]
					U2_new = np.append(U2_2,U2_1)
					Q2_new = np.append(Q2_2,Q2_1)
					
					Utot = U1+U2_new
					Qtot = Q1+Q2_new

					PA=arctan2(Utot,Qtot)*90/pi+90



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
						#plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-o",color=col,markersize="1.0")
						#plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"-",color=col,markersize="1.0")
						if(PA[ipha] > 50.0*PA[ipha+1]):
							plotAc.plot(phase_new[ipha:ipha+2],[PA[ipha],180.0],"--",color=col,linewidth=lwidth,dashes=[2,2])
						else:
							plotAc.plot(phase_new[ipha:ipha+2],PA[ipha:ipha+2],"--",color=col,linewidth=lwidth,dashes=[2,2])
						#linePA.set_dashes([2, 2, 2, 2])
						
				if(ic == 0):
					PA_acm0 = PA#norm_obsF[:,ene]#PA
					F_acm0 = norm_obsF[:,ene]
					Q_acm0 = norm_obsF[:,ene+1]
					U_acm0 = norm_obsF[:,ene+2]
					phase_acm0 = phase_new
					phase_acm0[len(phase_acm0)-1]=1.0

				if(ic == 1):
					PA_acm = PA#norm_obsF[:,ene]#PA
					phase_acm = phase_new
					phase_acm[len(phase_acm)-1]=1.0



		if(plot_all):
			if(two_spots):

				datafile = "../arcmancer/out3/polar_f600_bb_r12_m1.4_d120_i40_x10_obl2.csv"# (copy).csv"
				input = open(datafile, 'r')
				lines = input.readlines()
				input.close()
				Nchain_size = sum(1 for line in open(datafile))
				print(Nchain_size)
				egrid = 6#5#3
				full_chain= [[] for ix in range(egrid+1)]
				for j in range(0,len(full_chain)):
					for i in range(c_lines,Nchain_size): #not reading comment lines
						parts = lines[i].split(",")
						#print parts
						full_chain[j].append(float(parts[j]))
					parts = lines[c_lines].split(",")
				full_chain = np.array(full_chain)
				fluxsp2 = full_chain[1,:]	

				fluxsp2_1 = fluxsp2[0:int(len(fluxsp2)/2)]
				fluxsp2_2 = fluxsp2[int(len(fluxsp2)/2):len(fluxsp2)]
				fluxsp2_new = np.append(fluxsp2_2,fluxsp2_1)

				print(fluxsp2_new)
				norm_obsF[:,ene] = (fluxsp1 + fluxsp2_new)/(np.max(fluxsp1 + fluxsp2_new))


			#for i in range(ene,ene+1):
			for i in range(0,ene+3): #this if plotting all stokes fluxes to same plot
				for ipha in range(0,len(phase_new)-1):
					if(phase_new[ipha+1] > phase_new[ipha]):
						#plotAF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],color=colors[i],marker="o",markersize="1.0")
						if(i == 0):
							plotAF.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],"--",color=colors[i],linewidth=lwidth,dashes=[2,2])
						if(i == 1):
							plotAp.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],"--",color=colors[i],linewidth=lwidth,dashes=[2,2])
						if(i == 2):
							plotAd.plot(phase_new[ipha:ipha+2],norm_obsF[ipha:ipha+2,i],"--",color=colors[i],linewidth=lwidth,dashes=[2,2])


        


#plot_PA_residuals = True
#if(plot_QU):
#	plot_PA_residuals = False
#if(plot_PA_residuals):
if(plot_resid):
	col = "black"
	plotAr.set_xlim(0,1)
	plotAr.set_ylim(-10.0,10.0)

	plotAr.tick_params(axis='both', which='major', labelsize=fontsize,direction='in',pad=lpad,top=True,right = True)


	#plotAd.set_yticks([0,30,60,90,120,150,180])
	#plotAr.tick_params(axis='both', which='major', labelsize=ticksize,direction='in')
	#plotAd.tick_params(axis='y', which='major', labelsize=8)
	#plotAd.tick_params(axis='x', which='major', labelsize=labelsize)
	#plotAd.set_ylabel(r'$\chi_{\mathrm{acm}}/\chi_{\mathrm{vp}}$',fontsize=fontsize)
	plotAr.set_ylabel(r'$\Delta \chi \ [\mathrm{deg}]$',fontsize=fontsize)
	#plotAd.set_ylabel(r'$\frac{\chi_{\mathrm{acm}}-\chi_{\mathrm{vp}}}{\chi_{\mathrm{max}}-\chi_{\mathrm{min}}}$',fontsize=fontsize)
	#plotAd.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

	plotAr.yaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
	plotAr.set_yticks([-5.0,0.0,5.0])
	plotAr.set_yticklabels(["-5","0","5"],fontstyle="normal")
        
	#print(len(phase_VP04),len(phase_acm))
	PA_VP04_interp = interp1d(phase_VP04,PA_VP04)#,"cubic")
	#PA0_VP04_interp = interp1d(phase0_VP04,PA0_VP04)#,"cubic")
	#print(shift_phase(np.array(phase0_VP04),phshift1),PA0_VP04)
	#quit()
	#PA0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift0),PA0_VP04)#,"cubic")
	PA0_VP04_interp = interp1d(shift_phase(np.array(phase0_VP04),phshift1),PA0_VP04,fill_value='extrapolate')# workd with newer scipy
	PA0_VP04_2_interp = interp1d(shift_phase(np.array(phase0_VP04_2),phshift2),PA0_VP04_2,fill_value='extrapolate')# workd with newer scipy

	norm = 1.0#abs(np.max(PA_acm0)-np.min(PA_acm0))

	#print(phase_acm0)
	res_PA = (PA_acm0-PA0_VP04_interp(phase_acm0))/norm
	plotAr.plot(phase_acm0[abs(res_PA) < 10.0],res_PA[abs(res_PA) < 10.0],"-",markersize=5,color="red")
	res_PA_2 = (PA_acm0-PA0_VP04_2_interp(phase_acm0))/norm
	plotAr.plot(phase_acm0[abs(res_PA_2) < 10.0],res_PA_2[abs(res_PA_2) < 10.0],color="darkorange")
	#for ipha in range(0,len(phase_acm0)-1):
	#	#if(phase_new[ipha+1] > phase_new[ipha]):
	#	if(phase_acm0[ipha+1] > phase_acm0[ipha]):
	#		plotAd.plot(phase_acm0[ipha:ipha+2],(PA_acm0[ipha:ipha+2]-PA0_VP04_interp(phase_acm0[ipha:ipha+2]))/norm,color="blue")
	#		#plotAd.plot(phase_acm0[ipha:ipha+2],(PA_acm0[ipha:ipha+2]-PA_VP04_interp(phase_acm0[ipha:ipha+2]))/norm,color="green")


	#plotAp.set_ylim(0.0,0.02)

#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))



#Fine tuning:
#plotAF.axes.get_xaxis().set_ticks([])
#plotAp.axes.get_xaxis().set_ticks([])
#plotAc.axes.get_xaxis().set_ticks([])

plotAF.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAp.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
#plotAc.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAd.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
#plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)
if not (plot_resid):
	plotAc.set_xlabel(r'$\varphi / (2\pi)$',fontsize=fontsize)
else:
	plotAc.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
	plotAr.set_xlabel(r'$\varphi / (2\pi)$',fontsize=fontsize)

plotAc.set_yticks([0,50,100,150])
plotAc.set_yticklabels(["0","50","100","150"],fontstyle="normal")


#align manuallyl ylabels:
labelx = -0.16
plotAF.yaxis.set_label_coords(labelx, 0.5)
plotAd.yaxis.set_label_coords(labelx, 0.5)
plotAp.yaxis.set_label_coords(labelx, 0.5)
plotAc.yaxis.set_label_coords(labelx, 0.5)

plotAF.margins(x =0, y=0)
plotAd.margins(x =0, y=0)
plotAp.margins(x =0, y=0)
plotAc.margins(x =0, y=0)

#Adjusting the ticks and their labels ....
plotAp.set_ylim(-0.05,0.03)
plotAp.set_yticks([-0.04,-0.02,0.0,0.02])
plotAp.set_yticklabels(["-0.04","-0.02","0","0.02"],fontstyle="normal")

plotAd.set_ylim(-0.04,0.04)
plotAd.set_yticks([-0.02,0.0,0.02])
plotAd.set_yticklabels(["-0.02","0","0.02"],fontstyle="normal")

if(plot_resid):
	plotAr.margins(x =0)
	plotAr.yaxis.set_label_coords(labelx, 0.5)


#figA.tight_layout()
figA.subplots_adjust(wspace=0, hspace=0)
figA.subplots_adjust(left=0.15)

#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))
figA.savefig('res/B/plot.pdf',bbox_inches='tight')#.format(e))
figA.clf()



