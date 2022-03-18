import numpy
import matplotlib
matplotlib.use('agg')
from numpy import zeros, logspace, linspace, log
from pylab import *
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy import integrate
#from ixpeobssim.irfgen.toy import toy_modf
from polpulse import compf

def shift_signal(phase,signal,shift):
	phi = phase+shift
	for t in range(len(phi)):
		if(phi[t] > 1.0):
			phi[t] = phi[t]-1.0
		if(phi[t] < 0.0):
			phi[t] = phi[t]+1.0
	dim2 = False
	if(signal.ndim == 2):
		dim2 = True
		nene = len(signal[0,:])
		sig_shift = np.zeros((len(signal[:,0]),len(signal[0,:])))
	else:
		nene = 1
		sig_shift = np.zeros((len(signal)))
	for ie in range(0,nene):
		if(dim2):
			sint = interp1d(phase,signal[:,ie],fill_value='extrapolate',kind="linear")
			sig_shift[:,ie] = sint(phi)
		else:
			sint = interp1d(phase,signal[:],fill_value='extrapolate',kind="linear")
			sig_shift[:] = sint(phi)
	return sig_shift


def find_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


# Save a new model where two non-antipodal spots combined
def savePyt(temppytI,temppytQ,temppytU,Filename_new="pulse_test2",Filename_old="pulse_test"):
	#Filename = "pulse_test_ts1"
	Flux = open(Filename_old + 'FF.bin')
	phi = open(Filename_old + 'ff.bin')
	Flux1 = fromfile(Flux)
	phase = fromfile(phi)
	NEnergy =281  #n of energy bins in bin

	#modify Flux to correspond the new content:
	for i in range(0,NEnergy):
		Flux1[0+i*3:len(Flux1):3*NEnergy] = np.append(temppytI[:,i],temppytI[len(temppytI[:,i])-1,i]) 
		Flux1[1+i*3:len(Flux1):3*NEnergy] = np.append(temppytQ[:,i],temppytQ[len(temppytQ[:,i])-1,i]) 
		Flux1[2+i*3:len(Flux1):3*NEnergy] = np.append(temppytU[:,i],temppytI[len(temppytU[:,i])-1,i]) 

	outF = open(Filename_new + 'FF.bin','w')
	outf = open(Filename_new + 'ff.bin','w')
	Flux1.tofile(outF,format="%e")
	phase.tofile(outf,format="%e")


			
# READS MODEL MADE BY drive.py FROM .bin
def readPyt(Filename="pulse_test"):
	#Filename = "pulse_test_ts1"
	Flux = open(Filename + 'FF.bin')
	phi = open(Filename + 'ff.bin')
	Flux1 = fromfile(Flux)
	phi1 = fromfile(phi)
	NEnergy =281 # 500 #50 # n of energy bins in bin
	NPhase = 150 # n of phase bins in bin
	x_l, x_u = -3.7 , -1.2
	evere=.5109989e6 
	IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
	pytE,x_weight=IntEnergy
	pytkeV = (pytE*evere)/1e3
	
	emin = 1 # keV values used in integrating (although integrated values not used at the moment)
	emax = 12

	temppytI = np.zeros((NPhase,NEnergy))
	temppytQ = np.zeros((NPhase,NEnergy))
	temppytU = np.zeros((NPhase,NEnergy))

	pytI = np.zeros((NPhase-1,NEnergy))
	pytQ = np.zeros((NPhase-1,NEnergy))
	pytU = np.zeros((NPhase-1,NEnergy))

	for i in range(0,NEnergy):
		temppytI[:,i] = Flux1[0+i*3:len(Flux1):3*NEnergy]
		temppytQ[:,i] = Flux1[1+i*3:len(Flux1):3*NEnergy]
		temppytU[:,i] = Flux1[2+i*3:len(Flux1):3*NEnergy]
		for j in range(0,NPhase):
			if (np.isnan(temppytI[j,i])):
				temppytI[j,i] = 0.0
				temppytQ[j,i] = 0.0
				temppytU[j,i] = 0.0

	for i in range(NPhase-1):
		for j in range(NEnergy):
			pytI[i,j]=temppytI[i,j]
			pytQ[i,j]=temppytQ[i,j]
			pytU[i,j]=temppytU[i,j]


	pytI2 = np.zeros(NPhase-1)
	
	for p in range(NPhase-1):
		for e in range(find_idx(pytkeV,emin),find_idx(pytkeV,emax)):
			pytI2[p] = pytI2[p] + (pytI[p,e]+pytI[p,e+1])*(pytkeV[e+1]-pytkeV[e])

	pytI2 = 1/2*pytI2


	pytQ2 = np.zeros(NPhase-1)
	
	for p in range(NPhase-1):
		for e in range(find_idx(pytkeV,emin),find_idx(pytkeV,emax)):
			pytQ2[p] = pytQ2[p] + (pytQ[p,e]+pytQ[p,e+1])*(pytkeV[e+1]-pytkeV[e])

	pytQ2 = 1/2*pytQ2

	pytU2 = np.zeros(NPhase-1)
	
	for p in range(NPhase-1):
		for e in range(find_idx(pytkeV,emin),find_idx(pytkeV,emax)):
			pytU2[p] = pytU2[p] + (pytU[p,e]+pytU[p,e+1])*(pytkeV[e+1]-pytkeV[e])
			

	pytU2 = 1/2*pytU2

	pytPD = (numpy.sqrt(pytQ2**2+pytU2**2))/(pytI2)

	pytU2 = pytU2/pytI2
	pytQ2 = pytQ2/pytI2
	
	pytI2 = pytI2/np.max(pytI2)

	#print(len(pytI2),pytI2)
	#print(len(phi1),phi1)
	#exit()

	return pytI, pytQ, pytU, pytkeV, pytI2, pytQ2, pytU2, pytPD, phi1[0:NPhase-1]


def plot2():
	NEnergy = 281 #500 #50
	NPhase = 150
	phase1 = numpy.linspace(0,1,NPhase)
	x_l, x_u = -3.7 , -1.2
	evere=.5109989e6 
	IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
	E,x_weight=IntEnergy
	keV = (E*evere)/1e3
	params = [40.0,60.0,0.1171,-0.0,0.0]
	#phase = np.linspace(0.,1.,149)
	phase = phase1[0:NPhase-1]

	mass = 1.4
	rad = 12.0
	incl = params[0]
	theta = params[1]
	rho = 1.0
	pol = params[2]
	chi = params[3]
	pshift = params[4]


	antipodal= True
	
	if not antipodal:
		I1, Q1, U1, pytkeV, pytI2, pytQ2, pytU2, pytPD, phi = readPyt(Filename="pulses/pulse_compt2_r1t20i60")
		I2, Q2, U2, pytkeV, pytI2, pytQ2, pytU2, pytPD, phi = readPyt(Filename="pulses/pulse_compt1_r1t20i60")
		I3, Q3, U3 = I1-I2, Q1-Q2, U1-U2
	else:
		I2, Q2, U2, pytkeV, pytI2, pytQ2, pytU2, pytPD, phi = readPyt(Filename="pulses/pulse_compt1_r1t20i60")
		I3p, Q3p, U3p, pytkeV, pytI3p, pytQ3p, pytU3p, pytPDp, phi  = readPyt(Filename="pulses/pulse_compt1_r1t120i60")
		pshift = 0.5 
		I3 = shift_signal(phi,I3p,pshift)
		Q3 = shift_signal(phi,Q3p,pshift)
		U3 = shift_signal(phi,U3p,pshift)
		I1, Q1, U1 = I2+I3, Q2+Q3, U2+U3


	save_combined=False #True
	if save_combined:
		savePyt(I1,Q1,U1,Filename_new="pulses/pd_yes0_b_no0_xfc/pulse_phdcorr4_thom2_r1t20t120i60p0",Filename_old="pulses/pd_yes0_b_no0_xfc/pulse_thom1_r1t100i60p0")


	for ene in (find_idx(keV,2),find_idx(keV,6),find_idx(keV,10)):
		PD1 = (numpy.sqrt(Q1[:,ene]**2+U1[:,ene]**2))/(I1[:,ene])
		PD2 = (numpy.sqrt(Q2[:,ene]**2+U2[:,ene]**2))/(I2[:,ene])
		#PD3 = (numpy.sqrt(Q3[:,ene]**2+U3[:,ene]**2))/(I3[:,ene])
		Nphase=150
		PD3 = numpy.zeros((Nphase-1))
		for ip in range(0,Nphase-1):
			if(abs(I3[ip,ene]) > 1e-10):
				PD3[ip] = (numpy.sqrt(Q3[ip,ene]**2+U3[ip,ene]**2))/(I3[ip,ene])
			else:
				PD3[ip] = -100.0

		PA1=arctan2(-U1[:,ene],-Q1[:,ene])*90/pi+90
		PA2=arctan2(-U2[:,ene],-Q2[:,ene])*90/pi+90
		PA3=arctan2(-U3[:,ene],-Q3[:,ene])*90/pi+90

		nI1 = I1[:,ene]/np.max(I1[:,ene])
		#nI2 = I2[:,ene]/np.max(I2[:,ene])
		#nI3 = I3[:,ene]/np.max(I3[:,ene])
		nI2 = I2[:,ene]/np.max(I1[:,ene])
		nI3 = I3[:,ene]/np.max(I1[:,ene])

		fig = plt.figure(figsize=(10,10))
		fig.subplots_adjust(hspace=0)
		plt.rcParams.update({'font.size': 25})
		ax1 = plt.subplot(3,1,1)
		ax2 = plt.subplot(3,1,2)
		ax3 = plt.subplot(3,1,3)
		
		col1 = "black"
		col2 = "blue"
		col3 = "red"
		msz = 7
		ltick=12

		ax1.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
		ax1.plot(phase, nI1, label="both spots", color=col1)
		ax1.plot(phase, nI2, '--', label="spot1", color=col2)
		ax1.plot(phase, nI3, '.', label="spot2", color=col3,markersize=msz)
		ax1.set_xlim(0., 1.)
		ax1.set_ylim(-0.1, 1.1)
		ax1.set_ylabel("$I$/$I_{\mathrm{max}}$")

		ax2.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
		ax2.plot(phase, PD1, label="both spots", color=col1)
		ax2.plot(phase, PD2, '--', label="spot1", color=col2)
		ax2.plot(phase, PD3, '.', label="spot2", color=col3,markersize=msz)
		ax2.set_xlim(0., 1.)
		ax2.set_ylim(-0.01, 0.09)
		#ax2.set_ylim(-0.01, 0.1)
		ax2.set_ylabel("$p$")

		ax3.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
		ax3.plot(phase, PA1, label="both spots", color=col1)
		ax3.plot(phase, PA2, '--', label="spot1", color=col2)
		ax3.plot(phase, PA3, '.', label="spot2", color=col3,markersize=msz)
		ax3.set_xlim(0., 1.)
		ax3.set_ylim(0, 180)
		ax3.set_ylabel("$PA$ (deg)")
		ax3.set_xlabel("Phase $\phi/2\pi$")

		#ax1leg = ax1.legend(loc="lower right")

		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.setp(ax2.get_xticklabels(), visible=False)

		#ax1.tick_params(axis="x", which="both",direction="in",bottom=False,top=True)
		#ax1.tick_params(axis="y", which="both",direction="in",right=True,left=True)
		#ax2.tick_params(axis="y", which="both",direction="in",right=True,left=True)
		#ax3.tick_params(axis="y", which="both",direction="in",right=True,left=True)
		#ax3.tick_params(axis="x", which="both",direction="in",right=True,left=False)


		#plt.tight_layout()
		#fig.savefig("Comparison" + str(np.round(keV[ene],1)) + "keV.pdf")		



	###########################################################################################
	#NEW PART BEGINS:...
	#add evrything to a one plot:

	fig = plt.figure(figsize=(30,30))

	rc("text", usetex=True)
	rc("font", family="serif")
	rc("font",serif="Times")
	lbfontsz = 45#30 
	lwidth= 3.0 #3.5 #2.5 #2.0
	lwidth_dot= 5.0 
	msz = 8
	ltick=20#18#12

	#rc("xtick", labelsize=lbfontsz)
	#rc("ytick", labelsize=lbfontsz)
	#matplotlib.pyplot.rcParams.update({'font.size': lbfontsz})

	matplotlib.pyplot.rcParams.update({'ytick.major.width': lwidth})
	matplotlib.pyplot.rcParams.update({'xtick.major.width': lwidth})
	rc("axes", linewidth=lwidth)

	fig.subplots_adjust(hspace=0)
	fig.subplots_adjust(wspace=0.15)
	plt.rcParams.update({'font.size': lbfontsz})
	plt.rcParams['xtick.major.pad']='20'
	plt.rcParams['ytick.major.pad']='10'
	ax1a = plt.subplot(3,3,1)
	ax2a = plt.subplot(3,3,4)
	ax3a = plt.subplot(3,3,7)

	ax1b = plt.subplot(3,3,2)
	ax2b = plt.subplot(3,3,5)
	ax3b = plt.subplot(3,3,8)

	ax1c = plt.subplot(3,3,3)
	ax2c = plt.subplot(3,3,6)
	ax3c = plt.subplot(3,3,9)

	it = 0
	for ene in (find_idx(keV,2),find_idx(keV,5),find_idx(keV,8)):
		PD1 = (numpy.sqrt(Q1[:,ene]**2+U1[:,ene]**2))/(I1[:,ene])
		PD2 = (numpy.sqrt(Q2[:,ene]**2+U2[:,ene]**2))/(I2[:,ene])


		PA1=arctan2(-U1[:,ene],-Q1[:,ene])*90/pi+90
		PA2=arctan2(-U2[:,ene],-Q2[:,ene])*90/pi+90
		#PA3=arctan2(-U3[:,ene],-Q3[:,ene])*90/pi+90


		Nphase=150
		PD3 = numpy.zeros((Nphase-1))
		PA3=numpy.zeros((Nphase-1))

		for ip in range(0,Nphase-1):
			#print("I3 = ",I3[ip,ene])
			if(abs(I3[ip,ene]) > 1e-10):
				PD3[ip] = (numpy.sqrt(Q3[ip,ene]**2+U3[ip,ene]**2))/(I3[ip,ene])
			else:
				PD3[ip] = -100.0
			#if(it == 2):
			#	print(PD3[ip]," ",I3[ip,ene]," ",Q3[ip,ene]," ",U3[ip,ene])
			if(abs(Q3[ip,ene]) > 1e-10):
				PA3[ip]=arctan2(-U3[ip,ene],-Q3[ip,ene])*90/pi+90
			else:
				PA3[ip] = -100.0

		nI1 = I1[:,ene]/np.max(I1[:,ene])
		#nI2 = I2[:,ene]/np.max(I2[:,ene])
		#nI3 = I3[:,ene]/np.max(I3[:,ene])
		nI2 = I2[:,ene]/np.max(I1[:,ene])
		nI3 = I3[:,ene]/np.max(I1[:,ene])


		dash_le = 4.0
		dash_sp = 4.0

		dot_le = 1.0
		dot_sp = 2.0

		if(it == 0):
			ax1a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1a.plot(phase, nI1, label="both spots", color=col1,linewidth=lwidth)
			ax1a.plot(phase, nI2, '--', label="spot1", color=col2,linewidth=lwidth,dashes=(dash_le,dash_sp)) #same with linestyle='dashed'
			ax1a.plot(phase, nI3, linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))#markersize=msz)#,dashes=(dot_le,dot_sp))
			#ax1a.plot(phase, nI3, label="spot2", color=col3,linewidth=lwidth)
			ax1a.set_xlim(0., 1.)
			ax1a.set_ylim(-0.1, 1.1)
			#ax1a.set_ylabel("$F_{\mathrm{I}}/F_{\mathrm{I}}^{\mathrm{max}}$")#("$I$/$I_{\mathrm{max}}$")
			ax1a.set_ylabel("$F_{I}/F_{I}^{\mathrm{max}}$")#("$I$/$I_{\mathrm{max}}$")

			ax1a.set_yticks([0.0,0.25,0.50,0.75,1.0])
			ax1a.set_yticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])
			ax1a.set_xticks([0.0,0.25,0.50,0.75,1.0])

			ax2a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2a.plot(phase, PD1, label="both spots", color=col1, linewidth=lwidth)
			ax2a.plot(phase, PD2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax2a.plot(phase[PD3!=-100], PD3[PD3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))#, markersize=msz)
			ax2a.set_xlim(0., 1.)
			ax2a.set_ylim(-0.01, 0.25)
			#ax2a.set_ylim(-0.01, 0.1)
			ax2a.set_ylabel("$P_{\mathrm{obs}}$")
			ax2a.set_yticks([0.0,0.05,0.10,0.15,0.20,0.25])
			ax2a.set_yticklabels(["$0.0$","$0.05$","$0.10$","$0.15$","$0.20$","$0.25$"])
			ax2a.set_xticks([0.0,0.25,0.50,0.75,1.0])


			ax3a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3a.plot(phase, PA1, label="both spots", color=col1, linewidth=lwidth)
			ax3a.plot(phase, PA2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax3a.plot(phase[PA3!=-100], PA3[PA3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax3a.set_xlim(0., 1.)
			ax3a.set_ylim(0, 180)
			ax3a.set_ylabel("$PA$ (deg)")
			ax3a.set_xlabel("Phase $\phi/2\pi$")
			ax3a.set_xticks([0.0,0.25,0.50,0.75,1.0])
			ax3a.set_xticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])


			plt.setp(ax1a.get_xticklabels(), visible=False)
			plt.setp(ax2a.get_xticklabels(), visible=False)

		if(it == 1):
			ax1b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1b.plot(phase, nI1, label="both spots", color=col1, linewidth=lwidth)
			ax1b.plot(phase, nI2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax1b.plot(phase, nI3, linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax1b.set_xlim(0., 1.)
			ax1b.set_ylim(-0.1, 1.1)
			#ax1b.set_ylabel("$I$/$I_{\mathrm{max}}$")
			ax1b.set_yticks([0.0,0.25,0.50,0.75,1.0])
			#ax1b.set_yticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])
			ax1b.set_yticklabels(["","","","",""])
			ax1b.set_xticks([0.0,0.25,0.50,0.75,1.0])


			ax2b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2b.plot(phase, PD1, label="both spots", color=col1, linewidth=lwidth)
			ax2b.plot(phase, PD2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax2b.plot(phase[PD3!=-100], PD3[PD3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax2b.set_xlim(0., 1.)
			ax2b.set_ylim(-0.01, 0.25)
			#ax2b.set_ylim(-0.01, 0.1)
			#ax2b.set_ylabel("$p$")
			ax2b.set_yticks([0.0,0.05,0.10,0.15,0.20,0.25])
			#ax2b.set_yticklabels(["$0.0$","$0.05$","$0.10$","$0.15$","$0.20$","$0.25$"])
			ax2b.set_yticklabels(["","","","","",""])
			ax2b.set_xticks([0.0,0.25,0.50,0.75,1.0])

			ax3b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3b.plot(phase, PA1, label="both spots", color=col1, linewidth=lwidth)
			ax3b.plot(phase, PA2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax3b.plot(phase[PA3!=-100], PA3[PA3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax3b.set_xlim(0., 1.)
			ax3b.set_ylim(0, 180)
			#ax3b.set_ylabel("$PA$ (deg)")
			ax3b.set_xlabel("Phase $\phi/2\pi$")
			ax3b.set_xticks([0.0,0.25,0.50,0.75,1.0])
			ax3b.set_xticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])
			ax3b.set_yticks([50.0,100.0,150.0])
			ax3b.set_yticklabels(["","",""])

			plt.setp(ax1b.get_xticklabels(), visible=False)
			plt.setp(ax2b.get_xticklabels(), visible=False)

		if(it == 2):
			ax1c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1c.plot(phase, nI1, label="both spots", color=col1, linewidth=lwidth)
			ax1c.plot(phase, nI2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax1c.plot(phase, nI3, linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax1c.set_xlim(0., 1.)
			ax1c.set_ylim(-0.1, 1.1)
			#ax1c.set_ylabel("$I$/$I_{\mathrm{max}}$")
			ax1c.set_yticks([0.0,0.25,0.50,0.75,1.0])
			#ax1c.set_yticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])
			ax1c.set_yticklabels(["","","","",""])
			ax1c.set_xticks([0.0,0.25,0.50,0.75,1.0])


			ax2c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2c.plot(phase, PD1, label="both spots", color=col1, linewidth=lwidth)
			ax2c.plot(phase, PD2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax2c.plot(phase[PD3!=-100], PD3[PD3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax2c.set_xlim(0., 1.)
			ax2c.set_ylim(-0.01, 0.25)
			#ax2c.set_ylim(-0.01, 0.1)
			#ax2c.set_ylabel("$p$")
			ax2c.set_yticks([0.0,0.05,0.10,0.15,0.20,0.25])
			#ax2c.set_yticklabels(["$0.0$","$0.05$","$0.10$","$0.15$","$0.20$","$0.25$"])
			ax2c.set_yticklabels(["","","","","",""])
			ax2c.set_xticks([0.0,0.25,0.50,0.75,1.0])

			ax3c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3c.plot(phase, PA1, label="both spots", color=col1, linewidth=lwidth)
			ax3c.plot(phase, PA2, '--', label="spot1", color=col2, linewidth=lwidth,dashes=(dash_le,dash_sp))
			ax3c.plot(phase[PA3!=-100], PA3[PA3!=-100], linestyle='dotted', label="spot2", color=col3,linewidth=lwidth_dot,dashes=(dot_le,dot_sp))
			ax3c.set_xlim(0., 1.)
			ax3c.set_ylim(0, 180)
			#ax3c.set_ylabel("$PA$ (deg)")
			ax3c.set_xlabel("Phase $\phi/2\pi$")
			ax3c.set_xticks([0.0,0.25,0.50,0.75,1.0])
			ax3c.set_xticklabels(["$0.0$","$0.25$","$0.50$","$0.75$","$1.0$"])
			ax3c.set_yticks([50.0,100.0,150.0])
			ax3c.set_yticklabels(["","",""])


			plt.setp(ax1c.get_xticklabels(), visible=False)
			plt.setp(ax2c.get_xticklabels(), visible=False)

			
			#way to set the ticks by hand:
			#axr.set_yticks([1e18,1e20,1e22])
			#axr.set_yticklabels(["\$ 10^{18}\$","\$ 10^{20}\$","\$ 10^{22}\$"],fontstyle="normal")

		it = it+1


	ax3a.text(0.4,25.0, '2 keV')
	ax3b.text(0.4,25.0, '5 keV')
	ax3c.text(0.4,25.0, '8 keV')
	fig.savefig("figs/Comparison_allkeV_ap"+str(int(antipodal))+".pdf",bbox_inches='tight')		

	#Use this to fix font issues when adding figs to an article:
	#spath2 = "figs/Comparison_allkeV.pdf"
	#import subprocess
	#subprocess.check_output(['gs','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-dPDFSETTINGS=/prepress','-dEmbedAllFonts=true','-sOutputFile='+spath2[0:len(spath2)-4]+'_emb.pdf','-f',spath2])
        
	#add all Stokes to a one plot:

	fig = plt.figure(figsize=(30,30))

	fig.subplots_adjust(hspace=0)
	plt.rcParams.update({'font.size': 40})
	plt.rcParams['xtick.major.pad']='20'
	ax1a = plt.subplot(3,3,1)
	ax2a = plt.subplot(3,3,4)
	ax3a = plt.subplot(3,3,7)

	ax1b = plt.subplot(3,3,2)
	ax2b = plt.subplot(3,3,5)
	ax3b = plt.subplot(3,3,8)

	ax1c = plt.subplot(3,3,3)
	ax2c = plt.subplot(3,3,6)
	ax3c = plt.subplot(3,3,9)

	it = 0
	for ene in (find_idx(keV,2),find_idx(keV,5),find_idx(keV,8)):

		nI1 = I1[:,ene]/np.max(I1[:,ene])
		#nI2 = I2[:,ene]/np.max(I2[:,ene])
		#nI3 = I3[:,ene]/np.max(I3[:,ene])
		nI2 = I2[:,ene]/np.max(I1[:,ene])
		nI3 = I3[:,ene]/np.max(I1[:,ene])

		if(it == 0):
			ax1a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1a.plot(phase, nI1, label="both spots", color=col1)
			ax1a.plot(phase, nI2, '--', label="spot1", color=col2)
			ax1a.plot(phase, nI3, '.', label="spot2", color=col3,markersize=msz)
			ax1a.set_xlim(0., 1.)
			ax1a.set_ylim(-0.1, 1.1)
			ax1a.set_ylabel("$I$/$I_{\mathrm{max}}$")

			ax2a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2a.plot(phase, Q1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax2a.plot(phase, Q2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax2a.plot(phase, Q3[:,ene]/I3[:,ene], '.', label="spot2", color=col3,markersize=msz)
			ax2a.set_xlim(0., 1.)
			#ax2a.set_ylim(-0.01, 0.09)
			#ax2a.set_ylim(-0.01, 0.1)
			ax2a.set_ylabel("$Q/I$")

			ax3a.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3a.plot(phase, U1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax3a.plot(phase, U2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax3a.plot(phase, U3[:,ene]/I3[:,ene], '.', label="spot2", color=col3,markersize=msz)
			ax3a.set_xlim(0., 1.)
			#ax3a.set_ylim(0, 180)
			ax3a.set_ylabel("$U/I$")
			ax3a.set_xlabel("Phase $\phi/2\pi$")

			plt.setp(ax1a.get_xticklabels(), visible=False)
			plt.setp(ax2a.get_xticklabels(), visible=False)

		if(it == 1):
			ax1b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1b.plot(phase, nI1, label="both spots", color=col1)
			ax1b.plot(phase, nI2, '--', label="spot1", color=col2)
			ax1b.plot(phase, nI3, '.', label="spot2", color=col3,markersize=msz)
			ax1b.set_xlim(0., 1.)
			ax1b.set_ylim(-0.1, 1.1)
			#ax1b.set_ylabel("$I$/$I_{\mathrm{max}}$")

			Q3perI = np.zeros((Nphase-1))
			U3perI = np.zeros((Nphase-1))
			for ip in range(0,Nphase-1):
				if(abs(I3[ip,ene]) >  2164180067701):
					Q3perI[ip] = Q3[ip,ene]/I3[ip,ene]
					U3perI[ip] = U3[ip,ene]/I3[ip,ene]
				

			ax2b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2b.plot(phase, Q1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax2b.plot(phase, Q2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax2b.plot(phase, Q3perI, '.', label="spot2", color=col3,markersize=msz)
			ax2b.set_xlim(0., 1.)
			#ax2b.set_ylabel("$Q/I$")

			ax3b.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3b.plot(phase, U1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax3b.plot(phase, U2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax3b.plot(phase, U3perI, '.', label="spot2", color=col3,markersize=msz)
			ax3b.set_xlim(0., 1.)
			#ax3b.set_ylabel("$U/I$")
			ax3b.set_xlabel("Phase $\phi/2\pi$")

			plt.setp(ax1b.get_xticklabels(), visible=False)
			plt.setp(ax2b.get_xticklabels(), visible=False)

		if(it == 2):
			ax1c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax1c.plot(phase, nI1, label="both spots", color=col1)
			ax1c.plot(phase, nI2, '--', label="spot1", color=col2)
			ax1c.plot(phase, nI3, '.', label="spot2", color=col3,markersize=msz)
			ax1c.set_xlim(0., 1.)
			ax1c.set_ylim(-0.1, 1.1)
			#ax1c.set_ylabel("$I$/$I_{\mathrm{max}}$")

			ax2c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax2c.plot(phase, Q1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax2c.plot(phase, Q2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax2c.plot(phase, Q3[:,ene]/I3[:,ene], '.', label="spot2", color=col3,markersize=msz)
			ax2c.set_xlim(0., 1.)
			#ax2c.set_ylabel("$Q/I$")

			ax3c.tick_params(bottom=True, top=True, left=True, right=True, direction="in", length=ltick)
			ax3c.plot(phase, U1[:,ene]/I1[:,ene], label="both spots", color=col1)
			ax3c.plot(phase, U2[:,ene]/I2[:,ene], '--', label="spot1", color=col2)
			ax3c.plot(phase, U3[:,ene]/I3[:,ene], '.', label="spot2", color=col3,markersize=msz)
			ax3c.set_xlim(0., 1.)
			#ax3c.set_ylabel("$U/I$")
			ax3c.set_xlabel("Phase $\phi/2\pi$")


			plt.setp(ax1c.get_xticklabels(), visible=False)
			plt.setp(ax2c.get_xticklabels(), visible=False)

			
			#way to set the ticks by hand:
			#axr.set_yticks([1e18,1e20,1e22])
			#axr.set_yticklabels(["\$ 10^{18}\$","\$ 10^{20}\$","\$ 10^{22}\$"],fontstyle="normal")

		it = it+1

	fig.tight_layout()
	fig.savefig("figs/Stokes_allkeV_ap"+str(int(antipodal))+".pdf")		


if __name__ == "__main__":
	plot2()
