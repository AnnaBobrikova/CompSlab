
#import:
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
# m=fromfile(AtmName+'m.bin')
# print(m)
# exit()
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil, argmin
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d#,CubicSpline 
from scipy.interpolate import CubicSpline 
from scipy.interpolate import interp2d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect


#AtmName='../../../CompSlab/res/ixpe_spec_final/tau10_te01/grid_res_pbi_NN_no0/atmos_thom_p0'#

#AtmName='Comp'#

print("Four command line arguments required: source file, poldeg_intr, ang/ene, spath")



if len(sys.argv) < 4:
	print("Not enough arguments!!!")
	exit()

AtmName = str(sys.argv[1])
plot_spec = str(sys.argv[2]) #ang or ene
spath = str(sys.argv[3])

#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

NEnergy =  100 #20 #281 #20 #281 #50 #281 # 50# 101 # number of energy points (x)
NMu = 9 #3 #22 #3 #22 # 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
x_l, x_u = -3.7 , .3 #lower and upper bounds of the log_10 energy span
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function 
IntZenith = leggauss(NZenith) #  sample points and weights for integrations over zenith angle in positive and negative directions together

mu,mu_weight=IntZenith
x,x_weight=IntEnergy

keV = (x*evere)/1e3


T = 0.002  #~  1.0219978 keV, for comparison to Simpl

def find_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def poldegfromfile(x, e, d):
	inI = open(AtmName+'I.bin')
	inx = open(AtmName+'x.bin')
	inm = open(AtmName+'m.bin')
	x2=fromfile(inx)
	keV2=(x2*evere)/1e3
	mu=fromfile(inm)
	#print(mu)
	#exit()
	NEnergy=len(x2)
	NZenith=len(mu)
	NMu=NZenith/2
	Intensity=fromfile(inI).reshape((NEnergy,NZenith,2))

	Intensity[:,:,0] = np.nan_to_num(Intensity[:,:,0])
	Intensity[:,:,1] = np.nan_to_num(Intensity[:,:,1])

	fI = interp1d(x2, Intensity[:,d,0], kind="linear")
	I = fI(x)
	fQ = interp1d(x2, Intensity[:,d,1], kind="linear")
	Q = fQ(x)

	
	if(I[e]<1e-10):
		return 0

	p = Q[e]/I[e]

	return p
    

PDs = np.zeros((NEnergy,NZenith))
for e in range(NEnergy):
	for d in range(NZenith):
		#print(poldegfromfile(x, e, d))
		PDs[e,d] = poldegfromfile(x, e, d)*100.0



fig = figure(figsize=(10,8), dpi=300)
rc("text", usetex=True)
rc("font", family="serif")
rc("font",serif="Times")
plt.figure(1)

lbfontsz = 30#25#16#9#5#7#8
lwidth=2.0#3.0#1.5
dotsize=2.5

rc("xtick", labelsize=lbfontsz)
rc("ytick", labelsize=lbfontsz)
rc("axes", linewidth=lwidth)


plt.rcParams.update({'font.size': lbfontsz})
plt.rcParams.update({'axes.linewidth': lwidth})
plt.rcParams.update({'axes.labelsize': lbfontsz})
plt.rcParams.update({'axes.titlesize': lbfontsz})
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'xtick.labelsize': lbfontsz})
plt.rcParams.update({'ytick.labelsize': lbfontsz})
plt.rcParams.update({'lines.linewidth': lwidth})

plt.tight_layout()
plt.gcf().subplots_adjust(left=0.15)
plt.gcf().subplots_adjust(bottom=0.15)
ax = plt.subplot(1,1,1)
frame1 = plt.gca()

ax.xaxis.set_tick_params(width=lwidth,length=20.0)
ax.yaxis.set_tick_params(width=lwidth,length=20.0)
ax.tick_params(axis="y", which="both",direction="in",right=True,pad=8)
ax.tick_params(axis="x", which="both",direction="in",top=True,pad=10)
ax.tick_params(axis="both", which="minor",width=lwidth,length=10.0)


#cols = ["blue","magenta","red","pink","black", "grey"]
cols = ["black","blue","green","darkorange","red", "magenta"]

mus = IntZenith[0]

shift=1.0
        
id1 = find_idx(keV,2.0/shift)
id2 = find_idx(keV,5.0/shift)
id3 = find_idx(keV,8.0/shift)
id4 = find_idx(keV,12.0/shift)
id5 = find_idx(keV,18.0/shift)

#id1 = find_idx(keV,ene_scs[1])
#id2 = find_idx(keV,ene_scs[2])
#id3 = find_idx(keV,ene_scs[3])
#id4 = find_idx(keV,ene_scs[5])
#id5 = find_idx(keV,ene_scs[10])


print("plotting poldeg(\mu) for energies (keV):")
print(keV[id1]," ", keV[id2]," ",keV[id3])#," ",keV[id4]," ",keV[id5])


print("id1=",id1,", id2=",id2,", id3=",id3)

#plot_spec="ang"#"ang"#"neither" #"ene" #"ang" #"ene"
if(plot_spec=="ang"):
	#ax.plot(mus[NMu-1:],PDs[id1,NMu-1:],"-",color=cols[0],label="(1+z)*2.0 keV")
	ax.plot(mus[NMu:],PDs[id1,NMu:],"-",color=cols[0],label="2.0 keV",linewidth=lwidth)
	ax.plot(mus[NMu:],PDs[id2,NMu:],"-",color=cols[1],label="5.0 keV",linewidth=lwidth)
	ax.plot(mus[NMu:],PDs[id3,NMu:],"-",color=cols[2],label="8.0 keV",linewidth=lwidth)
	ax.plot(mus[NMu:],PDs[id4,NMu:],"-",color=cols[3],label="12.0 keV",linewidth=lwidth)
	ax.plot(mus[NMu:],PDs[id5,NMu:],"-",color=cols[4],label="18.0 keV",linewidth=lwidth)
	ax.plot(mus,np.zeros((len(mus))),"--",color="black",linewidth=lwidth) 

        
        
ide0 = find_idx(mus,0.01)        
ide1 = find_idx(mus,0.2)
ide2 = find_idx(mus,0.4)
ide3 = find_idx(mus,0.6)
ide4 = find_idx(mus,0.8)
ide5 = find_idx(mus,1.0)


if(plot_spec=="ene"):
	ax.plot(keV[:],PDs[:,ide0],"-",color=cols[0],label="$\mu$ = 0.0")
	ax.plot(keV[:],PDs[:,ide1],"-",color=cols[1],label="$\mu$ = 0.2")
	ax.plot(keV[:],PDs[:,ide2],"-",color=cols[2],label="$\mu$ = 0.4")
	ax.plot(keV[:],PDs[:,ide3],"-",color=cols[3],label="$\mu$ = 0.6")
	ax.plot(keV[:],PDs[:,ide4],"-",color=cols[4],label="$\mu$ = 0.8")
	ax.plot(keV[:],PDs[:,ide5],"-",color=cols[5],label="$\mu$ = 1.0")





PD_burst = -11.71*(1.0-mus[NMu:])/(1+3.582*mus[NMu:])
#ax.plot(mus[NMu:],PD_burst,"-.",color=cols[0],label="Eq for bursts",linewidth=lwidth)
#ax.legend()
#ax.tick_params(axis="y", which="both",direction="in",right=True)
#ax.tick_params(axis="x", which="both",direction="in",top=True)


if(plot_spec=="ang"):
	ax.set_xlabel('$\mu$',fontsize=lbfontsz)
	ax.set_ylabel('$P(\mu) \%$',fontsize=lbfontsz)
	ax.set_xlim(0.0,1.0)#(0.01,1.0)#(0.01, 1.)
	ax.plot(mus[NMu:],PD_burst,"-.",color=cols[0],label="Eq for bursts",linewidth=lwidth)
        
if(plot_spec=="ene"):
	ax.set_xlabel('$E$ (keV)',fontsize=lbfontsz)
	ax.set_ylabel('$P(E) \%$',fontsize=lbfontsz)
	#ax.set_xlim(0.1,50.0)
	ax.set_xlim(1.0,25.0)
	ax.set_ylim(-35.0,15.0)
	ax.set_xscale('log')
	ax.set_xticks([1.0,10.0])
	ax.set_xticklabels(["$1$","$10$"])


savefig(spath)
plt.close()

#USE THE FOLLOWING TO FIX FONT ISSUES WHEN ATTACHING THE FIGS THE A PAPER:
#import subprocess
#subprocess.check_output(['gs','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-dPDFSETTINGS=/prepress','-dEmbedAllFonts=true','-sOutputFile='+spath[0:len(spath)-4]+'_emb.pdf','-f',spath])
