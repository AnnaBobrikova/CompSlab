import numpy as np
from matplotlib.pyplot import *

def read_pas(datafile):
	Nchain_size = sum(1 for line in open(datafile))
	c_lines = 0
	egrid = 5
	full_chain= [[] for x in range(egrid+1)]
	input = open(datafile, 'r')
	lines = input.readlines()
	input.close()
	for j in range(0,len(full_chain)):
		for i in range(c_lines,Nchain_size): #not reading comment lines                                                                
			parts = lines[i].split(",")                                                                             
			full_chain[j].append(parts[j])
		parts = lines[c_lines].split(",")
	full_chain = np.array(full_chain)

	#phi = full_chain[0]
	#chi_tot = full_chain[1]
	#chi_0 = full_chain[2]
	#chi_1 = full_chain[3]
	#chi_prime = full_chain[4]
	#chi_prime_sph = full_chain[5]

	#chis = [chi_tot,chi_0,chi_1,chi_prime,chi_prime_sph]

	#print(phi)
	#print(chi_tot)
	#exit()

	return full_chain

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t60i40nu600.txt")
phi = alldata[0]
chi_tot, chi_0, chi_1, chi_prime, chi_prime_sph = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphTrue_t60i40nu600.txt")
chi_primeT = alldata[4]
chi_prime_sphT = alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t60i40nu300.txt")
chi_totf3, chi_0f3, chi_1f3, chi_primef3, chi_prime_sphf3 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t60i40nu100.txt")
chi_totf1, chi_0f1, chi_1f1, chi_primef1, chi_prime_sphf1 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t40i60nu600.txt")
chi_tot_it46, chi_0_it46, chi_1_it46, chi_prime_it46, chi_prime_sph_it46 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t20i80nu600.txt")
chi_tot_it28, chi_0_it28, chi_1_it28, chi_prime_it28, chi_prime_sph_it28 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]


#to compare with Vlad:
alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t10i50nu600.txt")
chi_tot_it51, chi_0_it51, chi_1_it51, chi_prime_it51, chi_prime_sph_it51 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t30i50nu600.txt")
chi_tot_it53, chi_0_it53, chi_1_it53, chi_prime_it53, chi_prime_sph_it53 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t60i50nu600.txt")
chi_tot_it56, chi_0_it56, chi_1_it56, chi_prime_it56, chi_prime_sph_it56 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t90i50nu600.txt")
chi_tot_it59, chi_0_it59, chi_1_it59, chi_prime_it59, chi_prime_sph_it59 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]
phi_it59 = alldata[0]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphTrue_t90i50nu600.txt")
chi_tot_it59T, chi_0_it59T, chi_1_it59T, chi_prime_it59T, chi_prime_sph_it59T = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]
phi_it59T = alldata[0]

alldata = read_pas("res/B/pa_eqs_figs/pa_SphTrue_t45i40nu600.txt")
chi_tot_it445T, chi_0_it445T, chi_1_it445T, chi_prime_it445T, chi_prime_sph_it445T = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]


alldata = read_pas("res/B/pa_eqs_figs/pa_SphFalse_t45i40nu600_dmr.txt")
chi_tot_it445, chi_0_it445, chi_1_it445, chi_prime_it445, chi_prime_sph_it445 = alldata[1], alldata[2], alldata[3], alldata[4], alldata[5]


print(chi_tot)

figA = figure(figsize=(14,12), dpi=300) #8,6
#rc("font", family="serif")
#rc("font",serif="Times")
matplotlib.pyplot.figure(1)
lbfontsz = 35 
lwidth= 3.5 #2.5#2.0#1.5 

#labelsize=35#30#20
#fontsize=35#50#35#25
#ticksize=25

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

plotAp=figA.add_subplot(1,1,1)

#if(plot_all):
#	plotAF=figA.add_subplot(4,1,1,yscale='linear') 
#	plotAp=figA.add_subplot(4,1,2)      #
#	plotAc=figA.add_subplot(4,1,3)      #
#	plotAd=figA.add_subplot(4,1,4)      #
#	#plotAd=figA.add_subplot(4,1,3)      #
#	#plotAc=figA.add_subplot(4,1,4)      #
#else:	
#	plotAc=figA.add_subplot(2,1,1)      #
#	plotAd=figA.add_subplot(2,1,2)      #



#figA.tight_layout()
figA.subplots_adjust(left=0.15)
figA.subplots_adjust(bottom=0.15)
#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))

plotAp.tick_params(axis="both", which="both", pad=14)#,length=6, width=lwidth)

cols = ["blue","green","red","black","darkorange","magenta"]

#plotAp.set_ylim(-15.0,50.0)
plotAp.set_ylabel(r'$\chi\,[\degree]$',fontsize=lbfontsz)
plotAp.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=lbfontsz)


chi_totn = np.zeros((len(chi_tot)))
chi_totsph = np.zeros((len(chi_tot)))
chi_totsph_it51 = np.zeros((len(chi_tot)))
chi_totsph_it53 = np.zeros((len(chi_tot)))
chi_totsph_it56 = np.zeros((len(chi_tot)))
chi_totsph_it59 = np.zeros((len(chi_tot_it59)))
chi_totsph_it59T = np.zeros((len(chi_tot_it59T)))
chi_totsph_it445T = np.zeros((len(chi_tot_it445T)))
print(len(phi_it59))
print(len(phi))

#phi_float = np.zeros((len(phi)))

for ij in range(0,len(chi_tot)):
	#phi_float[ij] = float(phi[ij])
	#print(phi[ij])
	#chi_totn[ij] = float(chi_tot[ij])# (float(chi_tot[ij])/2.0)- 90.0
	#if(chi_totn[ij] > 180.0):
	#	chi_totn[ij] = chi_totn[ij] - 180.0
	#if(chi_totn[ij] < 0.0):
	#	chi_totn[ij] = chi_totn[ij] + 180.0
	chi_totsph[ij] = float(chi_0[ij])+float(chi_prime_sph[ij])
	chi_totsph_it51[ij] = float(chi_0_it51[ij])+float(chi_prime_sph_it51[ij])
	chi_totsph_it53[ij] = float(chi_0_it53[ij])+float(chi_prime_sph_it53[ij])
	chi_totsph_it56[ij] = float(chi_0_it56[ij])+float(chi_prime_sph_it56[ij])
	#chi_totsph_it59[ij] = float(chi_0_it59[ij])+float(chi_prime_sph_it59[ij])
	chi_totsph_it445T[ij] = float(chi_0_it445T[ij])+float(chi_prime_sph_it445T[ij])

for ij in range(0,len(chi_tot_it59)):
	print(phi_it59[ij])
	chi_totsph_it59[ij] = float(chi_0_it59[ij])+float(chi_prime_sph_it59[ij])
	chi_totsph_it59T[ij] = float(chi_0_it59T[ij])+float(chi_prime_sph_it59T[ij])

#plotAp.plot(phi,chi_tot,"-",color=cols[5])
#plotAp.plot(phi,chi_totsph,"--",color=cols[5])
#plotAp.plot(phi,chi_0,"-",color=cols[1])

#plotAp.plot(phi,chi_1,"-",color=cols[3])
#plotAp.plot(phi,chi_prime,"-",color=cols[0])
#plotAp.plot(phi,chi_prime_sph,"-",color=cols[2])

#plotAp.plot(phi,chi_1f3,"--",color=cols[3])
#plotAp.plot(phi,chi_primef3,"--",color=cols[0])
#plotAp.plot(phi,chi_prime_sphf3,"--",color=cols[2])

#plotAp.plot(phi,chi_1f1,"--",color=cols[3])
#plotAp.plot(phi,chi_primef1,"--",color=cols[0])
#plotAp.plot(phi,chi_prime_sphf1,"--",color=cols[2])



#plotAp.plot(phi,chi_1_it46,"--",color=cols[3])
#plotAp.plot(phi,chi_prime_it46,"--",color=cols[0])
#plotAp.plot(phi,chi_prime_sph_it46,"--",color=cols[2])

dlen = 10
dspc = 10
dpnt = 4

#plotAp.plot(phi,chi_1_it28,"-.",color=cols[3],dashes=([dlen, dspc, dpnt, dspc]))
#plotAp.plot(phi,chi_prime_it28,"-.",color=cols[0],dashes=([dlen, dspc, dpnt, dspc]))
#plotAp.plot(phi,chi_prime_sph_it28,"-.",color=cols[2],dashes=([dlen, dspc, dpnt, dspc]))



#Compare with Vlad:
cols=['brown',
        'red',
        # 'xkcd:orange',
        # 'xkcd:dark yellow',
        # 'xkcd:dark yellow green',
        'green',
        # 'xkcd:dark cyan',
        # 'xkcd:blue',
        'purple'   
]

#plotAp.plot(phi_it59,chi_tot_it59,"-",color=cols[0])
#plotAp.plot(phi_it59,chi_totsph_it59T,"--",color=cols[0])
#plotAp.plot(phi,chi_tot_it53,"-",color=cols[1])
#plotAp.plot(phi,chi_totsph_it53,"--",color=cols[1])
#plotAp.plot(phi,chi_tot_it56,"-",color=cols[2])
#plotAp.plot(phi,chi_totsph_it56,"--",color=cols[2])
#plotAp.plot(phi,chi_tot_it51,"-",color=cols[3])
#plotAp.plot(phi,chi_totsph_it51,"--",color=cols[3])




def foldchi(c):
	#return (c+90)%180-90
	print(len(c))
	print(c)
	fc = np.zeros((len(c)))
	for ic in range(len(c)):
		print(c[ic], fc[ic])
		if(c[ic] > 90.0):
			fc[ic] = c[ic]-180.0
		elif(c[ic] < -90.0): 
			fc[ic] = c[ic]+180.0
		else:
			fc[ic] = c[ic]
	return fc


phi_float = phi.astype(np.float)
chi_tot_it445 = chi_tot_it445.astype(np.float)

phi_shifted = phi_float-0.5

#plotAp.plot(phi_shifted,foldchi(chi_tot_it445),"-",color=cols[0])
#plotAp.plot(phi_shifted,foldchi(chi_totsph_it445T),"--",color=cols[0])

plotAp.plot(phi_shifted[0:len(phi_shifted)/2],foldchi(chi_tot_it445[len(phi_shifted)/2:len(phi_shifted)]),"-",color=cols[0])
plotAp.plot(phi_shifted[0:len(phi_shifted)/2],foldchi(chi_totsph_it445T[len(phi_shifted)/2:len(phi_shifted)]),"--",color=cols[0])
plotAp.plot(phi_shifted[len(phi_shifted)/2:len(phi_shifted)],foldchi(chi_tot_it445[0:len(phi_shifted)/2]),"-",color=cols[0])
plotAp.plot(phi_shifted[len(phi_shifted)/2:len(phi_shifted)],foldchi(chi_totsph_it445T[0:len(phi_shifted)/2]),"--",color=cols[0])

plotAp.set_xlim(-0.5,0.5)
plotAp.set_ylim(-90.0,90.0)

#plotAp.plot(phi,chi_primeT,"--",color=cols[0])
#plotAp.plot(phi,chi_prime_sphT,"--",color=cols[2])

print("next saving the plot...")

figA.savefig('res/B/plot.pdf')#.format(e))
figA.clf()




