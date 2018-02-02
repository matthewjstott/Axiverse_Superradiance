#!/usr/bin/python
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import pyplot
plt.rcParams["figure.facecolor"] = 'w'
plt.rcParams["axes.facecolor"] = 'w'
plt.rcParams["savefig.facecolor"] = 'w'

colour_pallette1=['#045a8d','#2b8cbe','#74a9cf','#bdc9e1','#f1eef6']



def conf_legend():
	conf=[99,95,68]
	concolor=['orangered','red','maroon']
	Ep_handle={}
	Ep_label={}
	for i in range(0,3):
		Ep_handle[i]=[]
		Ep_label[i]=[]
		Ep_handle[i] = [mpatches.Patch(color=concolor[i], alpha=0.6, linewidth=0)]
		Ep_label[i] = [u'${0}\%$ CL'.format(conf[i])]

	handles2=[]
	labels2=[]
	for i in range(0,3):	
		handles2.extend(Ep_handle[i])
		labels2.extend(Ep_label[i])
		
	legend22 = plt.legend(handles2,labels2,loc='center right',bbox_to_anchor = [0.9325, 0.61],
	           ncol=2,prop={'size':12},numpoints=1)
	pyplot.gca().add_artist(legend22)
	plt.legend(loc='center right',bbox_to_anchor = [0.99, 0.27])
	
	
def quantum_levels_legend(l):
	p_handle={}
	p_label={}
	for i in range(0,5):
		p_handle[i]=[]
		p_label[i]=[]
		p_handle[i] = [mpatches.Patch(color=colour_pallette1[i], alpha=1.0, linewidth=1.5)]
		p_label[i] = [u'$l=m={0}$'.format(l[i])]

	plt.text(13.11, 0.34, r'$\mu_{\rm ax}=10^{-11}eV$', fontsize=15,bbox={'facecolor':'white', 'alpha':1.0, 'pad':12})
	#handle, label = ax.get_legend_handles_labels()
	handles=[]
	labels=[]
	for i in range(0,5):
		handles.extend(p_handle[i])
		labels.extend(p_label[i])
		
	legend2 = plt.legend(handles,labels,loc='lower right',
	           ncol=2,prop={'size':12},numpoints=1)
	pyplot.gca().add_artist(legend2)
	
	
def regge_plane_plot(x1,y1,sr_spins,sr_masses,sr_spin_up,sr_spin_low,sr_mass_up,sr_mass_low):
	fig, ax = plt.subplots(figsize=(10,6))
	for i in range(4,-1,-1):
		ax.fill_between(x1[i], y1[i], 1,facecolor=colour_pallette1[i],linewidth=2.0,zorder=2)

	labels=(r'$\rm Continuum\ Fit \ Black$'
	'\n' 
	r'$\rm Hole \ Data$') 
	ax.errorbar(sr_masses, sr_spins, yerr=[sr_spin_up,sr_spin_low], xerr=[sr_mass_up,sr_mass_low], fmt='o',color='k',label=labels)
	plt.legend(loc='lower right',prop={'size':12})
	plt.xlabel(r'$\rm Black \ Hole \ Mass \ \left(\rm{M_{\rm BH}} \ / M_{\odot} \right)$', ha='center', va='center',size=20,labelpad=15)
	plt.ylabel(r'$\rm Black \ Hole \ Spin \ \left( a_{*}\right)$',size=21)	
	plt.ylim(0,1)
	plt.xlim(0,x1[4].max())
	
def regge_region_plot(fx,fy,blackholes,rt,xtem,ytem,dytem,dxtem,regiesx,regiesy):
	plt.plot(fx,fy,linestyle='-',color='black')
	
	plt.fill_between(fx, fy,1, color='Gold',alpha=0.3)
	plt.xlim(fx.min(),fx.max())
	
	if rt == True:
		for i in range(len(regiesx)):
			plt.plot(regiesx[i],regiesy[i],linestyle='--',color='black')
	if blackholes == True:
		plt.errorbar(xtem, ytem, yerr=dytem, xerr=dxtem, fmt='o',color='k')
	plt.xlabel(r'${\rm M_{BH}} \left( M_{\odot} \right)$', ha='center', va='center',size=20,labelpad=15)
	plt.ylabel(r'$ \left( a_{*}\right)$',size=21)	
	plt.ylim(0,1)
	plt.xlim(2,20)
	
	
def intersection_plot(nx,ny,indx,indx2):
	plt.plot(nx[4][indx2[3]], ny[4][indy2[3]], 'ro')
	plt.plot(nx[0][0:indx[0]],ny[0][0:indx[0]])
	plt.plot(nx[1][indx2[0]:indx[1]],ny[1][indx2[0]:indx[1]])
	plt.plot(nx[2][indx2[1]:indx[2]],ny[2][indx2[1]:indx[2]])
	plt.plot(nx[3][indx2[2]:indx[3]],ny[3][indx2[2]:indx[3]])
	plt.plot(nx[4][indx2[3]:-1],ny[4][indx2[3]:-1])	
	
	
def superradiance_rates_plot(alpha,rates,ratesum):
	for i in range(0,5):
		plt.plot(alpha[i],rates[i],linewidth=2)
		plt.plot(alpha[i],ratesum[i],linewidth=2)
	plt.yscale('log')
	plt.xlabel(r'$\mu_{\rm ax}  r_g$', size=24,labelpad=4.15)
	plt.ylabel(r'$ \log_{10}(M_{\rm BH} \ IM(\omega))$',size=21,labelpad=2)
	plt.xlim(0,2.55)
	plt.ylim(10**-16.5,10**-6.5)

		