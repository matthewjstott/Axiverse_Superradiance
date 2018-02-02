#!/usr/bin/python
from __future__ import division
import numpy as np 
from scipy.stats import norm, chi2
import scipy as sp
import scipy.interpolate
import matplotlib.pyplot as plt
import math
import pandas as pd
import ConfigParser
import json
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.special import erfc


def read_in_matrix():
	
	configParser = ConfigParser.RawConfigParser()   
	configParser.read('configuration_card.ini')
	
	n = configParser.getint('General','Number of Axions')
	mo = configParser.getint('Model_Selection','Model')
	c = configParser.getfloat('Model_Selection','Dimension')
	
	a0 = configParser.getfloat('Hyperparameter','a0')
	sa = configParser.getfloat('Hyperparameter','sigma_a')
	b0 = configParser.getfloat('Hyperparameter','b0')
	sb = configParser.getfloat('Hyperparameter','sigma_b')
	
	kmin = configParser.getfloat('Hyperparameter','kmin')
	kmax = configParser.getfloat('Hyperparameter','kmax')
	mmin = configParser.getfloat('Hyperparameter','mmin')
	mmax = configParser.getfloat('Hyperparameter','mmax')
	
	s1 = configParser.getfloat('Hyperparameter','s1')
	s2 = configParser.getfloat('Hyperparameter','s2')
	
	return(n,mo,c,a0,sa,b0,sb,kmin,kmax,mmin,mmax,s1,s2)



def read_in():
	
	configParser = ConfigParser.RawConfigParser()   
	configParser.read('configuration_card.ini')
	
	axm = configParser.getfloat('Input', 'Axion Mass')
	astar = configParser.getfloat('Input', 'Black Hole Spin')
	g = configParser.get('Input', 'Gravitational Constant')
	l = json.loads(configParser.get('Input', 'Quantum Levels l'))
	m = json.loads(configParser.get('Input', 'Quantum Levels m'))
	n = json.loads(configParser.get('Input', 'Overtone Modes')) 
	bhml = configParser.getfloat('Input', 'Lower Black Hole Mass')
	bhmu = configParser.getfloat('Input', 'Upper Black Hold Mass')
		
	axm=float(axm)
	g=float(g)
	
	return(axm,astar,g,l,m,n,bhml,bhmu)

def effective_zvalue(a,b,c,d,fx,fy):
	
	dy = a
	dx = b
	y = c
	x = d
	xtem=[]
	ytem=[]
	dytem=[]
	dxtem=[]
	
	for i in range(len(x)):
		if x[i]<fx[-1]:
			xtem.append(x[i])
			ytem.append(y[i])
			dytem.append(dy[i])
			dxtem.append(dx[i])	
	f=[]
	for i in range(len(xtem)):
		f.append(regge_function(fx,fy,xtem[i]))
	fdx=[]
	for i in range(len(xtem)):
		fdx.append(grad(xtem[i],fx,fy))
	effvar=[]

	
	for i in range(len(ytem)):
		effvar.append(dytem[i]**2+fdx[i]**2*dxtem[i]**2)
	effsd = np.sqrt(effvar)
	
	zscore=[]
	zupper=[]
	zlower=[]
	for i in range(len(ytem)):
		zscore.append((f[i]-ytem[i])/effsd[i])
		zupper.append((1-y[i])/effsd[i])
		zlower.append((0-y[i])/effsd[i])
		
	return(zscore,zupper,zlower,effsd,xtem,ytem,dytem,dxtem)

def probability(z,zupper,zlower):
	prob = 0.5+(0.5*erf(z/np.sqrt(2)))
	exprob = 1 - prob
	prob_upper = (0.5*erfc(zupper/np.sqrt(2)))
	prob_lower = (0.5+ (0.5*erf(zlower/np.sqrt(2))))
	return(prob,exprob,prob_upper,prob_lower)
	
def normalised_probability(prob,prob_upper,prob_lower):
	nprob=[]
	for i in range(len(prob)):
		print 1-prob_upper[i]-prob_lower[i]
		nprob.append(prob[i]/(1-prob_upper[i]-prob_lower[i]))
	return(nprob)	
	

def black_hole_function_map(masses,fx,fy):
	ind=[]
	for i in range(len(masses)):
		fx=np.array(fx)
		#print len(fx)
		tem=find_nearest(fx,3.5)
		#print fy[tem],fx[tem]
		#print tem
	return(ind)	
	
def parameters(bhml,bhmu,g,axm,astar,ma_array):
	alpha={}
	bhms = np.linspace(0.01,200,5000)
	bhm = bhms*2.*10**30*5.6095886*10**35
	for i in range (len(ma_array)):
		alpha[i]=[]
		alpha[i] = bhm*g*ma_array[i]
		
	alphax = bhm*g*axm
	rg = g*bhm
	rp = rg + rg*(1-astar**2)**0.5
	wp = (1/(2*rg))*(astar/(1+(1-astar**2)**0.5))
	a = np.logspace(-2,0,500)
	mm = np.linspace(0.1,30,500)*2.*10**30*5.6095886*10**35*g
	X,Y = np.meshgrid(mm,a)
	time = np.logspace(-18.03,-18.03,1)
	return(bhms,bhm,alpha,rg,rp,wp,X,Y,time)


def black_hole_data():
	df=pd.read_csv('Black_hole_spin_data.csv', sep=',',header=None,encoding='latin-1')
	df.values

	sr_spins=[]
	sr_masses=[]
	sr_spin_up=[]
	sr_spin_low=[]
	sr_mass_up=[]
	sr_mass_low=[]
	sm_spins=[]
	sm_masses=[]
	sm_spin_up=[]
	sm_spin_low=[]
	sm_mass_up=[]
	sm_mass_low=[]

	for i in range (len(df[1][:])):
		if df[1][i] == 'Solar':
			sr_spins.append(float(df[5][i]))
			sr_masses.append(float(df[2][i]))
			sr_spin_up.append(float(df[6][i]))
			sr_spin_low.append(float(df[7][i]))
			sr_mass_up.append(float(df[3][i]))
			sr_mass_low.append(float(df[4][i]))
		if df[1][i] == 'Supermassive':
			sm_spins.append(float(df[5][i]))
			sm_masses.append(float(df[2][i]))
			sm_spin_up.append(float(df[6][i]))
			sm_spin_low.append(float(df[7][i]))
			sm_mass_up.append(float(df[3][i]))
			sm_mass_low.append(float(df[4][i]))
			
	return(sr_spins,sr_masses,sr_spin_up,sr_spin_low,sr_mass_up,sr_mass_low,sm_spins,sm_masses,sm_spin_up,sm_spin_low,sm_mass_up,sm_mass_low)	



def cov_ellipse(cov, q=None, nsig=None, **kwargs):
	if q is not None:
		q = np.asarray(q)
	elif nsig is not None:
		q = 2 * norm.cdf(nsig) - 1
	else:
		raise ValueError()
	r2 = chi2.ppf(q, 2)
	val, vec = np.linalg.eigh(cov)
	width, height = 2 * np.sqrt(val[:, None] * r2)
	rotation = np.degrees(np.arctan2(*vec[::-1, 0]))
	return width, height, rotation

def _rect_inter_inner(x1,x2):
	n1=x1.shape[0]-1
	n2=x2.shape[0]-1
	X1=np.c_[x1[:-1],x1[1:]]
	X2=np.c_[x2[:-1],x2[1:]]
	S1=np.tile(X1.min(axis=1),(n2,1)).T
	S2=np.tile(X2.max(axis=1),(n1,1))
	S3=np.tile(X1.max(axis=1),(n2,1)).T
	S4=np.tile(X2.min(axis=1),(n1,1))
	return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
	S1,S2,S3,S4=_rect_inter_inner(x1,x2)
	S5,S6,S7,S8=_rect_inter_inner(y1,y2)

	C1=np.less_equal(S1,S2)
	C2=np.greater_equal(S3,S4)
	C3=np.less_equal(S5,S6)
	C4=np.greater_equal(S7,S8)

	ii,jj=np.nonzero(C1 & C2 & C3 & C4)
	return ii,jj

def intersection(x1,y1,x2,y2):
	"""
INTERSECTIONS Intersections of curves.
   Computes the (x,y) locations where two curves intersect.  The curves
   can be broken with NaNs or have vertical segments.
usage:
x,y=intersection(x1,y1,x2,y2)
	Example:
	a, b = 1, 2
	phi = np.linspace(3, 10, 100)
	x1 = a*phi - b*np.sin(phi)
	y1 = a - b*np.cos(phi)
	x2=phi
	y2=np.sin(phi)+2
	x,y=intersection(x1,y1,x2,y2)
	plt.plot(x1,y1,c='r')
	plt.plot(x2,y2,c='g')
	plt.plot(x,y,'*k')
	plt.show()
	"""
	ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
	n=len(ii)

	dxy1=np.diff(np.c_[x1,y1],axis=0)
	dxy2=np.diff(np.c_[x2,y2],axis=0)

	T=np.zeros((4,n))
	AA=np.zeros((4,4,n))
	AA[0:2,2,:]=-1
	AA[2:4,3,:]=-1
	AA[0::2,0,:]=dxy1[ii,:].T
	AA[1::2,1,:]=dxy2[jj,:].T

	BB=np.zeros((4,n))
	BB[0,:]=-x1[ii].ravel()
	BB[1,:]=-x2[jj].ravel()
	BB[2,:]=-y1[ii].ravel()
	BB[3,:]=-y2[jj].ravel()

	for i in range(n):
		try:
			T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
		except:
			T[:,i]=np.NaN


	in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

	xy0=T[2:,in_range]
	xy0=xy0.T
	return xy0[:,0],xy0[:,1]

def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return idx
	
'''	
def grad2(x,fx,fy):
	fx=np.array(fx)
	interp = 2000
	indtemp = find_nearest(fx,x)
	tux,tuy = fx[indtemp],fy[indtemp]
	tempx=[]
	tempy=[]
	tempx.extend(fx[indtemp-2:indtemp+2])
	tempy.extend(fy[indtemp-2:indtemp+2])
	tempx=np.array(tempx)
	ntempx = np.linspace(tempx.min(), tempx.max(), interp)
	ntempy = sp.interpolate.interp1d(tempx, tempy, kind='cubic')(ntempx)
	m = (ntempy[interp/2+100]-ntempy[interp/2-100])/(ntempx[interp/2+100]-ntempx[interp/2-100])
	return(m)
'''

def grad(x,fx,fy):
	f = interp1d( fx, fy )
	a=0.01
	xtem = np.linspace(x-a,x+a,99)
	ytem = f(xtem)
	m = (ytem[-1]-ytem[0])/(xtem[-1]-xtem[0])
	return(m)	


	
def regge_region(nx,ny,indx,indx2):
	fx=[]
	fy=[]
	regiesx={}
	regiesy={}
	
	for i in range(len(ny)-1):
		regiesx[i]=[]
		regiesy[i] = []
		regiesx[i] =  nx[i][indx[i]:-1]
		regiesy[i] =  ny[i][indx[i]:-1]
		
	fx.extend(nx[0][0:indx[0]])
	fx.extend(nx[1][indx2[0]:indx[1]])
	fx.extend(nx[2][indx2[1]:indx[2]])
	fx.extend(nx[3][indx2[2]:indx[3]])
	fx.extend(nx[4][indx2[3]:-1])

	fy.extend(ny[0][0:indx[0]])
	fy.extend(ny[1][indx2[0]:indx[1]])
	fy.extend(ny[2][indx2[1]:indx[2]])
	fy.extend(ny[3][indx2[2]:indx[3]])
	fy.extend(ny[4][indx2[3]:-1])	
	
	fx=np.array(fx)
	fy=np.array(fy)	
	return(fx,fy,regiesx,regiesy)
	
def regge_function(fx,fy,x):
	f = interp1d( fx, fy )	
	y=f(x)
	return(y)
	
	
def regge_contour_outline(x1,y1):
	nx={}
	ny={}
	new_length = 2000
	for i in range(0,5):
		nx[i]=[]
		ny[i]=[]
		nx[i] = np.linspace(x1[i].min(), x1[i].max(), new_length)
		ny[i] = sp.interpolate.interp1d(x1[i], y1[i], kind='cubic')(nx[i])
	xi={}
	yi={}
	for i in range(0,4):
		xi[i]=[]
		yi[i]=[]
		xi[i],yi[i]=intersection(nx[i],ny[i],nx[i+1],ny[i+1])
	indx={}
	indy={}
	indx2={}
	indy2={}
	for i in range(0,4):
		indx[i]=[]
		indx2[i]=[]	
		indx[i] = find_nearest(nx[i],xi[i][0])
		indx2[i] = find_nearest(nx[i+1],xi[i][0])
	return(nx,ny,indx,indx2)
	

def superradiance_rates_detweiler(l2,m2,n2,alpha,astar,axm,rp,X,Y,ma_array):
	#axm2=[10**-12,10**-11]
	rad1={}
	rad2={}
	Z={}
	j2=l2
	AA=[]
	for j in range (len(ma_array)):
		A={}
		for i in range(len(m2)):
			A[i]=[]
			A[i] = alpha[j]**(4*l2[i]+5)*(m2[i]*astar-2*ma_array[j]*rp)
		AA.append(A.copy())	
	for i in range(len(m2)):
		rad1[i]=[]
		rad1[i] = (2**(4*l2[i]+2)*math.factorial(2*l2[i]+n2[i]+1))/((l2[i]+n2[i]+1)**(2*l2[i]+4)*math.factorial(n2[i]))
		rad2[i]=[]
		rad2[i] = (math.factorial(l2[i])/(math.factorial(2.*l2[i])*math.factorial(2.*l2[i]+1.)))**2
	B={}
	for i in range (len(m2)):
		B[i]=[]
		B[i] = rad1[i]*rad2[i]
	CC=[]
	for j in range (len(ma_array)):
		C={}
		for i in range(len(m2)):
			C[i] = []
			prd = 1
			axmm=ma_array[j]
			for lj in range(j2[i]):
				prd = prd*((lj+1)**2*(1-astar**2)+((m2[i]*astar-2*axmm*rp)**2))
			C[i] = prd
		CC.append(C.copy())	


	rates=[]
	for j in range(len(ma_array)):
		ratetwo={}
		for i in range(len(m2)):
			ratetwo[i] = []
			ratetwo[i] = AA[j][i]*B[i]*CC[j][i]
		rates.append(ratetwo.copy())
	
	if len(ma_array) > 1:
		rates2=rates
		for j in range(len(rates2)-1):
			print j
			for i in range (0,1):
				rates2[j+1][i] = np.add(rates2[j][i],rates2[j+1][i])
		ratesum2={}
		for i in range(len(rates2[0])):
			ratesum2[i]=[]
			ratesum2[i] = rates2[len(rates2)-1][i]
		
	ZZ=[]

	for k in range(len(ma_array)):
		Z={}
		axm=ma_array[k]
		for i in range(len(l2)):
			Z[i]=[]
			prdc = 1
			for lj in range(l2[i]):
				prdc = prdc*((lj+1)**2*(1-Y**2)+((m2[i]*Y-2*axm*(X + X*(1-Y**2)**0.5))**2))
			Z[i]=(axm*X)**(4*l2[i]+5)*(m2[i]*Y-2*axm*((X + X*(1-Y**2)**0.5)))*B[i]*prdc	
		ZZ.append(Z.copy())	
	'''	
	if len(ma_array) > 1:
		ratessum=ZZ
		for j in range(len(ratessum)-1):
			print j
			for i in range (len(ratessum[j])):
				print i
				ratessum[j+1][i] = np.add(ratessum[j][i],ratessum[j+1][i])
				print ratessum[j+1][i][0]
		ratesum={}
		for i in range(len(ratessum[0])):
			ratesum[i]=[]
			ratesum[i] = ratessum[len(ratessum)-1][i]
	
	print ratesum[0][0][0]
	print ZZ[0][0][0][0]	
	print ZZ[1][0][0][0]	
	'''	
	return(rates,ZZ,ratesum2)

def time_limit(aa,rg,axm):
	ev = 1/(6.582119*10**-19)
	years = 3.154*10**7
	rr = aa/(rg*axm)
	rr = rr*ev
	rr = 1/rr
	rr = rr/years
	print 'Contours for',rr[0],'Years'

	
	
def regge_contour_limits(X,Y,Z,aa):
	g = 6.7071186*10**-57
	cs={}
	p1={}
	v1={}
	x1={}
	y1={}
	for i in range (0,5):
		cs[i]=[]
		cs[i] = plt.contour(X/(2.*10**30*5.6095886*10**35*g),Y,Z[0][i],aa,cmap=plt.cm.bone,alpha=0.0,linewidths = 2.5,zorder=2)
		p1[i]=[]
		v1[i]=[]
		x1[i]=[]
		y1[i]=[]
		p1[i] = cs[i].collections[0].get_paths()[0]
		v1[i] = p1[i].vertices	
		x1[i] = v1[i][:,0]
		y1[i] = v1[i][:,1]	
	return(x1,y1)	