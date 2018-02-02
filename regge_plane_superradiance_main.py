####################################################

from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
import scipy as sp
import functions
import plotting_functions
import diag
np.set_printoptions(threshold=np.nan)

####################################################

class axiverse_parameters(object):
	
	def __init__(self):
		self.n,self.mo,self.c,self.a0,self.sa,self.b0,self.sb,self.kmin,self.kmax,self.mmin,self.mmax,self.s1,self.s2 = functions.read_in_matrix()
	
	def diagonalisation(self):
		self.ma_array,self.fef = diag.diag(self.n,self.mo,self.c,self.a0,self.sa,self.b0,self.sb,self.kmin,self.kmax,self.mmin,self.mmax,self.s1,self.s2)
		global axms
		axms = [10**-11.0,10**-11.0,10**-11]
	

class superradiance_calculator(object):
	
	def __init__(self):
		
		self.axm,self.astar,self.g,self.l,self.m,self.n,self.bhml,self.bhmu = functions.read_in()
		
		self.sr_spins,self.sr_masses,self.sr_spin_up,self.sr_spin_low,self.sr_mass_up,self.sr_mass_low,self.sm_spins,self.sm_masses,self.sm_spin_up,self.sm_spin_low,self.sm_mass_up,self.sm_mass_low = functions.black_hole_data()
		self.bhms,self.bhm,self.alpha,self.rg,self.rp,self.wp,self.X,self.Y,self.time = functions.parameters(self.bhml,self.bhmu,self.g,self.axm,self.astar,axms)
		
	 
	def output(self):
		
		self.rates,self.Z,self.ratesum = functions.superradiance_rates_detweiler(self.l,self.m,self.n,self.alpha,self.astar,self.axm,self.rp,self.X,self.Y,axms)
		#time_limit(time,rg,axm)
		self.x1,self.y1=functions.regge_contour_limits(self.X,self.Y,self.Z,self.time)
		self.nx,self.ny,self.indx,self.indx2=functions.regge_contour_outline(self.x1,self.y1)
		self.fx,self.fy,self.regiesx,self.regiesy=functions.regge_region(self.nx,self.ny,self.indx,self.indx2)
		self.ind=functions.black_hole_function_map(self.sr_masses,self.fx,self.fy)
		
	def stats(self):
		
		self.zscore,self.zupper,self.zlower,self.effsd,self.xtem,self.ytem,self.dytem,self.dxtem = functions.effective_zvalue(self.sr_spin_low, self.sr_mass_low, self.sr_spins, self.sr_masses, self.fx, self.fy)
		self.prob,self.exprob,self.upper_prob,self.lower_prob = functions.probability(self.zscore,self.zupper,self.zlower)
		self.nprob = functions.normalised_probability(self.prob,self.upper_prob,self.lower_prob)
	
		
		
		
	def print_out(self):
		ratesplot=False
		if ratesplot == True:
			plotting_functions.superradiance_rates_plot(self.alpha,self.rates,self.ratesum)
			plt.show()

		regge_zone = True
		if regge_zone == True:
			blackholes=True
			reggetrajectories=True
			plotting_functions.regge_region_plot(self.fx,self.fy,blackholes,reggetrajectories,self.xtem,self.ytem,self.dytem,self.dxtem,self.regiesx,self.regiesy)
			plt.show()

		regge_final = True
		if regge_final == True:
			plotting_functions.regge_plane_plot(self.x1,self.y1,self.sr_spins,self.sr_masses,self.sr_spin_up,self.sr_spin_low,self.sr_mass_up,self.sr_mass_low)
			plotting_functions.quantum_levels_legend(self.l)
			plotting_functions.conf_legend()
			plt.show()			

def main():
	my_matrix = axiverse_parameters()
	my_matrix.diagonalisation()
	
	my_calculator = superradiance_calculator()
	my_calculator.output()
	my_calculator.stats()
	my_calculator.print_out()

if __name__ == "__main__":
	main()

#########################################################################################################
#########################################################################################################


