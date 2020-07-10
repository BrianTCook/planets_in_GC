#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:39:44 2020

@author: BrianTCook
"""

import os

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

#partially from gravity_minimal.py in AMUSE textbook
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from amuse.lab import *

def gravity_minimal(N, W0):

	Mmin, Mmax = 0.1|units.MSun, 100.|units.MSun
	Rvir = 10.|units.parsec

	masses = new_kroupa_mass_distribution(N, Mmax)
	Mtot_init = masses.sum()
	converter = nbody_system.nbody_to_si(Mtot_init, Rvir)
	bodies = new_king_model(N, W0, convert_nbody=converter)
	bodies.mass = masses

	gravity = Hermite(converter)
	gravity.particles.add_particles(bodies)

	tend, dt = 100.|units.Myr, 0.5|units.Myr

	sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
	sim_times = [ t|units.Myr for t in sim_times_unitless ]

	plt.rc('font', family = 'serif')
	cmap = cm.get_cmap('nipy_spectral')

	rvals, zvals = [], []

	for j, t in enumerate(sim_times):

		if j%5 == 0:
			print('t is ', t)

		xs = [ gravity.particles.position[i].value_in(units.parsec)[0] for i in range(len(gravity.particles)) ]
		ys = [ gravity.particles.position[i].value_in(units.parsec)[1] for i in range(len(gravity.particles)) ]

		plt.figure()
		
		colors = np.log10(gravity.particles.mass.value_in(units.MSun))

		sc = plt.scatter(xs, ys, c=colors, s=1, alpha=0.8)
		cbar = plt.colorbar(sc)
		cbar.set_label(r'$\log_{10}(M_{\star}/M_{\odot})$', fontsize=14)

		plt.annotate(r'$t = %.02f$ Myr'%(t.value_in(units.Myr)), xy=(0.75, 0.1), xycoords='axes fraction', fontsize=8)

		plt.xlim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.ylim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.gca().set_aspect('equal')
		plt.xlabel('$x$ (pc)', fontsize=14)
		plt.ylabel('$y$ (pc)', fontsize=14)
		plt.title(r'King Model, $W_{0} = 7$', fontsize=14)

		plt.savefig('kingmodel_%s.png'%(str(j).rjust(5, '0')))
		plt.close()

		gravity.evolve_model(t)

	gravity.stop()

	return 0

    
if __name__ in ("__main__"):
    
	N = 1000
	W0 = 7.0
	gravity_minimal(N, W0)
