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
from amuse.ext.protodisk import ProtoPlanetaryDisk

def disk_in_cluster(N_stars, N_gasParticles, W0):

	#set up protoplanetary disk
	Mdisk = 1.|units.MSun
	Rdisk = 500.|units.AU # HD 106906 b is at 700 AU
	Rmin, Rmax = 1., Rdisk.value_in(units.AU)
	converter_disk = nbody_system.nbody_to_si(Mdisk, Rdisk)

	disk_particles = ProtoPlanetaryDisk(N_gasParticles, convert_nbody=converter_disk, densitypower=1.5, Rmin=Rmin, Rmax=Rmax/Rmin, q_out=1.0, discfraction=1.0).result

	#hydro = Fi(converter_disk, mode='openmp')
	#hydro.gas_particles.add_particles(disk_particles)

	#hydro_to_framework = 
	#framework_to_hydro = 

	#set up star cluster
	Mmin, Mmax = 0.1|units.MSun, 100.|units.MSun
	Rvir = 10.|units.parsec

	masses = new_kroupa_mass_distribution(N_stars, Mmax)
	Mtot_init = masses.sum()
	converter_cluster = nbody_system.nbody_to_si(Mtot_init, Rvir)
	stars = new_king_model(N_stars, W0, convert_nbody=converter_cluster)
	stars.mass = masses

	mass_diff = 10. #|units.MSun
	magic_index = -1

	#put the disk around star closest to 1 solar mass
	for i, star in enumerate(stars):

		i_diff = np.abs(star.mass.value_in(units.MSun) - 1.)

		if i_diff < mass_diff:

			magic_index = i
			mass_diff = i_diff

	for gas_particle in disk_particles:

		gas_particle.position += stars[magic_index].position

	gravity = BHTree(converter_cluster)
	gravity.particles.add_particles(stars)
	gravity.particles.add_particles(disk_particles)

	#gravity_to_framework = gravity.particles.new_channel_to()
	#framework_to_gravity = 

	#combined = bridge.Bridge()
	#combined.add_system(gravity, (hydro,))
	#combined.add_system(hydro, (gravity,))

	tend, dt = 100.|units.Myr, 0.5|units.Myr

	sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
	sim_times = [ t|units.Myr for t in sim_times_unitless ]

	plt.rc('font', family = 'serif')
	cmap = cm.get_cmap('nipy_spectral')

	rvals, zvals = [], []

	for j, t in enumerate(sim_times):

		if j%5 == 0:
			print('t is ', t)

		stars = gravity.particles[:N_stars]
		gas = gravity.particles[N_stars:]

		stars_x = [ stars.position[i].value_in(units.parsec)[0] for i in range(len(stars)) ]
		stars_y = [ stars.position[i].value_in(units.parsec)[1] for i in range(len(stars)) ]

		gas_x = [ gas.position[i].value_in(units.parsec)[0] for i in range(len(gas)) ]
		gas_y = [ gas.position[i].value_in(units.parsec)[1] for i in range(len(gas)) ]

		plt.figure()
		
		colors = np.log10(stars.mass.value_in(units.MSun))

		sc = plt.scatter(stars_x, stars_y, c=colors, s=1, alpha=0.6)
		gas_sc = plt.scatter(gas_x, gas_y, c='r', s=0.1, alpha=0.9, label='Protoplanetary Disk')
		cbar = plt.colorbar(sc)
		cbar.set_label(r'$\log_{10}(M_{\star}/M_{\odot})$', fontsize=14)

		plt.annotate(r'$t = %.02f$ Myr'%(t.value_in(units.Myr)), xy=(0.75, 0.1), xycoords='axes fraction', fontsize=8)

		plt.xlim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.ylim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.gca().set_aspect('equal')
		plt.xlabel('$x$ (pc)', fontsize=14)
		plt.ylabel('$y$ (pc)', fontsize=14)
		plt.legend(loc='upper left', fontsize=8)
		plt.title(r'King Model, $W_{0} = 7$', fontsize=14)

		plt.savefig('kingmodel_%s.png'%(str(j).rjust(5, '0')))
		plt.close()

		gravity.evolve_model(t)

	gravity.stop()

	return 0

    
if __name__ in ("__main__"):
    
	N_stars, N_gasParticles = 5000, 500
	W0 = 7.0
	disk_in_cluster(N_stars, N_gasParticles, W0)
