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
from amuse.couple import bridge
from amuse.ext.protodisk import ProtoPlanetaryDisk

np.random.seed(42)

def disk_in_cluster(N_stars, N_gasParticles, W0):

	#set up protoplanetary disk
	Mdisk = 1.|units.MSun
	Rdisk = 100.|units.AU # HD 106906 b is at 700 AU
	Rmin, Rmax = 1.|units.AU, Rdisk
	converter_disk = nbody_system.nbody_to_si(Mdisk, Rdisk)

	disk_particles = ProtoPlanetaryDisk(N_gasParticles, convert_nbody=converter_disk, densitypower=1.5, Rmin=Rmin.value_in(units.AU), Rmax=Rmax.value_in(units.AU), q_out=1.0, discfraction=1.0).result

	hydro = Fi(converter_disk, mode='openmp')
	hydro.gas_particles.add_particles(disk_particles)

	hydro_to_framework = hydro.gas_particles.new_channel_to(disk_particles)
	framework_to_hydro = disk_particles.new_channel_to(hydro.gas_particles)

	#set up star cluster
	Mmin, Mmax = 0.1|units.MSun, 100.|units.MSun
	Rvir = 1.|units.parsec

	masses = new_kroupa_mass_distribution(N_stars, Mmax)
	Mtot_init = masses.sum()
	converter_cluster = nbody_system.nbody_to_si(Mtot_init, Rvir)
	stars = new_king_model(N_stars, W0, convert_nbody=converter_cluster)
	stars.mass = masses

	dist_diff = 100. #|units.parsec
	magic_index = -1

	#put the disk around star closest to 1 solar mass
	for i, star in enumerate(stars):

		star_x = star.position.value_in(units.parsec)[0]
		star_y = star.position.value_in(units.parsec)[1]
		star_z = star.position.value_in(units.parsec)[2]

		star_dist = np.sqrt(star_x**2 + star_y**2 + star_z**2)

		i_diff = np.abs(star_dist - 1.0)

		if i_diff < dist_diff:

			magic_index = i
			dist_diff = i_diff

	for gas_particle in disk_particles:

		gas_particle.position += stars[magic_index].position

	gravitating_bodies = ParticlesSuperset([stars, disk_particles])

	gravity = Hermite(converter_cluster)
	gravity.particles.add_particles(gravitating_bodies)

	gravity_to_framework = gravity.particles.new_channel_to(gravitating_bodies)
	framework_to_gravity = gravitating_bodies.new_channel_to(gravity.particles)

	combined = bridge.Bridge()
	combined.add_system(gravity, (hydro,))
	combined.add_system(hydro, (gravity,))

	tend, dt = 10.|units.Myr, 0.005|units.Myr

	sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
	sim_times = [ t|units.Myr for t in sim_times_unitless ]

	plt.rc('font', family = 'serif')
	cmap = cm.get_cmap('nipy_spectral')

	rvals, zvals = [], []

	t0 = time.time()

	for j, t in enumerate(sim_times):

		if j%1 == 0:
			print('')	
			print('simulation time: ', t)
			print('wall time: %.03f minutes'%((time.time()-t0)/60.))

		stars = gravity.particles[:N_stars]
		gas = gravity.particles[N_stars:]

		stars_x = [ stars.position[i].value_in(units.parsec)[0] for i in range(len(stars)) ]
		stars_y = [ stars.position[i].value_in(units.parsec)[1] for i in range(len(stars)) ]

		gas_x = [ gas.position[i].value_in(units.parsec)[0] for i in range(len(gas)) ]
		gas_y = [ gas.position[i].value_in(units.parsec)[1] for i in range(len(gas)) ]

		plt.figure()
		
		colors = np.log10(stars.mass.value_in(units.MSun))

		sc = plt.scatter(stars_x, stars_y, c=colors, s=1, alpha=0.5, label='Stars')
		gas_sc = plt.scatter(gas_x, gas_y, c='r', s=1, alpha=0.8, label='Protoplanetary Disk')
		cbar = plt.colorbar(sc)
		cbar.set_label(r'$\log_{10}(M_{\star}/M_{\odot})$', fontsize=14)

		plt.annotate(r'$M_{\mathrm{cluster}} = %.03f \, M_{\odot}$'%(Mtot_init.value_in(units.MSun)), xy=(0.65, 0.15), xycoords='axes fraction', fontsize=8)
		plt.annotate(r'$t = %.02f$ Myr'%(t.value_in(units.Myr)), xy=(0.65, 0.1), xycoords='axes fraction', fontsize=8)

		plt.xlim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.ylim(-4.*Rvir.value_in(units.parsec), 4.*Rvir.value_in(units.parsec))
		plt.gca().set_aspect('equal')
		plt.xlabel('$x$ (pc)', fontsize=14)
		plt.ylabel('$y$ (pc)', fontsize=14)
		plt.legend(loc='upper left', fontsize=8)
		plt.title(r'King Model, $W_{0} = %.01f$'%(W0), fontsize=14)

		plt.savefig('kingmodel_%s.png'%(str(j).rjust(5, '0')))
		plt.close()

		combined.evolve_model(t)

	combined.stop()

	return 0

    
if __name__ in ("__main__"):
    
	N_stars, N_gasParticles = 500, 200
	W0 = 7.0
	disk_in_cluster(N_stars, N_gasParticles, W0)
