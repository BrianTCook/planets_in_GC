#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 10:39:44 2020

@author: BrianTCook
"""

#Circumvent a problem with using too many threads on OpenMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "yes"

#partially from gravity_minimal.py in AMUSE textbook
import time
from amuse.lab import Hermite
from amuse.lab import nbody_system
from amuse.lab import new_king_model

def gravity_minimal(N, W0):

    Mmin, Mmax = 0.1|units.MSun, 100.|units.MSun
    Rvir = 10.|units.parsec
    
    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter = nbody_system.nbody_to_si(Mtot_init, Rvir)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = Hermite()
    gravity.particles.add_particles(bodies)
    
    tend, dt = 10000.|units.yr, 100.|units.yr
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    rvals, zvals = [], []
    
    for j, t in enumerate(sim_times):
        
        xs = gravity.particles.position[0].value_in(units.parsec)
        ys = gravity.particles.position[1].value_in(units.parsec)

        plt.figure()
        
        plt.plot(xs, ys, c='k', label='M13')
        plt.legend(loc='upper right')

        plt.xlabel('$x$ (pc)', fontsize=14)
        plt.ylabel('$y$ (pc)', fontsize=14)

        plt.savefig('kingmodel_%s.png'%(str(j).rjust(5, '0')))
        plt.close()

        gravity.evolve_model(t)
        
    gravity.stop()

    
if __name__ in ("__main__"):
    
    N = 1000
    W0 = 7.0
    gravity_minimal(N, W0)