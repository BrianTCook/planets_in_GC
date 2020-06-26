#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:13:51 2020

@author: BrianTCook
"""

from __future__ import division

from amuse.lab import *
from amuse.couple import bridge
from amuse.support import io

from galpy.potential import MWPotential2014, KuzminKutuzovStaeckelPotential, to_amuse
from galpy.util import bovy_conversion

from gravity_code import gravity_code_setup
from fw_entropy import get_entropy
from cluster_maker import sort_clusters_by_attribute

import numpy as np
import pandas as pd
import gzip
import time
import math

def print_diagnostics(time, simulation_bodies, E_dyn, dE_dyn):
    
    print('------------')
    print('time: ', time)
    print('simulation_bodies.center_of_mass(): ', simulation_bodies.center_of_mass().value_in(units.kpc))
    print('E_dyn: ', E_dyn)
    print('dE_dyn: ', dE_dyn)
    print('------------')

def simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, sepBinary, 
               rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, radii, Norbiters, tend, dt):
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, 5.|units.parsec) #masses list is in solar mass units
    
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)
    
    #first thing is all particles in one superset, gravity is the code,
    #third thing is the list of orbiter bodies s.t. we can compute COMs independently
    #and plot them with different colors
    
    simulation_bodies, gravity, orbiter_bodies_list, cluster_colors, stellar = gravity_code_setup(code_name, orbiter_name, Mgalaxy, Rgalaxy, galaxy_code, sepBinary, 
                                                                                                  rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, radii, Norbiters)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(simulation_bodies)
    channel_from_framework_to_gravity = simulation_bodies.new_channel_to(gravity.particles)
    
    channel_from_stellar_to_framework = stellar.particles.new_channel_to(simulation_bodies)
    
    Ntotal = len(simulation_bodies)
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    #np.savetxt('times_in_Myr_%s_%s_Norbiters_%i.txt'%(code_name, orbiter_name, Norbiters), sim_times_unitless)
    
    delta_energies, clock_times = [], []
    
    body_masses = gravity.particles.mass

    '''
    #for 3D numpy array storage
    Nsavetimes = 50
    all_data = np.zeros((Nsavetimes+1, Ntotal, 6))
    mass_data = np.zeros((Nsavetimes+1, Ntotal))    
    '''

    #for saving in write_set_to_file
    filename = 'data_temp.csv'
    attributes = ('mass', 'x', 'y', 'z', 'vx', 'vy', 'vz')
    
    print('len(sim_times) is', len(sim_times))
    #gadget_flag = int(math.floor(len(sim_times)/Nsavetimes))
    
    t0 = time.time()
    #j_like_index = 0
    
    for j, t in enumerate(sim_times):

        clock_times.append(time.time()-t0) #will be in seconds
    
        '''
        if j == 0:
            E_dyn_init = gravity.kinetic_energy + gravity.potential_energy
            
        E_dyn = gravity.kinetic_energy + gravity.potential_energy
        dE_dyn = (E_dyn/E_dyn_init) - 1.
        
        delta_energies.append(dE_dyn)
         
        if j%gadget_flag == 0:
            
            print_diagnostics(t, simulation_bodies, E_dyn, dE_dyn)
             
            io.write_set_to_file(gravity.particles, filename, 'csv',
                                 attribute_types = (units.MSun, units.kpc, units.kpc, units.kpc, units.kms, units.kms, units.kms),
                                 attribute_names = attributes)
            
            data_t = pd.read_csv(filename, names=list(attributes))
            data_t = data_t.drop([0, 1, 2]) #removes labels units, and unit names
            
            masses = data_t['mass'].tolist()
            mass_data[j_like_index, :len(data_t.index)] = masses #in solar masses
        
            j_like_index += 1 
        '''   
        
        io.write_set_to_file(gravity.particles, filename, 'csv',
                             attribute_types = (units.MSun, units.kpc, units.kpc, units.kpc, units.kms, units.kms, units.kms),
                             attribute_names = attributes)
            
        data_t = pd.read_csv(filename, names=list(attributes))
        
    
        stellar.evolve_model(t)
        channel_from_stellar_to_framework.copy()
        
        channel_from_framework_to_gravity.copy()
        gravity.evolve_model(t)
        channel_from_gravity_to_framework.copy()

  
    np.savetxt(code_name + '_' + orbiter_name + '_masses_Norbiters_' + str(Norbiters) + '_dt_' + str(dt.value_in(units.Myr)) + '.txt', mass_data)

    return 0
