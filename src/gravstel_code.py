#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:04:23 2020

@author: BrianTCook
"""

from amuse.lab import *
#from amuse.ext.bridge import bridge
from amuse.couple import bridge

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, KuzminKutuzovStaeckelPotential, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

from gc_initialization import globularCluster

import numpy as np
import time

def gravity_code_setup(code_name, orbiter_name, Mgalaxy, Rgalaxy, galaxy_code, sepBinary, 
                       rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, radii, Norbiters):
    
    '''
    will need to ask SPZ if he meant for field, orbiter to be separate in non
    Nemesis gravity solvers?
    '''
    
    gravity = bridge.Bridge(use_threading=False)
    
    converter_parent = nbody_system.nbody_to_si(Mgalaxy, Rgalaxy)
    converter_sub = nbody_system.nbody_to_si(np.median(masses)|units.MSun, np.median(radii)|units.parsec) #masses list is in solar mass units
    
    list_of_orbiters = [ orbiter(code_name, orbiter_name, Mgalaxy, Rgalaxy, sepBinary,
                                     rvals, phivals, zvals, vrvals, vphivals, vzvals, masses, i) for i in range(Norbiters) ]
    
    orbiter_bodies_list = [ list_of_orbiters[i][0] for i in range(Norbiters) ] 
    orbiter_codes_list = [ list_of_orbiters[i][1] for i in range(Norbiters) ]
    
    print(len(list_of_orbiters))
    
    cluster_colors = []
    
    for i, orbiter_code in enumerate(orbiter_codes_list):   

        stars = orbiter_code.particles.copy()
        cluster_color = np.random.rand(3,)
        
        cluster_colors.append([cluster_color]*len(stars))

    channel = stars.new_channel_to(orbiter_code.particles)
    channel.copy_attributes(['mass', 'x','y','z','vx','vy','vz'])

    stellar = SeBa()
    stellar.particles.add_particles(cluster_code.particles)

    #bridges each cluster with the bulge, not the other way around though
    gravity.add_system(cluster_code, other_things)  
    

    return gravity.particles, gravity, orbiter_bodies_list, cluster_colors, stellar
