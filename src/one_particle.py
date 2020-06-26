#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 13:32:07 2020

@author: BrianTCook
"""

from __future__ import division

from amuse.lab import *
from amuse.couple import bridge
from amuse.support import io

from galpy.potential import MWPotential2014, KuzminKutuzovStaeckelPotential, to_amuse
from galpy.util import bovy_conversion

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gzip
import time
import math

def simulation():
    
        Hercules = Particles(1)
    Hercules.mass = 6e5 | units.MSun
    
    RA, DEC = 250.422, 36.460 #degrees
    D = 7.1 #distance from Sun, kpc
    vlos = -244.49 # pm 0.43, line-of-sight velocity in km/s
    
    mu_RA = -106.6 #pm 1.6, proper motion in km/s
    mu_DEC = -87.2 #pm 1.6, proper motion in km/s
    
    #will need to fix this
    Hercules.position = (20., 0., 0.) | units.kpc
    Hercules.velocity = (0., 200., 0.) | units.kms
    
    gravity = bridge.Bridge(use_threading=False)
    
    gravity_Hercules = Brutus()
    gravity_Hercules.particles.add_particles(Hercules)
    
    potential = KuzminKutuzovStaeckelPotential
    galaxy_code = to_amuse(potential, t=0.0, tgalpy=0.0, reverse=False, ro=None, vo=None)

    gravity.add_system(gravity_Hercules, galaxy_code) 

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(Hercules)
    channel_from_framework_to_gravity = Hercules.new_channel_to(gravity.particles)

    tend, dt = 1000.|units.Myr, 5.|units.Myr
    
    sim_times_unitless = np.arange(0., tend.value_in(units.Myr)+dt.value_in(units.Myr), dt.value_in(units.Myr))
    sim_times = [ t|units.Myr for t in sim_times_unitless ]
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    rvals, zvals = [], []
    
    for j, t in enumerate(sim_times):
        
        x_gal = gravity.particles.position[0].value_in(units.kpc)
        y_gal = gravity.particles.position[0].value_in(units.kpc)
        z_gal = gravity.particles.position[0].value_in(units.kpc)
        
        r_gal = np.sqrt(x_gal**2 + y_gal**2)

        rvals.append(r_gal)
        zvals.append(z_gal)

        plt.figure()
        
        plt.plot(rvals, zvals, c='k', label='M13')
        plt.legend(loc='upper right')
        
        plt.xlim(0., 50.)
        plt.ylim(-20., 20.)
        plt.xlabel('$r_{\mathrm{gal}}$ (kpc)', fontsize=14)
        plt.ylabel('$z_{\mathrm{gal}}$ (kpc)', fontsize=14)

        plt.savefig('m13_trajectory_%s.png'%(str(j).rjust(5, '0')))
        plt.close()

        gravity.evolve_model(t)
        
    gravity.stop()

    return 0

if __name__ in '__main__':
    
    simulation()
    
    print('hello world!')