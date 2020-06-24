#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:57:38 2020

@author: BrianTCook
"""

from gc_initialization import globularCluster
from gravstel_code import code
from simulation_script import print_diagnostics, simulation

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:05:09 2020

@author: BrianTCook
"""

import os
import sys
import time

sys.path.append(os.getcwd())

from amuse.lab import *
from amuse.couple import bridge

import random
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#other scripts
from simulation_script import simulation

from galpy.df import quasiisothermaldf
from galpy.potential import MWPotential2014, KuzminKutuzovStaeckelPotential, to_amuse
from galpy.util import bovy_conversion
from galpy.actionAngle import actionAngleStaeckel

random.seed(73)

if __name__ in '__main__':
    
    potential = KuzminKutuzovStaeckelPotential #galpy
    
    tend, dt = 100.|units.Myr, 0.2|units.Myr
    
    #uses a galpy function to evaluate the enclosed mass
    Mgalaxy, Rgalaxy = float(6.8e10)|units.MSun, 2.6|units.kpc #disk mass for MWPotential2014, Bovy(2015)
    
    data_directory = '/home/brian/Desktop/second_project_gcs/data/' #'/home/s1780638/second_project_gcs/data/'
    
    rvals_all = np.loadtxt(data_directory+'ICs/dehnen_rvals.txt')
    phivals_all = np.loadtxt(data_directory+'ICs/dehnen_phivals.txt')
    zvals_all = np.loadtxt(data_directory+'ICs/dehnen_zvals.txt')
    
    vrvals_all = np.loadtxt(data_directory+'ICs/bovy_vrvals.txt')
    vphivals_all = np.loadtxt(data_directory+'ICs/bovy_vphivals.txt')
    vzvals_all = np.loadtxt(data_directory+'ICs/bovy_vzvals.txt')
    
    masses_all = np.loadtxt(data_directory+'ICs/cluster_masses_for_sampling.txt')
    radii_all = np.loadtxt(data_directory+'ICs/cluster_radii_for_sampling.txt')

    #orbiter_name = 'SingleCluster'
    code_name = 'Nbody'

    t0 = time.time()
                
    print('\\\\\\\\\\\\\\\\\\\\\\\\')
    print(code_name, orbiter_name)
    print('\\\\\\\\\\\\\\\\\\\\\\\\')
    
    t_init = time.time()

    simulation(code_name, orbiter_name, potential, Mgalaxy, Rgalaxy, 
               sepBinary, rvals_all, phivals_all, zvals_all, vrvals_all, vphivals_all, vzvals_all, 
               masses_all, radii_all, Norbiters, tend, dt)
    
    print('time is: %.03f minutes'%((time.time()-t0)/60.))
    print('time to run last simulation: %.03f minutes'%((time.time()-t_init)/60.))

sys.exit()
