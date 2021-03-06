#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 17:02:05 2020

takes globular cluster coordinates given in Harris (1996) and converts from Cartesian to cylindrical coordinates

@author: BrianTCook
"""

import numpy as np
import pandas as pd

def harris_coord_converter():
    
    '''
    harris_array loads in a .txt file
    '''
    
    #data_directory = '/home/brian/Desktop/second_project_gcs/data/'
    data_directory = '/Users/BrianTCook/Desktop/Thesis/second_project_gcs/data/'
    harris_df = pd.read_csv(data_directory+'harris_GC_data.csv', header=None)
    harris_df.columns = [ 'ID', 'Name', '?', 'RA', 'DEC', 'L', 'B', 'R_Sun', 'R_gc', 'X', 'Y', 'Z' ]
    
    xvals, yvals, zvals = harris_df['X'].tolist(), harris_df['Y'].tolist(), harris_df['Z'].tolist()
    
    rvals = [ np.sqrt(xvals[i]**2 + yvals[i]**2) for i in range(len(harris_df.index)) ]
    phivals = [ np.arctan(yvals[i]/xvals[i]) for i in range(len(harris_df.index)) ]
    
    rvals = [ float(r) for r in rvals ]
    phivals = [ float(phi) for phi in phivals ]
    zvals = [ float(z) for z in zvals ]
    
    np.savetxt('MW_GC_rvals.txt', rvals)
    np.savetxt('MW_GC_phivals.txt', phivals)
    np.savetxt('MW_GC_zvals.txt', zvals)
    
    return 1

if __name__ in '__main__':
    harris_coord_converter()
    ahh = np.loadtxt('MW_GC_zvals.txt')
    print(ahh[0])
    print('hello world!')

