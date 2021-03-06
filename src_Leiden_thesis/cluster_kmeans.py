#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:17:37 2020

@author: BrianTCook
"""

import numpy as np
import time
import glob
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pylab as pl
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit
from cluster_table import sort_clusters_by_attribute
pd.options.mode.chained_assignment = None  # default='warn'

def sklearn_mapper(true_labels, kmeans_labels):
    
    input_output_dict = {}

    pairings = [ (x, y) for x, y in zip(true_labels, kmeans_labels) ]
    pairings_counted = [ [x,pairings.count(x)] for x in set(pairings) ]

    for counted_pairing in pairings_counted:
        
        pairing, count = counted_pairing
        first, second = pairing
        
        if first not in input_output_dict.keys():
            
            input_output_dict.update( {first: second} )
        
        if count > pairings.count((first, input_output_dict[first])):
            
            input_output_dict.pop(first)
            input_output_dict.update( {first: second} )
    
    return input_output_dict

def get_kmeans_result(snapshots, Norbiters, initial_masses):
    
    t0 = time.time()
    
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    datadir_AMUSE = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/Enbid-2.0/AMUSE_data/'
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_radii = np.loadtxt(datadir + '/ICs/cluster_radii_for_sampling.txt') 
    
    indices_dict, cluster_masses, cluster_dists = sort_clusters_by_attribute('|r|')
            
    cluster_populations_sorted = [ cluster_populations[indices_dict[i]] for i in range(Norbiters) ]        
    cluster_radii_sorted = [ cluster_radii[indices_dict[i]] for i in range(Norbiters) ]
    cluster_masses_sorted = [ cluster_masses[indices_dict[i]] for i in range(Norbiters) ]
    cluster_dists_sorted = [ cluster_dists[indices_dict[i]] for i in range(Norbiters) ]

    true_labels = []
    
    for i in range(Norbiters):
        true_labels += [ i for j in range(int(cluster_populations_sorted[i])) ]
    
    delta_max = 0.9
    all_deltas = [ [] for i in range(Norbiters) ]
    endpoints = [ '' for i in range(Norbiters) ]
    active_orbiters = [ i for i in range(Norbiters) ]
    
    star_masses = np.loadtxt(datadir+'tree_data/tree_SingleCluster_masses_Norbiters_64_dt_0.2.txt').T
    
    deltaMs = np.zeros(star_masses[0,:].shape)
    
    for i in range(1, len(deltaMs)):
        deltaMs[i] = np.mean(star_masses[:, i]) - np.mean(star_masses[:, i-1])

    deltaM = np.mean(deltaMs)

    for k, snapshot in enumerate(snapshots):
        
        print(snapshot)
        print('time is: %.03f minutes'%((time.time()-t0)/60.))
        
        if len(active_orbiters) == 0:
            
            return all_deltas, endpoints
        
        data_filename = glob.glob(datadir_AMUSE+'*_%s_Norbiters_%i.ascii'%(snapshot, Norbiters))
        data_3D = np.loadtxt(data_filename[0])[:, :3]
    
        df = pd.DataFrame(data_3D, columns=['x', 'y', 'z'])

        df['labels'] = true_labels
        
        star_masses_truncated = star_masses[:len(df.index), k]
        
        df['masses'] = star_masses_truncated

        for cluster_label in active_orbiters:
            
            #print('cluster_label: ', cluster_label)
            #print('current time: %.04f minutes'%((time.time()-t0)/60.))
            
            df_cluster = df.loc[df['labels'] == cluster_label]
            df_cluster['separation distance'] = '' #in parsecs
            m_cluster_init = cluster_masses_sorted[cluster_label]
                    
            #sort by distance from COM of cluster
            xc, yc, zc = np.median(df_cluster['x'].tolist()), np.median(df_cluster['y'].tolist()), np.median(df_cluster['z'].tolist())
        
            for i in df_cluster.index:
                
                dist_sq = (df_cluster.at[i, 'x']-xc)**2 + (df_cluster.at[i, 'y']-yc)**2 + (df_cluster.at[i, 'z']-zc)**2 
                dist_in_pc = np.sqrt(dist_sq) * 1000.
                df_cluster.loc[i, 'separation distance'] = dist_in_pc
            
            #praint('median separation distance, initial radius: %.03f pc, %.03f pc'%(np.median(df_cluster['separation distance'].tolist()), cluster_radii_sorted[cluster_label] ))
            
            df_cluster = df_cluster[ df_cluster['separation distance'] <= 2.*cluster_radii_sorted[cluster_label] ] #outside of original radius
            df_cluster = df_cluster.reset_index(drop=True)
            m_cluster_t = np.sum(df_cluster['masses'].tolist())
            nstars = len(df_cluster.index)
            
            '''                
            add column saying current label for plot showing N_in birth cluster,
            N_field, and N_adopted by other cluster?
            '''     

            if len(df_cluster.index) != 0:

                m_cluster_t = np.sum(df_cluster['masses'].tolist())
                eps = nstars * deltaM / m_cluster_t
                
                delta = 1. - m_cluster_t / m_cluster_init * (1 + eps)
                
                galactocentric_dist = np.sqrt( np.median(df_cluster['x'].tolist())**2. + np.median(df_cluster['y'].tolist())**2. + np.median(df_cluster['z'].tolist())**2. )
                endpoints[cluster_label] = 'final galactocentric distance: %.03f kpc'%(galactocentric_dist)
                
            else:
                
                delta = 1.
                
            all_deltas[cluster_label].append(delta)
            
            if delta >= delta_max:
                
                endpoints[cluster_label] = 'disruption time: %.00f Myr'%(k*2.)
                active_orbiters.remove(cluster_label) #removes orbiter label from active ones
                
                
    return all_deltas, endpoints, cluster_dists_sorted

def func_powerlaw(x, m, b):
    
    #log scaled fit
    
    return m*x + b

if __name__ in '__main__':

    snapshots = [ str(j*10).rjust(5, '0') for j in range(51) ]
    
    initial_masses = 0.
    
    '''
    Norbiters = 64
    all_deltas, endpoints, cluster_dists_sorted = get_kmeans_result(snapshots, Norbiters, initial_masses) #endpoints
    dt = 2. #Myr
    
    plt.rc('text', usetex = True)
    plt.rc('font', family = 'serif')
    
    colors = pl.cm.rainbow(np.linspace(0,1,Norbiters))
    '''
    
    xvals_powerlaw = np.loadtxt('xvals_powerlaw.txt')
    yvals_powerlaw = np.loadtxt('yvals_powerlaw.txt')
    
    xs = np.log10(xvals_powerlaw)
    ys = np.log10(yvals_powerlaw)
    
    popt, pcov = curve_fit(func_powerlaw, xs, ys)
    m_optimal, b_optimal = popt
    
    print('m: %.03f \pm %.03f'%(m_optimal, np.sqrt(pcov[0,0])))
    print('b: %.03f \pm %.03f'%(b_optimal, np.sqrt(pcov[1,1])))
    
    xs_fit = np.linspace(min(xvals_powerlaw), max(xvals_powerlaw), 100)
    ys_fit = [ 10**(b_optimal) * x**(m_optimal) for x in xs_fit ]
    
    xs_fit_log = np.linspace(min(xs), max(xs), 100)
    ys_fit_log = [ func_powerlaw(x, m_optimal, b_optimal) for x in xs_fit_log ]

    '''
    xvals_powerlaw, yvals_powerlaw = [], []
    
    for cluster, deltas in enumerate(all_deltas):
        
        print('cluster: %i, fate: %s'%(cluster, endpoints[cluster]))
        
        xvals_powerlaw.append(cluster_dists_sorted[cluster])
        yvals_powerlaw.append( np.mean( [ (deltas[i+1]-deltas[i])/dt for i in range(len(deltas)-1) ] ) )
    '''
        
    plt.scatter(xs, ys, c='r', s=4, label='Star Clusters')
    plt.plot(xs_fit_log, ys_fit_log, c='k', linewidth=1, label=r'Power law, $\alpha = %.03f \pm %.03f$'%(popt[0], np.sqrt(pcov[0,0])))

    #plt.gca().set_xscale('log')

    #np.savetxt('xvals_powerlaw.txt', xvals_powerlaw)
    #np.savetxt('yvals_powerlaw.txt', yvals_powerlaw)

    plt.legend(loc='lower left', fontsize=10)
    plt.xlabel(r'$\log \big(r_{0} \, (\mathrm{kpc})\big)$', fontsize=16)
    plt.ylabel(r'$\log \big(\langle \dot{\delta} \rangle \, (\mathrm{Myr}^{-1}) \big)$', fontsize=16)
    plt.gca().tick_params(labelsize='large')
    plt.tight_layout()
    plt.savefig('mass_loss_Norbiters_%i_powerlaw.pdf'%(Norbiters))    
    
    print('hello world!')
    
'''
success plotting

    hits_2D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_2D) ]
    hits_3D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_3D) ]
    hits_6D = [ 1 if y1 == y2 else 0 for y1, y2 in zip(true_labels, y_compare_6D) ]
    
    if hits_2D.count(1) > hits_2D.count(0):
        
        success_2D = np.sum(hits_2D)/len(hits_2D)
        
    if hits_2D.count(1) <= hits_2D.count(0):
        
        success_2D = 1. - np.sum(hits_2D)/len(hits_2D)
        
    if hits_3D.count(1) > hits_3D.count(0):
        
        success_3D = np.sum(hits_3D)/len(hits_3D)
        
    if hits_3D.count(1) <= hits_3D.count(0):
        
        success_3D = 1. - np.sum(hits_3D)/len(hits_3D)
        
    if hits_6D.count(1) > hits_6D.count(0):
        
        success_6D = np.sum(hits_6D)/len(hits_6D)
        
    if hits_6D.count(1) <= hits_6D.count(0):
        
        success_6D = 1. - np.sum(hits_6D)/len(hits_2D)


    if plot_2D == True:
        
        plt.figure()    
        plt.scatter(data_3D[:, 1], data_3D[:, 2], c=y_kmeans_2D, s=5, cmap='Dark2')
        plt.annotate('2D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_2D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_2D.jpg'%(Norbiters))
        plt.close()    
    
    if plot_3D == True:
        
        plt.figure()    
        plt.scatter(data_3D[:, 1], data_3D[:, 2], c=y_kmeans_3D, s=5, cmap='Dark2')
        plt.annotate('3D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_3D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_3D.jpg'%(Norbiters))
        plt.close()
        
    if plot_6D == True:
        
        plt.figure()    
        plt.scatter(data_6D[:, 1], data_6D[:, 2], c=y_kmeans_6D, s=5, cmap='Dark2')
        plt.annotate('6D coordinates', xy=(0.6, 0.2), xycoords='axes fraction', fontsize=12)
        plt.annotate('Success rate: %.05f'%(success_6D), xy=(0.6, 0.1), xycoords='axes fraction', fontsize=12)
        plt.xlabel(r'$y$ (kpc)', fontsize=20)
        plt.ylabel(r'$z$ (kpc)', fontsize=20)
        
        plt.xlim(-2, 2)
        plt.ylim(-1, 1)
        
        plt.title(r'$k$-means, $N_{clusters} = %i$'%(Norbiters), fontsize=16)
        plt.tight_layout()
        plt.savefig('kmeans_Norbiters_%i_6D.jpg'%(Norbiters))
        plt.close()
    
    return success_2D, success_3D, success_6D
   
for logN in range(logN_max+1):
    
    print('N = %i'%(2**logN))
    if plot_2D == True or plot_3D == True or plot_6D == True:
        
        plt.rc('text', usetex = True)
        plt.rc('font', family = 'serif')
    
    s_2, s_3, s_6 = get_kmeans_result(2**logN, plot_3D, plot_6D)
    s2_vals.append(s_2)
    s3_vals.append(s_3)
    s6_vals.append(s_6)
    
plt.figure()
plt.plot(logN_vals, s2_vals, label='2D distance metric')
plt.plot(logN_vals, s3_vals, label='3D distance metric')
plt.plot(logN_vals, s6_vals, label='6D distance metric')
plt.ylim(0, 1)
plt.legend(loc='best')
plt.title(r'$k$-means Cluster Identification', fontsize=16)
plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=14)
plt.ylabel('Labelling Accuracy', fontsize=14)
plt.tight_layout()
plt.savefig('accuracies_kmeans.pdf')
    
plt.figure()
plt.plot(logN_vals, s2_vals, label='2D distance metric')
plt.plot(logN_vals, s3_vals, label='3D distance metric')
plt.plot(logN_vals, s6_vals, label='6D distance metric')
plt.ylim(0, 1)
plt.legend(loc='best')
plt.title(r'$k$-means Cluster Identification', fontsize=16)
plt.xlabel(r'$\log_{2} N_{\mathrm{clusters}}$', fontsize=14)
plt.ylabel('Labelling Accuracy', fontsize=14)
plt.tight_layout()
plt.savefig('accuracies_kmeans.pdf')
'''
