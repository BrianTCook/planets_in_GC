import gzip
import numpy as np

def convert_numpy(code_names, orbiter_names, Norbiters_list):

    for code_name in code_names:
        for orbiter_name in orbiter_names:
            
            if code_name == 'nemesis' and orbiter_name == 'SingleStar':
                continue
            
            for Norbiters in Norbiters_list:
    
                f_all = gzip.GzipFile('all_data_%s_%s_Norbiters_%s.npy.gz'%(code_name, orbiter_name, str(Norbiters)), 'r')
                all_data = np.load(f_all)
                N_timesteps = len(all_data[:,0,0])
                
                for i in range(N_timesteps):
                    
                    data_to_keep = all_data[i, :, 1:] #gets rid of mass
                    np.savetxt('for_enbid_%s_%s_frame_%s_Norbiters_%s.ascii'%(code_name, orbiter_name,
                                                                              str(i).rjust(5, '0'), str(Norbiters)), data_to_keep)
                    
    return 0