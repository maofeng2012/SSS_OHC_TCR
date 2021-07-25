## load libraries
from netCDF4 import Dataset
import numpy as np
from scipy import integrate

#####################################################################################
############################# Function of Computing AMOC ############################
#####################################################################################
def AMOC40N(input_dir, output_dir, sim_type, start_yr, end_yr): 
    fid = Dataset(input_dir+'3d_ocean_data_'+sim_type+'_'+np.str(start_yr)+'_'+np.str(end_yr)+'.nc')
    levs = fid.variables['lev'][:]
    lats = fid.variables['lat'][:]
    lons = fid.variables['lon'][:]
    lon_2d, lat_2d = np.meshgrid(lons, lats)
    lon_2d_0_360 = lon_2d.copy() 
    lon_2d_0_360[lon_2d_0_360<0] = lon_2d_0_360[lon_2d_0_360<0] + 360.  
    rho = fid.variables['rho'][:]
    rho[rho.mask==True] = np.nan
    rho_5m = rho[0,:,:] 
    rho_int_1000m = integrate.simps(rho[levs<=1000,:,:], levs[levs<=1000], axis=0)/858.4215
    rho_int_300m = integrate.simps(rho[levs<=300,:,:], levs[levs<=300], axis=0)/298.30493
    ######## compute AMOC strength at lat=40 ########
    tmp = fid.variables['v'][:, (lats==latc-0.5) | (lats==latc+0.5), (lons>=-80.5) & (lons<=0.5)]
    tmp[tmp.mask==True] = np.nan
    v = np.nanmean(tmp, axis = 1)
    del tmp
    lons_amoc = lons[(lons>=-80.5) & (lons<=0.5)] * 110*10**3 * np.cos(latc*np.pi/180) # degree to m
    amoc = np.zeros((np.shape(levs)[0],)) * np.nan
    for j in np.arange(0,np.shape(amoc)[0]-1):
        v_i = v[:j,:]
        levs_i = levs[:j]
        data = np.zeros((np.shape(lons_amoc)[0],)) * np.nan
        for i in np.arange(0,np.shape(lons_amoc)[0]):
            v_f = v_i[:,i][~np.isnan(v_i[:,i])]
            levs_f = levs_i[~np.isnan(v_i[:,i])]
            if  levs_f.size > 1:
                data[i] = integrate.simps(v_f, levs_f, axis=0)
        dataf = data[~np.isnan(data)]
        lons_amoc_f = lons_amoc[~np.isnan(data)]
        if  dataf.size > 1:
            amoc[j] = np.nansum(dataf) * (lons_amoc[1] - lons_amoc[0])
            # amoc[j] = integrate.simps(dataf, lons_amoc_f)             
    fid.close()
    return amoc
    np.savetxt(output_dir+'AMOC@40N_'+sim_type+'_'+str(start_yr)+'_'+str(end_yr)+'_pattern.txt', amoc, fmt='%2.3f', delimiter=' ')




