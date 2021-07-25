## load libraries
from netCDF4 import Dataset
import numpy as np
from scipy import integrate


#####################################################################################
################################ Define OHC function  ###############################
#####################################################################################
def OHC_compute(input_dir, output_dir, sim_type, yr_start, yr_end): 
    fid = Dataset(input_dir+'3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
    levs = fid.variables['lev'][:]
    lats = fid.variables['lat'][:]
    lons = fid.variables['lon'][:]
    lon_2d, lat_2d = np.meshgrid(lons, lats).copy()
    lon_2d_0_360 = lon_2d.copy() 
    lon_2d_0_360[lon_2d_0_360<0] = lon_2d_0_360[lon_2d_0_360<0] + 360.  
    rho = fid.variables['rho'][:]
    rho[rho.mask==True] = np.nan
    rho_5m = rho[0,:,:].copy() 
#    rho_int_300m = integrate.simps(rho[levs<=300,:,:], levs[levs<=300], axis=0)/298.30493
#    rho_int_1000m = integrate.simps(rho[levs<=1008,:,:], levs[levs<=1008], axis=0)/1.0072571e+03
#    rho_int_2000m = integrate.simps(rho[levs<=2049,:,:], levs[levs<=2049], axis=0)/2.0488286e+03
    
    net_surf_heat = fid.variables['net_surf_heat'][:] # W/m2
    net_surf_heat[net_surf_heat.mask==True] = np.nan
    temp_yflux_adv_int_z = fid.variables['temp_yflux_adv_int_z'][:] # W
    temp_yflux_adv_int_z[temp_yflux_adv_int_z.mask==True] = np.nan
    wt = fid.variables['wt'][:] # m/s
    wt[wt.mask==True] = np.nan
    
    ######## compute OHC ######## 
    rho = fid.variables['rho'][:]
    rho[rho.mask==True] = np.nan
    temp = fid.variables['temp'][:] #  oC
    temp[temp.mask==True] = np.nan
    v = fid.variables['v'][:] #  m/s
    v[v.mask==True] = np.nan
    salt_prc = fid.variables['salt'][:] # psu
    z_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lon_2d)[0], np.shape(lon_2d)[1]))
    lat_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lon_2d)[0], np.shape(lon_2d)[1]))
    lon_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lon_2d)[0], np.shape(lon_2d)[1]))
    for i in np.arange(0, np.shape(levs)[0]):
        z_3d[i,:,:] = np.ones((np.shape(lon_2d)[0], np.shape(lon_2d)[1])) * (-levs[i])
        lat_3d[i,:,:] = lat_2d
        lon_3d[i,:,:] = lon_2d_0_360
    p_3d = gsw.p_from_z(z_3d,lat_3d)
    salt_abs_3d = gsw.SA_from_SP(salt_prc, p_3d, lon_3d, lat_3d)
    salt_abs_3d[salt_abs_3d.mask == True] = np.nan
    cp_3d = gsw.cp_t_exact(salt_abs_3d,temp,p_3d)
    cp_3d[cp_3d.mask==True] = np.nan
    ohc_lt100 = lon_2d.copy() * np.nan
    ohc_lt300 = lon_2d.copy() * np.nan
    ohc_lt700 = lon_2d.copy() * np.nan
    ohc_lt1000 = lon_2d.copy() * np.nan
    ohc_lt1500 = lon_2d.copy() * np.nan
    ohc_lt2000 = lon_2d.copy() * np.nan
    ohc_gt2000 = lon_2d.copy() * np.nan
    ohc_700_2000 = lon_2d.copy() * np.nan
    ohc_all = lon_2d.copy() * np.nan
    ohtv_lt100 = lon_2d.copy() * np.nan
    ohtv_lt300 = lon_2d.copy() * np.nan
    ohtv_lt700 = lon_2d.copy() * np.nan
    ohtv_lt1000 = lon_2d.copy() * np.nan
    ohtv_lt1500 = lon_2d.copy() * np.nan
    ohtv_lt2000 = lon_2d.copy() * np.nan
    ohtv_gt2000 = lon_2d.copy() * np.nan
    ohtv_700_2000 = lon_2d.copy() * np.nan
    ohtv_all = lon_2d.copy() * np.nan
    for i in np.arange(0,np.shape(ohc_all)[0],1):
        print(i, np.shape(ohc_all)[0])
        for j in np.arange(0,np.shape(ohc_all)[1],1):
            tmp = cp_3d[:,i,j]*rho[:,i,j]*(temp[:,i,j]+273.15)
            tmp_lt100 = tmp[levs<=106]
            levs_lt100 = levs[levs<=106]
            tmp_lt300 = tmp[levs<=300]
            levs_lt300 = levs[levs<=300]
            tmp_lt700 = tmp[levs<=729]
            levs_lt700 = levs[levs<=729]
            tmp_lt1000 = tmp[levs<=1008]
            levs_lt1000 = levs[levs<=1008]
            tmp_lt1500 = tmp[levs<=1573]
            levs_lt1500 = levs[levs<=1573]
            tmp_lt2000 = tmp[levs<=2049]
            levs_lt2000 = levs[levs<=2049]
            if  np.size(tmp[~np.isnan(tmp)]) > 1: 
                ohc_all[i,j] = integrate.simps(tmp[~np.isnan(tmp)],levs[~np.isnan(tmp)],axis = 0)/10.**9
            if  np.size(tmp_lt100[~np.isnan(tmp_lt100)]) > 1: 
                ohc_lt100[i,j] = integrate.simps(tmp_lt100[~np.isnan(tmp_lt100)],levs_lt100[~np.isnan(tmp_lt100)],axis = 0)/10.**9
            if  np.size(tmp_lt300[~np.isnan(tmp_lt300)]) > 1: 
                ohc_lt300[i,j] = integrate.simps(tmp_lt300[~np.isnan(tmp_lt300)],levs_lt300[~np.isnan(tmp_lt300)],axis = 0)/10.**9
            if  np.size(tmp_lt700[~np.isnan(tmp_lt700)]) > 1: 
                ohc_lt700[i,j] = integrate.simps(tmp_lt700[~np.isnan(tmp_lt700)],levs_lt700[~np.isnan(tmp_lt700)],axis = 0)/10.**9
            if  np.size(tmp_lt1000[~np.isnan(tmp_lt1000)]) > 1: 
                ohc_lt1000[i,j] = integrate.simps(tmp_lt1000[~np.isnan(tmp_lt1000)],levs_lt1000[~np.isnan(tmp_lt1000)],axis = 0)/10.**9
            if  np.size(tmp_lt1500[~np.isnan(tmp_lt1500)]) > 1: 
                ohc_lt1500[i,j] = integrate.simps(tmp_lt1500[~np.isnan(tmp_lt1500)],levs_lt1500[~np.isnan(tmp_lt1500)],axis = 0)/10.**9
            if  np.size(tmp_lt2000[~np.isnan(tmp_lt2000)]) > 1: 
                ohc_lt2000[i,j] = integrate.simps(tmp_lt2000[~np.isnan(tmp_lt2000)],levs_lt2000[~np.isnan(tmp_lt2000)],axis = 0)/10.**9
            del tmp
            # oceanic heat transport
            tmp = cp_3d[:,i,j]*rho[:,i,j]*(temp[:,i,j]+273.15) * v[:,i,j]
            tmp_lt100 = tmp[levs<=106]
            levs_lt100 = levs[levs<=106]
            tmp_lt300 = tmp[levs<=300]
            levs_lt300 = levs[levs<=300]
            tmp_lt700 = tmp[levs<=729]
            levs_lt700 = levs[levs<=729]
            tmp_lt1000 = tmp[levs<=1008]
            levs_lt1000 = levs[levs<=1008]
            tmp_lt1500 = tmp[levs<=1573]
            levs_lt1500 = levs[levs<=1573]
            tmp_lt2000 = tmp[levs<=2049]
            levs_lt2000 = levs[levs<=2049]
            if  np.size(tmp[~np.isnan(tmp)]) > 1: 
                ohtv_all[i,j] = integrate.simps(tmp[~np.isnan(tmp)],levs[~np.isnan(tmp)],axis = 0)/10.**9
            if  np.size(tmp_lt100[~np.isnan(tmp_lt100)]) > 1: 
                ohtv_lt100[i,j] = integrate.simps(tmp_lt100[~np.isnan(tmp_lt100)],levs_lt100[~np.isnan(tmp_lt100)],axis = 0)/10.**9
            if  np.size(tmp_lt300[~np.isnan(tmp_lt300)]) > 1: 
                ohtv_lt300[i,j] = integrate.simps(tmp_lt300[~np.isnan(tmp_lt300)],levs_lt300[~np.isnan(tmp_lt300)],axis = 0)/10.**9
            if  np.size(tmp_lt700[~np.isnan(tmp_lt700)]) > 1: 
                ohtv_lt700[i,j] = integrate.simps(tmp_lt700[~np.isnan(tmp_lt700)],levs_lt700[~np.isnan(tmp_lt700)],axis = 0)/10.**9
            if  np.size(tmp_lt1000[~np.isnan(tmp_lt1000)]) > 1: 
                ohtv_lt1000[i,j] = integrate.simps(tmp_lt1000[~np.isnan(tmp_lt1000)],levs_lt1000[~np.isnan(tmp_lt1000)],axis = 0)/10.**9
            if  np.size(tmp_lt1500[~np.isnan(tmp_lt1500)]) > 1: 
                ohtv_lt1500[i,j] = integrate.simps(tmp_lt1500[~np.isnan(tmp_lt1500)],levs_lt1500[~np.isnan(tmp_lt1500)],axis = 0)/10.**9
            if  np.size(tmp_lt2000[~np.isnan(tmp_lt2000)]) > 1: 
                ohtv_lt2000[i,j] = integrate.simps(tmp_lt2000[~np.isnan(tmp_lt2000)],levs_lt2000[~np.isnan(tmp_lt2000)],axis = 0)/10.**9
    
    ## write data to a netcdf file
    ncout = Dataset(output_dir+'OHC_density_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc', 'w', format='NETCDF4')
    # define axis size
    ncout.createDimension('lat',  np.shape(lats)[0])
    ncout.createDimension('lon',  np.shape(lons)[0])
    ncout.createDimension('lev',  np.shape(levs)[0])
    # create latitude axis
    lat = ncout.createVariable('lat', 'f4', ('lat'))
    lat.standard_name = 'lat'
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    lat[:] = lats[:]
    # create longitude axis
    lon = ncout.createVariable('lon', 'f4', ('lon'))
    lon.standard_name = 'lon'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    lon[:] = lons[:]
    # create vertical level axis
    lev = ncout.createVariable('lev', 'f4', ('lev'))
    lev.standard_name = 'lev'
    lev.long_name = 'vertical depths from the ocean surface'
    lev.units = 'm'
    lev[:] = levs
    
    # create variables
    vout = ncout.createVariable('rho_5m', 'f4', ('lat', 'lon'))
    vout.long_name = '5-m sea density'
    vout.units = 'kg m-3'
    vout[:] = rho_5m[:]
#    # create variables
#    vout = ncout.createVariable('rho_int_300m', 'f4', ('lat', 'lon'))
#    vout.long_name = '0-300-m sea density'
#    vout.units = 'kg m-3'
#    vout[:] = rho_int_300m[:]
#    # create variables
#    vout = ncout.createVariable('rho_int_1000m', 'f4', ('lat', 'lon'))
#    vout.long_name = '0-1000-m sea density'
#    vout.units = 'kg m-3'
#    vout[:] = rho_int_1000m[:]
#    # create variables
#    vout = ncout.createVariable('rho_int_2000m', 'f4', ('lat', 'lon'))
#    vout.long_name = '0-2000-m sea density'
#    vout.units = 'kg m-3'
#    vout[:] = rho_int_2000m[:]
    # create variables
    vout = ncout.createVariable('ohc_all', 'f4', ('lat', 'lon'))
    vout.long_name = 'All-depth OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_all[:]
    # create variables
    vout = ncout.createVariable('ohc_lt100', 'f4', ('lat', 'lon'))
    vout.long_name = '0-100-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt100[:]
    # create variables
    vout = ncout.createVariable('ohc_lt300', 'f4', ('lat', 'lon'))
    vout.long_name = '0-300-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt300[:]
    # create variables
    vout = ncout.createVariable('ohc_lt700', 'f4', ('lat', 'lon'))
    vout.long_name = '0-700-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt700[:]
    # create variables
    vout = ncout.createVariable('ohc_lt1000', 'f4', ('lat', 'lon'))
    vout.long_name = '0-1000-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt1000[:]
    # create variables
    vout = ncout.createVariable('ohc_lt1500', 'f4', ('lat', 'lon'))
    vout.long_name = '0-1500-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt1500[:]
    # create variables
    vout = ncout.createVariable('ohc_lt2000', 'f4', ('lat', 'lon'))
    vout.long_name = '0-2000-m OHC'
    vout.units = 'GJ m-2'
    vout[:] = ohc_lt2000[:]
    
    # create variables
    vout = ncout.createVariable('ohtv_all', 'f4', ('lat', 'lon'))
    vout.long_name = 'All-depth OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_all[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt100', 'f4', ('lat', 'lon'))
    vout.long_name = '0-100-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt100[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt300', 'f4', ('lat', 'lon'))
    vout.long_name = '0-300-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt300[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt700', 'f4', ('lat', 'lon'))
    vout.long_name = '0-700-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt700[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt1000', 'f4', ('lat', 'lon'))
    vout.long_name = '0-1000-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt1000[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt1500', 'f4', ('lat', 'lon'))
    vout.long_name = '0-1500-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt1500[:]
    # create variables
    vout = ncout.createVariable('ohtv_lt2000', 'f4', ('lat', 'lon'))
    vout.long_name = '0-2000-m OHT along V direction'
    vout.units = '10^9 W m-1'
    vout[:] = ohtv_lt2000[:]
    
    # create variables
    vout = ncout.createVariable('net_surf_heat', 'f4', ('lat', 'lon'))
    vout.long_name = 'surface ocean heat flux coming through coupler and mass transfer'
    vout.units = 'W m-2'
    vout[:] = net_surf_heat[:]
    # create variables
    vout = ncout.createVariable('temp_yflux_adv_int_z', 'f4', ('lat', 'lon'))
    vout.long_name = 'z-integral of cp*rho*dxt*v*temp with units of Watts'
    vout.units = 'W'
    vout[:] = temp_yflux_adv_int_z[:]
    
    # create variables for u 
    vout = ncout.createVariable('wt', 'f4', ('lev', 'lat', 'lon'))
    vout.long_name = 'vertical velocity'
    vout.units = 'm/sec'
    vout[:] = wt[:]
    # close files
    ncout.close()

