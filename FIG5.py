## load libraries
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

#######################################################################################            
############################# Read Observation data for NCEI ##########################
#######################################################################################
yr_s = 1968
yr_e = 2017
yr = np.arange(1968,2018,1) # NCEI
yr_f = np.arange(yr_s,yr_e+1,1)
fid = Dataset('salt_temp_rho_alpha_beta_yr_depth_lat_ATL_PAC_IND_obs_Levitus.nc')
lats = fid.variables['lat'][:]
depth = fid.variables['depth'][:]
lats_2d_NCEI, depth_2d_NCEI = np.meshgrid(lats, depth)

temp_nonatl_depth_lat = fid.variables['temp_nonatl_depth_lat'][(yr>=yr_s) & (yr<=yr_e),:,:]
temp_atl_depth_lat = fid.variables['temp_atl_depth_lat'][(yr>=yr_s) & (yr<=yr_e),:,:]

salt_nonatl_depth_lat = fid.variables['salt_nonatl_depth_lat'][(yr>=yr_s) & (yr<=yr_e),:,:]
salt_atl_depth_lat = fid.variables['salt_atl_depth_lat'][(yr>=yr_s) & (yr<=yr_e),:,:]
fid.close()

slope_temp_all = np.nan * np.zeros((np.shape(lats_2d_NCEI)[0], np.shape(lats_2d_NCEI)[1]))
pval_temp_all = np.nan * np.zeros((np.shape(lats_2d_NCEI)[0], np.shape(lats_2d_NCEI)[1]))
slope_salt_all = np.nan * np.zeros((np.shape(lats_2d_NCEI)[0], np.shape(lats_2d_NCEI)[1]))
pval_salt_all = np.nan * np.zeros((np.shape(lats_2d_NCEI)[0], np.shape(lats_2d_NCEI)[1]))

slope_salt_nonatl = slope_salt_all.copy()
slope_salt_atl = slope_salt_all.copy()
pval_salt_nonatl = slope_salt_all.copy()
pval_salt_atl = slope_salt_all.copy()

slope_temp_nonatl = slope_salt_all.copy()
slope_temp_atl = slope_salt_all.copy()
pval_temp_nonatl = slope_salt_all.copy()
pval_temp_atl = slope_salt_all.copy()


for i in np.arange(0,np.shape(lats_2d_NCEI)[0],1):
    for j in np.arange(0,np.shape(lats_2d_NCEI)[1],1):
        # salt
        slope_salt_nonatl[i,j], intcept, rval, pval_salt_nonatl[i,j], std_err = stats.linregress(yr_f, salt_nonatl_depth_lat[:, i,j])
        slope_salt_atl[i,j], intcept, rval, pval_salt_atl[i,j], std_err = stats.linregress(yr_f, salt_atl_depth_lat[:, i,j])
        # temp
        slope_temp_nonatl[i,j], intcept, rval, pval_temp_nonatl[i,j], std_err = stats.linregress(yr_f, temp_nonatl_depth_lat[:, i,j])
        slope_temp_atl[i,j], intcept, rval, pval_temp_atl[i,j], std_err = stats.linregress(yr_f, temp_atl_depth_lat[:, i,j])

slope_temp_nonatl_NCEI = slope_temp_nonatl.copy()
slope_salt_nonatl_NCEI = slope_salt_nonatl.copy()
slope_temp_atl_NCEI = slope_temp_atl.copy()
slope_salt_atl_NCEI = slope_salt_atl.copy()

pval_temp_nonatl_NCEI = pval_temp_nonatl.copy()
pval_salt_nonatl_NCEI = pval_salt_nonatl.copy()
pval_temp_atl_NCEI = pval_temp_atl.copy()
pval_salt_atl_NCEI = pval_salt_atl.copy()


################################# read ocean mask #############################
yr_start = 101
yr_end = 200
sim_type = 'ctl'
fid = Dataset('3d_ocean_data_' + sim_type + '_' + np.str(yr_start) + '_'+np.str(yr_end) + '.nc')
levs = fid.variables['lev'][:]
fid.close()

fid = Dataset('ocean_index_with_ATL_PAC_IND_mod.nc')
mask_i = fid.variables['ocean_index'][:]
mask_2d = mask_i.copy()
lat_mask = fid.variables['lat'][:]
lon_mask = fid.variables['lon'][:]
fid.close()

mask_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lat_mask)[0], np.shape(lon_mask)[0]))
for i in np.arange(0,np.shape(mask_3d)[0]):
    mask_3d[i,:,:] = mask_i
del mask_i

#####################################################################################
################################### read CTL data ###################################
#####################################################################################
yr_start = 101
yr_end = 200
sim_type = 'ctl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
levs = fid.variables['lev'][:]
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]
lats_2d, levs_2d = np.meshgrid(lats, levs)
rho = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho[rho.mask==True] = np.nan
temp = fid.variables['temp'][:] #  oC
temp[temp.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()


integral_params = np.zeros((np.shape(temp)[0], np.shape(temp)[1], np.shape(temp)[2]))
lons_2d_intergral, lats_2d_intergral = np.meshgrid(lons, lats)
for i in np.arange(0, np.shape(integral_params)[0], 1):
    integral_params[i,:,:] = 110 * 1000 * np.cos(lats_2d_intergral * np.pi / 180.)


temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_ctl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_ctl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

salt_tmp = salt.copy()
salt_tmp[mask_3d == 1.0] = np.nan
salt_cs_nonatl_ctl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp
salt_tmp = salt.copy()
salt_tmp[mask_3d != 1.0] = np.nan
salt_cs_atl_ctl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp


sim_type = 'ctlgl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
levs = fid.variables['lev'][:]
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:] 
rho = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho[rho.mask==True] = np.nan
temp = fid.variables['temp'][:] #  oC
temp[temp.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()

temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_ctlgl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_ctlgl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

salt_tmp = salt.copy()
salt_tmp[mask_3d == 1.0] = np.nan
salt_cs_nonatl_ctlgl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp
salt_tmp = salt.copy()
salt_tmp[mask_3d != 1.0] = np.nan
salt_cs_atl_ctlgl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp


#####################################################################################
################################### read 2CO2 data ##################################
#####################################################################################
yr_start = 161
yr_end = 180
sim_type = '2CO2std'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
levs = fid.variables['lev'][:]
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]
rho = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho[rho.mask==True] = np.nan
temp = fid.variables['temp'][:] #  oC
temp[temp.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()

temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_2co2std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_2co2std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

salt_tmp = salt.copy()
salt_tmp[mask_3d == 1.0] = np.nan
salt_cs_nonatl_2co2std = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp
salt_tmp = salt.copy()
salt_tmp[mask_3d != 1.0] = np.nan
salt_cs_atl_2co2std = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp


yr_start = 161
yr_end = 180
sim_type = '2CO2gl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
levs = fid.variables['lev'][:]
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]
#lon_2d, lat_2d = np.meshgrid(lons, lats)
#lon_2d_0_360 = lon_2d.copy() 
#lon_2d_0_360[lon_2d_0_360<0] = lon_2d_0_360[lon_2d_0_360<0] + 360.  
rho = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho[rho.mask==True] = np.nan
temp = fid.variables['temp'][:] #  oC
temp[temp.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()

temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_2co2gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_2co2gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

salt_tmp = salt.copy()
salt_tmp[mask_3d == 1.0] = np.nan
salt_cs_nonatl_2co2gl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp
salt_tmp = salt.copy()
salt_tmp[mask_3d != 1.0] = np.nan
salt_cs_atl_2co2gl = np.nansum(salt_tmp * integral_params, axis = 2) / 10.**6
del salt_tmp


####################################################################################
################################### compute diff ###################################
#####################################################################################
############ temp 
chg_temp_cs_nonatl_2co2std = temp_cs_nonatl_2co2std - temp_cs_nonatl_ctl
chg_temp_cs_nonatl_2co2gl = temp_cs_nonatl_2co2gl - temp_cs_nonatl_ctlgl

chg_temp_cs_atl_2co2std = temp_cs_atl_2co2std - temp_cs_atl_ctl
chg_temp_cs_atl_2co2gl = temp_cs_atl_2co2gl - temp_cs_atl_ctlgl

############ salt
chg_salt_cs_nonatl_2co2std = salt_cs_nonatl_2co2std - salt_cs_nonatl_ctl
chg_salt_cs_nonatl_2co2gl = salt_cs_nonatl_2co2gl - salt_cs_nonatl_ctlgl

chg_salt_cs_atl_2co2std = salt_cs_atl_2co2std - salt_cs_atl_ctl
chg_salt_cs_atl_2co2gl = salt_cs_atl_2co2gl - salt_cs_atl_ctlgl
############################################################################################            
############################### read OHC data for all runs #################################
############################################################################################ 
lon_2d_i = np.loadtxt('lon_ocean.txt', delimiter = ' ')
lat_2d = np.loadtxt('lat_ocean.txt', delimiter = ' ')
lon_2d = np.concatenate((lon_2d_i[:,lon_2d_i[0,:]>30]-360., lon_2d_i[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()
lat_2d = np.concatenate((lat_2d[:,lon_2d_i[0,:]>30], lat_2d[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()   

fid = Dataset('OHC_density_ctlgl_101_200.nc')
ohc_ctlgl_all = fid.variables['ohc_all'][:]
ohc_ctlgl_lt300 = fid.variables['ohc_lt300'][:]
ohc_ctlgl_lt700 = fid.variables['ohc_lt700'][:]
ohc_ctlgl_lt2000 = fid.variables['ohc_lt2000'][:]
fid.close()

fid = Dataset('OHC_density_2CO2gl_161_180.nc')
ohc_co2gl_all = fid.variables['ohc_all'][:]
ohc_co2gl_lt300 = fid.variables['ohc_lt300'][:]
ohc_co2gl_lt700 = fid.variables['ohc_lt700'][:]
ohc_co2gl_lt2000 = fid.variables['ohc_lt2000'][:]
fid.close()

fid = Dataset('OHC_density_ctl_101_200.nc')
ohc_ctl_all = fid.variables['ohc_all'][:]
ohc_ctl_lt300 = fid.variables['ohc_lt300'][:]
ohc_ctl_lt700 = fid.variables['ohc_lt700'][:]
ohc_ctl_lt2000 = fid.variables['ohc_lt2000'][:]
fid.close()

fid = Dataset('OHC_density_2CO2std_161_180.nc')
ohc_co2std_all = fid.variables['ohc_all'][:]
ohc_co2std_lt300 = fid.variables['ohc_lt300'][:]
ohc_co2std_lt700 = fid.variables['ohc_lt700'][:]
ohc_co2std_lt2000 = fid.variables['ohc_lt2000'][:]
fid.close()

chg_ohc_std_all = ohc_co2std_all - ohc_ctl_all
chg_ohc_std_lt300 = ohc_co2std_lt300 - ohc_ctl_lt300
chg_ohc_std_lt700 = ohc_co2std_lt700 - ohc_ctl_lt700
chg_ohc_std_300_700 = (ohc_co2std_lt700 - ohc_co2std_lt300) - (ohc_ctl_lt700 - ohc_ctl_lt300)
chg_ohc_std_700_2000 = (ohc_co2std_lt2000 - ohc_co2std_lt700) - (ohc_ctl_lt2000 - ohc_ctl_lt700)
chg_ohc_std_gt2000 = (ohc_co2std_all - ohc_co2std_lt2000) - (ohc_ctl_all - ohc_ctl_lt2000)

chg_ohc_gl_all = ohc_co2gl_all - ohc_ctlgl_all
chg_ohc_gl_lt300 = ohc_co2gl_lt300 - ohc_ctlgl_lt300
chg_ohc_gl_lt700 = ohc_co2gl_lt700 - ohc_ctlgl_lt700
chg_ohc_gl_300_700 = (ohc_co2gl_lt700 - ohc_co2gl_lt300) - (ohc_ctlgl_lt700 - ohc_ctlgl_lt300)
chg_ohc_gl_700_2000 = (ohc_co2gl_lt2000 - ohc_co2gl_lt700) - (ohc_ctlgl_lt2000 - ohc_ctlgl_lt700)
chg_ohc_gl_gt2000 = (ohc_co2gl_all - ohc_co2gl_lt2000) - (ohc_ctlgl_all - ohc_ctlgl_lt2000)

diff_ohc_std_gl_all = chg_ohc_std_all - chg_ohc_gl_all
diff_ohc_std_gl_lt300 = chg_ohc_std_lt300 - chg_ohc_gl_lt300
diff_ohc_std_gl_300_700 = chg_ohc_std_300_700 - chg_ohc_gl_300_700
diff_ohc_std_gl_700_2000 = chg_ohc_std_700_2000 - chg_ohc_gl_700_2000
diff_ohc_std_gl_gt2000 = chg_ohc_std_gt2000 - chg_ohc_gl_gt2000


################## read ocean mask index ###############
fid = Dataset('/Users/jiexue/Documents/Projects/Salinity/Work/OBS_MODEL/Data/ocean_index_with_ATL_PAC_IND_mod.nc')
ocean_mask = fid.variables['ocean_index'][:]
lat_mask = fid.variables['lat'][:]
lon_mask = fid.variables['lon'][:]
fid.close()
lons_2d, lats_2d = np.meshgrid(lon_mask, lat_mask)

ocean_mask[lats_2d <= -33] = 3.
#ocean_mask[(lats_2d > 0) & (ocean_mask==2)] = 25. # north Pacific
#ocean_mask[(lats_2d < 0) & (ocean_mask==2)] = 30. # south Pacific

# compute area for each grid
lon = lons_2d.copy()
lat = lats_2d.copy()
lon_diff = np.zeros((np.shape(lon)[0], np.shape(lon)[1]))
lat_diff = np.zeros((np.shape(lat)[0], np.shape(lat)[1]))
for i in np.arange(1,np.shape(lon)[1]-1,1):
    lon_diff[:,i] = (lon[:,i+1] - lon[:,i-1]) / 2.
for i in np.arange(1,np.shape(lat)[0]-1,1):   
    lat_diff[i,:] = (lat[i+1,:] - lat[i-1,:]) / 2.

lon_diff[:,0] = lon[:,1] - lon[:,0]
lon_diff[:,-1] = lon[:,-1] - lon[:,-2]
lat_diff[0,:] = lat[1,:] - lat[0,:]
lat_diff[-1,:] = lat[-1,:] - lat[-2,:] 
cos_lat = np.cos(lat*np.pi/180.)
area = (110*10**3)**2 * lon_diff * lat_diff * cos_lat


temp = (diff_ohc_std_gl_all * area).copy()
mean_diff_ohc_std_gl_all = np.array([np.nansum(temp[ocean_mask==1]),np.nansum(temp[ocean_mask!=1]), np.nansum(temp[ocean_mask==2]), np.nansum(temp[ocean_mask==3]), np.nansum(temp[ocean_mask==4]), np.nansum(temp)]) / 10.**24 * 10.**9
temp = (diff_ohc_std_gl_lt300 * area).copy()
mean_diff_ohc_std_gl_lt300 = np.array([np.nansum(temp[ocean_mask==1]), np.nansum(temp[ocean_mask!=1]), np.nansum(temp[ocean_mask==2]), np.nansum(temp[ocean_mask==3]), np.nansum(temp[ocean_mask==4]), np.nansum(temp)]) / 10.**24 * 10.**9
temp = (diff_ohc_std_gl_300_700 * area).copy()
mean_diff_ohc_std_gl_300_700 = np.array([np.nansum(temp[ocean_mask==1]), np.nansum(temp[ocean_mask!=1]), np.nansum(temp[ocean_mask==2]), np.nansum(temp[ocean_mask==3]), np.nansum(temp[ocean_mask==4]), np.nansum(temp)]) / 10.**24 * 10.**9
temp = (diff_ohc_std_gl_700_2000 * area).copy()
mean_diff_ohc_std_gl_700_2000 = np.array([np.nansum(temp[ocean_mask==1]), np.nansum(temp[ocean_mask!=1]), np.nansum(temp[ocean_mask==2]), np.nansum(temp[ocean_mask==3]), np.nansum(temp[ocean_mask==4]), np.nansum(temp)]) / 10.**24 * 10.**9
temp = (diff_ohc_std_gl_gt2000 * area).copy()
mean_diff_ohc_std_gl_gt2000 = np.array([np.nansum(temp[ocean_mask==1]), np.nansum(temp[ocean_mask!=1]), np.nansum(temp[ocean_mask==2]), np.nansum(temp[ocean_mask==3]), np.nansum(temp[ocean_mask==4]), np.nansum(temp)]) / 10.**24 * 10.**9


OHC_all = np.array([mean_diff_ohc_std_gl_all[-1], mean_diff_ohc_std_gl_lt300[-1], \
                    mean_diff_ohc_std_gl_300_700[-1], mean_diff_ohc_std_gl_700_2000[-1], \
                    mean_diff_ohc_std_gl_gt2000[-1]])
OHC_atl = np.array([mean_diff_ohc_std_gl_all[0], mean_diff_ohc_std_gl_lt300[0], \
                    mean_diff_ohc_std_gl_300_700[0], mean_diff_ohc_std_gl_700_2000[0], \
                    mean_diff_ohc_std_gl_gt2000[0]])
OHC_nonatl = np.array([mean_diff_ohc_std_gl_all[1], mean_diff_ohc_std_gl_lt300[1], \
                    mean_diff_ohc_std_gl_300_700[1], mean_diff_ohc_std_gl_700_2000[1], \
                    mean_diff_ohc_std_gl_gt2000[1]])

#######################################################################################            
################### plot cross sections of temp and salt GL - STD #####################
#######################################################################################
def new_cmap(levels):
    n = np.int((np.shape(levels)[0] - 1) / 2)
    x = 0.5
    cmap = plt.cm.RdBu_r
    lower = cmap(np.linspace(0, x, n))
    white = np.ones(((np.shape(levels)[0] - 1)-2*n,4))
    upper = cmap(np.linspace(1-x, 1, n))
    colors = np.vstack((lower, white, upper))
    cmap_f = matplotlib.colors.LinearSegmentedColormap.from_list('map_white', colors)
    return cmap_f

fig, axes = plt.subplots(nrows=4, ncols=2, figsize = (11,9.5))
#fig = plt.subplots(figsize = (12,9.5))
## ax 0, 0
lats_2d, levs_2d = np.meshgrid(lats, levs)
#ax = plt.subplot2grid((3, 2), (1, 0), colspan=1, rowspan = 1)
ax = axes[0, 0]
levels = np.arange(-4,4.5,0.5)
contours = ax.contourf(lats_2d, levs_2d, chg_salt_cs_atl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
#ax.plot(lats_2d[0,:], mld_atl_zonal, '-k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('S (psu; $\Delta$(standard) - $\Delta$(nudging))', fontsize=12, loc = 'left')
ax.set_title('S (10$^\mathrm{6}$ psu·m; Atlantic; FLOR)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'a', transform=ax.transAxes, fontsize=12, weight='bold')
ax.text(0.10, 1.25, 'Model response to CO2 doubling', transform=ax.transAxes, fontsize=14)
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-4, 6., 2.))
cbar.ax.tick_params()


## ax 1, 0
#ax = plt.subplot2grid((3, 2), (1, 1), colspan=1, rowspan = 1)
ax = axes[1, 0]
levels = np.arange(-20, 21, 1) 
cmap_f = new_cmap(levels)
contours = ax.contourf(lats_2d, levs_2d, chg_temp_cs_atl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
lats_2d_cont = lats_2d[:,lats_2d[0,:]>40]
levs_2d_cont = levs_2d[:,lats_2d[0,:]>40]
chg_temp_cs_atl_2co2std_cont = chg_temp_cs_atl_2co2std[:,lats_2d[0,:]>40]
contourx = ax.contour(lats_2d_cont, levs_2d_cont, chg_temp_cs_atl_2co2std_cont, levels = [0], colors='k',linestyles='-',linewidths=1)
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard) - $\Delta$(nudging))', fontsize=12, loc = 'left')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m; Atlantic; FLOR)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'b', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-20, 30, 10))
cbar.ax.tick_params()
# ax.grid(True)
ax.clabel(contourx, contourx.levels, inline=True, fmt='0', fontsize=10, manual = [(55, 1000)])


## ax 2, 0
lats_2d, levs_2d = np.meshgrid(lats, levs)
#ax = plt.subplot2grid((3, 2), (2, 0), colspan=1, rowspan = 1)
ax = axes[2, 0]
levels = np.arange(-4,4.5,0.5)
contours = ax.contourf(lats_2d, levs_2d, chg_salt_cs_nonatl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
#ax.plot(lats_2d[0,:], mld_nonatl_zonal, '-k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('S (psu;  $\Delta$(standard))', fontsize=12, loc = 'left')
ax.set_title('S (10$^\mathrm{6}$ psu·m; Non-Atlantic; FLOR)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'c', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-4,6.,2.))
cbar.ax.tick_params()


## ax 3, 0
#ax = plt.subplot2grid((3, 2), (2, 1), colspan=1, rowspan = 1)
ax = axes[3, 0]
levels = np.arange(-30, 31, 1)
contours = ax.contourf(lats_2d, levs_2d, chg_temp_cs_nonatl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard))', fontsize=12, loc = 'left')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m; Non-Atlantic; FLOR)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'd', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-30, 40, 10))
cbar.ax.tick_params()


# ax 0, 1
levels = np.arange(-4,4.5,0.5)
#ax = plt.subplot2grid((4, 2), (3, 0), colspan=1, rowspan = 1)
ax = axes[0, 1]
contours = ax.contourf(lats_2d_NCEI, depth_2d_NCEI, slope_salt_atl_NCEI*50*1.74, levels = levels, cmap='RdBu_r',extend='both')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_title('S (10$^\mathrm{6}$ psu·m/50yr; Atlantic; NCEI)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'e', transform=ax.transAxes, fontsize=12, weight='bold')
ax.text(0.35, 1.25, 'Observation', transform=ax.transAxes, fontsize=14)
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-4,6,2))
cbar.ax.tick_params()


# ax 1, 1
levels = np.arange(-20, 21, 1) 
ax = axes[1, 1]
#ax = plt.subplot2grid((4, 2), (3, 1), colspan=1, rowspan = 1)
contours = ax.contourf(lats_2d_NCEI, depth_2d_NCEI, slope_temp_atl_NCEI*50*1.74, levels = levels, cmap='RdBu_r',extend='both')
lats_2d_NCEI_cont = lats_2d_NCEI[:,lats_2d_NCEI[0,:]>30]
depth_2d_NCEI_cont = depth_2d_NCEI[:,lats_2d_NCEI[0,:]>30]
slope = slope_temp_atl_NCEI*50*1.74
slope_cont = slope[:,lats_2d_NCEI[0,:]>30]
contourx = ax.contour(lats_2d_NCEI_cont, depth_2d_NCEI_cont, slope_cont, levels = [0], colors='k',linestyles='-',linewidths=1)
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
#ax.set_ylabel('Depth (m)')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m/50yr; Atlantic; NCEI)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'f', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-20, 30, 10))
cbar.ax.tick_params()
ax.clabel(contourx, contourx.levels, inline=True, fmt='0', fontsize=10, manual = [(50, 1500)])


# ax 2, 1
levels = np.arange(-4,4.5,0.5)
#ax = plt.subplot2grid((4, 2), (3, 0), colspan=1, rowspan = 1)
ax = axes[2, 1]
contours = ax.contourf(lats_2d_NCEI, depth_2d_NCEI, slope_salt_nonatl_NCEI*50*1.74, levels = levels, cmap='RdBu_r',extend='both')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_title('S (10$^\mathrm{6}$ psu·m/50yr; Non-Atlantic; NCEI)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'g', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-4,6.,2.))
cbar.ax.tick_params()


# ax 3, 1
levels = np.arange(-30, 31, 1)
ax = axes[3, 1]
#ax = plt.subplot2grid((4, 2), (3, 1), colspan=1, rowspan = 1)
contours = ax.contourf(lats_2d_NCEI, depth_2d_NCEI, slope_temp_nonatl_NCEI*50*1.74, levels = levels, cmap='RdBu_r',extend='both')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
#ax.set_ylabel('Depth (m)')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m/50yr; Non-Atlantic; NCEI)', fontsize=12, loc = 'left')
ax.text(-0.12, 1.05, 'h', transform=ax.transAxes, fontsize=12, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-30, 40, 10))
cbar.ax.tick_params()

plt.subplots_adjust(wspace = 0.39)
plt.subplots_adjust(hspace = 0.42)
plt.savefig('FIG5_1st_review.pdf', bbox_inches='tight')


 
            
            
            