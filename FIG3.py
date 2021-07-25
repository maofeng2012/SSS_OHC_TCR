## load libraries
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

################################# read ocean mask #############################
yr_start = 101
yr_end = 200
sim_type = 'ctl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
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

sim_type = '2CO2gl'
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

############################################################################ 
#################################### mld ###################################
############################################################################ 
fid = Dataset('OHC_density_ctl_101_200.nc')
ohc_ctl_all = fid.variables['ohc_all'][:]
lat_mld = fid.variables['lat'][:]
lon_mld = fid.variables['lon'][:]
#lon_2d, lat_2d = np.meshgrid(lon, lat)
fid.close()

mld_ctl = np.loadtxt('mld_winter_ctl_101_200_pattern.txt', delimiter=' ') # lat, lon
mld_ctlgl = np.loadtxt('mld_winter_ctlgl_101_200_pattern.txt', delimiter=' ')
mld_co2std = np.loadtxt('mld_winter_2co2std_161_180_pattern.txt', delimiter=' ')
mld_co2gl = np.loadtxt('mld_winter_2co2gl_161_180_pattern.txt', delimiter=' ')

temp = mld_ctl.copy()
temp[mask_2d != 1] = np.nan
mld_atl_ctl_zonal = np.nanmean(temp, axis = 1)
temp = mld_ctlgl.copy()
temp[mask_2d != 1] = np.nan
mld_atl_ctlgl_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2std.copy()
temp[mask_2d != 1] = np.nan
mld_atl_co2std_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2gl.copy()
temp[mask_2d != 1] = np.nan
mld_atl_co2gl_zonal = np.nanmean(temp, axis = 1)

temp = mld_ctl.copy()
temp[mask_2d == 1] = np.nan
mld_nonatl_ctl_zonal = np.nanmean(temp, axis = 1)
temp = mld_ctlgl.copy()
temp[mask_2d == 1] = np.nan
mld_nonatl_ctlgl_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2std.copy()
temp[mask_2d == 1] = np.nan
mld_nonatl_co2std_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2gl.copy()
temp[mask_2d == 1] = np.nan
mld_nonatl_co2gl_zonal = np.nanmean(temp, axis = 1)

temp = mld_ctl.copy()
mld_all_ctl_zonal = np.nanmean(temp, axis = 1)
temp = mld_ctlgl.copy()
mld_all_ctlgl_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2std.copy()
mld_all_co2std_zonal = np.nanmean(temp, axis = 1)
temp = mld_co2gl.copy()
mld_all_co2gl_zonal = np.nanmean(temp, axis = 1)

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
fid = Dataset('ocean_index_with_ATL_PAC_IND_mod.nc')
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
#fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (11,5))
fig = plt.subplots(figsize = (12,7.5))
plt.rcParams.update({'font.family':'Arial'})
## ax 0, 0-1
ax = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan = 1)
x = np.arange(0, 5, 1)
bar_plot = ax.bar(x-0.15, OHC_all, 0.15, label='Globe', color = '#00a088')
bar_plot = ax.bar(x, OHC_atl, 0.15, label='Atlantic', color = '#4dbbd6')
bar_plot = ax.bar(x+0.15, OHC_nonatl, 0.15, label='Non-Atlantic', color = '#e64b34')
ax.legend(fontsize = 10,frameon=False)
ax.set_ylabel('10$^\mathrm{24}$ J', fontsize=10)
ax.set_title('OHC', loc = 'left', fontsize=10)
ax.text(0.9, 0.15, 'The impact of fixed SSS on model response to CO$_2$ doubling', fontsize=10, weight='bold')
ax.set_xticks(x)
#ax.set_yticks(np.arange(-1, 6, 1))
ax.set_ylim(-0.06,0.11)
ax.set_xticklabels(np.array(['All', '<300 m', '300-700 m', '700-2000 m', '>2000 m']))
x1 = np.arange(-0.4, 4.7, 0.1)
y1 = 0 * np.arange(-0.4, 4.7, 0.1)
ax.plot(x1, y1, '-', color = 'gray', linewidth = 0.5)
ax.set_xlim((-0.4, 4.4))
ax.tick_params(labelsize = 10)
ax.text(-0.05, 1.05, 'a', transform=ax.transAxes, fontsize=10, weight='bold')

## plot Atlantic mask
axin = inset_axes(ax,
                   width="30%",  # width = 10% of parent_bbox width
                   height="40%",  # height : 50%
                   loc='center',
                   bbox_to_anchor=(-0.24, 0.26, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
#                   loc='lower left',
m = Basemap(projection='cea', llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-329.5,urcrnrlon=29.5,resolution='c')
lons_2d, lats_2d = np.meshgrid(lon_mask, lat_mask)
lons_2d_plot = np.concatenate((lons_2d[:,lons_2d[0,:]>30]-360., lons_2d[:,lons_2d[0,:]<=30]),\
                          axis = 1).copy()
lats_2d_plot = np.concatenate((lats_2d[:,lons_2d[0,:]>30], lats_2d[:,lons_2d[0,:]<=30]),\
                          axis = 1).copy()
ocean_mask_plot = np.concatenate((ocean_mask[:,lons_2d[0,:]>30], ocean_mask[:,lons_2d[0,:]<=30]),\
                          axis = 1).copy()
ocean_mask_plot[ocean_mask_plot != 1] = 2
ocean_mask_plot[ocean_mask_plot == 2] = 3
ocean_mask_plot[ocean_mask_plot == 1] = 2 # Atlantic
ocean_mask_plot[ocean_mask_plot == 3] = 1 # Non-ATtlantic

#NonATL_mask_plot = ocean_mask_plot.copy()
#NonATL_mask_plot[NonATL_mask_plot == 1] = np.nan
#
#NonATL_mask_plot[NonATL_mask_plot != 1] = 2

x, y = m(lons_2d_plot, lats_2d_plot)
levels = np.arange(0,3,1) 
cmap = ListedColormap(['#4dbbd6', '#e64b34'])

contours = m.contourf(x, y,ocean_mask_plot,cmap=cmap, levels = levels, ax = axin)
m.drawcoastlines(ax = axin, linewidth = 0.2, color = 'black')
#m.drawparallels(np.arange(int(-60),int(90),30),labels=[1,0,0,0], linewidth = 0, ax = axin)
#m.drawparallels(np.array([0]),labels=[1,0,0,0], linewidth = 1, ax = axin)
#m.drawmeridians(np.arange(int(0),int(360),60),labels=[0,0,0,1], linewidth = 0, ax = axin)
m.fillcontinents(color='grey', ax = axin)
axins = inset_axes(axin,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
                   bbox_transform=axin.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=axin, cax = axins)
cbar.set_ticks([])
cbar.ax.text(1.2, 0.15, 'Non-ATL')
cbar.ax.text(1.2, 0.65, 'ATL')


## ax 1, 0
lats_2d, levs_2d = np.meshgrid(lats, levs)
ax = plt.subplot2grid((3, 2), (1, 0), colspan=1, rowspan = 1)
levels = np.arange(-4,4.5,0.5) 
#contours_1 = ax.contourf(lats_2d[:,lats_2d[0,:]<-33], levs_2d[:,lats_2d[0,:]<-33], -chg_salt_cs_all_2co2gl[:,lats_2d[0,:]<-33] + chg_salt_cs_all_2co2std[:,lats_2d[0,:]<-33], levels = levels, cmap='RdBu_r',extend='both', alpha = 0.6)
#contours = ax.contourf(lats_2d[:,lats_2d[0,:]>=-34], levs_2d[:,lats_2d[0,:]>=-34], -chg_salt_cs_all_2co2gl[:,lats_2d[0,:]>=-34] + chg_salt_cs_all_2co2std[:,lats_2d[0,:]>=-34], levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
contours = ax.contourf(lats_2d, levs_2d, -chg_salt_cs_atl_2co2gl + chg_salt_cs_atl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
#ax.plot(lats_2d[0,:], mld_atl_zonal_ctl, '-k')
#ax.plot(lats_2d[0,:], mld_atl_zonal_co2, '--k')
ax.plot(lats_2d[0,:], mld_atl_ctl_zonal, '-k')
ax.plot(lats_2d[0,:], mld_atl_co2std_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('S (psu; $\Delta$(standard) - $\Delta$(nudging))', fontsize=10, loc = 'left')
ax.set_title('S (10$^\mathrm{6}$ psu·m; Atlantic)', fontsize=10, loc = 'left')
ax.text(-0.12, 1.05, 'b', transform=ax.transAxes, fontsize=10, weight='bold')
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


## ax 1, 1
ax = plt.subplot2grid((3, 2), (1, 1), colspan=1, rowspan = 1)
levels = np.arange(-8, 8.5, 0.5) 
contours = ax.contourf(lats_2d, levs_2d, -chg_temp_cs_atl_2co2gl + chg_temp_cs_atl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
ax.plot(lats_2d[0,:], mld_atl_ctlgl_zonal, '-k')
ax.plot(lats_2d[0,:], mld_atl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard) - $\Delta$(nudging))', fontsize=10, loc = 'left')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m; Atlantic)', fontsize=10, loc = 'left')
ax.text(-0.12, 1.05, 'c', transform=ax.transAxes, fontsize=10, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-8, 12, 4))
cbar.ax.tick_params()


## ax 2, 0
lats_2d, levs_2d = np.meshgrid(lats, levs)
ax = plt.subplot2grid((3, 2), (2, 0), colspan=1, rowspan = 1)
levels = np.arange(-4,4.5,0.5) 
#contours_1 = ax.contourf(lats_2d[:,lats_2d[0,:]<-33], levs_2d[:,lats_2d[0,:]<-33], -chg_salt_cs_all_2co2gl[:,lats_2d[0,:]<-33] + chg_salt_cs_all_2co2std[:,lats_2d[0,:]<-33], levels = levels, cmap='RdBu_r',extend='both', alpha = 0.6)
#contours = ax.contourf(lats_2d[:,lats_2d[0,:]>=-34], levs_2d[:,lats_2d[0,:]>=-34], -chg_salt_cs_all_2co2gl[:,lats_2d[0,:]>=-34] + chg_salt_cs_all_2co2std[:,lats_2d[0,:]>=-34], levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
contours = ax.contourf(lats_2d, levs_2d, -chg_salt_cs_nonatl_2co2gl + chg_salt_cs_nonatl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
mld_atl_zonal_ctl = (mld_atl_ctl_zonal + mld_atl_ctlgl_zonal) / 2.
mld_atl_zonal_co2 = (mld_atl_co2std_zonal + mld_atl_co2gl_zonal) / 2.
mld_nonatl_zonal_ctl = (mld_nonatl_ctl_zonal + mld_nonatl_ctlgl_zonal) / 2.
mld_nonatl_zonal_co2 = (mld_nonatl_co2std_zonal + mld_nonatl_co2gl_zonal) / 2.
#ax.plot(lats_2d[0,:], mld_nonatl_zonal_ctl, '-k')
#ax.plot(lats_2d[0,:], mld_nonatl_zonal_co2, '--k')
ax.plot(lats_2d[0,:], mld_nonatl_ctl_zonal, '-k')
ax.plot(lats_2d[0,:], mld_nonatl_co2std_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0,2000))
ax.invert_yaxis()
ax.set_ylabel('Depth (m)')
#ax.set_title('S (psu;  $\Delta$(standard))', fontsize=10, loc = 'left')
ax.set_title('S (10$^\mathrm{6}$ psu·m;  Non-Atlantic)', fontsize=10, loc = 'left')
ax.text(-0.12, 1.05, 'd', transform=ax.transAxes, fontsize=10, weight='bold')
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


## ax 2, 1
ax = plt.subplot2grid((3, 2), (2, 1), colspan=1, rowspan = 1)
levels = np.arange(-8, 8.5, 0.5)  
contours = ax.contourf(lats_2d, levs_2d, -chg_temp_cs_nonatl_2co2gl + chg_temp_cs_nonatl_2co2std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
ax.plot(lats_2d[0,:], mld_nonatl_ctlgl_zonal, '-k')
ax.plot(lats_2d[0,:], mld_nonatl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard))', fontsize=10, loc = 'left')
ax.set_title('T (10$^\mathrm{6}$ $^\mathrm{o}$C·m; Non-Atlantic)', fontsize=10, loc = 'left')
ax.text(-0.12, 1.05, 'e', transform=ax.transAxes, fontsize=10, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-8, 12, 4))
cbar.ax.tick_params()

plt.subplots_adjust(wspace = 0.32)
plt.subplots_adjust(hspace = 0.42)
plt.savefig('FIG3.pdf', bbox_inches='tight')



 
            
            
            
