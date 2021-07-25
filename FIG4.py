## load libraries
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import gsw as gsw

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
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]
fid.close()

mask_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lats)[0], np.shape(lons)[0]))
for i in np.arange(0,np.shape(mask_3d)[0]):
    mask_3d[i,:,:] = mask_i
del mask_i

lon_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lats)[0], np.shape(lons)[0]))
lat_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lats)[0], np.shape(lons)[0]))
depth_3d = np.nan * np.zeros((np.shape(levs)[0], np.shape(lats)[0], np.shape(lons)[0]))
lon_2d_i, lat_2d_i = np.meshgrid(lons, lats)
for i in np.arange(0, np.shape(lon_3d)[0], 1):
    lon_3d[i,:,:] = lon_2d_i.copy()
    lat_3d[i,:,:] = lat_2d_i.copy()
    depth_3d[i,:,:] = -levs[i].copy()
p = gsw.p_from_z(depth_3d,lat_3d)

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
rho_ctlstd = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho_ctlstd[rho_ctlstd.mask==True] = np.nan
temp1 = fid.variables['temp'][:] #  oC
temp1[temp1.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()

salt_abs_ctlstd = gsw.SA_from_SP(salt,p,lon_3d,lat_3d)
CT = gsw.CT_from_t(salt_abs_ctlstd, temp1, p)   
[rho_i, alpha_ctlstd, beta_ctlstd] = gsw.rho_alpha_beta(salt_abs_ctlstd, CT, p)


yr_start = 161
yr_end = 180
sim_type = '2CO2std'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
rho_co2std = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho_co2std[rho_co2std.mask==True] = np.nan
temp2 = fid.variables['temp'][:] #  oC
temp2[temp2.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()
salt_abs_co2std = gsw.SA_from_SP(salt,p,lon_3d,lat_3d)
CT = gsw.CT_from_t(salt_abs_co2std, temp2, p)   
[rho_i, alpha_co2std, beta_co2std] = gsw.rho_alpha_beta(salt_abs_co2std, CT, p)

rho_beta_chg_std = (beta_ctlstd + beta_co2std) / 2. * (salt_abs_co2std - salt_abs_ctlstd) * rho_ctlstd
rho_alpha_chg_std = -(alpha_ctlstd + alpha_co2std) / 2. * (temp2 - temp1) * rho_ctlstd

rho_chg_std = rho_co2std - rho_ctlstd

integral_params = np.zeros((np.shape(temp1)[0], np.shape(temp1)[1], np.shape(temp1)[2]))
lons_2d_intergral, lats_2d_intergral = np.meshgrid(lons, lats)
for i in np.arange(0, np.shape(integral_params)[0], 1):
    integral_params[i,:,:] = 110 * 1000 * np.cos(lats_2d_intergral * np.pi / 180.)
    
temp = rho_beta_chg_std.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_beta_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_beta_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_alpha_chg_std.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_alpha_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_alpha_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_chg_std.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_ctlstd.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_ctlstd = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_ctlstd = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_co2std.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_co2std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_co2std = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

chg_temp_cs_atl_std = temp_cs_atl_co2std - temp_cs_atl_ctlstd
chg_temp_cs_nonatl_std = temp_cs_nonatl_co2std - temp_cs_nonatl_ctlstd

sim_type = 'ctlgl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
rho_ctlgl = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho_ctlgl[rho_ctlgl.mask==True] = np.nan
temp1 = fid.variables['temp'][:] #  oC
temp1[temp1.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()
salt_abs_ctlgl = gsw.SA_from_SP(salt,p,lon_3d,lat_3d)
CT = gsw.CT_from_t(salt_abs_ctlgl, temp1, p)   
[rho_i, alpha_ctlgl, beta_ctlgl] = gsw.rho_alpha_beta(salt_abs_ctlgl, CT, p)


yr_start = 161
yr_end = 180
sim_type = '2CO2gl'
fid = Dataset('3d_ocean_data_'+sim_type+'_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
levs = fid.variables['lev'][:]
lats = fid.variables['lat'][:]
lons = fid.variables['lon'][:]
rho_co2gl = fid.variables['rho'][:]  # kg/m3  (lev lat lon)
rho_co2gl[rho_co2gl.mask==True] = np.nan
temp2 = fid.variables['temp'][:] #  oC
temp2[temp2.mask==True] = np.nan
salt = fid.variables['salt'][:] # psu
salt[salt.mask==True] = np.nan
v = fid.variables['v'][:] # m/s
v[v.mask==True] = np.nan
fid.close()
salt_abs_co2gl = gsw.SA_from_SP(salt,p,lon_3d,lat_3d)
CT = gsw.CT_from_t(salt_abs_co2gl, temp2, p)   
[rho_i, alpha_co2gl, beta_co2gl] = gsw.rho_alpha_beta(salt_abs_co2gl, CT, p)


rho_beta_chg_gl = (beta_ctlgl + beta_co2gl) / 2. * (salt_abs_co2gl - salt_abs_ctlgl) * rho_ctlgl
rho_alpha_chg_gl = -(alpha_ctlgl + alpha_co2gl) / 2. * (temp2 - temp1) * rho_ctlgl

rho_chg_gl = rho_co2gl - rho_ctlgl

temp = rho_beta_chg_gl.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_beta_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_beta_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_alpha_chg_gl.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_alpha_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_alpha_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_chg_gl.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_ctlgl.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_ctlgl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_ctlgl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

temp = rho_co2gl.copy()
temp_tmp = temp.copy()
temp_tmp[mask_3d == 1.0] = np.nan
temp_cs_nonatl_co2gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp
temp_tmp = temp.copy()
temp_tmp[mask_3d != 1.0] = np.nan
temp_cs_atl_co2gl = np.nansum(temp_tmp * integral_params, axis = 2) / 10.**6
del temp_tmp

chg_temp_cs_atl_gl = temp_cs_atl_co2gl - temp_cs_atl_ctlgl
chg_temp_cs_nonatl_gl = temp_cs_nonatl_co2gl - temp_cs_nonatl_ctlgl


#######################################################################################            
################### plot cross sections of temp and salt GL - STD #####################
#######################################################################################
fig, axes = plt.subplots(nrows=3, ncols=2, figsize = (12,7.5))
#fig = plt.subplots(figsize = (12,4.5))
plt.rcParams.update({'font.family':'Arial'})
lats_2d, levs_2d = np.meshgrid(lats, levs)

## ax 0, 0
ax = axes[0, 0]
levels = np.arange(-2, 2.1, 0.1) 
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_atl_gl + temp_cs_atl_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_atl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_atl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.set_ylabel('Depth (m)')
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'a', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard) - $\Delta$(nudging))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'c', transform=ax.transAxes, fontsize=10, weight='bold')

ax = axes[0, 1]
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_nonatl_gl + temp_cs_nonatl_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_nonatl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_nonatl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'b', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Non-Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'e', transform=ax.transAxes, fontsize=10, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-8, 10, 1))
cbar.ax.tick_params()


## ax 1, 0
ax = axes[1, 0]
levels = np.arange(-2, 2.1, 0.1) 
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_atl_beta_gl + temp_cs_atl_beta_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_atl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_atl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.set_ylabel('Depth (m)')
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'c', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard) - $\Delta$(nudging))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1_{S}$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'c', transform=ax.transAxes, fontsize=10, weight='bold')

ax = axes[1, 1]
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_nonatl_beta_gl + temp_cs_nonatl_beta_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_nonatl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_nonatl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'd', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1_{S}$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Non-Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'e', transform=ax.transAxes, fontsize=10, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-8, 10, 1))
cbar.ax.tick_params()



## ax 2, 0
ax = axes[2, 0]
levels = np.arange(-2, 2.1, 0.1) 
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_atl_alpha_gl + temp_cs_atl_alpha_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_atl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_atl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.set_ylabel('Depth (m)')
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'e', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard) - $\Delta$(nudging))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1_{T}$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'c', transform=ax.transAxes, fontsize=10, weight='bold')

ax = axes[2, 1]
contours = ax.contourf(lats_2d, levs_2d, -temp_cs_nonatl_alpha_gl + temp_cs_nonatl_alpha_std, levels = levels, cmap='RdBu_r',extend='both', alpha = 1)
# ax.plot(lats_2d[0,:], mld_nonatl_ctlgl_zonal, '-k')
# ax.plot(lats_2d[0,:], mld_nonatl_co2gl_zonal, '--k')
ax.set_xlim((-60,60))
ax.set_xticks(np.arange(-60,80,20))
ax.set_xticklabels(['60$^\mathrm{o}$S','40$^\mathrm{o}$S','20$^\mathrm{o}$S','0$^\mathrm{o}$',\
                   '20$^\mathrm{o}$N','40$^\mathrm{o}$N','60$^\mathrm{o}$N'])
ax.set_ylim((0, 2000))
ax.invert_yaxis()
ax.text(-0.1, 1.1, 'f', transform=ax.transAxes, fontsize=11, weight='bold')
#ax.set_title('T ($^\mathrm{o}$C; $\Delta$(standard))', fontsize=10, loc = 'left')
ax.set_title('$\u03C1_{T}$ (10$^\mathrm{6}$ kg m$^\mathrm{-2}$; Non-Atlantic)', fontsize=10, loc = 'left')
#ax.text(-0.12, 1.05, 'e', transform=ax.transAxes, fontsize=10, weight='bold')
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.03, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-8, 10, 1))
cbar.ax.tick_params()

plt.subplots_adjust(wspace = 0.18)
plt.subplots_adjust(hspace = 0.37)
plt.savefig('FIG4.pdf', bbox_inches='tight')



 
            
            
            
