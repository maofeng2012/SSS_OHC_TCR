## load libraries
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def moving_average(a, n=10) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


################################################################################
################################# Read P-E data ################################
################################################################################
lon_2d_i = np.loadtxt('lon_atmos.txt', delimiter = ' ')
lat_2d = np.loadtxt('lat_atmos.txt', delimiter = ' ')

lon_2d_p_e = np.concatenate((lon_2d_i[:,lon_2d_i[0,:]>30]-360., lon_2d_i[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()
lat_2d_p_e = np.concatenate((lat_2d[:,lon_2d_i[0,:]>30], lat_2d[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy() 


# read data
p_e_ctl = np.loadtxt('p_e_ctl_101_200_pattern.txt', delimiter = ' ')
p_e_ctlgl = np.loadtxt('p_e_ctlgl_101_200_pattern.txt', delimiter = ' ')
p_e_2co2std = np.loadtxt('p_e_2co2std_161_180_pattern.txt', delimiter = ' ')
p_e_2co2gl = np.loadtxt('p_e_2co2gl_161_180_pattern.txt', delimiter = ' ')

chg_p_e_2co2std = p_e_2co2std - p_e_ctl

p_e_ctl = np.concatenate((p_e_ctl[:,lon_2d_i[0,:]>30], \
                              p_e_ctl[:,lon_2d_i[0,:]<=30]), axis = 1).copy()
chg_p_e_2co2std = np.concatenate((chg_p_e_2co2std[:,lon_2d_i[0,:]>30], \
                              chg_p_e_2co2std[:,lon_2d_i[0,:]<=30]), axis = 1).copy()


################################################################################
################################# Read data ####################################
################################################################################
yr_start = 101
yr_end = 219

lon = np.loadtxt('lon_ocean.txt', delimiter = ' ')
lat = np.loadtxt('lat_ocean.txt', delimiter = ' ')
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


######################## read TOA NETRAD data ######################## 
netrad_toa_ctl = np.loadtxt('netrad_toa_ctl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
netrad_toa_ctlgl = np.loadtxt('netrad_toa_ctlgl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
netrad_toa_2co2std = np.loadtxt('netrad_toa_2CO2std_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
netrad_toa_2co2gl = np.loadtxt('netrad_toa_2CO2gl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')

netrad_toa_ctl = moving_average(netrad_toa_ctl, n=10)
netrad_toa_ctlgl = moving_average(netrad_toa_ctlgl, n=10)

netrad_toa_2co2std = moving_average(netrad_toa_2co2std, n=10)
netrad_toa_2co2gl = moving_average(netrad_toa_2co2gl, n=10)

######################## read all-depth OHC data ######################## 
fid = Dataset('Annual_OHC_ctl_lat_lon_year_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
ohc_ctl = fid.variables['ohc'][:,:,:]
fid.close()
fid = Dataset('Annual_OHC_ctlgl_lat_lon_year_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
ohc_ctlgl = fid.variables['ohc'][:,:,:]
fid.close()
fid = Dataset('Annual_OHC_2CO2std_lat_lon_year_'+np.str(yr_start)+'_'+np.str(yr_end)+'.nc')
ohc_2co2std = fid.variables['ohc'][:,:,:]
fid.close()


fid = Dataset('Annual_OHC_2CO2gl_lat_lon_year_101_160.nc')
ohc_2co2gl_1 = fid.variables['ohc'][:,:,:]
fid.close()
fid = Dataset('Annual_OHC_2CO2gl_lat_lon_year_161_219.nc')
ohc_2co2gl_2 = fid.variables['ohc'][:,:,:]
fid.close()
ohc_2co2gl = np.concatenate((ohc_2co2gl_1, ohc_2co2gl_2), axis = 0)


chg_ohc_std = np.zeros((np.shape(ohc_ctl)[0],))
chg_ohc_gl = np.zeros((np.shape(ohc_ctl)[0],))

for yr in np.arange(0, np.shape(ohc_ctl)[0], 1):
    chg_ohc_std[yr] = np.nansum((ohc_2co2std[yr,:,:] - ohc_ctl[yr,:,:]) * area) / 10.**24 * 10.**9
    chg_ohc_gl[yr] = np.nansum((ohc_2co2gl[yr,:,:] - ohc_ctlgl[yr,:,:]) * area) / 10.**24 * 10.**9


############################################################################################            
############################### read OHC data for all runs #################################
############################################################################################ 
lon_2d_i = np.loadtxt('lon_ocean.txt', delimiter = ' ')
lat_2d = np.loadtxt('lat_ocean.txt', delimiter = ' ')
lon_2d = np.concatenate((lon_2d_i[:,lon_2d_i[0,:]>30]-360., lon_2d_i[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()
lat_2d = np.concatenate((lat_2d[:,lon_2d_i[0,:]>30], lat_2d[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()   

fid = Dataset('OHC_density_ctl_101_200.nc')
ohc_ctl_all = fid.variables['ohc_all'][:]
fid.close()
fid = Dataset('OHC_density_ctlgl_101_200.nc')
ohc_ctlgl_all = fid.variables['ohc_all'][:]
fid.close()
fid = Dataset('OHC_density_2CO2std_161_180.nc')
ohc_co2std_all = fid.variables['ohc_all'][:]
fid.close()
fid = Dataset('OHC_density_2CO2gl_161_180.nc')
ohc_co2gl_all = fid.variables['ohc_all'][:]
fid.close()


chg_ohc_co2std_all = ohc_co2std_all - ohc_ctl_all
chg_ohc_co2gl_all = ohc_co2gl_all - ohc_ctlgl_all

ohc_ctl_all = np.concatenate((ohc_ctl_all[:,lon_2d_i[0,:]>30], \
                              ohc_ctl_all[:,lon_2d_i[0,:]<=30]), axis = 1).copy()
ohc_co2std_all = np.concatenate((ohc_co2std_all[:,lon_2d_i[0,:]>30], \
                              ohc_co2std_all[:,lon_2d_i[0,:]<=30]), axis = 1).copy()
ohc_co2gl_all = np.concatenate((ohc_co2gl_all[:,lon_2d_i[0,:]>30], \
                              ohc_co2gl_all[:,lon_2d_i[0,:]<=30]), axis = 1).copy()

chg_ohc_co2std_all = np.concatenate((chg_ohc_co2std_all[:,lon_2d_i[0,:]>30], \
                              chg_ohc_co2std_all[:,lon_2d_i[0,:]<=30]), axis = 1).copy()
chg_ohc_co2gl_all = np.concatenate((chg_ohc_co2gl_all[:,lon_2d_i[0,:]>30], \
                              chg_ohc_co2gl_all[:,lon_2d_i[0,:]<=30]), axis = 1).copy()


############################################################################################            
############################### read SSS data for all runs #################################
############################################################################################ 
input_dir = '/Users/jiexue/Documents/Projects/Salinity/Work/Basics/Data_new/'
yr_start = '161'
yr_end = '180'

# read relative sst change
sss_ctl = np.loadtxt('5m_salinity_ctl_101_200_pattern.txt', delimiter = ' ')
sss_ctlgl = np.loadtxt('5m_salinity_ctlgl_101_200_pattern.txt', delimiter = ' ')

sss_co2std = np.loadtxt('5m_salinity_2CO2std_' + yr_start+'_'+yr_end+'_pattern.txt', delimiter = ' ')
sss_co2gl = np.loadtxt('5m_salinity_2CO2gl_' + yr_start+'_'+yr_end+'_pattern.txt', delimiter = ' ')

chg_sss_co2std = sss_co2std - sss_ctl
chg_sss_co2gl = sss_co2gl - sss_ctlgl


chg_sss_co2std = np.concatenate((chg_sss_co2std[:,lon_2d_i[0,:]>30], chg_sss_co2std[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()  
chg_sss_co2gl = np.concatenate((chg_sss_co2gl[:,lon_2d_i[0,:]>30], chg_sss_co2gl[:,lon_2d_i[0,:]<=30]),\
                          axis = 1).copy()  


#######################################################################################            
############################## plot SSS change CO2 - CTL ##############################
#######################################################################################
m = Basemap(projection='cea', llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-329.5,urcrnrlon=29.5,resolution='c')
plt.rcParams.update({'font.family':'Arial'})
fig = plt.figure(figsize = (12,4.8))
x, y = m(lon_2d, lat_2d)
levels_sss = np.arange(-0.8,0.9,0.1) 
levels_ohc_0 = np.arange(-12,13,1) 
levels_ohc_123 = np.arange(-4,4.5,0.5) 

# ax 0, 1 
ax = plt.subplot2grid((3, 2), (0, 1),colspan=1, rowspan = 1)
contours = m.contourf(x, y, chg_ohc_co2std_all - chg_ohc_co2gl_all,levels=levels_ohc_123,cmap='RdBu_r',extend='both',ax = ax)
m.drawcoastlines(ax = ax, linewidth = 0.2, color = 'black')
m.drawparallels(np.arange(int(-60),int(90),30),labels=[1,0,0,0], linewidth = 0, ax = ax)
m.drawparallels(np.array([0]),labels=[0,0,0,0], linewidth = 1, ax = ax)
m.drawmeridians(np.arange(int(0),int(360),60),labels=[0,0,0,0], linewidth = 0, ax = ax)
m.fillcontinents(color='grey', ax = ax)
ax.text(-0.1, 1.1, 'b', transform=ax.transAxes, fontsize=10, weight='bold')
ax.set_title('OHC (10$^{9}$ J m$^{-2}$)', fontsize = 10, loc = 'left')  
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax = ax, cax = axins, ticks = np.arange(-4.,6.,2.))

# ax 1, 1 
ax = plt.subplot2grid((3, 2), (1, 1),colspan=1, rowspan = 1)
contours = m.contourf(x, y, chg_sss_co2std - chg_sss_co2gl,levels=levels_sss,cmap='RdBu_r',extend='both',ax = ax)
m.drawcoastlines(ax = ax, linewidth = 0.2, color = 'black')
m.drawparallels(np.arange(int(-60),int(90),30),labels=[1,0,0,0], linewidth = 0, ax = ax)
m.drawparallels(np.array([0]),labels=[1,0,0,0], linewidth = 1, ax = ax)
m.drawmeridians(np.arange(int(0),int(360),60),labels=[0,0,0,0], linewidth = 0, ax = ax)
m.fillcontinents(color='grey', ax = ax)
ax.text(-0.1, 1.1, 'c', transform=ax.transAxes, fontsize=10, weight='bold')
ax.set_title('SSS (psu)', fontsize = 10, loc = 'left') 
# colorbar
axins = inset_axes(ax,
                  width="3.5%",  # width = 10% of parent_bbox width
                  height="100%",  # height : 50%
                  loc='lower left',
                  bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
                  bbox_transform=ax.transAxes,
                  borderpad=0,
                  )
cbar = plt.colorbar(contours, ax = ax, cax = axins, ticks = np.arange(-0.8,1.2,0.4))
cbar.ax.tick_params() 

# ax 2, 1

x, y = m(lon_2d_p_e, lat_2d_p_e)
ax = plt.subplot2grid((3, 2), (2, 1),colspan=1, rowspan = 1)
levels_p_e = np.arange(-1.5,1.75,0.25) 
contours = m.contourf(x, y, chg_p_e_2co2std, levels = levels_p_e, cmap='BrBG',extend='both', ax = ax)
m.drawcoastlines(ax = ax, linewidth = 0.2, color = 'black')
m.drawparallels(np.arange(int(-60),int(90),30),labels=[1,0,0,0], linewidth = 0, ax = ax)
m.drawparallels(np.array([0]),labels=[0,0,0,0], linewidth = 1, ax = ax)
m.drawmeridians(np.arange(int(0),int(360),60),labels=[0,0,0,1], linewidth = 0, ax = ax)
# m.fillcontinents(color='grey', ax = ax)
ax.text(-0.1, 1.1, 'd', transform=ax.transAxes, fontsize=10, weight='bold')
ax.set_title('P-E (mm d$^{-1}$)', fontsize = 10, loc = 'left')  
# colorbar
axins = inset_axes(ax,
                   width="3.5%",  # width = 10% of parent_bbox width
                   height="100%",  # height : 50%
                   loc='lower left',
                   bbox_to_anchor=(1.02, 0.0, 1.0, 1.0),
                   bbox_transform=ax.transAxes,
                   borderpad=0,
                   )
cbar = plt.colorbar(contours, ax=ax, cax = axins, ticks = np.arange(-1.5,2.,0.5))

######################## plot all-depth OHC ######################## 
# ax 0, 1 
ax = plt.subplot2grid((3, 2), (0, 0), colspan=1, rowspan = 3)
ax.spines['right'].set_color('blue')
ax.spines['left'].set_color('red')
#ax.spines['right'].set_linestyle((0, (5, 5)))
#ax.spines['left'].set_linestyle('-')
#ax = axes[0, 1]
# plot year 170 line
y = np.arange(0.,1.9,0.1)
x = 170* np.ones((np.shape(y)[0],1))
ax.plot(x, y, marker=None,linestyle = '-', color='grey')
yr_start = 110
yr_end = 215
years = np.arange(yr_start, yr_end, 1)
ax.plot(years, (netrad_toa_2co2std - netrad_toa_ctl)[5:], marker=None,linestyle = '-',color='blue',label = 'Net Rad (STD)')
ax.plot(years, (netrad_toa_2co2gl - netrad_toa_ctlgl)[5:], marker=None,linestyle = '--',color='blue',label = 'Net Rad (nudging)')
#ax.plot(years, netrad_toa_2co2atl - netrad_toa_ctlatl, marker=None,color='red',label = 'ATL nudging')
#ax.plot(years, netrad_toa_2co2nonatl - netrad_toa_ctlnonatl, marker=None,color='tab:brown',label = 'NonATL nudging')
ax.tick_params(axis="y", labelsize=10, colors = 'blue')
ax.tick_params(axis="x", labelsize=10)
#ax2.spines['right'].set_linestyle((0, (5, 5)))
#ax2.spines['left'].set_linestyle('-')
label = ax.set_ylabel('W m$^{-2}$', fontsize = 10)
label.set_color('blue')
ax.set_ylim((0.2, 1.37))
ax.set_xlim((100., 200))
ax.text(166, 0.3, 'CO$_{2}$ doubles', color = 'grey', rotation = 90)

ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
# plot OHC
yr_start = 100
yr_end = 200
years = np.arange(yr_start, yr_end, 1)
ax2.plot(years, chg_ohc_std[0:100], marker=None, linestyle = '-', color='red',label = 'OHC (STD)')
ax2.plot(years, chg_ohc_gl[0:100], marker=None, linestyle = '--', color='red',label = 'OHC (Nudging)')
#ax.plot(years, chg_ohc_atl[0:100], marker=None, color='red',label = 'ATL nudging')
#ax.plot(years, chg_ohc_nonatl[0:100], marker=None, color='tab:brown',label = 'NonATL nudging')
#ax.set_title('(c)', fontsize = 11, loc = 'left')
ax2.text(-0.12, 1.025, 'a', transform=ax.transAxes, fontsize=10, weight='bold')

#ax.legend(fontsize = 11)
ax.set_xlabel('Year', fontsize = 10)
label = ax2.set_ylabel('10$^{24}$ J', fontsize = 10)
label.set_color('red')
ax2.tick_params(axis="y", labelsize=10, colors='red')
ax2.set_ylim((0., 1.4))
ax2.spines['right'].set_color('red')
ax2.spines['left'].set_color('blue')

ax2.text(106, 1.1, 'OHC', color = 'red', fontsize = 10)
ax2.text(106, 1.2, 'TOA Net Radiation', color = 'blue', fontsize = 10)

plt.subplots_adjust(wspace= 0.10)
plt.subplots_adjust(hspace= 0.34)

plt.savefig('FIG2.pdf', bbox_inches='tight')







