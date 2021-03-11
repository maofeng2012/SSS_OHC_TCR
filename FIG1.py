## load libraries
import numpy as np
import matplotlib.pyplot as plt

########################## running mean funtion ######################
def moving_average(a, n=20) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

################################################################################
#################### Read global mean surf temp time series ####################
################################################################################
yr_start = 101
yr_end = 219
yrs = np.arange(yr_start, yr_end+1, 1)

######################## read Surf Temp data ######################## 
surfT_ctl = np.loadtxt('t_surf_ctl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
surfT_ctlgl = np.loadtxt('t_surf_ctlgl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
surfT_2co2std = np.loadtxt('t_surf_2CO2std_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')
surfT_2co2gl = np.loadtxt('t_surf_2CO2gl_'+np.str(yr_start)+'_'+np.str(yr_end)+'_series_raw.txt', delimiter=' ')

surfT_ctl = moving_average(surfT_ctl, n=20)
surfT_ctlgl = moving_average(surfT_ctlgl, n=20)

surfT_2co2std = moving_average(surfT_2co2std, n=20)
surfT_2co2gl = moving_average(surfT_2co2gl, n=20)

#######################################################################################            
######################################### plot ########################################
#######################################################################################
fig, axes = plt.subplots(1, 1, figsize = (5,4))
# ax 0
ax = axes
y = np.arange(0,2.9,0.1)
x = 170* np.ones((np.shape(y)[0],1))
ax.plot(x, y, marker=None,linestyle = '-', color='grey')
# plot surf temp
yr_start = 110
yr_end = 210
years = np.arange(yr_start, yr_end, 1)
ax.plot(years,(surfT_2co2std - surfT_ctl), marker=None,color='blue',linestyle = '-', label = 'STD')
ax.plot(years,(surfT_2co2gl - surfT_ctlgl), marker=None,color='blue',linestyle = '--', label = 'GL nudging')
#ax.legend(fontsize = 11, loc='lower left', bbox_to_anchor=(-0.1, -0.38), ncol = 4)
#ax.set_title('Surf Temp', fontsize = 11, loc = 'left')
ax.set_xlabel('Year', fontsize = 11)
ax.set_ylabel('K', fontsize = 11)
ax.tick_params(axis="x", labelsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.set_ylim((0., 2.5))
ax.set_xlim((100., 200))

tcr_2co2std = np.str(np.round((surfT_2co2std - np.nanmean(surfT_ctl))[years == 170][0], 1)) 
tcr_2co2gl = np.str(np.round((surfT_2co2gl - np.nanmean(surfT_ctlgl))[years == 170][0], 1)) 

ax.text(105, 2.1, 'TCR (STD) = '+tcr_2co2std+' K', color = 'black')
ax.text(105, 1.8, 'TCR (fixed-SSS-GL) = '+tcr_2co2gl+' K', color = 'black')
ax.text(163.6, 0.3, 'CO$_{2}$ doubles', color = 'grey', rotation = 90)
#ax.text(-0.17, 1.03, 'a', transform=ax.transAxes, fontsize=11, weight='bold')
plt.subplots_adjust(wspace= 0.27)
#plt.subplots_adjust(hspace= 0.5)
plt.savefig('FIG1_Surf_Temp_changes_{2CO2STD-2CO2GL}.pdf', bbox_inches='tight')
















