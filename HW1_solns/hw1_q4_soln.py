# This program solves problem 4 of homework 1 for Valpo's MET 430.
# It creates a time height plot of lapse rate from the Valparaiso Univ. profiling radiometer.
# Christopher Phillips

##### START OPTIONS #####

# Location of radiometer file
rfile = 'radiometer_2025-01-24_00-04-09_lv2.nc'

# Location to save plot
spath = 'hw1_q4_lapse_rate.png'

#####  END OPTIONS  #####

# Import required modules
import matplotlib.pyplot as pp
import netCDF4 as nc
import numpy as np

# Read in the radiometer file
# and extract variables
fn = nc.Dataset(rfile, 'r')
z = fn.variables['height'][:]*1000.0 # km -> m, shape z
time = fn.variables['time'][:] # hours since midnight UTC, shape t
temp = fn.variables['TEMP'][:] # K, shape t, z
fn.close()

# Compute centered difference
# Use forward and reverse differencing on grid edges
dTdz = np.zeros(temp.shape)

# The centered difference
for i in range(1, temp.shape[1]-1):
    dTdz[:,i] = (temp[:,i+1]-temp[:,i-1])/(z[i+1]-z[i-1])

# Forard and reverse differences
dTdz[:,0] = (temp[:,1]-temp[:,0])/(z[1]-z[0]) # Forward
dTdz[:,-1] = (temp[:,-1]-temp[:,-2])/(z[-1]-z[-2]) # Reverse

# Make the plot
fig, (ax1, ax2) = pp.subplots(ncols=2, constrained_layout=True, figsize=(8,6))

# First axis is temperature gradient through whole atmosphere
cont = ax1.contourf(time, z, dTdz.transpose(), cmap='bwr')
ax1.set_xlabel('Time\n(hours since midnight UTC)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Height (m)', fontsize=14, fontweight='bold')

# Second axis is temperature gradient through just the lower atmosphere
ax2.contourf(time, z, dTdz.transpose(), cmap='bwr')
ax2.set_xlabel('Time\n(hours since midnight UTC)', fontsize=14, fontweight='bold')
ax2.set_ylim(0, 2000)

# Colorbar can be shared since both plots have the same data limits
cb = fig.colorbar(cont, ax=[ax1,ax2], orientation='vertical')
cb.set_label('dT/dz', fontsize=14, fontweight='bold')

pp.savefig(spath)
pp.close()