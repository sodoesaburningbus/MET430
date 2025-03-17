# This script makes the plots for the shallow water model output
# Christopher Phillips

# Switch matplotlib backend
import matplotlib
matplotlib.use('agg')

# Import modules
import matplotlib.pyplot as pp
import netCDF4 as nc
import numpy as np

##### START OPTIONS ######

# Experiment name
exp = '100m'

# Location of the output file
fpath = f'./sw_output_{exp}.nc'

# Location to save outputs
sdir = f'images_{exp}/'

#####  END OPTIONS  #####

# Read in the output file
fn = nc.Dataset(fpath)
u = fn.variables['U'][:,:]
v = fn.variables['V'][:,:]
z = fn.variables['Z'][:,:]
time = fn.variables['time'][:]
x = fn.variables['X'][:]

# Make the plots
fs = 14
fw = 'bold'
for i in range(0, u.shape[0]):

    fig, (ax1, ax2, ax3) = pp.subplots(nrows=3, constrained_layout=True)
    ax1.set_title(f'Time {time[i]:.2f} s', loc='left', ha='left', fontsize=fs, fontweight=fw)

    ax1.plot(x, z[i,:], color='black')
    ax1.set_ylabel('Depth (m)', fontsize=fs, fontweight=fw)
    ax1.grid()
    ax1.set_ylim(500, 1500)

    ax2.plot(x, u[i,:], color='black')
    ax2.set_xlabel('X (m)', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('U (m/s)', fontsize=fs, fontweight=fw)
    ax2.grid()
    ax2.set_ylim(-50, 50)

    ax3.plot(x, v[i,:], color='black')
    ax3.set_xlabel('X (m)', fontsize=fs, fontweight=fw)
    ax3.set_ylabel('V (m/s)', fontsize=fs, fontweight=fw)
    ax3.grid()
    ax3.set_ylim(-30, 30)

    pp.savefig(f'{sdir}/frame_{i:04d}.png')
    pp.close()
