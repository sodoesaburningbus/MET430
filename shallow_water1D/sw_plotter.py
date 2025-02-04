# This script makes the plots for the shallow water model output
# Christopher Phillips

# Import modules
import matplotlib.pyplot as pp
import numpy as np

##### START OPTIONS ######

# Location of the output file
fpath = './sw_output.csv'

# Location to save outputs
sdir = 'images/'

# Interval between images
interval = 50

# Number of grid points in the model
nx = 500

#####  END OPTIONS  #####

# Lists to hold the data
time = []
u = []
z = []
x =  []

# Read in the output file
fn = list(open(fpath))
for line in fn[1:]:

    # Split the data
    dummy = line.split(',')
    time.append(dummy[0])
    u.append(dummy[1:nx+1])
    z.append(dummy[nx+1:2*nx+1])
    x.append(dummy[2*nx+1:3*nx+1])

# Convert to numpy arrays
time = np.array(time, dtype='float')
u = np.array(u, dtype='float')
z = np.array(z, dtype='float')
x = np.array(x, dtype='float')

# Make the plots
fs = 14
fw = 'bold'
for i in range(0, u.shape[0])[::interval]:

    fig, (ax1, ax2) = pp.subplots(nrows=2, constrained_layout=True)
    ax1.set_title(f'Time {time[i]:.2f} s', loc='left', ha='left', fontsize=fs, fontweight=fw)

    ax1.plot(x[i,:], z[i,:], color='black')
    ax1.set_ylabel('Depth (m)', fontsize=fs, fontweight=fw)
    ax1.grid()
    ax1.set_ylim(500, 1500)

    ax2.plot(x[i,:], u[i,:], color='black')
    ax2.set_xlabel('X (m)', fontsize=fs, fontweight=fw)
    ax2.set_ylabel('U (m/s)', fontsize=fs, fontweight=fw)
    ax2.grid()
    ax2.set_ylim(-35, 35)

    pp.savefig(f'{sdir}/frame_{i:04d}.png')
    pp.close()
