# This program solves problem 5 of homework 1 for Valpo's MET 430.
# It interpolates a sounding to a grid and plots it versus the original
# to ensure the interoplation worked properly.
# Christopher Phillips

##### START OPTIONS #####

# Location of radiometer file
ifile = 'input_sounding.txt'

# Location to save plot
spath = 'hw1_q5_soundings.png'

# Bounds of plot (m)
zlims = (0, 12000)

#####  END OPTIONS  #####

# Import required modules
import matplotlib.pyplot as pp
import numpy as np

# Dictionary to hold lists of the variables
vars = ['pres', 'height', 'temp', 'dewpt', 'rh', 'mixr', 'wdir', 'wspd', 'theta', 'thetaE', 'thetaV'] # Same order as columns in input_sounding.txt
sounding = {}
for v in vars:
    sounding[v] = []

# Read in the data file
fn = open(ifile, 'r')
for line in fn:

    # Split into columns
    cols = line.split()

    # Pull the variables and strip of any extra whitespace or end of line characters
    for i, v in enumerate(vars):
        sounding[v].append(cols[i].strip())

fn.close()

# Convert lists in sounding to numpy arrays
# Change data type from strings to floats at the same time
for v in vars:
    sounding[v] = np.array(sounding[v], dtype='float')

    # Subtract first height to make sounding above ground level
    # This is because model is from 0 to 12 km but soundings measure from above sea level.
    if (v == 'height'):
        sounding[v] -= sounding[v][0]

# Once all variables are in numpy arrays, can interpolate
isounding = {} # Interpolated sounding dictionary
iheights = np.arange(0.0, 12050.0, 50.0)
for v in vars:
    isounding[v] = np.interp(iheights, sounding['height'], sounding[v])

# Make the comparison plot
# Just check a few variables
fig, (ax1, ax2, ax3) = pp.subplots(ncols=3, figsize=(10,8), constrained_layout=True)

# Temperature
ax1.plot(sounding['temp'], sounding['height'], color='black', label='Original')
ax1.plot(isounding['temp'], isounding['height'], color='dodgerblue', label='Interpolated')
ax1.grid()
ax1.set_xlabel("Temperature ('C)", fontsize=14, fontweight='bold')
ax1.set_ylabel("Height (m)", fontsize=14, fontweight='bold')
ax1.legend()
ax1.set_ylim(zlims)

# Dewpoint
ax2.plot(sounding['dewpt'], sounding['height'], color='black', label='Original')
ax2.plot(isounding['dewpt'], isounding['height'], color='dodgerblue', label='Interpolated')
ax2.grid()
ax2.set_xlabel("Dewpoint ('C)", fontsize=14, fontweight='bold')
ax2.set_ylim(zlims)

# Wind Speed
ax3.plot(sounding['wspd'], sounding['height'], color='black', label='Original')
ax3.plot(isounding['wspd'], isounding['height'], color='dodgerblue', label='Interpolated')
ax3.grid()
ax3.set_xlabel("Wind Speed (kts)", fontsize=14, fontweight='bold')
ax3.set_ylim(zlims)

pp.savefig(spath)
pp.close()


