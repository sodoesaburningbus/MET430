# This script makes the plots for the shallow water model output
# Christopher Phillips

# Import modules
import matplotlib.pyplot as pp
import numpy as np

# Location of the output file
fpath = './sw_output.csv'

# Number of grid points in the model
nx = 100

# Lists to hold the data
time = []
u = []
z = []
x =  []

# Read in the output file
fn = open(fpath)
for line in fn:

    # Split the data
    dummy = line.split(',')
    time.append(dummy[0])
    u.append(dummy[1:nx+1])
    z.append(dummy[nx+1:2*nx+1])
    x.append(dummy[2*nx+1:3*nx+1])

# Convert to numpy arrays
u = np.array(u)
z = np.array(z)
x = np.array(x)

# Make the plots
for i in u.shape[0]:

    fig, ax = pp.subplots()
    ax.set_title(f'Time {time:.2f} s', loc='left', ha='left', fontsize=14, fontweight='bold')