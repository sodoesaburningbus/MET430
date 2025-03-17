### This program models a silpe oscillator (i.e. spring)
### and compares the results of two finite differencign schemes to
### the analytical solution.
###
### Christopher Phillips
### Valparaiso University

##### START OPTIONS #####

# Time step
dt = 1.0 # seconds

# Maximum integration time
tmax = 300.0 # seconds

# Initial spring perturbation
# I.E. how far is the spring pulled out to start
x0 = 0.50 # meters

# Spring constant
k = 0.05 # kg/s2

# Bob mass (the mass of the oscillator)
m = 1.0 # kg

# Location to save the outputs
spath_plot = './spring_plot_1s.png'

# Font options
fs = 14
fw = 'bold'

#####  END OPTIONS  #####

# Switch plotting backend
import matplotlib
matplotlib.use('agg')

# Import required modules
import matplotlib.pyplot as pp
import numpy as np

### First, make the analytical solution
time_a = np.arange(0, tmax+dt, dt)
x_a = x0*np.cos(np.sqrt(k/m)*time_a)
v_a = -x0*np.sqrt(k/m)*np.sin(np.sqrt(k/m)*time_a)


### Euler's Method
x_e = [x0]
v_e = [0]
time_e = [0]

# Can't use the proper x formulation on first step,
# so approximate with an Euler's step on the velocity too
v_new = v_e[-1]-dt*k/m*x_e[-1]
x_new = x_e[-1]+v_e[-1]*dt
x_e.append(x_new)
v_e.append(v_new)
time_e.append(time_e[-1]+dt)

# Now do the loop
while time_e[-1] < tmax:

    # Do the time step
    x_new = 2.0*x_e[-1]-x_e[-2]-(dt**2)*k/m*x_e[-1]
    v_new = v_e[-1]-dt*k/m*x_e[-1]

    # Store variables
    x_e.append(x_new)
    v_e.append(v_new)
    time_e.append(time_e[-1]+dt)


### Adams-Bashforth Scheme
# Need to do the first two time steps using Euler's method
x_ab = [x0]
v_ab = [0]
time_ab = [0]

for i in range(2):
    v_new = v_ab[-1]-dt*k/m*x_ab[-1]
    x_new = x_ab[-1]+v_ab[-1]*dt
    x_ab.append(x_new)
    v_ab.append(v_new)
    time_ab.append(time_ab[-1]+dt)

# Now go into Adams-Bashforth
step = dt/12.0*k/m
while time_ab[-1] < tmax:

    # Time steps
    v_new = v_ab[-1]-step*(23.0*x_ab[-1]-16.0*x_ab[-2]+5.0*x_ab[-3])
    x_new = 2.0*x_ab[-1]-x_ab[-2]-(dt**2)*k/m*x_ab[-1]

    # Store variables
    x_ab.append(x_new)
    v_ab.append(v_new)
    time_ab.append(time_ab[-1]+dt)

### Plot the results
fig, axes = pp.subplots(nrows=3, constrained_layout=True, figsize=(12,8))

# Analytical solutions
axes[0].set_title('Analytical', loc='left', ha='left', fontsize=12, fontweight='bold')
ax0b = axes[0].twinx()

axes[0].plot(time_a, x_a, color='black')
ax0b.plot(time_a, v_a, color='dodgerblue')

axes[0].set_ylabel('Position (m)', fontsize=fs, fontweight=fw)
ax0b.set_ylabel('Velocity (m/s)', fontsize=fs, fontweight=fw, color='dodgerblue')

# Euler's Method
axes[1].set_title("Euler's", loc='left', ha='left', fontsize=12, fontweight='bold')
ax1b = axes[1].twinx()

axes[1].plot(time_e, x_e, color='black')
ax1b.plot(time_e, v_e, color='dodgerblue')

axes[1].set_ylabel('Position (m)', fontsize=fs, fontweight=fw)
ax1b.set_ylabel('Velocity (m/s)', fontsize=fs, fontweight=fw, color='dodgerblue')

# Adams-Bashforth Scheme
axes[2].set_title('A-B', loc='left', ha='left', fontsize=12, fontweight='bold')
ax2b = axes[2].twinx()

axes[2].plot(time_ab, x_ab, color='black')
ax2b.plot(time_ab, v_ab, color='dodgerblue')

axes[2].set_ylabel('Position (m)', fontsize=fs, fontweight=fw)
ax2b.set_ylabel('Velocity (m/s)', fontsize=fs, fontweight=fw, color='dodgerblue')

axes[2].set_xlabel('Time (s)')

for ax in axes:
    ax.grid()
    ax.set_xlim(0, tmax)

pp.savefig(spath_plot)
pp.close()