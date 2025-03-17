### This program integrates the 1D shallow water equations
### within a pond (i.e. hard boundaries on the edges)
### Christopher Phillips

##### START OPTIONS #####

# Location to save the output file
spath = './sw_output_100m.nc'

# Grid information
dx = 100.0 # Grid spacing (meters)
nx = 500 # Number of grid points

# Time information
dt = 5.0 # Time stepping (seconds)
output_interval = 10.0 # Time step interval (seconds)
tmax = 3600.0*3.0 # Max time of integration (seconds)

# Other settings
g = 9.81 # Gravity (m/s2)
H = 1000.0 # Base height of the water surface (the rest state, meters)
eta0 = 400.0 # Initial perturbation depth
u0 = 0.0 # Initial current (m/s)
v0 = 0.0 # Initial v current (m/s)
gamma = 0.2 # A-R filter gamma
f = 10**-4 # Coriolis parameter (1/s)

#####  END OPTIONS  #####

# Import required modules
import netCDF4 as nc
import numpy as np

# Create output file to which model output will be saved
fn = nc.Dataset(spath, 'w')

# NetCDF4 dimensions
fn.createDimension('time', None)  # Unlimited time dimension
fn.createDimension('x', nx)  # Fixed x dimension

# NetCDF4 variables
times_out = fn.createVariable('time', np.float32, ('time',))
x_out = fn.createVariable('X', np.float32, ('x',))
Z_out = fn.createVariable('Z', np.float32, ('time', 'x'))
U_out = fn.createVariable('U', np.float32, ('time', 'x'))
V_out = fn.createVariable('V', np.float32, ('time', 'x'))

# Variable attributes
times_out.units = 's'
x_out.units = 'm'
Z_out.units = 'm'
U_out.units = 'm s-1'
V_out.units = 'm s-1'

### Helper functions for computing spatial derivatives
### Numpy gradient is faster, but this is illustrative
def spatial_difference(A):

    # Interior points
    dA = np.zeros(A.size)
    for i in range(1,A.size-1):
        dA[i] = A[i+1]-A[i-1]
    
    # Edges
    dA[0] = A[1]-0
    dA[-1] = 0.0-A[-2]

    return dA

# Initialize the grid
x = np.arange(nx)*dx
u_now = np.ones(x.size)*u0
v_now = np.ones(x.size)*v0
z_now = np.zeros(x.size)
z_now[int(x.size/2)-10:int(x.size/2)+10] += eta0*np.sin(np.linspace(0,np.pi,20))

# Do a single-forward step with Euler's
step = dt/dx
time = 0
du = np.zeros(x.size)
dv = np.zeros(x.size)
dz = np.zeros(z_now.size)

# Write out initial conditions
oc = 0 # Counter for tracking the number of output steps
x_out[:] = x
times_out[oc] = 0.0
U_out[oc,:] = u_now
V_out[oc,:] = v_now
Z_out[oc,:] = z_now+H

# Compute gradients
# See how we're using upwind for computing the gradient.
for i in range(1,du.size-1):
    if (u_now[i] < 0):
        du[i] = u_now[i+1]-u_now[i]
        dv[i] = v_now[i+1]-v_now[i]
        dz[i] = z_now[i+1]-z_now[i]
    else:
        du[i] = u_now[i]-u_now[i-1]
        dv[i] = v_now[i]-v_now[i-1]
        dz[i] = z_now[i]-z_now[i-1]

# Handle the edges (hard wall)
du[0] = u_now[1]-0.0
du[-1] = 0.0-u_now[-2]
dv[0] = v_now[1]-0.0
dv[-1] = 0.0-v_now[-2]
dz[0] = z_now[1]-0.0
dz[-1] = 0.0-z_now[-2]

# Do the time step
u_new = u_now-step*(g*dz+u_now*du-f*v_now*dx)
v_new = v_now-step*(u_now*dv+f*u_now*dx)
z_new = z_now-step*((H+z_now)*du+u_now*dz)

# Cycle the variables
u_old = u_now
v_old = v_now
z_old = z_now
u_now = u_new
v_now = v_new
z_now = z_new
time += dt

# Convert output time interval to an integer time step
output_step = int(output_interval/dt)

# Now do the leap-frog iterations
tc = 1  # Number of time steps done
while (time <= tmax):

    # Check if an output time and output if so
    if (tc % output_step == 0):
        oc += 1
        times_out[oc] = time
        U_out[oc,:] = u_now
        V_out[oc,:] = v_now
        Z_out[oc,:] = z_now+H

    # Compute gradients
    du = spatial_difference(u_now)
    dv = spatial_difference(v_now)
    dz = spatial_difference(z_now)

    # Do the time step
    u_new = u_now-step*(g*dz+u_now*du-f*v_now*dx)
    v_new = v_now-step*(u_now*dv+f*u_now*dx)
    z_new = z_old-step*((H+z_now)*du+u_now*dz)

    # Do the filters
    u_now = u_now + gamma*(u_new-2.0*u_now+u_old)
    v_now = v_now + gamma*(v_new-2.0*v_now+v_old)

    # Cycle variables
    u_old = u_now
    v_old = v_now
    z_old = z_now
    u_now = u_new
    v_now = v_new
    z_now = z_new
    time += dt
    tc += 1
    
# Close the output file
fn.close()