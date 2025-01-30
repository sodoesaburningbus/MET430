### This program integrates the 1D shallow water equations
### within a pond (i.e. hard boundaries on the edges)
### Christopher Phillips

##### START OPTIONS #####

# Location to save the output file
spath = './sw_output.csv'

# Grid information
dx = 10000.0 # Grid spacing (meters)
nx = 500 # Number of grid points

# Time information
dt = 5.0 # Time stepping (seconds)
output_interval = 10.0 # Time step interval (seconds)
tmax = 10800.0 # Max time of integration (seconds)

# Other settings
g = 9.81 # Gravity (m/s2)
H = 1000.0 # Base height of the water surface (the rest state, meters)
eta0 = 500 # Initial perturbation depth
u0 = 0.0 # Initial current (m/s)
gamma = 0.2 # A-R filter gamma

#####  END OPTIONS  #####

# Import required modules
import numpy as np

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
z_now = np.zeros(x.size)
z_now[int(x.size/2)-10:int(x.size/2)+10] += eta0*np.sin(np.linspace(0,np.pi,20))

# Do a single-forward step with Euler's
step = dt/dx
time = 0
du = np.zeros(x.size)
dz = np.zeros(z_now.size)

# Write out initial conditions
fn = open(spath,'w')
fn.write('Time'+(',u{}'*u_now.size).format(*range(u_now.size))
         +(',z{}'*z_now.size).format(*range(z_now.size))
         +(',x{}'*x.size).format(*range(x.size)))
fn.write('\n{:.3f}'.format(time)+(',{:.3f}'*u_now.size).format(*u_now)
         +(',{:.3f}'*z_now.size).format(*(z_now+H))
         +(',{:.3f}'*x.size).format(*x))

# Compute gradients
for i in range(1,du.size-1):
    if (u_now[i] < 0):
        du[i] = u_now[i+1]-u_now[i]
        dz[i] = z_now[i+1]-z_now[i]
    else:
        du[i] = u_now[i]-u_now[i-1]
        dz[i] = z_now[i]-z_now[i-1]


dz[0] = z_now[1]-0.0
dz[-1] = 0.0-z_now[-2]

# Do the time step
u_new = u_now+step*(g*dz-u_now*du)
z_new = z_now-step*((H+z_now)*du+u_now*dz)

# Cycle the variables
u_old = u_now
z_old = z_now
u_now = u_new
z_now = z_new
time += dt

# Now do the leap-frog iterations
while (time <= tmax):

    # Check if an output time
    fn.write('\n{:.3f}'.format(time)+(',{:.3f}'*u_now.size).format(*u_now)
        +(',{:.3f}'*z_now.size).format(*(z_now+H))
        +(',{:.3f}'*x.size).format(*x))

    """Euler's
    # Compute gradients
    for i in range(1,du.size-1):
        if (u_now[i] < 0):
            du[i] = u_now[i+1]-u_now[i]
            dz[i] = z_now[i+1]-z_now[i]
        else:
            du[i] = u_now[i]-u_now[i-1]
            dz[i] = z_now[i]-z_now[i-1]

    du[0] = u_now[0]-0
    du[-1] = 0.0-u_now[-1]
    dz[0] = z_now[0]-0.0
    dz[-1] = 0.0-z_now[-1]

    # Do the time step
    u_new = u_now+step*(g*dz-u_now*du)
    z_new = z_now-step*((H+z_now)*du+u_now*dz)
    """
    # Compute gradients
    du = spatial_difference(u_now)
    dz = spatial_difference(z_now)

    # Do the time step
    u_new = u_old-step*(g*dz-u_now*du)
    z_new = z_old-step*((H+z_now)*du+u_now*dz)

    # Do the filters
    u_now = u_now + gamma*(u_new-2.0*u_now+u_old)

    # Cycle variables
    u_old = u_now
    z_old = z_now
    u_now = u_new
    z_now = z_new
    time += dt
    
# Close the output file
fn.close()