import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from openpmd_viewer import OpenPMDTimeSeries, ParticleTracker
from scipy.constants import c, eV

mpl.rcParams.update({'font.size': 12})

# some useful constants from the input file
MeV = 1e6*eV
d1 = 0.1  
d2 = 0.19
d3 = 0.12
d4 = 0.2
xmin = -0.08
xmax =  0.08
ymin = -0.1
ymax =  0.3
zmax = d1+d2+d3+d4
zmin= -1e-3  

# get the data 
series = OpenPMDTimeSeries('./diags/diag1')

it = series.iterations
time = series.t
N_iterations = len(it)
species = series.avail_species
N_species = len(species)
N_particles = np.empty( N_species ) 

# arrays of all the trajectories 
x_trajectories = [None] * N_species
y_trajectories = [None] * N_species
z_trajectories = [None] * N_species

# array of initial energies 
E0 = [None] * N_species

for s in range( N_species ):

    # track all the particles of species s
    pt = ParticleTracker(series, iteration=it[0],  species=species[s], preserve_particle_index=True)
    ids0 = series.get_particle( ['id'], iteration=0, species=species[s], select=pt)
    
    N_particles[s] = pt.N_selected
    x_trajectories[s] = np.full( ( N_iterations, int(N_particles[s]) ), np.nan )
    y_trajectories[s] = np.full( ( N_iterations, int(N_particles[s]) ), np.nan )
    z_trajectories[s] = np.full( ( N_iterations, int(N_particles[s]) ), np.nan ) 
    E0[s] = np.full( (int(N_particles[s])), np.nan ) 
    
    # get the trajectories 
    for i in range( N_iterations ):
        x, y, z, ids = series.get_particle( ['x', 'y', 'z', 'id'], select=pt, iteration=it[i], species=species[s] )
        
        if np.shape(x) == (0,):
            continue 
            
        indeces = (ids==ids0)[0]                
        x_trajectories[s][i, indeces] = x[:]
        y_trajectories[s][i, indeces] = y[:]
        z_trajectories[s][i, indeces] = z[:]
        
        # get the initial energy
        if i==0:
           uz, m, ids = series.get_particle( ['uz', 'mass', 'id'], select=pt, iteration=it[i], species=species[s] )
           E0[s][indeces] = 0.5*m[:]*(uz[:]*c)**2 / MeV

# canvas 
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,10), dpi=300, subplot_kw={'projection': '3d'})

# draw relevant planes
X =  [xmin, xmax, xmax, xmin]
Y =  [ymin, ymin, ymax, ymax]
zE = [d1, d1, d1, d1]
zEB = [d1+d2, d1+d2, d1+d2, d1+d2]
zB = [d1+d2+d3, d1+d2+d3, d1+d2+d3, d1+d2+d3]
z0 = [0, 0, 0, 0]
z_detector = [zmax, zmax, zmax, zmax]
verts = [list(zip(X, Y, z0))]
ax.add_collection3d(Poly3DCollection(verts, linewidths=1, color='blue', alpha=0.1))
verts = [list(zip(X, Y, zE))]
ax.add_collection3d(Poly3DCollection(verts, linewidths=1, color='red', alpha=0.05))
verts = [list(zip(X, Y, zEB))]
ax.add_collection3d(Poly3DCollection(verts, linewidths=1, color='red', alpha=0.05))
verts = [list(zip(X, Y, zB))]
ax.add_collection3d(Poly3DCollection(verts, linewidths=1, color='red', alpha=0.05))
verts = [list(zip(X, Y, z_detector))]
ax.add_collection3d(Poly3DCollection(verts, linewidths=1, color='lime', alpha=0.1))

# plot trajectories in colors ordered by initial condition
cmap = ['winter', 'cool']
for s in range( N_species ):
    print(species[s])

    colors = plt.get_cmap(cmap[s])(np.linspace(0,1,int(N_particles[s])))
    indeces = E0[s].argsort()

    E0[s] = E0[s][indeces]
    x_trajectories[s][:,:] = x_trajectories[s][:,indeces] 
    y_trajectories[s][:,:] = y_trajectories[s][:,indeces] 
    z_trajectories[s][:,:] = z_trajectories[s][:,indeces] 

    for p in range( int(N_particles[s]) ):             
        j = indeces[p]
        ax.plot(x_trajectories[s][:,j], y_trajectories[s][:,j], z_trajectories[s][:,j], color=colors[j], lw=1) 

# plot initial position       
ax.scatter(0, 0, 0, marker='*', color='black', s=80)

ax.view_init(vertical_axis='x', azim=20, elev=25)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
ax.set_aspect('equal', 'box')
ax.set_box_aspect((1, 1, 1) , zoom=1.1)

# save
fig.savefig(f'trajs2species.png', dpi=300)
plt.close()
