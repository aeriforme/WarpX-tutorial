import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from openpmd_viewer import OpenPMDTimeSeries
from openpmd_viewer import ParticleTracker
from scipy.constants import eV, m_p, c 

mpl.use('Agg')
mpl.rcParams.update({'font.size': 18})

MeV=1e6*eV
Emax_hydrogen1_1 = 50*MeV
Emax_carbon12_6 = 10*MeV

# Workaround to get the correct particles' pids 
def flipHiBit(id):
    return id ^ (1 << 63)

def getParticleNumber(id):
    return (id ^ (1 << 63)) >> np.uint64(24)

def getScrapedParticleNumber(id):
    return id >> np.uint64(24)

series = OpenPMDTimeSeries('./diags/diag2/particles_at_zhi/')
series2 = OpenPMDTimeSeries('./diags/diag1/')

it = series.iterations
time = series.t
N_iterations = len(it)
species = series.avail_species
N_species = len(species)
cmap = ['rainbow', 'rainbow']

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,8), dpi=300)

for s in range(N_species):
    print(species[s])

    X, Y, E = [], [], []
    for i in range(N_iterations):
        x,y,z,ids = series.get_particle( ['x','y','z','id'], 
                                         iteration=it[i], 
                                         species=species[s], 
                                         plot=False )
        correct_ids =  getParticleNumber(ids)    
        
        uz0, ids0, m = series2.get_particle( ['uz','id','mass'],
                                             iteration=series2.iterations[0], 
                                             species=species[s], 
                                             plot=False )
        correct_ids0 = getScrapedParticleNumber(ids0)    
            
        indeces = np.where(np.in1d(correct_ids0,correct_ids))[0]

        E = np.append(E, 0.5*m[indeces]*(uz0[indeces]*c)**2/MeV)
        X = np.append(X, x) 
        Y = np.append(Y, y)
        
    sorted_indeces = np.argsort(E)   
    im = ax.scatter(X[sorted_indeces], Y[sorted_indeces], c=E[sorted_indeces], cmap=cmap[s], vmin=5, vmax=50)

fig.colorbar(im, label='E [MeV]')        
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_aspect('equal')

plt.tight_layout()
fig.savefig(f'screen.png', dpi=300)
plt.close()

