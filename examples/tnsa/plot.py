
import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import use, rcParams, cm  
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm, Normalize
from openpmd_viewer import OpenPMDTimeSeries
from scipy.constants import e as q_e, m_e, c, pi, epsilon_0, m_p

use('Agg')
rcParams.update({'font.size': 12})

plot_dir = 'plots'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
   
laser_wavelength=0.8e-6
n_c = epsilon_0 * m_e * (2*pi*c)**2 /(laser_wavelength*q_e)**2 

series = OpenPMDTimeSeries('./diags/diag1')

iterations = series.iterations
time = series.t
N_iterations = len(iterations)

for i, it in enumerate( iterations ):
    print(i)
    Ex,info = series.get_field(field='E', coord='x', iteration=it, plot=False)
    rhoe,info = series.get_field(field='rho_ele_targ', iteration=it, plot=False)
    rhoi,info = series.get_field(field='rho_ion_cont', iteration=it, plot=False)
        
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,8), dpi=300)
    
    imshow_extent = [info.z.min(), info.z.max(), info.x.min(), info.x.max()]
    n = (-rhoe+rhoi)/q_e/n_c
    im=ax.imshow(np.transpose(n),alpha=1.0, vmin=0.1, vmax=40,  extent=imshow_extent, cmap=cm.Greys)
    
    data = np.transpose(Ex)/(4e12)
    emax = 5
    emin = -5
    alphas = Normalize(0,emax, clip=True)(np.abs(data))
    alphas = np.clip(alphas, 0.1, 0.4) 
    colors = Normalize(emin, emax)(data)
    cmap = LinearSegmentedColormap.from_list("mycmap",  ["red", "white" ,"blue"]) # cm.PiYG
    colors = cmap(colors)
    colors[..., -1] = alphas
    im=ax.imshow(colors,extent=imshow_extent, cmap=cmap)

    ax.set_xlabel('z [m]')
    ax.set_ylabel('x [m]')
    
    plt.tight_layout()
    figname = os.path.join(plot_dir, f'T_{i:04d}.png')
    fig.savefig(figname, dpi=300)

    plt.close()
    
    
