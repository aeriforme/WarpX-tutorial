import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import use 
from matplotlib import rcParams 

from openpmd_viewer import OpenPMDTimeSeries

rcParams.update({'font.size': 12})

plot_dir = 'plots'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
   

series = OpenPMDTimeSeries('./diags/diag1')
iterations = series.iterations
time = series.t

for i, it in enumerate(iterations):
    print(i)

    
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10,8), dpi=300)

    # extract data and plot at once
    #rho, info = series.get_field(field='rho', slice_across='x', iteration=it, plot=True, cmap='rainbow')     

    # or 
    
    # extract the data 
    rho, info = series.get_field(field='rho', slice_across='x', iteration=it, plot=False) 
    
    # and then plot it 
    ax.imshow(rho, extent=info.imshow_extent, cmap='rainbow')
    ax.set_xlabel(info.axes[0]+' [m]')
    ax.set_ylabel(info.axes[1]+' [m]')
    ax.set_title(f"t = {series.t[i]} s, step = {it}")    
    
    figname = os.path.join(plot_dir, f'T_{i:04d}.png')
    fig.savefig(figname, dpi=300)
    plt.close()
    
