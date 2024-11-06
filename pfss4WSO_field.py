import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sunpy.map
from astropy.coordinates import SkyCoord
from matplotlib.colors import SymLogNorm

import pfsspy
from pfsspy import coords, tracing
from pfsspy.sample_data import get_gong_map

# Referring to raw data, from CR1642 to CR2287
cr_beg = 2209
cr_end = 2287
data_dir = 'E:/Research/Work/flux_transport_on_ss/Br_interp/data/'
save_dir = 'E:/Research/Work/flux_transport_on_ss/Br_ss/'
save_or_not = 0

# Importing GONG .fits as the Sunpy.map's meta
GONG_file = 'E:/Research/Data/GONG/fits/mrzqs_c2100.fits'
GONG_map = sunpy.map.Map(GONG_file)

# Defining plot settings
def set_axes_lims(ax):
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    
# Performing PFSS model
for i_cr in range(cr_beg, cr_end+1):
    
    # Importing WSO .csv (already interpolated to GONG grid) as the Sunpy.map's data
    WSO_file = 'CR' + str(i_cr) + '_interp.csv'
    df_Br_interp = pd.read_csv(data_dir+WSO_file, header=None)
    Br_interp  = df_Br_interp.values

    # Removing the mean, so that curl{B} = 0
    GONG_map = sunpy.map.Map(Br_interp - np.mean(Br_interp), GONG_map.meta)

    # Setting the grid of model
    nrho = 35
    rss = 2.5
    pfss_in = pfsspy.Input(GONG_map, nrho, rss)

    # Plotting the input synoptic map
    m = pfss_in.map
    fig = plt.figure()
    ax = plt.subplot(projection=m)
    m.plot()
    plt.colorbar()
    ax.set_title('Input photospheric field')
    set_axes_lims(ax)

    # Performing the PFSS model
    pfss_out = pfsspy.pfss(pfss_in)
    Br_ss = pfss_out.source_surface_br
    
    # Saving the source surface field
    if save_or_not == 1:
        save_file = 'CR' + str(i_cr) + '_ss.csv'
        df_Br_ss = pd.DataFrame(Br_ss.data)
        df_Br_ss.to_csv(save_dir+'data/'+save_file, header=False, index=False)

    # Plotting the source surface map
    fig = plt.figure()
    ax = plt.subplot(projection=Br_ss)

    Br_ss.plot()
    ax.plot_coord(pfss_out.source_surface_pils[0])
    plt.colorbar()
    ax.set_title('Source surface magnetic field')
    set_axes_lims(ax)
    
    if save_or_not == 1: 
        save_file = 'CR' + str(i_cr) + '_ss.png'
        plt.savefig(save_dir+'figures/'+save_file)

    # plt.show()
    plt.close()

    # Ploting magnetic fields on the interlayer
    interlayer = 0

    if interlayer != 0:
        # Get the radial magnetic field at a given height
        ridx = 15
        Br_il = pfss_out.bc[0][:, :, ridx]
        # Create a sunpy Map object using output WCS
        Br_il = sunpy.map.Map(Br_il.T, pfss_out.source_surface_Br_il.wcs)
        # Get the radial coordinate
        r = np.exp(pfss_out.grid.rc[ridx])

        # Create the figure and axes
        fig = plt.figure()
        ax = plt.subplot(projection=Br_il)

        # Plot the source surface map
        Br_il.plot(cmap='RdBu')
        # Plot formatting
        plt.colorbar()
        ax.set_title('$B_{r}$ ' + f'at r={r:.2f}' + '$r_{\\odot}$')
        set_axes_lims(ax)

        plt.show()

db