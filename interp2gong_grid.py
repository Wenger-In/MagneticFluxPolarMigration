import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import griddata

# Referring to raw data, from CR1642 to CR2287
cr_beg = 1642
cr_end = 2287
data_dir = 'E:/Research/Data/WSO/field/'
save_dir = 'E:/Research/Work/flux_transport_on_ss/Br_interp/'
save_or_not = 0

# Importing WSO grid
WSO_lat_file = 'lat_arr.dat'
WSO_lon_file = 'lon_arr.dat'
WSO_lat_arr = np.loadtxt(data_dir+WSO_lat_file)
WSO_lon_arr = np.loadtxt(data_dir+WSO_lon_file)

# Importing GONG grid
GONG_lat_sin_arr = np.linspace(-1+1/180, 1-1/180, 180)
GONG_lat_arr = np.rad2deg(np.arcsin(GONG_lat_sin_arr))
GONG_lon_arr = np.linspace(0.5, 359.5, 360)

# Constructing WSO and GONG grid
WSO_lonn, WSO_latt = np.meshgrid(WSO_lon_arr, WSO_lat_arr)
GONG_lonn, GONG_latt = np.meshgrid(GONG_lon_arr, GONG_lat_arr)

# Defining interpolation function
def InterpFunc(WSO_lon_arr, WSO_lat_arr, Br, GONG_lon_arr, GONG_lat_arr):
    f = interp2d(WSO_lon_arr, WSO_lat_arr, Br, kind='cubic')
    Br_interp = f(GONG_lon_arr, GONG_lat_arr)
    return Br_interp

# Interpolating WSO field to GONG grid
for i_cr in range(cr_beg, cr_end+1):
    
    # Importing WSO field
    Br_file = 'cr' + str(i_cr) + '.dat'
    Br = np.loadtxt(data_dir+Br_file)
    
    # Interpolating
    Br_interp = InterpFunc(WSO_lon_arr, WSO_lat_arr, Br, GONG_lon_arr, GONG_lat_arr)
    
    # Saving interpolated field
    if save_or_not == 1:
        save_file = 'CR' + str(i_cr) + '_interp.csv'
        df_Br_interp = pd.DataFrame(Br_interp)
        df_Br_interp.to_csv(save_dir+'data/'+save_file, header=False, index=False)
        np.savetxt(save_dir+'data/'+'lat_interp.txt', GONG_lat_arr)
        np.savetxt(save_dir+'data/'+'lon_interp.txt', GONG_lon_arr)
    
    # Plotting figures
    plt.figure(figsize=(8,8))
    
    plt.subplot(2,1,1)
    abs_max_WSO = np.max(Br)
    plt.pcolormesh(WSO_lonn, WSO_latt, Br, cmap='RdBu', vmin=-abs_max_WSO, vmax=abs_max_WSO)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.xlim(0, 360)
    plt.ylim(-90, 90)
    plt.ylabel('Carr Lat.')
    
    plt.subplot(2,1,2)
    abs_max_GONG = np.max(Br_interp)
    plt.pcolormesh(GONG_lonn, GONG_latt, Br_interp, cmap='RdBu', vmin=-abs_max_GONG, vmax=abs_max_GONG)
    plt.colorbar()
    plt.gca().set_aspect('equal')
    plt.xlim(0, 360)
    plt.ylim(-90, 90)
    plt.xlabel('Carr Lon.')
    plt.ylabel('Carr Lat.')
    
    plt.suptitle('Photospheric field on WSO and GONG grid')
    
    if save_or_not == 1: 
        save_file = 'CR' + str(i_cr) + '_interp.png'
        plt.savefig(save_dir+'figures/'+save_file)
    
    # plt.show()
    plt.close()
    
db