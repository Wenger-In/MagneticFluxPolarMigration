import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Referring to raw data, from CR1642 to CR2287
cr_beg = 1642
cr_end = 2287
phs_dir = 'E:/Research/Work/flux_transport_on_ss/Br_interp/data/'
ss_dir = 'E:/Research/Work/flux_transport_on_ss/Br_ss/data/'
save_dir = 'E:/Research/Work/flux_transport_on_ss/butterfly/'
save_or_not = 0

# Importing WSO grid
lat_file = 'lat_interp.txt'
lon_file = 'lon_interp.txt'
lat_arr = np.loadtxt(phs_dir+lat_file)
lon_arr = np.loadtxt(phs_dir+lon_file)

# Averaging among all longtitudes
Br_phs_ts = np.empty((180, 0))
Br_ss_ts = np.empty((180, 0))
for i_cr in range(cr_beg, cr_end+1):
    
    # Importing photospheric fields
    phs_file = 'CR' + str(i_cr) + '_interp.csv'
    df_Br_phs = pd.read_csv(phs_dir+phs_file, header=None)
    Br_phs  = df_Br_phs.values
    Br_phs_avg = np.nanmean(Br_phs, axis=1)
    
    # Importing source surface fields
    ss_file = 'CR' + str(i_cr) + '_ss.csv'
    df_Br_ss = pd.read_csv(ss_dir+ss_file, header=None)
    Br_ss  = df_Br_ss.values
    Br_ss_avg = np.nanmean(Br_ss, axis=1)
    
    # Combining to obtain the time series of the averaged fields
    Br_phs_ts = np.column_stack((Br_phs_ts, Br_phs_avg))
    Br_ss_ts = np.column_stack((Br_ss_ts, Br_ss_avg))

plt.figure(figsize=(12,6))

cr_arr = np.arange(cr_beg, cr_end+1)
crr, latt = np.meshgrid(cr_arr, lat_arr)

plt.subplot(2,1,1)
abs_max_phs = np.max(Br_phs_ts)
plt.pcolormesh(crr, latt, Br_phs_ts, cmap='RdBu', vmin=-abs_max_phs, vmax=abs_max_phs)
cbr = plt.colorbar()
cbr.set_label('Br (uT)')
plt.ylim(-90, 90)
plt.ylabel('Carr Lat.')
plt.title('Br on the photosphere')

plt.subplot(2,1,2)
abs_max_ss = np.max(Br_ss_ts)
plt.pcolormesh(crr, latt, Br_ss_ts, cmap='RdBu', vmin=-abs_max_ss, vmax=abs_max_ss)
cbr = plt.colorbar()
cbr.set_label('Br (uT)')
plt.ylim(-90, 90)
plt.xlabel('Carr Rotation')
plt.ylabel('Carr Lat.')
plt.title('Br on the source surface')

plt.show()

db