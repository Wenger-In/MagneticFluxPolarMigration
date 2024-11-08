# MagneticFluxPolarMigration

This repository is designed for the research of the poleward migration of solar magnetic fields. This research is based on the observations of Wilcox Solar Observatory, WSO, which provide the longest record of synoptic maps, and the PFSS model, which provides the connection between the photosphere and the source surface. We are interested in the 'butterfly diagram' on the source surface, and whether it differs from the 'butterfly diagram' on the photosphere. These may reflect the poleward migration of magnetic fields on different levels in the solar corona. 

This repository includes the following steps: 

1. Download the WSO synoptic maps from [WSO website](http://wso.stanford.edu/synopticl.html), where the *filled_synop* takes the grid on equal steps of longitude (from 360 to 0) and sine latitude (from +14.5/15 to -14.5/15).

   --> This step can be done by 'batch_download.py', and the original data is saved as 'CR****.txt'.

2. Neaten the original WSO data for easier reading, and save the corresponding grids. 

   --> This step can be done by 'neaten_WSO_photospheric_field.m', and the neatened data is saved as 'cr****.dat', together with the grids as 'lat_arr.dat' and 'lon_arr.dat'

3. Interpolate the WSO maps to GONG grids, and save the interpolated fields.
   
   --> This step can be done by 'interp2gong_grid.py', and the interpolated data is saved as 'CR****_interp.dat'. 

4. Construct Sunpy.map based on the interpolated data and corresponding meta, use PFSS model to obtain the magnetic field on source surface.

   --> This step can be done by 'pfss4WSO_field.py', and the source surface fields are saved as 'CR****_ss.dat'

5. Average the magnetic fields on the photosphere and the source surface among all longitudes, respectively. Then, plot the averaged fields on the coordinates of 'Time-Latitude', and find the features of the 'butterfly-diagrams'.

   --> This step can be done by 'plot_butterfly_diagram.py', and the averaged fields are saved as 'avg_Br_phs.dat' and 'avg_Br_ss.dat'.

**Note**: Here for the input of PFSS model, meta is required, however, WSO data does not have the meta data intrinsically. So, we will check whether we can obtain the same source surface fields with any meta. If the answer is YES, we can keep the same meta (from GONG) for the PFSS model.
