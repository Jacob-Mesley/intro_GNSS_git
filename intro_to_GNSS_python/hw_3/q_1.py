# code for question 1 of homework 3

# imports
import xarray as xr
import georinex as gr
import file
import plot
import matplotlib.pyplot as plt
import numpy as np


# Load the full RINEX once ------------------------------------------------------------

# print("before read")
# obs_data = gr.load("provided_data/NIST00USA_R_20252230000_01D_30S_MO.rnx")
# print("after read")
# # Keep only GPS data
# gps_data = obs_data.sel(sv=[sv for sv in obs_data.sv.values if sv.startswith("G")])
# # Save to a faster format (NetCDF or Zarr)
# gps_data.to_netcdf("gps_only.nc")

# import the saved data ----------------------------------------------------------------

# read in saved data
gps_data = xr.load_dataset("gps_only.nc")

# plot pseudorange vs time for PRN7
C1C_PRN07 = gps_data.C1C.sel(sv="G07").values
times = gps_data.time.to_index().to_pydatetime()
plot.plot_arrays(times, C1C_PRN07, "Time [Date]", "C1C Pseudorange [m]", "C1C vs Time for PRN7")

# plot signal-to-noise ratio vs time for PRN7
sig_to_noise_PRN07 = gps_data.S1C.sel(sv="G07").values
plot.scatter_arrays(times, sig_to_noise_PRN07,
                    "Time [Date]", "Signal to Noise Ratio", "Signal to Noise Ratio vs Time for PRN7")

# plot a new psudorange
C2W_PRN07 = gps_data.C2W.sel(sv="G07").values
plot.plot_arrays(times, C2W_PRN07, "Time [Date]", "C2W Pseudorange [m]", "C2W vs Time for PRN7")
plot.plot_arrays(times, C1C_PRN07 - C2W_PRN07, "Time [Date]", "(C1C-C2W) [m]", "(C1C-C2W) vs Time for PRN7")
plt.show()
