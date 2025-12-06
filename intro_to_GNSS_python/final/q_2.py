# work for question 1 of the final

# imports
import numpy as np
import convert
import calculate
import xarray as xr
import file
import georinex as gr
import plot
import matplotlib.pyplot as plt


# reading in the provided rnx files, and converting them to a formal that works with my code -------------
# print("before read")
# obs_data = gr.load("NIST00USA_R_20252570000_SINGLE_EPOCH.rnx")
# print("after read")
# # keep only GPS data
# gps_data = obs_data.sel(sv=[sv for sv in obs_data.sv.values if sv.startswith("G")])
# # save
# gps_data.to_netcdf("obs_data_NIST.nc")

# loading in the data ------------------------------------------
obs_data = xr.load_dataset("obs_data_NIST.nc")
# read in the broadcast ephem data
ephem_data, _ = file.read_clean_GPSbroadcast("brdc2570.25n")

# Initial Guess -----------------------------------------------
initial_guess_LLA = [40, -105, 1600]  # [Lat (deg), Lon (deg), Alt (m)]
initial_guess_ECEF = np.array(convert.lla2ecef(initial_guess_LLA))
print("Initial Guess [m] -> "+str(np.round(initial_guess_ECEF, 1)))

# do solution --------------------------------------------------
final_sol, clock_bias, prefits, postfits, error_ENU, HDOP, VDOP, num_sats = calculate.gps_pos_solution(
    initial_guess_ECEF, obs_data, ephem_data, 0, do_output=True)

# print final results ------------------------------------------
print("\n----- Final Iteration Info Expanded: -----")
decimals = 1
print("Final Solution LLA [deg, deg, m] -> "+str(np.round(convert.ecef2lla(final_sol), decimals)))
print("Clock Bias [m] -> "+str(np.round(clock_bias, decimals)))
print("GPS Sats Used to find Solution -> "+str(num_sats))
print("Pre-Fits [m] -> "+str(np.round(prefits, decimals)))
print("Post-Fits [m] -> "+str(np.round(postfits, decimals)))
print("Error in the ENU Frame [m] -> "+str(np.round(error_ENU, decimals)))
print("HDOP [m] -> "+str(np.round(HDOP, 3)))
print("VDOP [m] -> "+str(np.round(VDOP, 3)))


# for validating code against NIST
# print("\n---- NIST Validation ----")
# truth_pos_NIST_ECEF = np.array([-1288398.567, -4721696.932, 4078625.350])  # Location of NIST in Boulder, CO [m]
# print("initial error -> ", truth_pos_NIST_ECEF - initial_guess_ECEF)
# print("final error -> ", truth_pos_NIST_ECEF - final_sol)
