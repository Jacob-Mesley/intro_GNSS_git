# code for question 3 of the homework

# imports
import file
import sequencing as seq
import convert
import calculate
import numpy as np
import plot
import xarray as xr
import matplotlib.pyplot as plt

# Values to tweak to perform analysis
PRN_num = 7
only_do_when_visible = True

# QUESTION 3 --------------------------------------------------------------------------------------

# PART A --------------------------------------

# read in the OBS file
gps_data = xr.load_dataset("gps_only.nc")
times = gps_data.time.to_index().to_pydatetime()
week_num, obs_time_of_week_PRN = convert.datetime_to_gpsweek_tow(times)
obs_wn_tow_PRN = np.hstack((np.array([week_num]).T, np.array([obs_time_of_week_PRN]).T))

# calculate the ephem pos
ephem_data, _ = file.read_clean_GPSbroadcast("provided_data/brdc2230.25n")
pvt_data_PRN = convert.eph2pvt(ephem_data, obs_wn_tow_PRN, PRN_num)
pvt_pos_PRN = pvt_data_PRN[1]
pvt_week_num = pvt_data_PRN

# calculate the az, el and range
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
NIST_ECEF_array = np.array(NIST_ECEF)
for _ in range(pvt_pos_PRN.shape[0] - 1):
    NIST_ECEF_array = np.vstack((NIST_ECEF_array, NIST_ECEF))
az_el_range_PRN = calculate.az_el_range(NIST_ECEF_array, pvt_pos_PRN, only_do_when_visible)
# plot results
# plot.subplot_arrays(obs_time_of_week_PRN/3600, az_el_range_PRN, "Time [Hours]", ["Az [Deg]", "El [Deg]", "Range [M]"],
#                     "Measurements at NIST for PRN"+str(PRN_num))


# PART B ----------------------------------------

# compute Time of Transmission (Tt) ------------

# constants
c = 299792458  # [m/s]
omega = 7.2921151467e-5  # [rad/sec] for earth
clock_corr = pvt_data_PRN[3]  # [m]
# initialize the loop
R = az_el_range_PRN[:, 2]
Tr = obs_time_of_week_PRN
for j in range(5):
    # compute the time of transmission
    Tt = Tr - (R/c)
    # compute the new sat pos at new Tt
    wn_and_Tr = np.hstack((np.array([week_num]).T, np.array([Tt]).T))
    pvt_data_PRN = convert.eph2pvt(ephem_data, wn_and_Tr, PRN_num)
    ECEF_at_Tt = pvt_data_PRN[1]
    # rotate the sat position
    theta = omega*(Tr - Tt)
    rot_matrix = np.zeros((len(theta), 3, 3))
    for i in range(len(theta)):
        rot_matrix[i] = np.array([[np.cos(theta[i]), np.sin(theta[i]), 0],
                                  [-np.sin(theta[i]), np.cos(theta[i]), 0],
                                  [0, 0, 1]])
    ECEF_at_Tr = np.zeros((len(theta), 3))
    for i in range(len(theta)):
        ECEF_at_Tr[i] = rot_matrix[i] @ ECEF_at_Tt[i]
    # compute new range
    new_az_el_range = calculate.az_el_range(NIST_ECEF_array, ECEF_at_Tr, only_do_when_visible)
    R_new = new_az_el_range[:, 2]
    # test for convergence (take the difference, and use the max value)
    diff = np.nanmax(np.abs(R - R_new))
    if diff <= 1e-8:
        print("Converged after "+str(j)+" iterations")
        break
    else:
        R = R_new

# calculate the difference
expected_range = R
old_range = az_el_range_PRN[:, 2]
R_diff = old_range - expected_range
# plot all results
# array_of_arrays = np.hstack((np.array([old_range]).T, np.array([expected_range]).T))
# plot.plot_multiple_arrays(obs_time_of_week_PRN/3600, array_of_arrays, "Time [Hours]", "Range [m]",
#                           "Range Measurements at NIST for PRN"+str(PRN_num), ["Range from Ephemeris", "Expected Range"])
# plot.plot_arrays(obs_time_of_week_PRN/3600, R_diff, "Time [Hours]", "Range Difference [m]",
#                  "Ephem Range Minus Expected Range vs Time for PRN"+str(PRN_num))

# this was for testing if my data matches the provided examples to check with
# index = None
# start_comp_time = 86400
# for i in range(len(obs_time_of_week_PRN)):
#     if obs_time_of_week_PRN[i] == float(start_comp_time):
#         index = i
# print("The 4 Range Differences Starting at Time", start_comp_time, "for PRN", PRN_num, "->", R_diff[index:index+4])


# QUESTION 4 ----------------------------------------------------------------------------------------------

# PART A ---------------------------------------------------

# get data for PRN7
C1C_PRN07 = gps_data.C1C.sel(sv="G0"+str(PRN_num)).values  # CHANGE THIS IF YOU EVER CHANGE THE PRN_num VARIABLE
# plot the difference
range_diff_with_ephem = C1C_PRN07 - expected_range

array_of_arrays = np.hstack((np.array([C1C_PRN07]).T, np.array([expected_range]).T))
plot.plot_multiple_arrays(obs_time_of_week_PRN/3600, array_of_arrays, "Time [Hours]", "Range [m]",
                          "Range Measurements at NIST for PRN"+str(PRN_num), ["C1C", "Expected Range"])

plot.plot_arrays(obs_time_of_week_PRN/3600, range_diff_with_ephem, "Time [Hours]", "Range [m]",
                 "C1C Pseudorange Minus Expected Range vs Time for PRN"+str(PRN_num))
plt.show()
