# work for questions 1 through 6

# imports
import calculate
import xarray as xr
import plot
import matplotlib.pyplot as plt
import numpy as np
import convert
import file

# Variables to tweak
PRN_num = 14
show_plots_1_to_6 = False

# ----------------------------------------------
# PART 1 - Starting Point, No Error Corrections
# ----------------------------------------------

# read in data from the RNX file (read the saved gps file that was made in hw3)
rnx_data = xr.load_dataset("data/gps_only.nc")
C1C_PRN = rnx_data.C1C.sel(sv="G"+str(PRN_num)).values
times = rnx_data.time.to_index().to_pydatetime()
week_num, obs_time_of_week_PRN = convert.datetime_to_gpsweek_tow(times)
obs_wn_tow_PRN = np.hstack((np.array([week_num]).T, np.array([obs_time_of_week_PRN]).T))

# read in the broadcast ephem data
ephem_data, _ = file.read_clean_GPSbroadcast("data/brdc2230.25n")
pvt_data_PRN = convert.eph2pvt(ephem_data, obs_wn_tow_PRN, PRN_num)
pvt_pos_PRN = pvt_data_PRN[1]
pvt_week_num = pvt_data_PRN

# compute the range
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
NIST_ECEF_array = np.array(NIST_ECEF)
for _ in range(pvt_pos_PRN.shape[0] - 1):
    NIST_ECEF_array = np.vstack((NIST_ECEF_array, NIST_ECEF))
az_el_range_PRN = calculate.az_el_range(NIST_ECEF_array, pvt_pos_PRN, True)
range_PRN = az_el_range_PRN[:, 2]

# compute the expected range
expected_range = calculate.expected_range(range_PRN, NIST_ECEF, obs_time_of_week_PRN, week_num, ephem_data, PRN_num)

# calculate and plot the results
dPR0 = C1C_PRN - expected_range
time_hrs = obs_time_of_week_PRN/3600
# apply a mask to get rid of all values that are nan (aka, not visible, so no observation)
mask = ~np.isnan(dPR0)
dPR0 = dPR0[mask]
C1C_PRN = C1C_PRN[mask]
expected_range = expected_range[mask]
time_hrs = time_hrs[mask]
# plot the results
if show_plots_1_to_6:
    plot.plot_arrays(time_hrs, dPR0, "Time [Hours]", "(C1C-R) [m]",
                     "Pseudo-Range Minus Expected Range (dPR0) vs Time for PRN"+str(PRN_num))
# print out first and last values of dPR0
print("-------- PART 1, With No Corrections --------")
print(f"First Value in dPR0 -> {dPR0[0]:.3f} [m]")
print(f"Last Value in dPR0 -> {dPR0[-1]:.3f} [m]")


# ----------------------------------------------
# PART 2 - Clock Error
# ----------------------------------------------

# extract and plot the clock corrections
clock_corr = pvt_data_PRN[3][mask]
dPR1 = C1C_PRN - (expected_range - clock_corr)
y_data_matrix = np.hstack((np.array([clock_corr]).T, np.array([dPR1]).T))
if show_plots_1_to_6:
    plot.subplot_arrays(time_hrs, y_data_matrix, "Time [Hours]", ["Clock Correction [m]", "dPR1 [m]"],
                     "Clock Correction vs Time for PRN"+str(PRN_num))
# print out first and last values of dPR1
print("-------- PART 2, With Clock Correction --------")
print(f"First Value in dPR1 -> {dPR1[0]:.3f} [m]")
print(f"Last Value in dPR1 -> {dPR1[-1]:.3f} [m]")


# ----------------------------------------------
# PART 3 - Relativity
# ----------------------------------------------

# extract the relativity const
relCorr = pvt_data_PRN[4]
# apply mask to relCorr (and convert to meters)
c = 299792458  # [m/s]
relCorr = relCorr[mask]*c
dPR2 = C1C_PRN - (expected_range - clock_corr - relCorr)
# plot results
y_data_matrix = np.hstack((np.array([relCorr]).T, np.array([dPR2]).T))
if show_plots_1_to_6:
    plot.subplot_arrays(time_hrs, y_data_matrix, "Time [Hours]", ["Relativity Correction [m]", "dPR2 [m]"],
                     "Relativity Correction vs Time for PRN"+str(PRN_num))
# print first and last two values
print("-------- PART 3, With Relativity Correction --------")
print(f"First Value in dPR2 -> {dPR2[0]:.3f} [m]")
print(f"Last Value in dPR2 -> {dPR2[-1]:.3f} [m]")


# ----------------------------------------------
# PART 4 - Troposphere
# ----------------------------------------------

# calculate the tropospheric correction
zenith_delay = 2  # [m] (this was given/assumed)
el_rad = az_el_range_PRN[:, 1]*(np.pi/180)
tropoCorr = calculate.tropo_model(zenith_delay, el_rad)
tropoCorr = tropoCorr[mask]
# calculate the correction
dPR3 = C1C_PRN - (expected_range - clock_corr - relCorr + tropoCorr)
# plot the values
y_data_matrix = np.hstack((np.array([tropoCorr]).T, np.array([dPR3]).T))
if show_plots_1_to_6:
    plot.subplot_arrays(time_hrs, y_data_matrix, "Time [Hours]", ["Tropospheric Delay [m]", "dPR3 [m]"],
                     "Tropospheric Corrections vs Time for PRN"+str(PRN_num))
# print first and last two values
print("-------- PART 4, With Tropospheric Correction --------")
print(f"First Value in dPR3 -> {dPR3[0]:.3f} [m]")
print(f"Last Value in dPR3 -> {dPR3[-1]:.3f} [m]")


# ----------------------------------------------
# PART 5 - Ionosphere
# ----------------------------------------------

# calculate the ionospheric corrections
# common frequencies
f1 = 1575.42 * 10**6  # [Hz]
f2 = 1227.6 * 10**6  # [Hz]
f5 = 1176.45 * 10**6  # [Hz]
# getting the C2L pseudo-range
C2L_PRN = rnx_data.C2L.sel(sv="G"+str(PRN_num)).values
C2L_PRN = C2L_PRN[mask]
# calculate the ionospheric correction
PRIF, iono = calculate.iono_corr(C1C_PRN, f1, C2L_PRN, f2)
dPR4 = PRIF - (expected_range - clock_corr - relCorr + tropoCorr)
# plot the results
y_data_matrix = np.hstack((np.array([iono]).T, np.array([dPR4]).T))
if show_plots_1_to_6:
    plot.subplot_arrays(time_hrs, y_data_matrix, "Time [Hours]", ["Ionospheric Delay [m]", "dPR4 [m]"],
                     "Ionospheric Corrections vs Time for PRN"+str(PRN_num))
# print first and last two values
print("-------- PART 5, With Ionospheric Correction --------")
print(f"First Value in dPR4 -> {dPR4[1]:.3f} [m]")
print(f"Last Value in dPR4 -> {dPR4[-1]:.3f} [m]")


# ----------------------------------------------
# PART 6 - Comparison
# ----------------------------------------------

# collapse all data and plot
y_data_matrix = np.hstack((np.array([dPR1]).T, np.array([dPR2]).T, np.array([dPR3]).T, np.array([dPR4]).T))
if show_plots_1_to_6:
    plot.plot_multiple_arrays(time_hrs, y_data_matrix, "Time [Hours]", "Value [m]",
                     "GPS Corrections vs Time for PRN"+str(PRN_num), ["dPRN1, Clock Corr [m]", "dPR2, Relativity [m],",
                                                                      "dPRN3, Troposphere [m]", "dPR4, Ionosphere [m]"])
# plt.show()


# ----------------------------------------------
# PART 7 - Multipath
# ----------------------------------------------

# get all pseudo-ranges
C2W_PRN = rnx_data.C2W.sel(sv="G"+str(PRN_num)).values
C2W_PRN = C2W_PRN[mask]
C2L_PRN = rnx_data.C2L.sel(sv="G"+str(PRN_num)).values
C2L_PRN = C2L_PRN[mask]
C5Q_PRN = rnx_data.C5Q.sel(sv="G"+str(PRN_num)).values
C5Q_PRN = C5Q_PRN[mask]

# get all carrier phases
L1C_PRN = rnx_data.L1C.sel(sv="G"+str(PRN_num)).values
L1C_PRN = L1C_PRN[mask]
L2W_PRN = rnx_data.L2W.sel(sv="G"+str(PRN_num)).values
L2W_PRN = L2W_PRN[mask]
L2L_PRN = rnx_data.L2L.sel(sv="G"+str(PRN_num)).values
L2L_PRN = L2L_PRN[mask]
L5Q_PRN = rnx_data.L5Q.sel(sv="G"+str(PRN_num)).values
L5Q_PRN = L5Q_PRN[mask]

# calculate all data for the multipath (MP) and code-minus-carrier (CMC)
MP_C1C, CMC_C1C = calculate.mpath(C1C_PRN, L1C_PRN, f1, L2W_PRN, f2)
MP_C2W, CMC_C2W = calculate.mpath(C2W_PRN, L2W_PRN, f2, L1C_PRN, f1)
MP_C2L, CMC_C2L = calculate.mpath(C2L_PRN, L2L_PRN, f2, L1C_PRN, f1)
MP_C5Q, CMC_C5Q = calculate.mpath(C5Q_PRN, L5Q_PRN, f5, L1C_PRN, f1)

# plot CMC results
y_data_matrix = np.hstack((np.array([CMC_C1C]).T, np.array([CMC_C2W]).T, np.array([CMC_C2L]).T, np.array([CMC_C5Q]).T))
plot.plot_multiple_arrays(time_hrs, y_data_matrix, "Time [Hours]", "CMC [m]",
                 "Code Minus Carrier (CMC) for Multiple Pseudo-Ranges vs Time, PRN"+str(PRN_num), ["C1C", "C2W", "C2L", "C5Q"])

# get all signal-to-noise measurements
S1C_PRN = rnx_data.S1C.sel(sv="G"+str(PRN_num)).values
S1C_PRN = S1C_PRN[mask]
S2W_PRN = rnx_data.S2W.sel(sv="G"+str(PRN_num)).values
S2W_PRN = S2W_PRN[mask]
S2L_PRN = rnx_data.S2L.sel(sv="G"+str(PRN_num)).values
S2L_PRN = S2L_PRN[mask]
S5Q_PRN = rnx_data.S5Q.sel(sv="G"+str(PRN_num)).values
S5Q_PRN = S5Q_PRN[mask]

# plot signal-to-noise results
y_data_matrix = np.hstack((np.array([S1C_PRN]).T, np.array([S2W_PRN]).T, np.array([S2L_PRN]).T, np.array([S5Q_PRN]).T))
plot.plot_multiple_arrays(time_hrs, y_data_matrix, "Time [Hours]", "Signal-to-Noise Ratio [dB-Hz]",
                 "Signal-to-Noise vs Time, PRN"+str(PRN_num), ["S1C", "S2W", "S2L", "S5Q"])

# plot the multipath results
y_data_matrix = np.hstack((np.array([MP_C1C]).T, np.array([MP_C2W]).T, np.array([MP_C2L]).T, np.array([MP_C5Q]).T))
plot.plot_multiple_arrays(time_hrs, y_data_matrix, "Time [Hours]", "Multipath [m]",
                 "Multipath vs Time for Multiple Pseudo-Ranges, PRN"+str(PRN_num),
                          ["C1C (Second Signal L2W)", "C2W (Second Signal L1C)", "C2L (Second Signal L1C)", "C5Q (Second Signal L1C)"])

plt.show()
