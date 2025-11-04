# work for question 1 of the homework

# imports
import numpy as np
import convert
import calculate
import xarray as xr
import file
import sequencing as seq
import plot
import matplotlib.pyplot as plt

PRN_nums = [1, 3, 4, 6, 9, 16, 26, 28, 31]
number_of_PRNs = len(PRN_nums)

print("\n\n-------------------------- QUESTION 1 --------------------------")
print("\n---------- PART A ----------")
# find the error of the initial guess
truth_pos_ECEF = np.array([-1288398.567, -4721696.932, 4078625.350])  # Location of NIST in Boulder, CO [m]
initial_guess_LLA = [40, -105, 1600]  # [Lat (deg), Lon (deg), Alt (m)]
initial_guess_ECEF = np.array(convert.lla2ecef(initial_guess_LLA))
initial_error_ECEF = truth_pos_ECEF - initial_guess_ECEF
print("Initial Guess Error (Truth - Guess) in ECEF: \n"+str(np.round(initial_error_ECEF, 2))+" [m]")

# find the straight line distance between the guess and the truth
distance = calculate.dist_between_3d(truth_pos_ECEF, initial_guess_ECEF)
print("Magnitude of Distance Between Truth and Guess: \n["+str(np.round(distance, 2))+"] [m]")


print("\n---------- PART B ----------")
# read in data from the RNX (aka observation) file to know correct times
obs_data = xr.load_dataset("data/gps_only.nc")
# read in the broadcast ephem data
ephem_data, _ = file.read_clean_GPSbroadcast("data/brdc2230.25n")
# for disabling the loop
print_table = True
if print_table:
    range_corrs_all = np.zeros((number_of_PRNs, len(obs_data.C1C.sel(sv="G01").values), 7))
    # start table of data
    print(f"{'PRN':<5} {'PIF (m)':>10} {'R1 (m)':>10} {'EL (deg)':>10} {'AZ (deg)':>10} "
        f"{'BSV (m)':>10} {'RELAT (m)':>12} {'TROP (m)':>10}")
    print("-" * 87)
    # calculate all the data and print in the table (these are the PRNs saved to the data file)
    for PRN_num_in_array in range(number_of_PRNs):
        PRN_num = PRN_nums[PRN_num_in_array]
        # PRN 32 has a 2hr (7200 sec) offset wrt it's observations in the obs file (for whatever reason)
        if PRN_num == 32:
            offset = 7200
        else:
            offset = 0
        # get data
        index = 0
        PIF, expected_range, az_el_range, clkCorr, relCorr, tropoCorr = calculate.range_corrs(PRN_num,
                                            initial_guess_ECEF, obs_data, ephem_data, index, offset)
        # format all values
        PIF_form = np.round(PIF, 2)
        exp_rng_form = np.round(expected_range, 2)
        el_form = np.round(az_el_range[1], 1)
        az_form = np.round(az_el_range[0], 1)
        clckCorr_form = np.round(clkCorr, 2)
        relCorr_form = np.round(relCorr, 2)
        tropoCorr_form = np.round(tropoCorr, 2)
        # print formatted row
        print(
            f"{PRN_num:<5} {PIF_form:>10.1f} {exp_rng_form:>10.1f} {el_form:>10.1f} {az_form:>10.1f} "
            f"{clckCorr_form:>10.3f} {relCorr_form:>12.3f} {tropoCorr_form:>10.3f}")
        # write data to the data structure to be saved after the loop
        temp_array = np.column_stack((PIF, expected_range, az_el_range[1], az_el_range[0], clkCorr, relCorr, tropoCorr))
        # if save data, store results in the array
        range_corrs_all[PRN_num_in_array] = temp_array


print("\n\n-------------------------- QUESTION 2 --------------------------")
print("\n---------- PART A ----------")

# calculate pre-fit residuals based on data
_, rows, cols = range_corrs_all.shape
prefit_resids = np.zeros((number_of_PRNs, 1))
for PRN_num_in_array in range(number_of_PRNs):
    data = range_corrs_all[PRN_num_in_array]
    for i in range(rows):
        PIF = data[i, 0]
        expected_range = data[i, 1]
        clkCorr = data[i, 4]
        relCorr = data[i, 5]
        tropoCorr = data[i, 6]
        model = expected_range - clkCorr - relCorr + tropoCorr
        prefit_resids[PRN_num_in_array] = PIF - model
# print first pre-fits
print("First Pre-fit Residuals")
for PRN_num_in_array in range(number_of_PRNs):
    PRN_num = PRN_nums[PRN_num_in_array]
    first_prefit = np.round(prefit_resids[PRN_num_in_array], 2)
    # print result
    print("PRN "+str(PRN_num)+" -> "+str(first_prefit)+" [m]")


print("\n---------- PART B ----------")

# grabs correct time vector to pass into pvt data func
pseudo_range_temp = obs_data.C2W.sel(sv="G01").values
times = obs_data.time.to_index().to_pydatetime()
week_num, obs_tow = convert.datetime_to_gpsweek_tow(times)

# construct sensitivity (aka geometry) matrix
index = 0
obs_wn_tow_PRN = np.column_stack(([week_num[index]], [obs_tow[index]]))
g_matrix = np.zeros((number_of_PRNs, 4))
for PRN_num_in_array in range(number_of_PRNs):
    PRN_num = PRN_nums[PRN_num_in_array]
    pvt_data = convert.eph2pvt(ephem_data, obs_wn_tow_PRN, PRN_num)
    ECEF_pos = pvt_data[1][index]
    los_vector = calculate.direction_to_unit_vector(ECEF_pos, initial_guess_ECEF)
    g_matrix[PRN_num_in_array] = np.concatenate((los_vector, [1]))
print("Sensitivity Matrix (aka Geometry Matrix): ")
print(np.round(g_matrix, 5))

print("\n---------- PART C, D, E, F AND G ----------")
final_sol, clock_bias, prefits, postfits, error_ENU, HDOP, VDOP, num_sats = calculate.gps_pos_solution(
    initial_guess_ECEF, obs_data, ephem_data, 0, do_output=True)
final_diff = truth_pos_ECEF - final_sol
first_sol = [-1288397.1219354, -4721689.96271919, 4078621.49494259]
first_diff = truth_pos_ECEF - first_sol
print("\n----- Final Results -----")
print("First Iteration Error [m] -> "+str(np.round(first_diff, 2)))

print("\n----- Final Iteration Errors: -----")
print("Clock Bias [m] -> "+str(np.round(clock_bias, 2)))
print("Final Iteration Error [m] -> "+str(np.round(final_diff, 2)))
print("Pre-Fits [m] -> "+str(np.round(prefits, 2)))
print("Post-Fits [m] -> "+str(np.round(postfits, 2)))
print("Error in the ENU Frame [m] -> "+str(np.round(error_ENU, 2)))
print("HDOP [m] -> "+str(np.round(HDOP, 3)))
print("VDOP [m] -> "+str(np.round(VDOP, 3)))
print("Number of GPS Sats Used to find Solution -> "+str(len(num_sats)))

print("\n\n-------------------------- QUESTION 3 --------------------------")

# get useful sizing/plotting vars
data_len = len(obs_data.C1C.sel(sv="G01").values)
test_len = data_len
times = obs_data.time.to_index().to_pydatetime()
_, obs_tow = convert.datetime_to_gpsweek_tow(times)
time_hrs = obs_tow/3600


# GENERATE DATA ---------------------------------
generate_and_save_data = False
# initialize the data storage vars
if generate_and_save_data:
    results = []
    # generate all the data
    for i in range(test_len):
        # call function to find data
        result = calculate.gps_pos_solution(initial_guess_ECEF, obs_data, ephem_data, i)
        # if data can successfully be found
        if result:
            results.append(result)
        # if not
        else:
            results.append([np.nan])
        print("percent complete -> ", np.round(i/test_len * 100, 2), "%")
    # save data
    np.save("results.npy", np.array(results, dtype=object))


# DATA ANALYSIS -----------------------------------
# initialize arrays to plot
postfits = np.zeros(test_len)
prefits = np.zeros(test_len)
error_ENU_clk_bias = np.zeros((test_len, 4))
num_sats_HDOP_VDOP = np.zeros((test_len, 3))
results = np.load("results.npy", allow_pickle=True)
for i in range(test_len):
    result = results[i]
    if len(result) == 1:  # ie, an array with just a nan
        postfits[i] = np.nan
        prefits[i] = np.nan
        error_ENU_clk_bias[i] = [np.nan, np.nan, np.nan, np.nan]
        num_sats_HDOP_VDOP[i] = [np.nan, np.nan, np.nan]
    else:  # if thhe array is longer, assume it is full of valid numbers
        final_sol_temp, clock_bias_temp, prefits_temp, postfits_temp, error_ENU_temp, HDOP, VDOP, PRN_nums = result
        postfits[i] = np.mean(postfits_temp)
        prefits[i] = np.mean(prefits_temp)
        error_ENU_clk_bias[i] = np.hstack((error_ENU_temp, clock_bias_temp))
        num_sats_HDOP_VDOP[i] = np.hstack((len(PRN_nums), HDOP, VDOP))
# plot results
time = time_hrs[0:test_len]
plot.subplot_arrays(time, error_ENU_clk_bias, "Time [hrs]", ["East [m]", "North [m]", "Up [m]", "Clock Bias [m]"],
                    "ENU Error and Clock Bias vs Time")  # PART B
plot.subplot_arrays(time, num_sats_HDOP_VDOP, "Time [hrs]", ["Number of Visible Satellites", "HDOP", "VDOP"],
                    "Number of GPS Visible, HDOP and VDOP vs Time")  # PART C
fits = np.column_stack((prefits, postfits))
plot.subplot_arrays(time, fits, "Time [hrs]", ["Pre-Fits", "Post-Fits [m]"], "Pre and Post-Fit Average vs Time")  # PART D
plt.show()
