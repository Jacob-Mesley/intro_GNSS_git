# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant sequencing functions

# imports
import numpy as np

import calculate
import convert
import matplotlib.pyplot as plt


# for getting simplified measurements of an aircraft overhead pass
def simplified_overhead_pass_measurements(x_a, v_a, h):
    """
    Calculating the range, range rate, and zenith angle of an overhead pass at a snapshot in time.
    This assumes the aircraft moves in a straight line tangent/parallel to the observer.

    :param x_a: sat horizontal position, wrt the observer
    :param v_a: sat velocity, wrt the observer
    :param h: sat altitude
    :return:
    """
    ac_range = np.sqrt((x_a**2 + h**2))
    ac_zenith_ang = np.arctan((x_a/h))
    ac_range_rate = (x_a * v_a) / ac_range
    return ac_range, ac_range_rate, ac_zenith_ang


def ECEF2ENU_transformation_matrix(ref_lat_deg, ref_lon_deg):
    """
    Returns the rotation matrix that transforms ECEF coordinates to ENU coordinates
    at a given reference latitude and longitude.

    Input:
        ref_lat_deg : reference latitude in degrees
        ref_lon_deg : reference longitude in degrees
    Output:
        C_ECEF2ENU : 3x3 rotation matrix
    """
    # convert to radians
    lat = np.radians(ref_lat_deg)
    lon = np.radians(ref_lon_deg)
    # define rotation matrix
    C_ECEF2ENU = np.array([
        [-np.sin(lon),              np.cos(lon),              0],
        [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon),  np.cos(lat)],
        [ np.cos(lat)*np.cos(lon),  np.cos(lat)*np.sin(lon),  np.sin(lat)]
    ])
    return C_ECEF2ENU


def az_el_range(user_ECEF, sat_ECEF, only_when_visible=True, min_el_deg=0.0):
    # number of rows to calculate for
    rows = user_ECEF.shape[0]
    results = np.zeros((rows, 3))
    # calculating the angles and range
    for i in range(rows):
        # get data for this loop
        user = user_ECEF[i, :]
        sat = sat_ECEF[i, :]
        # convert user ECEF to LLA
        lat, lon, _ = convert.ecef2lla(user)
        # rotation matrix ECEF -> ENU
        C = ECEF2ENU_transformation_matrix(lat, lon)
        # relative vector and range (satellite w.r.t. user, in ECEF)
        rho_ecef = sat - user
        rng = np.linalg.norm(rho_ecef)
        # transform into ENU to calculate az and el
        rho_enu = C @ rho_ecef
        east, north, up = rho_enu
        # azimuth (0° = North, increasing clockwise)
        az = np.degrees(np.arctan2(east, north))
        if az < 0:
            az += 360.0
        # elevation
        el = np.degrees(np.arcsin(up / rng))
        # append the calculations to the results matrix
        # if the el (as measured in the ENU frame) is ever < 0 (meaning it's now below the horizon), store Nones
        if only_when_visible and el < min_el_deg:
            results[i] = np.array((np.nan, np.nan, np.nan))
        # if is visible (or function is told not to care)
        else:
            results[i] = [az, el, rng]
    # return the results matrix when finished calculating all rows
    return results


def expected_range(initial_range_array, gnd_pos_ECEF, initial_obs_tow, week_num,
                   read_clean_BD_ephem_data, PRN_num, only_do_when_visible=True):
    # constants
    c = 299792458  # [m/s]
    omega = 7.2921151467e-5  # [rad/sec] for earth
    # initialize the loop
    R = initial_range_array
    Tr = initial_obs_tow
    # convert the ground pos into a large array to work inside the loop
    # gnd_pos_ECEF_array = np.array(gnd_pos_ECEF)
    if initial_range_array.shape[0] - 1 > 0:
        for _ in range(initial_range_array.shape[0] - 1):
            gnd_pos_ECEF_array = np.vstack((gnd_pos_ECEF_array, gnd_pos_ECEF))
    else:
        gnd_pos_ECEF_array = np.array([np.array(gnd_pos_ECEF)])

    number_of_itters = 10
    for j in range(number_of_itters):
        # compute the time of transmission
        Tt = Tr - (R / c)
        # compute the new sat pos at new Tt
        wn_and_Tr = np.hstack((np.array([week_num]).T, np.array([Tt]).T))
        pvt_data_PRN = convert.eph2pvt(read_clean_BD_ephem_data, wn_and_Tr, PRN_num)
        ECEF_at_Tt = pvt_data_PRN[1]
        # rotate the sat position
        theta = omega * (Tr - Tt)
        rot_matrix = np.zeros((len(theta), 3, 3))
        for i in range(len(theta)):
            rot_matrix[i] = np.array([[np.cos(theta[i]), np.sin(theta[i]), 0],
                                      [-np.sin(theta[i]), np.cos(theta[i]), 0],
                                      [0, 0, 1]])
        ECEF_at_Tr = np.zeros((len(theta), 3))
        for i in range(len(theta)):
            ECEF_at_Tr[i] = rot_matrix[i] @ ECEF_at_Tt[i]

        # compute new range
        new_az_el_range = az_el_range(gnd_pos_ECEF_array, ECEF_at_Tr, only_when_visible=only_do_when_visible)
        R_new = new_az_el_range[:, 2]
        # test for convergence (take the difference, and use the max value)
        diff = np.nanmax(np.abs(R - R_new))
        if diff <= 1e-8:
            # print("Expected range calculation converged after " + str(j) + " iterations")
            return R_new
        else:
            R = R_new
    print("expected range exceeded "+str(number_of_itters)+" iterations without converging! quit loop")
    return


def tropo_model(zenith_delay, elevation):
    # to protect against lower elevations blowing up the simplified model
    min_elev = np.deg2rad(5.0)
    elevation = np.maximum(elevation, min_elev)
    # predefine the size of tropoCorr
    tropoCorr = np.zeros(len(elevation))
    # calculate a tropoCorr for each elevation
    for i in range(len(elevation)):
        tropoCorr[i] = zenith_delay / np.sin(elevation[i])
    # return entire tropoCorr vector
    return tropoCorr


def iono_corr(range_1, freq_1, range_2, freq_2):
    # calcuate the iono
    iono = (freq_2**2 / (freq_1**2 - freq_2**2)) * (range_2 - range_1)
    # calculate the pif
    pif = (freq_1**2 / (freq_1**2 - freq_2**2)) * range_1 - (freq_2**2 / (freq_1**2 - freq_2**2)) * range_2
    # return results
    return pif, iono


def mpath(p_range_1, carrier_freq_1, freq_1, carrier_freq_2, freq_2):

    # calculate multipath
    c = 299792458  # [m/s]
    multipath_1 = np.zeros(len(p_range_1))
    freq_frac_1 = (freq_1 ** 2 + freq_2 ** 2) / (freq_1 ** 2 - freq_2 ** 2)
    wavelength_1 = c / freq_1
    freq_frac_2 = (2 * freq_2 ** 2) / (freq_1 ** 2 - freq_2 ** 2)
    wavelength_2 = c / freq_2
    for i in range(len(p_range_1)):
        multipath_1[i] = p_range_1[i] - (freq_frac_1*carrier_freq_1[i]*wavelength_1) + (freq_frac_2*carrier_freq_2[i]*wavelength_2)

    # calculate the code minus carrier
    code_minus_carrier = np.zeros(len(p_range_1))
    for i in range(len(p_range_1)):
        code_minus_carrier[i] = p_range_1[i] - (carrier_freq_1[i]*wavelength_1)

    # return results
    return multipath_1, code_minus_carrier


def dist_between_3d(p1, p2):
    return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)


def range_corrs(PRN_num_int, user_ECEF, obs_file_data, ephem_file_data, index_in_obs_file, add_offset_sec=None):
    # format the passed PRN info
    if PRN_num_int < 10:
        PRN_num_str = "0"+str(PRN_num_int)
    else:
        PRN_num_str = str(PRN_num_int)

    # decompose the observation file data -------------
    pseudo_range_1 = obs_file_data.C1C.sel(sv="G" + PRN_num_str).values
    pseudo_range_1 = pseudo_range_1[index_in_obs_file]
    pseudo_range_2 = obs_file_data.C2W.sel(sv="G" + PRN_num_str).values
    pseudo_range_2 = pseudo_range_2[index_in_obs_file]
    times = obs_file_data.time.to_index().to_pydatetime()
    week_num, tow = convert.datetime_to_gpsweek_tow(times)
    # obs_wn_tow_PRN = np.hstack((np.array([week_num]).T, np.array([obs_time_of_week_PRN]).T))

    # consider offset (if there even is one)
    if add_offset_sec:
        for i in range(len(tow)):
            tow[i] = tow[i] + add_offset_sec
    obs_wn_tow_PRN = np.array([[week_num[index_in_obs_file]], [tow[index_in_obs_file]]]).T

    # get the iono-free pseudo-range ------------------
    f1 = 1575.42 * 10 ** 6  # [Hz]
    f2 = 1227.6 * 10 ** 6  # [Hz]
    f5 = 1176.45 * 10 ** 6  # [Hz]
    iono_free_range, _ = iono_corr(pseudo_range_1, f1, pseudo_range_2, f2)
    if np.isnan(iono_free_range):
        return None


    # get the PVT data --------------------------------
    pvt_data = convert.eph2pvt(ephem_file_data, obs_wn_tow_PRN, PRN_num_int)
    sat_pos_ECEF = pvt_data[1]

    # get the az, el, and range -----------------------
    user_ECEF_array = np.array([np.array(user_ECEF)])
    az_el_range_final = az_el_range(user_ECEF_array, sat_pos_ECEF, only_when_visible=False)
    initial_range = az_el_range_final[:, 2]
    # mask = ~np.isnan(initial_range)

    # get the expected range ---------------------------
    expected_range_final = expected_range(initial_range, user_ECEF, [tow[index_in_obs_file]],
            [week_num[index_in_obs_file]], ephem_file_data, PRN_num_int, only_do_when_visible=False)
    if expected_range_final is None:
        return None
    else:
        expected_range_final = expected_range_final[0]

    # clock correction ---------------------------------
    clock_corr_final = pvt_data[3]

    # relativity correction ----------------------------
    relCorr = pvt_data[4]
    c = 299792458  # [m/s]
    rel_corr_final = relCorr * c

    # tropospheric corrections -------------------------
    zenith_delay = 2  # [m] (this was given/assumed)
    el_rad = az_el_range_final[:, 1] * (np.pi / 180)
    tropo_corr_final = tropo_model(zenith_delay, el_rad)

    # return all results
    return iono_free_range, expected_range_final, az_el_range_final[0], clock_corr_final[0], rel_corr_final[0], tropo_corr_final[0]


def direction_to_unit_vector(from_pos, to_pos):
    direction = np.array(to_pos) - np.array(from_pos)
    norm = np.linalg.norm(direction)
    if norm == 0:
        return np.zeros_like(direction)
    return direction / norm


def gps_pos_solution(pos_guess_ECEF, obs_data, ephem_data, index_to_calc_solution, do_output=False):
    # initialize the loop
    maximum_iterations = 20
    iterations = 0
    min_sats_visible = 5
    index = index_to_calc_solution
    current_guess_ECEF = np.array(pos_guess_ECEF)

    # read in observation file
    times = obs_data.time.to_index().to_pydatetime()
    week_num, obs_tow = convert.datetime_to_gpsweek_tow(times)
    obs_wn_tow_PRN = np.column_stack(([week_num[index]], [obs_tow[index]]))

    # do a check to see which satellites should be visible (and therefore should be used)
    # also store pvt data
    PRN_nums = []
    ECEF_sat_positions = []
    # exclude PRN 32, because of it's weird offset issue with the specific data
    for i in range(1, 32):
        pvt_data = convert.eph2pvt(ephem_data, obs_wn_tow_PRN, i)
        ECEF_pos = np.array(pvt_data[1])
        temp = np.array([current_guess_ECEF])
        is_visible = calculate.az_el_range(temp, ECEF_pos, only_when_visible=True, min_el_deg=5.0)
        if np.isnan(is_visible[0, 0]):
            pass
        else:
            PRN_nums.append(i)
            ECEF_pos = ECEF_pos.flatten()
            ECEF_sat_positions.append(ECEF_pos)
    number_of_PRNs = len(PRN_nums)

    # enforce minimum sat visible
    if number_of_PRNs < min_sats_visible:
        return
    # loop through iterations
    for j in range(maximum_iterations):
        # calculate all the data
        range_corrs_all = np.zeros((number_of_PRNs, len(obs_data.C1C.sel(sv="G01").values), 7))
        for PRN_num_in_array in range(number_of_PRNs):
            PRN_num = PRN_nums[PRN_num_in_array]
            corrs = range_corrs(PRN_num, current_guess_ECEF, obs_data, ephem_data, index)
            if corrs is None:
                return
            PIF, expected_range, az_el_range, clkCorr, relCorr, tropoCorr = corrs
            # write data to the data structure to be saved after the loop
            temp_array = np.column_stack(
                (PIF, expected_range, az_el_range[1], az_el_range[0], clkCorr, relCorr, tropoCorr))
            range_corrs_all[PRN_num_in_array] = temp_array

        # calculate pre-fit residuals based on calculated data
        _, rows, cols = range_corrs_all.shape
        prefit_resids = np.zeros((number_of_PRNs, 1))
        for PRN_num_in_array in range(number_of_PRNs):
            data = range_corrs_all[PRN_num_in_array]
            PIF = data[0, 0]
            expected_range = data[0, 1]
            clkCorr = data[0, 4]
            relCorr = data[0, 5]
            tropoCorr = data[0, 6]
            model = expected_range - clkCorr - relCorr + tropoCorr
            prefit_resids[PRN_num_in_array] = PIF - model

        # construct sensitivity (aka geometry) matrix
        g_matrix = np.zeros((number_of_PRNs, 4))
        for PRN_num_in_array in range(number_of_PRNs):
            ECEF_pos = ECEF_sat_positions[PRN_num_in_array]
            los_vector = direction_to_unit_vector(ECEF_pos, current_guess_ECEF)
            g_matrix[PRN_num_in_array] = np.concatenate((los_vector, [1]))
        # calculate the deviation and new position guess
        x_hat = np.linalg.inv(g_matrix.T @ g_matrix) @ g_matrix.T @ prefit_resids
        new_guess_ECEF = (current_guess_ECEF + x_hat[0:3].T).flatten()
        # for counting number of iterations
        iterations = iterations+1
        state_vector_corr_mag = np.linalg.norm(x_hat[0:3])
        if do_output:
            decimals = 2
            print("----------- Iteration: "+str(j+1)+" -----------")
            print("State Deviation (x_hat) [m] -> "+str(np.round(x_hat.T.flatten(), decimals)))
            print("State Deviation (x_hat) Magnitude [m] -> "+str(np.round(state_vector_corr_mag, 4)))
            print("Updated State Estimate ECEF [m] -> "+str(np.round(new_guess_ECEF, decimals)))

        if state_vector_corr_mag < 0.01:  # only stop iteration once deviation is below 1 cm (ie, 0.01 meters)
            # set the final solution
            final_sol = new_guess_ECEF
            # calc final post fit resids
            post_fit_resids = prefit_resids - g_matrix @ x_hat
            # calc error in the ENU frame (this is unique to the problem i was doing)
            truth_pos_ECEF = np.array([-1288398.567, -4721696.932, 4078625.350])  # Location of NIST in Boulder, CO [m]
            ref_lla = convert.ecef2lla(truth_pos_ECEF)
            ENU_t_matrix = ECEF2ENU_transformation_matrix(ref_lla[0], ref_lla[1])
            final_diff = truth_pos_ECEF - final_sol
            final_error_ENU = ENU_t_matrix @ final_diff
            # calculate the HDOP and VDOP
            h = np.linalg.inv(g_matrix.T @ g_matrix)
            h_pos = h[:3, :3]
            h_tilde = ENU_t_matrix @ h_pos @ ENU_t_matrix.T
            HDOP = np.sqrt(h_tilde[0, 0] + h_tilde[1, 1])
            VDOP = np.sqrt(h_tilde[2, 2])
            # grab clock bias
            clock_bias = x_hat[3]
            # return final results
            return final_sol, clock_bias, prefit_resids.flatten(), post_fit_resids.flatten(), final_error_ENU, HDOP, VDOP, PRN_nums
        else:
            current_guess_ECEF = new_guess_ECEF
    print("After -> ("+str(maximum_iterations)+") Iterations, the gps solution did not converge")




