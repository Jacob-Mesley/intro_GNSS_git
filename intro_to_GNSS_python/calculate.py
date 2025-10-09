# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant sequencing functions

# imports
import numpy as np
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


def az_el_range(user_ECEF, sat_ECEF, only_when_visible=True):
    # number of rows to calculate for
    rows = user_ECEF.shape[0]
    results = np.zeros((rows, 3))
    # calculating the angles and range
    for i in range(rows):
        # get data for this loop
        user = user_ECEF[i, :]  #.flatten()
        sat = sat_ECEF[i, :]  #.flatten()
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
        if only_when_visible and el < 0:
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
    # gnd_pos_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
    gnd_pos_ECEF_array = np.array(gnd_pos_ECEF)
    for _ in range(initial_range_array.shape[0] - 1):
        gnd_pos_ECEF_array = np.vstack((gnd_pos_ECEF_array, gnd_pos_ECEF))
    for j in range(5):
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
        new_az_el_range = az_el_range(gnd_pos_ECEF_array, ECEF_at_Tr, only_do_when_visible)
        R_new = new_az_el_range[:, 2]
        # test for convergence (take the difference, and use the max value)
        diff = np.nanmax(np.abs(R - R_new))
        if diff <= 1e-8:
            print("Expected range calculation converged after " + str(j) + " iterations")
            return R_new
        else:
            R = R_new
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
