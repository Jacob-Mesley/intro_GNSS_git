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
