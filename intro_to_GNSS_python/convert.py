# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant conversion operations

# imports
import numpy as np
from datetime import datetime, timezone


def binary_array_to_hex(bit_array):
    # make sure it's a flat array of ints
    bit_array = np.array(bit_array, dtype=int).flatten()
    # join into a string of bits
    bit_string = ''.join(str(b) for b in bit_array)
    # convert binary string -> int -> hex
    hex_str = hex(int(bit_string, 2))[2:].upper()
    return hex_str


# converts 0 -> 1 and 1 -> -1 in an array of binary
def binary_array_to_1s(binary_array):
    converted_binary = []
    for i in range(len(binary_array)):
        if binary_array[i] == 0:
            converted_binary.append(1)
        elif binary_array[i] == 1:
            converted_binary.append(-1)
    return converted_binary


# converting the date and time to GPS week, and seconds into the week (provided by prof. P. Axelrad)
def cal2gps(dt):
    """
    Convert datetime to GPS week and time-of-week (seconds).
    GPS epoch: 1980-01-06 00:00:00 UTC
    """
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)
    delta = dt - gps_epoch
    wn = delta.days // 7
    tow = (delta.days % 7) * 86400 + delta.seconds + delta.microseconds * 1e-6
    return wn, tow


import numpy as np


# Helper: solve Kepler's Equation (mean -> eccentric anomaly)
def mean2eccentric(M, ecc, tol=1e-12, max_iter=50):
    """
    Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E.
    """
    # Normalize M into [-pi, pi] for stability
    M = np.mod(M + np.pi, 2 * np.pi) - np.pi
    # Initial guess
    if ecc < 0.8:
        E = M
    else:
        E = np.pi
    for _ in range(max_iter):
        f = E - ecc * np.sin(E) - M
        fprime = 1 - ecc * np.cos(E)
        dE = -f / fprime
        E = E + dE
        if abs(dE) < tol:
            break
    return E


def alm2pos(almanac, t_input, prn):
    """
    Calculate GPS satellite position and clock correction from a YUMA almanac.

    Args:
        almanac (np.ndarray): Almanac matrix (n x 25) from read_GPSyuma.
        t_input (np.ndarray): Times [WN, TOW] shape (k x 2).
        prn (int): PRN number of satellite to compute.

    Returns:
        health (np.ndarray): Health status (k x 1).
        satPos (np.ndarray): Satellite ECEF position (k x 3) in meters.
        satClkCorr (np.ndarray): Satellite clock correction (k x 1) in meters.
    """
    # WGS-84 constants
    muE = 3.986005e14  # m^3/s^2
    wE = 7.2921151467e-5  # rad/s
    c = 2.99792458e8  # m/s

    sz = t_input.shape[0]
    satPos = np.full((sz, 3), np.nan)
    health = np.full((sz, 1), np.nan)
    satClkCorr = np.full((sz, 1), np.nan)

    # Find almanac row for this PRN
    mask = almanac[:, 0] == prn
    if not np.any(mask):
        return health, satPos, satClkCorr  # no data

    sat_alm = almanac[mask][0]

    for tt in range(sz):
        Toe = sat_alm[16]  # col17 in MATLAB (0-based index 16)
        gps_wk = sat_alm[18]  # col19 in MATLAB (0-based index 18)

        dt = (t_input[tt, 0] - gps_wk) * 604800 + (t_input[tt, 1] - Toe)

        a = sat_alm[4] ** 2  # col5 sqrt_a
        ecc = sat_alm[3]  # col4
        n0 = np.sqrt(muE / a ** 3)
        n = n0 + sat_alm[2]  # col3 delta_n
        M = sat_alm[1] + n * dt  # col2 M0
        inc = sat_alm[6]  # col7
        perigee = sat_alm[7]  # col8

        # Eccentric anomaly
        E = mean2eccentric(M, ecc)
        cosE, sinE = np.cos(E), np.sin(E)

        # True anomaly
        nu = np.arctan2(np.sqrt(1 - ecc ** 2) * sinE, cosE - ecc)

        # Argument of latitude
        u = nu + perigee

        # Radius
        r = a * (1 - ecc * cosE)

        # Orbital plane coordinates
        xo = r * np.cos(u)
        yo = r * np.sin(u)

        # Corrected longitude of ascending node
        node = sat_alm[5] + (sat_alm[8] - wE) * dt - (wE * Toe)

        cosi, sini = np.cos(inc), np.sin(inc)
        coso, sino = np.cos(node), np.sin(node)

        # ECEF position (m)
        satPos[tt, 0] = xo * coso - yo * cosi * sino
        satPos[tt, 1] = xo * sino + yo * cosi * coso
        satPos[tt, 2] = yo * sini

        # Clock correction (meters)
        satClkCorr[tt, 0] = c * ((sat_alm[22] * dt + sat_alm[21]) * dt + sat_alm[20])

        # Health
        health[tt, 0] = sat_alm[24]

    return health, satPos, satClkCorr


def ecef2lla(ecef_xyz):
    """
    Convert ECEF coordinates to latitude, longitude, altitude (LLA).
    Input:
        ecef_xyz : ECEF position as an array [x, y, z] in meters
    Output:
        lat (deg), lon (deg), alt (m)
    """
    # for running iterations
    tol = 1e-8
    max_iter = 100
    # WGS84 constants
    R_earth = 6378137.0  # semi-major axis (meters)
    f = 1 / 298.257223563  # ellipsoid reciprocal flattening
    e_squared = 2*f - f**2  # eccentricity squared
    # decompose the ECEF inputs
    x, y, z = ecef_xyz
    # calculate longitude
    lon = np.arctan2(y, x)
    # distance from Z-axis
    p = np.sqrt(x**2 + y**2)
    # calculate the radius
    r = np.sqrt((x**2 + y**2 + z**2))
    # initialize the latitude approximation with this first guess
    lat = np.arcsin(z/r)
    # iterative latitude
    for _ in range(max_iter):
        c_plus = R_earth / np.sqrt(1 - e_squared*np.sin(lat)**2)
        lat_old = lat
        lat = np.arctan2((z + c_plus*e_squared*np.sin(lat)), p)
        # if tolerance has been met, break out of the loop
        if abs(lat_old - lat) < tol:
            break
    # use final lat to calculate the altitude
    c_plus = R_earth / np.sqrt(1 - e_squared*np.sin(lat)**2)
    alt = p / np.cos(lat) - c_plus
    # convert to degrees
    lat = np.degrees(lat)
    lon = np.degrees(lon)
    # return the result
    return [lat, lon, alt]


def lla2ecef(lla):
    """
    Convert Latitude, Longitude, Altitude (LLA) to ECEF coordinates.
    Input:
        lla : [lat (deg), lon (deg), alt (m)]
    Output:
        [x, y, z] in meters (ECEF)
    """
    # WGS84 constants
    R_earth = 6378137.0  # semi-major axis (meters)
    f = 1 / 298.257223563  # ellipsoid reciprocal flattening
    e_squared = 2 * f - f ** 2  # eccentricity squared
    # unpack input
    lat, lon, alt = lla
    # convert degrees to radians
    lat = np.radians(lat)
    lon = np.radians(lon)
    # prime vertical radius of curvature
    c_plus = R_earth / np.sqrt(1 - e_squared * np.sin(lat)**2)
    # compute coordinates
    x = (c_plus + alt) * np.cos(lat) * np.cos(lon)
    y = (c_plus + alt) * np.cos(lat) * np.sin(lon)
    z = (c_plus * (1 - e_squared) + alt) * np.sin(lat)

    return [x, y, z]


def eph2pvt(ephemeris, t_input, prn):
    """
    Compute GPS satellite position, velocity, and clock correction
    from broadcast ephemeris.

    Parameters
    ----------
    ephemeris : ndarray (n, 25)
        Ephemeris matrix from read_clean_GPSbroadcast
    t_input : ndarray (m, 2)
        GPS times to calculate at, as [WeekNumber, TimeOfWeek] per row
    prn : int
        Satellite PRN to compute for

    Returns
    -------
    health : ndarray (m,)
        Satellite health (0 = good)
    satPos : ndarray (m, 3)
        Satellite ECEF position [m]
    satVel : ndarray (m, 3)
        Satellite ECEF velocity [m/s]
    satClkCorr : ndarray (m,)
        Satellite clock correction [m]
    junk : int
        Placeholder (always 0 for now)
    tgd : ndarray (m,)
        Group delay [m]
    """
    # WGS-84 constants
    muE = 3.986005e14       # m^3/s^2
    wE  = 7.2921151467e-5   # rad/s
    c   = 2.99792458e8      # m/s

    m = t_input.shape[0]
    satPos     = np.full((m, 3), np.nan)
    satVel     = np.full((m, 3), np.nan)
    satClkCorr = np.full(m, np.nan)
    relCorr    = np.full(m, np.nan)
    tgd        = np.full(m, np.nan)
    health     = np.full(m, np.nan)

    # Get all ephemerides for this PRN
    mask = ephemeris[:, 0] == prn
    sat_eph0 = ephemeris[mask, :]
    if sat_eph0.size == 0:
        return health, satPos, satVel, satClkCorr, 0, tgd

    for tt in range(m):
        wk, tow = t_input[tt]

        # Time difference search
        dt_search = (wk - sat_eph0[:, 18]) * 604800 + (tow - sat_eph0[:, 16])
        dt_search = dt_search[dt_search >= -10]

        if dt_search.size == 0:
            continue
        dt_min_index = np.argmin(np.abs(dt_search))
        if np.abs(dt_search[dt_min_index]) > 2 * 3600:
            print("Warning: Ephemeris expired!")

        sat_eph = sat_eph0[dt_min_index, :]

        Toe     = sat_eph[16]
        gps_wk  = sat_eph[18]
        dt      = (wk - gps_wk) * 604800 + (tow - Toe)
        a       = sat_eph[4]**2
        ecc     = sat_eph[3]
        n0      = np.sqrt(muE / a**3)
        n       = n0 + sat_eph[2]
        M       = sat_eph[1] + n * dt
        inc     = sat_eph[6]
        perigee = sat_eph[7]

        # Eccentric anomaly
        E = mean2eccentric(M, ecc)
        cosE, sinE = np.cos(E), np.sin(E)
        nu = np.arctan2(np.sqrt(1 - ecc**2) * sinE, cosE - ecc)
        u = nu + perigee

        # Harmonic corrections
        Cuc, Cus, Crc, Crs, Cic, Cis = sat_eph[10:16]
        du = Cus * np.sin(2 * u) + Cuc * np.cos(2 * u)
        dr = Crs * np.sin(2 * u) + Crc * np.cos(2 * u)
        di = Cis * np.sin(2 * u) + Cic * np.cos(2 * u)

        r = a * (1 - ecc * cosE)
        u = u + du
        r = r + dr
        inc = inc + di + dt * sat_eph[9]

        cosu, sinu = np.cos(u), np.sin(u)
        cos2u, sin2u = np.cos(2 * u), np.sin(2 * u)

        xo, yo = r * cosu, r * sinu
        node = sat_eph[5] + (sat_eph[8] - wE) * dt - wE * Toe

        cosi, sini = np.cos(inc), np.sin(inc)
        coso, sino = np.cos(node), np.sin(node)

        # Position in ECEF
        satPos[tt, 0] = xo * coso - yo * cosi * sino
        satPos[tt, 1] = xo * sino + yo * cosi * coso
        satPos[tt, 2] = yo * sini

        # Velocities
        E_dot = n / (1 - ecc * cosE)
        nu_dot = E_dot * np.sqrt(1 - ecc**2) / (1 - ecc * cosE)
        i_dot = sat_eph[9] + 2 * nu_dot * (Cis * cos2u + Cic * sin2u)
        u_dot = nu_dot + 2 * nu_dot * (Cus * cos2u - Cuc * sin2u)
        r_dot = a * ecc * sinE * n / (1 - ecc * cosE) + \
                2 * nu_dot * (Crs * cos2u - Crc * sin2u)
        node_dot = sat_eph[8] - wE

        xo_dot = r_dot * cosu - yo * u_dot
        yo_dot = r_dot * sinu + xo * u_dot

        satVel[tt, 0] = (xo_dot - yo * cosi * node_dot) * coso - \
                        (xo * node_dot + yo_dot * cosi - yo * sini * i_dot) * sino
        satVel[tt, 1] = (xo_dot - yo * cosi * node_dot) * sino + \
                        (xo * node_dot + yo_dot * cosi - yo * sini * i_dot) * coso
        satVel[tt, 2] = yo_dot * sini + yo * cosi * i_dot

        # Clock correction (meters)
        Af0, Af1, Af2 = sat_eph[20:23]
        satClkCorr[tt] = c * ((Af2 * dt + Af1) * dt + Af0)

        health[tt] = sat_eph[24]
        relCorr[tt] = 0
        tgd[tt] = c * sat_eph[23]  # col 24

    junk = 0
    return health, satPos, satVel, satClkCorr, junk, tgd


# def datetime_to_gpsweek_tow(dt):
#     """Convert datetime to (gpsweek, seconds-of-week)."""
#     # GPS epoch: Jan 6, 1980
#     gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc)
#     dt = dt.replace(tzinfo=timezone.utc)
#     delta = (dt - gps_epoch).total_seconds()
#     gpsweek = int(delta // 604800)
#     tow = delta - gpsweek * 604800
#     return gpsweek, tow


def datetime_to_gpsweek_tow(date_time):
    """
    Convert array of datetimes to GPS week number and time-of-week (seconds).
    """
    GPS_EPOCH = datetime(1980, 1, 6, 0, 0, 0)
    gps_weeks = []
    tows = []
    for dt in date_time:
        delta = dt - GPS_EPOCH
        gps_week = delta.days // 7
        tow = (delta.days % 7) * 86400 + delta.seconds + dt.microsecond * 1e-6
        gps_weeks.append(gps_week)
        tows.append(tow)
    return np.array(gps_weeks), np.array(tows)

