# Author: Jacob Mesley
# File Created: 9/8/2025
# Last Edit: 9/8/2025
# Description: relevant file read functions

# imports
import numpy as np
from datetime import datetime, timedelta
import convert


# for reading sp3 files (converted from MATLAB code provided by prof. P. Axelrad)
def read_sp3(fname):
    """
    Read an SP3 orbit file.

    Returns:
        numpy.ndarray with columns:
        [week, TOW(s), PRN, X(km), Y(km), Z(km), clock_bias(us), constellation_id]
        constellation_id: 1=GPS, 2=GLO, 3=GAL, 4=BDS, 5=QZSS
    """
    constellation_map = {'G': 1, 'R': 2, 'E': 3, 'C': 4, 'J': 5}
    sp3 = []
    with open(fname, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            msg = line[0]
            if msg == '*':
                # Example line: *  2023  9  5  0 15  0.00000000
                temp = line[3:].split()
                year, month, day, hour, minute, sec = map(float, temp)
                dt = datetime(int(year), int(month), int(day), int(hour), int(minute), int(sec))
                wn, tow = convert.cal2gps(dt)
            elif msg == 'P':
                # Satellite line
                const_char = line[1]
                prn_num = int(line[2:4])
                temp = list(map(float, line[5:].split()))
                # temp contains [X, Y, Z, clock_bias, ...]
                X, Y, Z, clk = temp[:4]
                constellation = constellation_map.get(const_char, 0)
                sp3.append([wn, tow, prn_num, X, Y, Z, clk, constellation])
    return np.array(sp3)


def read_GPSyuma(yuma_filename, rollovers=0):
    """
    Reads a GPS YUMA almanac file and constructs a matrix of all ephemeris values.

    Args:
        yuma_filename (str): Path to YUMA almanac file.
        rollovers (int): Number of 1024-week rollovers (default=0).

    Returns:
        gps_alm (np.ndarray): Ephemeris matrix (n x 25).
        gps_alm_cell (list of dict): List of dictionaries with ephemeris values.
    """
    gps_alm_cell = []
    gps_alm = []

    try:
        with open(yuma_filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"Unable to find YUMA almanac: {yuma_filename}")

    j = 0
    i = 0
    while i < len(lines):
        tline = lines[i].strip()
        i += 1
        if not tline:
            continue
        if not tline.startswith("**"):  # header lines like "***Week ..."
            continue

        # skip junk line (***Week...)
        if i >= len(lines):
            break
        tline = lines[i].strip()
        i += 1
        if not tline:
            break

        eph = {}

        # ---- First line: PRN
        eph["prn"] = int(tline[26:].strip())

        # ---- Second line: Health
        eph["health"] = int(lines[i].strip()[26:]); i += 1

        # ---- Third line: eccentricity
        eph["ecc"] = float(lines[i].strip()[26:]); i += 1

        # ---- Fourth line: Toe
        eph["Toe"] = float(lines[i].strip()[26:]); i += 1

        # ---- Fifth line: inclination
        eph["incl"] = float(lines[i].strip()[26:]); i += 1

        # ---- Sixth line: rate of RAAN
        eph["ra_rate"] = float(lines[i].strip()[26:]); i += 1

        # ---- Seventh line: sqrt(a)
        eph["sqrt_a"] = float(lines[i].strip()[26:]); i += 1

        # ---- Eighth line: LoA
        eph["Loa"] = float(lines[i].strip()[26:]); i += 1

        # ---- Ninth line: argument of perigee
        eph["perigee"] = float(lines[i].strip()[26:]); i += 1

        # ---- Tenth line: mean anomaly
        eph["M0"] = float(lines[i].strip()[26:]); i += 1

        # ---- Eleventh line: Af0
        eph["Af0"] = float(lines[i].strip()[26:]); i += 1

        # ---- Twelfth line: Af1
        eph["Af1"] = float(lines[i].strip()[26:]); i += 1

        # ---- Thirteenth line: GPS week
        GPS_week = int(lines[i].strip()[26:]); i += 1
        eph["GPS_week"] = GPS_week + rollovers * 1024

        # Fill broadcast fields not in YUMA
        eph.update({
            "delta_n": 0, "i_rate": 0,
            "Cuc": 0, "Cus": 0, "Crc": 0, "Crs": 0, "Cic": 0, "Cis": 0,
            "IODE": 0, "Toc": 0, "Af2": 0
        })

        row = [
            eph["prn"], eph["M0"], eph["delta_n"], eph["ecc"], eph["sqrt_a"],
            eph["Loa"], eph["incl"], eph["perigee"], eph["ra_rate"], eph["i_rate"],
            eph["Cuc"], eph["Cus"], eph["Crc"], eph["Crs"], eph["Cic"], eph["Cis"],
            eph["Toe"], eph["IODE"], eph["GPS_week"], eph["Toc"], eph["Af0"], eph["Af1"],
            eph["Af2"], 0, eph["health"]
        ]

        gps_alm.append(row)
        gps_alm_cell.append(eph)

        j += 1

        # skip possible blank line
        if i < len(lines) and not lines[i].strip():
            i += 1

    gps_alm = np.array(gps_alm)
    return gps_alm, gps_alm_cell

