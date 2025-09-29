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


def read_clean_GPSbroadcast(navfilename, clean=False):
    """
    Reads an IGS GPS Broadcast Ephemeris (RINEX nav file) and constructs a matrix
    of ephemeris values. Optionally cleans bad ephemerides.

    Parameters
    ----------
    navfilename : str
        Path to the IGS broadcast ephemeris file (.nav).
    clean : bool, optional
        Whether to clean incorrect entries. Default = False.

    Returns
    -------
    gps_ephem : ndarray (n, 25)
        Broadcast ephemeris values:
          col1:  prn
          col2:  M0
          col3:  delta_n
          col4:  ecc
          col5:  sqrt_a
          col6:  Loa
          col7:  incl
          col8:  perigee
          col9:  ra_rate
          col10: i_rate
          col11: Cuc
          col12: Cus
          col13: Crc
          col14: Crs
          col15: Cic
          col16: Cis
          col17: Toe
          col18: IODE
          col19: GPS_week
          col20: Toc
          col21: Af0
          col22: Af1
          col23: Af2
          col24: TGD
          col25: health

    ionoparams : ndarray (8,)
        Klobuchar ionospheric parameters [A0 A1 A2 A3 B0 B1 B2 B3].
    """
    gps_ephem = []
    ALPHA, BETA = None, None

    with open(navfilename, 'r') as f:
        # --- Parse header ---
        while True:
            line = f.readline()
            if not line:
                raise ValueError("File ended before END OF HEADER.")
            if "ION ALPHA" in line:
                parts = line.split()[:4]
                ALPHA = [float(x.replace('D', 'E')) for x in parts]
            if "ION BETA" in line:
                parts = line.split()[:4]
                BETA = [float(x.replace('D', 'E')) for x in parts]
            if "END OF HEADER" in line:
                break

        ionoparams = np.array(ALPHA + BETA if (ALPHA and BETA) else [])

        # --- Parse ephemeris blocks ---
        while True:
            line = f.readline()
            if not line:
                break

            if not line.strip():
                continue

            prn = int(line[0:2])
            Af0 = float(line[22:41].replace('D', 'E'))
            Af1 = float(line[41:60].replace('D', 'E'))
            Af2 = float(line[60:79].replace('D', 'E'))

            # line 2
            line = f.readline()
            IODE = float(line[3:22].replace('D', 'E'))
            Crs = float(line[22:41].replace('D', 'E'))
            delta_n = float(line[41:60].replace('D', 'E'))
            M0 = float(line[60:79].replace('D', 'E'))

            # line 3
            line = f.readline()
            Cuc = float(line[3:22].replace('D', 'E'))
            ecc = float(line[22:41].replace('D', 'E'))
            Cus = float(line[41:60].replace('D', 'E'))
            sqrt_a = float(line[60:79].replace('D', 'E'))

            # line 4
            line = f.readline()
            Toe = float(line[3:22].replace('D', 'E'))
            Toc = Toe
            Cic = float(line[22:41].replace('D', 'E'))
            Loa = float(line[41:60].replace('D', 'E'))
            Cis = float(line[60:79].replace('D', 'E'))

            # line 5
            line = f.readline()
            incl = float(line[3:22].replace('D', 'E'))
            Crc = float(line[22:41].replace('D', 'E'))
            perigee = float(line[41:60].replace('D', 'E'))
            ra_rate = float(line[60:79].replace('D', 'E'))

            # line 6
            line = f.readline()
            i_rate = float(line[3:22].replace('D', 'E'))
            GPS_week = int(float(line[41:60].replace('D', 'E')))
            if GPS_week < 1024:
                GPS_week += 2048

            # line 7
            line = f.readline()
            health = float(line[22:41].replace('D', 'E'))
            TGD = float(line[41:60].replace('D', 'E'))

            # line 8 (skip)
            f.readline()

            gps_ephem.append([
                prn, M0, delta_n, ecc, sqrt_a, Loa, incl, perigee,
                ra_rate, i_rate, Cuc, Cus, Crc, Crs, Cic, Cis, Toe,
                IODE, GPS_week, Toc, Af0, Af1, Af2, TGD, health
            ])

    gps_ephem = np.array(gps_ephem)

    # --- Optional cleaning ---
    if clean and gps_ephem.size > 0:
        badrows = []
        for prn in np.unique(gps_ephem[:, 0]):
            sat_ephem_idx = np.where(gps_ephem[:, 0] == prn)[0]
            sat_ephem = gps_ephem[sat_ephem_idx, :]
            newuploads = np.where(sat_ephem[:, 16] % 3600 != 0)[0]  # Toe in col 17 -> index 16
            for idx in newuploads:
                if idx < len(sat_ephem) - 1:
                    if (sat_ephem[idx + 1, 16] - sat_ephem[idx, 16]) < 240:
                        badrows.append(sat_ephem_idx[idx + 1])
        gps_ephem = np.delete(gps_ephem, badrows, axis=0)

    return gps_ephem, ionoparams

