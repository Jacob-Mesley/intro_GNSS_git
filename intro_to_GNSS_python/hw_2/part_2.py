# part 2 of the homework

# imports
import numpy as np
import convert
import calculate


# given/ found user locations
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
aero_building_backyard_LLA = [40.010279, -105.243965, 1655]  # lat long from Google, alt = altitude of Boulder, CO


# PART A ----------------------------------------------------------------------------------------------------------

# convert NIST_ECEF to LLA
NIST_LLA = convert.ecef2lla(NIST_ECEF)
# define EQUA_LLA
EQUA_LLA = [0, NIST_LLA[1], 0]  # at the equator, same lon at NIST, and assume at sea level (no altitude)
# convert EQUA_LLA to ECEF
EQUA_ECEF = convert.lla2ecef(EQUA_LLA)
# convert aero_building_backyard_LLA to ECEF
aero_building_backyard_ECEF = convert.lla2ecef(aero_building_backyard_LLA)

# print out table with results
locations = [
    {"Name": "NIST",
     "LLA": NIST_LLA,
     "ECEF": NIST_ECEF},
    {"Name": "Equator (at NIST lon)",
     "LLA": EQUA_LLA,
     "ECEF": EQUA_ECEF},
    {"Name": "Aero Building Backyard",
     "LLA": aero_building_backyard_LLA,
     "ECEF": aero_building_backyard_ECEF},
]
# print formatted table
print("\nResults:")
print("-" * 110)
print(f"{'Location':<25} {'Latitude (deg)':>15} {'Longitude (deg)':>15} {'Altitude (m)':>15} "
      f"{'X (m)':>15} {'Y (m)':>15} {'Z (m)':>15}")
print("-" * 110)
for loc in locations:
    lat, lon, alt = loc["LLA"]
    x, y, z = loc["ECEF"]
    print(f"{loc['Name']:<25} {lat:15.6f} {lon:15.6f} {alt:15.0f} {x:15.0f} {y:15.0f} {z:15.0f}")
print("-" * 110)


# PART B -------------------------------------------------------------------------------------------------------------

# checking function for NIST
NIST_tmatrix = calculate.ECEF2ENU_transformation_matrix(NIST_LLA[0], NIST_LLA[1])
NIST_ECEF_vector = np.array(NIST_ECEF)
NIST_ENU_vector = NIST_tmatrix @ NIST_ECEF_vector
NIST_ENU_r_norm = NIST_ENU_vector / np.linalg.norm(NIST_ENU_vector)
NIST_ecef_r_norm = NIST_ECEF_vector / np.linalg.norm(NIST_ECEF_vector)
print("---------------------------------------------------------------------------------")
print("NIST ECEF Unit Vector -> ", NIST_ecef_r_norm)
print("NIST ENU Unit Vector -> ", NIST_ENU_r_norm)

# checking function for equator
EQUA_tmatrix = calculate.ECEF2ENU_transformation_matrix(EQUA_LLA[0], EQUA_LLA[1])
EQUA_ECEF_vector = np.array(EQUA_ECEF)
EQUA_ENU_vector = EQUA_tmatrix @ EQUA_ECEF_vector
EQUA_ENU_r_norm = EQUA_ENU_vector / np.linalg.norm(EQUA_ENU_vector)
EQUA_ecef_r_norm = EQUA_ECEF_vector / np.linalg.norm(EQUA_ECEF_vector)
print("---------------------------------------------------------------------------------")
print("EQUA ECEF Unit Vector -> ", EQUA_ecef_r_norm)
print("EQUA ENU Unit Vector -> ", EQUA_ENU_r_norm)
