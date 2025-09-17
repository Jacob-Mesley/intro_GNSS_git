# part 4 of the homework

# imports
import plot
import matplotlib.pyplot as plt
import numpy as np
import file
import convert
import calculate


# PART A ------------------------------------------------------

# given/ found user locations
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
aero_building_backyard_LLA = [40.010279, -105.243965, 1655]  # lat long from Google, alt = altitude of Boulder, CO
NIST_LLA = convert.ecef2lla(NIST_ECEF)
EQUA_LLA = [0, NIST_LLA[1], 0]  # at the equator, same lon at NIST, and assume at sea level (no altitude)
EQUA_ECEF = convert.lla2ecef(EQUA_LLA)
aero_building_backyard_ECEF = convert.lla2ecef(aero_building_backyard_LLA)

# read in and extract all GPA data from the sp3 file
all_GPS_data_sp3 = file.read_sp3("provided_data/IGS0OPSFIN_20252230000_01D_15M_ORB.SP3")
all_GPS_pos = all_GPS_data_sp3[:, 3:6]*1000
week_num, time_of_week = all_GPS_pos[:, 0], all_GPS_pos[:, 1]

# find lat and long wrt NIST and plot
NIST_ECEF_array = np.array(NIST_ECEF)
for _ in range(all_GPS_data_sp3.shape[0] - 1):
    NIST_ECEF_array = np.vstack((NIST_ECEF_array, NIST_ECEF))
user_GPS_lla_NIST = calculate.az_el_range(NIST_ECEF_array, all_GPS_pos)
plot.plot_az_el(user_GPS_lla_NIST[:, 0], user_GPS_lla_NIST[:, 1], title="GPS Coverage from NIST")

# find lat long wrt equator
user_GPS_lla_EQUA = np.array(EQUA_ECEF)
for _ in range(all_GPS_data_sp3.shape[0] - 1):
    user_GPS_lla_EQUA = np.vstack((user_GPS_lla_EQUA, EQUA_ECEF))
user_GPS_lla_EQUA = calculate.az_el_range(user_GPS_lla_EQUA, all_GPS_pos)
plot.plot_az_el(user_GPS_lla_EQUA[:, 0], user_GPS_lla_EQUA[:, 1], title="GPS Coverage from EQUA")

# find lat long wrt aero building backyard
are0_building_backyard_ECEF_array = np.array(aero_building_backyard_ECEF)
for _ in range(all_GPS_data_sp3.shape[0] - 1):
    are0_building_backyard_ECEF_array = np.vstack((are0_building_backyard_ECEF_array, aero_building_backyard_ECEF))
user_GPS_lla = calculate.az_el_range(are0_building_backyard_ECEF_array, all_GPS_pos)
plot.plot_az_el(user_GPS_lla[:, 0], user_GPS_lla[:, 1], title="GPS Coverage from Aero Building Backyard")

# show all plots
plt.show()
