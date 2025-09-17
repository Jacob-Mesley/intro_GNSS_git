# part 3 of the homework

# imports
import numpy as np
import calculate
import file
import sequencing as seq
import plot
import matplotlib.pyplot as plt

# PART A ------------------------------------------------------------------------------

# sanity check for the az_el_range function
R_earth = 6378137.0  # [m]
test_user_pos_ECEF = np.array([[R_earth, 0, 0]])  # on the equator
test_GPS_pos_ECEF = np.array([[R_earth + 2e6, 0, 0]])  # on the equator, directly overhead of the user
test_az_el_range = calculate.az_el_range(test_user_pos_ECEF, test_GPS_pos_ECEF)
print("TEST POS: GPS directly overhead at altitude 2e6 meters, the \n"
      "[az (deg), el (deg), range (m)] is ->", test_az_el_range.flatten())


# PART B -------------------------------------------------------------------------------

# should be the same as the PRN chosen in part 1 of the homework
PRN_num = 10
# read in GPS data from the sp3 file
all_GPS_data_sp3 = file.read_sp3("provided_data/IGS0OPSFIN_20252230000_01D_15M_ORB.SP3")
# extract the PRN10 sat info
GPS_PRN_data = seq.extract_rows_by_index(all_GPS_data_sp3, PRN_num-1, 31)
# extract the pos and time from the sat
GPS_PRN_pos = GPS_PRN_data[:, 3:6]*1000
week_num, time_of_week = GPS_PRN_data[:, 0], GPS_PRN_data[:, 1]

# creating the NIST_ECEF array to pair with the GPS_PRN_pos array
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
NIST_ECEF_array = np.array(NIST_ECEF)
for _ in range(GPS_PRN_pos.shape[0] - 1):
    NIST_ECEF_array = np.vstack((NIST_ECEF_array, NIST_ECEF))

# calculating the az, el, and range at NIST for PRN10, and plotting results
measurements = calculate.az_el_range(NIST_ECEF_array, GPS_PRN_pos, only_when_visible=False)
plot.subplot_arrays(time_of_week/3600, measurements, "Time [Hrs]", ["Az [deg]", "El [deg]", "Range [m]"], "GPS-PRN" +
                    str(PRN_num) + " Measurements as seen by NIST (Assuming Always Visible)")


# PART C -----------------------------------------------------------------------------------

# repeat of the earlier plot, but only when the sat is visible
measurements = calculate.az_el_range(NIST_ECEF_array, GPS_PRN_pos, only_when_visible=True)
plot.subplot_arrays(time_of_week/3600, measurements, "Time [Hrs]", ["Az [deg]", "El [deg]", "Range [m]"], "GPS-PRN" +
                    str(PRN_num) + " Measurements as seen by NIST (Adjusted for Visibility)")
# plt.show()


# PART D -----------------------------------------------------------------------------------

PRN_num = 2
# read in GPS data from the sp3 file
all_GPS_data_sp3 = file.read_sp3("provided_data/IGS0OPSFIN_20252230000_01D_15M_ORB.SP3")
# extract the PRN10 sat info
GPS_PRN_data = seq.extract_rows_by_index(all_GPS_data_sp3, PRN_num-1, 31)
# extract the pos and time from the sat
GPS_PRN_pos = GPS_PRN_data[:, 3:6]*1000
week_num, time_of_week = GPS_PRN_data[:, 0], GPS_PRN_data[:, 1]

# creating the NIST_ECEF array to pair with the GPS_PRN_pos array
NIST_ECEF = [-1288398.567, -4721696.932, 4078625.350]  # provided by prof (in [m])
NIST_ECEF_array = np.array(NIST_ECEF)
for _ in range(GPS_PRN_pos.shape[0] - 1):
    NIST_ECEF_array = np.vstack((NIST_ECEF_array, NIST_ECEF))

measurements = calculate.az_el_range(NIST_ECEF_array, GPS_PRN_pos, only_when_visible=True)
plot.subplot_arrays(time_of_week/3600, measurements, "Time [Hrs]", ["Az [deg]", "El [deg]", "Range [m]"], "GPS-PRN" +
                    str(PRN_num) + " Measurements as seen by NIST (Adjusted for Visibility)")
plt.show()




