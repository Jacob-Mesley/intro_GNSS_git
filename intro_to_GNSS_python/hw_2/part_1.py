# part 1 of the homework

# imports
import file
import plot
import matplotlib.pyplot as plt
import sequencing as seq
import convert
import numpy as np

# choose a PRN number
PRN_num = 10


# PART A ---------------------------------------------------------------------------------------

# read in GPS data from the sp3 file
all_GPS_data_sp3 = file.read_sp3("provided_data/IGS0OPSFIN_20252230000_01D_15M_ORB.SP3")
# extract the PRN10 sat info
GPS_PRN_data = seq.extract_rows_by_index(all_GPS_data_sp3, PRN_num-1, 31)
# extract the pos and time from the sat
GPS_PRN_pos = GPS_PRN_data[:, 3:6]*1000
week_num, time_of_week = GPS_PRN_data[:, 0], GPS_PRN_data[:, 1]


# PART B ----------------------------------------------------------------------------------------

# read the yuma data
all_GPS_data_yuma, _ = file.read_GPSyuma("provided_data/YUMA223.alm")
# get the almanac week number, and combine with the time of week from the sp3 (assuming PRN10)
PRN_data_yuma = all_GPS_data_yuma[PRN_num-1]
week_num_array_yuma = np.full(GPS_PRN_data.shape[0], PRN_data_yuma[18])
yuma_time_matrix = np.vstack((week_num_array_yuma, GPS_PRN_data[:, 1]))
yuma_time_matrix = yuma_time_matrix.T
# convert the almanac data to ECEF pos
health, GPS_PRN_alm_pos, clock_corr = convert.alm2pos(all_GPS_data_yuma, yuma_time_matrix, PRN_num)
# plot results
plot.subplot_multiple_arrays(time_of_week/3600, np.array([GPS_PRN_pos, GPS_PRN_alm_pos]), "Time [Hrs]", ["X [m]",
"Y [m]", "Z [m]"], ["SP3", "Almanac"], "GPS-PRN" + str(PRN_num) + " Position vs Time")


# PART C -------------------------------------------------------------------------------------

# plotting difference between the almanac and sp3
pos_diff_PRN = GPS_PRN_alm_pos - GPS_PRN_pos
plot.subplot_arrays(time_of_week/3600, pos_diff_PRN, "Time [Hrs]", ["X [m]", "Y [m]", "Z [m]"], "GPS-PRN" +
                    str(PRN_num) + " Almanac Minus SP3 Position vs Time")
plt.show()
