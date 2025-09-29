# question 2 code

# imports
import file
import convert
import plot
import sequencing as seq
import matplotlib.pyplot as plt

# pick the PRN number to perform analysis on
PRN_num = 7

# PART A AND B ------------------------------------------------------------------------------

# read in the sp3 file from homework 3, and get PRN data
sp3_data = file.read_sp3("provided_data/IGS0OPSFIN_20252230000_01D_15M_ORB.SP3")
sp3_data_PRN = seq.extract_PRN_from_sp3_data(sp3_data, PRN_num)
sp3_pos_PRN = sp3_data_PRN[:, 3:6]*1000
sp3_week_num_PRN, sp3_time_of_week_PRN = sp3_data_PRN[:, 0], sp3_data_PRN[:, 1]
sp3_wn_tow_PRN = sp3_data_PRN[:, 0:2]
# plot the pos vs time
plot.subplot_arrays(sp3_time_of_week_PRN/3600, sp3_pos_PRN,
                    "Time [Hours]", ["X [m]", "Y [m]", "Z [m]"], "SP3 Data vs Time for PRN"+str(PRN_num))
# read in and clean up the ephemeris data
ephem_data, _ = file.read_clean_GPSbroadcast("provided_data/brdc2230.25n")
pvt_data_PRN = convert.eph2pvt(ephem_data, sp3_wn_tow_PRN, PRN_num)
pvt_pos_PRN = pvt_data_PRN[1]
# plot the pos vs time
plot.subplot_arrays(sp3_time_of_week_PRN/3600, pvt_pos_PRN,
                    "Time [Hours]", ["X [m]", "Y [m]", "Z [m]"], "PVT Data vs Time for PRN"+str(PRN_num))
# plot the difference
plot.subplot_arrays(sp3_time_of_week_PRN/3600, pvt_pos_PRN-sp3_pos_PRN,
                    "Time [Hours]", ["X [m]", "Y [m]", "Z [m]"], "(PVT-SP3) vs Time for PRN"+str(PRN_num))
# plt.show()


# PART C ------------------------------------------------------------------------------------

# the clock bias is already calculated internal to the function, it uses the a0 and a1 parameters to find this
# value. it is redundant to calculate it twice, so extract this bias form the output of the function called earlier.
clock_bias_from_function = pvt_data_PRN[3]
# plot the clock corr
plot.plot_arrays(sp3_time_of_week_PRN/3600, clock_bias_from_function,
                 "Time [Hours]", "Clock Bias [m]", "Clock Bias vs Time for PRN"+str(PRN_num))
plt.show()
