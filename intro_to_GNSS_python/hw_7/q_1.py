# work for homework 7, q1

# imports
from scipy.io import loadmat
import matplotlib.pyplot as plt
import plot
import calculate
import file
import numpy as np
import convert

# QUESTION 1 ------------------------------------
# read in the almanac for Aug 24th, 2025 using the alamanc
alm_data, _ = file.read_GPSyuma("current_yuma.alm")
prn_nums = []
for i in range(len(alm_data)):
    prn_nums.append(alm_data[i][0])

# convert almanac to position data for the provided time
wn_tow = np.array([[333.0, (22*3600 + 39*60 + 12)]])
prn_pos_data_ECEF = []
for i in range(len(prn_nums)):
    _, pos_data_ECEF, _ = convert.alm2pos(alm_data, wn_tow, prn_nums[i])
    prn_pos_data_ECEF.append(pos_data_ECEF)

# calculate the az-el-range of visible gps wrt to the user
user_lla = [40.0, -105.244, 1620]  # (lat [deg N], long [deg W], alt [m])
user_ECEF = convert.lla2ecef(user_lla)
user_ECEF_formatted = np.array([user_ECEF])
GPS_az_el_range = []
visible_prns = []
for i in range(len(prn_nums)):
    prn_pos_data_formatted = np.array(prn_pos_data_ECEF[i])
    temp = calculate.az_el_range(user_ECEF_formatted, prn_pos_data_formatted)
    temp = temp.flatten()
    if not np.isnan(temp[0]):
        visible_prns.append(prn_nums[i])
        GPS_az_el_range.append(temp.flatten())

# generate and show skyplot
az = []
el = []
for i in range(len(visible_prns)):
    az.append(GPS_az_el_range[i][0])
    el.append(GPS_az_el_range[i][1])
plot.plot_az_el(az, el, "GPS Sky Plot - August 24th 2025, 22:39:12", visible_prns)
plt.show()
