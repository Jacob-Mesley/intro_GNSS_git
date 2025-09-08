# work for question 2 of the homework
import calculate
# import useful files
import plot
import numpy as np
import matplotlib.pyplot as plt

# PART B ------------------------------------------------------------------------------------

x_limit = 300  # bounds of the positions to calculate measurements for
# x_values to calculate measurements for
x_values = np.linspace(-x_limit, x_limit, x_limit*2+1)  # [m]
# given constant parameters
v_0 = 50  # [m/s^2]
h_0 = 100  # [m]

# initialize measurement variables
x_range = np.zeros(x_limit*2+1)
x_range_rate = np.zeros(x_limit*2+1)
x_zenith_ang = np.zeros(x_limit*2+1)

# calculate the measurements for each value of x
for i in range(len(x_values)):
    rng, rng_rate, znth_ang = calculate.simplified_overhead_pass_measurements(x_values[i], v_0, h_0)
    x_range[i] = rng
    x_range_rate[i] = rng_rate
    x_zenith_ang[i] = znth_ang

# plotting the results
plot.plot_arrays(x_values, x_range, "Position [m]", "Range [m]", "Position vs Range")
plot.plot_arrays(x_values, x_range_rate, "Position [m]", "Range Rate [m/s]", "Position vs Range Rate")
plot.plot_arrays(x_values, x_zenith_ang, "Position [m]", "Zenith [Rad]", "Position vs Zenith Angle")
plt.show()


