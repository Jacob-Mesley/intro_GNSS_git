# work for question 4 of the homework
import numpy as np
import matplotlib.pyplot as plt
import plot
import convert
import sequencing as seq
import calculate
import generate

# PART A -------------------------------------------------------

# generate the PRN19 CA code
PRN19_CA_code_binary = seq.create_CA_code(1023, [3, 6])
# convert the binary CA code to 1 and -1s
PRN19_CA_code_1s = convert.binary_array_to_1s(PRN19_CA_code_binary)
# perform the auto correlation on itself and plot
Rxy_A, lag_A = seq.cyc_corr_basic(PRN19_CA_code_1s, PRN19_CA_code_1s)
plot.plot_arrays_as_step(lag_A, Rxy_A, "Lag", "Correlation", "PRN19 vs PRN19 Auto-Corr")
# plt.show()

# PART B -----------------------------------------------------------

# shift the previously calculated PRN19 CA code by 200 bit delay
PRN19_CA_code_1s_shifted = seq.shift_and_wrap_array(PRN19_CA_code_1s, 200)
# perform autocorr and plot
Rxy_B, lag_B = seq.cyc_corr_basic(PRN19_CA_code_1s_shifted, PRN19_CA_code_1s)
plot.plot_arrays_as_step(lag_B, Rxy_B, "Lag", "Correlation", "PRN19 vs PRN19 Delayed (By 200 Bits) Auto-Corr")
# plt.show()

# PART C -----------------------------------------------------------

PRN25_CA_code_binary = seq.create_CA_code(1023, [5, 7])
PRN25_CA_code_1s = convert.binary_array_to_1s(PRN25_CA_code_binary)
Rxy_C, lag_C = seq.cyc_corr_basic(PRN19_CA_code_1s, PRN25_CA_code_1s)
plot.plot_arrays_as_step(lag_C, Rxy_C, "Lag", "Correlation", "PRN19 vs PRN25 Auto-Corr")
# plt.show()

# PART D -----------------------------------------------------------

PRN5_CA_code_binary = seq.create_CA_code(1023, [1, 9])
PRN5_CA_code_1s = convert.binary_array_to_1s(PRN5_CA_code_binary)
Rxy_D, lag_D = seq.cyc_corr_basic(PRN19_CA_code_1s, PRN5_CA_code_1s)
plot.plot_arrays_as_step(lag_D, Rxy_D, "Lag", "Correlation", "PRN19 vs PRN5 Auto-Corr")
# plt.show()

# PART E --------------------------------------------------------------------

# create the new PRNs and sum them
PRN19_delay_350 = seq.shift_and_wrap_array(PRN19_CA_code_1s, 350)
PRN25_delay_905 = seq.shift_and_wrap_array(PRN25_CA_code_1s, 905)
PRN5_delay_75 = seq.shift_and_wrap_array(PRN5_CA_code_1s, 75)
PRN_addition = PRN19_delay_350 + PRN25_delay_905 + PRN5_delay_75

# perform autocorr and plot
Rxy_E, lag_E = seq.cyc_corr_basic(PRN_addition, PRN19_CA_code_1s)
plot.plot_arrays_as_step(lag_E, Rxy_E, "Lag", "Correlation", "PRN19 vs PRN Addition (19, 25, 5) Auto-Corr")
# plt.show()

# PART F ----------------------------------------------------------------------

noise_vector = generate.noise_vector(1023, 4)
plotted_matrix = np.vstack((PRN19_delay_350, PRN25_delay_905, PRN5_delay_75, noise_vector))
plot.subplot_arrays_as_step(plotted_matrix, "Bits", "Bit Value", title="X_1, X_2, X_3, and Noise Plots, Respectively", y_limits=(-15, 15))
# plt.show()

# PART G -----------------------------------------------------------------------

PRN_with_noise = PRN_addition + noise_vector
Rxy_G, lag_G = seq.cyc_corr_basic(PRN_with_noise, PRN19_CA_code_1s)
plot.plot_arrays_as_step(lag_G, Rxy_G, "Lag", "Correlation", "PRN19 vs PRN Addition (19, 25, 5, Random Noise) Auto-Corr")
plt.show()
