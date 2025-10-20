# work for question 3 of the midterm

# imports
import sequencing as seq
import convert
import plot
import numpy as np
import matplotlib.pyplot as plt


# PART A -------------------------------------------------------

# provided codes
A_code = [1, 1, 1, -1, -1, 1, -1]
B_code = [1, -1, -1, 1, 1, -1, 1]

# do auto-correlation to both codeS
Rxy_A, lag_A = seq.cyc_corr_basic(A_code, A_code)
Rxy_B, lag_B = seq.cyc_corr_basic(B_code, B_code)

# plot results of both auto-correlations
# plot.plot_arrays(lag_A, Rxy_A, "Lag", "Correlation", "Auto-Corr for A")
# plot.plot_arrays(lag_B, Rxy_B, "Lag", "Correlation", "Auto-Corr for B")
y_data_matrix = np.hstack((np.array([Rxy_A]).T, np.array([Rxy_B]).T))
plot.plot_multiple_arrays(lag_A, y_data_matrix, "Lag", "Correlation", "Auto-Corr", ["Code A", "Code B"])
plt.show()
