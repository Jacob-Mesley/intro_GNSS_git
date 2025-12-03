# work for homework 7, q2 - q5

# imports
from scipy.io import loadmat
import matplotlib.pyplot as plt
import plot
import calculate
import file
import numpy as np
import convert
import sequencing as seq

# QUESTION 2 ------------------------------------
# load in the provided mat file
mat = loadmat("HW7data.mat")
signal = mat['signal'].flatten()  # 1x4e6 sized matrix containing complex values
# creating PRN3 C/A code
PRN_tap = [4, 8]  # PRN-03
# PRN_tap = [3, 8]  # PRN-31
# PRN_tap = [6, 8]  # PRN-26
# PRN_tap = [9, 10]  # PRN-16
PRN_num_string = "PRN03"
PRN_CA_code = seq.create_CA_code(1023, PRN_tap)
PRN_CA_code_1s = 1 - 2 * PRN_CA_code

# creating time vector
samp_freq = 8e6  # [Hz] 8 MHz
duration = 2e-3  # [sec] 1 millisecond
num_samples = int(samp_freq * duration)
time = np.arange(num_samples) / samp_freq
# corresponding PRN CA code it to the time vector
chips_per_ms = 1023
chip_rate = chips_per_ms / 1e-3
chip_duration = 1 / chip_rate  # duration of one chip
chip_indices = ((time // chip_duration) % chips_per_ms).astype(int)
# new PRN code
PRN_code_formatted = PRN_CA_code_1s[chip_indices]

# testing the function
DOP_freq = [0, 3500, 1000]
signal_shift = [9, 294, 294]
for i in range(len(DOP_freq)):
    S = calculate.complex_correlator(signal, DOP_freq[i], signal_shift[i], PRN_code_formatted, time)
    print("DOP Freq: ("+str(DOP_freq[i])+"), Sample Shift: ("+str(signal_shift[i])+"), S = "+str(np.round(S, 1)))


# QUESTION 3 -----------------------------------------------------------

# do the search for PRN-XX
generate_data = False
file_name = "complex_corr_2ms_"+PRN_num_string+".npy"
delay_span = np.arange(0, 1023, 0.5)  # 8185
DOP_freq_span = np.arange(0, 10500, 500)
if generate_data:
    data = np.zeros((len(delay_span), len(DOP_freq_span)))
    for i in range(len(delay_span)):
        for j in range(len(DOP_freq_span)):
            S = calculate.complex_correlator(signal, DOP_freq_span[j], delay_span[i], PRN_code_formatted, time)
            data[i, j] = abs(S)
        print("Progress "+str(np.round(i/len(delay_span)*100, 2))+" %")
    np.save(file_name, data)
else:
    data = np.load(file_name)

# normalize data to peak of 1
max_val = np.max(data)
if max_val != 0:
    data_norm = data / max_val
else:
    data_norm = data

# convert code shift to actual chip shift
delay_span = delay_span * (chip_rate/samp_freq)
# find index of maximum value
max_idx = np.unravel_index(np.argmax(data), data.shape)
max_code_shift = delay_span[max_idx[0]]
max_doppler = DOP_freq_span[max_idx[1]]
print(f"Max occurs at: code_shift = {max_code_shift}, doppler = {max_doppler} Hz")

# extract the entire bin for the strongest peak
strongest_codeshift_bin = data[max_idx[0], :]   # row -> all Doppler bins at this code shift
strongest_doppler_bin = data[:, max_idx[1]]   # column -> all code-shift bins at this Doppler

# make plots
plot.plot_arrays(DOP_freq_span, strongest_codeshift_bin, "Doppler Range [Hz]", "Magnitude of S", PRN_num_string+" - Peak Bin: Sample Shift")
plot.plot_arrays(delay_span, strongest_doppler_bin, "Code Delay [Chips]", "Magnitude of S", PRN_num_string+" - Peak Bin: Doppler Shift")
plot.complex_corr(DOP_freq_span, delay_span, data_norm, title=PRN_num_string+" - Complex Correlation of Doppler and Sample Shift")
plt.show()
