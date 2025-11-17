# work for all homework 6 questions (1-8)

# imports
import numpy as np
from scipy import signal
import sequencing as seq
import matplotlib.pyplot as plt
import plot

print("|---------------- PART 1 ----------------|\n")
# creating the time array
f_sample_rate = 50e6  # sampling frequency [Hz]
duration = 10e-3  # duration [sec]
t = np.arange(0, duration, 1/f_sample_rate)
print("Sample Rate ->", f_sample_rate/1e6, "MHz, Duration ->", duration*1e3, "msec")
# nyquist frequency
f_nyquist = f_sample_rate/2
print("Nyquist Frequency ->", f_nyquist/1e6, "MHz")
# frequency resolution
f_resolution = f_sample_rate/len(t)
print("Frequency Resolution ->", f_resolution, "Hz")
# duration to produce an f_resolution of 5 Hz
desired_f_res = 5
desired_len = f_sample_rate / desired_f_res
desired_duration = desired_len / f_sample_rate
print("Duration required for a Frequency Resolution of 5 Hz ->", desired_duration, "sec")


print("\n|---------------- PART 2 ----------------|\n")
# creating a sine wave and plotting on scopes
f_sine = 40e3  # [Hz]
amplitude = 1.0  # [V]
sine_wave = amplitude * np.sin(2 * np.pi * f_sine * t)
# plot.oscope_and_spectrum_analyser(t, sine_wave, str(f_sine/1e3)+" kHz Sine Wave", oscope_xlim_microsec=(0, 200), freq_xlim_kHz=(0, 100))


print("\n|---------------- PART 3 ----------------|\n")
# creating a square wave and plotting scopes
f_square = 40e3  # [Hz]
amplitude = 1.0  # [V]
square_wave = amplitude * signal.square(2 * np.pi * f_square * t)
# plot.oscope_and_spectrum_analyser(t, square_wave, str(f_square/1e3)+" kHz Square Wave", oscope_xlim_microsec=(0, 200), freq_xlim_kHz=(0, 500))


print("\n|---------------- PART 4 ----------------|\n")
# create maximal length G1 code
g1_code = seq.create_G1_code(1023)
# create the "stream" we would observe if this code was received at the provided
# chipping_rate, and we sampled at the provided f_sample_rate
chipping_rate = 1.023*1e6  # [Hz]
chip_duration = 1 / chipping_rate
chip_indices = np.floor(t / chip_duration).astype(int) % len(g1_code)
g1_code_stream = g1_code[chip_indices]
# printing required elements from the array
element_checks = np.linspace(485, 495, 11)
print("[Time (micro-sec), Value (From Maximal Length 1023 G-Code)] for elements (485) to (495), where each "
      "\narray starts at element 1 and indexes up")
for element in element_checks:
    i = int(element-1)
    print("("+str(int(element))+") -> ["+str(np.round(t[i]*1e6, 2))+", "+str(g1_code_stream[i])+"]")
# how long to repeat a PRN code
duration_to_repeat = 1023 / chipping_rate
print("Duration for one G1 code to repeat is -> "+str(np.round(duration_to_repeat*1e6, 2))+" micro-seconds")
# plot results
plot.oscope_and_spectrum_analyser(t, g1_code_stream, "Signal Analysis for G1 Code", oscope_xlim_microsec=(0, 100), freq_xlim_kHz=(0, 5000))


print("\n|---------------- PART 5 ----------------|\n")
# create maximal length PRN9 code (tapped at 3 and 10)
PRN_code = seq.create_CA_code(1023, [3, 10])
# create the "stream" we would observe if this code was received at the provided
# chipping_rate, and we sampled at the provided f_sample_rate
chipping_rate = 1.023*1e6  # [Hz]
chip_duration = 1 / chipping_rate
chip_indices = np.floor(t / chip_duration).astype(int) % len(PRN_code)
PRN_stream = PRN_code[chip_indices]
# printing required elements from the array
element_checks = np.linspace(485, 495, 11)
print("[Time (micro-sec), Value (From PRN9)] for elements (485) to (495), where each "
      "\narray starts at element 1 and indexes up")
for element in element_checks:
    i = int(element-1)
    print("("+str(int(element))+") -> ["+str(np.round(t[i]*1e6, 2))+", "+str(PRN_stream[i])+"]")
# how long to repeat a PRN code
duration_to_repeat = 1023 / chipping_rate
print("Duration for one PRN code to repeat is -> "+str(np.round(duration_to_repeat*1e6, 2))+" micro-seconds")
# plot results
# plot.oscope_and_spectrum_analyser(t, PRN_stream, "Signal Analysis for PRN9", oscope_xlim_microsec=(0, 100), freq_xlim_kHz=(0, 5000))


print("\n|---------------- PART 6 ----------------|\n")
# provided
chipping_rate = 1.023*1e6  # [Hz]
chip_duration = 1 / chipping_rate
chip_indices = np.floor(t / chip_duration).astype(int) % len(g1_code)
g1_code_stream = g1_code[chip_indices]
f_carrier = 5.115 * 1e6  # 5.115 MHz
# convert the g1_code to +1 and -1s, create the carrier signal, and bpsk modulate
g1_bipolar = 2 * g1_code_stream - 1
carrier = np.sin(2 * np.pi * f_carrier * t)
bpsk_signal = g1_bipolar * carrier
# plot results
# plot.oscope_and_spectrum_analyser(t, bpsk_signal, str(f_carrier/1e6)+" MHz Carrier Wave", oscope_xlim_microsec=(10, 20), freq_xlim_kHz=(0, 10000))


print("\n|---------------- PART 7 ----------------|\n")
# create noisy bpsk signal
noise_std = 1
noise = np.random.normal(0, noise_std, size=bpsk_signal.shape)
bpsk_noisy = bpsk_signal + noise
# plot results
# plot.oscope_and_spectrum_analyser(t, bpsk_noisy, "("+str(f_carrier/1e6)+" MHz Carrier Wave + Noise)", oscope_xlim_microsec=(10, 20), freq_xlim_kHz=(0, 10000))
# multiply noisy signal by clean expected signal and plot
bpsk_multiply = bpsk_noisy * bpsk_signal
# plot.oscope_and_spectrum_analyser(t, bpsk_multiply, "("+str(f_carrier/1e6)+" MHz Carrier Wave + Noise) * "+str(f_carrier/1e6)+" MHz Carrier Wave", oscope_xlim_microsec=(10, 20))
plt.show()


