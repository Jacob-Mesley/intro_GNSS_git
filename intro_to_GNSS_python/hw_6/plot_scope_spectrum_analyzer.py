import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def plot_scope_spectrum_analyzer(t, y, zoomt=None, zoomf=None, zoomp=None):
    """
    Dawson Beatty, based on matlab script by Penina Axelrad
    Generates plots in the time domain and frequency domain of a signal
    Both are generated at the full scale and at a zoomed in time scale
    specified by the user, for a total of 4 subplots

    Inputs:
    t - time vector (assumed to be evenly spaced at the sample time)
    y - signal (amplitude) that you are analyzing
    zoomt - is max value of time for the second zoomed-in plot
    zoomf - [min, max] value of frequency for the second zoomed in plot
    zoomp - [min, max] value of power for plots

    Example call:
    t = np.arange(0, 10, 0.001)
    y = np.sin(2 * np.pi * 3 * t)
    plot_scope_spectrum_analyzer(t, y, 2, [0, 5], [-100, 0])
    """

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 10),
        constrained_layout=True,
    )

    # Time Plots - Signal Amplitude vs Time
    # Plots on the left are at full scale, right is zoomed in
    axes[0, 0].plot(t, y)
    axes[0, 0].set(
        xlabel="TIME (SEC)",
        ylabel="AMPLITUDE (V)",
        xlim=[np.min(t), np.max(t)],
    )
    axes[0, 0].grid()

    axes[0, 1].plot(t, y)
    axes[0, 1].set(
        xlabel="TIME (SEC)",
        xlim=[np.min(t), np.max(t)],
    )
    axes[0, 1].grid()
    if zoomt:
        axes[0, 1].set_xlim([0, zoomt])

    # Power Plots - Signal Power vs Frequency
    N = len(t)
    fs = 1 / (t[1] - t[0])
    f = fs * np.arange(0, N) / N

    py = 20 * np.log10(np.sqrt(2) * np.abs(np.fft.fft(y) / N))

    axes[1, 0].plot(f, py, marker=".")
    axes[1, 0].set(xlabel="FREQ (Hz)", ylabel="POWER (dBW)", xlim=[0, np.max(f) / 2])
    axes[1, 0].grid()

    axes[1, 1].plot(f, py, marker=".")
    axes[1, 1].set(xlabel="FREQ (HZ)", xlim=[0, np.max(f) / 2])
    axes[1, 1].grid()

    if zoomf:
        axes[1, 1].set_xlim(zoomf)

    if zoomp:
        axes[1, 1].set_ylim(zoomp)

    plt.show()
