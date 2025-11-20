# Author: Jacob Mesley
# File Created: 9/2/2025
# Last Edit: 9/2/2025
# Description: file that stores all relevant plotting helper functions

# imports
import matplotlib.pyplot as plt
import numpy as np
import itertools

# simple plotting function given two arrays
def plot_arrays(x_data, y_data, x_axis_title="x-axis", y_axis_title="y-axis", title="plot title"):
    plt.figure(figsize=(12, 4))
    plt.plot(x_data, y_data)
    plt.grid(True)
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(title)
    plt.tight_layout()



def plot_multiple_arrays(x_data, y_data_matrix, x_axis_title="x-axis", y_axis_title="y-axis", title="plot title", labels=None):
    rows, cols = y_data_matrix.shape
    # Make figure wide
    plt.figure(figsize=(12, 4))
    # Define line styles to cycle through
    line_styles = ['-', '--', '-.', ':']
    style_cycle = itertools.cycle(line_styles)
    for i in range(cols):
        style = next(style_cycle)
        if labels is not None and i < len(labels):
            plt.plot(x_data, y_data_matrix[:, i], linestyle=style, label=labels[i])
        else:
            plt.plot(x_data, y_data_matrix[:, i], linestyle=style)
    plt.grid(True)
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(title)
    if labels is not None:
        plt.legend()
    plt.tight_layout()


def scatter_arrays(x_data, y_data, x_axis_title="x-axis", y_axis_title="y-axis", title="plot title"):
    plt.figure(figsize=(12, 2))
    plt.figure()
    plt.scatter(x_data, y_data)
    plt.grid(True)
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(title)
    plt.tight_layout()


def subplot_arrays(x_data_array, y_data_matrix, x_axis_label, y_axis_labels, title="plot title"):
    rows, cols = y_data_matrix.shape
    # create the subplots
    fig, axs = plt.subplots(cols, 1, figsize=(12, 2.0 * cols), sharex=True)
    axs = np.atleast_1d(axs)  # ensure axs is iterable even if rows==1
    for i in range(cols):
        axs[i].scatter(x_data_array, y_data_matrix[:, i])
        axs[i].set_ylabel(y_axis_labels[i])
        axs[i].grid(True)
    # Final formatting
    axs[0].set_title(title)
    axs[-1].set_xlabel(x_axis_label)
    plt.tight_layout()
    return fig, axs


def subplot_multiple_arrays(x_data_array, y_matrix_of_matrices, x_axis_label, y_axis_labels, line_labels=None, title="plot title"):
    """
    Create subplots with multiple lines in each subplot.

    Parameters
    ----------
    x_data_array : array-like
        Shared x-axis values.
    y_matrix_of_matrices : ndarray
        Shape (num_lines, N, num_subplots).
        - num_lines: number of lines to plot in each subplot
        - N: number of x points
        - num_subplots: number of subplots (vertical stack)
    x_axis_label : str
        Label for shared x-axis.
    y_axis_labels : list of str
        List of labels for each subplot's y-axis.
    line_labels : list of str, optional
        Labels for each line (shared across subplots).
    title : str
        Title for the figure.
    """
    rows, cols = y_matrix_of_matrices[0].shape
    num_lines = y_matrix_of_matrices.shape[0]

    # linestyles to cycle through
    linestyles = ['-', '--', '-.', ':']

    fig, axs = plt.subplots(cols, 1, figsize=(8, 2 * cols), sharex=True)
    axs = np.atleast_1d(axs)

    for i in range(cols):
        for j in range(num_lines):
            label = line_labels[j] if line_labels is not None else None
            style = linestyles[j % len(linestyles)]
            axs[i].plot(x_data_array, y_matrix_of_matrices[j, :, i], linestyle=style, label=label)
        axs[i].set_ylabel(y_axis_labels[i])
        axs[i].grid(True)
        if line_labels is not None:
            axs[i].legend()

    axs[0].set_title(title)
    axs[-1].set_xlabel(x_axis_label)
    plt.tight_layout()
    return fig, axs



# simple plotting function given two arrays
def plot_arrays_as_step(x_data, y_data, x_axis_title="x-axis", y_axis_title="y-axis", title="plot title"):
    plt.figure()
    n = len(x_data)
    plt.step(x_data, y_data)
    plt.grid(True)
    plt.xlim([0, n])
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(title)
    plt.tight_layout()


# simple plotting function given on array
def plot_array_as_step(array_to_plot, x_axis_label="x-axis", y_axis_label="y-axis", title="plot title"):
    plt.figure(figsize=(12, 3))
    plt.step(range(1, len(array_to_plot)+1), array_to_plot, where="post")
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()


def subplot_arrays_as_step(data_matrix, x_axis_label="x-axis", y_axis_label="y-axis", title: str = "plot title", y_limits=None):
    # determine number of measurements and steps
    rows, cols = data_matrix.shape
    x_axis_data = np.arange(1, cols + 1)
    # Create the subplots
    fig, axs = plt.subplots(rows, 1, figsize=(8, 3 * rows), sharex=True)
    axs = np.atleast_1d(axs)  # ensure axs is iterable even if rows==1
    for i in range(rows):
        array_to_plot = data_matrix[i, :]
        axs[i].step(range(1, len(array_to_plot)+1), array_to_plot, where="post")
        axs[i].set_ylabel(y_axis_label)
        axs[i].grid(True)
        if y_limits is not None:
            axs[i].set_ylim(y_limits)
    # Final formatting
    axs[0].set_title(title)
    axs[-1].set_xlabel(x_axis_label)
    plt.tight_layout()
    return fig, axs


def plot_az_el(az, el, title="plot title", svs=None, ax=None, ):
    """
    Creates an az-el sky plot of satellites.

    Parameters
    ----------
    az : array-like
        Vector of azimuth angles, in degrees.
    el : array-like
        Vector of elevation angles, in degrees.
    svs : array-like
        Vector of satellite PRN numbers.
        Use zeros to avoid printing PRN labels.
    ax : matplotlib Axes, optional
        Axes handle. If None, creates a new polar plot.

    Returns
    -------
    ax : matplotlib Axes
        The polar plot axes handle.
    """
    az = np.asarray(az)
    el = np.asarray(el)
    svs = np.asarray(svs)
    if az.shape != el.shape:
        raise ValueError("AZ and EL must be the same size.")
    if ax is None:
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.set_theta_zero_location("N")   # 0° at North
        ax.set_theta_direction(-1)        # clockwise azimuth
        ax.set_ylim(0, 90)
        ax.set_yticks(range(0, 91, 15))
        ax.set_yticklabels([str(90 - ang) for ang in range(0, 91, 15)])
        ax.grid(True)
    ax.set_title(title)
    # Convert elevation to polar radius (90° overhead → 0, horizon → 90)
    r = 90 - el
    theta = np.deg2rad(az)
    # Plot satellite positions
    ax.plot(theta, r, '.k', markersize=4)
    # Add PRN labels next to each satellite (if svs provided)
    if svs is not None:
        for th, rr, prn in zip(theta, r, svs):
            if prn != 0:
                # Slight offset so the text doesn't overlap the point
                ax.text(th, rr + 2, str(prn),
                        fontsize=8,
                        ha='center',
                        va='bottom')


def oscope_and_spectrum_analyser(time, signal_values, title="plot title", oscope_xlim_microsec=None, freq_xlim_kHz=None):
    # do spectrum analyser (do FFT)
    N = len(time)
    fs = 1 / (time[1] - time[0])
    freqs = fs * np.arange(0, N) / N
    py = 20 * np.log10(np.sqrt(2) * np.abs(np.fft.fft(signal_values) / N) + 1e-12)
    # create subplot
    fig, axs = plt.subplots(2, 1, figsize=(10, 6))
    # plot o-scope
    axs[0].plot(time*1e6, signal_values)
    axs[0].set_ylabel("Magnitude [V]")
    axs[0].set_xlabel("Time [micro-sec]")
    axs[0].set_title("Magnitude vs Time (Oscilloscope View)")
    if oscope_xlim_microsec:
        axs[0].set_xlim(oscope_xlim_microsec)
    axs[0].grid(True)
    # plot spectrum analyzer
    axs[1].plot(freqs/1e3, py)
    axs[1].set_ylabel("Magnitude [dBW]")
    axs[1].set_xlabel("Frequency [kHz]")
    axs[1].set_title("Magnitude vs Frequency (Spectrum Analyzer View, Sampling Rate of "+str(fs/1e6)+" MHz)")
    if freq_xlim_kHz:
        axs[1].set_xlim(freq_xlim_kHz)
    axs[1].grid(True)
    # finish plot
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    return fig, axs

def complex_corr(DOP_freq_span, delay_span, data, title="3D Correlation Surface", DOP_lims=None, delay_lims=None):
    # Mesh
    X, Y = np.meshgrid(DOP_freq_span, delay_span)

    # Figure
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Surface plot
    surf = ax.plot_surface(
        X, Y, data,
        cmap="viridis",
        edgecolor="none",
        linewidth=0,
        antialiased=True
    )
    # axis labels
    ax.set_xlabel("Doppler Frequency (Hz)")
    ax.set_ylabel("Code Delay (chips)")
    ax.set_zlabel("Correlation Magnitude")
    ax.set_title(title)

    # Axis limits (optional)
    if DOP_lims is not None:
        ax.set_xlim(DOP_lims)
    if delay_lims is not None:
        ax.set_ylim(delay_lims)
    # Colorbar
    fig.colorbar(surf, shrink=0.6, label="Correlation Value")
    plt.tight_layout()

    return fig, ax
