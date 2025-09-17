# Author: Jacob Mesley
# File Created: 9/2/2025
# Last Edit: 9/2/2025
# Description: file that stores all relevant plotting helper functions

# imports
import matplotlib.pyplot as plt
import numpy as np


# simple plotting function given two arrays
def plot_arrays(x_data, y_data, x_axis_title="x-axis", y_axis_title="y-axis", title="plot title"):
    plt.figure()
    plt.plot(x_data, y_data)
    plt.grid(True)
    plt.xlabel(x_axis_title)
    plt.ylabel(y_axis_title)
    plt.title(title)
    plt.tight_layout()


def subplot_arrays(x_data_array, y_data_matrix, x_axis_label, y_axis_labels, title="plot title"):
    rows, cols = y_data_matrix.shape
    # create the subplots
    fig, axs = plt.subplots(cols, 1, figsize=(8, 2 * cols), sharex=True)
    axs = np.atleast_1d(axs)  # ensure axs is iterable even if rows==1
    for i in range(cols):
        axs[i].scatter(x_data_array, y_data_matrix[:, i])
        axs[i].set_ylabel(y_axis_labels[i])
        axs[i].grid(True)
    # Final formatting
    axs[0].set_title(title)
    axs[-1].set_xlabel(x_axis_label)
    axs[i].set_xlim(np.min(x_data_array)-1, np.max(x_data_array)+1)
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

    # # Add PRN labels
    # if svs is not None:
    #     for i, prn in enumerate(svs):
    #         if prn != 0:
    #             ax.text(theta[i], r[i], str(prn), fontsize=8,
    #                     ha='left', va='center')
    #
    #     return ax
