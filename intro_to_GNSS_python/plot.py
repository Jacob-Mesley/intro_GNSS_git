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