# Author: Jacob Mesley
# File Created: 9/5/2025
# Last Edit: 9/5/2025
# Description: file that stores all relevant data/measurement generation functions

# imports
import numpy as np
import calculate

def noise_vector(vector_length, std_dev):
    """
    Generate a noise vector of length n, normally distributed
    with zero mean and given standard deviation.
    """
    return np.random.normal(loc=0.0, scale=std_dev, size=vector_length)