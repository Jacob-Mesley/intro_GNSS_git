# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant sequencing functions

# imports
import numpy as np
import matplotlib.pyplot as plt


# for getting simplified measurements of an aircraft overhead pass
def simplified_overhead_pass_measurements(x_a, v_a, h):
    """
    Calculating the range, range rate, and zenith angle of an overhead pass at a snapshot in time.
    This assumes the aircraft moves in a straight line tangent/parallel to the observer.

    :param x_a: sat horizontal position, wrt the observer
    :param v_a: sat velocity, wrt the observer
    :param h: sat altitude
    :return:
    """
    ac_range = np.sqrt((x_a**2 + h**2))
    ac_zenith_ang = np.arctan((x_a/h))
    ac_range_rate = (x_a * v_a) / ac_range
    return ac_range, ac_range_rate, ac_zenith_ang
