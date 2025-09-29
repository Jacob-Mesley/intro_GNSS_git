# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant sequencing functions

# imports
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
import operator

# bitwise xor function
def xor(array_of_bits):
    return reduce(operator.xor, array_of_bits)


# one step of a shift register
def shift_register_step(reg, taps_zero_indexed):
    output = reg[-1]
    # xor the register at the tap locations
    bits_to_xor = []
    for tap in taps_zero_indexed:
        bits_to_xor.append(reg[tap])
    xor_result = xor(bits_to_xor)
    # shift the register (using numpy functions)
    reg = np.roll(reg, 1)
    reg[0] = xor_result
    # return the new register, and the "output" of the shift
    return reg, output


# for creating a CA code using a specific PRN tap
def create_CA_code(code_length, PRN_taps_one_indexed):
    # decompose the PRN taps
    PRN_tap1 = PRN_taps_one_indexed[0] - 1
    PRN_tap2 = PRN_taps_one_indexed[1] - 1
    # creating the G1 10 bit shift register
    G1 = np.ones(10, dtype=int)
    G1_taps = [2, 9]  # this is zero index. One indexed would be [3, 10]
    # creating the G2 10 bit shift register
    G2 = np.ones(10, dtype=int)
    G2_taps = [1, 2, 5, 7, 8, 9]  # this is zero index. One indexed would be [2, 3, 6, 8, 9, 10]
    # initialize the CA code length
    CA_code = np.zeros(code_length, dtype=int)

    # generate the CA code
    for i in range(code_length):
        # generate the new bit in the CA code
        G1_output = G1[-1]
        G2i = xor([G2[PRN_tap1], G2[PRN_tap2]])
        new_bit = xor([G2i, G1_output])
        CA_code[i] = new_bit
        # shift the registers
        G1, _ = shift_register_step(G1, G1_taps)
        G2, _ = shift_register_step(G2, G2_taps)
    # return resulting CA code
    return CA_code


# for comparing the difference between two arrays
def get_mismatching_indices(code_1, code_2):
    mismatch_index_array = []
    for i in range(len(code_1)):
        if code_1[i] != code_2[i]:
            mismatch_index_array.append(i)
    return mismatch_index_array


# given matlab function converted to python
def cyc_corr_basic(x, y):
    """
    Function to compute the cyclic correlation
    First input is the "received" signal (x)
    Second input is the "replica" signal (y) - the one that shifts
    """
    n = len(x)
    x = np.reshape(x, (1, n))
    yshifted = np.reshape(y, (n, 1))
    lag = np.arange(0, n+1)
    Rxy = np.full_like(lag, np.nan, dtype=float)
    for i in range(len(lag)):
        # dot product acts as multiplication of every element, then summation of every element
        # this matches the given equation in hw 1, prob 2-2
        Rxy[i] = np.dot(x, yshifted)
        # circular shift down by 1
        yshifted = np.roll(yshifted, 1, axis=0)
    return Rxy, lag


# function for shifting and wrapping arrays
# a pos number means a delay (shift right)
# negative number means an advance (shift left)
def shift_and_wrap_array(array_to_shift, number_to_shift_by):
    return np.roll(array_to_shift, number_to_shift_by)


def extract_rows_by_index(matrix, index, chunk_size=32):
    """
    Extracts the same row index from each chunk of rows in a 2D array.

    Args:
        matrix (np.ndarray): Input 2D array of shape (N, M).
        index (int): Row index to extract within each chunk.
        chunk_size (int): Number of rows per chunk (default=32).

    Returns:
        np.ndarray: New 2D array containing the extracted rows.
    """
    n_chunks = matrix.shape[0] // chunk_size
    extracted = [matrix[i*chunk_size + index] for i in range(n_chunks)]
    return np.array(extracted)


def extract_PRN_from_sp3_data(sp3_data, PRN_num):
    PRN_num = float(PRN_num)
    mask = sp3_data[:, 2].astype(str) == str(PRN_num)
    return sp3_data[mask]

def extract_PRN_from_ephem(ephem_data, PRN_num):
    PRN_num = float(PRN_num)
    mask = ephem_data[:, 0].astype(str) == str(PRN_num)
    return ephem_data[mask]
