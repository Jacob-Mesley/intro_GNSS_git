# Author: Jacob Mesley
# File Created: 9/4/2025
# Last Edit: 9/4/2025
# Description: file that stores all relevant conversion operations

# imports
import numpy as np


def binary_array_to_hex(bit_array):
    # make sure it's a flat array of ints
    bit_array = np.array(bit_array, dtype=int).flatten()
    # join into a string of bits
    bit_string = ''.join(str(b) for b in bit_array)
    # convert binary string -> int -> hex
    hex_str = hex(int(bit_string, 2))[2:].upper()
    return hex_str


# converts 0 -> 1 and 1 -> -1 in an array of binary
def binary_array_to_1s(binary_array):
    converted_binary = []
    for i in range(len(binary_array)):
        if binary_array[i] == 0:
            converted_binary.append(1)
        elif binary_array[i] == 1:
            converted_binary.append(-1)
    return converted_binary
