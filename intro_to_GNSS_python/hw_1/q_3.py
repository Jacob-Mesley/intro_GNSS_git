# work for question 3 of the homework
import numpy as np
import matplotlib.pyplot as plt
import plot
import convert
import sequencing as seq


# for formatting and calculating other functions within my code
def calculate_and_print_CA_code_for_PRN(number_of_bits, PRN_taps, PRN_number):
    print("------------------------")
    # creating a PRN CA code and converting to hexadecimal
    PRN_CA_code = seq.create_CA_code(number_of_bits, [PRN_taps[0], PRN_taps[1]])
    PRN_first_16_hexa = convert.binary_array_to_hex(PRN_CA_code[0:16])
    PRN_last_16_hexa = convert.binary_array_to_hex(PRN_CA_code[-16:])
    print("First and Last 16 bits of the CA code PRN", PRN_number, ", 1023 bits long \nconverted to Hexadecimal =",
          PRN_first_16_hexa, ",", PRN_last_16_hexa)
    # plotting the first and last 16 of the generated sequence
    first_16_bits = PRN_CA_code[0:16]
    last_16_bits = PRN_CA_code[-16:]
    bit_matrix = np.array([first_16_bits.T, last_16_bits.T])
    plot.subplot_arrays_as_step(bit_matrix, "Bits", "Bit Value",
                                "First 16 (Top) and Last 16 (Bottom) Bits in 1023 Bit Long C/A Code (PRN-" + PRN_number + ")")

# PART A --------------------------------------------------------------------------------------------

# creating a PRN 19 CA code and converting to hexadecimal
calculate_and_print_CA_code_for_PRN(1023, [3, 6], "19")


# PART B -------------------------------------------------------------------------------------------

print("------------------------")
# create standard CA code
PRN19_CA_code = seq.create_CA_code(1023, [3, 6])
# create a C/A code that is twice as long
PRN19_CA_code_long = seq.create_CA_code(1023*2, [3, 6])
# take last half of newly generated code
PRN19_CA_code_second_half = PRN19_CA_code_long[-1023:]
# search for the number of indices that are not matching
mismatching_indices = seq.get_mismatching_indices(PRN19_CA_code, PRN19_CA_code_second_half)
number_of_mismatches = len(mismatching_indices)
print("Number of mis-matching indices between the 0-1023 C/A code and \nthe 1024-2046 is", number_of_mismatches)


# PART C -------------------------------------------------------------------------------------------

# creating a PRN 25 CA code and converting to hexadecimal
calculate_and_print_CA_code_for_PRN(1023, [5, 7], "25")


# PART D -------------------------------------------------------------------------------------------

# creating a PRN 5 CA code and converting to hexadecimal
calculate_and_print_CA_code_for_PRN(1023, [1, 9], "5")

# show all plots
plt.show()


