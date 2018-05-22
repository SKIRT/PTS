#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log
from pts.core.units.unit import PhotometricUnit
from pts.core.units.parsing import parse_unit as u
from pts.core.tools import numbers
from pts.core.units.parsing import parse_quantity
from pts.core.tools import parsing
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

description = "testing binary real conversions"

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W") # 3.828×10^26 W

# -----------------------------------------------------------------

class BinaryTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BinaryTest, self).__init__(*args, **kwargs)
        
    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Test dust mass
        self.test_dust_mass()

        # 3.
        self.test_dust_mass_with_rounding()

        # 4.
        self.test_dust_mass_second()

        # 5.
        self.test_speed_of_light()

        # 6. Experimental
        self.test_speed_of_light_experimental()

        # 6.
        self.test_gray_generation()

        # 7.
        self.test_gray_conversion_dynamic()

        # 8.
        self.test_gray_conversion_fixed()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BinaryTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def test_dust_mass(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing the dust mass ...")

        ndigits = 4

        nbits = numbers.binary_digits_for_significant_figures(ndigits)

        print(str(nbits) + " bits for " + str(ndigits) + " digits")

        mass_range = parsing.quantity_range("1500000.0 Msun > 300000000.0 Msun")

        unit = "Msun"

        minimum = mass_range.min
        maximum = mass_range.max

        low = minimum.to(unit).value
        high = maximum.to(unit).value

        value = parse_quantity("1.5e7 Msun").to(unit).value

        # Test : ROUNDING IN TEST BUT NOT IN BETWEEN CONVERSION!!
        if light_test_from_number_rounding(value, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_dust_mass_with_rounding(self):

        """
        This fucntion ...
        :return: 
        """

        # Inform the user
        log.info("Testing dust mass with rounding ...")

        ndigits = 4

        nbits = numbers.binary_digits_for_significant_figures(ndigits)

        #print(str(nbits) + " bits for " + str(ndigits) + " digits")

        mass_range = parsing.quantity_range("1500000.0 Msun > 300000000.0 Msun")

        unit = "Msun"

        minimum = mass_range.min
        maximum = mass_range.max

        low = minimum.to(unit).value
        high = maximum.to(unit).value

        # (FROM ERROR RESULTS IN MODELING:)

        # WITH ROUNDING:

        binary_string_a = [1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1]
        binary_string_b = [1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0]

        if test_from_binary_string_rounding(binary_string_a, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        if test_from_binary_string_rounding(binary_string_b, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        # Convert to numbers
        number_a = numbers.binary_string_to_float(binary_string_a, low, high, nbits)
        number_b = numbers.binary_string_to_float(binary_string_b, low, high, nbits)

        if test_from_number_rounding(number_a, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        if test_from_number_rounding(number_b, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_dust_mass_second(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Second dust mass test ...")

        ndigits = 4

        genome = [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]
        check = [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        nbits = len(genome)

        print(nbits)

        unit = "Msun"

        mass_range = parsing.quantity_range("1500000.0 Msun > 300000000.0 Msun")

        minimum = mass_range.min
        maximum = mass_range.max

        low = minimum.to(unit).value
        high = maximum.to(unit).value

        if test_from_binary_string_rounding(genome, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        if test_from_binary_string_rounding(check, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        genome_mass = numbers.binary_string_to_float(genome, low, high, nbits)

        rounded_genome_mass = numbers.round_to_n_significant_digits(genome_mass, ndigits)

        genome_mass_to_binary = numbers.float_to_binary_string(genome_mass, low, high, nbits)
        rounded_genome_mass_to_binary = numbers.float_to_binary_string(rounded_genome_mass, low, high, nbits)

        print("original:", genome_mass)
        print("rounded:", rounded_genome_mass)
        print("genome:", genome)
        print("original to binary:", genome_mass_to_binary)
        print("rounded to binary: ", rounded_genome_mass_to_binary)
        print("check:             ", check)

    # -----------------------------------------------------------------

    def test_speed_of_light(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing speed of light ...")

        ndigits = 3

        nbits = numbers.nbits_for_ndigits(ndigits)

        print(str(nbits) + " bits for " + str(ndigits) + " digits")

        minimum = 0.01 * speed_of_light
        maximum = speed_of_light

        unit = "km/s"

        # Convert value to binary
        value = (0.333 * speed_of_light).to(unit).value

        low = minimum.to(unit).value
        high = maximum.to(unit).value

        #print(value)
        #print(low)
        #print(high)

        # Test without rounding IN CONVERSION, BUT FOR COMPARISON THERE IS ROUNDING
        if light_test_from_number_rounding(value, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

        # With rounding
        if test_from_number_rounding(value, low, high, ndigits): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_speed_of_light_experimental(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing speed of light ...")

        ndigits = 3

        #nbits = numbers.nbits_for_ndigits(ndigits)

        minimum = 0.01 * speed_of_light
        maximum = speed_of_light

        unit = "km/s"

        # Convert value to binary
        value = (0.333 * speed_of_light).to(unit).value

        low = minimum.to(unit).value
        high = maximum.to(unit).value

        nbits = numbers.nbits_for_ndigits_experimental(ndigits, low, high)

        print(str(nbits) + " bits for " + str(ndigits) + " digits")

        # Test without rounding IN CONVERSION, BUT FOR COMPARISON THERE IS ROUNDING
        if light_test_from_number_rounding(value, low, high, ndigits, experimental=True): log.success("Test succeeded")
        else: log.error("Test failed")

        # With rounding
        if test_from_number_rounding(value, low, high, ndigits, experimental=True): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_gray_generation(self):

        """
        This function ...
        :return: 
        """

        #n = int(raw_input())

        # Inform the user
        log.info("Testing the generation of Gray code ...")

        n = 14

        if check_gray_generation(n): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_gray_conversion_dynamic(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing the Gray code conversions (dynamic nbits) ...")

        max_number = 15

        if check_gray_conversion(max_number): log.success("Test succeeded")
        else: log.error("Test failed")

    # -----------------------------------------------------------------

    def test_gray_conversion_fixed(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Testing the Gray code conversions (fixed nbits) ...")

        max_number = 15

        nbits = numbers.nbits_for_integer(max_number)

        if check_gray_conversion(max_number, nbits): log.success("Test succeeded")
        else: log.error("Test failed")

# -----------------------------------------------------------------

def test_from_binary_string_rounding(binary_string, low, high, ndigits):

    """
    This function ...
    :param binary_string: 
    :param low: 
    :param high: 
    :param ndigits: 
    :return: 
    """

    # Detrmine number of bits
    nbits = numbers.binary_digits_for_significant_figures(ndigits)

    # To number
    number = numbers.binary_string_to_float(binary_string, low, high, nbits)

    # Test from number
    return test_from_number_rounding(number, low, high, ndigits)

# -----------------------------------------------------------------

def test_from_binary_rounding(binary, low, high, ndigits):

    """
    This function ...
    :return: 
    """

    # Detrmine number of bits
    nbits = numbers.binary_digits_for_significant_figures(ndigits)

    # To number
    number = numbers.binary_to_float(binary, low, high, nbits)

    # Test from number
    return test_from_number_rounding(number, low, high, ndigits)

# -----------------------------------------------------------------

def test_from_number_rounding(number, low, high, ndigits, experimental=False):

    """
    This function ...
    :param number:
    :param low:
    :param high:
    :param ndigits:
    :param experimental:
    :return: 
    """

    # Detrmine number of bits
    if experimental: nbits = numbers.nbits_for_ndigits_experimental(ndigits, low, high)
    else: nbits = numbers.binary_digits_for_significant_figures(ndigits)

    # Number a to binary
    number_binary_string = numbers.float_to_binary_string(number, low, high, nbits)

    # back to number
    number_back = numbers.binary_string_to_float(number_binary_string, low, high, nbits)

    #print("")
    #print("SHOULD BE OK:")
    #print(number)
    #print(number_back)
    #print("")

    # NOW ROUND TO THE NUMBER OF DIGITS
    number_rounded = numbers.round_to_n_significant_digits(number, ndigits)

    #print("FAILS?")
    #print(number)
    ##print(number_rounded)

    # NOW CONVERT
    number_rounded_binary_string = numbers.float_to_binary_string(number_rounded, low, high, nbits)

    # Back to number
    number_rounded_back = numbers.binary_string_to_float(number_rounded_binary_string, low, high, nbits)

    #print(number_rounded_back)
    #print("")

    # Round
    number_rounded = numbers.round_to_n_significant_digits(number, ndigits)
    number_rounded_back_rounded = numbers.round_to_n_significant_digits(number_rounded_back, ndigits)

    # Check
    #return number_rounded_back == number

    return number_rounded == number_rounded_back_rounded

# -----------------------------------------------------------------

def light_test_from_number_no_rounding(number, low, high, ndigits):

    """
    This function ...
    :param number: 
    :param low: 
    :param high: 
    :param ndigits: 
    :return: 
    """

    # ACTUALLY THIS FUNCTION CAN BE EXPECTED TO FAIL!!

    # Detrmine number of bits
    nbits = numbers.binary_digits_for_significant_figures(ndigits)

    number_back = number_to_binary_and_back(number, low, high, nbits)

    largest_error = 0.5 * (high - low) / numbers.nintegers_for_nbits(nbits)

    if number_back != number:
        print(str(number) + " != " + str(number_back))
        abserror = abs(number_back - number)
        rel = abs((number_back - number) / number)
        print("rel error:", rel * 100, "%")
        print("abs error:", abserror)
        print("largest error:", largest_error)
    else: print(str(number) + " == " + str(number_back))

    return number_back == number

# -----------------------------------------------------------------

def light_test_from_number_rounding(number, low, high, ndigits, experimental=False):

    """
    This fucntion ...
    :param number: 
    :param low: 
    :param high: 
    :param ndigits: 
    :param experimental:
    :return: 
    """

    # Detrmine number of bits
    if experimental: nbits = numbers.nbits_for_ndigits_experimental(ndigits, low, high)
    else: nbits = numbers.binary_digits_for_significant_figures(ndigits)

    number_back = number_to_binary_and_back(number, low, high, nbits)

    largest_error = 0.5 * (high - low) / numbers.nintegers_for_nbits(nbits)

    number_rounded = numbers.round_to_n_significant_digits(number, ndigits)
    number_back_rounded = numbers.round_to_n_significant_digits(number_back, ndigits)

    if number_back_rounded != number_rounded:
        print(str(number_rounded) + " != " + str(number_back_rounded))
        abserror = abs(number_back_rounded - number_rounded)
        rel = abserror / number_rounded
        print("rel error:", rel * 100, "%")
        print("abs error:", abserror)
        print("largest error:", largest_error)
    else: print(str(number_rounded) + " == " + str(number_back_rounded))

    # Check
    return number_back_rounded == number_rounded

# -----------------------------------------------------------------

def number_to_binary_and_back(number, low, high, nbits):

    """
    This function ...
    :param number:
    :param low:
    :param high:
    :param nbits:
    :return: 
    """

    # Number a to binary
    number_a_binary_string = numbers.float_to_binary_string(number, low, high, nbits)
    return numbers.binary_string_to_float(number_a_binary_string, low, high, nbits)

# -----------------------------------------------------------------

def check_gray_generation(dimension):

    """
    This function ...
    :param dimension: 
    :return: 
    """

    nbits = dimension

    codes = numbers.binary_gray_code(dimension)

    #if n >= 1:
    #    for i in range(len(g)): print(g[i])

    # Try some integers

    # nintegers_for_nbits
    nintegers = numbers.nintegers_for_nbits(nbits)

    # Generate some random integers
    nrandom = 100
    for integer in list(np.random.randint(0, nintegers, size=nrandom)):

        # Get binary string representation
        binary = numbers.integer_to_binary_string(integer, nbits=nbits)

        # Get Gray code
        gray = codes[integer]

        #print(gray)

        # Check whether Gray to binary gives same
        binary_via_gray = numbers.gray_binary_string_to_binary_string(gray)

        #print(binary, binary_via_gray)

        if binary_via_gray != binary: return False

    # Test succeeded
    return True

# -----------------------------------------------------------------

def check_gray_conversion(max_number, nbits=None):

    """
    This function ...
    :param max_number: 
    :param nbits:
    :return: 
    """

    for i in range(max_number + 1):

        binary = numbers.integer_to_binary_string(i, nbits=nbits)
        gray = numbers.binary_string_to_gray_binary_string(binary)
        binary2 = numbers.gray_binary_string_to_binary_string(gray)
        integer = numbers.binary_string_to_integer(binary2)

        #print('int:%2i -> bin:%12r -> gray:%12r -> bin:%12r -> int:%2i' % (i, binary, gray, binary2, integer))

        # OUTPUT:
        # int: 0 -> bin:         [0] -> gray:         [0] -> bin:         [0] -> int: 0
        # int: 1 -> bin:         [1] -> gray:         [1] -> bin:         [1] -> int: 1
        # int: 2 -> bin:      [1, 0] -> gray:      [1, 1] -> bin:      [1, 0] -> int: 2
        # int: 3 -> bin:      [1, 1] -> gray:      [1, 0] -> bin:      [1, 1] -> int: 3
        # int: 4 -> bin:   [1, 0, 0] -> gray:   [1, 1, 0] -> bin:   [1, 0, 0] -> int: 4
        # int: 5 -> bin:   [1, 0, 1] -> gray:   [1, 1, 1] -> bin:   [1, 0, 1] -> int: 5
        # int: 6 -> bin:   [1, 1, 0] -> gray:   [1, 0, 1] -> bin:   [1, 1, 0] -> int: 6
        # int: 7 -> bin:   [1, 1, 1] -> gray:   [1, 0, 0] -> bin:   [1, 1, 1] -> int: 7
        # int: 8 -> bin:[1, 0, 0, 0] -> gray:[1, 1, 0, 0] -> bin:[1, 0, 0, 0] -> int: 8
        # int: 9 -> bin:[1, 0, 0, 1] -> gray:[1, 1, 0, 1] -> bin:[1, 0, 0, 1] -> int: 9
        # int:10 -> bin:[1, 0, 1, 0] -> gray:[1, 1, 1, 1] -> bin:[1, 0, 1, 0] -> int:10
        # int:11 -> bin:[1, 0, 1, 1] -> gray:[1, 1, 1, 0] -> bin:[1, 0, 1, 1] -> int:11
        # int:12 -> bin:[1, 1, 0, 0] -> gray:[1, 0, 1, 0] -> bin:[1, 1, 0, 0] -> int:12
        # int:13 -> bin:[1, 1, 0, 1] -> gray:[1, 0, 1, 1] -> bin:[1, 1, 0, 1] -> int:13
        # int:14 -> bin:[1, 1, 1, 0] -> gray:[1, 0, 0, 1] -> bin:[1, 1, 1, 0] -> int:14
        # int:15 -> bin:[1, 1, 1, 1] -> gray:[1, 0, 0, 0] -> bin:[1, 1, 1, 1] -> int:15

        if integer != i: return False

    return True

# -----------------------------------------------------------------
