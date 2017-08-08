#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.parameters

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ...core.tools import types, numbers
from ...core.basics.log import log

# -----------------------------------------------------------------

def is_binary_values(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    for element in sequence:
        if element != 0 and element != 1: return False
    return True

# -----------------------------------------------------------------

def is_binary_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return is_binary_values(genome)

# -----------------------------------------------------------------

def is_real_values(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    for element in sequence:
        if not types.is_real_type(element): return False
    return True

# -----------------------------------------------------------------

def is_real_genome(genome):

    """
    This function ...
    :param genome:
    :return:
    """

    return is_real_values(genome)

# -----------------------------------------------------------------

def equal_individuals(individual_a, individual_b, rtol=1e-5, atol=1e-8):

    """
    This function ...
    :param individual_a:
    :param individual_b:
    :param rtol:
    :param atol:
    :return:
    """

    return equal_genomes(individual_a.genomeList, individual_b.genomeList, rtol=rtol, atol=atol)

# -----------------------------------------------------------------

def equal_genomes(genome_a, genome_b, rtol=1e-5, atol=1e-8, binary_parameters=None, return_comparison=False):

    """
    This function ...
    :param genome_a: 
    :param genome_b: 
    :param rtol:
    :param atol:
    :param binary_parameters:
    :param return_comparison:
    :return: 
    """

    # Convert into arrays
    genome_a = np.array(genome_a)
    genome_b = np.array(genome_b)

    # Binary: check exact
    if is_binary_values(genome_a):

        # Parameters are specified for doing binary to real conversion
        if binary_parameters is not None:

            minima = binary_parameters.minima
            maxima = binary_parameters.maxima
            nbits = binary_parameters.nbits
            gray = binary_parameters.gray
            ndigits = binary_parameters.ndigits
            scales = binary_parameters.scales

            # Compare the binary genomes with conversion to real parameters
            return equal_binary_genomes_with_conversion(genome_a, genome_b, minima, maxima, nbits, gray, ndigits, scales, return_comparison=return_comparison)

        # Check whether the binary genomes are
        else:
            if return_comparison: return equal_binary_genomes_exact(genome_a, genome_b), genome_a, genome_b
            else: return equal_binary_genomes_exact(genome_a, genome_b)

    # Real: check with certain tolerance
    elif is_real_values(genome_b):

        if return_comparison: return np.isclose(genome_a, genome_b, rtol=rtol, atol=atol), genome_a, genome_b
        else: return np.isclose(genome_a, genome_b, rtol=rtol, atol=atol)

    # Unrecognized 1D genome list
    else: raise ValueError("Genome list not recognized: " + str(genome_a))

# -----------------------------------------------------------------

def equal_binary_genomes_exact(genome_a, genome_b):

    """
    This fucntion ...
    :param genome_a: 
    :param genome_b: 
    :return: 
    """

    return np.all(genome_a == genome_b)

# -----------------------------------------------------------------

def equal_binary_genomes_with_conversion(genome_a, genome_b, minima, maxima, nbits, gray, ndigits, scales, return_comparison=False):

    """
    This function ...
    :param genome_a: 
    :param genome_b: 
    :param minima:
    :param maxima:
    :param nbits:
    :param gray:
    :param ndigits:
    :param scales:
    :param return_comparison:
    :return: 
    """

    parameters_a_rounded = get_parameters_from_genome_rounded(genome_a, minima, maxima, nbits, scales, gray, ndigits)
    parameters_b_rounded = get_parameters_from_genome_rounded(genome_b, minima, maxima, nbits, scales, gray, ndigits)

    # Loop over the parameters
    # get_parameter_values_from_genome
    #for value_a_rounded, value_b_rounded in get_parameters_from_genomes_rounded(genome_a, genome_b, minima, maxima, nbits, gray, ndigits, scales):
    for value_a_rounded, value_b_rounded in zip(parameters_a_rounded, parameters_b_rounded):

        # Fail!
        if value_a_rounded != value_b_rounded:
            if return_comparison: return False, parameters_a_rounded, parameters_b_rounded
            else: return False

    # All checks passed
    if return_comparison: return True, parameters_a_rounded, parameters_b_rounded
    else: return True, parameters_a_rounded, parameters_b_rounded

# -----------------------------------------------------------------

def get_parameters_from_genome(genome, minima, maxima, nbits, parameter_scales, gray=False):

    """
    This function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :param parameter_scales:
    :param gray:
    :return:
    """

    # Scale the minima and maxima!
    minima, maxima = scale_extrema(minima, maxima, parameter_scales)

    # Get scaled form
    scaled = get_parameters_from_genome_scaled(genome, minima, maxima, nbits, gray=gray)

    # Unscale if necessary
    parameters = unscale_parameters(scaled, parameter_scales)

    # Return the parameters
    return parameters

# -----------------------------------------------------------------

def get_parameters_from_genome_rounded(genome, minima, maxima, nbits, parameter_scales, gray, ndigits):

    """
    THis function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :param parameter_scales:
    :param gray:
    :param ndigits:
    :return:
    """

    parameters = get_parameters_from_genome(genome, minima, maxima, nbits, parameter_scales, gray)

    # Loop over the parameters
    for index in range(len(parameters)):

        value = parameters[index]
        rounded = numbers.round_to_n_significant_digits(value, ndigits[index])

        # Replace the value
        parameters[index] = rounded

    # Return the list of rounded values
    return parameters

# -----------------------------------------------------------------

def get_parameters_from_genomes(genome_a, genome_b, minima, maxima, nbits, scales, gray):

    """
    This function ...
    :param genome_a:
    :param genome_b:
    :param minima:
    :param maxima:
    :param nbits:
    :param scales:
    :param gray:
    :return:
    """

    parameters_a = get_parameters_from_genome(genome_a, minima, maxima, nbits, scales, gray)
    parameters_b = get_parameters_from_genome(genome_b, minima, maxima, nbits, scales, gray)

    # Return the pairs
    return zip(parameters_a, parameters_b)

# -----------------------------------------------------------------

def get_parameters_from_genomes_rounded(genome_a, genome_b, minima, maxima, nbits, scales, gray, ndigits):

    """
    This function ...
    :param genome_a:
    :param genome_b:
    :param minima:
    :param maxima:
    :param nbits:
    :param scales:
    :param gray:
    :param ndigits:
    :return:
    """

    parameters_a = get_parameters_from_genome_rounded(genome_a, minima, maxima, nbits, scales, gray, ndigits=ndigits)
    parameters_b = get_parameters_from_genome_rounded(genome_b, minima, maxima, nbits, scales, gray, ndigits=ndigits)

    # Return the pairs
    return zip(parameters_a, parameters_b)

# -----------------------------------------------------------------

#def get_parameters_from_genomes(genome_a, genome_b, minima, maxima, nbits, gray):

    #"""
    #This function ...
    #:param genome_a:
    #:param genome_b:
    #:param minima:
    #:param maxima:
    #:param nbits:
    #:param gray:
    #:return:
    #"""

    #nparameters = len(minima)

    #bit_slices = numbers.generate_bit_slices(nbits)

    #pairs = []

    # Loop over the parameters
    #for index in range(nparameters):

        # Get the first n bits
        #bits_a = genome_a[bit_slices[index]]
        #bits_b = genome_b[bit_slices[index]]

        # Convert into real value
        #if gray:
        #    value_a = numbers.gray_binary_string_to_float(bits_a, low=minima[index], high=maxima[index], nbits=nbits[index])
        #    value_b = numbers.gray_binary_string_to_float(bits_b, low=minima[index], high=maxima[index], nbits=nbits[index])
        #else:
        #    value_a = numbers.binary_string_to_float(bits_a, low=minima[index], high=maxima[index], nbits=nbits[index])
        #    value_b = numbers.binary_string_to_float(bits_b, low=minima[index], high=maxima[index], nbits=nbits[index])

        #pairs.append((value_a, value_b))

    # Return the pairs
    #return pairs

# -----------------------------------------------------------------

#def get_parameters_from_genomes_rounded(genome_a, genome_b, minima, maxima, nbits, gray, ndigits):

    #"""
    #This function ...
    #:param genome_a:
    #:param genome_b:
    #:param minima:
    #:param maxima:
    #:param nbits:
    #:param gray:
    #:param ndigits:
    #:return:
    #"""

    #nparameters = len(minima)

    #bit_slices = numbers.generate_bit_slices(nbits)

    #pairs = []

    # Loop over the parameters
    #for index in range(nparameters):

        # Get the first n bits
        #bits_a = genome_a[bit_slices[index]]
        #bits_b = genome_b[bit_slices[index]]

        # Convert into real value
        #if gray:
        #    value_a = numbers.gray_binary_string_to_float(bits_a, low=minima[index], high=maxima[index], nbits=nbits[index])
        #    value_b = numbers.gray_binary_string_to_float(bits_b, low=minima[index], high=maxima[index], nbits=nbits[index])
        #else:
        #    value_a = numbers.binary_string_to_float(bits_a, low=minima[index], high=maxima[index], nbits=nbits[index])
        #    value_b = numbers.binary_string_to_float(bits_b, low=minima[index], high=maxima[index], nbits=nbits[index])

        #get_parameters_from_genomes()

        #value_a_rounded = numbers.round_to_n_significant_digits(value_a, ndigits[index])
        #value_b_rounded = numbers.round_to_n_significant_digits(value_b, ndigits[index])

        #pairs.append((value_a_rounded, value_b_rounded))

    # Return the pairs
    #return pairs

# -----------------------------------------------------------------

def parameters_to_binary_string(parameters, minima, maxima, nbits):

    """
    This function ...
    :param parameters:
    :param minima:
    :param maxima
    :param nbits:
    :return:
    """

    binary_string = []

    # Convert parameters into binary strings
    for index in range(len(parameters)):
        # Convert floating point value into binary string with specific number of bits
        value = parameters[index]
        binary_string_parameter = numbers.float_to_binary_string(value, low=minima[index], high=maxima[index],
                                                                 nbits=nbits[index])

        # Extend the complete string
        binary_string.extend(binary_string_parameter)

    # Return the binary string
    return binary_string

# -----------------------------------------------------------------

def parameters_to_gray_binary_string(parameters, minima, maxima, nbits):

    """
    This function ...
    :param parameters:
    :param minima:
    :param maxima:
    :param nbits:
    :return:
    """

    binary_string = []

    # Convert parameters into binary strings
    for index in range(len(parameters)):

        # Convert floating point value into Gray binary string with specific number of bits
        value = parameters[index]
        gray_binary_string_parameter = numbers.float_to_gray_binary_string(value, low=minima[index], high=maxima[index],
                                                                           nbits=nbits[index])

        # Extend the complete string
        binary_string.extend(gray_binary_string_parameter)

    # Return the binary string
    return binary_string

# -----------------------------------------------------------------

def binary_string_to_parameters(genome, minima, maxima, nbits):

    """
    This function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :return:
    """

    nparameters = len(minima)
    assert nparameters == len(maxima) == len(nbits)

    # Initialize list for the parameters
    parameters = []

    # Generate bit slices
    bit_slices = numbers.generate_bit_slices(nbits)

    # Loop over the parameters
    for index in range(nparameters):

        # Get the first n bits
        bits = genome[bit_slices[index]]

        # Convert into real value
        value = numbers.binary_string_to_float(bits, low=minima[index], high=maxima[index], nbits=nbits[index])

        # Add the value
        parameters.append(value)

    # Return the parameter values
    return parameters

# -----------------------------------------------------------------

def gray_binary_string_to_parameters(genome, minima, maxima, nbits):

    """
    This function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :return:
    """

    nparameters = len(minima)
    assert nparameters == len(maxima) == len(nbits)

    # Initialize list for the parameters
    parameters = []

    # Generate bit slices
    bit_slices = numbers.generate_bit_slices(nbits)

    # Loop over the parameters
    for index in range(nparameters):

        # Get the first n bits
        bits = genome[bit_slices[index]]

        # Convert into real value
        value = numbers.gray_binary_string_to_float(bits, low=minima[index], high=maxima[index], nbits=nbits[index])

        # Add the value
        parameters.append(value)

    # Return the parameter values
    return parameters

# -----------------------------------------------------------------

def get_parameters_from_binary_genome_scaled(genome, minima, maxima, nbits, gray=False):

    """
    This function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :param gray:
    :return:
    """

    # Convert
    if gray: parameters = gray_binary_string_to_parameters(genome, minima, maxima, nbits)
    else: parameters = binary_string_to_parameters(genome, minima, maxima, nbits)

    # Return
    return parameters

# -----------------------------------------------------------------

def get_parameters_from_genome_scaled(genome, minima, maxima, nbits, gray=False):

    """
    This function ...
    :param genome:
    :param minima:
    :param maxima:
    :param nbits:
    :param gray:
    :return:
    """

    # If binary genome, convert binary individual into actual parameters list
    if is_binary_genome(genome): parameters = get_parameters_from_binary_genome_scaled(genome, minima, maxima, nbits, gray=gray)
    elif is_real_genome(genome): parameters = list(genome)
    else: raise ValueError("Unrecognized genome: " + str(genome))

    # Return the parameters
    return parameters

# -----------------------------------------------------------------

def get_binary_genome_from_scaled_parameters(parameters, minima, maxima, nbits, gray=False):

    """
    This function ...
    :param parameters:
    :param minima:
    :param maxima:
    :param nbits:
    :param gray:
    :return:
    """

    # # Convert
    if gray: binary_string = parameters_to_gray_binary_string(parameters, minima, maxima, nbits)
    else: binary_string = parameters_to_binary_string(parameters, minima, maxima, nbits)

    # Return the binary string genome
    return binary_string

# -----------------------------------------------------------------

def get_binary_genome_from_parameters(parameters, minima, maxima, nbits, parameter_scales, gray=False):

    """
    This function ...
    :param parameters:
    :param minima:
    :param maxima:
    :param nbits:
    :param parameter_scales:
    :param gray:
    :return:
    """

    # Scale
    scaled_parameters = scale_parameters(parameters, parameter_scales)

    # Scale the minima and maxima
    minima, maxima = scale_extrema(minima, maxima, parameter_scales)

    # Return the binary genome
    return get_binary_genome_from_scaled_parameters(scaled_parameters, minima, maxima, nbits, gray=gray)

# -----------------------------------------------------------------

def scale_parameters(parameters, scales):

    """
    This function ...
    :param parameters:
    :param scales:
    :return:
    """

    new_parameters = []

    # Loop over the scales and parameter values
    for scale, parameter_value in zip(scales, parameters):

        # Convert to scale
        value = numbers.to_scale(parameter_value, scale)

        # Add the value to the parameters
        new_parameters.append(value)

    # Return the new parameter set
    return new_parameters

# -----------------------------------------------------------------

def unscale_parameters(parameters, scales):

    """
    This function ...
    :param parameters:
    :param scales:
    :return:
    """

    new_parameters = []

    # Loop over the scales and scaled parameter values
    for scale, scaled_value in zip(scales, parameters):

        # Convert to scale
        value = numbers.unscale(scaled_value, scale)

        # Add the value to the parameters
        new_parameters.append(value)

    # Return the new parameter set
    return new_parameters

# -----------------------------------------------------------------

def scale_extrema(minima, maxima, scales):

    """
    This function ...
    :param minima:
    :param maxima:
    :param scales:
    :return:
    """

    scaled_minima = []
    scaled_maxima = []

    # Loop over the values
    for minimum, maximum, scale in zip(minima, maxima, scales):

        # Convert to scale
        scaled_minimum = numbers.to_scale(minimum, scale)
        scaled_maximum = numbers.to_scale(maximum, scale)

        # Add to the lists
        scaled_minima.append(scaled_minimum)
        scaled_maxima.append(scaled_maximum)

    # Return the scaled values
    return scaled_minima, scaled_maxima

# -----------------------------------------------------------------

def unscale_extrema(scaled_minima, scaled_maxima, scales):

    """
    This function ...
    :param scaled_minima:
    :param scaled_maxima:
    :param scales:
    :return:
    """

    minima = []
    maxima = []

    # Loop over the values
    for scaled_minimum, scaled_maximum, scale in zip(scaled_minima, scaled_maxima, scales):

        # Convert to scale
        minimum = numbers.unscale(scaled_minima, scale)
        maximum = numbers.unscale(scaled_maxima, scale)

        # Add to the list
        minima.append(minimum)
        maxima.append(maximum)

    # Return the unscaled values
    return minima, maxima

# -----------------------------------------------------------------

def scale_parameter_sets(sets, scales):

    """
    This function:
    :param sets:
    :param scales
    :return:
    """

    # Initialize a new parameter set list
    parameter_sets = []

    # Loop over the parameter sets
    for original_parameters in sets:

        # Scale
        parameters = scale_parameters(original_parameters, scales)

        # Add the parameter set to the list of parameter sets
        parameter_sets.append(parameters)

    # Return the new parameter sets
    return parameter_sets

# -----------------------------------------------------------------

def unscale_parameter_sets(sets, scales):

    """
    This function ...
    :param sets:
    :param scales:
    :return:
    """

    # Initialize a new parameter set list
    parameter_sets = []

    # Loop over the parameter sets
    for scaled_parameters in sets:
        # Scale
        parameters = unscale_parameters(scaled_parameters, scales)

        # Add the parameter set to the list of parameter sets
        parameter_sets.append(parameters)

    # Return the new parameter sets
    return parameter_sets

# -----------------------------------------------------------------

def scale_extrema_sets(sets, scales):

    """
    This function ...
    :param sets:
    :param scales:
    :return:
    """

    extrema_sets = []

    # Loop over
    for minima, maxima in sets:

        # Scale
        scaled_minima, scaled_maxima = scale_extrema(minima, maxima, scales)

        # Add the extrema
        extrema_sets.append((scaled_minima, scaled_maxima))

    # Return the new extrema
    return extrema_sets

# -----------------------------------------------------------------

def unscale_extrema_sets(sets, scales):

    """
    This function ...
    :param sets:
    :param scales:
    :return:
    """

    extrema_sets = []

    # Loop over
    for scaled_minima, scaled_maxima in sets:

        # Unscale
        minima, maxima = unscale_extrema(scaled_minima, scaled_maxima, scales)

        # Add the extrema
        extrema_sets.append((minima, maxima))

    # Return the new extrema
    return extrema_sets

# -----------------------------------------------------------------

def round_parameters(parameters, ndigits):

    """
    This function ...
    :param parameters:
    :param ndigits:
    :return:
    """

    rounded_parameters = []

    # Loop over the parameters
    for index in range(len(parameters)):
        # Get the value
        value = parameters[index]

        # Convert to relevant number of digits
        value = numbers.round_to_n_significant_digits(value, ndigits[index])

        # Add to the list
        rounded_parameters.append(value)

    # Return
    return rounded_parameters

# -----------------------------------------------------------------

# Convert into arrays
def show_genome_differences(genome_a, genome_b, rtol=1e-5, atol=1e-8, binary_parameters=None):

    """
    This function ...
    :param genome_a: 
    :param genome_b:
    :param binary_parameters:
    :return: 
    """

    genome_a = np.array(genome_a)
    genome_b = np.array(genome_b)

    # Binary: check exact
    if is_binary_values(genome_a):

        # Parameters are specified for doing binary to real conversion
        if binary_parameters is not None:

            minima = binary_parameters.minima
            maxima = binary_parameters.maxima
            nbits = binary_parameters.nbits
            gray = binary_parameters.gray
            ndigits = binary_parameters.ndigits

            parameter_index = 0
            for value_a, value_b in get_parameters_from_genomes(genome_a, genome_b, minima, maxima, nbits, gray):

                # Round both values
                value_a_rounded = numbers.round_to_n_significant_digits(value_a, ndigits[parameter_index])
                value_b_rounded = numbers.round_to_n_significant_digits(value_b, ndigits[parameter_index])

                adiff = abs(value_a_rounded - value_b_rounded)
                rdiff = adiff / value_a_rounded

                if adiff > atol or rdiff > rtol:

                    log.error("Different (rounded) values for parameter #" + str(parameter_index) + ":")
                    log.error(" - value_a: " + str(value_a))
                    log.error(" - value_b: " + str(value_b))
                    log.error(" - value_a_rounded: " + str(value_a_rounded))
                    log.error(" - value_b_rounded: " + str(value_b_rounded))
                    log.error(" - absolute difference: " + str(adiff))
                    log.error(" - relative difference: " + str(rdiff * 100.) + "%")
                    log.error(" - absolute tolerance: " + str(atol))
                    log.error(" - relative tolerance: " + str(rtol * 100.) + "%")

                parameter_index += 1

        else:

            # Show the positions for which there is a mismatch
            indices = np.where(genome_a != genome_b)
            log.error(" - different bits: " + ",".join(str(index) for index in indices))

    # Real: check with certain tolerance
    elif is_real_values(genome_b):

        adiff = abs(genome_a - genome_b) # numpy arrays
        rdiff = adiff / genome_a

        log.error(" - absolute differences: " + str(adiff))
        log.error(" - relative differences: " + str(rdiff * 100.) + "%")
        log.error(" - absolute tolerance: " + str(atol))
        log.error(" - relative tolerance: " + str(rtol * 100.) + "%")

    # Unrecognized 1D genome list
    else: raise ValueError("Genome list not recognized: " + str(genome_a))

# -----------------------------------------------------------------
