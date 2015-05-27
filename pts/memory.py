#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# This ad-hoc script estimates the amount of memory needed by a particular SKIRT simulation.
# The current implementation assumes that
#  - the stellar system does not consume a lot of memory (i.e. no SPH or AMR)
#  - there is a dust system

# -----------------------------------------------------------------

import os.path

# -----------------------------------------------------------------

## Return the estimated memory usage for a given ski file (in gigabytes)
def estimate_memory(skifile, inputpath=None):

    # Get the number of wavelengths (force all calculations to floating point, avoiding integer overrun)
    Nlambda = None
    try:

        Nlambda = float(skifile.nwavelengths())

    except ValueError:

        # Get the name of the wavelength data file
        filename = skifile.wavelengthsfile()

        # Open the wavelengths file
        filepath = os.path.join(inputpath, filename)
        with open(filepath) as file:

            first = file.readlines()[0]

        Nlambda = float(first.split(" ")[0])

    # Get the number of dust cells
    try:
        Ncells = skifile.ncells()
    except ValueError:

        raise ValueError("Disabled for now")
        #Ncells = float(raw_input("\033[91m" + "\033[1m" + "?  How many dust cells do you expect for this simulation? ... \033[0m"))

    # Get the number of dust components
    Ncomps = skifile.ncomponents()

    # Get the number of items in the dust library
    Nitems = skifile.nlibitems()

    # Get the number of dust populations (all dust mixes combined)
    Npops = skifile.npopulations()

    # Get other simulation properties
    dustEmission = skifile.dustemission()
    selfAbsorption = skifile.dustselfabsorption()
    transientHeating = skifile.transientheating()

    # Overhead
    Ndoubles = 50e6 + (Nlambda + Ncells + Ncomps + Npops) * 10

    # Instruments
    for instrument in skifile.npixels(Nlambda):

        Ndoubles += instrument[2]

    # Dust system
    Ndoubles += (Ncomps+1) * Ncells
    if dustEmission:
        Ndoubles += Nlambda * Ncells
        if selfAbsorption:
            Ndoubles += Nlambda * Ncells

    # Dust grid tree
    if True:

        Ndoubles += Ncells * 1.2 * 50

    # Dust mixes
    Ndoubles += Nlambda * Npops * 3
    Ndoubles += (Nlambda + Npops) * 10

    # Dust library
    if dustEmission:
        Ndoubles += Ncells + Nlambda*Nitems

    # Transient heating
    if dustEmission and transientHeating:

        NT = 2250.
        Ndoubles += Ncomps * (Nlambda+1) * NT
        Ndoubles += Npops * 5./8.*NT*NT

    # 8 bytes in a double
    Nbytes = Ndoubles * 8

    # Return the number of gigabytes
    return Nbytes/1e9

# -----------------------------------------------------------------