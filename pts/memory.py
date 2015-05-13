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

## Return the estimated memory usage for a given ski file (in gigabytes)
def estimate_memory(skifile):

    # Get the number of wavelengths (force all calculations to floating point, avoiding integer overrun)
    Nlambda = float(skifile.nwavelengths())

    # Get the number of dust cells
    Ncells = skifile.ncells()

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
    #Ndoubles = 50e6 + (Nlambda + Ncells + Ncomps + Npops) * 10
    #print 1, Ndoubles*8/1e9

    Ndoubles = 0

    # Instruments
    for instrument in skifile.npixels():

        Ndoubles += instrument[2]

    # Dust system
    Ndoubles += (Ncomps+1) * Ncells
    if dustEmission:
        Ndoubles += Nlambda * Ncells
        if selfAbsorption:
            Ndoubles += Nlambda * Ncells

    # Dust grid
    #if treeGrid:
    #    Ndoubles += Ncells * 1.2 * 50
    #print 4, Ndoubles*8/1e9

    # Dust mixes
    Ndoubles += Nlambda * Npops * 3
    Ndoubles += (Nlambda + Npops) * 10

    # dust library
    if dustEmission:
        Ndoubles += Ncells + Nlambda*Nitems

    # transient heating
    if dustEmission and transientHeating:
        NT = 2250.
        Ndoubles += Ncomps * (Nlambda+1) * NT
        Ndoubles += Npops * 5./8.*NT*NT

    # 8 bytes in a double
    Nbytes = Ndoubles * 8

    #print "Estimated memory use: {:.1f} GB".format(Nbytes/1e9)

    return Nbytes/1e9

# -----------------------------------------------------------------