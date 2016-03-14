#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('text', usetex=True)

# -----------------------------------------------------------------

def main():

    #outpath  = "modelChecks/"
    #inpath   = "SKIRTOutput/iteration5_J14/"
    #refSED   = "Files/M31skirtSED.dat"
    plotting = True

    # Constants
    parsec = 3.08567758e16  # in meters
    Lsun = 3.846e26  # in Watt
    
    # Distance of M81 in Mpc
    D = 3.6 # 3.6 Mpc
    
    D *= 1e6 * parsec

    # scan for sed files
    #sedfiles = []
    #os.chdir(inpath)
    #for file in glob.glob("M31*full_i0_sed.dat"):
    #    sedfiles.append(file)
    #print "Found "+str(len(sedfiles))+" sed files."
    #os.chdir('../../')

    print "Reading in data..."

    # Load the observed SED
    observed_sed_path = "SED_simple.dat"
    input   = np.loadtxt(observed_sed_path)
    
    # Get the wavelength (in micron), the observed flux (in Jansky) and the error (also in Jansky)
    obsWls  = input[:,0]  # in micron
    obsFlux = input[:,1]
    obsErr  = input[:,2]

    # Convert to solar luminosities
    obsFlux = obsFlux * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun
    obsErr  = obsErr * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun

    # Load the simulated SED
    skirt_sed_path = "M81_m1_1e4_i59_sed.dat"
    input      = np.loadtxt(skirt_sed_path)
    
    # Get the wavelength (in micron) and the simulated flux
    modWls     = input[:,0]
    modFlux    = input[:,1]
    
    #modDirect  = input[:,2]
    #modStellarScatter = input[:,3]
    #modDust    = input[:,4]
    #modDustScatter = input[:,5]
    #modTrans   = input[:,6]
    
    # Convert to solar luminosities
    #modFlux    = modFlux    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    modFlux    = modFlux   * 4.*np.pi*D**2 / Lsun  # W/m2 to Lsun
    
    #modDirect  = modDirect  * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    #modStellarScatter = modStellarScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    #modDust    = modDust    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    #modTrans   = modTrans   * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    #modDustScatter = modDustScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
    
    # Create the plot
    plt.figure(figsize=(10,6))
    plt.ylabel('Luminosity/L$_\odot$',fontsize=24)
    plt.xlabel('$\lambda/\mu\mathrm{m}$',fontsize=24)
    #plt.xlim(0.1,1.e3)
    #plt.ylim(1e7,1.e11)
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(labelsize=20)
    #plt.subplots_adjust(bottom=0.1)
    
    #plt.plot(modWls,modDirect, 'c-', label="Direct")
    #plt.plot(modWls,modStellarScatter, 'm-', label="Scatter")
    #plt.plot(modWls,modDust, 'r-', label="Dust")
    #plt.plot(modWls,modTrans, 'b-', label="Transparent")
    
    # Plot the simulated and observed fluxes
    plt.plot(modWls, modFlux, 'k-', linewidth=2, label="Total")
    plt.errorbar(obsWls, obsFlux, yerr=obsErr, color='g', fmt='ko', markersize=8, label="Observed")
    
    plt.tight_layout()
    plt.legend(loc='upper right',numpoints=1,markerscale=1.5,fontsize=12)
    
    output_path = "sed.pdf"
    
    plt.savefig(output_path, format='pdf')

    #plt.show()
    plt.close()

# -----------------------------------------------------------------
