#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_seds

# -----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import glob
import os
from scipy import interpolate
from scipy import integrate

def main():

    outpath  = "modelChecks/"
    inpath   = "SKIRTOutput/"
    refSED   = "Files/M31skirtSED_weight.dat"
    plotting = True
    OutputChi2List = False
    # Bands in refSED for which use=1
    bands  = ['FUV','NUV', 'u','g','r','i','z', 'W1','3.6','4.5','W2', '5.8','8',
              'W3','W4', '24','70', '100','160', '250','350','500']

    #bands  = ['FUV', '3.6', '24', '100','160', '250','350','500']

    D = 2.4222569e22 # 0.785 Mpc in m
    Lsun = 3.846e26
    

    # scan for sed files
    sedfiles = []
    os.chdir(inpath)
    for file in glob.glob("M31_*_i77.5_sed.dat"):
        sedfiles.append(file)
    print "Found "+str(len(sedfiles))+" sed files."
    os.chdir('../')

    print "Reading in data..."

    # Observed SED
    input   = np.loadtxt(refSED)
    obsWls  = input[:,0]
    obsFlux = input[:,1]
    obsErr  = input[:,2]
    chi2weight = input[:,3]

    obsFlux = obsFlux * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun
    obsErr  = obsErr * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(obsWls*1.e-6) / Lsun

    compareFlux = obsFlux[chi2weight>=0]
    compareErr  = obsErr[chi2weight>=0]
    

    # write chi2 values to file
    if OutputChi2List:
        file = open(outpath+"chi2list_weighted.dat",'w')
        file.write('#model    chi2\n')

    chi2list = []
    sedlist  = []
    
    for sed in sedfiles:
        print 'Processing '+sed
        # skirt SED
        input      = np.loadtxt(inpath+sed)
        modWls     = input[:,0]
        modFlux    = input[:,1]
        modDirect  = input[:,2]
        modStellarScatter = input[:,3]
        modDust    = input[:,4]
        modDustScatter = input[:,5]
        modTrans   = input[:,6]
        
        modFlux    = modFlux    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDirect  = modDirect  * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modStellarScatter = modStellarScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDust    = modDust    * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modTrans   = modTrans   * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun
        modDustScatter = modDustScatter * 1.e-26 * 4.*np.pi*D**2 * 3.e8/(modWls*1.e-6) / Lsun

        modBands = np.array([])
        for band in bands:
            filter = "Files/Filters/transmission_"+band+".dat"
            modBands = np.append(modBands,convolveFilter(modFlux,modWls,filter))

        chi2 = np.sum(chi2weight[chi2weight>=0]*(compareFlux - modBands)**2/compareErr**2)
        chi2list.append(chi2)
        sedlist.append(sed.replace('M31_','').replace('_i77.5_sed.dat',''))
    
        ####################################################
        if plotting:
            # SED
            plt.figure(figsize=(10,6))
            plt.ylabel('Luminosity/L$_\odot$',fontsize=24)
            plt.xlabel('$\lambda/\mu\mathrm{m}$',fontsize=24)
            plt.xlim(0.1,1.e3)
            plt.ylim(1e7,1.e11)
            #plt.ylim(0.1,1.e5)
            plt.xscale('log')
            plt.yscale('log')
            plt.tick_params(labelsize=20)
            #plt.subplots_adjust(bottom=0.1)
            plt.plot(modWls,modDirect, 'c-', label="Direct")
            plt.plot(modWls,modStellarScatter, 'm-', label="Scatter")
            plt.plot(modWls,modDust, 'r-', label="Dust")
            plt.plot(modWls,modTrans, 'b-', label="Transparent")
            plt.plot(modWls,modFlux, 'k-', linewidth=2, label="Total")
            plt.errorbar(obsWls,obsFlux,yerr=obsErr,color='g', fmt='ko',markersize=8, label="Observed")
            plt.tight_layout()
            plt.legend(loc='upper right',numpoints=1,markerscale=1.5,fontsize=12)
            sed = sed.replace('dat','pdf')
            plt.savefig(outpath+sed, format='pdf')

            #plt.show()
            plt.close()

    # Sort according to model number
    if OutputChi2List:
        sedlist = np.array(sedlist, dtype=int)
        
        idx = np.argsort(sedlist)
        rankedSEDList = sedlist[idx]
        rankedChi2List = np.array(chi2list)[idx]
    # Sort according to chi2 value
    idx = np.argsort(chi2list)

    sortedSEDs = []
    sortedChi2 = []
    for i in range(0,len(chi2list)):
        sortedSEDs.append(sedfiles[idx[i]])
        sortedChi2.append(chi2list[idx[i]])
        if OutputChi2List:
            file.write(str(rankedSEDList[i])+'    '+str(rankedChi2List[i])+'\n')

    if OutputChi2List:
        file.close()


    print 'Best fit chi2 values: '
    i=0
    while i < len(sortedSEDs) and i < 5:
        print 'Chi2 = '+str(sortedChi2[i])+', SED file = '+sortedSEDs[i]
        i += 1

def convolveFilter(flux,wl,filter):
    input = np.loadtxt(filter)
    response_wl     = input[:,0] * 1.e-4
    response_trans  = input[:,1]
    
    minwl = response_wl[0]
    maxwl = response_wl[len(response_wl)-1]
    intwl = np.copy(wl)
    intwl[intwl > maxwl] = -1
    intwl[intwl < minwl] = -1
    wlrange = intwl > 0
    intwl = intwl[wlrange]
    transmission = np.zeros(len(flux))
    
    interpfunc = interpolate.interp1d(response_wl,response_trans, kind='linear')
    transmission[wlrange] = interpfunc(intwl)
    tot_trans = integrate.simps(transmission,wl)
    tot_flux  = integrate.simps(transmission*flux,wl)


    return tot_flux/tot_trans

    #plt.plot(wl,transmission,'bo')
    #plt.plot(wl,transmission,'k-')
    #plt.plot(response_wl,response_trans,'r+')
    #plt.xscale('log')
    #plt.show()


if __name__ == '__main__':
    main()