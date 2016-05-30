#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_attenuation

# -----------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import glob
import os
from scipy import interpolate

def main():

    outpath         = "modelChecks/iteration5_J14/"
    inpath          = "SKIRTOutput/iteration5_J14/"
    refAtt_SMC      = "Files/AttenuationLawSMC.dat"
    refAtt_MWRv5    = "Files/AttenuationLawMWRv5.dat"
    refAtt_Calzetti = "Files/AttenuationLawCalzetti.dat"
    refAtt_MAPPINGS = "Files/AttenuationLawMAPPINGS.dat"
    refSED_MAPPINGS = "SKIRTrun/models/testHeating/MappingsHeating/mappings_i0_sed.dat"

    # scan for sed files
    sedfiles = []
    os.chdir(inpath)
    for file in glob.glob("M31_212full_i77.5_sed.dat"):
        sedfiles.append(file)
    print "Found "+str(len(sedfiles))+" sed files."
    os.chdir('../../')

    print "Reading in data..."

    # Observed SED
    #input   = np.loadtxt(refAtts)
    #obsWls     = input[:,0]
    #obsFlux = input[:,1]
    #obsErr  = input[:,2]
    #useChi2  = input[:,3]

    input = np.loadtxt(refAtt_MAPPINGS)
    wl_mappings  = input[:,0] # wl in micron from long to short wl.
    abs_mappings = input[:,1] # factor
    # V band attenuation = total dust mass / total surface of M31 / 1e5 * 0.67 see data file.
    Av_mappings = 2143279.799/ (0.137**2*15230) / 1e5 * 0.67
    # Attenuation law (not normalized). abs for V band is 2.23403e-01
    A_mappings = abs_mappings / 2.23403e-01 * Av_mappings
    # Create interpolation function
    A_interpfunc = interpolate.interp1d(wl_mappings,A_mappings, kind='linear')
    
    # load the MAPPINGS SED
    input = np.loadtxt(refSED_MAPPINGS)
    flux_mappings = input[:,1] # in Jy


    input = np.loadtxt(refAtt_SMC)
    wl_SMC     = input[:,0]*0.0001 # wl in micron
    n_atts_SMC = input[:,1]

    input = np.loadtxt(refAtt_Calzetti)
    wl_Calzetti     = input[:,0]*0.0001 # wl in micron
    n_atts_Calzetti = input[:,1] #/ 1.0124025 # properly normalize function at 0.550 micron

    input = np.loadtxt(refAtt_MWRv5)
    wl_MWRv5     = input[:,0]*0.0001 # wl in micron
    n_atts_MWRv5 = input[:,1] # for Rv = 5.

    wl_MW, n_atts_MW = makeMilkyWayAtt(0.1,1000,1000)

    # From Dong et al 2014
    wl_M31 = np.array([1928, 2246, 2600, 2704, 3355, 3897, 4319, 4747, 5447 , 6656, 8057 , 11534, 15369])*1.e-4
    n_atts_M31  = [2.87, 3.92, 2.61, 2.80, 1.92, 1.58, 1.37, 1.21, 1.0,  0.73, 0.61, 0.46, 0.43]
    err_att_M31 = [ 0.38, 0.502, 0.091, 0.261, 0.029, 0.019, 0.018, 0.006, 0., 0.013, 0.031, 0.094, 0.099]

    # From Battisit et al. 2016
    wl_B16 = np.arange(0.125,0.832,0.01)
    x = 1./wl_B16
    Qfit_B16 = -2.488 + 1.803*x -0.261*x**2 + 0.0145*x**3
    interpfunc = interpolate.interp1d(wl_B16,Qfit_B16, kind='linear')
    Qfit_B16_V = interpfunc(0.55) # Interpolate to find attenuation at V band central wavelengths
    n_atts_B16 = Qfit_B16/Qfit_B16_V


    for sed in sedfiles:
        print "Processing "+sed
        # skirt SED
        input      = np.loadtxt(inpath+sed)
        modWls     = input[:,0]
        modFlux    = input[:,1]
        modDirect  = input[:,2]
        modStellarScatter = input[:,3]
        modDust    = input[:,4]
        modDustScatter = input[:,5]
        modTrans   = input[:,6]
        
        atts_mappings = A_interpfunc(modWls)    # Attenuation from dust in star forming regions
        delta_flux_mappings = flux_mappings * (10**(atts_mappings/2.5)-1) # additional flux from mappings attenuation

        atts_diff = -2.5*np.log10(modFlux/modTrans) # Attenuation from diffuse dust
        interpfunc = interpolate.interp1d(modWls,atts_diff, kind='linear')
        att_diff_V      = interpfunc(0.55) # Interpolate to find attenuation at V band central wavelengths
        n_atts_diff = atts_diff/att_diff_V

        atts_tot = -2.5*np.log10(modFlux/(modTrans+delta_flux_mappings))
        interpfunc = interpolate.interp1d(modWls,atts_tot, kind='linear')
        att_tot_V      = interpfunc(0.55) # Interpolate to find attenuation at V band central wavelengths
        n_atts_tot = atts_tot/att_tot_V

        interpfunc = interpolate.interp1d(modWls,atts_mappings, kind='linear')
        att_mappings_V      = interpfunc(0.55) # Interpolate to find attenuation at V band central wavelengths
        n_atts_mappings = atts_mappings/att_mappings_V

        f = open(outpath + sed.replace('_sed.dat','_att.dat'),'w')
        f.write('# wavelength_micron    diffuse_attenuation    Mappings_attenuation \n')
        for i in range(0,len(atts_diff)):
            f.write(str(modWls[i])+'    '+str(atts_tot[i])+'    '+str(atts_diff[i])+'    '+str(atts_mappings[i])+'\n')
        f.close()

        print ' V band attenuation from SF regions:   '+str(att_mappings_V)
        print ' V band attenuation from diffuse dust: '+str(att_diff_V)
        print ' V band attenuation in total:          '+str(att_tot_V)

        proportion = delta_flux_mappings / (modTrans-modFlux)
        print ' V band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls-0.55))])
        print ' NUV band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls-0.227))])
        print ' FUV band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls-0.153))])
        print ' Mean energy fraction of mappings attenuation to total:   ' + str(np.mean(proportion[modWls < 5 ]))
        print ' Median energy fraction of mappings attenuation to total: ' + str(np.median(proportion[modWls < 5 ]))

        #plt.plot(modWls, proportion)
        #plt.xscale('log')
        #plt.show()
                    

        ####################################################

        # Plot the attenuation curves
        plt.figure(figsize=(10,10))
        plt.ylabel('$A_\lambda/A_V$',fontsize=28)
        plt.xlabel('$\lambda/\mu$m',fontsize=28)
        plt.xlim(0.1,2)
        plt.ylim(0,8)
        plt.xscale('log')
        x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2]
        plt.xticks(x,x)
        #plt.yscale('log')
        plt.tick_params(labelsize=17)
        #plt.subplots_adjust(bottom=0.1)
        
        plt.plot(wl_MW,n_atts_MW, 'r-', label="MW",linewidth=2)
        plt.plot(wl_SMC,n_atts_SMC, 'b--', label="SMC",linewidth=2)
        plt.plot(wl_B16,Qfit_B16+1, 'g-.', label="Battisti+16",linewidth=3)
        #plt.plot(wl_Calzetti,n_atts_Calzetti, 'g-', label="Cal+94",linewidth=2)

        plt.plot(modWls,n_atts_tot, 'k-', label="M31 tot",linewidth=2)
        plt.plot(modWls,n_atts_diff, 'k--', label="M31 diffuse",linewidth=2)
        plt.plot(modWls,n_atts_mappings, 'k:', label="M31 SF regions",linewidth=2)
 
        plt.errorbar(wl_M31,n_atts_M31,yerr=err_att_M31,color='m', fmt='ko',markersize=8, label="M31- Dong+14")
        #plt.plot(1./wl_MWRv5,n_atts_MWRv5, 'm-', label="MW Rv=5")
        #plt.plot([1./0.55,1./0.55],[0,15], 'k-')

        plt.tight_layout()
        plt.legend(loc='upper right',numpoints=1,markerscale=1.5,fontsize=24)
        sed = sed.replace('_sed.dat','_att.pdf')
        plt.savefig(outpath+sed, format='pdf')

        #plt.show()
        plt.close()

def makeMilkyWayAtt2(minWl, maxWl,Nsamp):
    
    # Parameter values from Fitzpatrick & Massa 2007, table 5.
    x0 = 4.592
    gamma = 0.922
    c1 = -0.175
    c2 = 0.807
    c3 = 2.991
    c4 = 0.319
    c5 = 6.097
    k_ir = 1.057
    Rv = 3.001
    #Rv = 5.5
    
    wl = np.logspace(np.log10(minWl),np.log10(maxWl),Nsamp)
    idx_V = (np.abs(wl-0.550)).argmin() # index closest to V band = 0.55 micron

    x = 1./wl
    D = Lorentzian(x,x0,gamma)
    k = np.zeros(Nsamp)
    
    for i in range(0,len(x)):
        if x[i] <= c5:
            k[i] = c1 + c2*x[i] + c3*D[i]
        else:
            k[i] = c1 + c2*x[i] + c3*D[i] + c4*(x[i]-c5)**2

    Al_Av = k/Rv + 1.
    # print Al_Av[idx]
    return wl, Al_Av - (Al_Av[idx]-1.)
    #return wl, Al_Av

def makeMilkyWayAtt(minWl, maxWl,Nsamp):
    
    # Parameter values from Fitzpatrick & Massa 2007, table 5.
    x0 = 4.592
    gamma = 0.922
    c1 = -0.175
    c2 = 0.807
    c3 = 2.991
    c4 = 0.319
    c5 = 6.097
    O1 = 2.055
    O2 = 1.322
    O3 = 0.0
    k_ir = 1.057
    Rv = 3.001
    #Rv = 5.5
    
    wl_UV = np.logspace(np.log10(minWl),np.log10(0.2700),Nsamp/2) # UV part stops at 0.27 micron
    wl_ir = np.logspace(np.log10(0.2700),np.log10(maxWl),Nsamp/2) # optical-IR part starts at 0.27 micron
    idx = (np.abs(wl_ir-0.550)).argmin() # index closest to V band = 0.55 micron
    idx_U2 = (np.abs(wl_UV-0.2700)).argmin() # index closest to U2 band = 0.27 micron
    idx_U1 = (np.abs(wl_UV-0.2600)).argmin() # index closest to U1 band = 0.26 micron

    # construct UV attenuation curve
    x = 1./wl_UV
    D = Lorentzian(x,x0,gamma)
    k_UV = np.zeros(Nsamp/2)
    
    for i in range(0,len(x)):
        if x[i] <= c5:
            k_UV[i] = c1 + c2*x[i] + c3*D[i]
        else:
            k_UV[i] = c1 + c2*x[i] + c3*D[i] + c4*(x[i]-c5)**2

    # construct ir attenuation curve
    sample_wl = np.array([10000., 4., 2., 1.3333, 0.5530, 0.4000, 0.3300, 0.2700, 0.2600])
    sample_k = np.append( k_ir*sample_wl[0:4]**-1.84 - Rv, [O3,O2,O1, k_UV[idx_U2], k_UV[idx_U1]])

    spline = interpolate.UnivariateSpline(1./sample_wl,sample_k)
    k_ir = spline(1./wl_ir)

    wl    = np.append(wl_UV,wl_ir)
    Al_Av = np.append(k_UV,k_ir)/Rv + 1.
    return wl, Al_Av

def Lorentzian(x, x0, gamma):
    return x*x / ((x*x-x0*x0)**2 + (x*gamma)**2)

if __name__ == '__main__':
    main()