#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.starformation Post processing data from EAGLE simulation results.
#
# This module incorporates post processing functions to handle EAGLE output, for use in PTS.
#
# James Trayford - 26/01/2015

# -----------------------------------------------------------------

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

## Reads a hdf5 field's attributes into a python dictionary.
def readAttrs(hdf, field='Header'):
     fieldobj = hdf[field]
     headkeys = list(fieldobj.attrs)
     values = []
     for i in headkeys:
         values.append(fieldobj.attrs[str(i)])
     return dict(zip(headkeys, values))

# -----------------------------------------------------------------

## Reads the Schmidt parameters into a python structure.
def schmidtParameters(consts, runpars):

    # extract relevent unit conversions
    CM_PER_MPC   = consts['CM_PER_MPC']
    GAMMA        = consts['GAMMA']
    GRAV         = consts['GRAVITY']
    K_B          = consts['BOLTZMANN']
    M_PROTON     = consts['PROTONMASS']
    M_SUN        = consts['SOLAR_MASS']
    SEC_PER_YEAR = consts['SEC_PER_YEAR']

     # extract relevent runtime parameters used to create EAGLE snapshot
    GammaEff     = runpars['EOS_Jeans_GammaEffective']
    InitH        = runpars['InitAbundance_Hydrogen']
    RhoHi        = runpars['SF_SchmidtLawHighDensThresh_HpCM3']
    RhoNorm      = runpars['EOS_NormPhysDens_HpCM3']
    SchmidtCoeff = runpars['SF_SchmidtLawCoeff_MSUNpYRpKPC2']
    SchmidtExp   = runpars['SF_SchmidtLawExponent']
    SchmidtExpHi = runpars['SF_SchmidtLawHighDensExponent']
    T_JeansNorm  = runpars['EOS_Jeans_TempNorm_K']

    # Normalisation in cgs units
    Norm_cgs     = SchmidtCoeff * pow(pow(CM_PER_MPC / 1.e6, 2) / M_SUN , SchmidtExp - 1) / (1.e6 * SEC_PER_YEAR)

    # High density Threshold
    RhoHi_cgs    = RhoHi * M_PROTON / InitH

    # Density normalisation in cgs
    RhoNorm_cgs  = RhoNorm * M_PROTON / InitH

    # Min total Pressure
    P_totc       = RhoNorm * T_JeansNorm * K_B / (InitH * 1.22)

    # Pressure at high density Schmidt law break
    PBreak_cgs   = P_totc * (RhoHi/RhoNorm) ** GammaEff

    # Assume f_g = 1
    NormHi_cgs   = Norm_cgs * (GAMMA * PBreak_cgs / GRAV) ** ((SchmidtExp - SchmidtExpHi) * 0.5)

    # tuple of universal SF parameters
    sfparams     = RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff

    # tuples of high and low pressure SF parameters
    sf_lo        = Norm_cgs, GAMMA/GRAV, SchmidtExp
    sf_hi        = NormHi_cgs, GAMMA/GRAV, SchmidtExpHi

    return sfparams, sf_lo, sf_hi

# -----------------------------------------------------------------

## Function to obtain SFR of gas from which star particles formed.
#
# Inputs:
#  - rho_form: gas density at formation of star particle
#  - mass: mass of star particle
#  - schmidtpars: parameters for implementing Schmidt law from schmidtParameters()
#
# Outputs:
#  - SFR = Star formation rate for gas particle in input mass units per year
#
def getSFR(rho_form, mass, schmidtpars):

    # unpack universal SF law parameters
    RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff = schmidtpars[0]

    # Pressure at star formation
    P_form = P_totc * (rho_form / RhoNorm_cgs) ** GammaEff

    # unpack high and low pressure SF law parameters
    sf_lo, sf_hi = schmidtpars[1:]

    # calculate SFR
    if type(rho_form) == np.ndarray:
        hidx = rho_form > RhoHi_cgs
        SFR  = np.zeros(rho_form.size)
        if np.any(hidx):
            SFR[hidx] = mass[hidx] * sf_hi[0] * (sf_hi[1] * P_form[hidx]) ** ((sf_hi[2] - 1) * 0.5)
        if np.any(-hidx):
            SFR[-hidx] = mass[-hidx] * sf_lo[0] * (sf_lo[1] * P_form[-hidx]) ** ((sf_lo[2] - 1) * 0.5)
    else:
        if rho_form > RhoHi_cgs:
            SFR = mass * sf_hi[0] * (sf_hi[1] * P_form) ** ((sf_hi[2] - 1) * 0.5)
        else:
            SFR = mass * sf_lo[0] * (sf_lo[1] * P_form) ** ((sf_lo[2] - 1) * 0.5)

     # return SFR converted to input mass units per year from per second
    return np.array(SFR) * 3.15569e7

# -----------------------------------------------------------------

## Function to obtain ambient pressure of gas from which star particles formed.
#
# Inputs:
#  - rho: gas density of star forming particle
#  - schmidtpars: parameters for implementing Schmidt law from schmidtParameters()
#
# Outputs:
#  - P_tot: Ambient pressure from polytropic effective EoS (Schaye & Dalla Vecchia (2004))
#
def getPtot(rho, schmidtpars):
    RhoNorm_cgs, RhoHi_cgs, P_totc, PBreak_cgs, GammaEff = schmidtpars[0]
    P_form = P_totc * (rho / RhoNorm_cgs) ** GammaEff
    return P_form

# -----------------------------------------------------------------

## Function to sample a star forming gas particles into a number of sub-particles.
#
# Inputs:
#  - sfr: star formation rate in solar masses per yr
#  - m_gas: particle mass in solar masses
#
# Outputs:
#  - nested arrays with a list of subparticles for each parent input particle:
#     - ms: sub-particle stellar masses in solar masses
#     - ts: lookback times of sub-particle formation
#     - idxs: index of the sub-particle's parent particle in input array
#  - mdiffs: mass of parent particles locked up in new stars; this can be subtracted from the parent gas
#            particles for mass conservation
#
def stochResamp(sfr, m_gas):

    # mass resampling parameters (see Kennicutt & Evans 2012 section 2.5)
    m_min = 700         # minimum mass of sub-particle in M_solar
    m_max = 1e6         # maximum mass of sub-particle in M_solar
    alpha = 1.8         # exponent of power-law mass function
    alpha1 = 1. - alpha

    # age resampling parameters
    thresh_age = 1e8    # period over which to resample in yr (100 Myr)

    # initialise lists for output
    ms   = [[]]
    ts   = [[]]
    idxs = [[]]
    mdiffs = []

    # for each parent particle, determine the star-forming sub-particles
    for i in range(sfr.size):
        sfri = sfr[i]
        mi = m_gas[i]

        # determine the maximum number of sub-particles based on the minimum sub-particle mass
        N = int(max(1,np.ceil(mi/m_min)))

        # generate random sub-particle masses from a power-law distribution between min and max values
        X = np.random.random(N)
        m = (m_min**alpha1 + X*(m_max**alpha1-m_min**alpha1))**(1./alpha1)

        # limit and normalize the list of sub-particles to the total mass of the parent
        mlim = m[np.cumsum(m)<=mi]
        if len(mlim)<1: mlim = m[:1]
        m = mi/mlim.sum() * mlim
        N = len(m)

        # generate random decay lookback time for each sub-particle
        X = np.random.random(N)               # X in range (0,1]
        t = thresh_age + mi/sfri * np.log(1-X)

        # determine mask for sub-particles that form stars by present day
        issf = t > 0.

        # add star-forming sub-particles to the output lists
        ms.append(m[issf])
        ts.append(t[issf])
        idxs.append([i]*np.count_nonzero(issf))
        mdiffs.append(m[issf].sum())

    # convert sub-particle lists into numpy arrays
    ms     = np.hstack(ms)
    ts     = np.hstack(ts)
    idxs   = np.hstack(idxs).astype(int)
    mdiffs = np.array(mdiffs)

    return ms, ts, idxs, mdiffs

# -----------------------------------------------------------------

## Function to randomly shift the positions of HII region sub-particles within the smoothing sphere of their parent
#
# Arguments:
#  - r: parent positions; updated by this function to the shifted positions
#  - h: the smoothing lengths of the parents
#  - h_mapp: the smoothing lengths of the sub-particles
#
def stochShiftPos(r, h, h_mapp):
    # the offset sampling smoothing length is determined so that in the limit of infinite particles,
    # the light distribution is the same as the parent particle kernel;
    # assuming Gaussian kernels this means h_sampling**2 + h_mapp**2 = h**2.
    h_sampling = np.sqrt(np.maximum(0,h*h - h_mapp*h_mapp))

    # sample the offset from a scaled gaussian that resembles a cubic spline kernel
    # (see the documentation of the SPHDustDistribution class in SKIRT)
    r[:,0] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)
    r[:,1] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)
    r[:,2] += h_sampling * np.random.normal(scale=0.29, size=h_sampling.shape)

# -----------------------------------------------------------------
