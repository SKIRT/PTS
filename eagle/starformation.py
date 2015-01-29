#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.starformation Post processing data from EAGLE simulation results.
#
# This module incorporates post processing functions to handle EAGLE output, for use in PTS.
#
# James Trayford - 26/01/2015

# -----------------------------------------------------------------

from __future__ import division
import numpy as np

# -----------------------------------------------------------------

## Reads a hdf5 field's attributes into a python dictionary
def readAttrs(hdf, field='Header'):
     fieldobj = hdf[field]
     headkeys = list(fieldobj.attrs)
     values = []
     for i in headkeys:
         values.append(fieldobj.attrs[str(i)])
     return dict(zip(headkeys, values))

# -----------------------------------------------------------------

## Reads the Schmidt parameters into a python structure
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

# Function to sample a star forming gas particles into a number of sub particles at a set mass resolution.
#
# Inputs:
#  - sfr: star formation rate in solar masses per yr
#  - m_gas: particle mass in solar masses
#  - del_m: mass dicretisation in solar masses
#  - thresh_age: period over which to resample in yr (100 Myr by default)
#
# Outputs:
#  - ms: sub-particle stellar masses in solar masses
#  - ts: lookback times of sub-particle formation
#  - sfrs: the birth star formation rate of sub-particles
#  - idxs: index array to index inherited sub-particle properties from particle arrays
#  - mdiffs: mass of parent particles locked up in new stars; this can be subtracted from the parent gas
#            particles for mass conservation
#
def stochResamp(sfr, m_gas, delm=1.e4, thresh_age = 1e8):

    # Number of sub-particles for resampling at resolution of delm
    Ns = np.round(m_gas / delm).astype(int)

    # find masses and SF time scales of sub-particles
    sspm = m_gas / Ns
    taus = m_gas / sfr

    # initialise lists of lists
    ms   = [[]]
    ts   = [[]]
    Zs   = [[]]
    sfrs = [[]]
    idxs = [[]]

    # calculate mass converted to stars in each particle
    mdiffs = []
    for i in range(Ns.size):
        X = np.random.random(Ns[i])           # generate Ns[i] uniform random variables in range (0,1]
        t = thresh_age + taus[i]*np.log(1-X)  # calculate decay lookback time

        isCnv = t > 0.                        # check if star forms by present day
        ages  = t[isCnv]                      # remove stars forming in the future
        N_cnv = isCnv.sum()                   # calculate number converted

        ts.append(ages)
        ms.append([sspm[i]]*N_cnv)
        sfrs.append([sfr[i]/Ns[i]]*N_cnv)
        idxs.append([i]*N_cnv)
        mdiffs.append(sspm[i]*N_cnv)

    # convert sub-particle lists into numpy arrays
    ms     = np.hstack(ms)
    ts     = np.hstack(ts)
    sfrs   = np.hstack(sfrs)
    idxs   = np.hstack(idxs).astype(int)
    mdiffs = np.array(mdiffs)

    return ms, ts, sfrs, idxs, mdiffs

# -----------------------------------------------------------------
