#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.mappings Contains the Mappings class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np

# Import the relevant PTS classes and modules
from ...core.data.sed import SED
from ...core.tools import introspection, tables
from ...core.tools import filesystem as fs
from ...core.units.parsing import parse_unit as u
from pts.core.tools.utils import lazyproperty
from ...core.basics.log import log

# -----------------------------------------------------------------

# Determine the path to the Bruzual-Charlot SED directory
bc_path = fs.join(introspection.skirt_repo_dir, "dat", "SED", "BruzualCharlot")
chabrier_path = fs.join(bc_path, "chabrier") # chabrier IMF

# -----------------------------------------------------------------

class BruzualCharlot(object):
    
    """
    This class ...
    """

    def __init__(self, metallicity, age):

        """
        The constructor ...
        :param metallicity:
        :param age:
        """

        # Set properties
        self.metallicity = metallicity
        self.age = age

    # -----------------------------------------------------------------

    @lazyproperty
    def sed(self):

        """
        This function ...
        :return:
        """

        # Create the intrinsic SED and return it
        return create_bruzual_charlot_sed(self.metallicity, self.age)

    # -----------------------------------------------------------------

    @classmethod
    def dust_mass_for_luminosity(cls, metallicity, compactness, pressure, covering_factor, luminosity, wavelength):

        """
        This function ...
        :param metallicity:
        :param compactness:
        :param pressure:
        :param covering_factor:
        :param luminosity:
        :param wavelength:
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def luminosity_at(self, wavelength, unit="W/micron"):

        """
        This function ...
        :param wavelength:
        :param unit:
        :return:
        """

        return self.sed.photometry_at(wavelength).to(unit)

    # -----------------------------------------------------------------

    def luminosity_for_filter(self, fltr, unit="W/micron"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        luminosity = fltr.convolve(self.sed.wavelengths(unit="micron", asarray=True), self.sed.photometry(unit="W/micron", asarray=True)) # also in W/micron
        luminosity = luminosity * u("W/micron")

        # Return the luminosity in the desired unit
        return luminosity.to(unit)

# -----------------------------------------------------------------

def create_bruzual_charlot_sed(metallicity, age, mass=1.):

    """
    This function ...
    :param metallicity:
    :param age:
    :param mass:
    :return:
    """

    #seconds_in_year = 31557600.

    #// find the appropriate SED from interpolating in the BC library

    #// get the luminosities for arbitrary mass and for the appropriate Z and t
    #setLuminosities(bc.luminosities(1., _metallicity, _age / Constants::year()));

    # double M, double Z, double t, double z=0

    # Load the data
    _lambdav, _tv, _zv, j_dict = load_bruzual_charlot_data()

    nlambda = 1221
    nz = 6
    nt = 221
    hZ = 0.0

    if metallicity <= _zv[0]: mL = mR = 0
    elif metallicity >= _zv[nz-1]: mL = mR = nz-1
    else:

        mL = locate_clip(_zv,metallicity)
        mR = mL+1
        #double ZL = _Zv[mL];
        #double ZR = _Zv[mR];
        zl = _zv[mL]
        zr = _zv[mR]
        #hZ = (Z-ZL)/(ZR-ZL);
        hz = (metallicity-zl)/(zr-zl)

    ht = 0.0

    age_years = age.to("yr").value

    #print(age, age_years)

    #tl = tr = None
    #print(_tv[0])
    #print(_tv[-1])

    if age_years <= _tv[0]: pL = pR = 0
    elif age_years >= _tv[nt-1]: pL = pR = nt - 1
    else:

        #pL = NR::locate_clip(_tv,t);
        pL = locate_clip(_tv, age_years)
        pR = pL + 1
        #double tL = _tv[pL];
        #double tR = _tv[pR];
        tl = _tv[pL]
        tr = _tv[pR]
        #ht = (t-tL)/(tR-tL);
        ht = (age_years - tl) / (tr - tl)

    #print(tl, tr, pL, pR, ht)

    #const Array& jLLv = _jvv(pL,mL);
    #const Array& jLRv = _jvv(pL,mR);
    #const Array& jRLv = _jvv(pR,mL);
    #const Array& jRRv = _jvv(pR,mR);

    #Array jv(Nlambda);

    jv = [None] * nlambda

    jLLv = j_dict[(pL, mL)]
    jLRv = j_dict[(pL, mR)]
    jRLv = j_dict[(pR, mL)]
    jRRv = j_dict[(pR, mR)]

    #for (int k=0; k<Nlambda; k++)
    for k in range(nlambda):

        jv[k] = (1.0-ht)*(1.0-hZ) * jLLv[k] + (1.0-ht)*hZ*jLRv[k] + ht*(1.0-hZ)*jRLv[k] + ht * hZ * jRRv[k]

    #// resample to the possibly redshifted simulation wavelength grid,
    #// convert emissivities to luminosities (i.e. multiply by the wavelength bins),
    #// multiply by the mass of the population (in solar masses),
    #// and return the result

    #return NR::resample<NR::interpolate_loglog>(_lambdagrid->lambdav()*(1-z), _lambdav, jv)
    #                                * _lambdagrid->dlambdav() * M;

    # Return
    #return jv

    #return NR::resample < NR::interpolateLogLog > (_lambdagrid->lambdav() * (1 - z), _lambdav, jv) * _lambdagrid->dlambdav() * M

    # Array dlambdav(Nlambda);
    #setWavelengthBins(lambdav, dlambdav);   // set the grid with temporary widths so that lambdamin/max work
    #for (int ell=0; ell<Nlambda; ell++)
    #{
    #    dlambdav[ell] = lambdamax(ell)-lambdamin(ell);
    #}

    # In micrometers
    wavelength_column = [wavelength_meter * 1e6 for wavelength_meter in _lambdav]

    dlambdas = []
    for index in range(nlambda):
        dlambda = lambdamax(wavelength_column, index) - lambdamin(wavelength_column, index)
        dlambdas.append(dlambda)

    if hasattr(mass, "unit"): mass_in_solar = mass.to("Msun").value
    else: mass_in_solar = mass

    # Create luminosities
    luminosity_column = [None] * nlambda
    for index in range(nlambda):
        j = jv[index]
        dlambda = dlambdas[index]
        lum = j * dlambda * mass_in_solar
        luminosity_column[index] = lum

    # Create the SED
    sed = SED.from_arrays(wavelength_column, luminosity_column, wavelength_unit="micron", photometry_unit="W/micron")

    # Return the SED
    return sed

# -----------------------------------------------------------------

def lambdamin(wavelengths, ell):

    """
    This function ...
    :param wavelengths:
    :param ell:
    :return:
    """

    return wavelengths[0] if ell == 0 else math.sqrt(wavelengths[ell-1]*wavelengths[ell])

# -----------------------------------------------------------------

def lambdamax(wavelengths, ell):

    """
    This function ...
    :param wavelengths:
    :param ell:
    :return:
    """

    Nlambda = len(wavelengths)
    return wavelengths[-1] if ell == Nlambda - 1 else math.sqrt(wavelengths[ell]*wavelengths[ell+1])

# -----------------------------------------------------------------

def load_bruzual_charlot_data():

    """
    This function ...
    :return:
    """

    # number of items in the library
    nlambda = 1221
    nt = 221
    nz = 6

    # local constants for units
    Angstrom = 1e-10
    solar_luminosity = 3.846e26 * u("W")
    Lsun = solar_luminosity.to("W").value

    #// Prepare the vectors for the Bruzual & Charlot library SEDs
    #_lambdav.resize(Nlambda)
    #_tv.resize(Nt)
    #_Zv.resize(NZ)
    #_jvv.resize(Nt, NZ, Nlambda)

    _lambdav = [None] * nlambda
    _tv = [None] * nt
    #_zv = [None] * nz
    #_jvv = [[None] * ]

    _zv = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
    zcodev = ["m22", "m32", "m42", "m52", "m62", "m72"]

    j_dict = dict()

    # Read the wavelength, age and emissivity vectors from the Bruzual & Charlot library
    for m in range(nz):

        filepath = fs.join(chabrier_path, "bc2003_lr_" + zcodev[m] + "_chab_ssp.ised_ASCII")
        #print(filepath)

        # Inform the user
        log.info("Reading SED data from file " + filepath + "...")

        fh = fs.read_words(filepath, newlines=True)

        first = next(fh)
        nt_file = int(first)
        if nt_file != nt: raise IOError("iNt is not equal to Nt")

        p = 0
        #for p in range(nt):
        while p < nt:

            word = next(fh)
            if word == "\n": continue
            t = float(word)
            _tv[p] = t    # // age in file in yr, we want in yr
            p += 1

        #print(_tv)

        # skip six lines...
        nlines = 0
        while nlines < 6:
            word = next(fh)
            if word == "\n": nlines += 1

        line = next(fh)
        inlambda = int(line)
        if inlambda != nlambda: raise IOError("iNlambda is not equal to Nlambda")

        k = 0
        #for k in range(nlambda):
        while k < nlambda:
            #double lambda;
            #bcfile >> lambda;
            word = next(fh)
            #print(list(word))
            if word == "\n": continue
            lam = float(word)
            #if lam > 1: print(lam)
            _lambdav[k] = lam * Angstrom   #;   // lambda in file in A, we want in m
            #if _lambdav[k] > 1: print(lam)
            k += 1

        #from ...magic.tools import plotting
        #plotting.plot()

        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(range(nlambda), _lambdav)
        #ax = plt.gca()
        #ax.set_yscale("log")
        #plt.show()

        #print(_lambdav)

        p = 0
        while p < nt:
        #for p in range(nt):

            #Array& jv = _jvv(p,m);
            #bcfile >> iNlambda;
            #if (iNlambda != Nlambda)
            #    throw FATALERROR("iNlambda is not equal to Nlambda");

            word = next(fh)
            if word == "\n": word = next(fh)
            #print(word)
            inlambda = int(word)
            #print(inlambda, nlambda)
            if inlambda != nlambda: raise IOError("iNlambda is not equal to Nlambda")

            jv = [None] * nlambda

            k = 0
            #for k in range(nlambda):
            while k < nlambda:

                #double j;
                #bcfile >> j;
                word = next(fh)
                if word == "\n": continue
                j = float(word)
                jv[k] = j * Lsun / Angstrom  # // emissivity in file in Lsun/A, we want in W/m.
                k += 1

            # Set the data
            j_dict[(p,m)] = jv

            #print(jv)

            #int idummy;
            #bcfile >> idummy;
            word = next(fh)
            if word == "\n": word = next(fh)
            idummy = int(word)
            #print(idummy)
            k = 0
            #for k in range(idummy):
            while k < idummy:

                #double dummy;
                #bcfile >> dummy;
                dummy = next(fh)
                if dummy == "\n": continue
                k += 1

            p += 1

    # Return the data
    return _lambdav, _tv, _zv, j_dict

# -----------------------------------------------------------------

def find_nearest(array, value):

    """
    This function ...
    :param array:
    :param value:
    :return:
    """

    idx = (np.abs(array-value)).argmin()
    return array[idx]

# -----------------------------------------------------------------

def locate_clip(array, value):

    """
    This function ...
    :param array: actually a list
    :param value:
    :return:
    """

    #n = array.size
    n = len(array)
    if value < array[0]: return 0
    return locate_basic_impl(array, value, n-1)

# -----------------------------------------------------------------

def locate_basic_impl(xv, x, n):

    """
    This function ...
    :param xv:
    :param x:
    :param n:
    :return:
    """

    jl = -1
    ju = n

    while ju - jl > 1:

        jm = (ju + jl) >> 1
        if x < xv[jm]: ju = jm
        else: jl = jm

    return jl

# -----------------------------------------------------------------
