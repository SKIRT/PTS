#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.data.extinction Contains the ExtinctionCurve classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy import interpolate

# Import astronomical modules
from astropy.units import Unit, spectral
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools import tables, introspection, arrays
from ...core.tools import filesystem as fs
from ..basics.curve import Curve
from ...magic.tools.wavelengths import extinction_wavelength_range

# -----------------------------------------------------------------

class ExtinctionCurve(Curve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Units
        x_unit = "micron"

        # Names
        x_name = "Wavelength"
        y_name = "Extinction"

        # Descriptions
        x_description = "Wavelength"
        y_description = "Extinction"

        kwargs["x_unit"] = x_unit
        kwargs["y_unit"] = None
        kwargs["x_name"] = x_name
        kwargs["y_name"] = y_name
        kwargs["x_description"] = x_description
        kwargs["y_description"] = y_description

        # If data is passed
        if "wavelengths" in kwargs and "extinctions" in kwargs:

            wavelengths = kwargs.pop("wavelengths")
            extinctions = kwargs.pop("extinctions")

        else: wavelengths = extinctions = None

        # Call the constructor of the base class
        super(ExtinctionCurve, self).__init__(*args, **kwargs)

        # Add the data
        if wavelengths is not None:
            for index in range(len(wavelengths)): self.add_row([wavelengths[index], extinctions[index]])

    # -----------------------------------------------------------------

    @classmethod
    def from_range(cls, wavelength_range=extinction_wavelength_range, npoints=100, scale="logarithmic", rv=3.1):

        """
        This function ...
        :param wavelength_range:
        :param npoints:
        :param scale
        :param rv:
        :return:
        """

        # Generate
        if scale == "logarithmic": wavelengths = wavelength_range.log(npoints)
        elif scale == "linear": wavelengths = wavelength_range.linear(npoints)
        else: raise ValueError("Invalid scale: " + scale)

        # Create and return
        return cls(wavelengths=wavelengths, rv=rv)

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :param unit:
        :param asarray:
        :param add_unit:
        :return:
        """

        if asarray: return arrays.plain_array(self["Wavelength"], unit=unit, array_unit=self.column_unit("Wavelength"))
        else: return arrays.array_as_list(self["Wavelength"], unit=unit, add_unit=add_unit, array_unit=self.column_unit("Wavelength"))

    # -----------------------------------------------------------------

    def extinctions(self, asarray=False):

        """
        This function ...
        :param asarray:
        :return:
        """

        if asarray: return arrays.plain_array(self["Extinction"])
        else: return arrays.array_as_list(self["Extinction"])

    # -----------------------------------------------------------------

    def extinction_at(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        interpolated = interpolate.interp1d(self.wavelengths(unit="micron", asarray=True), self.extinctions(asarray=True), kind='linear')
        return interpolated(wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def normalize_at(self, wavelength, value=1.):

        """
        This function ...
        :param wavelength:
        :param value:
        :return:
        """

        extinction_wavelength = self.extinction_at(wavelength)
        self["Extinction"] /= extinction_wavelength * value

# -----------------------------------------------------------------

extinction_path = fs.join(introspection.pts_subproject_dir("core"), "data", "extinction_functions.pyx")

import pyximport

#USE_CYTHON = True
fname = extinction_path

core_data_directory_path = fs.join(introspection.pts_subproject_dir("core"), "data")

#extinction_module_name = "extinction"
extinction_module_name = "pts.core.data.extinction_functions"

extern_path = fs.join(introspection.pts_subproject_dir("core"), "data", "extern")
bs_c_path = fs.join(extern_path, "bs.c")
bs_h_path = fs.join(extern_path, "bs.h")
bsplines_path = fs.join(extern_path, "bsplines.pxi")

sourcefiles = [fname, bs_c_path]
dependsfiles = [bs_h_path, bsplines_path]
include_dirs = [np.get_include(), extern_path]
#extensions = [Extension(extinction_module_name, sourcefiles, include_dirs=include_dirs,
#                        depends=dependsfiles, extra_compile_args=['-std=c99'])]

pyximport.install(build_dir=core_data_directory_path, setup_args={"include_dirs":include_dirs}, reload_support=True, pyimport=True)
#pyximport.install(setup_args={"include_dirs":include_dirs}, reload_support=True, pyimport=True)
from pts.core.data import extinction_functions

# -----------------------------------------------------------------

class CardelliClaytonMathisExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, wavelengths, rv=3.1):

        """
        This function ...
        :param wavelengths:
        :param rv:
        """

        # Generate
        extinctions = extinction_functions.ccm89(wavelengths, 1.0, rv)

        # Call the constructor of the base class
        super(CardelliClaytonMathisExtinctionCurve, self).__init__(wavelengths=wavelengths, extinctions=extinctions)

# -----------------------------------------------------------------

class ODonnellExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, wavelengths, rv=3.1):

        """
        This function ...
        :param wavelengths:
        """

        # Generate
        extinctions = extinction_functions.odonnell94(wavelengths, 1.0, rv)

        # Call the constructor of the base class
        super(ODonnellExtinctionCurve, self).__init__(wavelengths=wavelengths, extinctions=extinctions)

# -----------------------------------------------------------------

class FitzpatrickExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, wavelengths, rv=3.1):

        """
        Thisf unction ...
        :param wavelengths:
        :param rv:
        """

        # Generate
        extinctions = extinction_functions.fitzpatrick99(wavelengths, 1.0, rv)

        # Call the constructor of the base class
        super(FitzpatrickExtinctionCurve, self).__init__(wavelengths=wavelengths, extinctions=extinctions)

# -----------------------------------------------------------------

class FitzpatrickMassaExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, wavelengths):

        """
        Thisf unction ...
        :param wavelengths:
        """

        # Generate
        extinctions = extinction_functions.fm07(wavelengths, 1.0)

        # Call the constructor of the base class
        super(FitzpatrickMassaExtinctionCurve, self).__init__(wavelengths=wavelengths, extinctions=extinctions)

# -----------------------------------------------------------------

class CalzettiExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, wavelengths, rv=3.1):

        """
        This function ...
        :param wavelengths:
        :param rv:
        """

        # Generate
        extinctions = extinction_functions.calzetti00(wavelengths, 1.0, rv)

        # Call the constructor of the base class
        super(CalzettiExtinctionCurve, self).__init__(wavelengths=wavelengths, extinctions=extinctions)

# -----------------------------------------------------------------

# class MilkyWayAttenuationCurve(AttenuationCurve):
#
#     """
#     This class ...
#     """
#
#     def __init__(self, *args, **kwargs):
#
#         """
#         This function ...
#         :param rv: the value of R(V). If None, the mean Milky Way attenuation curve is generated. (See Fitzpatrick & Massa, 2007)
#         """
#
#         # Determine the path to the data file
#         #path = fs.join(attenuation_data_path, "AttenuationLawMWRv5.dat")
#
#         # Load the attenuation data
#         #wavelengths_angstrom, alambda_av = np.loadtxt(path, unpack=True)
#
#         # Convert wavelengths into micron
#         #wavelengths = wavelengths_angstrom * 0.0001
#
#         # Attenuations
#         #attenuations = alambda_av
#
#         # Get data
#         wavelengths, attenuations = generate_milky_way_attenuations(0.1, 1000, 1000)
#
#         kwargs["wavelengths"] = wavelengths
#         kwargs["attenuations"] = attenuations
#
#         # Call the constructor of the base class
#         super(MilkyWayAttenuationCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

# class SMCAttenuationCurve(AttenuationCurve):
#
#     """
#     This class ...
#     """
#
#     # Determine the path to the data file
#     dat_path = fs.join(attenuation_data_path, "AttenuationLawSMC.dat")
#
#     # -----------------------------------------------------------------
#
#     def __init__(self, *args, **kwargs):
#
#         """
#         This function ...
#         :param args:
#         :param kwargs:
#         """
#
#         # Load the attenuation data
#         wavelengths_angstrom, alambda_av = np.loadtxt(self.dat_path, unpack=True)
#
#         # Convert wavelengths into micron
#         wavelengths = wavelengths_angstrom * 0.0001
#
#         # Attenuations
#         attenuations = alambda_av
#
#         kwargs["wavelengths"] = wavelengths
#         kwargs["attenuations"] = attenuations
#
#         # Call the constructor of the base class
#         super(SMCAttenuationCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

# class CalzettiAttenuationCurve(AttenuationCurve):
#
#     """
#     This class ...
#     """
#
#     # Determine the path to the data file
#     path = fs.join(attenuation_data_path, "AttenuationLawCalzetti.dat")
#
#     # -----------------------------------------------------------------
#
#     def __init__(self, *args, **kwargs):
#
#         """
#         This function ...
#         :param args:
#         :param kwargs:
#         """
#
#         # Load the attenuation data
#         wavelengths_angstrom, alambda_av = np.loadtxt(self.path, unpack=True)
#
#         # Convert wavelengths into micron
#         wavelengths = wavelengths_angstrom * 0.0001
#
#         # Attenuations
#         attenuations = alambda_av
#
#         kwargs["wavelengths"] = wavelengths
#         kwargs["attenuations"] = attenuations
#
#         # Call the constructor of the base class
#         super(CalzettiAttenuationCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

class BattistiExtinctionCurve(ExtinctionCurve):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # From Battisit et al. 2016

        wl_B16 = np.arange(0.125, 0.832, 0.01)
        x = 1. / wl_B16
        Qfit_B16 = -2.488 + 1.803 * x - 0.261 * x ** 2 + 0.0145 * x ** 3
        interpfunc = interpolate.interp1d(wl_B16, Qfit_B16, kind='linear')
        Qfit_B16_V = interpfunc(0.55)  # Interpolate to find attenuation at V band central wavelengths
        n_atts_B16 = Qfit_B16 / Qfit_B16_V

        wavelengths = wl_B16
        attenuations = Qfit_B16 + 1. # TODO: understand this more ...

        kwargs["wavelengths"] = wavelengths
        kwargs["attenuations"] = attenuations

        # Call the constructor of the base class
        super(BattistiExtinctionCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

# class MappingsAttenuationCurve(AttenuationCurve):
#
#     """
#     This class ...
#     """
#
#     path = fs.join(attenuation_data_path, "AttenuationLawMAPPINGS.dat")
#
#     # -----------------------------------------------------------------
#
#     def __init__(self, *args, **kwargs):
#
#         """
#         This function ...
#         :param args:
#         :param kwargs:
#         """
#
#         attenuation = kwargs.pop("attenuation")
#         wavelength = kwargs.pop("wavelength")
#
#         # Load the data
#         # wl in micron from long to short wl.
#         # ABS attenuations (see header of data file)
#         wavelengths, abs_attenuations = np.loadtxt(self.path, unpack=True)
#
#         # CREATE A TABLE SO WE CAN EASILY SORT THE COLUMNS FOR INCREASING WAVELENGTH
#         names = ["Wavelength", "ABS attenuation"]
#         # Create the table
#         abs_table = tables.new([wavelengths, abs_attenuations], names)
#         abs_table["Wavelength"].unit = Unit("micron")
#         # Sort the table on wavelength
#         abs_table.sort("Wavelength")
#
#         wavelengths = np.array(list(abs_table["Wavelength"]))
#         abs_attenuations = np.array(list(abs_table["ABS attenuation"]))
#
#         # Find the ABS attenuation at the specified wavelength
#         interpolated = interpolate.interp1d(wavelengths, abs_attenuations, kind='linear')
#         abs_wavelength = interpolated(wavelength.to("micron").value)
#
#         # 'Real' attenuations
#         attenuations_mappings = abs_attenuations / abs_wavelength * attenuation
#
#         kwargs["wavelengths"] = wavelengths
#         kwargs["attenuations"] = attenuations_mappings
#
#         # Call the constructor of the base class
#         super(MappingsAttenuationCurve, self).__init__(*args, **kwargs)

# -----------------------------------------------------------------

def generate_milky_way_attenuations(wavelength_min, wavelength_max, Nsamp):

    """
    This function ...
    :param wavelength_min:
    :param wavelength_max:
    :param Nsamp:
    :return:
    """

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

    wl_UV = np.logspace(np.log10(wavelength_min),np.log10(0.2700), Nsamp/2) # UV part stops at 0.27 micron
    wl_ir = np.logspace(np.log10(0.2700),np.log10(wavelength_max), Nsamp/2) # optical-IR part starts at 0.27 micron
    idx = (np.abs(wl_ir-0.550)).argmin() # index closest to V band = 0.55 micron
    idx_U2 = (np.abs(wl_UV-0.2700)).argmin() # index closest to U2 band = 0.27 micron
    idx_U1 = (np.abs(wl_UV-0.2600)).argmin() # index closest to U1 band = 0.26 micron

    # construct UV attenuation curve
    x = 1./wl_UV
    D = Lorentzian(x, x0, gamma)
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

# -----------------------------------------------------------------

def Lorentzian(x, x0, gamma):

    """
    This function ...
    :param x:
    :param x0:
    :param gamma:
    :return:
    """

    return x*x / ((x*x-x0*x0)**2 + (x*gamma)**2)

# -----------------------------------------------------------------
