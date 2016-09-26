#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.filter Working with wavelength filters (or bands).
#
# An instance of the Filter class in this module represents a particular wavelength band, including its response or
# transmission curve, and allows integrating a given spectrum over the band. A filter instance can be constructed from
# one of the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import sys
import types
import numpy as np
from scipy.interpolate import interp1d
from lxml import etree

# -----------------------------------------------------------------

# Filtername , Frequency, filename
# Planck_350    857 GHz   857   HFI
# Planck_550    545 GHz   545   HFI
# Planck_850    352 GHz   353   HFI
# Planck_1380   222 GHz   217   HFI
# Planck_2100   142 GHz   143   HFI
# Planck_3000   100 GHz   100   HFI
# Planck_4260   70 GHz    070   LFI
# Planck_6810   44 GHz    044   LFI
# Planck_10600  28 GHz    030   LFI

planck_info = dict()
planck_info["Planck_350"] = ("857", "HFI")
planck_info["Planck_550"] = ("545", "HFI")
planck_info["Planck_850"] = ("353", "HFI")
planck_info["Planck_1380"] = ("217", "HFI")
planck_info["Planck_2100"] = ("143", "HFI")
planck_info["Planck_3000"] = ("100", "HFI")
planck_info["Planck_4260"] = ("070", "LFI")
planck_info["Planck_6810"] = ("044", "LFI")
planck_info["Planck_10600"] = ("030", "LFI")

# -----------------------------------------------------------------
#  Filter class
# -----------------------------------------------------------------

## An instance of the Filter class represents a particular wavelength band, including its response or transmission
# curve and some basic properties such as its mean and pivot wavelengths. The class provides a function to
# integrate a given spectrum over the band. A filter instance can be constructed by name from one of
# the provided resource files describing specific instruments or standard wavelength bands. Alternatively,
# a filter can be created with a uniform transmission curve over a certain wavelength range.
#
# The precise formalae involved in the integration over the filter and the calculation of the pivot wavelength
# depend on whether the instrument counts photons (photon counter) or measures energy (bolometer).
#
class Filter:
    # ---------- Constructing -------------------------------------

    ## The constructor constructs a new Filter instance in one of two ways, depending on the type of the argument.
    #
    # If \em filterspec is a tuple with two numbers, a bolometer-type filter is constructed with a uniform
    # transmission curve in the indicated range. The (min,max) wavelength values must be expressed in micron.
    #
    # If \em filterspec is a string, the constructor locates a "$PYTHONPATH/dat/filter" directory, and searches
    # in that directory for a VOTable resource file with a name that contains the specified string. If there is
    # exactly one such file, its contents is loaded into the new Filter instance. If there is no such file, or
    # if the specified string matches multiple file names, the constructor raises an error.
    #
    # The VOTable resource files do not properly specify the type of filter (photon counter or bolometer).
    # The constructor classifies filters based on the instrument name in the filter identifier.
    # While this simple heuristic works for now, it is not guaranteed to be correct for future cases.
    #
    # The tables below list the filters available at the time of writing. More filter definitions can be
    # downloaded from the filter profile service web site http://svo2.cab.inta-csic.es/theory/fps.
    #
    #| Filter-spec | \f$\lambda_\text{ctr}\f$ | \f$\lambda_\text{piv}\f$ | Type | Description
    #|------- -----|--------------------------|--------------------------|------|------------
    #| 2MASS/2MASS.J | 1.239 | 1.239 | Pho | 2MASS J
    #| 2MASS/2MASS.H | 1.649 | 1.649 | Pho | 2MASS H
    #| 2MASS/2MASS.Ks | 2.163 | 2.164 | Pho | 2MASS Ks
    #| GALEX/GALEX.FUV | 0.1526 | 0.1535 | Pho | GALEX FUV
    #| GALEX/GALEX.NUV | 0.2329 | 0.2301 | Pho | GALEX NUV
    #| Herschel/Pacs.blue | 71.33 | 70.77 | Bol | Herschel Pacs blue filter
    #| Herschel/Pacs.green | 102.3 | 100.8 | Bol | Herschel Pacs green filter
    #| Herschel/Pacs.red | 166.1 | 161.9 | Bol | Herschel Pacs red filter
    #| Herschel/SPIRE.PSW_ext | 257.6 | 252.5 | Bol | Herschel SPIRE PSW filter (extended sources)
    #| Herschel/SPIRE.PMW_ext | 357.5 | 354.3 | Bol | Herschel SPIRE PMW filter (extended sources)
    #| Herschel/SPIRE.PLW_ext | 518.4 | 515.4 | Bol | Herschel SPIRE PLW filter (extended sources)
    #| IRAS/IRAS.12mu | 11.43 | 11.41 | Pho | IRAS 12 micron
    #| IRAS/IRAS.25mu | 23.97 | 23.61 | Pho | IRAS 25 micron
    #| IRAS/IRAS.60mu | 61.88 | 60.41 | Pho | IRAS 60 micron
    #| IRAS/IRAS.100mu | 101 | 101.1 | Pho | IRAS 100 micron
    #| Misc/MCPS.U | 0.3666 | 0.3646 | Pho | MCPS Johnson U
    #| Misc/MCPS.B | 0.4406 | 0.4434 | Pho | MCPS Johnson B
    #| Misc/MCPS.V | 0.5426 | 0.5493 | Pho | MCPS Johnson V
    #| Misc/MCPS.I | 0.8643 | 0.8739 | Pho | MCPS Johnson I
    #| SLOAN/SDSS.u | 0.3565 | 0.3557 | Pho | SDSS u
    #| SLOAN/SDSS.g | 0.47 | 0.4702 | Pho | SDSS g
    #| SLOAN/SDSS.r | 0.6174 | 0.6176 | Pho | SDSS r
    #| SLOAN/SDSS.i | 0.7534 | 0.749 | Pho | SDSS i
    #| SLOAN/SDSS.z | 0.8782 | 0.8947 | Pho | SDSS z
    #| Spitzer/IRAC.I1 | 3.547 | 3.551 | Pho | IRAC I1
    #| Spitzer/IRAC.I2 | 4.502 | 4.496 | Pho | IRAC I2
    #| Spitzer/IRAC.I3 | 5.716 | 5.724 | Pho | IRAC I3
    #| Spitzer/IRAC.I4 | 7.856 | 7.884 | Pho | IRAC I4
    #| Spitzer/MIPS.24mu | 23.47 | 23.59 | Bol | MIPS 24 microns
    #| Spitzer/MIPS.70mu | 70.52 | 70.89 | Bol | MIPS 70 microns
    #| Spitzer/MIPS.160mu | 156.9 | 155.4 | Bol | MIPS 160 microns
    #| UKIRT/UKIDSS.Z | 0.8825 | 0.8826 | Pho | UKIDSS Z
    #| UKIRT/UKIDSS.Y | 1.03 | 1.031 | Pho | UKIDSS Y
    #| UKIRT/UKIDSS.J | 1.249 | 1.25 | Pho | UKIDSS J
    #| UKIRT/UKIDSS.H | 1.638 | 1.635 | Pho | UKIDSS H
    #| UKIRT/UKIDSS.K | 2.206 | 2.206 | Pho | UKIDSS K
    #| WISE/WISE.W1 | 3.466 | 3.39 | Pho | WISE W1 filter
    #| WISE/WISE.W2 | 4.644 | 4.641 | Pho | WISE W2 filter
    #| WISE/WISE.W3 | 13.22 | 12.57 | Pho | WISE W3 filter
    #| WISE/WISE.W4 | 22.22 | 22.31 | Pho | WISE W4 filter
    #| Generic/Johnson.U | Johnson U
    #| Generic/Johnson.B | Johnson B
    #| Generic/Johnson.V | Johnson V
    #| Generic/Johnson.R | Johnson R
    #| Generic/Johnson.I | Johnson I
    #| Generic/Johnson.J | Johnson J
    #| Generic/Johnson.M | Johnson M
    def __init__(self, filterspec, name=None):

        # Planck filters have to be handled seperately
        if isinstance(filterspec, types.StringTypes) and "planck" in filterspec.lower():

            # Import Astropy stuff
            from astropy.units import Unit, spectral

            this_path = os.path.dirname(os.path.abspath(__file__))
            core_path = os.path.dirname(this_path)
            planck_transmissions_path = os.path.join(core_path, "dat", "filters", "Planck")

            # Get info
            frequency_string = planck_info[filterspec][0]
            instrument = planck_info[filterspec][1]

            filename = instrument + "_BANDPASS_F" + frequency_string + ".txt"
            filepath = os.path.join(planck_transmissions_path, filename)

            wavenumbers, transmissions = np.loadtxt(filepath, unpack=True, skiprows=2, usecols=(0,1))

            # Remove zero from the wavenumbers
            if wavenumbers[0] == 0:
                wavenumbers = wavenumbers[1:]
                transmissions = transmissions[1:]

            # HFI
            if instrument == "HFI":
                wavenumbers = wavenumbers * Unit("1/cm")
                wavelengths = (1.0/wavenumbers).to("micron").value

            # LFI
            else:
                frequencies = wavenumbers * Unit("GHz")
                wavelengths = frequencies.to("micron", equivalencies=spectral()).value

            # REVERSE
            wavelengths = np.flipud(wavelengths)
            transmissions = np.flipud(transmissions)

            self._WavelengthMin = np.min(wavelengths)
            self._WavelengthMax = np.max(wavelengths)
            self._WavelengthCen = None

            self._Wavelengths = wavelengths
            self._Transmission = transmissions

            # Calculate integrals
            integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission)
            integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission / (self._Wavelengths ** 2))

            integral3 = np.trapz(x=self._Wavelengths, y=self._Transmission * self._Wavelengths)

            # Mean wavelength, ID and description
            self._WavelengthMean = integral3 / integral1
            self._WavelengthEff = None
            self._FilterID = "Planck/" + instrument + "." + frequency_string
            self._Description = instrument + " " + frequency_string

            # Not photon counter
            self._PhotonCounter = False

            # Calculate integrated transmission and pivot wavelength
            self._IntegratedTransmission = integral1
            self._WavelengthPivot = np.sqrt(integral1 / integral2)

            self._EffWidth = None

            self.true_filter = True

            # Success
            return

        # string --> load from appropriate resource file
        elif isinstance(filterspec, types.StringTypes):
            for pythondir in sys.path:
                filterdir = os.path.join(pythondir, "pts", "core", "dat", "filters")
                if os.path.isdir(filterdir):

                    filterfiles = filter(lambda fn: fn.endswith(".xml") and filterspec in fn, os.listdir(filterdir))
                    if len(filterfiles) > 1: raise ValueError("filter spec " + filterspec + " is ambiguous")
                    if len(filterfiles) < 1: raise ValueError("no filter found with spec " + filterspec)

                    # load the XML tree
                    tree = etree.parse(open(os.path.join(filterdir, filterfiles[0]), 'r'))

                    # verify the wavelength unit to be Angstrom
                    unit = tree.xpath("//RESOURCE/PARAM[@name='WavelengthUnit'][1]/@value")[0]
                    if unit!='Angstrom': raise ValueError("VOTable uses unsupported unit: " + unit)

                    # load some basic properties (converting from Angstrom to micron)
                    self._WavelengthMin = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMin'][1]/@value")[0])
                    self._WavelengthMax = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMax'][1]/@value")[0])
                    self._WavelengthCen = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthCen'][1]/@value")[0])
                    self._WavelengthMean = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthMean'][1]/@value")[0])
                    self._WavelengthEff = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WavelengthEff'][1]/@value")[0])
                    self._FilterID = tree.xpath("//RESOURCE/PARAM[@name='filterID'][1]/@value")[0]
                    self._Description = tree.xpath("//RESOURCE/PARAM[@name='Description'][1]/@value")[0]
                    self._Description = self._Description.replace("&#956;m", "micron")
                    self._EffWidth = 1e-4*float(tree.xpath("//RESOURCE/PARAM[@name='WidthEff'][1]/@value")[0])

                    # load the transmission table (converting wavelengths from Angstrom to micron)
                    values = np.array(tree.xpath("//RESOURCE/TABLE/DATA/TABLEDATA[1]/TR/TD/text()"), dtype=float)
                    if len(values)<4: raise ValueError("transmission table not found in filter definition")
                    self._Wavelengths,self._Transmission = np.reshape(values, (-1,2)).T
                    self._Wavelengths *= 1e-4

                    # determine the filter type (there seems to be no better heuristic than using the instrument name)
                    self._PhotonCounter = not any(["/"+x in self._FilterID.lower() for x in ("mips","pacs","spire")])

                    # calculate the pivot wavelength
                    if self._PhotonCounter:
                        integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission*self._Wavelengths)
                        integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission/self._Wavelengths)
                    else:
                        integral1 = np.trapz(x=self._Wavelengths, y=self._Transmission)
                        integral2 = np.trapz(x=self._Wavelengths, y=self._Transmission/(self._Wavelengths**2))
                    self._IntegratedTransmission = integral1
                    self._WavelengthPivot = np.sqrt(integral1/integral2)

                    self.true_filter = True

                    # report success
                    return

            raise ValueError("filter resource path not found")

        # range --> construct ad hoc uniform bolometer
        else:
            if name is None: name = "Uniform"
            self._WavelengthMin, self._WavelengthMax = map(float,filterspec)
            self._WavelengthCen = 0.5 * (self._WavelengthMin + self._WavelengthMax)
            self._WavelengthMean = self._WavelengthCen
            self._WavelengthEff = self._WavelengthCen
            self._FilterID = name + "_[{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._Description = name + " filter in range [{},{}]".format(self._WavelengthMin,self._WavelengthMax)
            self._Wavelengths = np.array((self._WavelengthMin,self._WavelengthMax))
            self._Transmission = np.ones((2,))
            self._PhotonCounter = False
            self._IntegratedTransmission = self._WavelengthMax - self._WavelengthMin
            self._WavelengthPivot = np.sqrt(self._WavelengthMin * self._WavelengthMax)
            self._EffWidth = self._WavelengthMax - self._WavelengthMin
            self.true_filter = False

    # ---------- Alternative construction -------------------------------

    @classmethod
    def from_string(cls, name):

        """
        This function allows constructing a Filter object from a less strict name.
        :param name:
        :return:
        """

        # First try an exact match
        #try: return cls(name)
        #except ValueError: pass

        # Define the different possible names for each of the filters
        galex_fuv_names = ["GALEX.FUV", "GALEX FUV", "FUV", "the FUV-band", "GALEX_FUV"]
        galex_nuv_names = ["GALEX.NUV", "GALEX NUV", "NUV", "the NUV-band", "GALEX_NUV"]
        sdss_u_names = ["SDSS.u", "SDSS u", "u", "SDSSu", "the u-band", "SDSS_u"]
        sdss_g_names = ["SDSS.g", "SDSS g", "g", "SDSSg", "the g-band", "SDSS_g"]
        sdss_r_names = ["SDSS.r", "SDSS r", "r", "SDSSr", "the r-band", "SDSS_r"]
        sdss_i_names = ["SDSS.i", "SDSS i", "i", "SDSSi", "the i-band", "SDSS_i"]
        sdss_z_names = ["SDSS.z", "SDSS z", "z", "SDSSz", "the z-band", "SDSS_z"]
        mass_h_names = ["2MASS.H", "2MASS H", "H", "SDSSh", "the H-band", "2MASS_H"]
        mass_j_names = ["2MASS.J", "2MASS J", "J", "SDSSj", "the J-band", "2MASS_J"]
        mass_k_names = ["2MASS.K", "2MASS.Ks", "2MASS K", "2MASS Ks", "K", "Ks", "the K-band", "2MASS_K", "2MASS_Ks"]
        irac_i1_names = ["IRAC.I1", "IRAC I1", "IRAC 3.6", "IRAC 3.6um", "IRAC 3.6mu", "I1", "IRAC1", "IRAC-1", "the IRAC-1 band", "Spitzer 3.6", "IRAC_I1", "Spitzer_3.6"]
        irac_i2_names = ["IRAC.I2", "IRAC I2", "IRAC 4.5", "IRAC 4.5um", "IRAC 4.5mu", "I2", "IRAC2", "IRAC-2", "the IRAC-2 band", "Spitzer 4.5", "IRAC_I2", "Spitzer_4.5"]
        irac_i3_names = ["IRAC.I3", "IRAC I3", "IRAC 5.8", "IRAC 5.8um", "IRAC 5.8mu", "I3", "IRAC3", "IRAC-3", "the IRAC-3 band", "Spitzer 5.8", "IRAC_I3", "Spitzer_5.8"]
        irac_i4_names = ["IRAC.I4", "IRAC I4", "IRAC 8.0", "IRAC 8.0um", "IRAC 8.0mu", "I4", "IRAC4", "IRAC-4", "the IRAC-4 band", "Spitzer 8.0", "IRAC_I4", "Spitzer_8.0"]
        wise_w1_names = ["WISE.W1", "WISE W1", "W1", "WISE1", "WISE-1", "the WISE-1 band", "w1", "WISE 3.4", "WISE_W1", "WISE_3.4"]
        wise_w2_names = ["WISE.W2", "WISE W2", "W2", "WISE2", "WISE-2", "the WISE-2 band", "w2", "WISE 4.6", "WISE_W2", "WISE_4.6"]
        wise_w3_names = ["WISE.W3", "WISE W3", "W3", "WISE3", "WISE-3", "the WISE-3 band", "w3", "WISE 12", "WISE_W3", "WISE_12"]
        wise_w4_names = ["WISE.W4", "WISE W4", "W4", "WISE4", "WISE-4", "the WISE-4 band", "w4", "WISE 22", "WISE_W4", "WISE_22"]
        mips_24_names = ["MIPS.24mu", "MIPS.24um", "MIPS.24", "MIPS 24mu", "MIPS 24um", "MIPS 24", "24mu", "24um", "MIPS-24", "the MIPS-24 band", "Spitzer 24", "MIPS_24mu", "MIPS_24um", "Spitzer_24"]
        mips_70_names = ["MIPS.70mu", "MIPS.70um", "MIPS.70", "MIPS 70mu", "MIPS 70um", "MIPS 70", "MIPS-70", "the MIPS-70 band", "Spitzer 70", "MIPS_70mu", "MIPS_70um", "Spitzer_70"]
        mips_160_names = ["MIPS.160mu", "MIPS.160um", "MIPS.160", "MIPS 160mu", "MIPS 160um", "MIPS 160", "MIPS-160", "the MIPS-160 band", "Spitzer 160", "MIPS_160mu", "MIPS_160um", "Spitzer_160"]
        pacs_blue_names = ["Pacs.blue", "PACS.BLUE", "PACS blue", "PACS BLUE", "Pacs 70mu", "Pacs 70um", "PACS 70mu", "PACS 70um", "PACS-70", "the PACS-70 band", "Pacs 70", "Pacs blue", "PACS 70", "Pacs_blue", "PACS_BLUE", "PACS_70"]
        pacs_green_names = ["Pacs.green", "PACS.GREEN", "PACS green", "PACS GREEN", "Pacs 100mu", "Pacs 100um", "PACS 100mu", "PACS 100um", "PACS-100", "the PACS-100 band", "Pacs 100", "Pacs green", "PACS 100", "Pacs_green", "PACS_GREEN", "PACS_100"]
        pacs_red_names = ["Pacs.red", "PACS.RED", "PACS red", "PACS RED", "Pacs 160mu", "Pacs 160um", "PACS 160mu", "PACS 160um", "PACS-160", "the PACS-160 band", "Pacs 160", "Pacs red", "PACS 160", "Pacs_red", "PACS_RED", "PACS_160"]
        spire_psw_names = ["SPIRE.PSW", "SPIRE PSW", "SPIRE 250mu", "SPIRE 250um", "SPIRE-250", "the SPIRE-250 band", "SPIRE.PSW_ext", "SPIRE 250", "SPIRE PSW_ext", "SPIRE_PSW", "SPIRE_250"]
        spire_pmw_names = ["SPIRE.PMW", "SPIRE PMW", "SPIRE 350mu", "SPIRE 350um", "SPIRE-350", "the SPIRE-250 band", "SPIRE.PMW_ext", "SPIRE 350", "SPIRE PMW_ext", "SPIRE_PMW", "SPIRE_350"]
        spire_plw_names = ["SPIRE.PLW", "SPIRE PLW", "SPIRE 500mu", "SPIRE 500um", "SPIRE-500", "the SPIRE-500 band", "SPIRE.PLW_ext", "SPIRE 500", "SPIRE PLW_ext", "SPIRE_PLW", "SPIRE_500"]

        planck_350_names = ["Planck 350", "Planck_350", "Planck 857 GHz"]
        planck_550_names = ["Planck 550", "Planck_550", "Planck 545 GHz"]
        planck_850_names = ["Planck 850", "Planck_850", "Planck 353 GHz"]
        planck_1380_names = ["Planck 1380", "Planck_1380", "Planck 217 GHz"]
        planck_2100_names = ["Planck 2100", "Planck_2100", "Planck 143 GHz"]
        planck_3000_names = ["Planck 3000", "Planck_3000", "Planck 100 GHz"]
        planck_4260_names = ["Planck 4260", "Planck_4260", "Planck 70 GHz"]
        planck_6810_names = ["Planck 6810", "Planck_6810", "Planck 44 GHz"]
        planck_10600_names = ["Planck 10600", "Planck_10600", "Planck 30 GHz"]

        # Generic filters
        johnson_u_names = ["Johnson U", "U"]
        johnson_b_names = ["Johnson B", "B"]
        johnson_v_names = ["Johnson V", "V"]
        johnson_r_names = ["Johnson R", "R"]
        johnson_i_names = ["Johnson I", "I"]

        # IRAS filters
        iras_12_names = ["IRAS.12", "IRAS 12", "IRAS.12um", "IRAS 12um", "IRAS.12mu", "IRAS 12mu", "IRAS_12"]
        iras_25_names = ["IRAS.25", "IRAS 25", "IRAS.25um", "IRAS 25um", "IRAS.25mu", "IRAS 25mu", "IRAS_25"]
        iras_60_names = ["IRAS 60", "IRAS 60", "IRAS.60um", "IRAS 60um", "IRAS.60mu", "IRAS 60mu", "IRAS_60"]
        iras_100_names = ["IRAS 100", "IRAS 100", "IRAS.100um", "IRAS 100um", "IRAS.100mu", "IRAS 100mu", "IRAS_100"]

        # Ha
        if name == "656_1": return cls(name)
        ha_names = ["Ha", "H alpha", "H-alpha", "H-a", "Halpha", "H_alpha", "Mosaic Halpha", "H_alph"]

        # SWIFT UVOT filters
        uvot_u_names = ["UVOT U", "Swift U", "Swift_U"]
        uvot_b_names = ["UVOT B", "Swift B", "Swift_B"]
        uvot_v_names = ["UVOT V", "Swift V", "Swift_V"]
        uvot_uvw2_names = ["UVOT W2", "UVW2", "W2", "Swift W2", "SWIFT W2", "UVOT_W2", "Swift_"]
        uvot_uvm2_names = ["UVOT M2", "UVM2", "M2", "Swift M2", "SWIFT M2", "UVOT_M2"]
        uvot_uvw1_names = ["UVOT W1", "UVW1", "W1", "Swift W1", "SWIFT W1", "UVOT_W1"]

        # Select the right filter
        if name in galex_fuv_names: return cls("GALEX.FUV")
        elif name in galex_nuv_names: return cls("GALEX.NUV")
        elif name in sdss_u_names: return cls("SDSS.u")
        elif name in sdss_g_names: return cls("SDSS.g")
        elif name in sdss_r_names: return cls("SDSS.r")
        elif name in sdss_i_names: return cls("SDSS.i")
        elif name in sdss_z_names: return cls("SDSS.z")
        elif name in mass_h_names: return cls("2MASS.H")
        elif name in mass_j_names: return cls("2MASS.J")
        elif name in mass_k_names: return cls("2MASS.Ks")
        elif name in irac_i1_names: return cls("IRAC.I1")
        elif name in irac_i2_names: return cls("IRAC.I2")
        elif name in irac_i3_names: return cls("IRAC.I3")
        elif name in irac_i4_names: return cls("IRAC.I4")
        elif name in wise_w1_names: return cls("WISE.W1")
        elif name in wise_w2_names: return cls("WISE.W2")
        elif name in wise_w3_names: return cls("WISE.W3")
        elif name in wise_w4_names: return cls("WISE.W4")
        elif name in mips_24_names: return cls("MIPS.24mu")
        elif name in mips_70_names: return cls("MIPS.70mu")
        elif name in mips_160_names: return cls("MIPS.160mu")
        elif name in pacs_blue_names: return cls("Pacs.blue")
        elif name in pacs_green_names: return cls("Pacs.green")
        elif name in pacs_red_names: return cls("Pacs.red")
        elif name in spire_psw_names: return cls("SPIRE.PSW_ext")
        elif name in spire_pmw_names: return cls("SPIRE.PMW_ext")
        elif name in spire_plw_names: return cls("SPIRE.PLW_ext")
        elif name in johnson_u_names: return cls("Johnson.U")
        elif name in johnson_b_names: return cls("Johnson.B")
        elif name in johnson_v_names: return cls("Johnson.V")
        elif name in johnson_r_names: return cls("Johnson.R")
        elif name in johnson_i_names: return cls("Johnson.I")
        elif name in iras_12_names: return cls("IRAS.12mu")
        elif name in iras_25_names: return cls("IRAS.25mu")
        elif name in iras_60_names: return cls("IRAS.60mu")
        elif name in iras_100_names: return cls("IRAS.100mu")
        elif name in ha_names: return cls("Halpha")
        elif name in uvot_uvw2_names: return cls("UVOT.UVW2")
        elif name in uvot_uvm2_names: return cls("UVOT.UVM2")
        elif name in uvot_uvw1_names: return cls("UVOT.UVW1")
        elif name in planck_350_names: return cls("Planck_350")
        elif name in planck_550_names: return cls("Planck_550")
        elif name in planck_850_names: return cls("Planck_850")
        elif name in planck_1380_names: return cls("Planck_1380")
        elif name in planck_2100_names: return cls("Planck_2100")
        elif name in planck_3000_names: return cls("Planck_3000")
        elif name in planck_4260_names: return cls("Planck_4260")
        elif name in planck_6810_names: return cls("Planck_6810")
        elif name in planck_10600_names: return cls("Planck_10600")
        else: raise ValueError("No corresponding filter found (name = " + name + ")")

    @classmethod
    def from_instrument_and_band(cls, instrument, band):
        return cls.from_string(instrument + "." + band)

    # ---------- Retrieving information -------------------------------

    @property
    def skirt_description(self): # returns the name as defined in the SKIRT LuminosityStellarCompNormalization class

        if self.name == "GALEX.FUV": return "FUV"
        elif self.name == "GALEX.NUV": return "NUV"
        elif self.name == "Johnson.U": return "U"
        elif self.name == "Johnson.B": return "B"
        elif self.name == "Johnson.V": return "V"
        elif self.name == "Johnson.R": return "R"
        elif self.name == "Johnson.I": return "I"
        elif self.name == "2MASS.J": return "J"
        elif self.name == "2MASS.H": return "H"
        elif self.name == "2MASS.Ks": return "K"
        elif self.name == "SDSS.u": return "SDSSu"
        elif self.name == "SDSS.g": return "SDSSg"
        elif self.name == "SDSS.r": return "SDSSr"
        elif self.name == "SDSS.i": return "SDSSi"
        elif self.name == "SDSS.z": return "SDSSz"
        elif self.name == "IRAC.I1": return "IRAC1"
        elif self.name == "IRAC.I2": return "IRAC2"
        elif self.name == "WISE.W1": return "WISE1"
        elif self.name == "WISE.W2": return "WISE2"
        else: raise ValueError("The band " +  self.name + " is not defined in SKIRT")

    @property
    def aniano_name(self): # Returns the name as appearing in the Aniano kernel and psf FITS files

        if self.name == "GALEX.FUV": return "GALEX_FUV"
        elif self.name == "GALEX.NUV": return "GALEX_NUV"
        elif self.name == "IRAC.I1": return "IRAC_3.6"
        elif self.name == "IRAC.I2": return "IRAC_4.5"
        elif self.name == "IRAC.I3": return "IRAC_5.8"
        elif self.name == "IRAC.I4": return "IRAC_8.0"
        elif self.name == "WISE.W1": return "WISE_FRAME_3.4"
        elif self.name == "WISE.W2": return "WISE_FRAME_4.6"
        elif self.name == "WISE.W3": return "WISE_FRAME_11.6"
        elif self.name == "WISE.W4": return "WISE_FRAME_22.1"
        elif self.name == "MIPS.24mu": return "MIPS_24"
        elif self.name == "MIPS.70mu": return "MIPS_70"
        elif self.name == "MIPS.160mu": return "MIPS_160"
        elif self.name == "Pacs.blue": return "PACS_70"
        elif self.name == "Pacs.green": return "PACS_100"
        elif self.name == "Pacs.red": return "PACS_160"
        elif self.name == "SPIRE.PSW_ext": return "SPIRE_250"
        elif self.name == "SPIRE.PMW_ext": return "SPIRE_350"
        elif self.name == "SPIRE.PLW_ext": return "SPIRE_500"
        else: raise ValueError("The band " + self.name + " is not defined for the Aniano set of kernels")

    @property
    def name(self):
        if self.true_filter: return self._FilterID.split("/")[1]
        else: return self._FilterID.split("_")[0]

    @property
    def observatory(self):
        if self.true_filter: return self._FilterID.split("/")[0]
        else: return None

    @property
    def instrument(self):
        if self.true_filter: return self._FilterID.split("/")[1].split(".")[0]
        else: return None

    @property
    def band(self):
        if self.true_filter: return self._FilterID.split("/")[1].split(".")[1].replace("_ext", "")
        elif self._FilterID.startswith("Uniform"): return None
        else: return self._FilterID.split("_")[0]

    ## This function returns the filter in string format ('instrument' 'band')
    def __str__(self):
        if self.true_filter: return self.instrument + " " + self.band
        else: return self.name

    ## This function returns a unique identifier for the filter.
    def filterID(self):
        return self._FilterID

    ## This function returns a human-readable description for the filter.
    def description(self):
        return self._Description

    ## This function returns the mean wavelength for the filter, in micron.
    def meanwavelength(self):
        return self._WavelengthMean

    ## This function returns the mean wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def mean(self):
        from astropy.units import Unit
        return self.meanwavelength() * Unit("micron")

    ## This function returns the effective wavelength for the filter, in micron.
    def effectivewavelength(self):
        return self._WavelengthEff

    ## This function returns the effective wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def effective(self):
        from astropy.units import Unit
        return self.effectivewavelength() * Unit("micron")

    ## This function returns the minimum wavelength for the filter, in micron.
    def minwavelength(self):
        return self._WavelengthMin

    ## This function returns the minimum wavelength for the filter as a quantity. Astropy is imported inside this
    # function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def min(self):
        from astropy.units import Unit
        return self.minwavelength() * Unit("micron")

    ## This function returns the maximum wavelength for the filter, in micron.
    def maxwavelength(self):
        return self._WavelengthMax

    ## This function returns the maximum wavelength for the filter as a quantity. Astropy is imported inside this
    #  function to support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def max(self):
        from astropy.units import Unit
        return self.maxwavelength() * Unit("micron")

    ## This function returns the center wavelength for the filter, in micron. The center wavelength is
    # defined as the wavelength halfway between the two points for which filter response or transmission
    # (depending on the filter type) is half maximum.
    def centerwavelength(self):
        return self._WavelengthCen

    ## This function returns the center wavelength as a quantity. Astropy is imported inside this function to
    #  support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def center(self):
        from astropy.units import Unit
        return self.centerwavelength() * Unit("micron")

    ## This function returns the pivot wavelength for the filter, in micron. The pivot wavelength is defined
    # as the wavelength that connects the filter-averaged wavelength and frequency-style fluxes through
    # \f$\left<F_\nu\right> = \left<F_\lambda\right>\lambda_\text{pivot}^2/c\f$. The value depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int\lambda R(\lambda) \,\mathrm{d}\lambda }
    #     {  \int R(\lambda) \,\mathrm{d}\lambda/\lambda } }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \lambda_\mathrm{pivot} = \sqrt{ \frac{ \int T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda/\lambda^2 } }. \f]
    def pivotwavelength(self):
        return self._WavelengthPivot

    ## This function returns the pivot wavelength as a quantity. Astropy is imported inside this function to
    #  support PTS installations without Astropy for users that don't use this (new) function.
    @property
    def pivot(self):
        from astropy.units import Unit
        return self.pivotwavelength() * Unit("micron")

    ## This function returns the effective bandwith, in micron.
    def effective_bandwidth(self):
        return self._EffWidth

    # ---------- Integrating --------------------------------------

    ## This function calculates and returns the filter-averaged value \f$\left<F_\lambda\right>\f$ for a given
    # spectral energy distribution \f$F_\lambda(\lambda)\f$. The calculation depends
    # on the filter type. For a photon counter with response curve \f$R(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int\lambda F_\lambda(\lambda)R(\lambda) \,\mathrm{d}\lambda }
    #     { \int\lambda R(\lambda) \,\mathrm{d}\lambda }. \f]
    # For a bolometer with transmission curve \f$T(\lambda)\f$,
    # \f[ \left<F_\lambda\right> = \frac{ \int F_\lambda(\lambda)T(\lambda) \,\mathrm{d}\lambda }
    #     { \int T(\lambda) \,\mathrm{d}\lambda }. \f]
    #
    # The quantities \f$F_\lambda(\lambda)\f$ must be expressed per unit of wavelength (and \em not, for example,
    # per unit of frequency). The resulting \f$\left<F_\lambda\right>\f$ has the same units as the input distribition.
    # \f$F_\lambda(\lambda)\f$ can be expressed in any units (as long as it is per unit of wavelength) and it can
    # represent various quantities; for example a flux density, a surface density, or a luminosity density.
    #
    # The function accepts two arguments:
    # - \em wavelengths: a numpy array specifying the wavelengths \f$\lambda_\ell\f$, in micron, in increasing order,
    #   on which the spectral energy distribution is sampled. The integration is performed on a wavelength grid that
    #   combines the grid points given here with the grid points on which the filter response or transmission curve
    #   is defined.
    # - \em densities: a numpy array specifying the spectral energy distribution(s) \f$F_\lambda(\lambda_\ell)\f$
    #   per unit of wavelength. This can be an array with the same length as \em wavelengths, or a multi-dimensional
    #   array where the last dimension has the same length as \em wavelengths.
    #   The returned result will have the shape of \em densities minus the last (or only) dimension.
    def convolve(self, wavelengths, densities):
        # define short names for the involved wavelength grids
        wa = wavelengths
        wb = self._Wavelengths

        # create a combined wavelength grid, restricted to the overlapping interval
        w1 = wa[ (wa>=wb[0]) & (wa<=wb[-1]) ]
        w2 = wb[ (wb>=wa[0]) & (wb<=wa[-1]) ]
        w = np.unique(np.hstack((w1,w2)))
        if len(w) < 2: return 0

        # log-log interpolate SED and transmission on the combined wavelength grid
        # (use scipy interpolation function for SED because np.interp does not support broadcasting)
        F = np.exp(interp1d(np.log(wa), _log(densities), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._Transmission), left=0., right=0.))

        # perform the integration
        if self._PhotonCounter:
            return np.trapz(x=w, y=w*F*T) / self._IntegratedTransmission
        else:
            return np.trapz(x=w, y=F*T) / self._IntegratedTransmission

    ## This function calculates and returns the integrated value for a given spectral energy distribution over the
    #  filter's wavelength range,
    def integrate(self, wavelengths, densities):
        # define short names for the involved wavelength grids
        wa = wavelengths
        wb = self._Wavelengths

        # create a combined wavelength grid, restricted to the overlapping interval
        w1 = wa[(wa >= wb[0]) & (wa <= wb[-1])]
        w2 = wb[(wb >= wa[0]) & (wb <= wa[-1])]
        w = np.unique(np.hstack((w1, w2)))
        if len(w) < 2: return 0

        # log-log interpolate SED and transmission on the combined wavelength grid
        # (use scipy interpolation function for SED because np.interp does not support broadcasting)
        F = np.exp(interp1d(np.log(wa), _log(densities), copy=False, bounds_error=False, fill_value=0.)(np.log(w)))
        T = np.exp(np.interp(np.log(w), np.log(wb), _log(self._Transmission), left=0., right=0.))

        # perform the integration
        if self._PhotonCounter: return np.trapz(x=w, y=w * F * T)
        else: return np.trapz(x=w, y=F * T)

## This private helper function returns the natural logarithm for positive values, and a large negative number
# (but not infinity) for zero or negative values.
def _log(X):
    zeromask = X<=0
    logX = np.empty(X.shape)
    logX[zeromask] = -750.  # the smallest (in magnitude) negative value x for which np.exp(x) returns zero
    logX[~zeromask] = np.log(X[~zeromask])
    return logX

# -----------------------------------------------------------------
