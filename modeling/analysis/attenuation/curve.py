#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.attenuation Contains the AttenuationCurveAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AttenuationAnalysisComponent
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....core.data.sed import SED
from ....core.data.attenuation import AttenuationCurve, SMCAttenuationCurve, MilkyWayAttenuationCurve, CalzettiAttenuationCurve, BattistiAttenuationCurve, MappingsAttenuationCurve
from ....core.plot.attenuation import AttenuationPlotter

# -----------------------------------------------------------------

# IR-through-UV extinction curve (Fitzpatrick+, 2007):
# http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJ/663/320
#  - J/ApJ/663/320/stars: Survey stars (tables 1, 3 and 4 of paper) (328 rows)

# UVOT imaging of M81 and Holmberg IX (Hoversten+, 2011):
#  - J/AJ/141/205/table3: Source photometry and fitted parameters
#    http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/AJ/141/205/table3
#  -> has "AV": Best fitting internal extinction A_V

# -----------------------------------------------------------------

class AttenuationCurveAnalyser(AttenuationAnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationCurveAnalyser, self).__init__(config)

        # -- Attributes --

        # The SEDs
        self.total_sed = None
        self.transparent_sed = None

        # The attenuation curves
        self.attenuation_sfr = None
        self.attenuation_diffuse = None
        self.attenuation_total = None

        # The reference attenuation curves
        self.references = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the simulated SEDs
        self.load_seds()

        # 3. Calculate the attenuation curve
        self.calculate_attenuation()

        # 4. Get the reference attenuation curves
        self.get_references()

        # 5. Writing
        self.write()

        # 6. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AttenuationCurveAnalyser, self).setup()

    # -----------------------------------------------------------------

    def load_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated total and transparent SED ...")

        # Determine the path to the SED file corresponding to the 'earth' instrument
        sed_path = fs.join(self.analysis_run.out_path, self.galaxy_name + "_earth_sed.dat")

        # Load the total SED
        self.total_sed = SED.from_skirt(sed_path, contribution="total")

        # Load the transparent SED
        self.transparent_sed = SED.from_skirt(sed_path, contribution="transparent")

    # -----------------------------------------------------------------

    def calculate_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the attenuations ...")

        # Calculate the diffuse attenuation
        self.calculate_diffuse_attenuation()

        # Calculate the SFR attenuation
        self.calculate_sfr_attenuation()

        # Calculate the total attenuation
        self.calculate_total_attenuation()

    # -----------------------------------------------------------------

    def calculate_diffuse_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the diffuse stellar contribution to the attenuation curve ...")

        # Calculate the attenuations of the diffuse dust component
        self.attenuation_diffuse = AttenuationCurve.from_seds(self.total_sed, self.transparent_sed)

        # Find the V-band attenuation for the diffuse dust component
        v_band_attenuation_diffuse = self.attenuation_diffuse.attenuation_at(self.v_band_wavelength)

        # Debugging
        log.debug("The V-band attenuation from diffuse dust is " + str(v_band_attenuation_diffuse))

        # Normalize the attenuation curve to the V-band attenuation
        self.attenuation_diffuse.normalize_at(self.v_band_wavelength)

    # -----------------------------------------------------------------

    def calculate_sfr_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the SFR contribution to the attenuation curve ...")

        # Get the dust mass in SFR clouds for the model
        dust_mass = self.analysis_run.model.sfr_dust_mass

        # Get the dust mass column density
        average_column_density = (dust_mass / self.truncation_area).to("Msun / kpc2").value

        # V band attenuation = total dust mass in SFR clouds / total surface of M31 / 1e5 * 0.67
        a_v_mappings = 0.67 * average_column_density / 1e5

        # Creat the attenuation curve for the SFR dust
        self.attenuation_sfr = MappingsAttenuationCurve(attenuation=a_v_mappings, wavelength=self.v_band_wavelength)

        # Find the V-band attenuation for the SFR dust component
        v_band_attenuation_sfr = self.attenuation_sfr.attenuation_at(self.v_band_wavelength)

        # Debugging
        log.debug("The V-band attenuation from dust in SF clouds is " + str(v_band_attenuation_sfr))

        # Normalize the attenuation curve to the V-band attenuation
        self.attenuation_sfr.normalize_at(self.v_band_wavelength)

    # -----------------------------------------------------------------

    def calculate_total_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the total attenuation curve ...")

        wavelengths = None

        # load the MAPPINGS SED
        input = np.loadtxt(refSED_MAPPINGS)
        flux_mappings = input[:, 1]  # in Jy

        #atts_mappings = A_interpfunc(modWls)  # Attenuation from dust in star forming regions

        atts_mappings = []
        for wavelength in wavelengths:
            att_mappings = self.attenuation_sfr.attenuation_at(wavelength)
            atts_mappings.append(att_mappings)
        atts_mappings = np.array(atts_mappings)

        delta_flux_mappings = flux_mappings * (10 ** (atts_mappings / 2.5) - 1)  # additional flux from mappings attenuation

        # SED of transparent + delta flux mappings
        transparent_sfr_sed = self.transparent_sed + delta_flux_mappings

        atts_tot = -2.5 * np.log10(modFlux / (modTrans + delta_flux_mappings))
        interpfunc = interpolate.interp1d(modWls, atts_tot, kind='linear')
        att_tot_V = interpfunc(0.55)  # Interpolate to find attenuation at V band central wavelengths
        n_atts_tot = atts_tot / att_tot_V

        # Calculate the attenuations for all the dust
        # self.attenuation_total = AttenuationCurve.from_seds(self.total_sed, transparent_sfr_sed)
        # Find the total V-band attenuation
        # v_band_attenuation_total = self.attenuation_total.attenuation_at(v_band_wavelength)


        # print(' V band attenuation from star-forming regions:   ' + str(att_mappings_V))

        # print(' total V band attenuation:          ' + str(att_tot_V))

        # proportion = delta_flux_mappings / (modTrans - modFlux)

        # print(' V band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls - 0.55))]))
        # print(' NUV band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls - 0.227))]))
        # print(' FUV band energy fraction of mappings attenuation to total: ' + str(proportion[np.argmin(np.abs(modWls - 0.153))]))
        # print(' Mean energy fraction of mappings attenuation to total:   ' + str(np.mean(proportion[modWls < 5])))
        # print(' Median energy fraction of mappings attenuation to total: ' + str(np.median(proportion[modWls < 5])))

    # -----------------------------------------------------------------

    def get_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the reference attenuation curves ...")

        # Load the Milky Way, SMC and Calzetti attenuation curves
        self.references["Milky Way"] = MilkyWayAttenuationCurve()
        self.references["SMC"] = SMCAttenuationCurve()
        self.references["Calzetti"] = CalzettiAttenuationCurve()
        self.references["Battisti"] = BattistiAttenuationCurve()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the attenuation curves
        self.write_attenuations()

    # -----------------------------------------------------------------

    def write_attenuations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the attenuations ...")

        # Determine the path to the diffuse dust attenuation file
        diffuse_path = fs.join(self.attenuation_curve_path, "diffuse.dat")

        # Write the diffuse dust attenuations
        self.attenuation_diffuse.saveto(diffuse_path)

        # Determine the path to the SFR attenuation file
        sfr_path = fs.join(self.attenuation_curve_path, "sfr.dat")

        # Write the SFR attenuation curve
        self.attenuation_sfr.saveto(sfr_path)

        # Determine the path to the total dust attenuation file
        total_path = fs.join(self.attenuation_curve_path, "total.dat")

        # Write the total attenuation curve
        self.attenuation_total.saveto(total_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Create a new AttenuationPlotter instance
        plotter = AttenuationPlotter()

        # Add the diffuse attenuation curve
        plotter.add_attenuation_curve(self.attenuation_diffuse, "diffuse")

        # Add the SFR attenuation curve
        plotter.add_attenuation_curve(self.attenuation_sfr, "SFR")

        # Add the total attenuation curve
        plotter.add_attenuation_curve(self.attenuation_total, "total")

        # Add the reference attenuation curves
        for name in self.references: plotter.add_attenuation_curve(self.references[name], name)

        # Determine the path to the attenuation plot file
        path = fs.join(self.attenuation_curve_path, "attenuation.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
