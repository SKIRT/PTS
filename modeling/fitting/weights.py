#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.weights Contains the WeightsCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from .tables import WeightsTable
from ...core.tools import sequences
from ...magic.tools import wavelengths
from ...core.tools import formatting as fmt

# -----------------------------------------------------------------

# Define original regime names
uv_name = "uv"
optical_name = "optical"
nir_name = "nir"
mir_name = "mir"
fir_name = "fir"
submm_microwave_name = "submm-microwave"

# List of all original regime names
original_regime_names = [uv_name, optical_name, nir_name, mir_name, fir_name, submm_microwave_name]

# -----------------------------------------------------------------

# Define physical regime names
ionizing_name = wavelengths.ionizing
young_name = wavelengths.young
evolved_name = wavelengths.evolved
mix_name = wavelengths.mix
aromatic_name = wavelengths.aromatic
thermal_name = wavelengths.thermal

# List of all physical regime names
physical_regime_names = wavelengths.physical_regimes

# -----------------------------------------------------------------

# ALL REGIME NAMES
all_regime_names = original_regime_names + physical_regime_names

# -----------------------------------------------------------------

class WeightsCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(WeightsCalculator, self).__init__(*args, **kwargs)

        # The filters for which to calculate the weight
        self.filters = None

        # Create the table to contain the weights
        self.table = WeightsTable()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Calculate weights
        self.calculate()

        # 3. Write
        if self.config.write: self.write()

        # 4. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    @property
    def nfilters(self):

        """
        This function ...
        :return:
        """

        return len(self.filters)

    # -----------------------------------------------------------------

    @property
    def has_filters(self):

        """
        This function ...
        :return:
        """

        return self.nfilters > 0

    # -----------------------------------------------------------------

    @property
    def no_filters(self):

        """
        This function ...
        :return:
        """

        return self.nfilters == 0

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(WeightsCalculator, self).setup(**kwargs)

        # Get the filters
        if kwargs.get("filters", None) is not None: self.filters = kwargs.pop("filters")
        else: self.filters = self.config.filters

        # Check that there are filters
        if self.no_filters: raise ValueError("No filters")

        # Check flags that should not be enabled in physical mode
        if self.config.physical:

            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")

            if self.config.no_uv: raise ValueError("Error")
            if self.config.no_optical: raise ValueError("Error")
            if self.config.no_nir: raise ValueError("Error")
            if self.config.no_mir: raise ValueError("Error")
            if self.config.no_fir: raise ValueError("Error")
            if self.config.no_submm_microwave: raise ValueError("Error")

        # Check flags that should not be enabled in original mode
        else:

            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")

            if self.config.no_ionizing: raise ValueError("Error")
            if self.config.no_young: raise ValueError("Error")
            if self.config.no_evolved: raise ValueError("Error")
            if self.config.no_mix: raise ValueError("Error")
            if self.config.no_aromatic: raise ValueError("Error")
            if self.config.no_thermal: raise ValueError("Error")

    # -----------------------------------------------------------------

    @lazyproperty
    def regimes(self):

        """
        This function ...
        :return:
        """

        # Return the regime names
        if self.config.physical: regimes = self.get_regimes_physical()
        else: regimes = self.get_regimes_original()

        # Check number of regimes
        if len(regimes) == 0: raise ValueError("No regimes")

        # Return
        return regimes

    # -----------------------------------------------------------------

    def get_regimes_original(self):

        """
        This function ...
        :return:
        """

        # Only UV
        if self.config.only_uv:

            if self.config.no_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = [uv_name]

        # Only optical
        elif self.config.only_optical:

            if self.config.no_optical: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = [optical_name]

        # Only NIR
        elif self.config.only_nir:

            if self.config.no_nir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = [nir_name]

        # Only MIR
        elif self.config.only_mir:

            if self.config.no_mir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = [mir_name]

        # Only FIR
        elif self.config.only_fir:

            if self.config.no_fir: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_submm_microwave: raise ValueError("Error")
            regimes = [fir_name]

        # Only submm/microwave
        elif self.config.only_submm:

            if self.config.no_submm: raise ValueError("Error")
            if self.config.only_uv: raise ValueError("Error")
            if self.config.only_optical: raise ValueError("Error")
            if self.config.only_nir: raise ValueError("Error")
            if self.config.only_mir: raise ValueError("Error")
            if self.config.only_fir: raise ValueError("Error")
            regimes = [submm_microwave_name]

        # Regimes are defined
        elif self.config.regimes is not None: regimes = self.config.regimes[:]

        # Regimes are not defined
        else: regimes = original_regime_names[:]

        # Ignore certain regimes?
        if self.config.no_uv: regimes = sequences.removed_item(regimes, uv_name)
        if self.config.no_optical: regimes = sequences.removed_item(regimes, optical_name)
        if self.config.no_nir: regimes = sequences.removed_item(regimes, nir_name)
        if self.config.no_mir: regimes = sequences.removed_item(regimes, mir_name)
        if self.config.no_fir: regimes = sequences.removed_item(regimes, fir_name)
        if self.config.no_submm_microwave: regimes = sequences.removed_item(regimes, submm_microwave_name)

        # Return the regimes
        return regimes

    # -----------------------------------------------------------------

    def get_regimes_physical(self):

        """
        This function ...
        :return:
        """

        # Only ionizing
        if self.config.only_ionizing:

            if self.config.no_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")
            regimes = [ionizing_name]

        # Only young
        elif self.config.only_young:

            if self.config.no_young: raise ValueError("Error")
            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")
            regimes = [young_name]

        # Only evolved
        elif self.config.only_evolved:

            if self.config.no_evolved: raise ValueError("Error")
            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")
            regimes = [evolved_name]

        # Only mix
        elif self.config.only_mix:

            if self.config.no_mix: raise ValueError("Error")
            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")
            regimes = [mix_name]

        # Only aromatic
        elif self.config.only_aromatic:

            if self.config.no_aromatic: raise ValueError("Error")
            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_thermal: raise ValueError("Error")
            regimes = [aromatic_name]

        # Only thermal
        elif self.config.only_thermal:

            if self.config.no_thermal: raise ValueError("Error")
            if self.config.only_ionizing: raise ValueError("Error")
            if self.config.only_young: raise ValueError("Error")
            if self.config.only_evolved: raise ValueError("Error")
            if self.config.only_mix: raise ValueError("Error")
            if self.config.only_aromatic: raise ValueError("Error")
            regimes = [thermal_name]

        # Regimes are defined
        elif self.config.regimes is not None: regimes = self.config.regimes[:]

        # Regimes are not defined
        else: regimes = physical_regime_names[:]

        # Ignore certain regimes?
        if self.config.no_ionizing: regimes = sequences.removed_item(regimes, ionizing_name)
        if self.config.no_young: regimes = sequences.removed_item(regimes, young_name)
        if self.config.no_evolved: regimes = sequences.removed_item(regimes, evolved_name)
        if self.config.no_mix: regimes = sequences.removed_item(regimes, mix_name)
        if self.config.no_aromatic: regimes = sequences.removed_item(regimes, aromatic_name)
        if self.config.no_thermal: regimes = sequences.removed_item(regimes, thermal_name)

        # Return the regimes
        return regimes

    # -----------------------------------------------------------------

    @property
    def include_uv(self):

        """
        This function ...
        :return:
        """

        return uv_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def uv_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_uv: return self.config.uv
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_optical(self):

        """
        This function ...
        :return:
        """

        return optical_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def optical_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_optical: return self.config.optical
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_nir(self):

        """
        This function ...
        :return:
        """

        return nir_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def nir_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_nir: return self.config.nir
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_mir(self):

        """
        This function ...
        :return:
        """

        return mir_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def mir_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_mir: return self.config.mir
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_fir(self):

        """
        This function ...
        :return:
        """

        return fir_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def fir_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_fir: return self.config.fir
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_submm_microwave(self):

        """
        This function ...
        :return:
        """

        return submm_microwave_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def submm_microwave_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_submm_microwave: return self.config.submm_microwave
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_ionizing(self):

        """
        This function ...
        :return:
        """

        return ionizing_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_ionizing: return self.config.ionizing
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_young(self):

        """
        This function ...
        :return:
        """

        return young_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def young_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_young: return self.config.young
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_evolved(self):

        """
        This function ...
        :return:
        """

        return evolved_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def evolved_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_evolved: return self.config.evolved
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_mix(self):

        """
        This function ...
        :return:
        """

        return mix_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def mix_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_mix: return self.config.mix
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_aromatic(self):

        """
        This function ...
        :return:
        """

        return aromatic_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def aromatic_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_aromatic: return self.config.aromatic
        else: return 0.

    # -----------------------------------------------------------------

    @property
    def include_thermal(self):

        """
        This function ...
        :return:
        """

        return thermal_name in self.regimes

    # -----------------------------------------------------------------

    @lazyproperty
    def thermal_weight(self):

        """
        This function ...
        :return:
        """

        if self.include_thermal: return self.config.thermal
        else: return 0.

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Get the weights
        if self.config.physical: weights = self.calculate_weights_physical()
        else: weights = self.calculate_weights_original()

        # Add to weights table
        for fltr in weights: self.table.add_point(fltr, weights[fltr])

    # -----------------------------------------------------------------

    def calculate_weights_original(self):

        """
        This function ...
        :return:
        """

        # Get the weights
        return calculate_weights_filters(self.filters, uv=self.uv_weight, optical=self.optical_weight,
                                            nir=self.nir_weight, mir=self.mir_weight, fir=self.fir_weight,
                                            submm_microwave=self.submm_microwave_weight)

    # -----------------------------------------------------------------

    def calculate_weights_physical(self):

        """
        This function ...
        :return:
        """

        # Get the weights
        return calculate_weights_filters_physical(self.filters, ionizing=self.ionizing_weight, young=self.young_weight,
                                                  evolved=self.evolved_weight, mix=self.mix_weight, aromatic=self.aromatic_weight,
                                                  thermal=self.thermal_weight)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the weights table
        self.write_table()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        path = self.output_path_file("weights.dat")

        # Inform the user
        log.info("Writing the weights table to '" + path + "' ...")

        # Write
        self.table.saveto(path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the regime weights
        self.show_regimes()

        # Table
        self.show_table()

    # -----------------------------------------------------------------

    def show_regimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the regime weights ...")

        # Physical
        if self.config.physical: self.show_regimes_physical()

        # Original
        else: self.show_regimes_original()

    # -----------------------------------------------------------------

    def show_regimes_physical(self):

        """
        This function ...
        :return:
        """

        # Show weights per regime
        print("")
        print(fmt.underlined + fmt.yellow + "REGIMES:" + fmt.reset)
        print("")

        # Ionizing
        if self.include_ionizing: print(fmt.green + " - " + fmt.bold + "Ionizing: " + fmt.reset_bold + str(self.ionizing_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Ionizing: " + fmt.reset_bold + "--" + fmt.reset)

        # Young
        if self.include_young: print(fmt.green + " - " + fmt.bold + "Young: " + fmt.reset_bold + str(self.young_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Young: " + fmt.reset_bold + "--" + fmt.reset)

        # Evolved
        if self.include_evolved: print(fmt.green + " - " + fmt.bold + "Evolved: " + fmt.reset_bold + str(self.evolved_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Evolved: " + fmt.reset_bold + "--" + fmt.reset)

        # Mix
        if self.include_mix: print(fmt.green + " - " + fmt.bold + "Mix: " + fmt.reset_bold + str(self.mix_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Mix: " + fmt.reset_bold + "--" + fmt.reset)

        # Aromatic
        if self.include_aromatic: print(fmt.green + " - " + fmt.bold + "Aromatic: " + fmt.reset_bold + str(self.aromatic_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Aromatic: " + fmt.reset_bold + "--" + fmt.reset)

        # Thermal
        if self.include_thermal: print(fmt.green + " - " + fmt.bold + "Thermal: " + fmt.reset_bold + str(self.thermal_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Thermal: " + fmt.reset_bold + "--" + fmt.reset)

        print("")

    # -----------------------------------------------------------------

    def show_regimes_original(self):

        """
        This function ...
        :return:
        """

        # Show weight per regime
        print("")
        print(fmt.underlined + fmt.yellow + "REGIMES:" + fmt.reset)
        print("")

        # UV
        if self.include_uv: print(fmt.green + " - " + fmt.bold + "UV: " + fmt.reset_bold + str(self.uv_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "UV: " + fmt.reset_bold + "--" + fmt.reset)

        # Optical
        if self.include_optical: print(fmt.green + " - " + fmt.bold + "Optical: " + fmt.reset_bold + str(self.optical_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Optical: " + fmt.reset_bold + "--" + fmt.reset)

        # NIR
        if self.include_nir: print(fmt.green + " - " + fmt.bold + "NIR: " + fmt.reset_bold + str(self.nir_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "NIR: " + fmt.reset_bold + "--" + fmt.reset)

        # MIR
        if self.include_mir: print(fmt.green + " - " + fmt.bold + "MIR: " + fmt.reset_bold + str(self.mir_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "MIR: " + fmt.reset_bold + "--" + fmt.reset)

        # FIR
        if self.include_fir: print(fmt.green + " - " + fmt.bold + "FIR: " + fmt.reset_bold + str(self.fir_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "FIR: " + fmt.reset_bold + "--" + fmt.reset)

        # Submm-microwave
        if self.include_submm_microwave: print(fmt.green + " - " + fmt.bold + "Submm-microwave: " + fmt.reset_bold + str(self.submm_microwave_weight) + fmt.reset)
        else: print(fmt.red + " - " + fmt.bold + "Submm-microwave: " + fmt.reset_bold + "--" + fmt.reset)

        print("")

    # -----------------------------------------------------------------

    def show_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the weights table ...")

        # Show the weights
        print("")
        print(self.table)

# -----------------------------------------------------------------

def calculate_weights_filters(filters, uv=1, optical=1, nir=1, mir=1, fir=1, submm_microwave=1):

    """
    This function ...
    :param filters:
    :param uv:
    :param optical:
    :param nir:
    :param mir:
    :param fir:
    :param submm_microwave:
    :return:
    """

    # Get bands per regime
    uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands = split_filters_regimes(filters)

    # Get nbands per regime
    nuv = len(uv_bands)
    noptical = len(optical_bands)
    nnir = len(nir_bands)
    nmir = len(mir_bands)
    nfir = len(fir_bands)
    nsubmm_microwave = len(submm_microwave_bands)

    # Determine regime weights
    uv_weight, optical_weight, nir_weight, mir_weight, fir_weight, submm_radio_weight = calculate_weights(nuv, noptical,
                                                                                                          nnir, nmir,
                                                                                                          nfir,
                                                                                                          nsubmm_microwave,
                                                                                                          uv=uv,
                                                                                                          optical=optical,
                                                                                                          nir=nir,
                                                                                                          mir=mir,
                                                                                                          fir=fir,
                                                                                                          submm_microwave=submm_microwave)

    # Initialize dictionary for the weights per filter
    weights = OrderedDict()

    # Loop over the bands in each group and set the weight in the weights table
    for fltr in uv_bands: weights[fltr] = uv_weight
    for fltr in optical_bands: weights[fltr] = optical_weight
    for fltr in nir_bands: weights[fltr] = nir_weight
    for fltr in mir_bands: weights[fltr] = mir_weight
    for fltr in fir_bands: weights[fltr] = fir_weight
    for fltr in submm_microwave_bands: weights[fltr] = submm_radio_weight

    # Return the weights
    return weights

# -----------------------------------------------------------------

def calculate_weights_filters_physical(filters, ionizing=1, young=1, evolved=1, mix=1, aromatic=1, thermal=1):

    """
    This function ...
    :param filters:
    :param ionizing:
    :param young:
    :param evolved:
    :param mix:
    :param aromatic:
    :param thermal:
    :return:
    """

    # Get bands per physical regime
    ionizing_bands, young_bands, evolved_bands, mix_bands, aromatic_bands, thermal_bands = split_filters_regimes_physical(filters)

    # Get nbands per regime
    nionizing = len(ionizing_bands)
    nyoung = len(young_bands)
    nevolved = len(evolved_bands)
    nmix = len(mix_bands)
    naromatic = len(aromatic_bands)
    nthermal = len(thermal_bands)

    # Determine regime weights
    ionizing_weight, young_weight, evolved_weight, mix_weight, aromatic_weight, thermal_weight = calculate_weights_physical(nionizing, nyoung, nevolved, nmix, naromatic, nthermal, ionizing, young, evolved, mix, aromatic, thermal)

    # Initialize dictionary for the weights per filter
    weights = OrderedDict()

    # Loop over the bands in each group and set the weight in the weights table
    for fltr in ionizing_bands: weights[fltr] = ionizing_weight
    for fltr in young_bands: weights[fltr] = young_weight
    for fltr in evolved_bands: weights[fltr] = evolved_weight
    for fltr in mix_bands: weights[fltr] = mix_weight
    for fltr in aromatic_bands: weights[fltr] = aromatic_weight
    for fltr in thermal_bands: weights[fltr] = thermal_weight

    # Return the weights
    return weights

# -----------------------------------------------------------------

def split_filters_regimes(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    # Initialize lists to contain the filters of the different wavelength ranges
    uv_bands = []
    optical_bands = []
    nir_bands = []
    mir_bands = []
    fir_bands = []
    submm_microwave_bands = []

    # Loop over the filters
    for fltr in filters:

        # Get the central wavelength
        wavelength = fltr.wavelength

        # Get a string identifying which portion of the wavelength spectrum this wavelength belongs to
        spectrum = wavelengths.name_in_spectrum(wavelength)

        # Determine to which group
        if spectrum[0] == "UV":
            uv_bands.append(fltr)
        elif spectrum[0] == "Optical":
            optical_bands.append(fltr)
        elif spectrum[0] == "Optical/IR":
            optical_bands.append(fltr)
        elif spectrum[0] == "IR":
            if spectrum[1] == "NIR": nir_bands.append(fltr)
            elif spectrum[1] == "MIR": mir_bands.append(fltr)
            elif spectrum[1] == "FIR": fir_bands.append(fltr)
            else: raise RuntimeError("Unknown IR range")
        elif spectrum[0] == "Submm": submm_microwave_bands.append(fltr)
        elif spectrum[0] == "Radio" and spectrum[1] == "Microwave": submm_microwave_bands.append(fltr)
        else: raise RuntimeError("Unknown wavelength range: " + str(spectrum))

    # Return
    return uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands

# -----------------------------------------------------------------

def split_filters_regimes_physical(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    # Initialize lists to contain the filters of the different physical regimes
    ionizing_bands = []
    young_bands = []
    evolved_bands = []
    mix_bands = []
    aromatic_bands = []
    thermal_bands = []

    # Loop over the filters
    for fltr in filters:

        # Get the physical regime
        regime = wavelengths.physical_regime_for_filter(fltr)

        # Determine in which group
        if regime == wavelengths.ionizing: ionizing_bands.append(fltr)
        elif regime == wavelengths.young: young_bands.append(fltr)
        elif regime == wavelengths.evolved: evolved_bands.append(fltr)
        elif regime == wavelengths.mix: mix_bands.append(fltr)
        elif regime == wavelengths.aromatic: aromatic_bands.append(fltr)
        elif regime == wavelengths.thermal: thermal_bands.append(fltr)
        else: raise RuntimeError("Unknown physical regime: " + str(regime))

    # Return the filters
    return ionizing_bands, young_bands, evolved_bands, mix_bands, aromatic_bands, thermal_bands

# -----------------------------------------------------------------

def get_nbands_per_regime(filters):

    """
    This function ...
    :param filters:
    :return:
    """

    uv_bands, optical_bands, nir_bands, mir_bands, fir_bands, submm_microwave_bands = split_filters_regimes(filters)
    return len(uv_bands), len(optical_bands), len(nir_bands), len(mir_bands), len(fir_bands), len(submm_microwave_bands)

# -----------------------------------------------------------------

def get_nbands_per_regime_physical(filters):

    """
    This fuction ...
    :param filters:
    :return:
    """

    ionizing_bands, young_bands, evolved_bands, mix_bands, aromatic_bands, thermal_bands = split_filters_regimes_physical(filters)
    return len(ionizing_bands), len(young_bands), len(evolved_bands), len(mix_bands), len(aromatic_bands), len(thermal_bands)

# -----------------------------------------------------------------

def calculate_weights(nuv, noptical, nnir, nmir, nfir, nsubmm_microwave, uv=1, optical=1, nir=1, mir=1, fir=1, submm_microwave=1):

    """
    This function ...
    :param nuv:
    :param noptical:
    :param nnir:
    :param nmir:
    :param nfir:
    :param nsubmm_microwave:
    :param uv:
    :param optical:
    :param nir:
    :param mir:
    :param fir:
    :param submm_microwave:
    :return:
    """

    # Initialize the number of groups
    number_of_groups = 0

    # Check which groups are present
    has_uv = nuv > 0
    has_optical = noptical > 0
    has_nir = nnir > 0
    has_mir = nmir > 0
    has_fir = nfir > 0
    has_submm_microwave = nsubmm_microwave > 0

    # Check how many groups
    if has_uv: number_of_groups += 1
    if has_optical: number_of_groups += 1
    if has_nir: number_of_groups += 1
    if has_mir: number_of_groups += 1
    if has_fir: number_of_groups += 1
    if has_submm_microwave: number_of_groups += 1
    nall_groups = 6

    # Determine total number of data points
    number_of_data_points = nuv + noptical + nnir + nmir + nfir + nsubmm_microwave

    # Determine normalizations
    total_normalization = uv + optical + nir + mir + fir + submm_microwave
    uv = float(uv) / total_normalization * nall_groups
    optical = float(optical) / total_normalization * nall_groups
    nir = float(nir) / total_normalization * nall_groups
    mir = float(mir) / total_normalization * nall_groups
    fir = float(fir) / total_normalization * nall_groups
    submm = float(submm_microwave) / total_normalization * nall_groups

    # Determine the weight for each group of filters
    uv_weight = uv / (nuv * number_of_groups) * number_of_data_points if has_uv else 0.0
    optical_weight = optical / (noptical * number_of_groups) * number_of_data_points if has_optical else 0.0
    nir_weight = nir / (nnir * number_of_groups) * number_of_data_points if has_nir else 0.0
    mir_weight = mir / (nmir * number_of_groups) * number_of_data_points if has_mir else 0.0
    fir_weight = fir / (nfir * number_of_groups) * number_of_data_points if has_fir else 0.0
    submm_microwave_weight = submm / (nsubmm_microwave * number_of_groups) * number_of_data_points if has_submm_microwave else 0.0

    # Debugging
    if has_uv: log.debug("UV: number of bands = " + str(nuv) + ", weight = " + str(uv_weight))
    if has_optical: log.debug("Optical: number of bands = " + str(noptical) + ", weight = " + str(optical_weight))
    if has_nir: log.debug("NIR: number of bands = " + str(nnir) + ", weight = " + str(nir_weight))
    if has_mir: log.debug("MIR: number of bands = " + str(nmir) + ", weight = " + str(mir_weight))
    if has_fir: log.debug("FIR: number of bands = " + str(nfir) + ", weight = " + str(fir_weight))
    if has_submm_microwave: log.debug("Submm/microwave: number of bands = " + str(nsubmm_microwave) + ", weight = " + str(submm_microwave_weight))

    # Return the weights
    return uv_weight, optical_weight, nir_weight, mir_weight, fir_weight, submm_microwave_weight

# -----------------------------------------------------------------

def calculate_weights_physical(nionizing, nyoung, nevolved, nmix, naromatic, nthermal, ionizing=1, young=1, evolved=1,
                               mix=1, aromatic=1, thermal=1):
    """
    This function ...
    :param nionizing:
    :param nyoung:
    :param nevolved:
    :param nmix:
    :param naromatic:
    :param nthermal:
    :param ionizing:
    :param young:
    :param evolved:
    :param mix:
    :param aromatic:
    :param thermal:
    :return:
    """

    # Initialize the number of groups
    number_of_groups = 0

    # Check which groups are present
    has_ionizing = nionizing > 0
    has_young = nyoung > 0
    has_evolved = nevolved > 0
    has_mix = nmix > 0
    has_aromatic = naromatic > 0
    has_thermal = nthermal > 0

    # Check how many groups
    if has_ionizing: number_of_groups += 1
    if has_young: number_of_groups += 1
    if has_evolved: number_of_groups += 1
    if has_mix: number_of_groups += 1
    if has_aromatic: number_of_groups += 1
    if has_thermal: number_of_groups += 1
    nall_groups = 6

    # Determine total number of data points
    number_of_data_points = nionizing + nyoung + nevolved + nmix + naromatic + nthermal

    # Determine normalizations
    total_normalization = ionizing + young + evolved + mix + aromatic + thermal
    ionizing = float(ionizing) / total_normalization * nall_groups
    young = float(young) / total_normalization * nall_groups
    evolved = float(evolved) / total_normalization * nall_groups
    mix = float(mix) / total_normalization * nall_groups
    aromatic = float(aromatic) / total_normalization * nall_groups
    thermal = float(thermal) / total_normalization * nall_groups

    # Determine the weight for each group of filters
    ionizing_weight = ionizing / (nionizing * number_of_groups) * number_of_data_points if has_ionizing else 0.0
    young_weight = young / (nyoung * number_of_groups) * number_of_data_points if has_young else 0.0
    evolved_weight = evolved / (nevolved * number_of_groups) * number_of_data_points if has_evolved else 0.0
    mix_weight = mix / (nmix * number_of_groups) * number_of_data_points if has_mix else 0.0
    aromatic_weight = aromatic / (naromatic * number_of_groups) * number_of_data_points if has_aromatic else 0.0
    thermal_weight = thermal / (nthermal * number_of_groups) * number_of_data_points if has_thermal else 0.0

    # Debugging
    if has_ionizing: log.debug("Ionizing: number of bands = " + str(nionizing) + ", weight = " + str(ionizing_weight))
    if has_young: log.debug("Young: number of bands = " + str(nyoung) + ", weight = " + str(young_weight))
    if has_evolved: log.debug("Evolved: number of bands = " + str(nevolved) + ", weight = " + str(evolved_weight))
    if has_mix: log.debug("Mix: number of bands = " + str(nmix) + ", weight = " + str(mix_weight))
    if has_aromatic: log.debug("Aromatic: number of bands = " + str(naromatic) + ", weight = " + str(aromatic_weight))
    if has_thermal: log.debug("Thermal: number of bands = " + str(nthermal) + ", weight = " + str(thermal_weight))

    # Return the weights
    return ionizing_weight, young_weight, evolved_weight, mix_weight, aromatic_weight, thermal_weight

# -----------------------------------------------------------------
