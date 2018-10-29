#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.images Contains the ImagesAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.misc.images import ObservedImageMaker
from ...core.units.parsing import parse_unit as u
from ...core.tools.utils import lazyproperty
from ...core.tools import sequences
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

earth_name = "earth"

# -----------------------------------------------------------------

class ImagesAnalyser(AnalysisRunComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        Thisn function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AnalysisRunComponent, self).__init__(*args, **kwargs)

        # The mock observed images
        self.images = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Make the images
        self.make_images()

        # Writing
        self.write()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisRunComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------
    # ANALYSIS RUN
    # -----------------------------------------------------------------

    @property
    def simulation_output_path(self):
        return self.analysis_run.total_output_path

    # -----------------------------------------------------------------

    @property
    def images_path(self):
        return self.analysis_run.images_path

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):
        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------
    # FILTERS
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range(self):
        filters = []
        for fltr in self.observed_filters:
            if not self.wavelength_grid.covers(fltr.wavelength):
                log.warning("The '" + str(fltr) + "' filter is not covered by the wavelength range of the analysis simulations: not making observations for this filter")
                continue
            filters.append(fltr)
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range_without_iras(self):
        return [fltr for fltr in self.observed_filters_in_range if fltr not in self.iras_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):
        return sequences.sorted_by_attribute(self.observed_filters_in_range_without_iras, "wavelength")

    # -----------------------------------------------------------------

    @lazyproperty
    def present_image_names(self):
        return fs.files_in_path(self.images_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def present_filters(self):
        return [parse_filter(name) for name in self.present_image_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def image_maker_filters(self):
        return sequences.elements_not_in_other(self.filters, self.present_filters)

    # -----------------------------------------------------------------
    # INPUT
    # -----------------------------------------------------------------

    @lazyproperty
    def image_maker_input(self):

        """
        This function ...
        :return:
        """

        # Set input
        input_dict = dict()

        # The output path of the simulation
        input_dict["simulation_output_path"] = self.simulation_output_path

        # The output path for the images
        input_dict["output_path"] = self.images_path
        input_dict["output_paths_instruments"] = {earth_name: self.images_path}

        # Filters and instruments
        input_dict["filters"] = self.image_maker_filters
        input_dict["instrument_names"] = [earth_name]

        # Set coordinate system of the datacube
        input_dict["wcs_path"] = self.analysis_run.reference_map_path
        input_dict["wcs_instrument"] = earth_name

        # Unit conversion
        input_dict["unit"] = u("Jy")  # self.misc_options.images_unit

        # Convolution
        if self.config.convolve:
            input_dict["auto_psfs"] = True
            # input_dict["kernel_paths"] = self.misc_options.images_kernels
            input_dict["fwhms_dataset"] = self.static_photometry_dataset  # self.misc_options.fwhms_dataset # path or dataset is possible

        # Set dataset for rebinning
        input_dict["rebin_dataset"] = self.static_photometry_dataset  # path or dataset is possible
        input_dict["rebin_instrument"] = earth_name

        # NO SPECTRAL CONVOLUTION FOR CERTAIN IMAGES?
        input_dict["no_spectral_convolution_filters"] = self.planck_filters

        # Return
        return input_dict

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :param filters:
        :return:
        """

        # Inform the user
        log.info("Making the images ...")

        # Create the maker
        maker = ObservedImageMaker()

        # Set options
        maker.config.spectral_convolution = True

        # Write intermediate results
        maker.config.write_intermediate = self.config.intermediate
        maker.config.write_kernels = self.config.kernels

        # Set number of processes to one
        maker.config.nprocesses_local = 1

        # Settings for convolution
        maker.config.check_wavelengths = True
        maker.config.ignore_bad = True
        maker.config.skip_ignored_bad_convolution = False
        maker.config.skip_ignored_bad_closest = False

        # No plotting
        maker.config.plot = False

        # Group (don't add instrument name prefixes), but write earth instrument output into the main directory (see input)
        maker.config.group = True

        # Don't convolve unless auto_psfs is enabled in input
        maker.config.convolve = False

        # Run
        maker.run(**self.image_maker_input)

        # Update the images
        # Check for if all images were already made
        if earth_name in maker.images: self.images.update(maker.images[earth_name])

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
