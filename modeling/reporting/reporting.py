#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.reporting.reporting Contains the Reporter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import ReportingComponent
from ...core.basics.log import log
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...magic.core.image import Image

# -----------------------------------------------------------------

steps = ["data", "preparation_initialization", "preparation", "decomposition", "photometry", "maps", "input_initialization", "exploration", "fitting", "analysis"]

# -----------------------------------------------------------------

class Reporter(ReportingComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(Reporter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make a report for the data
        if self.config.step == steps[0]: self.report_data()

        # 3. Make a report for the preparation initialization step
        elif self.config.step == steps[1]: self.report_preparation_initialization()

        # 4. Make a report for the data preparation step
        elif self.config.steps == steps[2]: self.report_preparation()

        # 5. Make a report for the decomposition step
        elif self.config.steps == steps[3]: self.report_decomposition()

        # 6. Make a report for the photometry step
        elif self.config.steps == steps[4]: self.report_photometry()

        # 7. Make a report for the map making step
        elif self.config.steps == steps[5]: self.report_map_making()

        # 8. Make a report for the input initialization step
        elif self.config.steps == steps[6]: self.report_input_initialization()

        # 9. Make a report for the parameter exploration step
        elif self.config.steps == steps[7]: self.report_exploration()

        # 10. Make a report for the SED fitting step
        elif self.config.steps == steps[8]: self.report_fitting()

        # 11. Make a report for the analysis step
        elif self.config.steps == steps[9]: self.report_analysis()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Reporter, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def report_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a report for the data ...")

        # Initialize the columns of the table
        names_column = []
        paths_column = []
        observatory_column = []
        instrument_column = []
        band_column = []
        unit_column = []
        pixelscale_column = []
        has_errors_column = []
        prep_names_column = []
        names = ["Image name", "Image path", "Observatory", "Instrument", "Band", "Unit", "Pixelscale", "Has errors",
                 "Preparation name"]

        # Loop over all subdirectories of the data directory
        for path, origin in fs.directories_in_path(self.data_images_path, returns=["path", "name"]):

            # Loop over all FITS files found in the current subdirectory
            for image_path, image_name in fs.files_in_path(path, extension="fits", not_contains="Error", returns=["path", "name"]):

                # Open the image
                image = Image.from_file(image_path)

                # Check if an error map is present (as one of the image frames or as a seperate FITS file)
                has_errors = "errors" in image.frames or fs.is_file(fs.join(path, image_name + "_Error.fits"))

                if image.filter is not None:
                    observatory = image.filter.observatory
                    instrument = image.filter.instrument
                    band = image.filter.band
                    prep_name = instrument + " " + band
                else:
                    observatory = None
                    instrument = None
                    band = None
                    prep_name = image_name

                names_column.append(image_name)
                paths_column.append(image_path)
                observatory_column.append(observatory)
                instrument_column.append(instrument)
                band_column.append(band)
                unit_column.append(str(image.unit))
                pixelscale_column.append(image.average_pixelscale.to("arcsec").value)
                has_errors_column.append(has_errors)
                prep_names_column.append(prep_name)

        # Create the table
        data = [names_column, paths_column, observatory_column, instrument_column, band_column, unit_column,
                pixelscale_column, has_errors_column, prep_names_column]
        table = tables.new(data, names)

        # Debugging
        log.debug("Writing the data report to '" + self.data_report_path + "'...")

        # Save the table
        tables.write(table, self.data_report_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def report_preparation_initialization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a report for the preparation initialization step ...")

        # Initialize the columns of the table
        names_column = []
        initialized_column = []
        errors_column = []
        galaxy_column = []
        star_column = []
        saturation_column = []
        other_column = []
        segments_column = []
        statistics_column = []
        fwhm_column = []

        names = ["Image name", "Initialized", "Errors", "Galaxy region", "Star region", "Saturation region",
                 "Other sources region", "Segmentation maps", "Statistics file", "FWHM"]

        # Loop over all subdirectories of the preparation directory
        for path, name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Determine the path to the initialized image
            image_path = fs.join(path, "initialized.fits")
            initialized = fs.is_file(image_path)

            # Determine the path to the sources directory
            sources_path = fs.join(path, "sources")
            sources_present = fs.is_directory(sources_path)

            # Check the presence of the region files
            has_galaxy_region = fs.is_file(fs.join(sources_path, "galaxies.reg")) if sources_present else False
            has_star_region = fs.is_file(fs.join(sources_path, "stars.reg")) if sources_present else False
            has_saturation_region = fs.is_file(fs.join(sources_path, "saturation.reg")) if sources_present else False
            has_other_region = fs.is_file(fs.join(sources_path, "other_sources.reg")) if sources_present else False

            # Check the presence of the segmentation image
            has_segments = fs.is_file(fs.join(sources_path, "segments.fits")) if sources_present else False

            # Check the presence of the statistics file
            has_statistics = fs.is_file(fs.join(sources_path, "statistics.dat")) if sources_present else False

            # Open the image
            image = Image.from_file(image_path)

            # Check the presence of an error frame
            has_errors = "errors" in image.frames.keys()

            # Get the FWHM
            fwhm = image.fwhm.to("arcsec").value

            # Fill in the columns
            names_column.append(name)
            initialized_column.append(initialized and sources_present)
            errors_column.append(has_errors)
            galaxy_column.append(has_galaxy_region)
            star_column.append(has_star_region)
            saturation_column.append(has_saturation_region)
            other_column.append(has_other_region)
            segments_column.append(has_segments)
            statistics_column.append(has_statistics)
            fwhm_column.append(fwhm)

        # Create the table
        data = [names_column, initialized_column, errors_column, galaxy_column, star_column, saturation_column,
                other_column, segments_column, statistics_column, fwhm_column]
        table = tables.new(data, names)

        # Debugging
        log.debug("Writing the preparation initialization report to '" + self.preparation_initialization_report_path + "'...")

        # Save the table
        tables.write(table, self.preparation_initialization_report_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def report_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a report from the data preparation step ...")

        # Initialize the columns of the table
        names_column = []
        extracted_column = []
        corrected_column = []
        converted_column = []
        convolved_column = []
        rebinned_column = []
        subtracted_column = []
        result_column = []
        unit_column = []
        pixelscale_column = []
        fwhm_column = []
        errors_column = []
        sky_column = []
        sources_column = []

        names = ["Image name", "Sources extracted", "Corrected for extinction", "Unit converted", "Convolved",
                 "Rebinned", "Sky subtracted", "Result", "Unit", "Pixelscale", "FWHM", "Has errors", "Has sky",
                 "Has sources mask"]

        # Loop over all subdirectories of the preparation directory
        for path, name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Debugging
            log.debug("Checking " + name + " image ...")

            # Determine the path to the extracted image
            extracted_path = fs.join(path, "extracted.fits")
            has_extracted = fs.is_file(extracted_path)

            # Determine the path to the extinction-corrected image
            corrected_path = fs.join(path, "corrected_for_extinction.fits")
            has_corrected = fs.is_file(corrected_path)

            # Determine the path to the unit-converted image
            converted_path = fs.join(path, "converted_unit.fits")
            has_converted = fs.is_file(converted_path)

            # Determine the path to the convolved image
            convolved_path = fs.join(path, "convolved.fits")
            has_convolved = fs.is_file(convolved_path)

            # Determine the path to the rebinned image
            rebinned_path = fs.join(path, "rebinned.fits")
            has_rebinned = fs.is_file(rebinned_path)

            # Determine the path to the sky-subtracted image
            subtracted_path = fs.join(path, "subtracted.fits")
            has_subtracted = fs.is_file(subtracted_path)

            # Determine the path to the prepared image
            result_path = fs.join(path, "result.fits")
            has_result = fs.is_file(result_path)

            # If the prepared image is present, open it and get some properties
            if has_result:

                result = Image.from_file(result_path)

                unit = str(result.unit)
                pixelscale = result.average_pixelscale.to("arcsec").value
                fwhm = result.fwhm.to("arcsec").value if result.fwhm is not None else None
                has_errors = "errors" in result.frames.keys() and not result.frames["errors"].all_zero
                has_sky = "sky" in result.frames.keys() and not result.frames["sky"].all_zero
                has_sources = "sources" in result.masks.keys()

            else:

                unit = None
                pixelscale = None
                fwhm = None
                has_errors = None
                has_sky = None
                has_sources = None

            # Fill in the columns
            names_column.append(name)
            extracted_column.append(has_extracted)
            corrected_column.append(has_corrected)
            converted_column.append(has_converted)
            convolved_column.append(has_convolved)
            rebinned_column.append(has_rebinned)
            subtracted_column.append(has_subtracted)
            result_column.append(has_result)
            unit_column.append(unit)
            pixelscale_column.append(pixelscale)
            fwhm_column.append(fwhm)
            errors_column.append(has_errors)
            sky_column.append(has_sky)
            sources_column.append(has_sources)

        # Create the table
        data = [names_column, extracted_column, corrected_column, converted_column, convolved_column, rebinned_column,
                subtracted_column, result_column, unit_column, pixelscale_column, fwhm_column, errors_column,
                sky_column, sources_column]
        table = tables.new(data, names)

        # Debugging
        log.info("Writing the preparation report to '" + self.preparation_report_path + "' ...")

        # Save the table
        tables.write(table, self.preparation_report_path, format="ascii.ecsv")

    # -----------------------------------------------------------------

    def report_decomposition(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_photometry(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_map_making(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_input_initialization(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_exploration(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_fitting(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def report_analysis(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
