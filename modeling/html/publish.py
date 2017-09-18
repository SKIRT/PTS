#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.publish Contains the Publisher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, table_class
from .all import AllPagesGenerator
from ...core.remote.host import load_host
from ...core.remote.mounter import RemoteMounter
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ..core.environment import index_filename
from ...core.tools import browser

# -----------------------------------------------------------------

modeling_name = "modelling"
images_name = "images"

# -----------------------------------------------------------------

base_url = "http://users.ugent.be/~sjversto"

# -----------------------------------------------------------------

class Publisher(HTMLPageComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Publisher, self).__init__(*args, **kwargs)

        # The mounter and mount path
        self.mounter = None
        self.mount_path = None

    # -----------------------------------------------------------------

    @lazyproperty
    def host(self):

        """
        This function ...
        :return:
        """

        return load_host(self.config.host_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_modeling_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.mount_path, modeling_name)
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_galaxy_path(self):

        """
        This function ..
        :return:
        """

        path = fs.join(self.remote_modeling_path, self.galaxy_name)
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_galaxy_images_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.remote_galaxy_path, images_name)
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_truncation_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.remote_galaxy_path, "truncation")
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_maps_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.remote_galaxy_path, "maps")
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_index_url(self):

        """
        This function ...
        :return:
        """

        return fs.join(base_url, modeling_name, self.galaxy_name, index_filename)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Setup
        self.setup(**kwargs)

        # 2. Generate the pages
        if self.config.regenerate or self.config.generate: self.generate()

        # 3. Upload
        self.upload()

        # 4. Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Publisher, self).setup(**kwargs)

        # Create the mounter and mount the host
        self.mounter = RemoteMounter()
        self.mount_path = self.mounter.mount(self.host)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the pages ...")

        # Create the generator
        generator = AllPagesGenerator()

        # Set modeling path
        generator.config.path = self.config.path

        # Set settings
        generator.config.regenerate = self.config.regenerate
        generator.config.replot = self.config.replot
        generator.config.details = self.config.details

        # Run the generator
        generator.run()

    # -----------------------------------------------------------------

    def upload(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading ...")

        # Upload pages
        self.upload_pages()

        # Upload images
        self.upload_images()

    # -----------------------------------------------------------------

    def upload_pages(self):

        """
        This function ...
        :return:
        """

        # 1. Upload the pages
        self.upload_main_pages()

        # 2. Upload truncation pages
        if self.config.details: self.upload_truncation_pages()

        # 3. Upload maps pages
        if self.config.details: self.upload_maps_pages()

    # -----------------------------------------------------------------

    def upload_main_pages(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the pages ...")

        # Index page
        if self.has_index_page: self.upload_index_page()

        # Status page
        if self.has_status_page: self.upload_status_page()

        # Data page
        if self.has_data_page: self.upload_data_page()

        # Preparation page
        if self.has_preparation_page: self.upload_preparation_page()

        # Components page
        if self.has_components_page: self.upload_components_page()

        # Photometry page
        if self.has_photometry_page: self.upload_photometry_page()

        # Maps page
        if self.has_maps_page: self.upload_maps_page()

        # Model page
        if self.has_model_page: self.upload_model_page()

        # Fitting page
        if self.has_fitting_page: self.upload_fitting_page()

        # Datacubes page
        if self.has_datacubes_page: self.upload_datacubes_page()

        # Fluxes page
        if self.has_fluxes_page: self.upload_fluxes_page()

        # Images
        if self.has_images_page: self.upload_images_page()

        # Colours page
        if self.has_colours_page: self.upload_colours_page()

        # Attenuation page
        if self.has_attenuation_page: self.upload_attenuation_page()

        # Heating page
        if self.has_heating_page: self.upload_heating_page()

    # -----------------------------------------------------------------

    def upload_index_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the index page ...")

        # Update
        updated = fs.update_file_in(self.index_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_status_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the status page ...")

        # Update
        updated = fs.update_file_in(self.status_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_data_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the data page ...")

        # Update
        updated = fs.update_file_in(self.data_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_preparation_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the preparation page ...")

        # Update
        updated = fs.update_file_in(self.preparation_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_components_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the components page ...")

        # Update
        updated = fs.update_file_in(self.components_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_photometry_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the photometry page ...")

        # Update
        updated = fs.update_file_in(self.photometry_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_maps_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the maps page ...")

        # Update
        updated = fs.update_file_in(self.maps_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_model_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the model page ...")

        # Update
        updated = fs.update_file_in(self.model_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_fitting_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the fitting page ...")

        # Update
        updated = fs.update_file_in(self.fitting_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_datacubes_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the datacubes page ...")

        # Update
        updated = fs.update_file_in(self.datacubes_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_fluxes_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the fluxes page ...")

        # Update
        updated = fs.update_file_in(self.fluxes_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_images_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the images page ...")

        # Update
        updated = fs.update_file_in(self.images_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Succesfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_colours_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the colours page ...")

        # Update
        updated = fs.update_file_in(self.colours_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_attenuation_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading attenuation page ...")

        # Update
        updated = fs.update_file_in(self.attenuation_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_heating_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the heating page ...")

        # Update
        updated = fs.update_file_in(self.heating_page_path, self.remote_galaxy_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_truncation_pages(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the truncation page ...")

        # Truncation ellispe
        if self.has_truncation_ellipse_page: self.upload_truncation_ellipse_page()

        # Significance levels
        if self.has_truncation_significance_page: self.upload_significance_levels_page()

    # -----------------------------------------------------------------

    def upload_truncation_ellipse_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading truncation ellipse page ...")

        # Update
        updated = fs.update_file_in(self.truncation_ellipse_page_path, self.remote_truncation_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_significance_levels_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading significance levels page ...")

        # Update
        updated = fs.update_file_in(self.truncation_significance_page_path, self.remote_truncation_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_maps_pages(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Uploading the maps pages ...")

        # All maps page
        if self.has_all_maps_page: self.upload_all_maps_page()

        # Maps summary page
        if self.has_maps_summary_page: self.upload_maps_summary_page()

        # Old stellar maps
        if self.has_old_maps_page: self.upload_old_maps_page()

        # Young stellar maps
        if self.has_young_maps_page: self.upload_young_maps_page()

        # Ionizing stellar maps
        if self.has_ionizing_maps_page: self.upload_ionizing_maps_page()

        # Dust maps
        if self.has_dust_maps_page: self.upload_dust_maps_page()

        # Clip maps
        if self.has_clip_maps_page: self.upload_clip_maps_page()

        # Maps selection
        #self.upload_maps_selection_page()

    # -----------------------------------------------------------------

    def upload_all_maps_page(self):

        """
        Thisf unction ...
        :return:
        """

        # Infomr the user
        log.info("Uploading all maps page ...")

        # Update
        updated = fs.update_file_in(self.all_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_maps_summary_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading maps summary page ...")

        # Update
        updated = fs.update_file_in(self.maps_summary_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_old_maps_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading old maps page ...")

        # Update
        updated = fs.update_file_in(self.old_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_young_maps_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading young maps page ...")

        # Update
        updated = fs.update_file_in(self.young_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_ionizing_maps_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading ionizing maps page ...")

        # Update
        updated = fs.update_file_in(self.ionizing_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_dust_maps_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading dust maps page ...")

        # Update
        updated = fs.update_file_in(self.dust_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_clip_maps_page(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Uploading clip maps page ...")

        # Update
        updated = fs.update_file_in(self.clip_maps_page_path, self.remote_maps_path, create=True, report=log.is_debug())

        # Report
        if updated: log.success("Successfully uploaded the page")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def upload_maps_selection_page(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Uploading maps selection page ...")

    # -----------------------------------------------------------------

    def upload_images(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Uploading the images ...")

        # Main pages images
        self.upload_main_images()

    # -----------------------------------------------------------------

    def upload_main_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Uploading the images of the main pages ...")

        # Synchronize
        updated = fs.update_directory(self.images_path, self.remote_galaxy_images_path, create=True, report=log.is_debug())

        if updated: log.success("Succesfully uploaded the images")
        else: log.info("Already up-to-date")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Open in browser
        browser.open_url(self.galaxy_index_url)

    # -----------------------------------------------------------------

    @property
    def page_path(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
