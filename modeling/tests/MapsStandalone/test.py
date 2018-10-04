#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.magic.maps.colour.colour import make_map as make_colour_map
from pts.magic.maps.ssfr.colours import make_map as make_ssfr_map
from pts.magic.maps.tir.multi import make_map as make_tir_map_multi
from pts.magic.maps.tir.single import make_map as make_tir_map_single
from pts.magic.maps.youngstars.young import make_map as make_young_stellar_map
from pts.magic.maps.oldstars.disk import make_map as make_old_stellar_map
from pts.magic.maps.attenuation.cortese import make_map as make_fuv_attenuation_map
from pts.magic.maps.attenuation.buat import make_map as make_nuv_or_fuv_attenuation_map
from pts.magic.maps.dust.attenuation import make_map as make_diffuse_dust_map
from pts.magic.maps.dust.hot import make_map as make_hot_dust_map
from pts.magic.maps.ionizingstars.ionizing import make_map as make_ionizing_stellar_map
from pts.core.test.implementation import TestImplementation
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import resolve_name, get_distance, get_inclination, get_center
from pts.core.tools import introspection
from pts.magic.core.list import FrameList
from pts.core.tools import types
from pts.core.filter.filter import parse_filter
from pts.dustpedia.core.properties import DustPediaProperties
from pts.magic.convolution.kernels import has_variable_fwhm, get_fwhm
from pts.magic.convolution.aniano import AnianoKernels
from pts.magic.core.frame import Frame
from pts.magic.services.s4g import get_properties, get_components
from pts.modeling.basics.models import SersicModel3D
from pts.core.simulation.skifile import SkiFile
from pts.core.launch.launcher import SingleImageSKIRTLauncher
from pts.modeling.basics.projection import GalaxyProjection
from pts.modeling.basics.instruments import SimpleInstrument

# -----------------------------------------------------------------

description = "testing the map making for a certain galaxy"

# -----------------------------------------------------------------

maps_commands = ["make_colours_maps", "make_ssfr_maps", "make_tir_maps", "make_attenuation_maps", "make_old_stellar_maps", "make_dust_map", "make_young_stellar_maps", "make_ionizing_stellar_maps"]

# -----------------------------------------------------------------

ssfr_colour = "FUV-H"

# -----------------------------------------------------------------

colour_filename = "FUV-H.fits"
ssfr_filename = "sSFR.fits"
tir_single_filename = "TIR_single.fits"
tir_multi_filename = "TIR_multi.fits"
fuv_attenuation_filename = "FUV_attenuation.fits"
nuv_attenuation_filename = "NUV_attenuation.fits"
old_bulge_filename = "old_bulge.fits"
old_disk_filename = "old_disk.fits"
diffuse_dust_filename = "diffuse_dust.fits"
hot_dust_filename = "hot_dust.fits"
young_stars_filename = "young_stars.fits"
ionizing_stars_filename = "ionizing_stars.fits"

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M77
m77_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "M77")

# Determine the path to the H-alpha image
halpha_path = fs.join(m77_data_path, "Halpha", "NGC1068_Halpha.fits")

# -----------------------------------------------------------------

# The path to the template ski files directory
template_path = fs.join(introspection.pts_dat_dir("modeling"), "ski")

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

class MapsStandaloneTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor
        :param kwargs:
        """

        # Call the constructor of the base class
        super(MapsStandaloneTest, self).__init__(*args, **kwargs)

        # The galaxy distance and inclination
        self.galaxy_center = None
        self.distance = None
        self.inclination = None

        # The data path
        self.data_path = None

        # The DustPedia database
        self.database = None

        # The frames
        self.frames = None

        # FWHMs of images
        self.fwhms = None

        # Paths for maps
        self.maps_path = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the data
        if not self.from_existing_reference: self.get_data()
        else: self.check_reference_data()

        # Load the frames
        self.load_frames()

        # Make colour maps
        self.make_colour_maps()

        # Make sSFR maps
        self.make_ssfr_maps()

        # Make TIR maps
        self.make_tir_maps()

        # Make attenuation maps
        self.make_attenuation_maps()

        # Make maps of old stars
        self.make_old_stars_maps()

        # Make dust maps
        self.make_dust_maps()

        # Make young stars maps
        self.make_young_stars_maps()

        # Make ionizing stars maps
        self.make_ionizing_stars_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsStandaloneTest, self).setup(**kwargs)

        # Create data path
        self.data_path = fs.create_directory_in(self.path, "data")

        # Create maps path
        self.maps_path = fs.create_directory_in(self.path, "maps")

        # Login to the DustPedia database
        self.database = DustPediaDatabase()
        username, password = get_account()
        self.database.login(username, password)

        # Get the galaxy center
        self.galaxy_center = get_center(self.config.galaxy)

        # Get the distance
        self.distance = get_distance(self.config.galaxy)

        # Get the inclination
        self.inclination = get_inclination(self.config.galaxy)

        # Get the DustPedia properties instance
        properties = DustPediaProperties()

        # Get FWHMs
        self.fwhms = properties.fwhms

    # -----------------------------------------------------------------

    @property
    def from_existing_reference(self):

        """
        This function ...
        :return:
        """

        return self.config.reference_path is not None or self.config.reference_test is not None

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the image data ...")

        # Resolve the name
        galaxy_name = resolve_name(self.config.galaxy)

        # Download the images
        self.database.download_images(galaxy_name, self.data_path, error_maps=False, not_instruments=["DSS"])

    # -----------------------------------------------------------------

    def check_reference_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the reference data ...")

        # Determine simulation directory
        if self.config.reference_path is not None: data_path = self.config.reference_path
        elif self.config.reference_test is not None: data_path = fs.join(introspection.pts_tests_dir, self.config.reference_test, "data")
        else: raise ValueError("Reference path and reference test settings are None")

        # Check whether directory exist and not empty
        if not fs.is_directory(data_path): raise ValueError("Directory does not exist: " + data_path)
        if fs.is_empty(data_path): raise ValueError("Empty directory: " + data_path)

        # Remove data directory for this test
        fs.remove_directory(self.data_path)

        # Set data path
        self.data_path = data_path

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the frames ...")

        # Load as frame list
        self.frames = FrameList.from_directory(self.data_path, recursive=True, not_contains="_DSS")

        # Load and add the H-alpha image
        self.frames.append_from_file(halpha_path, "Halpha")

        # Set the distance (to each frame)
        self.frames.distance = self.distance

    # -----------------------------------------------------------------

    def get_frame(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        frame = self.frames[fltr]

        # Set FWHM if necessary
        if frame.fwhm is None:
            if has_variable_fwhm(frame.filter): frame.fwhm = self.fwhms[frame.filter]
            else: frame.fwhm = get_fwhm(fltr)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    def get_wcs(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)

        return self.frames[fltr].wcs

    # -----------------------------------------------------------------

    def get_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine path
        path = fs.join(self.maps_path, name)
        return Frame.from_file(path)

    # -----------------------------------------------------------------

    def make_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making colour maps ...")

        # Make FUV-H
        fuv_h = make_colour_map(self.get_frame("FUV"), self.get_frame("H"))

        # Determine path
        path = fs.join(self.maps_path, colour_filename)

        # Save the map
        fuv_h.saveto(path)

    # -----------------------------------------------------------------

    def make_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making sSFR maps ...")

        # Make FUV-H sSFR (in this case this is trivial)
        ssfr = make_ssfr_map(fuv_h=self.get_map(colour_filename))

        # Determine path
        path = fs.join(self.maps_path, ssfr_filename)

        # Save the map
        ssfr.saveto(path)

    # -----------------------------------------------------------------

    def make_tir_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR maps ...")

        # Make TIR maps based on single bands
        self.make_tir_maps_single()

        # Make TIR maps based on multiple bands
        self.make_tir_maps_multi()

    # -----------------------------------------------------------------

    def make_tir_maps_single(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR maps based on single bands ...")

        # Make TIR map
        tir = make_tir_map_single(self.get_frame("Pacs blue"))

        # Determine path
        path = fs.join(self.maps_path, tir_single_filename)

        # Save
        tir.saveto(path)

    # -----------------------------------------------------------------

    def make_tir_maps_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR maps based on multiple bands ...")

        # Make TIR map
        tir = make_tir_map_multi(self.get_frame("MIPS 24mu"), self.get_frame("Pacs blue"), self.get_frame("Pacs red"))

        # Determine path
        path = fs.join(self.maps_path, tir_multi_filename)

        # Save
        tir.saveto(path)

    # -----------------------------------------------------------------

    def make_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making attenuation maps ...")

        # Make FUV attenuation map
        self.make_attenuation_maps_fuv()

        # Make NUV attenuation map
        self.make_attenuation_maps_nuv()

    # -----------------------------------------------------------------

    def make_attenuation_maps_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the FUV attenuation map ...")

        # Get input
        fuv = self.get_frame("FUV")
        tir = self.get_map(tir_multi_filename)
        ssfr = self.get_map(ssfr_filename)

        # Make FUV attenuation map
        fuv_attenuation = make_fuv_attenuation_map(fuv, tir, ssfr, ssfr_colour)

        # Determine path
        path = fs.join(self.maps_path, fuv_attenuation_filename)

        # Save
        fuv_attenuation.saveto(path)

    # -----------------------------------------------------------------

    def make_attenuation_maps_nuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the NUV attenuation map ...")

        # Get input
        nuv = self.get_frame("NUV")
        tir = self.get_map(tir_multi_filename)

        # Make NUV attenuation map
        nuv_attenuation = make_nuv_or_fuv_attenuation_map(nuv, tir)

        # Detemrine the path
        path = fs.join(self.maps_path, nuv_attenuation_filename)

        # Save
        nuv_attenuation.saveto(path)

    # -----------------------------------------------------------------

    def make_old_stars_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making old stellar maps ...")

        # Make bulge map
        self.make_old_stars_maps_bulge()

        # Make disk map
        self.make_old_stars_maps_disk()

    # -----------------------------------------------------------------

    def make_old_stars_maps_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making map of old stellar bulge ...")

        # Get decomposition properties
        #properties = get_properties(self.config.galaxy)
        components = get_components(self.config.galaxy)
        disk_pa = components["disk"].position_angle

        # Create a Sersic model for the bulge
        bulge_deprojection_method = "tilt"
        print("axial ratio:", components["bulge"].axial_ratio)
        print("inclination:", self.inclination)
        components["bulge"].axial_ratio = 1./components["bulge"].axial_ratio
        bulge = SersicModel3D.from_2d(components["bulge"], self.inclination, disk_pa, azimuth_or_tilt=bulge_deprojection_method)

        # CREATE INSTRUMENT

        # Create the 'earth' projection system
        azimuth = 0.0
        projection = GalaxyProjection.from_wcs(self.get_wcs("I2"), self.galaxy_center, self.distance, self.inclination, azimuth, disk_pa)

        # Create instrument from projection
        instrument = SimpleInstrument.from_projection(projection)

        # Inform the user
        log.debug("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = fs.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        npackages = 1e7

        # Set the number of photon packages
        ski.setpackages(npackages)

        # Set the bulge geometry
        ski.set_stellar_component_geometry(0, bulge)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        ski.add_instrument(instrument_name, instrument)

        # Determine the simulation path
        simulation_path = fs.create_directory_in(self.path, "bulge")

        # Determine the path to the ski file
        ski_path = fs.join(simulation_path, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory and create it
        out_path = fs.create_directory_in(simulation_path, "out")

        # Inform the user
        log.debug("Running the bulge simulation ...")

        # The SKIRT launching environment
        launcher = SingleImageSKIRTLauncher()

        # Load the PSF kernel and prepare
        aniano = AnianoKernels()
        psf = aniano.get_psf("I2")
        psf.prepare_for(self.get_wcs("I2"))

        # Simulate the bulge image
        fluxdensity = components["bulge"].fluxdensity
        bulge_image = launcher.run(ski_path, out_path, self.get_wcs("I2"), fluxdensity, psf, instrument_name=instrument_name, progress_bar=True)

        # Determine path for the bulge map
        path = fs.join(self.maps_path, old_bulge_filename)

        # Save the bulge map
        bulge_image.saveto(path)

    # -----------------------------------------------------------------

    def make_old_stars_maps_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making map of old stellar disk ...")

        # Get input
        total = self.get_frame("I2")
        bulge = self.get_map(old_bulge_filename)

        # Make the map
        old_disk = make_old_stellar_map(total, bulge)

        # Determine the path for the disk map
        path = fs.join(self.maps_path, old_disk_filename)

        # Save the disk map
        old_disk.saveto(path)

    # -----------------------------------------------------------------

    def make_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making dust maps ...")

        # Diffuse dust through attenuation
        self.make_diffuse_dust_maps()

        # Hot dust
        self.make_hot_dust_maps()

    # -----------------------------------------------------------------

    def make_diffuse_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making map of diffuse dust ...")

        # Make map (trivial in this case)
        dust = make_diffuse_dust_map(self.get_map(fuv_attenuation_filename))

        # Detemrine path
        path = fs.join(self.maps_path, diffuse_dust_filename)

        # Save the dust map
        dust.saveto(path)

    # -----------------------------------------------------------------

    def make_hot_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making map of hot dust ...")

        # To be determined for each galaxy
        correction_factor = 0.25

        # Make map
        dust = make_hot_dust_map(self.get_frame("MIPS 24mu"), self.get_map(old_disk_filename), correction_factor)

        # Detemrine path
        path = fs.join(self.maps_path, hot_dust_filename)

        # Save the dust map
        dust.saveto(path)

    # -----------------------------------------------------------------

    def make_young_stars_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps of young stars ...")

        # Get input
        fuv = self.get_frame("FUV")
        fuv_attenuation = self.get_map(fuv_attenuation_filename)
        old_disk = self.get_map(old_disk_filename)
        factor = 0.25 # to be determined for each galaxy

        # Make the map of young stars
        young_stars = make_young_stellar_map(self.get_frame("FUV"), self.get_map(fuv_attenuation_filename), old_disk, factor)

        # Determine the path
        path = fs.join(self.maps_path, young_stars_filename)

        # Save the map
        young_stars.saveto(path)

    # -----------------------------------------------------------------

    def make_ionizing_stars_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps of ionizing stars ...")

        # Make map of ionizing stars
        ionizing_stars = make_ionizing_stellar_map(self.get_frame("Halpha"), self.get_map(hot_dust_filename))

        # Detemrine the path
        path = fs.join(self.maps_path, )

        # Save the map
        ionizing_stars.saveto(path)

# -----------------------------------------------------------------
