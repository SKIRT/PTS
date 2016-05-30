#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.decomposition Contains the GalaxyDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astroquery.vizier import Vizier
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .component import DecompositionComponent
from ...core.basics.map import Map
from ...core.tools import inspection
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.skifile import SkiFile
from ..preparation import unitconversion
from ...core.basics.errorbar import ErrorBar
from ...core.simulation.arguments import SkirtArguments
from ...core.simulation.execute import SkirtExec
from ...magic.tools import catalogs
from ...magic.basics.vector import Extent, Position
from ...magic.basics.skygeometry import SkyEllipse, SkyCoordinate
from ...magic.basics.skyregion import SkyRegion
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.tools import tables
from ..basics.models import SersicModel, ExponentialDiskModel
from ..basics.instruments import SimpleInstrument
from ...magic.misc.kernels import AnianoKernels
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection

# -----------------------------------------------------------------

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = fs.join(inspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

# The path to the template ski files directory
template_path = fs.join(inspection.pts_dat_dir("modeling"), "ski")

# -----------------------------------------------------------------

class GalaxyDecomposer(DecompositionComponent):
    
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
        super(GalaxyDecomposer, self).__init__(config)

        # The NGC name of the galaxy
        self.ngc_id = None
        self.ngc_id_nospaces = None

        # The path to the disk and bulge directories
        self.bulge_directory = None
        self.disk_directory = None
        self.model_directory = None

        # The Vizier querying object
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1

        # The bulge and disk parameters
        self.parameters = Map()

        # The SKIRT execution context
        self.skirt = SkirtExec()

        # The bulge and disk model
        self.bulge = None
        self.disk = None

        # The bulge and disk image
        self.bulge_image = None
        self.disk_image = None
        self.model_image = None

        # The projection systems
        self.projections = dict()

        # The instruments
        self.instruments = dict()

        # The reference coordinate system
        self.reference_wcs = None

        # The PSF (of the reference image) for convolution with the simulated images
        self.psf = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new GalaxyDecomposer instance
        decomposer = cls(arguments.config)

        # Set the input and output path
        decomposer.config.path = arguments.path

        # Return the new instance
        return decomposer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Get the parameters
        self.get_parameters()

        # 3. Create the models
        self.create_models()

        # 4. Create the projection systems
        self.create_projections()

        # 4. Create the instruments
        self.create_instruments()

        # 2. Simulate the bulge and disk images
        self.simulate()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyDecomposer, self).setup()

        # TEMP: provide a cfg file for this class
        self.config.bulge_packages = 1e7
        self.config.disk_packages = 1e8

        # Get the NGC name of the galaxy
        self.ngc_id = catalogs.get_ngc_name(self.galaxy_name)
        self.ngc_id_nospaces = self.ngc_id.replace(" ", "")

        # Determine the path to the bulge and disk directories
        self.bulge_directory = self.full_output_path("bulge")
        self.disk_directory = self.full_output_path("disk")
        self.model_directory = self.full_output_path("model")

        # Create the bulge and disk directories
        fs.create_directories([self.bulge_directory, self.disk_directory, self.model_directory])

        # Get the coordinate system describing the pixel grid of the prepared images
        reference_path = fs.join(self.prep_path, self.reference_image, "result.fits")
        self.reference_wcs = CoordinateSystem.from_file(reference_path)

        # Load the PSF
        aniano = AnianoKernels()
        self.psf = aniano.get_psf("PACS_160")

    # -----------------------------------------------------------------

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the structural galaxy parameters from the S4G catalog ...")


        # Tracing spiral density waves in M81: (Kendall et. al, 2008)
        # - Sersic bulge:
        #    n = 2.62
        #    Re = 46.2 arcsec
        #    b/a = 0.71
        #    PA = −31.9°
        # - Exponential disk:
        #    Rs = 155.4 arcsec
        #    b/a = 0.52
        #    PA = −28.3°

        # Query the S4G catalog using Vizier for general parameters
        self.get_general_parameters()

        # Get the decomposition parameters
        #self.get_p4() currently (writing on 31 of march 2016) there is a problem with the effective radius values
        # (at least for M81) on Vizier as well as in the PDF version of table 7 (S4G models homepage).

        # Parse the S4G table 8 to get the decomposition parameters
        self.get_parameters_from_table()

        #self.parameters.bulge.n = 2.62
        #self.parameters.bulge.PA = Angle(-31.9 + 90., "deg")
        #self.parameters.bulge.q = 0.71
        #value = 46.2 * Unit("arcsec")
        #self.parameters.bulge.Re = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())

        #value = 155.4 * Unit("arcsec")
        #self.parameters.disk.hr = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())
        #self.parameters.disk.q = 0.52
        #self.parameters.disk.PA = Angle(-28.3 + 90., "deg")

    # -----------------------------------------------------------------

    def get_general_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the S4G catalog ...")

        # Get parameters from S4G catalog
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        # Galaxy name for S4G catalog
        self.parameters.galaxy_name = table["Name"][0]

        # Galaxy center from decomposition (?)
        ra_center = table["_RAJ2000"][0]
        dec_center = table["_DEJ2000"][0]
        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame='fk5')
        self.parameters.center = center

        # Distance
        self.parameters.distance = table["Dmean"][0] * Unit("Mpc")
        self.parameters.distance_error = table["e_Dmean"][0] * Unit("Mpc")

        # Major axis, ellipticity, position angle
        self.parameters.major_arcsec = table["amaj"][0] * Unit("arcsec")
        self.parameters.major = (self.parameters.distance * self.parameters.major_arcsec).to("pc", equivalencies=dimensionless_angles())

        # Ellipticity
        self.parameters.ellipticity = table["ell"][0]
        self.parameters.position_angle = Angle(table["PA"][0] + 90.0, Unit("deg"))

        # Magnitudes
        asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        self.parameters.i1_mag = asymptotic_ab_magnitude_i1
        self.parameters.i1_mag_error = asymptotic_ab_magnitude_i1_error
        self.parameters.i2_mag = asymptotic_ab_magnitude_i2
        self.parameters.i2_mag_error = asymptotic_ab_magnitude_i2_error

        self.parameters.i1_fluxdensity = unitconversion.ab_to_jansky(self.parameters.i1_mag) * Unit("Jy")
        i1_fluxdensity_lower = unitconversion.ab_to_jansky(self.parameters.i1_mag + self.parameters.i1_mag_error) * Unit("Jy")
        i1_fluxdensity_upper = unitconversion.ab_to_jansky(self.parameters.i1_mag - self.parameters.i1_mag_error) * Unit("Jy")
        i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=self.parameters.i1_fluxdensity)
        self.parameters.i1_error = i1_error.average

        self.parameters.i2_fluxdensity = unitconversion.ab_to_jansky(self.parameters.i2_mag) * Unit("Jy")
        i2_fluxdensity_lower = unitconversion.ab_to_jansky(self.parameters.i2_mag + self.parameters.i2_mag_error) * Unit("Jy")
        i2_fluxdensity_upper = unitconversion.ab_to_jansky(self.parameters.i2_mag - self.parameters.i2_mag_error) * Unit("Jy")
        i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=self.parameters.i2_fluxdensity)
        self.parameters.i2_error = i2_error.average

        # Other ...
        #absolute_magnitude_i1 = table["M3.6"][0]
        #absolute_magnitude_i2 = table["M4.5"][0]
        #stellar_mass = 10.0**table["logM_"][0] * u.Unit("Msun")

        # Inform the user
        log.info("Querying the catalog of radial profiles for 161 face-on spirals ...")

        # Radial profiles for 161 face-on spirals (Munoz-Mateos+, 2007)
        radial_profiles_result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/658/1006")

        distance = float(radial_profiles_result[0][0]["Dist"])
        inclination = Angle(float(radial_profiles_result[0][0]["i"]), "deg")

        # Set the inclination
        self.parameters.inclination = inclination

    # -----------------------------------------------------------------

    def get_p4(self):

        """
        This function ...
        :return:
        """

        #http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/ApJS/219/4

        # J/ApJS/219/4: S4G pipeline 4: multi-component decompositions (Salo+, 2015)

        #  - J/ApJS/219/4/galaxies: Parameters of the galaxies; center, outer orientation, sky background; and 1-component Sersic fits (tables 1 and 6) (2352 rows)
        #  - J/ApJS/219/4/table7: *Parameters of final multicomponent decompositions (Note) (4629 rows)

        # Inform the user
        log.info("Querying the S4G pipeline 4 catalog ...")

        # Get the "galaxies" table
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/ApJS/219/4/galaxies"])
        table = result[0]

        # PA: [0.2/180] Outer isophote position angle
        # e_PA: [0/63] Standard deviation in PA
        # Ell:  [0.008/1] Outer isophote ellipticity
        # e_Ell: [0/0.3] Standard deviation in Ell

        pa = Angle(table["PA"][0] - 90., "deg")
        pa_error = Angle(table["e_PA"][0], "deg")

        ellipticity = table["Ell"][0]
        ellipticity_error = table["e_Ell"][0]

        # Get the "table7" table
        result = self.vizier.get_catalogs("J/ApJS/219/4/table7")
        table = result[0]

        # Name: Galaxy name
        # Mod: Type of final decomposition model
        # Nc: [1/4] Number of components in the model (1-4)
        # Q: [3/5] Quality flag, 5=most reliable
        # C: Physical interpretation of the component
        # Fn: The GALFIT function used for the component (sersic, edgedisk, expdisk, ferrer2 or psf)
        # f1: [0.006/1] "sersic" fraction of the total model flux
        # mag1:  [7/19.4] "sersic" total 3.6um AB magnitude
        # q1: [0.1/1] "sersic" axis ratio
        # PA1: [0.08/180] "sersic" position angle [deg]
        # Re: [0.004/430] "sersic" effective radius (Re) [arcsec]
        # n: [0.01/20] "sersic" parameter n
        # f2: [0.02/1] "edgedisk" fraction of the total model flux
        # mu02: [11.8/24.6] "edgedisk" central surface face-on brightness (µ0) [mag/arcsec2]
        # PA2: [-90/90] "edgedisk" PA [deg]
        # hr2: [1/153] "edgedisk" exponential scale length (hr) [arcsec]
        # hz2: [0.003/39] "edgedisk" z-scale hz [arcsec]
        # f3: [0.02/1] "expdisk" fraction of the total model flux
        # mag3: [6.5/18.1] "expdisk" total 3.6um AB magnitude [mag]
        # q3: [0.1/1] "expdisk" axis ratio
        # PA3: [-90/90] "expdisk" position angle [deg]
        # hr3: [0.7/332] "expdisk" exponential scale length (hr) [arcsec]
        # mu03: [16.4/25.3] "expdisk" central surface face-on brightness (µ0) [mag/arcsec2]
        # f4: [0.003/0.6] "ferrer2" fraction of the total model flux
        # mu04: [16/24.8] "ferrer2" central surface sky brightness (µ0) [mag/arcsec2]
        # q4: [0.01/1] "ferrer2" axis ratio
        # PA4: [-90/90] "ferrer2" position angle [deg]
        # Rbar: [3.7/232.5] "ferrer2" outer truncation radius of the bar (Rbar) [arcsec]
        # f5: [0.001/0.4] "psf" fraction of the total model flux
        # mag5: [11.5/21.1] "psf" total 3.6um AB magnitude [mag]

        indices = tables.find_indices(table, self.ngc_id_nospaces, "Name")

        labels = {"sersic": 1, "edgedisk": 2, "expdisk": 3, "ferrer2": 4, "psf": 5}

        #units = {"f": None, "mag": "mag", "q": None, "PA": "deg", }

        # Loop over the indices
        for index in indices:

            model_type = table["Mod"][index]
            number_of_components = table["Nc"][index]
            quality = table["Q"][index]
            interpretation = table["C"][index]
            functionname = table["Fn"][index]

            component_parameters = Map()

            if self.parameters.model_type is not None: assert model_type == self.parameters.model_type
            if self.parameters.number_of_components is not None: assert number_of_components == self.parameters.number_of_components
            if self.parameters.quality is not None: assert quality == self.parameters.quality
            self.parameters.model_type = model_type
            self.parameters.number_of_components = number_of_components
            self.parameters.quality = quality

            for key in table.colnames:

                if not key.endswith(str(labels[functionname])): continue

                parameter = key[:-1]

                value = table[key][index]

                if parameter == "PA":
                    value = Angle(value + 90., "deg")
                    if quadrant(value) == 2: value = value - Angle(180., "deg")
                    elif quadrant(value) == 3: value = value + Angle(180., "deg")

                    if value.to("deg").value > 180.: value = value - Angle(360., "deg")
                    elif value.to("deg").value < -180.: value = value + Angle(360., "deg")
                elif parameter == "mag":
                    parameter = "fluxdensity"
                    value = unitconversion.ab_to_jansky(value) * Unit("Jy")
                elif parameter == "mu0": value = value * Unit("mag/arcsec2")
                elif parameter == "hr":
                    value = value * Unit("arcsec")
                    value = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())
                elif parameter == "hz":
                    value = value * Unit("arcsec")
                    value = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())

                component_parameters[parameter] = value

            if functionname == "sersic":

                re = table["Re"][index] * Unit("arcsec")
                component_parameters["Re"] = (self.parameters.distance * re).to("pc", equivalencies=dimensionless_angles())
                component_parameters["n"] = table["n"][index]

            elif functionname == "ferrer2":

                rbar = table["Rbar"][index] * Unit("arcsec")
                component_parameters["Rbar"] = (self.parameters.distance * rbar).to("pc", equivalencies=dimensionless_angles())

            if interpretation == "B": # bulge

                self.parameters.bulge = component_parameters

            elif interpretation == "D": # disk

                self.parameters.disk = component_parameters

            else: raise RuntimeError("Unrecognized component: " + interpretation)

         # Determine the full path to the parameters file
        path = self.full_output_path("parameters.dat")

        # Write the parameters to the specified location
        write_parameters(self.parameters, path)

    # -----------------------------------------------------------------

    def get_parameters_from_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Parsing S4G table 8 to get the decomposition parameters ...")

        #The table columns are:
        # (1) the running number (1-2352),
        # (2) the galaxy name,
        # (3) the type of final decomposition model,
        # (4) the number N of components in the final model , and
        # (5) the quality flag Q.

        # 6- 8 for unresolved central component ('psf'; phys, frel,  mag),
        # 9-15 for inner sersic-component ('sersic1'; phys,  frel,  mag,    re,    ar,      pa,    n),
        # 16-22 for inner disk-component ('expo1'; phys,  frel,  mag,    hr,    ar,      pa,   mu0),
        # 23-28 for inner ferrers-component ('ferrers1'; phys,  frel,  mu0,   rout,   ar,     pa),
        # 29-34 for inner edge-on disk component ('edgedisk1';  phys, frel, mu0,  rs,    hs,      pa ).
        # 35-41 for outer sersic-component ('sersic2'; phys,  frel,  mag,    re,    ar,      pa,    n),
        # 42-49 for outer disk-component ('expo2'; phys,  frel,  mag,    hr,    ar,      pa,   mu0),
        # 50-55 for outer ferrers-component ('ferrers2'; phys,  frel,  mu0,   rout,   ar,     pa),
        # 56-61 for outer edge-on disk component ('edgedisk2';  phys, frel, mu0,  rs,    hs,      pa ).

        sersic_1_index = 8
        disk_1_index = 15

        #For each function:

        #the first entry stands for the physical intepretation of the component:
        #'N' for a central source (or unresolved bulge), 'B' for a bulge (or elliptical), 'D' for a disk, 'BAR' for a bar, and 'Z' for an edge-on disk.

        # 'rel    =  the relative contribution of the component to the total model flux,
        #  mag    =  the component's total 3.6 micron AB magnitude,
        #  mu0    =  the central  surface brightness (mag/arcsec^2;  de-projected central surface brightness for expdisk and edgedisk, and
        #                                                          sky-plane central surface brightness for ferrer)
        #  ar     =  axial ratio
        #  pa     =  position angle (degrees ccw from North)
        #  n      =  sersic index
        #  hr     =  exponential scale lenght (arcsec)
        #  rs     =  radial scale lenght (arcsec)
        #  hs     =  vertical scale height (arcsec)
        #  rout   =  bar outer truncation radius (arcsec)

        # SERSIC: phys,  frel,  mag,    re,    ar,      pa,    n
        # DISK: phys,  frel,  mag,    hr,    ar,      pa,   mu0

        with open(local_table_path, 'r') as s4g_table:

            for line in s4g_table:

                splitted = line.split()

                if len(splitted) < 2: continue

                name = splitted[1]

                # Only look at the line corresponding to the galaxy
                if name != self.ngc_id_nospaces: continue

                self.parameters.model_type = splitted[2]
                self.parameters.number_of_components = splitted[3]
                self.parameters.quality = splitted[4]

                # BULGE
                self.parameters.bulge = Map()
                #self.parameters.bulge.interpretation = splitted[sersic_1_index].split("|")[1]
                self.parameters.bulge.f = float(splitted[sersic_1_index + 1])
                mag = float(splitted[sersic_1_index + 2])
                self.parameters.bulge.fluxdensity = unitconversion.ab_to_jansky(mag) * Unit("Jy")

                # Effective radius in pc
                re_arcsec = float(splitted[sersic_1_index + 3]) * Unit("arcsec")
                self.parameters.bulge.Re = (self.parameters.distance * re_arcsec).to("pc", equivalencies=dimensionless_angles())

                self.parameters.bulge.q = float(splitted[sersic_1_index + 4])
                self.parameters.bulge.PA = Angle(float(splitted[sersic_1_index + 5]) - 90., "deg")
                self.parameters.bulge.n = float(splitted[sersic_1_index + 6])

                # DISK
                self.parameters.disk = Map()
                #self.parameters.disk.interpretation = splitted[disk_1_index].split("|")[1]
                self.parameters.disk.f = float(splitted[disk_1_index + 1])
                mag = float(splitted[disk_1_index + 2])
                self.parameters.disk.fluxdensity = unitconversion.ab_to_jansky(mag) * Unit("Jy")

                # Scale length in pc
                hr_arcsec = float(splitted[disk_1_index + 3]) * Unit("arcsec")
                self.parameters.disk.hr = (self.parameters.distance * hr_arcsec).to("pc", equivalencies=dimensionless_angles())

                self.parameters.disk.q = float(splitted[disk_1_index + 4]) # axial ratio
                self.parameters.disk.PA = Angle(float(splitted[disk_1_index + 5]) - 90., "deg")
                self.parameters.disk.mu0 = float(splitted[disk_1_index + 6]) * Unit("mag/arcsec2")

    # -----------------------------------------------------------------

    def get_spiral_properties(self):

        """
        This function ...
        :return:
        """

        # http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/582/A86

        # J/A+A/582/A86: Catalogue of features in the S4G (Herrera-Endoqui+, 2015)

        # - J/A+A/582/A86/table2: Properties of bars, ring- and lens-structures in the S4G (2387 rows)
        # - J/A+A/582/A86/table3: Properties of spiral arms in the S4G (1854 rows)

        # Get table2
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/A+A/582/A86/table2"])
        table = result[0]

        # Name: Galaxy name
        # Class: Morphological classification
        # Type: Type of feature
        # sma: Semi-major axis [arcsec]
        # PA: Position angle [deg]
        # Ell: Ellipticity
        # smaEll: Semi-major axis from ellipticity [arcsec]
        # dsma: Deprojected semi-major axis [arcsec]
        # dPA: Deprojected position angle [deg]
        # dEll: Deprojected ellipticity
        # dsmaEll: Deprojexted semi-major axis from Ell [arcsec]
        # Qual: [1/3] Quality flag

        # Get table 3
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/A+A/582/A86/table3"])
        table = result[0]

        # Name: Galaxy name
        # Class: Morphological classification
        # Type: Type of arms
        # Segment: Segment
        # Pitchang: Pitch angle [deg]
        # ri: Inner radius [arcsec]
        # ro: Outer radius [arcsec]
        # Qual: [1/2] Quality flag

    # -----------------------------------------------------------------

    def create_models(self):

        """
        :return:
        """

        # Create the bulge model
        self.create_bulge_model()

        # Create the disk model
        self.create_disk_model()

    # -----------------------------------------------------------------

    def create_bulge_model(self):

        """
        :return:
        """

        # Create a Sersic model for the bulge
        self.bulge = SersicModel.from_galfit(self.parameters.bulge, self.parameters.inclination, self.parameters.disk.PA)

    # -----------------------------------------------------------------

    def create_disk_model(self):

        """
        :return:
        """

        # Create an exponential disk model for the disk
        self.disk = ExponentialDiskModel.from_galfit(self.parameters.disk, self.parameters.inclination, self.parameters.disk.PA)

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        # Create the 'earth' projection system
        azimuth = 0.0
        self.projections["earth"] = GalaxyProjection.from_wcs(self.reference_wcs, self.parameters.center, self.parameters.distance,
                                                              self.parameters.inclination, azimuth, self.parameters.disk.PA)

        # Create the face-on projection system
        self.projections["faceon"] = FaceOnProjection.from_wcs(self.reference_wcs, self.parameters.center, self.parameters.distance)

        # Create the edge-on projection system
        self.projections["edgeon"] = EdgeOnProjection.from_wcs(self.reference_wcs, self.parameters.center, self.parameters.distance)

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projection systems
        for name in self.projections:

            # Create the instrument from the projection system
            self.instruments[name] = SimpleInstrument.from_projection(self.projections[name])

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Simulate the stellar bulge without deprojection
        self.simulate_bulge_simple()

        # Simulate the stellar bulge
        self.simulate_bulge()

        # Simulate the stellar disk
        self.simulate_disk()

        # Simulate the bulge + disk model
        self.simulate_model()

    # -----------------------------------------------------------------

    def simulate_bulge_simple(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = fs.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.bulge_packages)

        # Change the ski file parameters
        # component_id, index, radius, y_flattening=1, z_flattening=1
        ski.set_stellar_component_sersic_geometry(0, self.parameters.bulge.n, self.parameters.bulge.Re, y_flattening=self.parameters.bulge.q)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Create the instrument
        distance = self.parameters.distance
        inclination = 0.0
        azimuth = Angle(90., "deg")
        #position_angle = self.parameters.bulge.PA + Angle(90., "deg") # + 90° because we can only do y_flattening and not x_flattening
        position_angle = self.parameters.bulge.PA
        pixels_x = self.reference_wcs.xsize
        pixels_y = self.reference_wcs.ysize
        pixel_center = self.parameters.center.to_pixel(self.reference_wcs)
        center = Position(0.5*pixels_x - pixel_center.x - 0.5, 0.5*pixels_y - pixel_center.y - 0.5)
        center_x = center.x * Unit("pix")
        center_y = center.y * Unit("pix")
        center_x = (center_x * self.reference_wcs.pixelscale.x.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * self.reference_wcs.pixelscale.y.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        field_x_angular = self.reference_wcs.pixelscale.x.to("deg/pix") * pixels_x * Unit("pix")
        field_y_angular = self.reference_wcs.pixelscale.y.to("deg/pix") * pixels_y * Unit("pix")
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())
        fake = SimpleInstrument(distance, inclination, azimuth, position_angle, field_x_physical, field_y_physical, pixels_x, pixels_y, center_x, center_y)

        # Add the instrument
        ski.add_instrument("earth", fake)

        # Create the directory to simulate the bulge
        simple_bulge_directory = self.full_output_path("bulge_simple")
        fs.create_directory(simple_bulge_directory)

        # Determine the path to the ski file
        ski_path = fs.join(simple_bulge_directory, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.join(simple_bulge_directory, "out")

        # Create the output directory
        fs.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        bulge_image_path = fs.join(out_path, "bulge_earth_total.fits")

        # Check if the output contains the "bulge_earth_total.fits" file
        if not fs.is_file(bulge_image_path):
            raise RuntimeError("Something went wrong with the simple bulge simulation: output FITS file missing")

    # -----------------------------------------------------------------

    def simulate_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = fs.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.bulge_packages)

        # Set the bulge geometry
        ski.set_stellar_component_geometry(0, self.bulge)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.bulge_directory, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.join(self.bulge_directory, "out")

        # Create the output directory
        fs.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        bulge_image_path = fs.join(out_path, "bulge_earth_total.fits")

        # Check if the output contains the "bulge_earth_total.fits" file
        if not fs.is_file(bulge_image_path):
            raise RuntimeError("Something went wrong with the bulge simulation: output FITS file missing")

        # Open the bulge image
        self.bulge_image = Frame.from_file(bulge_image_path)

        # Set the coordinate system of the bulge image
        self.bulge_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the bulge image to the bulge flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.parameters.bulge.fluxdensity
        self.bulge_image *= fluxdensity.to("Jy").value / np.sum(self.bulge_image)
        self.bulge_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the bulge image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.bulge_image = self.bulge_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def simulate_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the disk image ...")

        # Load the disk ski file template
        disk_template_path = fs.join(template_path, "disk.ski")
        ski = SkiFile(disk_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.disk_packages)

        # Change the ski file parameters
        ski.set_stellar_component_geometry(0, self.disk)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.disk_directory, "disk.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.join(self.disk_directory, "out")

        # Create the output directory
        fs.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the disk simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        disk_image_path = fs.join(out_path, "disk_earth_total.fits")

        # Check if the output contains the "disk_earth_total.fits" file
        if not fs.is_file(disk_image_path):
            raise RuntimeError("Something went wrong with the disk simulation: output FITS file missing")

        # Open the disk image
        self.disk_image = Frame.from_file(disk_image_path)

        # Set the coordinate system of the disk image
        self.disk_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the disk image to the disk flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.parameters.disk.fluxdensity
        self.disk_image *= fluxdensity.to("Jy").value / np.sum(self.disk_image)
        self.disk_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the disk image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.disk_image = self.disk_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def simulate_model(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge+disk model image ...")

        # Load the disk ski file template
        disk_template_path = fs.join(template_path, "model.ski")
        ski = SkiFile(disk_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.disk_packages)

        # Change the ski file parameters
        ski.set_stellar_component_geometry(0, self.disk)
        ski.set_stellar_component_geometry(1, self.bulge)

        # Set the luminosities of the two components
        ski.set_stellar_component_luminosities(0, [self.parameters.disk.f])
        ski.set_stellar_component_luminosities(1, [self.parameters.bulge.f])

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.model_directory, "model.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.join(self.model_directory, "out")

        # Create the output directory
        fs.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the disk+bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        model_image_path = fs.join(out_path, "model_earth_total.fits")

        # Check if the output contains the "model_earth_total.fits" file
        if not fs.is_file(model_image_path):
            raise RuntimeError("Something went wrong with the disk+bulge simulation: output FITS file missing")

        # Open the model image
        self.model_image = Frame.from_file(model_image_path)

        # Set the coordinate system of the model image
        self.model_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the model image to the bulge+disk flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.parameters.bulge.fluxdensity + self.parameters.disk.fluxdensity # sum of bulge and disk component flux density
        self.model_image *= fluxdensity.to("Jy").value / np.sum(self.model_image)
        self.model_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the model image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.model_image = self.model_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the final bulge and disk images
        self.write_images()

        # Write out the parameters in a data file
        self.write_parameters()

        # Write the projection systems
        self.write_projections()

        # Write out the disk ellipse
        self.write_disk_ellipse()

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the bulge image and save it
        final_bulge_path = fs.join(self.components_path, "bulge.fits")
        self.bulge_image.save(final_bulge_path)

        # Determine the path to the disk image and save it
        final_disk_path = fs.join(self.components_path, "disk.fits")
        self.disk_image.save(final_disk_path)

        # Determine the path to the model image and save it
        final_model_path = fs.join(self.components_path, "model.fits")
        self.model_image.save(final_model_path)

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing data file with parameters ...")

        # Determine the full path to the parameters file
        path = self.full_output_path("parameters.dat")

        # Write the parameters to the specified location
        write_parameters(self.parameters, path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Loop over the projection systems
        for name in self.projections:

            # Determine the path to the projection file
            path = self.full_output_path(name + ".proj")

            # Write the projection system
            self.projections[name].save(path)

    # -----------------------------------------------------------------

    def write_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing regions file with disk ellipse ...")

        minor = (1.0 - self.parameters.ellipticity) * self.parameters.major_arcsec

        # Ellipse radius
        radius = Extent(self.parameters.major_arcsec, minor)

        # Create sky ellipse
        sky_ellipse = SkyEllipse(self.parameters.center, radius, self.parameters.position_angle)

        # Create region
        region = SkyRegion()
        region.append(sky_ellipse)
        region_path = self.full_output_path("disk.reg")
        region.save(region_path)

# -----------------------------------------------------------------

def write_parameters(parameters, path):

    """
    This function ...
    :param parameters:
    :param path:
    :return:
    """

    # Create the parameters file
    with open(path, 'w') as parameter_file:

        # Add general info
        print("Name:", parameters.galaxy_name, file=parameter_file)
        print("Center RA:", str(parameters.center.ra.to("deg").value) + " deg", file=parameter_file)
        print("Center DEC:", str(parameters.center.dec.to("deg").value) + " deg", file=parameter_file)
        print("Major axis length:", str(parameters.major), file=parameter_file)
        print("Ellipticity:", parameters.ellipticity, file=parameter_file)
        print("Position angle:", str(parameters.position_angle.to("deg").value) + " deg", file=parameter_file)
        print("Distance:", str(parameters.distance), file=parameter_file)
        print("Distance error:", str(parameters.distance_error), file=parameter_file)
        print("Inclination:", str(parameters.inclination.to("deg").value) + " deg", file=parameter_file)
        print("IRAC 3.6um flux density:", parameters.i1_fluxdensity, file=parameter_file)
        print("IRAC 3.6um flux density error:", parameters.i1_error, file=parameter_file)
        print("IRAC 4.5um flux density:", parameters.i2_fluxdensity, file=parameter_file)
        print("IRAC 4.5um flux density error:", parameters.i2_error, file=parameter_file)

        #print("Model type:", parameters.model_type, file=parameter_file)
        #print("Number of components:", parameters.number_of_components, file=parameter_file)
        #print("Quality:", parameters.quality, file=parameter_file)

        # Add components parameters
        for component in ["bulge", "disk"]:

            #print(component.title() + ": Interpretation:", parameters[component].interpretation, file=parameter_file)
            print(component.title() + ": Relative contribution:", parameters[component].f, file=parameter_file)
            print(component.title() + ": IRAC 3.6um flux density:", parameters[component].fluxdensity, file=parameter_file)
            print(component.title() + ": Axial ratio:", parameters[component].q, file=parameter_file)
            print(component.title() + ": Position angle:", str(parameters[component].PA.to("deg").value) + " deg", file=parameter_file) # (degrees ccw from North)

            if component == "bulge":

                print(component.title() + ": Effective radius:", str(parameters[component].Re), file=parameter_file)
                print(component.title() + ": Sersic index:", parameters[component].n, file=parameter_file)

            elif component == "disk":

                print(component.title() + ": Central surface brightness:", parameters[component].mu0, file=parameter_file) # (mag/arcsec2)
                print(component.title() + ": Exponential scale length:", str(parameters[component].hr), file=parameter_file)

# -----------------------------------------------------------------

def load_parameters(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Create parameters structure
    parameters = Map()
    parameters.bulge = Map()
    parameters.disk = Map()

    ra = None
    dec = None

    # Read the parameter file
    with open(path, 'r') as parameter_file:

        # Loop over all lines in the file
        for line in parameter_file:

            splitted = line.split(": ")

            # Bulge parameters
            if splitted[0] == "Bulge":

                if splitted[1] == "Relative contribution": parameters.bulge.f = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.bulge.fluxdensity = get_quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.bulge.q = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.bulge.PA = get_angle(splitted[2])
                elif splitted[1] == "Effective radius": parameters.bulge.Re = get_quantity(splitted[2])
                elif splitted[1] == "Sersic index": parameters.bulge.n = float(splitted[2])

            # Disk parameters
            elif splitted[0] == "Disk":

                if splitted[1] == "Relative contribution": parameters.disk.f = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.disk.fluxdensity = get_quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.disk.q = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.disk.PA = get_angle(splitted[2])
                elif splitted[1] == "Central surface brightness": parameters.disk.mu0 = get_quantity(splitted[2])
                elif splitted[1] == "Exponential scale length": parameters.disk.hr = get_quantity(splitted[2])

            # Other parameters
            elif len(splitted) == 2:

                if splitted[0] == "Name": parameters.galaxy_name = splitted[1]
                elif splitted[0] == "Center RA": ra = get_quantity(splitted[1])
                elif splitted[0] == "Center DEC": dec = get_quantity(splitted[1])
                elif splitted[0] == "Major axis length": parameters.major = get_quantity(splitted[1])
                elif splitted[0] == "Ellipticity": parameters.ellipticity = float(splitted[1])
                elif splitted[0] == "Position angle": parameters.position_angle = get_angle(splitted[1])
                elif splitted[0] == "Distance": parameters.distance = get_quantity(splitted[1])
                elif splitted[0] == "Distance error": parameters.distance_error = get_quantity(splitted[1])
                elif splitted[0] == "Inclination": parameters.inclination = get_angle(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density": parameters.i1_fluxdensity = get_quantity(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density error": parameters.i1_error = get_quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density": parameters.i2_fluxdensity = get_quantity(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density error": parameters.i2_error = get_quantity(splitted[1])

    # Add the center coordinate
    parameters.center = SkyCoordinate(ra=ra, dec=dec)

    # Return the parameters
    return parameters

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object and return it
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def get_angle(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create an Angle object and return it
    if unit is not None: value = Angle(value, unit)
    return value

# -----------------------------------------------------------------

def quadrant(angle):

    """
    This function ...
    :param angle:
    :return:
    """

    if -180 <= angle.to("deg").value < -90.: return 3
    elif -90. <= angle.to("deg").value < 0.0: return 4
    elif 0.0 <= angle.to("deg").value < 90.: return 1
    elif 90. <= angle.to("deg").value < 180.: return 2
    elif 180. <= angle.to("deg").value < 270.: return 3
    elif 270. <= angle.to("deg").value <= 360.: return 4
    else: raise ValueError("Failed to determine quadrant for " + str(angle))

# -----------------------------------------------------------------
