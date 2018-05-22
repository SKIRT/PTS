#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.services.s4g Contains the S4G class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles
from astropy.coordinates import Angle
from astroquery.vizier import Vizier

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import tables, introspection
from ...core.basics.map import Map
from ...modeling.preparation import unitconversion
from ...modeling.basics.models import SersicModel2D, ExponentialDiskModel2D
from ...core.basics.configurable import Configurable
from ..tools import catalogs
from ..basics.coordinate import SkyCoordinate
from ...modeling.basics.properties import GalaxyProperties
from ...core.tools import formatting as fmt
from ...core.units.parsing import parse_unit as u
from ...core.tools.utils import lazyproperty
from ...core.tools.angles import quadrant

# -----------------------------------------------------------------

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

def get_properties(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    # Run the S4G service for the specified galaxy
    s4g = S4G()
    s4g.config.galaxy_name = galaxy_name
    s4g.run()

    # Return the properties
    return s4g.properties

# -----------------------------------------------------------------

def get_components(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    # Run the S4G service for the specified galaxy
    s4g = S4G()
    s4g.config.galaxy_name = galaxy_name
    s4g.run()

    # Return the components
    return s4g.components

# -----------------------------------------------------------------

def get_table_lines():

    """
    This function ...
    :return:
    """

    lines = []
    with open(local_table_path, 'r') as s4g_table:
        for line in s4g_table:
            splitted = line.split()
            if len(splitted) < 2: continue
            lines.append(line)
    return lines

# -----------------------------------------------------------------

def get_galaxy_names():

    """
    This function ...
    :return:
    """

    names = []
    for line in get_table_lines():
        splitted = line.split()
        name = splitted[1]
        names.append(name)
    return names

# -----------------------------------------------------------------

def has_galaxy(name):

    """
    This function ...
    :param name:
    :return:
    """

    # Get the NGC name of the galaxy
    ngc_name = catalogs.get_ngc_name(name)
    ngc_name_nospaces = ngc_name.replace(" ", "")

    # Check
    return ngc_name_nospaces in get_galaxy_names()

# -----------------------------------------------------------------

class S4G(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(S4G, self).__init__(*args, **kwargs)

        # The inclination
        self.inclination = None

        # The galaxy properties object
        self.properties = None

        # The dictionary of components
        self.components = dict()

        # The Vizier service
        self.vizier = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 2. Get galaxy properties
        self.get_properties()

        # 3. Get the components
        self.get_components()

        # 4. Show
        if self.config.show: self.show()

        # 5. Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(S4G, self).setup(**kwargs)

        # Get the inclination
        self.inclination = kwargs.pop("inclination", None)

        # Create the galaxy properties object
        self.properties = GalaxyProperties()
        self.properties.name = self.config.galaxy_name

        # Get the NGC name of the galaxy
        self.properties.ngc_name = self.ngc_name

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_name(self):

        """
        This function ...
        :return:
        """

        return catalogs.get_ngc_name(self.config.galaxy_name)

    # -----------------------------------------------------------------

    @property
    def ngc_name_nospaces(self):

        """
        This function ...
        :return:
        """

        return self.ngc_name.replace(" ", "")

    # -----------------------------------------------------------------

    def get_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy properties ...")

        # S4G
        self.get_s4g_properties()

        # To first order, cos i = b/a where a & b
        # are the observed major and minor
        # axes (assume disk is intrinsically
        # circular)

        # But this implies galaxies have zero
        # thickness at i = 90°! Better to assume
        # that spirals are oblate ellipsoids with
        # intrinsic axis ratios a:a:c. If q = c/a,
        # then after a bit a simple geometry

        # cos^2(i) = [ (b/a)^2 - q^2 ] / (1-q^2)


        # The best fit to observed spirals yields
        # q = 0.13 for Sc galaxies

        # Ellipticity (1-b/a)

        # b_to_a = 1. - self.properties.ellipticity
        #
        # # Calculate the inclination
        # inclination_deg = 90. - math.degrees(math.acos(b_to_a))
        # inclination = Angle(inclination_deg, "deg")
        #
        # # Check the inclination
        # if self.inclination is not None:
        #     difference = abs(self.inclination - inclination)
        #     rel_difference = difference / self.inclination
        #     if rel_difference > 0.1:
        #         log.warning("The inclination angle calculated based on the decomposition differs by " + str(rel_difference*100) + "% from the specified inclination")
        #
        # # Set the incliantion
        # self.properties.inclination = inclination

    # -----------------------------------------------------------------

    def get_s4g_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the S4G catalog ...")

        # The Vizier querying object, specifying the necessary columns for this class
        self.vizier = Vizier(columns=['Name', 'RAJ2000', 'DEJ2000', 'Dmean', 'e_Dmean', 'amaj', 'ell', '[3.6]', '[4.5]', 'e_[3.6]', 'e_[4.5]', 'PA'])
        self.vizier.ROW_LIMIT = -1

        # Get parameters from S4G catalog
        result = self.vizier.query_object(self.config.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        # Galaxy name for S4G catalog
        self.properties.name = table["Name"][0]
        # Galaxy center from decomposition (?)
        ra_center = table["RAJ2000"][0]
        dec_center = table["DEJ2000"][0]
        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame='fk5')
        self.properties.center = center

        # Center position
        #self.properties.center = SkyCoordinate(ra=self.info["RA"][0], dec=self.info["DEC"][0], unit="deg") # center position from DustPedia
        
        # Distance
        self.properties.distance = table["Dmean"][0] * u("Mpc")
        self.properties.distance_error = table["e_Dmean"][0] * u("Mpc")

        # Major axis, ellipticity, position angle
        self.properties.major_arcsec = table["amaj"][0] * u("arcsec")
        self.properties.major = (self.properties.distance * self.properties.major_arcsec).to("pc", equivalencies=dimensionless_angles())

        # Ellipticity
        self.properties.ellipticity = table["ell"][0]
        self.properties.position_angle = Angle(table["PA"][0] + 90.0, u("deg"))

        # Magnitudes
        asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        #self.properties.i1_mag = asymptotic_ab_magnitude_i1
        #self.properties.i1_mag_error = asymptotic_ab_magnitude_i1_error
        #self.properties.i2_mag = asymptotic_ab_magnitude_i2
        #self.properties.i2_mag_error = asymptotic_ab_magnitude_i2_error

        #self.properties.i1_fluxdensity = unitconversion.ab_to_jansky(self.properties.i1_mag) * Unit("Jy")
        #i1_fluxdensity_lower = unitconversion.ab_to_jansky(
        #    self.properties.i1_mag + self.properties.i1_mag_error) * Unit("Jy")
        #i1_fluxdensity_upper = unitconversion.ab_to_jansky(
        #    self.properties.i1_mag - self.properties.i1_mag_error) * Unit("Jy")
        #i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=self.properties.i1_fluxdensity)
        #self.properties.i1_error = i1_error.average

        #self.properties.i2_fluxdensity = unitconversion.ab_to_jansky(self.properties.i2_mag) * Unit("Jy")
        #i2_fluxdensity_lower = unitconversion.ab_to_jansky(
        #    self.properties.i2_mag + self.properties.i2_mag_error) * Unit("Jy")
        #i2_fluxdensity_upper = unitconversion.ab_to_jansky(
        #    self.properties.i2_mag - self.properties.i2_mag_error) * Unit("Jy")
        #i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=self.properties.i2_fluxdensity)
        #self.properties.i2_error = i2_error.average

        # Other ...
        # absolute_magnitude_i1 = table["M3.6"][0]
        # absolute_magnitude_i2 = table["M4.5"][0]
        # stellar_mass = 10.0**table["logM_"][0] * u.Unit("Msun")

    # -----------------------------------------------------------------

    def get_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the components ...")

        # self.get_p4()  # currently (writing on the 31th of march 2016) there is a problem with the effective radius values
        # (at least for M81) on Vizier as well as in the PDF version of table 7 (S4G models homepage).

        # Parse the S4G table 8 to get the decomposition parameters
        self.get_parameters_from_table()

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
        result = self.vizier.query_object(self.config.galaxy_name, catalog=["J/ApJS/219/4/galaxies"])
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

        indices = tables.find_indices(table, self.ngc_name_nospaces, "Name")

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
                    value = unitconversion.ab_to_jansky(value) * u("Jy")
                elif parameter == "mu0": value = value * u("mag/arcsec2")
                elif parameter == "hr":
                    value = value * u("arcsec")
                    value = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())
                elif parameter == "hz":
                    value = value * u("arcsec")
                    value = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())

                component_parameters[parameter] = value

            if functionname == "sersic":

                re = table["Re"][index] * u("arcsec")
                component_parameters["Re"] = (self.parameters.distance * re).to("pc", equivalencies=dimensionless_angles())
                component_parameters["n"] = table["n"][index]

            elif functionname == "ferrer2":

                rbar = table["Rbar"][index] * u("arcsec")
                component_parameters["Rbar"] = (self.parameters.distance * rbar).to("pc", equivalencies=dimensionless_angles())

            if interpretation == "B": # bulge

                self.parameters.bulge = component_parameters

            elif interpretation == "D": # disk

                self.parameters.disk = component_parameters

            else: raise RuntimeError("Unrecognized component: " + interpretation)

        # Determine the full path to the parameters file
        #path = fs.join(self.components_path, "parameters.dat")

        # Write the parameters to the specified location
        #write_parameters(self.parameters, path)

    # -----------------------------------------------------------------

    def get_parameters_from_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the structural galaxy parameters from the S4G catalog ...")

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

        edgedisk_1_index = 28

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

        # EDGEDISK: phys frel mu0    rs    hs      pa

        with open(local_table_path, 'r') as s4g_table:

            for line in s4g_table:

                splitted = line.split()

                if len(splitted) < 2: continue

                name = splitted[1]

                #print(list(name))

                # Only look at the line corresponding to the galaxy
                if name != self.ngc_name_nospaces: continue

                #if self.ngc_name_nospaces not in name: continue

                #self.parameters.model_type = splitted[2]
                #self.parameters.number_of_components = splitted[3]
                #self.parameters.quality = splitted[4]

                ## BULGE

                #self.parameters.bulge.interpretation = splitted[sersic_1_index].split("|")[1]
                bulge_f = float(splitted[sersic_1_index + 1])
                mag = float(splitted[sersic_1_index + 2])
                bulge_fluxdensity = unitconversion.ab_to_jansky(mag) * u("Jy")

                # Effective radius in pc
                re_arcsec = float(splitted[sersic_1_index + 3]) * u("arcsec")
                bulge_re = (self.properties.distance * re_arcsec).to("pc", equivalencies=dimensionless_angles())

                bulge_q = float(splitted[sersic_1_index + 4])
                bulge_pa = Angle(float(splitted[sersic_1_index + 5]) - 90., "deg")
                bulge_n = float(splitted[sersic_1_index + 6])

                # Create the bulge component
                bulge = SersicModel2D(rel_contribution=bulge_f, fluxdensity=bulge_fluxdensity, effective_radius=bulge_re,
                                      axial_ratio=bulge_q, position_angle=bulge_pa, index=bulge_n)

                # Add the bulge to the components dictionary
                self.components["bulge"] = bulge

                ## DISK

                if splitted[disk_1_index + 1] != "-":

                    #self.parameters.disk.interpretation = splitted[disk_1_index].split("|")[1]
                    disk_f = float(splitted[disk_1_index + 1])
                    mag = float(splitted[disk_1_index + 2])
                    disk_fluxdensity = unitconversion.ab_to_jansky(mag) * u("Jy")

                    # Scale length in pc
                    hr_arcsec = float(splitted[disk_1_index + 3]) * u("arcsec")
                    disk_hr = (self.properties.distance * hr_arcsec).to("pc", equivalencies=dimensionless_angles())

                    disk_q = float(splitted[disk_1_index + 4]) # axial ratio
                    disk_pa = Angle(float(splitted[disk_1_index + 5]) - 90., "deg")
                    disk_mu0 = float(splitted[disk_1_index + 6]) * u("mag/arcsec2")

                    # Create the disk component
                    disk = ExponentialDiskModel2D(rel_contribution=disk_f, fluxdensity=disk_fluxdensity, scalelength=disk_hr,
                                                  axial_ratio=disk_q, position_angle=disk_pa, mu0=disk_mu0)

                else:

                    # phys frel mu0    rs    hs      pa

                    disk_f = float(splitted[edgedisk_1_index + 1])
                    mag = None
                    fluxdensity = None

                    # frel mu0    rs    hs      pa

                    #print(disk_f)

                    disk_mu0 = float(splitted[edgedisk_1_index + 2]) * u("mag/arcsec2")

                    #  rs    hs      pa

                    #  rs     =  radial scale lenght (arcsec)
                    #  hs     =  vertical scale height (arcsec)

                    rs = float(splitted[edgedisk_1_index + 3])
                    hs = float(splitted[edgedisk_1_index + 4])

                    #print(rs)
                    #print(hs)

                    disk_pa = Angle(float(splitted[edgedisk_1_index + 5]) - 90., "deg")

                    #print(disk_pa)

                    disk_q = hs / rs

                    # Create the disk component
                    disk = ExponentialDiskModel2D(rel_contribution=disk_f, scalelength=rs, axial_ratio=disk_q, position_angle=disk_pa, mu0=disk_mu0)

                # Add the disk to the components dictionary
                self.components["disk"] = disk

    # -----------------------------------------------------------------

    @property
    def disk_pa(self):

        """
        This function ...
        :return:
        """

        return self.components["disk"].position_angle

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        print(fmt.green + "Galaxy properties" + fmt.reset + ":")
        print("")
        print(self.properties)
        print("")

        for name in self.components:

            print(fmt.green + name + fmt.reset + ":")
            print("")
            print(self.components[name])
            print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the properties
        self.write_properties()

        # Write the components
        self.write_components()

    # -----------------------------------------------------------------

    def write_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the properties ...")

        # Determine the path and save
        path = self.output_path_file("properties.dat")
        self.properties.saveto(path)

    # -----------------------------------------------------------------

    def write_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the components ...")

        # Loop over the components
        for name in self.components:

            # Determine the path
            path = self.output_path_file(name + ".mod")

            # Save the model
            self.components[name].saveto(path)

# -----------------------------------------------------------------
