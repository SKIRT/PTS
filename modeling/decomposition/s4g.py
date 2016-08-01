#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.s4g Contains the S4GDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .component import DecompositionComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import tables, introspection
from astroquery.vizier import Vizier
from ...core.basics.map import Map
from ..preparation import unitconversion
from ..basics.models import SersicModel2D, ExponentialDiskModel2D

# -----------------------------------------------------------------

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = fs.join(introspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

class S4GDecomposer(DecompositionComponent):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(S4GDecomposer, self).__init__(config)

        # The Vizier querying object
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1

        # The path for the S4G decomposition models
        self.components_2d_s4g_path = None

        # The dictionary of components
        self.components = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

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

        #self.get_p4()  # currently (writing on the 31th of march 2016) there is a problem with the effective radius values

        # (at least for M81) on Vizier as well as in the PDF version of table 7 (S4G models homepage).

        # Parse the S4G table 8 to get the decomposition parameters
        self.get_parameters_from_table()

        # self.parameters.bulge.n = 2.62
        # self.parameters.bulge.PA = Angle(-31.9 + 90., "deg")
        # self.parameters.bulge.q = 0.71
        # value = 46.2 * Unit("arcsec")
        # self.parameters.bulge.Re = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())

        # value = 155.4 * Unit("arcsec")
        # self.parameters.disk.hr = (self.parameters.distance * value).to("pc", equivalencies=dimensionless_angles())
        # self.parameters.disk.q = 0.52
        # self.parameters.disk.PA = Angle(-28.3 + 90., "deg")

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(S4GDecomposer, self).setup()

        # Set the path
        self.components_2d_s4g_path = fs.create_directory_in(self.components_2d_path, "S4G")

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

                #print(list(name))

                # Only look at the line corresponding to the galaxy
                if name != self.ngc_id_nospaces: continue

                #if self.ngc_id_nospaces not in name: continue



                #self.parameters.model_type = splitted[2]
                #self.parameters.number_of_components = splitted[3]
                #self.parameters.quality = splitted[4]

                ## BULGE

                #self.parameters.bulge.interpretation = splitted[sersic_1_index].split("|")[1]
                bulge_f = float(splitted[sersic_1_index + 1])
                mag = float(splitted[sersic_1_index + 2])
                bulge_fluxdensity = unitconversion.ab_to_jansky(mag) * Unit("Jy")

                # Effective radius in pc
                re_arcsec = float(splitted[sersic_1_index + 3]) * Unit("arcsec")
                bulge_re = (self.galaxy_properties.distance * re_arcsec).to("pc", equivalencies=dimensionless_angles())

                bulge_q = float(splitted[sersic_1_index + 4])
                bulge_pa = Angle(float(splitted[sersic_1_index + 5]) - 90., "deg")
                bulge_n = float(splitted[sersic_1_index + 6])

                # Create the bulge component
                bulge = SersicModel2D(rel_contribution=bulge_f, fluxdensity=bulge_fluxdensity, effective_radius=bulge_re,
                                      axial_ratio=bulge_q, position_angle=bulge_pa, index=bulge_n)

                # Add the bulge to the components dictionary
                self.components["bulge"] = bulge

                ## DISK

                #self.parameters.disk.interpretation = splitted[disk_1_index].split("|")[1]
                disk_f = float(splitted[disk_1_index + 1])
                mag = float(splitted[disk_1_index + 2])
                disk_fluxdensity = unitconversion.ab_to_jansky(mag) * Unit("Jy")

                # Scale length in pc
                hr_arcsec = float(splitted[disk_1_index + 3]) * Unit("arcsec")
                disk_hr = (self.galaxy_properties.distance * hr_arcsec).to("pc", equivalencies=dimensionless_angles())

                disk_q = float(splitted[disk_1_index + 4]) # axial ratio
                disk_pa = Angle(float(splitted[disk_1_index + 5]) - 90., "deg")
                disk_mu0 = float(splitted[disk_1_index + 6]) * Unit("mag/arcsec2")

                # Create the disk component
                disk = ExponentialDiskModel2D(relative_contribution=disk_f, fluxdensity=disk_fluxdensity, scalelength=disk_hr,
                                              axial_ratio=disk_q, position_angle=disk_pa, mu0=disk_mu0)

                # Add the disk to the components dictionary
                self.components["disk"] = disk

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the components
        self.write_components()

    # -----------------------------------------------------------------

    def write_components(self):

        """
        This function ...
        :return:
        """

        # Loop over the components
        for name in self.components:

            # Determine the path
            path = fs.join(self.components_2d_s4g_path, name + ".mod")

            # Save the model
            self.components[name].save(path)

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
