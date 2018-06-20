#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.simulation.simulations Contains the ComponentSimulations class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty

# Import the relevant PTS classes and modules
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.data.attenuation import AttenuationCurve
from ...magic.tools import extinction
from ...magic.core.frame import Frame
from ...magic.core.datacube import DataCube
from ...core.basics.log import log
from .simulation import earth_name, faceon_name, edgeon_name
from .simulation import total_contribution, scattered_contribution, direct_contribution, transparent_contribution
from .simulation import dust_contribution, dust_direct_contribution, dust_scattered_contribution
from ...core.units.parsing import parse_unit as u
from ...magic.core.list import uniformize

# -----------------------------------------------------------------

stellar_dust_sed_split_wavelength = 5. * u("micron")

# -----------------------------------------------------------------

class ComponentSimulations(object):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, name, observed, distance=None, earth_wcs=None):

        """
        The constructor ...
        :param name:
        :param observed:
        :param distance:
        :param earth_wcs:
        """

        # Set the name
        self.name = name

        # Set the observed simulation
        self.observed = observed

        # Set the distance
        if distance is not None: self.distance = distance

        # Set the coordinate system
        self.earth_wcs = earth_wcs

    # -----------------------------------------------------------------

    @property
    def has_earth_wcs(self):

        """
        This function ...
        :return:
        """

        return self.earth_wcs is not None

    # -----------------------------------------------------------------

    @property
    def has_observed(self):

        """
        This function ...
        :return:
        """

        return self.observed is not None

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.observed.distance

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        if not self.has_observed: log.warning("The simulation of the observed component does not exist: cannot set the distance")
        else: self.observed.distance = value

    # -----------------------------------------------------------------

    @property
    def observed_output_path(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed: return None
        else: return self.observed.output_path

    # -----------------------------------------------------------------

    @property
    def observed_output(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed: return None
        else: return self.observed.output

    # -----------------------------------------------------------------

    @property
    def observed_data(self):

        """
        This function ...
        :return:
        """

        if not self.has_observed: return None
        else: return self.observed.data

    # -----------------------------------------------------------------

    @property
    def has_observed_output(self):

        """
        This function ...
        :return:
        """

        return self.has_observed and self.observed.has_output

    # -----------------------------------------------------------------

    @property
    def has_observed_data(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_output and self.observed_data.has_any

    # -----------------------------------------------------------------

    @property
    def has_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_observed and self.has_observed_output and self.observed_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed and self.has_observed_output and self.observed_data.has_images

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_sed_faceon(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_sed_edgeon(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_cube_faceon(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def has_faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube_faceon

    # -----------------------------------------------------------------

    @abstractproperty
    def has_intrinsic_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def has_edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.wavelength_grid()

    # -----------------------------------------------------------------
    # SEDs
    # -----------------------------------------------------------------

    @property
    def observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_data.sed_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.observed_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_observed_sed_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.observed_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_sed(self):

        """
        Thisn function ...
        :return:
        """

        return self.has_edgeon_observed_sed_orientation

    # -----------------------------------------------------------------

    @property
    def faceon_observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.sed_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_sed(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.seds[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.seds[faceon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_sed_contributions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and len(self.observed_data.seds[edgeon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_transparent_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_sed_contributions

    # -----------------------------------------------------------------

    @property
    def observed_sed_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_sed_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[faceon_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_sed_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[edgeon_name][dust_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed_dust - self.observed_sed_dust_scattered

    # -----------------------------------------------------------------

    @property
    def observed_sed_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_sed_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.seds[earth_name][transparent_contribution]

    # -----------------------------------------------------------------
    # CUBES
    # -----------------------------------------------------------------

    @property
    def observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.observed_data.image_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_earth_observed_cube_orientation(self):

        """
        Thisfunction ...
        :return:
        """

        return earth_name in self.observed_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_earth_observed_cube(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_earth_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.observed_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.observed_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_edgeon_observed_cube_orientation

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.image_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_earth_observed_cube and len(self.observed_data.images[earth_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_faceon_observed_cube and len(self.observed_data.images[faceon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_other_observed_cube_contributions_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_edgeon_observed_cube and len(self.observed_data.images[edgeon_name]) > 1

    # -----------------------------------------------------------------

    @property
    def has_transparent_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_earth_observed_cube and self.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_transparent_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_faceon_observed_cube and self.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_transparent_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_data and self.has_edgeon_observed_cube and self.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    @property
    def observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][direct_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][direct_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_direct_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_direct.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_direct_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.faceon_observed_cube_direct.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_direct_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_observed_cube_direct.integrate()

    # -----------------------------------------------------------------

    @property
    def observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][scattered_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_scattered_frame(self):

        """
        Thins function ...
        :return:
        """

        return self.observed_cube_scattered.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_scattered_frame_faceon(self):

        """
        This function ...
        :return:
        """

        return self.faceon_observed_cube_scattered.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_scattered_frame_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_observed_cube_scattered.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        # kwargs["no_fwhm"] = "return"
        #         kwargs["no_pixelscale"] = "shape"
        #         kwargs["distance"] = self.galaxy_distance

        #intrinsic = self.intrinsic_stellar_cube
        #observed = self.observed_stellar_cube
        #return self.intrinsic_stellar_cube - self.observed_stellar_cube

        #from ...magic.core.list import convolve_rebin_and_convert
        #intrinsic, observed = convolve_rebin_and_convert(self.intrinsic_stellar_cube, self.observed_stellar_cube, no_fwhm="return", no_pixelscale="shape", distance=self.distance)
        intrinsic, observed = uniformize(self.intrinsic_stellar_cube, self.observed_stellar_cube, distance=self.distance)

        # Return
        return intrinsic - observed

    # -----------------------------------------------------------------

    @property
    def has_observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_cube and self.has_observed_stellar_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube_faceon - self.faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_cube_faceon and self.has_faceon_observed_stellar_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube_edgeon - self.edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_cube_absorbed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_stellar_cube_edgeon and self.has_edgeon_observed_stellar_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_cube_attenuated(self):

        """
        This function ...
        :return:
        """

        # attenuated = transparent - direct stellar
        return self.intrinsic_stellar_cube - self.observed_cube_direct

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_cube_attenuated(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube_faceon - self.faceon_observed_cube_direct

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_cube_attenuated(self):

        """
        Thsi function ...
        :return:
        """

        return self.intrinsic_stellar_cube_edgeon - self.edgeon_observed_cube_direct

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][dust_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_dust(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][dust_contribution]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube_dust - self.observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.faceon_observed_cube_dust - self.faceon_observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_cube_dust_direct(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_observed_cube_dust - self.edgeon_observed_cube_dust_scattered

    # -----------------------------------------------------------------

    @property
    def observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_dust_scattered(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][dust_scattered_contribution]

    # -----------------------------------------------------------------

    @property
    def observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[earth_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[faceon_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_observed_cube_transparent(self):

        """
        This function ...
        :return:
        """

        return self.observed_data.images[edgeon_name][transparent_contribution]

    # -----------------------------------------------------------------

    @property
    def has_cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_cell_properties

    # -----------------------------------------------------------------

    @property
    def cell_volumes(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_volumes

    # -----------------------------------------------------------------

    @property
    def cell_dust_densities(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_dust_densities

    # -----------------------------------------------------------------

    @property
    def cell_mass_fractions(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_mass_fractions

    # -----------------------------------------------------------------

    @property
    def cell_optical_depths(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_optical_depths

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_x_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_y_coordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_z_coordinates

    # -----------------------------------------------------------------

    @property
    def has_grid_files(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_grid_files

    # -----------------------------------------------------------------

    @property
    def grid_xy_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xy_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_yz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_yz_filepath

    # -----------------------------------------------------------------

    @property
    def grid_xyz_filepath(self):

        """
        This function ...
        :return:
        """

        return self.observed.grid_xyz_filepath

    # -----------------------------------------------------------------

    @property
    def has_cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.observed.has_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def cell_stellar_density(self):

        """
        This function ...
        :return:
        """

        return self.observed.cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def has_full_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_sed_contributions

    # -----------------------------------------------------------------

    @property
    def has_full_sed_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_sed_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_full_sed_edgeon(self):

        """
        This function ..
        :return:
        """

        return self.has_other_observed_sed_contributions_edgeon

    # -----------------------------------------------------------------

    @property
    def has_full_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_cube_contributions

    # -----------------------------------------------------------------

    @property
    def has_full_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_cube_contributions_faceon

    # -----------------------------------------------------------------

    @property
    def has_full_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_other_observed_cube_contributions_edgeon

    # -----------------------------------------------------------------

    def get_stellar_part(self, sed_or_cube, fit_stellar=False, fit_dust=False, full=False):

        """
        This function ...
        :param sed_or_cube:
        :param fit_stellar: improve the accuracy by fitting a smooth function to the stellar part of the SED (or cube pixels)
        :param fit_dust: improve the accuracy by fitting a smooth function (MBB) to the dust part of the SED (or cube pixels), and subtracting this from the total
        :param full: keep the full wavelength grid
        :return:
        """

        # Fit stuff
        if fit_stellar or fit_dust: raise NotImplementedError("Not yet implemented")

        # Don't fit, just split at a certain wavelength
        else:

            if full: return sed_or_cube.flattened_above(stellar_dust_sed_split_wavelength)
            else: return sed_or_cube.splice(max_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    def get_dust_part(self, sed_or_cube, fit_stellar=False, fit_dust=False, full=False):

        """
        This function ...
        :param sed_or_cube:
        :param fit_stellar: improve the accuracy by fitting a smooth function to the stellar part of the SED (or cube pixels), and subtracting this from the total
        :param fit_dust: improve the accuracy by fitting a smooth function (MBB) to the dust part of the SED (or cube pixels)
        :param full: keep the full wavelength grid
        :return:
        """

        # Fit stuff?
        if fit_stellar or fit_dust: raise NotImplementedError("Not yet implemented")

        # Don't fit, just split at a certain wavelength
        else:

            if full: return sed_or_cube.flattened_below(stellar_dust_sed_split_wavelength)
            else: return sed_or_cube.splice(min_wavelength=stellar_dust_sed_split_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        # Has full instrument
        if self.has_full_sed and self.has_observed_diffuse_dust_sed: return self.observed_sed - self.observed_diffuse_dust_sed

        elif self.has_full_sed and self.has_observed_dust_sed: return self.observed_sed - self.observed_dust_sed + self.intrinsic_dust_sed

        # No full instrument data: get the stellar part of the total SED
        else: return self.get_stellar_part(self.observed_sed, full=True) + self.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed and self.has_observed_diffuse_dust_sed) or (self.has_full_sed and self.has_observed_dust_sed and self.has_intrinsic_dust_sed) or self.has_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed_faceon and self.has_faceon_observed_diffuse_dust_sed: return self.faceon_observed_sed - self.faceon_observed_diffuse_dust_sed

        elif self.has_full_sed_faceon and self.has_faceon_observed_dust_sed: return self.faceon_observed_sed - self.faceon_observed_dust_sed + self.faceon_intrinsic_dust_sed

        else: return self.get_stellar_part(self.faceon_observed_sed, full=True) + self.faceon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed_faceon and self.has_faceon_observed_diffuse_dust_sed) or (self.has_full_sed_faceon and self.has_faceon_observed_dust_sed and self.has_faceon_intrinsic_dust_sed) or self.has_faceon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed_edgeon and self.has_edgeon_observed_diffuse_dust_sed: return self.edgeon_observed_sed - self.edgeon_observed_diffuse_dust_sed

        elif self.has_full_sed_edgeon and self.has_edgeon_observed_dust_sed: return self.edgeon_observed_sed - self.edgeon_observed_dust_sed + self.edgeon_intrinsic_dust_sed

        else: return self.get_stellar_part(self.edgeon_observed_sed, full=True) + self.edgeon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed_edgeon and self.has_edgeon_observed_diffuse_dust_sed) or (self.has_full_sed_edgeon and self.has_edgeon_observed_dust_sed and self.has_edgeon_intrinsic_dust_sed) or self.has_edgeon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        # Has full instrument and diffuse dust cube
        # stellar = observed (=> BUT THIS IS WITH MAPPINGS INTERNAL DUST EM) - dust DIFFUSE
        # => GOOD, EMITTED MAPPINGS IS ALSO ABSORBED INTRINSIC MAPPINGS STELLAR [ENERGY BALANCE IN MAPPINGS])
        #         = direct + scattered
        if self.has_full_cube and self.has_observed_diffuse_dust_cube: return self.observed_cube - self.observed_diffuse_dust_cube

        # Has full instrument, and has intrinsic dust cube
        elif self.has_full_cube and self.has_observed_dust_cube: return self.observed_cube - self.observed_dust_cube + self.intrinsic_dust_cube

        # No full instrument data: get the stellar part of the total spectral cube
        # + INTRINSIC INTERNAL DUST CUBE -> ASSUME THIS RADIATION IS NOT REPROCESSED
        else: return self._get_observed_stellar_cube_from_stellar_part_and_intrinsic_dust()

    # -----------------------------------------------------------------

    def _get_observed_stellar_cube_from_stellar_part_and_intrinsic_dust(self):

        """
        This function ...
        :return: 
        """

        stellar, intrinsic_dust = uniformize(self.get_stellar_part(self.observed_cube, full=True), self.intrinsic_dust_cube)
        return stellar + intrinsic_dust

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_cube and self.has_observed_diffuse_dust_cube) or (self.has_full_cube and self.has_observed_diffuse_dust_cube and self.has_intrinsic_dust_cube) or (self.has_intrinsic_dust_cube)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        # Has full instrument and diffuse dust cube
        if self.has_full_cube_faceon and self.has_faceon_observed_diffuse_dust_cube: return self.faceon_observed_cube - self.faceon_observed_diffuse_dust_cube

        # Has full instrument
        elif self.has_full_cube_faceon and self.has_faceon_observed_dust_cube: return self.faceon_observed_cube - self.faceon_observed_dust_cube + self.faceon_intrinsic_dust_cube

        # No full instrument data
        else: return self._get_faceon_observed_stellar_cube_from_stellar_part_and_intrinsic_dust()

    # -----------------------------------------------------------------

    def _get_faceon_observed_stellar_cube_from_stellar_part_and_intrinsic_dust(self):

        """
        This function ...
        :return:
        """

        stellar, dust = uniformize(self.get_stellar_part(self.faceon_observed_cube, full=True), self.faceon_intrinsic_dust_cube)
        return stellar + dust

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_cube_faceon and self.has_faceon_observed_diffuse_dust_cube) or (self.has_full_cube_faceon and self.has_faceon_observed_dust_cube and self.has_faceon_intrinsic_dust_cube) or self.has_faceon_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_stellar_cube(self):

        """
        This function ...
        :return:
        """

        # Has full instrument and diffuse dust cube
        if self.has_full_cube_edgeon and self.has_edgeon_observed_diffuse_dust_cube: return self.edgeon_observed_cube - self.edgeon_observed_diffuse_dust_cube

        # Has full instrument
        elif self.has_full_cube_edgeon and self.has_edgeon_observed_dust_cube: return self.edgeon_observed_cube - self.edgeon_observed_dust_cube + self.edgeon_intrinsic_dust_cube

        # No full instrument data
        else: return self._get_edgeon_observed_stellar_cube_from_stellar_part_and_intrinsic_dust()

    # -----------------------------------------------------------------

    def _get_edgeon_observed_stellar_cube_from_stellar_part_and_intrinsic_dust(self):

        """
        This function ...
        :return:
        """

        stellar, dust = uniformize(self.get_stellar_part(self.edgeon_observed_cube, full=True), self.edgeon_intrinsic_dust_cube)
        return stellar + dust

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_stellar_cube(self):

        """
        Thisn function ...
        :return:
        """

        return (self.has_full_cube_edgeon and self.has_edgeon_observed_diffuse_dust_cube) or (self.has_full_cube_edgeon and self.has_edgeon_observed_dust_cube and self.has_edgeon_intrinsic_dust_cube) or self.has_edgeon_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @abstractproperty
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def intrinsic_stellar_sed(self):

        """
        Thisf unction ...
        :return:
        """

        #return self.intrinsic_sed.splice(max_wavelength=stellar_dust_sed_split_wavelength)
        #return self.get_stellar_part(self.intrinsic_sed)
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        #return self.intrinsic_cube.splice(max_wavelength=stellar_dust_sed_split_wavelength)
        #return self.get_stellar_part(self.intrinsic_cube)
        return self.intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.faceon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_faceon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube_faceon

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_cube_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_edgeon_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_bolometric_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_stellar_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_stellar_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed: return self.observed_sed_dust + self.intrinsic_dust_sed
        else: return self.get_dust_part(self.observed_sed, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed: return self.observed_sed_dust
        else: return self.get_dust_part(self.observed_sed, full=True) - self.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube: return self.observed_cube_dust + self.intrinsic_dust_cube
        else: return self.get_dust_part(self.observed_cube, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube: return self.observed_cube_dust
        else: return self.get_dust_part(self.observed_cube, full=True) - self.intrinsic_dust_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed_faceon: return self.faceon_observed_sed_dust + self.faceon_intrinsic_dust_sed
        else: return self.get_dust_part(self.faceon_observed_sed, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed_faceon: return self.faceon_observed_sed_dust
        else: return self.get_dust_part(self.faceon_observed_sed, full=True) - self.faceon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube_faceon: return self.faceon_observed_cube_dust + self.faceon_intrinsic_dust_cube
        else: return self.get_dust_part(self.faceon_observed_cube, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube_faceon: return self.faceon_observed_cube_dust
        else: return self.get_dust_part(self.faceon_observed_cube, full=True) - self.faceon_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_full_sed_edgeon: return self.edgeon_observed_sed_dust + self.edgeon_intrinsic_dust_sed
        else: return self.get_dust_part(self.edgeon_observed_sed, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_diffuse_dust_sed(self):

        """
        Thisf unction ...
        :return:
        """

        if self.has_full_sed_edgeon: return self.edgeon_observed_sed_dust
        else: return self.get_dust_part(self.edgeon_observed_sed, full=True) - self.edgeon_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube_edgeon: return self.edgeon_observed_cube_dust + self.edgeon_intrinsic_dust_cube
        else: return self.get_dust_part(self.edgeon_observed_cube, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_full_cube_edgeon: return self.edgeon_observed_cube_dust
        else: return self.get_dust_part(self.edgeon_observed_cube, full=True) - self.edgeon_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed and self.has_intrinsic_dust_sed) or self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed_faceon and self.has_intrinsic_dust_sed) or self.has_faceon_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_dust_sed(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_sed_edgeon and self.has_intrinsic_dust_sed) or self.has_edgeon_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_full_sed or (self.has_observed_sed and self.has_intrinsic_dust_sed)

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_full_sed_faceon or (self.has_faceon_observed_sed and self.has_intrinsic_dust_sed)

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_diffuse_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_full_sed_edgeon or (self.has_edgeon_observed_sed and self.has_intrinsic_dust_sed)

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_cube and self.has_intrinsic_dust_cube) or self.has_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_cube_faceon and self.has_faceon_intrinsic_dust_cube) or self.has_faceon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_dust_cube(self):

        """
        This function ...
        :return:
        """

        return (self.has_full_cube_edgeon and self.has_edgeon_intrinsic_dust_cube) or self.has_edgeon_observed_cube

    # -----------------------------------------------------------------

    @property
    def has_observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_full_cube or (self.has_observed_cube and self.has_intrinsic_dust_cube)

    # -----------------------------------------------------------------

    @property
    def has_faceon_observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_full_cube_faceon or (self.has_faceon_observed_cube and self.has_faceon_intrinsic_dust_cube)

    # -----------------------------------------------------------------

    @property
    def has_edgeon_observed_diffuse_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_full_cube_edgeon or (self.has_edgeon_observed_cube and self.has_edgeon_intrinsic_dust_cube)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_part(self.intrinsic_sed, full=True)

    # ----------------------------------------------------------------

    @property
    def faceon_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_sed # INTRINSIC = ISOTROPIC

    # ----------------------------------------------------------------

    @property
    def edgeon_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_sed # INTRINSIC = ISOTROPIC

    # ----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_part(self.intrinsic_cube, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.get_dust_part(self.faceon_intrinsic_cube, full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_dust_cube(self):

        """
        Thisfunctino ...
        :return:
        """

        return self.get_dust_part(self.edgeon_intrinsic_cube, full=True)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_faceon_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_sed # INTRINSIC = ISOTROPIC

    # -----------------------------------------------------------------

    @property
    def has_edgeon_intrinsic_dust_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_sed  # INTRINSIC = ISOTROPIC

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @property
    def has_faceon_intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_edgeon_intrinsic_dust_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.observed_dust_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.observed_dust_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_dust_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_sed.integrate()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_dust_cube.integrate()

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_dust_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_dust_cube

    # -----------------------------------------------------------------

    @memoize_method
    def observed_photometry_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        return self.observed_sed.photometry_at(wavelength, interpolate=interpolate)

    # -----------------------------------------------------------------

    @memoize_method
    def observed_frame_at(self, wavelength, interpolate=False):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        if interpolate: raise NotImplementedError("Interpolating is not supported")
        return self.observed_cube.get_frame_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    @property
    def has_observed_photometry(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed

    # -----------------------------------------------------------------

    @property
    def has_observed_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube

    # -----------------------------------------------------------------

    @memoize_method
    def intrinsic_photometry_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        return self.intrinsic_sed.photometry_at(wavelength, interpolate=interpolate)

    # -----------------------------------------------------------------

    @memoize_method
    def intrinsic_frame_at(self, wavelength, interpolate=False):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        if interpolate: raise NotImplementedError("Interpolating is not supported")
        return self.intrinsic_cube.get_frame_for_wavelength(wavelength)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_photometry(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_cube

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_curve(self):

        """
        This function ...
        :return:
        """

        observed = self.observed_stellar_sed
        intrinsic = self.intrinsic_stellar_sed
        return AttenuationCurve.from_seds(observed, intrinsic)

    # -----------------------------------------------------------------

    @lazyproperty
    def attenuation_cube(self):

        """
        This function ...
        :return:
        """

        # Get the datacubes as arrays
        observed = self.observed_stellar_cube.asarray() # wavelength axis first
        intrinsic = self.intrinsic_stellar_cube.asarray()  # wavelength axis first

        # Calculate the attenuations
        attenuation = extinction.attenuation(observed, intrinsic)

        # Create and return the datacube of attenuations
        return DataCube.from_array(attenuation, wavelength_grid=self.wavelength_grid)

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        # Get the observed stellar luminosity
        observed = self.observed_stellar_luminosity.to("W", distance=self.distance).value

        # Get the intrinsic stellar luminosity
        #intrinsic = self.intrinsic_stellar_luminosity.to("W", distance=self.distance).value
        # includes intrinsic dust emission in template, so also ACTUAL intrinsic bolometric lum
        intrinsic = self.intrinsic_bolometric_luminosity.to("W", distance=self.distance).value

        # Calculate and return the extinction
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_attenuation_frame(self):

        """
        This function ...
        :return:
        """

        # Get the observed stellar frame
        observed = self.observed_stellar_frame
        observed.convert_to("W", distance=self.distance)

        # Get the intrinsic frame
        intrinsic = self.intrinsic_bolometric_frame
        intrinsic.convert_to("W", distance=self.distance)

        # Calculate and return the extinction
        attenuation = extinction.attenuation(observed.data, intrinsic.data)
        return Frame(attenuation, wcs=observed.wcs)

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_stellar_luminosity and self.has_intrinsic_bolometric_luminosity #self.has_intrinsic_stellar_luminosity

    # -----------------------------------------------------------------

    @property
    def has_bolometric_attenuation_frame(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_stellar_frame

    # -----------------------------------------------------------------

    @memoize_method
    def attenuation_at(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        # Get the observed luminosity
        observed = self.observed_photometry_at(wavelength, interpolate=interpolate).to("W/micron", wavelength=wavelength, distance=self.distance).value

        # Get the intrinsic luminosity
        intrinsic = self.intrinsic_photometry_at(wavelength, interpolate=interpolate).to("W/micron", wavelength=wavelength, distance=self.distance).value

        # Calculate and return the extinction
        return extinction.attenuation(observed, intrinsic)

    # -----------------------------------------------------------------

    @memoize_method
    def attenuation_at_frame(self, wavelength, interpolate=True):

        """
        This function ...
        :param wavelength:
        :param interpolate:
        :return:
        """

        # Get the observed frame
        observed = self.observed_frame_at(wavelength, interpolate=interpolate)
        observed.convert_to("W/micron", wavelength=wavelength, distance=self.distance)

        # Get the intrinsic frame
        intrinsic = self.intrinsic_frame_at(wavelength, interpolate=interpolate)
        intrinsic.convert_to("W/micron", wavelength=wavelength, distance=self.distance)

        # Calculate and return the extinction
        attenuation = extinction.attenuation(observed.data, intrinsic.data)
        return Frame(attenuation, wcs=observed.wcs, unit="W/micron", wavelength=wavelength)

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_sed and self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_attenuation_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_observed_cube and self.has_intrinsic_cube

# -----------------------------------------------------------------
