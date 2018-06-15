#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.model Contains the IntrinsicComponentSimulation, ObservedComponentSimulation,
#  and ComponentSimulation class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from .simulations import ComponentSimulations
from .simulation import ObservedComponentSimulation, IntrinsicComponentSimulation
from .simulation import earth_name, faceon_name, edgeon_name, total_contribution
from ...magic.core.frame import Frame
from ...magic.core.datacube import DataCube
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class SingleComponentSimulations(ComponentSimulations):
    
    """
    Objects of this class describe the simulation(s) of radiative transfer model of a certain stellar component.
    """

    def __init__(self, name, observed, intrinsic=None, distance=None, map_earth=None, map_faceon=None, map_edgeon=None,
                 earth_wcs=None):

        """
        This function ...
        :param name:
        :param intrinsic:
        :param observed:
        :param distance:
        :param map_earth:
        :param map_faceon:
        :param map_edgeon:
        :param earth_wcs:
        """

        # Call the constructor of the base class
        super(SingleComponentSimulations, self).__init__(name, observed, distance=distance, earth_wcs=earth_wcs)

        # Set the intrinsic simulation
        self.intrinsic = intrinsic

        # Set the maps
        self.map_earth = map_earth
        self.map_faceon = map_faceon
        self.map_edgeon = map_edgeon

    # -----------------------------------------------------------------

    @classmethod
    def from_output_paths(cls, name, observed, intrinsic=None, distance=None, map_earth_path=None, map_earth=None,
                          map_faceon_path=None, map_faceon=None, map_edgeon_path=None, map_edgeon=None, earth_wcs=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic:
        :param distance:
        :param map_earth_path:
        :param map_earth:
        :param map_faceon_path:
        :param map_faceon:
        :param map_edgeon_path:
        :param map_edgeon:
        :param earth_wcs:
        :return:
        """

        # Load observed simulation
        if not fs.is_directory(observed): raise ValueError("Observed simulation directory does not exist")
        if fs.has_files_in_path(observed): observed = ObservedComponentSimulation.from_output_path(observed, earth_wcs=earth_wcs)
        else:
            log.warning("Observed simulation has not been performed: no output files")
            observed = None

        # Load intrinsic simulation
        if intrinsic is not None:
            if not fs.is_directory(intrinsic): raise ValueError("Intrinsic simulation directory does not exist")
            if fs.has_files_in_path(intrinsic): intrinsic = IntrinsicComponentSimulation.from_output_path(intrinsic)
            else:
                log.warning("Intrinsic simulation has not been performed: no output files")
                intrinsic = None

        # Load earth map
        if map_earth is not None:
            if map_earth_path is not None: raise ValueError("Cannot specify both map_earth and map_earth_path")
        elif map_earth_path is not None: map_earth = Frame.from_file(map_earth_path, wcs=earth_wcs)
        else: map_earth = None

        # Load faceon map
        if map_faceon is not None:
            if map_faceon_path is not None: raise ValueError("Cannot specify both map_faceon and map_faceon_path")
        elif map_faceon_path is not None: map_faceon = Frame.from_file(map_faceon_path)
        else: map_faceon = None

        # Load edgeon map
        if map_edgeon is not None:
            if map_edgeon_path is not None: raise ValueError("Cannot specify both map_edgeon and map_edgeon_path")
        elif map_edgeon_path is not None: map_edgeon = Frame.from_file(map_edgeon_path)
        else: map_edgeon = None

        # Create and return
        return cls(name, observed, intrinsic=intrinsic, distance=distance, map_earth=map_earth, map_faceon=map_faceon,
                   map_edgeon=map_edgeon, earth_wcs=earth_wcs)

    # -----------------------------------------------------------------

    @property
    def has_map_earth(self):

        """
        This function ...
        :return:
        """

        return self.map_earth is not None

    # -----------------------------------------------------------------

    @property
    def has_map_faceon(self):

        """
        This function ...
        :return:
        """

        return self.map_faceon is not None

    # -----------------------------------------------------------------

    @property
    def has_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.map_edgeon is not None

    # -----------------------------------------------------------------

    @property
    def has_intrinsic(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic is not None

    # -----------------------------------------------------------------

    @property
    def intrinsic_output_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.output_path if self.has_intrinsic else None

    # -----------------------------------------------------------------

    @property
    def intrinsic_output(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.output

    # -----------------------------------------------------------------

    @property
    def intrinsic_data(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic.data

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_output(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic and self.intrinsic.has_output

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_data(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_output and self.intrinsic_data.has_any

    # -----------------------------------------------------------------

    @property
    def has_sed_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic and self.has_intrinsic_output and self.intrinsic_data.has_seds

    # -----------------------------------------------------------------

    @property
    def has_cube_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic and self.has_intrinsic_output and self.intrinsic_data.has_images

    # -----------------------------------------------------------------

    @property
    def has_faceon_cube_from_intrinsic_simulation(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_intrinsic and self.has_intrinsic_output and self.has_faceon_intrinsic_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_edgeon_cube_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic and self.has_intrinsic_output and self.has_edgeon_intrinsic_cube_orientation

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        return self.has_transparent_sed or self.has_sed_from_intrinsic_simulation

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_transparent_cube or self.has_cube_from_intrinsic_simulation or (self.has_map_earth and self.has_intrinsic_sed)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_transparent_cube_faceon or self.has_faceon_cube_from_intrinsic_simulation or (self.has_map_faceon and self.has_intrinsic_sed)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_edgeon(self):

        """
        This fnuction ...
        :return:
        """

        return self.has_transparent_cube_edgeon or self.has_edgeon_cube_from_intrinsic_simulation or (self.has_map_edgeon and self.has_intrinsic_sed)

    # -----------------------------------------------------------------
    # INTRINSIC SEDs
    # -----------------------------------------------------------------

    @property
    def has_other_intrinsic_sed_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_data.sed_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_intrinsic_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.intrinsic_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_intrinsic_sed_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.intrinsic_data.sed_paths_instruments

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.sed_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.sed_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_sed_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.sed_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def sed_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.seds[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_sed_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.seds[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_sed_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.seds[edgeon_name][total_contribution]

    # -----------------------------------------------------------------
    # INTRINSIC CUBES
    # -----------------------------------------------------------------

    @property
    def has_other_intrinsic_cube_orientations(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_data.image_paths_instruments) > 1

    # -----------------------------------------------------------------

    @property
    def has_faceon_intrinsic_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return faceon_name in self.intrinsic_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def has_edgeon_intrinsic_cube_orientation(self):

        """
        This function ...
        :return:
        """

        return edgeon_name in self.intrinsic_data.image_paths_instruments

    # -----------------------------------------------------------------

    @property
    def intrinsic_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.image_paths_instruments[earth_name]

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.image_paths_instruments[faceon_name]

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_cube_path(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.image_paths_instruments[edgeon_name]

    # -----------------------------------------------------------------

    @property
    def cube_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.images[earth_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def faceon_cube_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.images[faceon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def edgeon_cube_from_intrinsic_simulation(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_data.images[edgeon_name][total_contribution]

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # Transparent SED is written out
        if self.has_transparent_sed: return self.observed_sed_transparent

        # Has intrinsic simulation
        elif self.has_sed_from_intrinsic_simulation: return self.sed_from_intrinsic_simulation

        # Cannot be calculated
        else: raise ValueError("Intrinsic SED cannot be calculated")

    # -----------------------------------------------------------------

    @property
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube: return self.observed_cube_transparent

        # Has intrinsic simulation
        elif self.has_cube_from_intrinsic_simulation: return self.cube_from_intrinsic_simulation

        # Has map and intrinsic SED
        elif self.has_map_earth and self.has_intrinsic_sed: return self.cube_from_map_and_intrinsic_sed

        # Cannot be calculated
        else: raise ValueError("Intrinsic cube cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def cube_from_map_and_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        return DataCube.from_sed_and_map(self.intrinsic_sed, self.map_earth)

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION, SO SAME AS EARTH PROJECTION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube_faceon: return self.faceon_observed_cube_transparent

        # Has intrinsic simulation
        elif self.has_faceon_cube_from_intrinsic_simulation: return self.faceon_cube_from_intrinsic_simulation

        # Has map and intrinsic SED
        elif self.has_map_faceon and self.has_intrinsic_sed: return self.faceon_cube_from_map_and_intrinsic_sed

        # Cannot be calculated
        else: raise ValueError("Intrinsic cube cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_cube_from_map_and_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        return DataCube.from_sed_and_map(self.intrinsic_sed, self.map_faceon)

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION, SO SAME AS EARTH PROJECTION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube_edgeon: return self.edgeon_observed_cube_transparent

        # Has intrinsic simulation
        elif self.has_edgeon_cube_from_intrinsic_simulation: return self.edgeon_cube_from_intrinsic_simulation

        # Has map and intrinsic SED
        elif self.has_map_edgeon and self.has_intrinsic_sed: return self.edgeon_cube_from_map_and_intrinsic_sed

        # Cannot be calculated
        else: raise ValueError("Intrinsic cube cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_cube_from_map_and_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        return DataCube.from_sed_and_map(self.intrinsic_sed, self.map_edgeon)

# -----------------------------------------------------------------
