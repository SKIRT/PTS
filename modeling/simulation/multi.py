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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import sequences
from ...core.data.sed import load_sed
from ...magic.core.datacube import DataCube
from .simulations import ComponentSimulations
from .simulation import ObservedComponentSimulation

# -----------------------------------------------------------------

class MultiComponentSimulations(ComponentSimulations):

    """
    Objects of this class describe the simulation(s) of a radiative transfer model containing of multiple stellar components.
    """

    def __init__(self, name, observed, intrinsic_seds=None, intrinsic_cubes=None, intrinsic_cubes_faceon=None,
                 intrinsic_cubes_edgeon=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_seds:
        :param intrinsic_cubes:
        :param intrinsic_cubes_faceon:
        :param intrinsic_cubes_edgeon:
        :param distance:
        """

        # Call the constructor of the base class
        super(MultiComponentSimulations, self).__init__(name, observed, distance=distance)

        # Set the SEDs of the components
        self.intrinsic_seds = intrinsic_seds

        # Set the datacubes of the components
        self.intrinsic_cubes = intrinsic_cubes
        self.intrinsic_cubes_faceon = intrinsic_cubes_faceon
        self.intrinsic_cubes_edgeon = intrinsic_cubes_edgeon

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, name, observed, intrinsic_sed_paths=None, intrinsic_cube_paths=None,
                         intrinsic_cube_faceon_paths=None, intrinsic_cube_edgeon_paths=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic_sed_paths:
        :param intrinsic_cube_paths:
        :param intrinsic_cube_faceon_paths:
        :param intrinsic_cube_edgeon_paths:
        :param distance:
        :return:
        """

        # Load observed simulation
        observed = ObservedComponentSimulation.from_output_path(observed)

        # Load intrinsic SEDs
        if intrinsic_sed_paths is not None:
            intrinsic_seds = OrderedDict()
            for component_name in intrinsic_sed_paths: intrinsic_seds[component_name] = load_sed(intrinsic_sed_paths[component_name])
        else: intrinsic_seds = None

        # Load intrinsic cubes
        if intrinsic_cube_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            intrinsic_cubes = OrderedDict()
            for component_name in intrinsic_cube_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid from SED
                intrinsic_cubes[component_name] = DataCube.from_file(intrinsic_cube_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes = None

        # Load intrinsic faceon cubes
        if intrinsic_cube_faceon_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            intrinsic_cubes_faceon = OrderedDict()
            for component_name in intrinsic_cube_faceon_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid
                intrinsic_cubes_faceon[component_name] = DataCube.from_file(intrinsic_cube_faceon_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes_faceon = None

        # Load intrinsic edgeon cubes
        if intrinsic_cube_edgeon_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            intrinsic_cubes_edgeon = OrderedDict()
            for component_name in intrinsic_cube_edgeon_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid
                intrinsic_cubes_edgeon[component_name] = DataCube.from_file(intrinsic_cube_edgeon_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes_edgeon = None

        # Create and return
        return cls(name, observed, intrinsic_seds=intrinsic_seds, intrinsic_cubes=intrinsic_cubes,
                   intrinsic_cubes_faceon=intrinsic_cubes_faceon, intrinsic_cubes_edgeon=intrinsic_cubes_edgeon,
                   distance=distance)

    # -----------------------------------------------------------------

    @property
    def component_names(self):

        """
        This function ...
        :return:
        """

        if self.has_intrinsic_seds: return self.intrinsic_seds.keys()
        elif self.has_intrinsic_cubes: return self.intrinsic_cubes.keys()
        elif self.has_intrinsic_cubes_faceon: return self.intrinsic_cubes_faceon.keys()
        elif self.has_intrinsic_cubes_edgeon: return self.intrinsic_cubes_edgeon.keys()
        else: raise ValueError("Component names are not defined")

    # -----------------------------------------------------------------

    @property
    def nintrinsic_seds(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_seds)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_cubes)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        return len(self.intrinsic_cubes_faceon)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes_edgeon(self):

        """
        Thisfunction ...
        :return:
        """

        return len(self.intrinsic_cubes_edgeon)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_seds(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_seds is not None and self.nintrinsic_seds > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cubes is not None and self.nintrinsic_cubes > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cubes_faceon is not None and self.nintrinsic_cubes_faceon > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_cubes_edgeon is not None and self.nintrinsic_cubes_edgeon > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_sed: return True
        elif self.has_intrinsic_seds: return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_cube: return True
        elif self.has_intrinsic_cubes: return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_faceon(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_cube_faceon: return True
        elif self.has_intrinsic_cubes_faceon: return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_edgeon(self):

        """
        This function ...
        :return:
        """

        if self.has_transparent_cube_edgeon: return True
        elif self.has_intrinsic_cubes_edgeon: return True
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # Transparent SED is written out
        if self.has_transparent_sed: return self.observed_sed_transparent

        # Has intrinsic SEDs, add them
        elif self.has_intrinsic_seds: return sequences.sum(self.intrinsic_seds.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic SED cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube: return self.observed_cube_transparent

        # Has intrinsic cubes, add them
        elif self.has_intrinsic_cubes: return sequences.sum(self.intrinsic_cubes.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic datacube cannot be calculated")

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        #if self.has_faceon_transparent_cube: return self.faceon_observed_cube_transparent

        # Cannot be calculated
        #else: raise ValueError("Intrinsic datacube from face-on orientation cannot be calculated")

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        #if self.has_edgeon_transparent_cube: return self.edgeon_observed_cube_transparent

        # Cannot be calculated
        #else: raise ValueError("Intrinsic datacube from edge-on orientation cannot be calculated")

        raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------
