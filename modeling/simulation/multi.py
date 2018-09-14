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
from ...magic.core.list import uniformize
from ...core.tools.utils import LazyDictionary

# -----------------------------------------------------------------

class MultiComponentSimulations(ComponentSimulations):

    """
    Objects of this class describe the simulation(s) of a radiative transfer model containing of multiple stellar components.
    """

    def __init__(self, name, observed, component_simulations, intrinsic_seds=None, intrinsic_cubes=None, intrinsic_cubes_faceon=None,
                 intrinsic_cubes_edgeon=None, distance=None, earth_wcs=None):

        """
        This function ...
        :param name:
        :param observed:
        :param component_simulations:
        :param intrinsic_seds:
        :param intrinsic_cubes:
        :param intrinsic_cubes_faceon:
        :param intrinsic_cubes_edgeon:
        :param distance:
        :param earth_wcs:
        """

        # Call the constructor of the base class
        super(MultiComponentSimulations, self).__init__(name, observed, distance=distance, earth_wcs=earth_wcs)

        # The component simulations
        self.component_simulations = component_simulations

        # Set the SEDs of the components, if passed
        if intrinsic_seds is not None: self.intrinsic_seds = intrinsic_seds

        # Set the datacubes of the components, if passed
        if intrinsic_cubes is not None: self.intrinsic_cubes = intrinsic_cubes
        if intrinsic_cubes_faceon is not None: self.intrinsic_cubes_faceon = intrinsic_cubes_faceon
        if intrinsic_cubes_edgeon is not None: self.intrinsic_cubes_edgeon = intrinsic_cubes_edgeon

    # -----------------------------------------------------------------

    @classmethod
    def from_output_path(cls, name, observed_path, component_simulations, intrinsic_sed_paths=None, intrinsic_seds=None,
                         intrinsic_cube_paths=None, intrinsic_cubes=None, intrinsic_cube_faceon_paths=None,
                         intrinsic_cubes_faceon=None, intrinsic_cube_edgeon_paths=None, intrinsic_cubes_edgeon=None, distance=None, earth_wcs=None):

        """
        This function ...
        :param name:
        :param observed_path:
        :param component_simulations:
        :param intrinsic_sed_paths:
        :param intrinsic_seds:
        :param intrinsic_cube_paths:
        :param intrinsic_cubes:
        :param intrinsic_cube_faceon_paths:
        :param intrinsic_cubes_faceon:
        :param intrinsic_cube_edgeon_paths:
        :param intrinsic_cubes_edgeon:
        :param distance:
        :param earth_wcs:
        :return:
        """

        # Load observed simulation
        observed = ObservedComponentSimulation.from_output_path(observed_path, earth_wcs=earth_wcs)

        # Create and return
        return cls.from_observed(name, observed, component_simulations=component_simulations, intrinsic_sed_paths=intrinsic_sed_paths, intrinsic_seds=intrinsic_seds,
                         intrinsic_cube_paths=intrinsic_cube_paths, intrinsic_cubes=intrinsic_cubes, intrinsic_cube_faceon_paths=intrinsic_cube_faceon_paths,
                         intrinsic_cubes_faceon=intrinsic_cubes_faceon, intrinsic_cube_edgeon_paths=intrinsic_cube_edgeon_paths, intrinsic_cubes_edgeon=intrinsic_cubes_edgeon,
                         distance=distance, earth_wcs=earth_wcs)

    # -----------------------------------------------------------------

    @classmethod
    def from_output(cls, name, observed_output, component_simulations, intrinsic_sed_paths=None, intrinsic_seds=None,
                     intrinsic_cube_paths=None, intrinsic_cubes=None, intrinsic_cube_faceon_paths=None,
                     intrinsic_cubes_faceon=None, intrinsic_cube_edgeon_paths=None, intrinsic_cubes_edgeon=None, distance=None, earth_wcs=None):

        """
        This function ...
        :param name:
        :param observed_output:
        :param component_simulations:
        :param intrinsic_sed_paths:
        :param intrinsic_seds:
        :param intrinsic_cube_paths:
        :param intrinsic_cubes:
        :param intrinsic_cube_faceon_paths:
        :param intrinsic_cubes_faceon:
        :param intrinsic_cube_edgeon_paths:
        :param intrinsic_cubes_edgeon:
        :param distance:
        :param earth_wcs:
        :return:
        """

        # Load observed simulation
        observed = ObservedComponentSimulation.from_output(observed_output, earth_wcs=earth_wcs)

        # Create and return
        return cls.from_observed(name, observed, component_simulations=component_simulations, intrinsic_sed_paths=intrinsic_sed_paths, intrinsic_seds=intrinsic_seds,
                         intrinsic_cube_paths=intrinsic_cube_paths, intrinsic_cubes=intrinsic_cubes, intrinsic_cube_faceon_paths=intrinsic_cube_faceon_paths,
                         intrinsic_cubes_faceon=intrinsic_cubes_faceon, intrinsic_cube_edgeon_paths=intrinsic_cube_edgeon_paths, intrinsic_cubes_edgeon=intrinsic_cubes_edgeon,
                         distance=distance, earth_wcs=earth_wcs)

    # -----------------------------------------------------------------

    @classmethod
    def from_observed(cls, name, observed, component_simulations, intrinsic_sed_paths=None, intrinsic_seds=None,
                     intrinsic_cube_paths=None, intrinsic_cubes=None, intrinsic_cube_faceon_paths=None,
                     intrinsic_cubes_faceon=None, intrinsic_cube_edgeon_paths=None, intrinsic_cubes_edgeon=None, distance=None, earth_wcs=None):

        """
        This function ...
        :param name:
        :param observed:
        :param component_simulations:
        :param intrinsic_sed_paths:
        :param intrinsic_seds:
        :param intrinsic_cube_paths:
        :param intrinsic_cubes:
        :param intrinsic_cube_faceon_paths:
        :param intrinsic_cubes_faceon:
        :param intrinsic_cube_edgeon_paths:
        :param intrinsic_cubes_edgeon:
        :param distance:
        :param earth_wcs:
        :return:
        """

        # Load intrinsic SEDs
        if intrinsic_seds is not None:
            if intrinsic_sed_paths is not None: raise ValueError("Cannot specify both intrinsic SEDs and intrinsic SED paths")
        if intrinsic_sed_paths is not None:
            intrinsic_seds = OrderedDict()
            for component_name in intrinsic_sed_paths: intrinsic_seds[component_name] = load_sed(intrinsic_sed_paths[component_name])
        else: intrinsic_seds = None

        # Load intrinsic cubes (EARTH)
        if intrinsic_cubes is not None:
            if intrinsic_cube_paths is not None: raise ValueError("Cannot specify both intrinsic cubes and intrinsic cube paths")
        elif intrinsic_cube_paths is not None:
            #if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            #intrinsic_cubes = OrderedDict()
            intrinsic_cubes = LazyDictionary(DataCube.from_file, distance=distance, wcs=earth_wcs)
            for component_name in intrinsic_cube_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid from SED
                #intrinsic_cubes[component_name] = DataCube.from_file(intrinsic_cube_paths[component_name], wavelength_grid=wavelength_grid, wcs=earth_wcs)
                intrinsic_cubes.set(component_name, intrinsic_cube_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes = None

        # Load intrinsic faceon cubes
        if intrinsic_cubes_faceon is not None:
            if intrinsic_cube_faceon_paths is not None: raise ValueError("Cannot specify both intrinsic cubes and intrinsic cube paths")
        elif intrinsic_cube_faceon_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            #intrinsic_cubes_faceon = OrderedDict()
            intrinsic_cubes_faceon = LazyDictionary(DataCube.from_file, distance=distance)
            for component_name in intrinsic_cube_faceon_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid
                #intrinsic_cubes_faceon[component_name] = DataCube.from_file(intrinsic_cube_faceon_paths[component_name], wavelength_grid=wavelength_grid)
                intrinsic_cubes_faceon.set(component_name, intrinsic_cube_faceon_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes_faceon = None

        # Load intrinsic edgeon cubes
        if intrinsic_cubes_edgeon is not None:
            if intrinsic_cube_edgeon_paths is not None: raise ValueError("Cannot specify both intrinsic cubes and intrinsic cube paths")
        elif intrinsic_cube_edgeon_paths is not None:
            if intrinsic_seds is None: raise ValueError("When passing the filepaths of datacubes, the filepaths of corresponding simulated SEDs also have to be specified (for the wavelength grid)")
            #intrinsic_cubes_edgeon = OrderedDict()
            intrinsic_cubes_edgeon = LazyDictionary(DataCube.from_file, distance=distance)
            for component_name in intrinsic_cube_edgeon_paths:
                if component_name not in intrinsic_seds: raise ValueError("SED of component '" + component_name + "' is not loaded")
                wavelength_grid = intrinsic_seds[component_name].wavelength_grid() # create wavelength grid
                #intrinsic_cubes_edgeon[component_name] = DataCube.from_file(intrinsic_cube_edgeon_paths[component_name], wavelength_grid=wavelength_grid)
                intrinsic_cubes_edgeon.set(component_name, intrinsic_cube_edgeon_paths[component_name], wavelength_grid=wavelength_grid)
        else: intrinsic_cubes_edgeon = None

        # Create and return
        return cls(name, observed, component_simulations, intrinsic_seds=intrinsic_seds, intrinsic_cubes=intrinsic_cubes,
                   intrinsic_cubes_faceon=intrinsic_cubes_faceon, intrinsic_cubes_edgeon=intrinsic_cubes_edgeon,
                   distance=distance, earth_wcs=earth_wcs)

    # -----------------------------------------------------------------

    @property
    def component_names(self):
        #if self.has_intrinsic_seds: return self.intrinsic_seds.keys()
        #elif self.has_intrinsic_cubes: return self.intrinsic_cubes.keys()
        #elif self.has_intrinsic_cubes_faceon: return self.intrinsic_cubes_faceon.keys()
        #elif self.has_intrinsic_cubes_edgeon: return self.intrinsic_cubes_edgeon.keys()
        #else: raise ValueError("Component names are not defined")
        return self.component_simulations.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_seds(self):
        seds = OrderedDict()
        for name in self.component_names: seds[name] = self.component_simulations[name].intrinsic_sed
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_cubes(self):
        cubes = OrderedDict()
        for name in self.component_names: cubes[name] = self.component_simulations[name].intrinsic_cube
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_cubes_faceon(self):
        cubes = OrderedDict()
        for name in self.component_names: cubes[name] = self.component_simulations[name].faceon_intrinsic_cube
        return cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_cubes_edgeon(self):
        cubes = OrderedDict()
        for name in self.component_names: cubes[name] = self.component_simulations[name].edgeon_intrinsic_cube
        return cubes

    # -----------------------------------------------------------------

    @property
    def nintrinsic_seds(self):
        return len(self.intrinsic_seds)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes(self):
        return len(self.intrinsic_cubes)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes_faceon(self):
        return len(self.intrinsic_cubes_faceon)

    # -----------------------------------------------------------------

    @property
    def nintrinsic_cubes_edgeon(self):
        return len(self.intrinsic_cubes_edgeon)

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_seds(self):
        return self.intrinsic_seds is not None and self.nintrinsic_seds > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes(self):
        return self.intrinsic_cubes is not None and self.nintrinsic_cubes > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes_faceon(self):
        return self.intrinsic_cubes_faceon is not None and self.nintrinsic_cubes_faceon > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cubes_edgeon(self):
        return self.intrinsic_cubes_edgeon is not None and self.nintrinsic_cubes_edgeon > 0

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed(self):
        if self.has_transparent_sed: return True
        elif self.has_intrinsic_seds: return True
        else: return False

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_faceon(self):
        return self.has_intrinsic_sed # INTRINSIC = ISOTROPIC

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_sed_edgeon(self):
        return self.has_intrinsic_sed # INTRINSIC = ISOTROPIC

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
    def intrinsic_cubes_uniformized(self):
        return uniformize(*self.intrinsic_cubes.values())

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
        elif self.has_intrinsic_cubes: return sequences.sum(self.intrinsic_cubes_uniformized) #return sequences.sum(self.intrinsic_cubes.values())

        # Cannot be calculated
        else: raise ValueError("Intrinsic datacube cannot be calculated")

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):
        return self.has_transparent_cube or self.has_intrinsic_cubes

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_sed(self):
        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_cubes_uniformized(self):
        return uniformize(*self.intrinsic_cubes_faceon.values())

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube_faceon: return self.faceon_observed_cube_transparent

        # Has intrinsic cubes, add them
        elif self.has_intrinsic_cubes_faceon: return sequences.sum(self.faceon_intrinsic_cubes_uniformized)

        # Cannot be calculated
        else: raise ValueError("Intrinsic datacube from the face-on orientation cannot be calculated")

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_faceon(self):
        return self.has_transparent_cube_faceon or self.has_intrinsic_cubes_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_sed(self):
        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_cubes_uniformized(self):
        return uniformize(*self.intrinsic_cubes_edgeon.values())

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Transparent cube is written out
        if self.has_transparent_cube_edgeon: return self.edgeon_observed_cube_transparent

        # Has intrinsic cubes, add them
        elif self.has_intrinsic_cubes_edgeon: return sequences.sum(self.edgeon_intrinsic_cubes_uniformized)

        # Cannot be calculated
        else: raise ValueError("Intrinsic datacube from the edge-on orientation cannot be calculated")

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube_edgeon(self):
        return self.has_transparent_cube_edgeon or self.has_intrinsic_cubes_edgeon

# -----------------------------------------------------------------
