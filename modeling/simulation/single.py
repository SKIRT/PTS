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

# -----------------------------------------------------------------

class SingleComponentSimulations(ComponentSimulations):
    
    """
    Objects of this class describe the simulation(s) of radiative transfer model of a certain stellar component.
    """

    def __init__(self, name, observed, intrinsic=None, distance=None):

        """
        This function ...
        :param name:
        :param intrinsic:
        :param observed:
        :param distance:
        """

        # Call the constructor of the base class
        super(SingleComponentSimulations, self).__init__(name, observed, distance=distance)

        # Set the intrinsic simulation
        self.intrinsic = intrinsic

    # -----------------------------------------------------------------

    @classmethod
    def from_output_paths(cls, name, observed, intrinsic=None, distance=None):

        """
        This function ...
        :param name:
        :param observed:
        :param intrinsic:
        :param distance:
        :return:
        """

        # Load observed simulation
        if not fs.is_directory(observed): raise ValueError("Observed simulation directory does not exist")
        if fs.has_files_in_path(observed): observed = ObservedComponentSimulation.from_output_path(observed)
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

        # Create and return
        return cls(name, observed, intrinsic=intrinsic, distance=distance)

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
    def has_intrinsic_sed(self):

        """
        This function ..
        :return:
        """

        return self.has_transparent_sed

    # -----------------------------------------------------------------

    @property
    def has_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        return self.has_transparent_cube

    # -----------------------------------------------------------------

    @property
    def intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # Check whether transparent SED is created
        if not self.has_transparent_sed: raise ValueError("Intrinsic SED cannot be calculated")

        # Return
        return self.observed_sed_transparent

    # -----------------------------------------------------------------

    @property
    def intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check whether transparant SED is created
        if not self.has_transparent_cube: raise ValueError("Intrinsic cube cannot be calculated")

        # Return
        return self.observed_cube_transparent

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def faceon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check
        #if not self.has_faceon_transparent_cube: raise ValueError("Intrinsic cube from face-on orientation cannot be calculated")

        # Return
        #return self.faceon_observed_cube_transparent

        raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_sed(self):

        """
        This function ...
        :return:
        """

        # ISOTROPIC RADIATION
        return self.intrinsic_sed

    # -----------------------------------------------------------------

    @property
    def edgeon_intrinsic_cube(self):

        """
        This function ...
        :return:
        """

        # Check
        #if not self.has_edgeon_transparent_cube: raise ValueError("Intrinsic cube from edge-on orientation cannot be calculated")

        # Return
        #return self.edgeon_observed_cube_transparent

        raise NotImplementedError("Not implemented")

# -----------------------------------------------------------------
