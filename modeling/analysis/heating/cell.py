#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.cell Contains the CellDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty, lazyfileproperty
from ...core.data import Data3D
from ...projection.data import project_data, project_faceon, project_edgeon
from ....core.basics.distribution import Distribution, Distribution2D
from ....magic.core.frame import Frame
from ....magic.core.image import Image
from ....core.plot.distribution import plot_distribution, plot_distributions, plot_2d_distribution
from ....magic.tools.plotting import plot_map

# -----------------------------------------------------------------

old = "old"
young = "young"
ionizing = "ionizing"
dust = "dust"
extra = "extra"
disk_components = [old, young, ionizing, extra, dust]

# -----------------------------------------------------------------

class CellDustHeatingAnalyser(DustHeatingAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(CellDustHeatingAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Writing
        self.write()

        # Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CellDustHeatingAnalyser, self).setup(**kwargs)

        # Set recreate maps
        self.set_recreate_maps()

        # Set replot maps
        self.set_replot_maps()

    # -----------------------------------------------------------------

    @property
    def do_extra(self):
        return self.has_extra

    # -----------------------------------------------------------------

    def set_recreate_maps(self):

        """
        This function ...
        :return:
        """

        # Set recreate maps flags
        if self.config.recreate_maps:
            self.config.recreate_lowres_maps = True
            self.config.recreate_highres_maps = True

        # Set recreate lowres maps flags
        if self.config.recreate_lowres_maps:
            self.config.recreate_lowres_map_midplane = True
            self.config.recreate_lowres_map_faceon = True
            self.config.recreate_lowres_map_edgeon = True

        # Set recreate highres maps flags
        if self.config.recreate_highres_maps:
            self.config.recreate_highres_map_midplane = True
            self.config.recreate_highres_map_faceon = True
            self.config.recreate_highres_map_edgeon = True

        # Lowres midplane
        if self.config.recreate_lowres_map_midplane:
            self.config.replot_lowres_map_midplane = True
            if self.has_lowres_map_midplane: fs.remove_file(self.lowres_map_midplane_path)
            if self.has_lowres_frame_midplane_interpolated: fs.remove_file(self.lowres_frame_midplane_interpolated_path)

        # Lowres faceon
        if self.config.recreate_lowres_map_faceon:
            self.config.replot_lowres_map_faceon = True
            if self.has_lowres_map_faceon: fs.remove_file(self.lowres_map_faceon_path)
            if self.has_lowres_frame_faceon_interpolated: fs.remove_file(self.lowres_frame_faceon_interpolated_path)

        # Lowres edgeon
        if self.config.recreate_lowres_map_edgeon:
            self.config.replot_lowres_map_edgeon = True
            if self.has_lowres_map_edgeon: fs.remove_file(self.lowres_map_edgeon_path)
            if self.has_lowres_frame_edgeon_interpolated: fs.remove_file(self.lowres_frame_edgeon_interpolated_path)

        # Highres midplane
        if self.config.recreate_highres_map_midplane:
            self.config.replot_highres_map_midplane = True
            if self.has_highres_map_midplane: fs.remove_file(self.highres_map_midplane_path)
            if self.has_highres_frame_midplane_interpolated: fs.remove_file(self.highres_frame_midplane_interpolated_path)

        # Highres faceon
        if self.config.recreate_highres_map_faceon:
            self.config.replot_highres_map_faceon = True
            if self.has_highres_map_faceon: fs.remove_file(self.highres_map_faceon_path)
            if self.has_highres_frame_faceon_interpolated: fs.remove_file(self.highres_frame_faceon_interpolated_path)

        # Highres edgeon
        if self.config.recreate_highres_map_edgeon:
            self.config.replot_highres_map_edgeon = True
            if self.has_highres_map_edgeon: fs.remove_file(self.highres_map_edgeon_path)
            if self.has_highres_frame_edgeon_interpolated: fs.remove_file(self.highres_frame_edgeon_interpolated_path)

    # -----------------------------------------------------------------

    def set_replot_maps(self):

        """
        This function ...
        :return:
        """

        # Set replot maps flags
        if self.config.replot_maps:
            self.config.replot_lowres_maps = True
            self.config.replot_highres_maps = True

        # Set replot lowres maps flags
        if self.config.replot_lowres_maps:
            self.config.replot_lowres_map_midplane = True
            self.config.replot_lowres_map_faceon = True
            self.config.replot_lowres_map_edgeon = True

        # Set replot highres maps flags
        if self.config.replot_highres_maps:
            self.config.replot_highres_map_midplane = True
            self.config.replot_highres_map_faceon = True
            self.config.replot_highres_map_edgeon = True

        # Lowres midplane
        if self.config.replot_lowres_map_midplane:
            if self.has_lowres_map_midplane_plot: fs.remove_file(self.lowres_map_midplane_plot_path)

        # Lowres faceon
        if self.config.replot_lowres_map_faceon:
            if self.has_lowres_map_faceon_plot: fs.remove_file(self.lowres_map_faceon_plot_path)

        # Lowres edgeon
        if self.config.replot_lowres_map_edgeon:
            if self.has_lowres_map_edgeon_plot: fs.remove_file(self.lowres_map_edgeon_plot_path)

        # Highres midplane
        if self.config.replot_highres_map_midplane:
            if self.has_highres_map_midplane_plot: fs.remove_file(self.highres_map_midplane_plot_path)

        # Highres faceon
        if self.config.replot_highres_map_faceon:
            if self.has_highres_map_faceon_plot: fs.remove_file(self.highres_map_faceon_plot_path)

        # Highres edgeon
        if self.config.replot_highres_map_edgeon:
            if self.has_highres_map_edgeon_plot: fs.remove_file(self.highres_map_edgeon_plot_path)

    # -----------------------------------------------------------------
    # PROJECTIONS
    # -----------------------------------------------------------------

    @property
    def earth_projection(self):
        return self.analysis_run.earth_projection

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):
        return self.analysis_run.faceon_projection

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):
        return self.analysis_run.edgeon_projection

    # -----------------------------------------------------------------
    # UNITS
    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

    # -----------------------------------------------------------------

    @property
    def luminosity_unit(self):
        return "W"

    # -----------------------------------------------------------------
    # TOTAL
    # -----------------------------------------------------------------

    @property
    def total_absorptions_name(self):
        return "Labs_total"

    # -----------------------------------------------------------------

    @property
    def total_absorptions_description(self):
        return "Cell absorbed luminosities in the total simulation"

    # -----------------------------------------------------------------

    @property
    def total_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_total.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_absorptions(self):
        return fs.is_file(self.total_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_absorptions_path", True, write=False)
    def total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.total_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.total_contribution_absorption_luminosities * conversion_factor

        # Create the data with external xyz
        return Data3D.from_values(self.total_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.total_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # OLD
    # -----------------------------------------------------------------

    @property
    def old_absorptions_name(self):
        return "Labs_old"

    # -----------------------------------------------------------------

    @property
    def old_absorptions_description(self):
        return "Cell absorbed luminosities of the old stars"

    # -----------------------------------------------------------------

    @property
    def old_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_old.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_absorptions(self):
        return fs.is_file(self.old_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "old_absorptions_path", True, write=False)
    def old_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.old_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.old_contribution_absorption_luminosities * conversion_factor

        # Create the data with external xyz
        return Data3D.from_values(self.old_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.old_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # YOUNG
    # -----------------------------------------------------------------

    @property
    def young_absorptions_name(self):
        return "Labs_young"

    # -----------------------------------------------------------------

    @property
    def young_absorptions_description(self):
        return "Cell absorbed luminosities of the young stars"

    # -----------------------------------------------------------------

    @property
    def young_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_young.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_absorptions(self):
        return fs.is_file(self.young_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "young_absorptions_path", True, write=False)
    def young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.young_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.young_contribution_absorption_luminosities * conversion_factor

        # Create the data with external xyz
        return Data3D.from_values(self.young_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.young_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # IONIZING
    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_name(self):
        return "Labs_ionizing"

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_description(self):
        return "Cell absorbed luminosities of the ionizing stars"

    # -----------------------------------------------------------------

    @property
    def ionizing_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_ionizing.dat")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_absorptions(self):
        return fs.is_file(self.ionizing_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "ionizing_absorptions_path", True, write=False)
    def ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.ionizing_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.ionizing_contribution_absorption_luminosities * conversion_factor

        # Create the data with external xyz
        return Data3D.from_values(self.ionizing_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.ionizing_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # INTERNAL
    # -----------------------------------------------------------------

    @property
    def internal_absorptions_name(self):
        return "Labs_internal"

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_description(self):
        return "Cell absorbed luminosities in star formation regions dust"

    # -----------------------------------------------------------------

    @property
    def internal_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_internal.dat")

    # -----------------------------------------------------------------

    @property
    def has_internal_absorptions(self):
        return fs.is_file(self.internal_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def volumes(self):
        return self.model.cell_volumes

    # -----------------------------------------------------------------

    @property
    def ionizing_densities(self):
        return self.model.sfr_cell_stellar_density

    # -----------------------------------------------------------------

    @property
    def internal_absorbed_luminosity(self):
        return self.model.intrinsic_dust_luminosity_sfr.to(self.luminosity_unit, distance=self.galaxy_distance).value

    # -----------------------------------------------------------------

    @lazyproperty
    def relative_masses_ionizing(self):
        masses = self.volumes * self.ionizing_densities
        masses /= sum(masses) # normalized
        return masses

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "internal_absorptions_path", True, write=False)
    def internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # absorptions_ionizing_internal_watt = volumes * density * absorbed_energy
        luminosities_watt = self.relative_masses_ionizing * self.internal_absorbed_luminosity

        # Create the data with external xyz
        return Data3D.from_values(self.internal_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.internal_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # EXTRA
    # -----------------------------------------------------------------

    @property
    def extra_absorptions_name(self):
        return "Labs_extra"

    # -----------------------------------------------------------------

    @property
    def extra_absorptions_description(self):
        return "Cell absorbed luminosities of the extra component"

    # -----------------------------------------------------------------

    @property
    def extra_absorptions_path(self):
        return fs.join(self.cell_heating_path, "absorptions_extra.dat")

    # -----------------------------------------------------------------

    @property
    def has_extra_absorptions(self):
        return fs.is_file(self.extra_absorptions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "extra_absorptions_path", True, write=False)
    def extra_absorptions(self):

        """
        This function ...
        :return:
        """

        # Calculate luminosities
        conversion_factor = self.extra_contribution_absorption_unit.conversion_factor(self.luminosity_unit)
        luminosities_watt = self.extra_contribution_absorption_luminosities * conversion_factor

        # Create the data with external xyz
        return Data3D.from_values(self.extra_absorptions_name, luminosities_watt,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.extra_absorptions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, unit=self.luminosity_unit,
                                  distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # DIFFUSE HEATING FRACTIONS
    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_name(self):
        return "Funev_diffuse"

    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_description(self):
        return "Heating fraction by the unevolved stellar populations (diffuse dust)"

    # -----------------------------------------------------------------

    @property
    def young_absorption_values(self):
        return self.young_absorptions.values

    # -----------------------------------------------------------------

    @property
    def ionizing_absorption_values(self):
        return self.ionizing_absorptions.values

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_values_diffuse(self):
        return self.young_absorption_values + self.ionizing_absorption_values

    # -----------------------------------------------------------------

    @property
    def extra_absorption_values(self):
        return self.extra_absorptions.values

    # -----------------------------------------------------------------

    @property
    def total_absorption_values_diffuse(self):
        return self.total_absorptions.values

    # -----------------------------------------------------------------

    @property
    def diffuse_fractions_path(self):
        return fs.join(self.cell_heating_path, "diffuse_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_diffuse_fractions(self):
        return fs.is_file(self.diffuse_fractions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "diffuse_fractions_path", True, write=False)
    def diffuse_fractions(self):

        """
        This function ...
        :return:
        """

        # Calculate
        fractions = self.unevolved_absorption_values_diffuse / self.total_absorption_values_diffuse

        # Create the data with external xyz
        return Data3D.from_values(self.diffuse_fractions_name, fractions,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.diffuse_fractions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def diffuse_fraction_values(self):
        return self.diffuse_fractions.values

    # -----------------------------------------------------------------
    # TOTAL HEATING FRACTIONS
    # -----------------------------------------------------------------

    @property
    def total_fractions_name(self):
        return "Funev_total"

    # -----------------------------------------------------------------

    @property
    def total_fractions_description(self):
        return "Heating fraction by the unevolved stellar population (all dust)"

    # -----------------------------------------------------------------

    @property
    def internal_absorption_values(self):
        return self.internal_absorptions.values

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_values(self):
        return self.unevolved_absorption_values_diffuse + self.internal_absorption_values

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_values(self):
        return self.total_absorption_values_diffuse + self.internal_absorption_values

    # -----------------------------------------------------------------

    @property
    def total_fractions_path(self):
        return fs.join(self.cell_heating_path, "total_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_fractions(self):
        return fs.is_file(self.total_fractions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "total_fractions_path", True, write=False)
    def total_fractions(self):

        """
        This function ...
        :return:
        """

        # Calculate
        fractions = self.unevolved_absorption_values / self.total_absorption_values

        # Create the data with external xyz
        return Data3D.from_values(self.total_fractions_name, fractions,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.total_fractions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def total_fraction_values(self):
        return self.total_fractions.values

    # -----------------------------------------------------------------
    # EXTRA HEATING FRACTIONS
    # -----------------------------------------------------------------

    @property
    def extra_fractions_name(self):
        return "Fextra_total"

    # -----------------------------------------------------------------

    @property
    def extra_fractions_description(self):
        return "Heating fraction by the extra component (all dust)"

    # -----------------------------------------------------------------

    @property
    def extra_fractions_path(self):
        return fs.join(self.cell_heating_path, "extra_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_extra_fractions(self):
        return fs.is_file(self.extra_fractions_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Data3D, "extra_fractions_path", True, write=False)
    def extra_fractions(self):

        """
        This function ...
        :return:
        """

        # Calculate
        fractions = self.extra_absorption_values / self.total_absorption_values

        # Create the data with external xyz
        return Data3D.from_values(self.extra_fractions_name, fractions,
                                  self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                  length_unit=self.length_unit, description=self.extra_fractions_description,
                                  xyz_filepath=self.cell_coordinates_filepath, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def extra_fraction_values(self):
        return self.extra_fractions.values

    # -----------------------------------------------------------------
    # VALID CELLS
    # -----------------------------------------------------------------

    @lazyproperty
    def zero_absorption_mask(self):
        return self.total_absorptions.zeroes

    # -----------------------------------------------------------------

    @lazyproperty
    def diffuse_fractions_unphysical_mask(self):
        return self.diffuse_fractions.where_greater_than(1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def invalid_diffuse_fractions_mask(self):
        return self.diffuse_fractions.invalid + self.diffuse_fractions_unphysical_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_diffuse_fractions_mask(self):
        return np.logical_not(self.invalid_diffuse_fractions_mask)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_fractions_unphysical_mask(self):
        return self.total_fractions.where_greater_than(1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def invalid_total_fractions_mask(self):
        return self.total_fractions.invalid + self.total_fractions_unphysical_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_total_fractions_mask(self):
        return np.logical_not(self.invalid_total_fractions_mask)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_diffuse_fraction_values(self):
        return self.diffuse_fraction_values[self.valid_diffuse_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_total_fraction_values(self):
        return self.total_fraction_values[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_fractions_unphysical_mask(self):
        return self.extra_fractions.where_greater_than(1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def invalid_extra_fractions_mask(self):
        return self.extra_fractions.invalid + self.extra_fractions_unphysical_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_fractions_mask(self):
        return np.logical_not(self.invalid_extra_fractions_mask)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_fraction_values(self):
        return self.extra_fraction_values[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------

    @property
    def cell_weights(self):
        return self.cell_mass_fractions

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_diffuse_cell_weights(self):
        return self.cell_weights[self.valid_diffuse_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_total_cell_weights(self):
        return self.cell_weights[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_cell_weights(self):
        return self.cell_weights[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_x_coordinates(self):
        return self.total_fractions.x[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_y_coordinates(self):
        return self.total_fractions.y[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_z_coordinates(self):
        return self.total_fractions.z[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_radii(self):
        return self.total_fractions.radii[self.valid_total_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_x_coordinates(self):
        return self.extra_fractions.x[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_y_coordinates(self):
        return self.extra_fractions.y[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_z_coordinates(self):
        return self.extra_fractions.z[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_extra_radii(self):
        return self.extra_fractions.radii[self.valid_extra_fractions_mask]

    # -----------------------------------------------------------------
    # DISTRIBUTION DIFFUSE
    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_name(self):
        return "Heating fraction (diffuse dust)"

    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_path(self):
        return fs.join(self.cell_heating_path, "distribution_diffuse.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_diffuse(self):
        return fs.is_file(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution, "distribution_diffuse_path", True, write=False)
    def distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population (diffuse) ...")

        # Generate the distribution
        # Weights are dust mass fraction
        #print(self.valid_diffuse_fraction_values.shape)
        #print(self.valid_diffuse_cell_weights.shape)
        return Distribution.from_values(self.distribution_diffuse_name, self.valid_diffuse_fraction_values, nbins=self.config.nbins, weights=self.valid_diffuse_cell_weights)

    # -----------------------------------------------------------------
    # DISTRIBUTION TOTAL
    # -----------------------------------------------------------------

    @property
    def distribution_total_name(self):
        return "Heating fraction (all dust)"

    # -----------------------------------------------------------------

    @property
    def distribution_total_path(self):
        return fs.join(self.cell_heating_path, "distribution_total.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_total(self):
        return fs.is_file(self.distribution_total_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution, "distribution_total_path", True, write=False)
    def distribution_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population ...")

        # Generate the distribution
        # Weights are dust mass fraction
        return Distribution.from_values(self.distribution_total_name, self.valid_total_fraction_values, nbins=self.config.nbins, weights=self.valid_total_cell_weights)

    # -----------------------------------------------------------------
    # DISTRIBUTION EXTRA
    # -----------------------------------------------------------------

    @property
    def distribution_extra_name(self):
        return "Extra heating fraction (all dust)"

    # -----------------------------------------------------------------

    @property
    def distribution_extra_path(self):
        return fs.join(self.cell_heating_path, "distribution_extra.dat")

    # -----------------------------------------------------------------

    @property
    def has_distribution_extra(self):
        return fs.is_file(self.distribution_extra_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution, "distribution_extra_path", True, write=False)
    def distribution_extra(self):
        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the extra component ...")

        # Generate the distribution
        # Weights are dust mass fraction
        return Distribution.from_values(self.distribution_extra_name, self.valid_extra_fraction_values,
                                        nbins=self.config.nbins, weights=self.valid_extra_cell_weights)


    # -----------------------------------------------------------------
    # RADIAL DISTRIBUTION
    # -----------------------------------------------------------------

    @property
    def radial_distribution_name(self):
        #return "Heating fraction"
        return "Number of cells"

    # -----------------------------------------------------------------

    @property
    def radial_distribution_description(self):
        return "Radial distribution of the dust cell heating fraction"

    # -----------------------------------------------------------------

    @property
    def radial_distribution_path(self):
        return fs.join(self.cell_heating_path, "radial_distribution.fits")

    # -----------------------------------------------------------------

    @property
    def has_radial_distribution(self):
        return fs.is_file(self.radial_distribution_path)

    # -----------------------------------------------------------------

    @property
    def radial_distribution_x_name(self):
        return "Radius [" + str(self.length_unit) + "]"

    # -----------------------------------------------------------------

    @property
    def radial_distribution_y_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution2D, "radial_distribution_path", True, write=False)
    def radial_distribution(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the unevolved stellar population ...")

        # Generate the radial distribution
        # name, x, y, weights=None, nbins=200, x_name=None, y_name=None, x_unit=None, y_unit=None, description=None
        return Distribution2D.from_values(self.radial_distribution_name, self.valid_radii, self.valid_total_fraction_values,
                                          weights=self.valid_total_cell_weights, x_name=self.radial_distribution_x_name,
                                          y_name=self.radial_distribution_y_name, x_unit=self.length_unit,
                                          nbins=self.config.nradial_bins, description=self.radial_distribution_description)

    # -----------------------------------------------------------------
    # EXTRA RADIAL DISTRIBUTION
    # -----------------------------------------------------------------

    # -----------------------------------------------------------------

    @property
    def extra_radial_distribution_description(self):
        return "Radial distribution of the dust cell extra heating fraction"

    # -----------------------------------------------------------------

    @property
    def extra_radial_distribution_path(self):
        return fs.join(self.cell_heating_path, "extra_radial_distribution.fits")

    # -----------------------------------------------------------------

    @property
    def has_extra_radial_distribution(self):
        return fs.is_file(self.extra_radial_distribution_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Distribution2D, "extra_radial_distribution_path", True, write=False)
    def extra_radial_distribution(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the extra component ...")

        # Generate the radial distribution
        return Distribution2D.from_values(self.radial_distribution_name, self.valid_extra_radii,
                                          self.valid_extra_fraction_values,
                                          weights=self.valid_extra_cell_weights,
                                          x_name=self.radial_distribution_x_name,
                                          y_name=self.radial_distribution_y_name, x_unit=self.length_unit,
                                          nbins=self.config.nradial_bins,
                                          description=self.extra_radial_distribution_description)

    # -----------------------------------------------------------------
    # MIDPLANE HEIGHT
    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_component_scaleheight(self):

        """
        This function ...
        :return:
        """

        # Old stellar disk
        if self.config.midplane_component == old: return self.model.old_disk_scaleheight

        # Young stellar disk
        elif self.config.midplane_component == young: return self.model.young_scaleheight

        # Ionizing stellar disk
        elif self.config.midplane_component == ionizing: return self.model.sfr_scaleheight

        # Extra component
        elif self.config.midplane_component == extra: return self.model.extra_scaleheight

        # Dust disk
        elif self.config.midplane_component == dust: return self.model.dust_scaleheight

        # Invalid
        else: raise ValueError("Invalid midplane component: '" + self.config.midplane_component + "'")

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_height(self):
        return self.midplane_component_scaleheight * self.config.midplane_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def midplane_height_pc(self):
        return self.midplane_height.to("pc").value

    # -----------------------------------------------------------------
    # LOWRES MAPS
    # -----------------------------------------------------------------
    #   MIDPLANE
    # -----------------------------------------------------------------

    @property
    def lowres_map_midplane_path(self):
        return fs.join(self.cell_heating_path, "lowres_midplane.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_midplane(self):
        return fs.is_file(self.lowres_map_midplane_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "lowres_map_midplane_path", True, write=False)
    def lowres_map_midplane(self):

        """
        Thisf unction ...
        :return:
        """

        return project_faceon(self.total_fractions_name, self.total_fractions, spacing=self.config.lowres_map_spacing,
                              spacing_factor=self.config.lowres_map_spacing_factor, height=self.midplane_height,
                              return_stddev=True, return_ncells=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def lowres_frame_midplane(self):
        return self.lowres_map_midplane.primary

    # -----------------------------------------------------------------

    @property
    def lowres_stddev_midplane(self):
        return self.lowres_map_midplane.frames.stddev

    # -----------------------------------------------------------------

    @property
    def lowres_ncells_midplane(self):
        return self.lowres_map_midplane.frames.ncells

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def lowres_map_faceon_path(self):
        return fs.join(self.cell_heating_path, "lowres_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_faceon(self):
        return fs.is_file(self.lowres_map_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "lowres_map_faceon_path", True, write=False)
    def lowres_map_faceon(self):

        """
        This function ...
        :return:
        """

        return project_faceon(self.total_fractions_name, self.total_fractions, spacing=self.config.lowres_map_spacing,
                              spacing_factor=self.config.lowres_map_spacing_factor, return_stddev=True, return_ncells=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def lowres_frame_faceon(self):
        return self.lowres_map_faceon.primary

    # -----------------------------------------------------------------

    @property
    def lowres_stddev_faceon(self):
        return self.lowres_map_faceon.frames.stddev

    # -----------------------------------------------------------------

    @property
    def lowres_ncells_faceon(self):
        return self.lowres_map_faceon.frames.ncells

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def lowres_map_edgeon_path(self):
        return fs.join(self.cell_heating_path, "lowres_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_edgeon(self):
        return fs.is_file(self.lowres_map_edgeon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "lowres_map_edgeon_path", True, write=False)
    def lowres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return project_edgeon(self.total_fractions_name, self.total_fractions, spacing=self.config.lowres_map_spacing,
                              spacing_factor=self.config.lowres_map_spacing_factor, return_stddev=True, return_ncells=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def lowres_frame_edgeon(self):
        return self.lowres_map_edgeon.primary

    # -----------------------------------------------------------------

    @property
    def lowres_stddev_edgeon(self):
        return self.lowres_map_edgeon.frames.stddev

    # -----------------------------------------------------------------

    @property
    def lowres_ncells_edgeon(self):
        return self.lowres_map_edgeon.frames.ncells

    # -----------------------------------------------------------------
    # HIGHRES MAPS
    # -----------------------------------------------------------------
    #   MIDPLANE
    # -----------------------------------------------------------------

    @property
    def highres_map_midplane_path(self):
        return fs.join(self.cell_heating_path, "highres_midplane.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_midplane(self):
        return fs.is_file(self.highres_map_midplane_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "highres_map_midplane_path", True, write=False)
    def highres_map_midplane(self):

        """
        This function ...
        :return:
        """

        return project_data(self.total_fractions_name, self.total_fractions, self.faceon_projection,
                            return_stddev=True, return_ncells=True, as_image=True, height=self.midplane_height)

    # -----------------------------------------------------------------

    @property
    def highres_frame_midplane(self):
        return self.highres_map_midplane.primary

    # -----------------------------------------------------------------

    @property
    def highres_stddev_midplane(self):
        return self.highres_map_midplane.frames.stddev

    # -----------------------------------------------------------------

    @property
    def highres_ncells_midplane(self):
        return self.highres_map_midplane.frames.ncells

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def highres_map_faceon_path(self):
        return fs.join(self.cell_heating_path, "highres_faceon.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_faceon(self):
        return fs.is_file(self.highres_map_faceon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "highres_map_faceon_path", True, write=False)
    def highres_map_faceon(self):

        """
        This function ...
        :return:
        """

        return project_data(self.total_fractions_name, self.total_fractions, self.faceon_projection,
                            return_ncells=True, return_stddev=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def highres_frame_faceon(self):
        return self.highres_map_faceon.primary

    # -----------------------------------------------------------------

    @property
    def highres_stddev_faceon(self):
        return self.highres_map_faceon.frames.stddev

    # -----------------------------------------------------------------

    @property
    def highres_ncells_faceon(self):
        return self.highres_map_faceon.frames.ncells

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def highres_map_edgeon_path(self):
        return fs.join(self.cell_heating_path, "highres_edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_edgeon(self):
        return fs.is_file(self.highres_map_edgeon_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Image, "highres_map_edgeon_path", True, write=False)
    def highres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        return project_data(self.total_fractions_name, self.total_fractions, self.edgeon_projection,
                            return_ncells=True, return_stddev=True, as_image=True)

    # -----------------------------------------------------------------

    @property
    def highres_frame_edgeon(self):
        return self.highres_map_edgeon.primary

    # -----------------------------------------------------------------

    @property
    def highres_stddev_edgeon(self):
        return self.highres_map_edgeon.frames.stddev

    # -----------------------------------------------------------------

    @property
    def highres_ncells_edgeon(self):
        return self.highres_map_edgeon.frames.ncells

    # -----------------------------------------------------------------
    # INTERPOLATED MAPS
    # -----------------------------------------------------------------

    def interpolate_map(self, frame, ncells, replace_nans=True):

        """
        This function ...
        :param frame:
        :param ncells:
        :param replace_nans:
        :return:
        """

        # Inform the user
        #log.info("Creating the interpolated version of the map of the heating fraction ...")

        # Copy
        interpolated = frame.copy()

        # Get outside nans
        outside_nans = interpolated.nans.largest()
        # plotting.plot_mask(outside_nans, title="outside nans")
        #outside_nans.saveto(fs.join(self.cell_heating_path, "outside_nans.fits"))
        not_nans = outside_nans.inverse()
        not_nans.disk_dilate(radius=self.config.not_nans_dilation_radius)
        # not_nans.fill_holes()
        #not_nans.saveto(fs.join(self.cell_heating_path, "not_nans.fits"))
        do_nans = not_nans.largest().inverse()

        # Get mask
        where = ncells.where_smaller_than(self.config.min_ncells)
        where = where * do_nans.inverse()  # don't interpolate outside (where ncells = 0)

        # plotting.plot_mask(where, title="where smaller than " + str(self.config.min_ncells))
        # plotting.plot_mask(self.map_interpolated.nans, title="nans")

        # Replace NaNs to zero that have to stay NaNs (don't interpolate)
        if replace_nans:
            interpolated[do_nans] = 0.0
            do_nans.disk_dilate(radius=self.config.not_nans_dilation_radius)

        # Put pixels to NaN
        interpolated.replace_by_nans(where)

        # plotting.plot_mask(self.map_interpolated.nans, title="nans")
        # exit()

        # Interpolate nans
        interpolated.interpolate_nans(sigma=2., error_on_max=replace_nans)
        interpolated.replace_by_nans(do_nans)

        # Return the interpolated frame
        return interpolated

    # -----------------------------------------------------------------
    # LOWRES MIDPLANE
    # -----------------------------------------------------------------

    @property
    def lowres_frame_midplane_interpolated_path(self):
        return fs.join(self.cell_heating_path, "lowres_midplane_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_frame_midplane_interpolated(self):
        return fs.is_file(self.lowres_frame_midplane_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "lowres_frame_midplane_interpolated_path", True)
    def lowres_frame_midplane_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolate
        return self.interpolate_map(self.lowres_frame_midplane, self.lowres_ncells_midplane)

    # -----------------------------------------------------------------
    # LOWRES FACEON
    # -----------------------------------------------------------------

    @property
    def lowres_frame_faceon_interpolated_path(self):
        return fs.join(self.cell_heating_path, "lowres_faceon_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_frame_faceon_interpolated(self):
        return fs.is_file(self.lowres_frame_faceon_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "lowres_frame_faceon_interpolated_path", True)
    def lowres_frame_faceon_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolate
        return self.interpolate_map(self.lowres_frame_faceon, self.lowres_ncells_faceon)

    # -----------------------------------------------------------------
    # LOWRES EDGEON
    # -----------------------------------------------------------------

    @property
    def lowres_frame_edgeon_interpolated_path(self):
        return fs.join(self.cell_heating_path, "lowres_edgeon_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_lowres_frame_edgeon_interpolated(self):
        return fs.is_file(self.lowres_frame_edgeon_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "lowres_frame_edgeon_interpolated_path", True)
    def lowres_frame_edgeon_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolate
        return self.interpolate_map(self.lowres_frame_edgeon, self.lowres_ncells_edgeon)

    # -----------------------------------------------------------------
    # HIGHRES MIDPLANE
    # -----------------------------------------------------------------

    @property
    def highres_frame_midplane_interpolated_path(self):
        return fs.join(self.cell_heating_path, "highres_midplane_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_frame_midplane_interpolated(self):
        return fs.is_file(self.highres_frame_midplane_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "highres_frame_midplane_interpolated_path", True)
    def highres_frame_midplane_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolate
        return self.interpolate_map(self.highres_frame_midplane, self.highres_ncells_midplane)

    # -----------------------------------------------------------------
    # HIGHRES FACEON
    # -----------------------------------------------------------------

    @property
    def highres_frame_faceon_interpolated_path(self):
        return fs.join(self.cell_heating_path, "highres_faceon_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_frame_faceon_interpolated(self):
        return fs.is_file(self.highres_frame_faceon_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "highres_frame_faceon_interpolated_path", True)
    def highres_frame_faceon_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolated
        return self.interpolate_map(self.highres_frame_faceon, self.highres_ncells_faceon)

    # -----------------------------------------------------------------
    # HIGHRES EDGEON
    # -----------------------------------------------------------------

    @property
    def highres_frame_edgeon_interpolated_path(self):
        return fs.join(self.cell_heating_path, "highres_edgeon_interpolated.fits")

    # -----------------------------------------------------------------

    @property
    def has_highres_frame_edgeon_interpolated(self):
        return fs.is_file(self.highres_frame_edgeon_interpolated_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Frame, "highres_frame_edgeon_interpolated_path", True)
    def highres_frame_edgeon_interpolated(self):

        """
        This function ...
        :return:
        """

        # Interpolated
        return self.interpolate_map(self.highres_frame_edgeon, self.highres_ncells_edgeon)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write absorptions
        self.write_absorptions()

        # Write the heating fractions
        self.write_fractions()

        # Write the distributions
        self.write_distributions()

        # Write the radial distributions
        self.write_radial_distributions()

        # Write maps
        self.write_maps()

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption data ...")

        # Total
        if self.do_write_total_absorptions: self.write_total_absorptions()

        # Old
        if self.do_write_old_absorptions: self.write_old_absorptions()

        # Young
        if self.do_write_young_absorptions: self.write_young_absorptions()

        # Ionizing
        if self.do_write_ionizing_absorptions: self.write_ionizing_absorptions()

        # Internal
        if self.do_write_internal_absorptions: self.write_internal_absorptions()

        # Extra
        if self.do_write_extra_absorptions: self.write_extra_absorptions()

    # -----------------------------------------------------------------

    @property
    def do_write_total_absorptions(self):
        return not self.has_total_absorptions

    # -----------------------------------------------------------------

    def write_total_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.total_absorptions.saveto(self.total_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_old_absorptions(self):
        return not self.has_old_absorptions

    # -----------------------------------------------------------------

    def write_old_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.old_absorptions.saveto(self.old_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_young_absorptions(self):
        return not self.has_young_absorptions

    # -----------------------------------------------------------------

    def write_young_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.young_absorptions.saveto(self.young_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ionizing_absorptions(self):
        return not self.has_ionizing_absorptions

    # -----------------------------------------------------------------

    def write_ionizing_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.ionizing_absorptions.saveto(self.ionizing_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_internal_absorptions(self):
        return not self.has_internal_absorptions

    # -----------------------------------------------------------------

    def write_internal_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.internal_absorptions.saveto(self.internal_absorptions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_extra_absorptions(self):
        return self.do_extra and not self.has_extra_absorptions

    # -----------------------------------------------------------------

    def write_extra_absorptions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.extra_absorptions.saveto(self.extra_absorptions_path)

    # -----------------------------------------------------------------

    def write_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the heating fraction data ...")

        # Diffuse
        if self.do_write_diffuse_fractions: self.write_diffuse_fractions()

        # Total
        if self.do_write_total_fractions: self.write_total_fractions()

        # Extra
        if self.do_write_extra_fractions: self.write_extra_fractions()

    # -----------------------------------------------------------------

    @property
    def do_write_diffuse_fractions(self):
        return not self.has_diffuse_fractions

    # -----------------------------------------------------------------

    def write_diffuse_fractions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.diffuse_fractions.saveto(self.diffuse_fractions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_total_fractions(self):
        return not self.has_total_fractions

    # -----------------------------------------------------------------

    def write_total_fractions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.total_fractions.saveto(self.total_fractions_path)

    # -----------------------------------------------------------------

    @property
    def do_write_extra_fractions(self):
        return not self.has_extra_fractions

    # -----------------------------------------------------------------

    def write_extra_fractions(self):

        """
        This function ...
        :return:
        """

        # Save
        self.extra_fractions.saveto(self.extra_fractions_path)

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distributions of the heating fraction ...")

        # Diffuse
        if self.do_write_distribution_diffuse: self.write_distribution_diffuse()

        # Total
        if self.do_write_distribution_total: self.write_distribution_total()

        # Extra
        if self.do_write_distribution_extra: self.write_distribution_extra()

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_diffuse(self):
        return not self.has_distribution_diffuse

    # -----------------------------------------------------------------

    def write_distribution_diffuse(self):

        """
        This function ...
        :return:
        """

        # Save
        self.distribution_diffuse.saveto(self.distribution_diffuse_path)

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_total(self):
        return not self.has_distribution_total

    # -----------------------------------------------------------------

    def write_distribution_total(self):

        """
        This function ...
        :return:
        """

        # Save
        self.distribution_total.saveto(self.distribution_total_path)

    # -----------------------------------------------------------------

    @property
    def do_write_distribution_extra(self):
        return self.do_extra and not self.has_distribution_extra

    # -----------------------------------------------------------------

    def write_distribution_extra(self):

        """
        This function ...
        :return:
        """

        # Save
        self.distribution_extra.saveto(self.distribution_extra_path)

    # -----------------------------------------------------------------

    @property
    def do_write_radial_distribution(self):
        return not self.has_radial_distribution

    # -----------------------------------------------------------------

    def write_radial_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the radial distributions of the heating fraction ...")

        # Radial
        if self.do_write_radial_distribution: self.write_radial_distribution()

        # Extra radial
        if self.do_write_extra_radial_distribution: self.write_extra_radial_distribution()

    # -----------------------------------------------------------------

    def write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Save
        self.radial_distribution.saveto(self.radial_distribution_path)

    # -----------------------------------------------------------------

    @property
    def do_write_extra_radial_distribution(self):
        return self.do_extra and not self.has_extra_radial_distribution

    # -----------------------------------------------------------------

    def write_extra_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Save
        self.extra_radial_distribution.saveto(self.extra_radial_distribution_path)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Low-res
        if self.config.lowres: self.write_lowres_maps()

        # High-res
        if self.config.highres: self.write_highres_maps()

    # -----------------------------------------------------------------

    def write_lowres_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the low-resolution heating fraction maps ...")

        # Midplane
        if self.do_write_lowres_map_midplane: self.write_lowres_map_midplane()

        # Face-on
        if self.do_write_lowres_map_faceon: self.write_lowres_map_faceon()

        # Edge-on
        if self.do_write_lowres_map_edgeon: self.write_lowres_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_lowres_map_midplane(self):
        return not self.has_lowres_map_midplane

    # -----------------------------------------------------------------

    def write_lowres_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Save
        self.lowres_map_midplane.saveto(self.lowres_map_midplane_path)

    # -----------------------------------------------------------------

    @property
    def do_write_lowres_map_faceon(self):
        return not self.has_lowres_map_faceon

    # -----------------------------------------------------------------

    def write_lowres_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.lowres_map_faceon.saveto(self.lowres_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_lowres_map_edgeon(self):
        return not self.has_lowres_map_edgeon

    # -----------------------------------------------------------------

    def write_lowres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.lowres_map_edgeon.saveto(self.lowres_map_edgeon_path)

    # -----------------------------------------------------------------

    def write_highres_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the high-resolution heating fraction maps ...")

        # Midplane
        if self.do_write_highres_map_midplane: self.write_highres_map_midplane()

        # Face-on
        if self.do_write_highres_map_faceon: self.write_highres_map_faceon()

        # Edge-on
        if self.do_write_highres_map_edgeon: self.write_highres_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def do_write_highres_map_midplane(self):
        return not self.has_highres_map_midplane

    # -----------------------------------------------------------------

    def write_highres_map_midplane(self):

        """
        This function ...
        :return:
        """

        # Save
        self.highres_map_midplane.saveto(self.highres_map_midplane_path)

    # -----------------------------------------------------------------

    @property
    def do_write_highres_map_faceon(self):
        return not self.has_highres_map_faceon

    # -----------------------------------------------------------------

    def write_highres_map_faceon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.highres_map_faceon.saveto(self.highres_map_faceon_path)

    # -----------------------------------------------------------------

    @property
    def do_write_highres_map_edgeon(self):
        return not self.has_highres_map_edgeon

    # -----------------------------------------------------------------

    def write_highres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.highres_map_edgeon.saveto(self.highres_map_edgeon_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Distributions
        self.plot_distributions()

        # Radial distributions
        self.plot_radial_distributions()

        # Maps
        self.plot_maps()

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the distributions ...")

        # Diffuse
        if self.do_plot_distribution_diffuse: self.plot_distribution_diffuse()

        # Total
        if self.do_plot_distribution_total: self.plot_distribution_total()

        # Extra
        if self.do_plot_distribution_extra: self.plot_distribution_extra()

        # Diffuse and total
        if self.do_plot_distribution_diffuse_and_total: self.plot_distribution_diffuse_and_total()

    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_plot_path(self):
        return fs.join(self.cell_heating_path, "distribution_diffuse.pdf")

    # -----------------------------------------------------------------

    @property
    def has_distribution_diffuse_plot(self):
        return fs.is_file(self.distribution_diffuse_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_distribution_diffuse(self):
        return not self.has_distribution_diffuse_plot

    # -----------------------------------------------------------------

    def plot_distribution_diffuse(self):

        """
        Thisf unction ...
        :return:
        """

        # Plot
        plot_distribution(self.distribution_diffuse, path=self.distribution_diffuse_plot_path)

    # -----------------------------------------------------------------

    @property
    def distribution_total_plot_path(self):
        return fs.join(self.cell_heating_path, "distribution_total.pdf")

    # -----------------------------------------------------------------

    @property
    def has_distribution_total_plot(self):
        return fs.is_file(self.distribution_total_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_distribution_total(self):
        return not self.has_distribution_total_plot

    # -----------------------------------------------------------------

    def plot_distribution_total(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_distribution(self.distribution_total, path=self.distribution_total_plot_path)

    # -----------------------------------------------------------------


    @property
    def distribution_extra_plot_path(self):
        return fs.join(self.cell_heating_path, "distribution_extra.pdf")

    # -----------------------------------------------------------------

    @property
    def has_distribution_extra_plot(self):
        return fs.is_file(self.distribution_extra_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_distribution_extra(self):
        return self.do_extra and not self.has_distribution_extra_plot

    # -----------------------------------------------------------------

    def plot_distribution_extra(self):

        """
        This function ...
        :return:
        """

        # Plot
        plot_distribution(self.distribution_extra, path=self.distribution_extra_plot_path)

    # -----------------------------------------------------------------

    @property
    def distribution_diffuse_and_total_plot_path(self):
        return fs.join(self.cell_heating_path, "distribution_diffuse_and_total.pdf")

    # -----------------------------------------------------------------

    @property
    def has_distribution_diffuse_and_total_plot(self):
        return fs.is_file(self.distribution_diffuse_and_total_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_distribution_diffuse_and_total(self):
        return not self.has_distribution_diffuse_and_total_plot

    # -----------------------------------------------------------------

    def plot_distribution_diffuse_and_total(self):

        """
        This function ...
        :return:
        """

        # Set distributions
        distributions = OrderedDict()
        distributions["unevolved"] = self.distribution_total
        distributions["unevolved (diffuse)"] = self.distribution_diffuse

        # Plot
        plot_distributions(distributions, path=self.distribution_diffuse_and_total_plot_path, alpha=0.5)

    # -----------------------------------------------------------------

    def plot_radial_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the radial distributions ...")

        # Radial
        if self.do_plot_radial_distribution: self.plot_radial_distribution()

        # Extra radial
        if self.do_plot_extra_radial_distribution: self.plot_extra_radial_distribution()

    # -----------------------------------------------------------------

    @property
    def radial_distribution_plot_path(self):
        return fs.join(self.cell_heating_path, "radial_distribution.pdf")

    # -----------------------------------------------------------------

    @property
    def has_radial_distribution_plot(self):
        return fs.is_file(self.radial_distribution_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_radial_distribution(self):
        return not self.has_radial_distribution_plot

    # -----------------------------------------------------------------

    @lazyproperty
    def radial_distribution_plot_radii_pc(self):

        """
        This function ...
        :return:
        """

        if self.config.radial_distribution_radii is None: return []
        else: return [radius.to("pc").value for radius in self.config.radial_distribution_radii]

    # -----------------------------------------------------------------

    def plot_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Set title
        title = "Radial distribution of the heating fraction of the unevolved stellar population"

        # Create the plot file
        plot_2d_distribution(self.radial_distribution, x_lines=self.radial_distribution_plot_radii_pc, title=title, path=self.radial_distribution_plot_path)

    # -----------------------------------------------------------------

    @property
    def extra_radial_distribution_plot_path(self):
        return fs.join(self.cell_heating_path, "extra_radial_distribution.pdf")

    # -----------------------------------------------------------------

    @property
    def has_extra_radial_distribution_plot(self):
        return fs.is_file(self.extra_radial_distribution_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_extra_radial_distribution(self):
        return not self.has_extra_radial_distribution_plot

    # -----------------------------------------------------------------

    def plot_extra_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Set title
        title = "Radial distribution of the heating fraction of the extra component"

        # Create the plot file
        plot_2d_distribution(self.extra_radial_distribution, x_lines=self.radial_distribution_plot_radii_pc, title=title, path=self.extra_radial_distribution_plot_path)

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the heating fraction ...")

        # Lowres
        if self.config.lowres: self.plot_lowres_maps()

        # Highres
        if self.config.highres: self.plot_highres_maps()

    # -----------------------------------------------------------------

    def plot_heating_fraction_map(self, frame, path=None):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        # Plot
        plot_map(frame, path=path, interval=(0., 1.), contours=self.config.contours, ncontours=self.config.nlevels)

    # -----------------------------------------------------------------

    def plot_lowres_maps(self):

        """
        This function ...
        :return:
        """

        # Midplane
        if self.do_plot_lowres_map_midplane: self.plot_lowres_map_midplane()

        # Faceon
        if self.do_plot_lowres_map_faceon: self.plot_lowres_map_faceon()

        # Edgeon
        if self.do_plot_lowres_map_edgeon: self.plot_lowres_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def lowres_map_midplane_plot_path(self):
        return fs.join(self.cell_heating_path, "lowres_midplane.pdf")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_midplane_plot(self):
        return fs.is_file(self.lowres_map_midplane_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_lowres_map_midplane(self):
        return not self.has_lowres_map_midplane_plot

    # -----------------------------------------------------------------

    def plot_lowres_map_midplane(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.lowres_frame_midplane_interpolated, path=self.lowres_map_midplane_plot_path)

    # -----------------------------------------------------------------

    @property
    def lowres_map_faceon_plot_path(self):
        return fs.join(self.cell_heating_path, "lowres_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_faceon_plot(self):
        return fs.is_file(self.lowres_map_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_lowres_map_faceon(self):
        return not self.has_lowres_map_faceon_plot

    # -----------------------------------------------------------------

    def plot_lowres_map_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.lowres_frame_faceon_interpolated, path=self.lowres_map_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def lowres_map_edgeon_plot_path(self):
        return fs.join(self.cell_heating_path, "lowres_edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_lowres_map_edgeon_plot(self):
        return fs.is_file(self.lowres_map_edgeon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_lowres_map_edgeon(self):
        return not self.has_lowres_map_edgeon_plot

    # -----------------------------------------------------------------

    def plot_lowres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.lowres_frame_edgeon_interpolated, path=self.lowres_map_edgeon_plot_path)

    # -----------------------------------------------------------------

    def plot_highres_maps(self):

        """
        This function ...
        :return:
        """

        # Midplane
        if self.do_plot_highres_map_midplane: self.plot_highres_map_midplane()

        # Faceon
        if self.do_plot_highres_map_faceon: self.plot_highres_map_faceon()

        # Edgeon
        if self.do_plot_highres_map_edgeon: self.plot_highres_map_edgeon()

    # -----------------------------------------------------------------

    @property
    def highres_map_midplane_plot_path(self):
        return fs.join(self.cell_heating_path, "highres_midplane.pdf")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_midplane_plot(self):
        return fs.is_file(self.highres_map_midplane_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_highres_map_midplane(self):
        return not self.has_highres_map_midplane_plot

    # -----------------------------------------------------------------

    def plot_highres_map_midplane(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.highres_frame_midplane_interpolated, path=self.highres_map_midplane_plot_path)

    # -----------------------------------------------------------------

    @property
    def highres_map_faceon_plot_path(self):
        return fs.join(self.cell_heating_path, "highres_faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_faceon_plot(self):
        return fs.is_file(self.highres_map_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_highres_map_faceon(self):
        return not self.has_highres_map_faceon_plot

    # -----------------------------------------------------------------

    def plot_highres_map_faceon(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.highres_frame_faceon_interpolated, path=self.highres_map_faceon_plot_path)

    # -----------------------------------------------------------------

    @property
    def highres_map_edgeon_plot_path(self):
        return fs.join(self.cell_heating_path, "highres_edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_highres_map_edgeon_plot(self):
        return fs.is_file(self.highres_map_edgeon_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_highres_map_edgeon(self):
        return not self.has_highres_map_edgeon_plot

    # -----------------------------------------------------------------

    def plot_highres_map_edgeon(self):

        """
        This function ...
        :return:
        """

        self.plot_heating_fraction_map(self.highres_frame_edgeon_interpolated, path=self.highres_map_edgeon_plot_path)

# -----------------------------------------------------------------
