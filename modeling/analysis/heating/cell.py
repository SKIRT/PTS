#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.celldustheating Contains the CellDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.basics.distribution import Distribution, Distribution2D
from .tables import AbsorptionTable
from ....core.tools.utils import lazyproperty

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

        # -- Attributes --

        # The table with the absorbed luminosities
        self.absorptions = None

        # The heating fraction of the unevolved stellar population for each dust cell
        self.heating_fractions = None

        # The distribution of heating fractions
        self.distribution = None

        # The 2D distribution of heating fractions
        self.radial_distribution = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the absorption table
        self.get_absorption_table()

        # 3. Calculate the heating fraction of the unevolved stellar population
        self.calculate_heating_unevolved()

        # 4. Calculate the distribution of the heating fraction of the unevolved stellar population
        self.calculate_distribution()

        # 5. Calculate the distribution of the heating fraction of the unevolved stellar population as a function of radius
        self.calculate_radial_distribution()

        # 6. Writing
        self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(CellDustHeatingAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def ncells(self):

        """
        This function ...
        :return:
        """

        return len(self.total_contribution_absorption_data)

    # -----------------------------------------------------------------

    def get_absorption_table(self):

        """
        This function ...
        :return:
        """

        # Create
        if not self.has_absorptions: self.create_absorption_table()

        # Load the table
        else: self.load_absorption_table()

    # -----------------------------------------------------------------

    def create_absorption_table(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating the absorption table ...")

        # Get the coordinates
        x = self.total_contribution_absorption_data["X coordinate of cell center"]
        y = self.total_contribution_absorption_data["Y coordinate of cell center"]
        z = self.total_contribution_absorption_data["Z coordinate of cell center"]

        # Get luminosity for total stellar population
        total_absorptions = self.total_contribution_absorption_data["Absorbed bolometric luminosity"]

        # Get luminosity for old stellar population
        old_absorptions = self.old_contribution_absorption_data["Absorbed bolometric luminosity"]

        # Get luminosity for young stellar population
        young_absorptions = self.young_contribution_absorption_data["Absorbed bolometric luminosity"]

        # Get luminosity for ionizing stellar population
        ionizing_absorptions = self.ionizing_contribution_absorption_data["Absorbed bolometric luminosity"]

        # Create the table
        self.absorptions = AbsorptionTable.from_columns(x, y, z, total_absorptions, old_absorptions, young_absorptions, ionizing_absorptions)

    # -----------------------------------------------------------------

    def load_absorption_table(self):

        """
        This function ...
        :return:
        """

        # Success
        log.success("Absorption table has already been created: loading from file ...")

        # Load
        self.absorptions = AbsorptionTable.from_file(self.absorption_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def zero_absorption_mask(self):

        """
        This function ...
        :return:
        """

        return self.absorptions.zero_absorption_mask

    # -----------------------------------------------------------------

    def calculate_heating_unevolved(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the heating fraction of the unevolved stellar population ...")

        # SEBA:

        ## Total energy absorbed in the new component
        ## Derived from ../../SKIRTrun/models/testHeating/MappingsHeating/plotSEDs.py
        #Lnew = 176495776.676  # in Lsun

        #print ''reading data... '
        #input = np.loadtxt(path + totISRFfile, skiprows=1)
        #ID = input[:, 0]
        #volume = input[:, 1]
        #density = input[:, 2]
        #massFrac = input[:, 3]
        #density_new = input[:, 5]
        #x = input[:, 6]
        #y = input[:, 7]
        #tot = input[:, 9] / Lsun
        #old = input[:, 10] / Lsun
        #yng = input[:, 11] / Lsun
        #new = input[:, 12] / Lsun

        #energy_new = volume * density_new * Lnew
        #F_abs_yng = (yng + new + energy_new) / (old + yng + new + energy_new)

        #cell_properties = self.model.cell_properties
        #volumes = cell_properties["Volume"]
        volumes = self.model.cell_volumes
        absorbed_energy = self.model.intrinsic_sfr_dust_luminosity
        density = self.model.sfr_cell_stellar_density

        #absorptions_unevolved_diffuse = self.absorptions["Absorbed bolometric luminosity of the young stellar population"] + self.absorptions["Absorbed bolometric luminosity of the ionizing stellar population"]
        #absorptions_unevolved_diffuse = self.absorptions["young"] + self.absorptions["ionizing"]
        absorptions_unevolved_diffuse = self.absorptions.unevolved(unit="W", asarray=True)

        #absorptions_ionizing_internal = None # TODO !!

        #absorptions_total = self.absorptions["Absorbed bolometric luminosity of the total stellar population"]
        #absorptions_total = self.absorptions["total"]
        absorptions_total = self.absorptions.total(unit="W", asarray=True)
        #absorptions_total = absorptions_unevolved_diffuse + absorptions_ionizing_internal + absorptions_evolved # TODO !!

        # Calculate the heating fraction of the unevolved stellar population in each dust cell
        self.heating_fractions = absorptions_unevolved_diffuse / absorptions_total

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_nans(self):

        """
        Thisn function ...
        :return:
        """

        return np.isnan(self.heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_infs(self):

        """
        This function ...
        :return:
        """

        return np.isinf(self.heating_fractions)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fraction_unphysical(self):

        """
        Thisn function ...
        :return:
        """

        greater_than_one_mask = self.heating_fractions > 1.0
        ngreater_than_one = np.sum(greater_than_one_mask)
        relative_ngreater_than_one = float(ngreater_than_one) / len(self.heating_fractions)

        # Debugging
        log.debug(str(ngreater_than_one) + " pixels have a heating fraction greater than unity (" + str(relative_ngreater_than_one * 100) + "%)")

        # Return
        return greater_than_one_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_fractions_mask(self):

        """
        This function ...
        :return:
        """

        return self.heating_fraction_nans + self.heating_fraction_infs + self.heating_fraction_unphysical

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_heating_fractions(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.heating_fractions, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_weights(self):

        """
        This function ...
        :return:
        """

        return self.cell_properties["Mass fraction"]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_weights(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.cell_weights, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def radii(self):

        """
        This function ...
        :return:
        """

        x_coords = self.absorptions["x"]
        y_coords = self.absorptions["y"]
        return np.sqrt(x_coords**2 + y_coords**2)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_radii(self):

        """
        This fuction ...
        :return:
        """

        return np.ma.MaskedArray(self.radii, self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def x_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return self.absorptions["x"]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_x_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.absorptions["x"], mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def y_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return self.absorptions["y"]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_y_coordinates(self):

        """
        Thisn function ...
        :return:
        """

        return np.ma.MaskedArray(self.absorptions["y"], mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    @lazyproperty
    def z_coordinates(self):

        """
        Thisnf unction ...
        :return:
        """

        return self.absorptions["z"]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_z_coordinates(self):

        """
        This function ...
        :return:
        """

        return np.ma.MaskedArray(self.absorptions["z"], mask=self.heating_fractions_mask).compressed()

    # -----------------------------------------------------------------

    def calculate_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the distribution of heating fractions of the unevolved stellar population ...")

        # Generate the distribution
        self.distribution = Distribution.from_values("Heating fraction", self.valid_heating_fractions, nbins=self.config.nbins, weights=self.valid_cell_weights)

    # -----------------------------------------------------------------

    def calculate_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the radial distribution of heating fractions of the unevolved stellar population ...")

        # Generate the radial distribution
        self.radial_distribution = Distribution2D.from_values(self.valid_radii, self.valid_heating_fractions,
                                                              weights=self.valid_cell_weights, x_name="radius (pc)",
                                                              y_name="Heating fraction of unevolved stars", nbins=self.config.nradial_bins)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the absorption table
        self.write_absorptions()

        # Write the distribution of heating fractions
        self.write_distribution()

        # Write the radial distribution of heating fractions
        self.write_radial_distribution()

    # -----------------------------------------------------------------

    @property
    def absorption_table_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.cell_heating_path, "absorptions.dat")

    # -----------------------------------------------------------------

    @property
    def has_absorptions(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.absorption_table_path)

    # -----------------------------------------------------------------

    def write_absorptions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption table ...")

        # Save the table
        self.absorptions.saveto(self.absorption_table_path)

    # -----------------------------------------------------------------

    @property
    def distribution_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the distribution file
        return fs.join(self.cell_heating_path, "distribution.dat")

    # -----------------------------------------------------------------

    def write_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distribution of heating fractions of the unevolved stellar population ...")

        # Save the distribution
        self.distribution.saveto(self.distribution_path)

    # -----------------------------------------------------------------

    @property
    def radial_distribution_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the radial distribution file
        return fs.join(self.cell_heating_path, "radial_distribution.dat")

    # -----------------------------------------------------------------

    def write_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the radial distribution of heating fractions of the unevolved stellar population ...")

        # Save the radial distribution
        self.radial_distribution.saveto(self.radial_distribution_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the distribution of heating fractions
        self.plot_distribution()

        # Plot the radial distribution of heating fractions
        self.plot_radial_distribution()

        # Plot a map of the heating fraction for a face-on view of the galaxy
        self.plot_map()

    # -----------------------------------------------------------------

    def plot_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a histogram of the distribution of heating fractions of the unevolved stellar population ...")

        # Determine the path to the plot file
        path = fs.join(self.cell_heating_path, "distribution.pdf")

        # Create the plot file
        self.distribution.plot(title="Distribution of the heating fraction of the unevolved stellar population", path=path)

    # -----------------------------------------------------------------

    def plot_radial_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a 2D histogram of the radial distribution of the heating fractions of the unevolved stellar population ...")

        # Determine the path to the plot file
        path = fs.join(self.cell_heating_path, "radial_distribution.pdf")

        # Create the plot file
        self.radial_distribution.plot(title="Radial distribution of the heating fraction of the unevolved stellar population", path=path)

    # -----------------------------------------------------------------

    @property
    def heating_map_plot_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the plot file
        return fs.join(self.cell_heating_path, "map.pdf")

    # -----------------------------------------------------------------

    def plot_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a map of the heating fraction of the unevolved stellar population for a face-on view of the galaxy ...")

        # Create figure
        plt.figure()

        x = self.valid_x_coordinates
        y = self.valid_y_coordinates
        z = self.valid_heating_fractions

        #plt.pcolormesh(x, y, z, cmap='RdBu', vmin=0.0, vmax=1.0)

        from matplotlib import mlab

        x_ticks = x
        y_ticks = y

        z_grid = mlab.griddata(x, y, z, x, y)

        from matplotlib.backends import backend_agg as agg
        from matplotlib import cm

        # plot
        # fig = Figure()  # create the figure
        fig = plt.figure()
        agg.FigureCanvasAgg(fig)  # attach the rasterizer
        ax = fig.add_subplot(1, 1, 1)  # make axes to plot on
        ax.set_title("Interpolated Contour Plot of Experimental Data")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        cmap = cm.get_cmap("hot")  # get the "hot" color map
        contourset = ax.contourf(x_ticks, y_ticks, z_grid, 10, cmap=cmap)

        # Plot
        plt.savefig(self.heating_map_plot_path)
        plt.close()

# -----------------------------------------------------------------
