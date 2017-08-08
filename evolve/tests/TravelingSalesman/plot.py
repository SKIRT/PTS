#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from PIL import Image, ImageDraw, ImageFont

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log

# -----------------------------------------------------------------

class Plotter(object):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Get properties
        self.coordinates = kwargs.pop("coordinates")
        self.frequency = kwargs.pop("frequency")

        # The output path
        self.output_path = None

        self.counter = 0

    # -----------------------------------------------------------------

    def add_generation(self, engine):

        """
        This function ...
        :param engine:
        :return:
        """

        if engine.currentGeneration % self.frequency == 0:

            best = engine.bestIndividual()
            self.add_best(best)

    # -----------------------------------------------------------------

    def add_best(self, best):

        """
        This function ...
        :param best:
        :return:
        """

        self.write_tour_to_img(best)
        self.counter += 1

    # -----------------------------------------------------------------

    def write_tour_to_img(self, tour):

        """
        The function to plot the graph
        :param tour:
        """

        padding = 20
        coords = [(x+padding,y+padding) for (x,y) in self.coordinates]
        maxx, maxy = 0, 0

        for x, y in coords:

          maxx = max(x, maxx)
          maxy = max(y, maxy)

        maxx += padding
        maxy += padding

        img = Image.new("RGB",(int(maxx),int(maxy)),color=(255,255,255))

        font = ImageFont.load_default()
        d = ImageDraw.Draw(img)
        num_cities = len(tour)

        # Loop over the cities
        for i in range(num_cities):

          j = (i+1) % num_cities
          city_i = tour[i]
          city_j = tour[j]
          x1,y1 = coords[city_i]
          x2,y2 = coords[city_j]

          d.line((int(x1),int(y1),int(x2),int(y2)),fill=(0,0,0))
          d.text((int(x1)+7,int(y1)-5),str(i),font=font,fill=(32,32,32))

        for x, y in coords:

          x, y = int(x),int(y)
          d.ellipse((x-5,y-5,x+5,y+5),outline=(0,0,0),fill=(196,196,196))

        del d

        # Save
        if self.output_path is not None:

            # Determine the plot path
            path = fs.join(self.output_path, str(self.counter) + ".png")

            # Debugging
            log.debug("Saving the plot to '" + path + "' ...")

            # Save
            img.save(path, "PNG")

        else: raise RuntimeError("Cannot show the plot, specify an output path")

# -----------------------------------------------------------------
