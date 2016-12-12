#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.python Open a python session with PTS modules loaded, and with remote capabilities if requested.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from code import InteractiveConsole

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.remote.host import find_host_ids
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.remote.remote import Remote
from pts.core.tools.logging import log

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()
definition.add_positional_optional("host_id", "string", "remote host ID", choices=find_host_ids())

# Read the command line arguments
setter = ArgumentConfigurationSetter("python", "Open a python session with PTS modules loaded, and with remote capabilities if requested.")
config = setter.run(definition)

# -----------------------------------------------------------------

class PTSConsole(InteractiveConsole):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        """

        # Call the constructor of the base class
        InteractiveConsole.__init__(self, *args)

        # The remote
        #self.remote = kwargs.pop("remote", None)

        self.host_id = kwargs.pop("host_id", None)

    # -----------------------------------------------------------------

    def runsource(self, source, filename="<input>", symbol="single"):

        """
        THis function ...
        :param source:
        :param filename:
        :param symbol:
        :return:
        """

        #if self.remote is not None:
        if self.host_id is not None:

            source = source.replace("Frame.from_file(", "RemoteFrame")

            if "Frame.from_file" in source:

                source = source.split(" = Frame")[0] + " = RemoteFrame.from_file(" + source.split("Frame(")[1].split(",")[0] + ", '" + self.host_id + "')"

                print(source)

            #source = source.replace("Image.from_file(", "RemoteImage")
            #source = source.replace("DataCube.from_file", "RemoteDataCube")
            #print(source)

        # Run the code
        InteractiveConsole.runsource(self, source, filename=filename, symbol=symbol)

# -----------------------------------------------------------------

# Console
#if config.host_id is not None:
#    remote = Remote()
#    remote.setup(config.host_id)
#else: remote = None
#console = PTSConsole(remote=remote)
console = PTSConsole(host_id=config.host_id)

log.info("Importing modules ...")

# Load PTS modules
console.runcode("from pts.core.tools import filesystem as fs")
console.runcode("from pts.core.tools import introspection")
console.runcode("from pts.core.tools.logging import log")
console.runcode("from pts.magic.core.frame import Frame")
console.runcode("from pts.magic.core.image import Image")
console.runcode("from pts.magic.core.datacube import DataCube")
console.runcode("from pts.magic.core.kernel import ConvolutionKernel")
console.runcode("from pts.magic.remote import RemoteFrame, RemoteImage, RemoteDataCube")

# Do the interaction
console.interact()

# -----------------------------------------------------------------
