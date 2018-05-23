#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.analyser Contains the SimulationAnalyser class, used for analysing simulation output.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..launch.basicanalyser import BasicAnalyser, steps, extraction, plotting, misc
from ..launch.batchanalyser import BatchAnalyser
from ..test.scalinganalyser import ScalingAnalyser
from ..basics.log import log
from ..simulation.simulation import RemoteSimulation
from ..tools import filesystem as fs
from ..simulation.remote import get_simulation_for_host
from ..tools import sequences
from .options import extraction_names, plotting_names, misc_names
from ..tools import formatting as fmt
from ..tools import types
from ..tools.stringify import tostr

# -----------------------------------------------------------------

batch = "batch"
scaling = "scaling"
extra = "extra"

# -----------------------------------------------------------------

all_steps = steps + [batch, scaling]
all_steps_and_extra = all_steps + [extra]

# -----------------------------------------------------------------

def has_analysed(simulation, steps, features=None):

    """
    This function ...
    :param simulation:
    :param steps:
    :param features:
    :return:
    """

    # Check
    if types.is_string_type(steps): steps = [steps]
    elif types.is_string_sequence(steps): pass
    else: raise ValueError("Invalid value for 'steps'")

    # Only one step defined
    if len(steps) == 1:

        # Get the step
        step = steps[0]

        # Extraction
        if step == extraction:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, extraction_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(extraction_names) + "'")

                # Return
                return sequences.contains_any(simulation.analysed_extraction, features)

            # Features are not defined
            else: return sequences.has_any(simulation.analysed_extraction)

        # Plotting
        elif step == plotting:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, plotting_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(plotting_names) + "'")

                # Return
                return sequences.contains_any(simulation.analysed_plotting, features)

            # Features are not defined
            else: return sequences.has_any(simulation.analysed_plotting)

        # Misc
        elif step == misc:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, misc_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(misc_names) + "'")

                # Return
                return sequences.contains_any(simulation.analysed_misc, features)

            # Features are not defined
            else: return sequences.has_any(simulation.analysed_misc)

        # Batch
        elif step == batch:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

            # Return flag
            return simulation.analysed_batch

        # Scaling
        elif step == scaling:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

            # Return flag
            return simulation.analysed_scaling

        # Invalid
        else: raise ValueError("Invalid step: '" + step + "'")

    # Multiple steps defined
    else:

        # Check whether features are not defined
        if features is not None: raise ValueError("Features cannot be specified with multiple steps")

        # Reset extraction
        if extraction in steps and sequences.has_any(simulation.analysed_extraction): return True

        # Reset plotting
        if plotting in steps and sequences.has_any(simulation.analysed_plotting): return True

        # Reset misc
        if misc in steps and sequences.has_any(simulation.analysed_misc): return True

        # Reset batch
        if batch in steps and simulation.analysed_batch: return True

        # Reset scaling
        if scaling in steps and simulation.analysed_scaling: return True

        # Return False
        return False

# -----------------------------------------------------------------

def is_analysed(simulation, steps, features=None):

    """
    This function ...
    :param simulation:
    :param steps:
    :param features:
    :return:
    """

    # Check
    if types.is_string_type(steps): steps = [steps]
    elif types.is_string_sequence(steps): pass
    else: raise ValueError("Invalid value for 'steps'")

    # Only one step defined
    if len(steps) == 1:

        # Get the step
        step = steps[0]

        # Extraction
        if step == extraction:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, extraction_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(extraction_names) + "'")

                # Return
                return sequences.contains_all(simulation.analysed_extraction, features)

            # Features are not defined
            else: return sequences.contains_all(simulation.analysed_extraction, extraction_names)

        # Plotting
        elif step == plotting:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, plotting_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(plotting_names) + "'")

                # Return
                return sequences.contains_all(simulation.analysed_plotting, features)

            # Features are not defined
            else: return sequences.contains_all(simulation.analysed_plotting, plotting_names)

        # Misc
        elif step == misc:

            # Features are defined
            if features is not None:

                # Check features
                if not sequences.all_in(features, misc_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(misc_names) + "'")

                # Return
                return sequences.contains_all(simulation.analysed_misc, features)

            # Features are not defined
            else: return sequences.contains_all(simulation.analysed_misc, misc_names)

        # Batch
        elif step == batch:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

            # Return flag
            return simulation.analysed_batch

        # Scaling
        elif step == scaling:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

            # Return flag
            return simulation.analysed_scaling

        # Invalid
        else: raise ValueError("Invalid step: '" + step + "'")

    # Multiple steps defined
    else:

        # Check whether features are not defined
        if features is not None: raise ValueError("Features cannot be specified with multiple steps")

        # Reset extraction
        if extraction in steps and not sequences.contains_all(simulation.analysed_extraction, extraction_names): return False

        # Reset plotting
        if plotting in steps and not sequences.contains_all(simulation.analysed_plotting, plotting_names): return False

        # Reset misc
        if misc in steps and not sequences.contains_all(simulation.analysed_misc, plotting_names): return False

        # Reset batch
        if batch in steps and not simulation.analysed_batch: return False

        # Reset scaling
        if scaling in steps and not simulation.analysed_scaling: return False

        # Return True
        return True

# -----------------------------------------------------------------

def reanalyse_simulation(simulation, steps, features=None, not_steps=None, not_features=None, config=None, extra_configs=None):

    """
    This function ...
    :param simulation:
    :param steps:
    :param features:
    :param not_steps:
    :param not_features:
    :param config:
    :param extra_configs:
    :return:
    """

    # Check
    if types.is_string_type(steps): steps = [steps]
    elif types.is_string_sequence(steps): pass
    else: raise ValueError("Invalid value for 'steps'")

    # Only one step defined
    if len(steps) == 1:

        # Flag indicating whether this simulation has been analysed or not
        simulation.analysed = False

        # Get the step
        step = steps[0]

        # Extraction
        if step == extraction:

            # All features
            if features is None: simulation.analysed_extraction = []

            # Select features
            else:

                # Check features
                if not sequences.all_in(features, extraction_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(extraction_names) + "'")

                # Set features
                simulation.analysed_extraction = sequences.elements_not_in_other(extraction_names, features, check_existing=True)

        # Plotting
        elif step == plotting:

            # All features
            if features is None: simulation.analysed_plotting = []

            # Select features
            else:

                # Check features
                if not sequences.all_in(features, plotting_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(plotting_names) + "'")

                # Set features
                simulation.analysed_plotting = sequences.elements_not_in_other(plotting_names, features, check_existing=True)

        # Misc
        elif step == misc:

            # All features
            if features is None: simulation.analysed_misc = []

            # Select features
            else:

                # Check features
                if not sequences.all_in(features, misc_names): raise ValueError("Features contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(misc_names) + "'")

                # Set features
                simulation.analysed_misc = sequences.elements_not_in_other(misc_names, features, check_existing=True)

        # Batch
        elif step == batch:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

            # Set analysed_batch flag to False
            simulation.analysed_batch = False

        # Scaling
        elif step == scaling:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

            # Set analysed_scaling flag to False
            simulation.analysed_scaling = False

        # Extra
        elif step == extra:

            # All features (=classes)
            if features is None: simulation.unset_analysed_extra()

            # Select features
            else:

                # Check classes
                if not sequences.all_in(features, simulation.analyser_class_names): raise ValueError("Classes contains invalid value(s): '" + tostr(features) + "': can only contain '" + tostr(simulation.analyser_class_names) + "'")

                # Remove
                #for class_name in features: simulation.remove_
                simulation.analysed_extra = sequences.elements_not_in_other(simulation.analyser_class_names, features, check_existing=True)

        # Invalid
        else: raise ValueError("Invalid step: '" + step + "'")

    # Multiple steps defined
    else:

        # Check whether features are not defined
        if features is not None: raise ValueError("Features cannot be specified with multiple steps")

        # Flag indicating whether this simulation has been analysed or not
        simulation.analysed = False

        # Reset extraction
        if extraction in steps: simulation.analysed_extraction = []

        # Reset plotting
        if plotting in steps: simulation.analysed_plotting = []

        # Reset misc
        if misc in steps: simulation.analysed_misc = []

        # Reset batch
        if batch in steps: simulation.analysed_batch = False

        # Reset scaling
        if scaling in steps: simulation.analysed_scaling = False

    # Now analyse the simulation
    analyse_simulation(simulation, not_steps=not_steps, not_features=not_features, config=config, extra_configs=extra_configs)

# -----------------------------------------------------------------

def analyse_simulation(simulation, not_steps=None, not_features=None, config=None, extra_configs=None):

    """
    This function ...
    :param simulation:
    :param not_steps:
    :param not_features:
    :param config:
    :param extra_configs:
    :return:
    """

    # Steps are defined
    if not_steps is not None:

        # One step is defined
        if len(not_steps) == 1:

            # Get the step
            not_step = not_steps[0]

            # Extraction
            if not_step == extraction:

                # Features are not specified
                if not_features is None: simulation.analysed_extraction = extraction_names

                # Features are defined
                else:

                    # Check features
                    if not sequences.all_in(not_features, extraction_names): raise ValueError("Features contains invalid value(s): '" + tostr(not_features) + "': can only contain '" + tostr(extraction_names) + "'")

                    # Add features
                    sequences.extend_unique(simulation.analysed_extraction, not_features)

            # Plotting
            elif not_step == plotting:

                # Features are not specified
                if not_features is None: simulation.analysed_plotting = plotting_names

                # Features are defined
                else:

                    # Check features
                    if not sequences.all_in(not_features, plotting_names): raise ValueError("Features contains invalid value(s): '" + tostr(not_features) + "': can only contain '" + tostr(plotting_names) + "'")

                    # Add features
                    sequences.extend_unique(simulation.analysed_plotting, not_features)

            # Misc
            elif not_step == misc:

                # Features are not specified
                if not_features is None: simulation.analysed_misc = misc_names

                # Features are defined
                else:

                    # Check features
                    if not sequences.all_in(not_features, misc_names): raise ValueError("Features contains invalid value(s): '" + tostr(not_features) + "': can only contain '" + tostr(misc_names) + "'")

                    # Add features
                    sequences.extend_unique(simulation.analysed_misc, not_features)

            # Batch
            elif not_step == batch:

                # Check whether features are not defined
                if not_features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

                # Set flag
                simulation.analysed_batch = True

            # Scaling
            elif not_step == scaling:

                # Check whether features are not defined
                if not_features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

                # Set flag
                simulation.analysed_scaling = True

            # Invalid
            else: raise ValueError("Invalid step '" + not_step + "'")

        # Multiple steps are defined
        else:

            # Check whether features are not defined
            if not_features is not None: raise ValueError("Features cannot be specified with multiple steps")

            # Set all extracted
            if extraction in not_steps: simulation.analysed_extraction = extraction_names

            # Set all plotted
            if plotting in not_steps: simulation.analysed_plotting = plotting_names

            # Set all misc
            if misc in not_steps: simulation.analysed_misc = misc_names

            # Set batch
            if batch in not_steps: simulation.analysed_batch = True

            # Set scaling
            if scaling in not_steps: simulation.analysed_scaling = True

    # Steps are not defined
    elif not_features is not None: raise ValueError("Cannot specify features when step is not defined")

    # Create simulation analyser
    analyser = SimulationAnalyser(config=config)

    # Run the analyser on the simulation
    analyser.run(simulation=simulation, extra_configs=extra_configs)

# -----------------------------------------------------------------

def show_analysis_steps(simulation, do_basic=True, do_batch=True, do_scaling=True, do_extra=True):

    """
    This function ...
    :param simulation:
    :param do_basic:
    :param do_batch:
    :param do_scaling:
    :param do_extra:
    :return:
    """

    # Set flags
    basic = do_basic
    batch = do_batch and simulation.from_batch
    scaling = do_scaling and simulation.from_scaling_test
    extra = do_extra

    print("")

    # Show basic analysis steps
    if basic:

        print(" - basic:")
        print("")
        show_basic_analysis_steps(simulation)
        print("")

    # Show batch
    if batch:

        if simulation.analysed_batch: print(fmt.green + " - batch: analysed" + fmt.reset)
        else: print(fmt.yellow + " - batch: not yet analysed" + fmt.reset)
        print("")

    # Show scaling
    if scaling:

        if simulation.analysed_scaling: print(fmt.green + " - scaling: analysed" + fmt.reset)
        else: print(fmt.yellow + " - scaling: not yet analysed" + fmt.reset)
        print("")

    # Show extra analysis steps
    if extra:

        print(" - extra")
        print("")
        show_extra_analysis_steps(simulation)
        print("")

# -----------------------------------------------------------------

def show_basic_analysis_steps(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    # Set flags
    extraction = simulation.analysis.any_extraction
    plotting = simulation.analysis.any_plotting
    miscellaneous = simulation.analysis.any_misc

    if extraction:

        print("    * extraction:")
        print("")
        show_extraction_steps(simulation)
        print("")

    if plotting:

        print("    * plotting:")
        print("")
        show_plotting_steps(simulation)
        print("")

    if miscellaneous:

        print("    * miscellaneous:")
        print("")
        show_misc_steps(simulation)
        print("")

# -----------------------------------------------------------------

def show_extraction_steps(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    from .options import progress_name, timeline_name, memory_name

    # Get options
    extraction_options = simulation.analysis.extraction

    # Set flags
    progress_extraction = extraction_options.progress
    timeline_extraction = extraction_options.timeline
    memory_extraction = extraction_options.memory

    # Extracted?
    extracted_progress = progress_name in simulation.analysed_extraction
    extracted_timeline =  timeline_name in simulation.analysed_extraction
    extracted_memory = memory_name in simulation.analysed_extraction

    if progress_extraction:
        if extracted_progress: print(fmt.green + "       + progress: extracted" + fmt.reset)
        else: print(fmt.yellow + "       + progress: not yet extracted" + fmt.reset)

    if timeline_extraction:
        if extracted_timeline: print(fmt.green + "       + timeline: extracted" + fmt.reset)
        else: print(fmt.yellow + "       + timeline: not yet extracted" + fmt.reset)

    if memory_extraction:
        if extracted_memory: print(fmt.green + "       + memory: extracted" + fmt.reset)
        else: print(fmt.yellow + "       + memory: not yet extracted" + fmt.reset)

# -----------------------------------------------------------------

def show_plotting_steps(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    from .options import progress_name, timeline_name, memory_name, seds_name, grids_name

    # Get options
    plotting_options = simulation.analysis.plotting

    # Set flags
    seds_plotting = plotting_options.seds
    grids_plotting = plotting_options.grids
    progress_plotting = plotting_options.progress
    timeline_plotting = plotting_options.timeline
    memory_plotting = plotting_options.memory

    # Plotted?
    plotted_seds = seds_name in simulation.analysed_plotting
    plotted_grids = grids_name in simulation.analysed_plotting
    plotted_progress = progress_name in simulation.analysed_plotting
    plotted_timeline = timeline_name in simulation.analysed_plotting
    plotted_memory = memory_name in simulation.analysed_plotting

    if seds_plotting:
        if plotted_seds: print(fmt.green + "       + seds: plotted" + fmt.reset)
        else: print(fmt.yellow + "       + seds: not yet plotted" + fmt.reset)

    if grids_plotting:
        if plotted_grids: print(fmt.green + "       + grids: plotted" + fmt.reset)
        else: print(fmt.yellow + "       + grids: not yet plotted" + fmt.reset)

    if progress_plotting:
        if plotted_progress: print(fmt.green + "       + progress: plotted" + fmt.reset)
        else: print(fmt.yellow + "       + progress: not yet plotted" + fmt.reset)

    if timeline_plotting:
        if plotted_timeline: print(fmt.green + "       + timeline: plotted" + fmt.reset)
        else: print(fmt.yellow + "       + timeline: not yet plotted" + fmt.reset)

    if memory_plotting:
        if plotted_memory: print(fmt.green + "       + memory: plotted" + fmt.reset)
        else: print(fmt.yellow + "       + memory: not yet plotted" + fmt.reset)

# -----------------------------------------------------------------

def show_misc_steps(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    from .options import rgb_name, animations_name, fluxes_name, fluxes_from_images_name, images_name

    # Get options
    misc_options = simulation.analysis.misc

    # Set flags
    rgb = misc_options.rgb
    animations = misc_options.animations
    observed_fluxes = misc_options.fluxes
    observed_fluxes_from_images = misc_options.fluxes_from_images
    observed_images = misc_options.images

    # Analysed?
    has_rgb = rgb_name in simulation.analysed_misc
    has_animations = animations_name in simulation.analysed_misc
    has_fluxes = fluxes_name in simulation.analysed_misc
    has_fluxes_from_images = fluxes_from_images_name in simulation.analysed_misc
    has_images = images_name in simulation.analysed_misc

    # Show

    if rgb:
        if has_rgb: print(fmt.green + "       + rgb: created" + fmt.reset)
        else: print(fmt.yellow + "       + rgb: not yet created" + fmt.reset)

    if animations:
        if has_animations: print(fmt.green + "       + animations: created" + fmt.reset)
        else: print(fmt.yellow + "       + animations: not yet created" + fmt.reset)

    if observed_fluxes:
        if has_fluxes: print(fmt.green + "       + fluxes: calculated" + fmt.reset)
        else: print(fmt.yellow + "       + fluxes: not yet calculated" + fmt.reset)

    if observed_fluxes_from_images:
        if has_fluxes_from_images: print(fmt.green + "       + fluxes from images: calculated" + fmt.reset)
        else: print(fmt.yellow + "       + fluxes from images: not yet calculated" + fmt.reset)

    if observed_images:
        if has_images: print(fmt.green + "       + images: created" + fmt.reset)
        else: print(fmt.yellow + "       + images: not yet created" + fmt.reset)

# -----------------------------------------------------------------

def show_extra_analysis_steps(simulation):

    """
    This function ....
    :param simulation:
    :return:
    """

    # Loop over the 'extra' analyser classes that are defined for this simulation
    for analyser_class in simulation.analyser_classes:

        # Get name
        class_name = analyser_class.__name__

        # Check
        if class_name in simulation.analysed_extra: print(fmt.green + "    * " + class_name + ": analysed" + fmt.reset)
        else: print(fmt.yellow + "    * " + class_name + ": not yet analysed" + fmt.reset)

# -----------------------------------------------------------------

class SimulationAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SimulationAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = None
        self.batch_analyser = None
        self.scaling_analyser = None

        # Configs for the extra analysers
        self.extra_configs = None

    # -----------------------------------------------------------------

    @property
    def analysed_batch(self):

        """
        This function ...
        :return:
        """

        return self.simulation.analysed_batch

    # -----------------------------------------------------------------

    @property
    def analysed_scaling(self):

        """
        This function ...
        :return:
        """

        return self.simulation.analysed_scaling

    # -----------------------------------------------------------------

    @property
    def basic(self):

        """
        This function ...
        :return:
        """

        return self.config.do_basic

    # -----------------------------------------------------------------

    @property
    def batch(self):

        """
        This function ...
        :return:
        """

        return self.config.do_batch and self.simulation.from_batch and not self.analysed_batch

    # -----------------------------------------------------------------

    @property
    def scaling(self):

        """
        This function ...
        :return:
        """

        return self.config.do_scaling and self.simulation.from_scaling_test and not self.analysed_scaling

    # -----------------------------------------------------------------

    @property
    def extra(self):

        """
        This function ...
        :return:
        """

        return self.config.do_extra

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. If the simulation has no analysis options, finish the procedure right away
        if self.simulation.analysis is None: return

        # 3. Run the basic analysis
        if self.basic: self.analyse_basic()

        # 4. Run the batch analysis
        if self.batch: self.analyse_batch()

        # 5. Analyse the scaling
        if self.scaling: self.analyse_scaling()

        # 6. Perform extra analysis
        if self.extra: self.analyse_extra()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationAnalyser, self).setup(**kwargs)

        # Make a local reference to the simulation object
        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        elif self.config.remote is not None and self.config.id is not None: self.load_simulation()
        else: raise ValueError("No simulation is specified")

        # Set the analysers
        self.basic_analyser = BasicAnalyser(self.config.basic)
        self.batch_analyser = BatchAnalyser(self.config.batch)
        self.scaling_analyser = ScalingAnalyser()

        # Set flags
        self.basic_analyser.config.ignore_missing_data = self.config.ignore_missing_data
        self.batch_analyser.config.ignore_missing_data = self.config.ignore_missing_data

        # Get extra configs
        if kwargs.get("extra_configs", None) is not None: self.extra_configs = kwargs.pop("extra_configs")
        else: self.extra_configs = dict()

    # -----------------------------------------------------------------

    def load_simulation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the simulation ...")

        # Load simulation
        self.simulation = get_simulation_for_host(self.config.remote, self.config.id)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the simulation analyser ...")

        # Set everything to None
        self.simulation = None
        self.basic_analyser = None
        self.batch_analyser = None
        self.scaling_analyser = None

    # -----------------------------------------------------------------

    def analyse_basic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the simulation output ...")

        # Run the analyser on the simulation
        self.basic_analyser.run(simulation=self.simulation)

    # -----------------------------------------------------------------

    def analyse_batch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the properties relevant for the batch of simulations ...")

        # Run the batch analyser on the simulation
        self.batch_analyser.run(simulation=self.simulation, timeline=self.basic_analyser.timeline, memory=self.basic_analyser.memory)

        # Set flag
        self.simulation.analysed_batch = True
        self.simulation.save()

    # -----------------------------------------------------------------

    def analyse_scaling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the scaling results ...")

        # Run the scaling analyser
        self.scaling_analyser.run(simulation=self.simulation, timeline=self.basic_analyser.timeline, memory=self.basic_analyser.memory)

        # Set flag
        self.simulation.analysed_scaling = True
        self.simulation.save()

    # -----------------------------------------------------------------

    def get_config_for_analyser_class(self, class_name):

        """
        This function ...
        :param class_name:
        :return:
        """

        if class_name not in self.extra_configs: return None
        else: return self.extra_configs[class_name]

    # -----------------------------------------------------------------

    def analyse_extra(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing extra analysis on the simulation output ...")

        # Loop over the 'extra' analyser classes that are defined for this simulation
        for analyser_class in self.simulation.analyser_classes:

            # Get name
            class_name = analyser_class.__name__

            # Check
            if class_name in self.simulation.analysed_extra:
                log.debug("Analysis with the " + class_name + " class has already been performed")
                continue

            # Debugging
            log.debug("Running the " + class_name + " on the simulation ...")

            # Create an instance of the analyser class
            analyser = analyser_class.for_simulation(self.simulation, config=self.get_config_for_analyser_class(class_name), check_required=False)

            # Run the analyser, giving this simulation analyser instance as an argument
            analyser.run(simulation_analyser=self)

            # Add name to analysed_extra
            self.simulation.analysed_extra.append(class_name)
            self.simulation.save()

        # Indicate that this simulation has been analysed
        self.simulation.analysed = True
        self.simulation.save()

        # If requested, remove the local output directory
        if isinstance(self.simulation, RemoteSimulation) and self.simulation.remove_local_output: fs.remove_directory(self.simulation.output_path)

# -----------------------------------------------------------------
