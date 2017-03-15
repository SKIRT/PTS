#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.logging Provides functions for creating loggers and linking log files to them.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import sys
import types
import logging
import warnings

# -----------------------------------------------------------------

# Suppress Astropy warnings
# Import astronomical modules
try:
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
except ImportError: pass

# -----------------------------------------------------------------

# Add the 'START' log level
START = 24
logging.addLevelName(START, "START")
def start(self, message, *args, **kwargs):
    if self.isEnabledFor(START): self._log(START, message, args, **kwargs)
logging.Logger.start = start

# Add the 'SUCCESS' log level
SUCCESS = 25
logging.addLevelName(SUCCESS, "SUCCESS")
#Logger = logging.getLoggerClass()
def success(self, message, *args, **kwargs):
    if self.isEnabledFor(SUCCESS): self._log(SUCCESS, message, args, **kwargs)
logging.Logger.success = success

# -----------------------------------------------------------------

# The global log = a singleton of this module
log = None

# -----------------------------------------------------------------

#def exception_handler(type, value, tb):
#    log.exception(str(type.__name__) + ": {0}".format(str(value)))
#    print(tb)

# Set exception handler
#sys.excepthook = exception_handler

def handle_exception(exc_type, exc_value, exc_traceback):

    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    #log.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    log.error(str(exc_type.__name__) + ": {0}".format(str(exc_value)), exc_info=(exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

# -----------------------------------------------------------------

def init_log(level="INFO"):

    """
    This function ...
    :param level:
    :return:
    """

    global log

    # Create a new named logger
    log = logging.getLogger("pts")

    # Create formatter
    #formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(message)s (%(module)s)", "%d/%m/%Y %H:%M:%S")
    formatter = MyFormatter()

    # Create handler
    handler = ColorizingStreamHandler()
    handler.setFormatter(formatter)

    # Set level and stream handler
    log.setLevel(level)
    log.addHandler(handler)

    def is_debug(self): return self.level <= 10
    log.is_debug = types.MethodType(is_debug, log)

    # Show welcome message
    log.info("Welcome to PTS")

# -----------------------------------------------------------------

def setup_log(level="INFO", path=None, memory=False):

    """
    This function ...
    :param level:
    :param path:
    :param memory:
    :return:
    """

    log.setLevel(level)

    if path is not None:

        # Create file handler
        fh = logging.FileHandler(path)

        # Set the formatter
        fh.setFormatter(log.handlers[0].formatter)

        # Set the level
        fh.setLevel(level)

        # Add the handler to the log instance
        log.addHandler(fh)

    # Memory logging
    if memory: pass

    return log

# -----------------------------------------------------------------

# Custom formatter
class MyFormatter(logging.Formatter):

    timestamp = "%(asctime)s.%(msecs)03d"
    message = "%(message)s"
    module = "(%(module)s)"

    debug_char = "D"
    info_char = " "
    start_char = "-"
    success_char = "-"
    warning_char = "!"
    error_char = "*"

    def __init__(self, fmt="%(levelno)s: %(msg)s", datefmt="%d/%m/%Y %H:%M:%S"):
        logging.Formatter.__init__(self, fmt, datefmt)


    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.debug_char + " " + MyFormatter.message

        elif record.levelno == logging.INFO:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.info_char + " " + MyFormatter.message

        elif record.levelno == logging._levelNames["SUCCESS"]:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.success_char + " " + MyFormatter.message

        elif record.levelno == logging._levelNames["START"]:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.start_char + " " + MyFormatter.message

        elif record.levelno == logging.WARNING:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.warning_char + " " + MyFormatter.message

        elif record.levelno == logging.ERROR:
            self._fmt = MyFormatter.timestamp + " " + MyFormatter.error_char + " " + MyFormatter.message

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result

# -----------------------------------------------------------------

class ColorizingStreamHandler(logging.StreamHandler):

    @property
    def is_tty(self):
        isatty = getattr(self.stream, 'isatty', None)
        return isatty and isatty()

    # -----------------------------------------------------------------

    def emit(self, record):
        try:
            message = self.format(record)
            stream = self.stream
            if not self.is_tty:
                stream.write(message)
            else:
                self.output_colorized(message)
            stream.write(getattr(self, 'terminator', '\n'))
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

    # -----------------------------------------------------------------

    def output_colorized(self, message):
            self.stream.write(message)

    # -----------------------------------------------------------------

    def format(self, record):
        message = logging.StreamHandler.format(self, record)
        if self.is_tty:
            # Don't colorize any traceback
            parts = message.split('\n', 1)
            parts[0] = self.colorize(parts[0], record)
            message = '\n'.join(parts)
        return message

    # -----------------------------------------------------------------

    # color names to indices
    color_map = {
        'black': 0,
        'red': 1,
        'green': 2,
        'yellow': 3,
        'blue': 4,
        'magenta': 5,
        'cyan': 6,
        'white': 7,
    }

    #levels to (background, foreground, bold/intense)
    level_map = {
        logging.DEBUG: (None, 'blue', False),
        logging.INFO: (None, None, False),
        logging._levelNames["START"]: (None, None, True),
        logging._levelNames["SUCCESS"]: (None, 'green', False),
        logging.WARNING: (None, 'magenta', False),
        logging.ERROR: (None, 'red', False),
        logging.CRITICAL: ('red', 'white', True) }

    csi = '\x1b['
    reset = '\x1b[0m'

    def colorize(self, message, record):

        if record.levelno in self.level_map:
            bg, fg, bold = self.level_map[record.levelno]
            params = []
            if bg in self.color_map:
                params.append(str(self.color_map[bg] + 40))
            if fg in self.color_map:
                params.append(str(self.color_map[fg] + 30))
            if bold:
                params.append('1')
            if params:
                message = ''.join((self.csi, ';'.join(params),
                                   'm', message, self.reset))
        return message

# -----------------------------------------------------------------

class MemuseFilter(logging.Filter):

    def filter(self, record):
        """ This function overrides logging.Filter, adds memuse as a field
        """
        record.memuse = self.str_mem()
        return True

    # Following code from http://stackoverflow.com/a/938800/819110:
    _proc_status = '/proc/%d/status' % os.getpid()
    _scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
              'KB': 1024.0, 'MB': 1024.0*1024.0}

    def _VmB(self,VmKey):
        """Private.
        """
        # get pseudo file  /proc/<pid>/status
        try:
            t = open(self._proc_status)
            v = t.read()
            t.close()
        except:
            return 0.0  # non-Linux?
        # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # whitespace
        if len(v) < 3:
            return 0.0  # invalid format?
        # convert Vm value to bytes
        return float(v[1]) * self._scale[v[2]]

    def memory(self,since=0.0):
        """Return memory usage in bytes.
        """
        return self._VmB('VmSize:') - since

    def swapsize(self,since=0.0):
        """Return swap size in bytes.
        """
        return self._VmB('VmSwap:') - since

    def byte_to_mb(self,byte):
        """return size in MB (being lazy)
        """
        return byte/(1024*1024)

    def str_mem(self):
        """Return a string with the total memuse and swap size in MB
        """
        return "MemTotal:%.0fM,Swap:%.0fM"%(self.byte_to_mb(self.memory()),self.byte_to_mb(self.swapsize()) )

# -----------------------------------------------------------------

def new_memory_log():

    """
    This function ...
    :return:
    """

    #memory_log_conf_path = os.path.join(introspection.pts_root_dir, "memlogging.conf")
    #logging.basicConfig(memory_log_conf_path)

    #logging.basicConfig(format="%(asctime)-15s %(name)-5s %(levelname)-8s %(memuse)-22s %(message)s")

    log = logging.getLogger('')               # Get root logger

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)-15s %(name)-5s %(levelname)-8s %(memuse)-22s %(message)s")
    sh.setFormatter(formatter)

    log.addHandler(sh)

    f = MemuseFilter()                        # Create filter
    log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler

    return log

# -----------------------------------------------------------------

# Initialize a global logger
init_log()

# -----------------------------------------------------------------

def customwarn(message, category, filename, lineno, file=None, line=None):
    #sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
    log.warning(str(category.__name__) + ": " + str(message) + " [file:" + str(filename) + ", line:" + str(lineno) + "]")
warnings.showwarning = customwarn

# -----------------------------------------------------------------

class suppress_logging(object):

    def __enter__(self):

        self.original_level = log.level
        log.setLevel("WARNING")

    def __exit__(self, type, value, traceback):

        log.setLevel(self.original_level)

# -----------------------------------------------------------------

class no_debugging(object):

    def __enter__(self):

        self.original_level = log.level
        log.setLevel("INFO")

    def __exit__(self, type, value, traceback):

        log.setLevel(self.original_level)

# -----------------------------------------------------------------
