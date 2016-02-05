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
import logging

# -----------------------------------------------------------------

log = None

def setup_custom_logger(level="INFO"):

    global log

    log = logging.getLogger("pts")

    formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(message)s (%(module)s)", "%d/%m/%Y %H:%M:%S")

    #handler = logging.StreamHandler()
    handler = ColorizingStreamHandler()
    handler.setFormatter(formatter)

    #log = logging.getLogger(name)
    log.setLevel(level)
    log.addHandler(handler)

    return log

# -----------------------------------------------------------------

# Add the 'SUCCESS' log level
#SUCCESS = 25
#logging.addLevelName(SUCCESS, "SUCCESS")

#Logger = logging.getLoggerClass()

#def success(self, message, *args, **kwargs):
#    print(self.getEffectiveLevel())
#    if self.isEnabledFor(SUCCESS): self._log(SUCCESS, message, args, **kwargs)

#logging.Logger.success = success

#class PTSLogger(Logger):
    #def success(self, message, *args, **kwargs):
        #if self.isEnabledFor(SUCCESS): self._log(SUCCESS, message, args, **kwargs)

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
        #logging._levelNames["SUCCESS"]: (None, 'yellow', False),
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

# Global logger
#log = logging.getLogger("pts")
#log = None
#sh = ColorizingStreamHandler(sys.stdout)
#sh.setLevel(logging.INFO)
#formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s (%(origin)s)", "%d/%m/%Y %H:%M:%S")
#sh.setFormatter(formatter)
#log.addHandler(sh)

# -----------------------------------------------------------------

def init_log(name=None, level="INFO", path=None, memory=False):

    """
    Initializes the PTS log--in most circumstances this is called
    when a pts do script is launched. Parameters:
    :param name:
    :param level:
    :param path:
    """

    global log

    # Set the name to 'pts' if none was given
    if name is None: name = "pts"

    # Create a logger
    #log = logging.getLogger(name) # Does not provide a functioning log ??
    #log = logger.AstropyLogger(name)
    #log = logging.Logger(name)

    # Create a stream handler
    sh = ColorizingStreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)

    # Create and set the formatter
    formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(message)s (%(origin)s)", "%d/%m/%Y %H:%M:%S")
    sh.setFormatter(formatter)

    # Add the stream handler
    log.addHandler(sh)

    # Add file handler if requested
    if path is not None:

        # Create file handler
        fh = logging.FileHandler(path)

        # Set the formatter
        fh.setFormatter(formatter)

        # Set the level
        fh.setLevel(level)

        # Add the handler to the log instance
        log.addHandler(fh)

    # Return the logger so that the do script that calls this function can immediately use it
    return log

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

    #memory_log_conf_path = os.path.join(inspection.pts_root_dir, "memlogging.conf")
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
