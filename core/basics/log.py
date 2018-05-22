#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.log Contains the global logger.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import sys
import logging
import warnings
import inspect
import importlib
from contextlib import contextmanager

# -----------------------------------------------------------------

# The global log = a singleton of this module
log = None

# -----------------------------------------------------------------

# Suppress Astropy warnings
# Import astronomical modules
try:
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
except ImportError: pass

# -----------------------------------------------------------------

def customwarn(message, category, filename, lineno, file=None, line=None):

    #sys.stdout.write(warnings.formatwarning(message, category, filename, lineno))
    log.warning(str(category.__name__) + ": " + str(message) + " [file:" + str(filename) + ", line:" + str(lineno) + "]")

warnings.showwarning = customwarn

# -----------------------------------------------------------------

#logging_levels = ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'FATAL',]
logging_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'START', 'SUCCESS']

# -----------------------------------------------------------------

# Add the 'START' and 'SUCCESS' log level
START = 24
SUCCESS = 25
logging.addLevelName(START, "START")
logging.addLevelName(SUCCESS, "SUCCESS")

# -----------------------------------------------------------------

def handle_exception(exc_type, exc_value, exc_traceback):

    if issubclass(exc_type, KeyboardInterrupt) or issubclass(exc_type, SystemExit):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    #log.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))
    log.error(str(exc_type.__name__) + ": {0}".format(str(exc_value)), exc_info=(exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

# -----------------------------------------------------------------

def _init_log(level=None):

    """
    This function ...
    :param level:
    :return:
    """

    global log

    orig_logger_cls = logging.getLoggerClass()
    logging.setLoggerClass(PTSLogger)

    # Create a new named logger
    try:

        log = logging.getLogger("pts")
        log._set_defaults()
        log.welcome()

    finally: logging.setLoggerClass(orig_logger_cls)

    # Set level?
    if level is not None: log.setLevel(level)

    # Return the logger
    return log

# -----------------------------------------------------------------

def setup_log(level="INFO", path=None, memory=False, show_origins=False):

    """
    This function ...
    :param level:
    :param path:
    :param memory:
    :param show_origins:
    :return:
    """

    # Set the log level
    log.setLevel(level)

    # Add log file
    if path is not None: log.add_log_file(path)

    # Memory logging
    if memory: log.show_memuse()

    # Show origin
    if show_origins: log.show_origins()

    # Return the logger
    return log

# -----------------------------------------------------------------

class PTSFormatter(logging.Formatter):

    """
    This class ...
    """

    timestamp = "%(asctime)s.%(msecs)03d"
    message = "%(message)s"
    module = "(%(module)s)"

    debug_char = "D"
    info_char = " "
    start_char = "-"
    success_char = "-"
    warning_char = "!"
    error_char = "*"

    # -----------------------------------------------------------------

    def __init__(self, fmt="%(levelno)s: %(msg)s", datefmt="%d/%m/%Y %H:%M:%S"):

        """
        This function ...
        :param fmt:
        :param datefmt:
        """

        # Call the constructor of the base class
        logging.Formatter.__init__(self, fmt, datefmt)

    # -----------------------------------------------------------------

    def format(self, record):

        """
        This function ...
        :param record:
        :return:
        """

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = self.timestamp + " " + self.debug_char + " " + self.message

        elif record.levelno == logging.INFO:
            self._fmt = self.timestamp + " " + self.info_char + " " + self.message

        elif record.levelno == logging._levelNames["SUCCESS"]:
            self._fmt = self.timestamp + " " + self.success_char + " " + self.message

        elif record.levelno == logging._levelNames["START"]:
            self._fmt = self.timestamp + " " + self.start_char + " " + self.message

        elif record.levelno == logging.WARNING:
            self._fmt = self.timestamp + " " + self.warning_char + " " + self.message

        elif record.levelno == logging.ERROR:
            self._fmt = self.timestamp + " " + self.error_char + " " + self.message

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        # Return the formatted result
        return result

# -----------------------------------------------------------------

class ColorizingStreamHandler(logging.StreamHandler):

    """
    This class ...
    """

    @property
    def is_tty(self):

        isatty = getattr(self.stream, 'isatty', None)
        return isatty and isatty()

    # -----------------------------------------------------------------

    def emit(self, record):

        """
        This function ...
        :param record:
        :return:
        """

        if record.levelno <= logging.INFO: stream = sys.stdout
        else: stream = sys.stderr

        try:

            message = self.format(record)
            #stream = self.stream

            #if not self.is_tty: stream.write(message)
            #else: self.output_colorized(message)

            #record.message = "{0} [{1:s}]".format(record.msg, record.origin)
            #print(": " + record.message, file=stream)
            #stream.write(record.getMessage())
            stream.write(message)
            stream.write(getattr(self, 'terminator', '\n'))
            self.flush()

            #stream.write(getattr(self, 'terminator', '\n'))
            #self.flush()

        except (KeyboardInterrupt, SystemExit): raise
        except: self.handleError(record)

    # -----------------------------------------------------------------

    def format(self, record):

        """
        This function ...
        :param record:
        :return:
        """

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

# Get the original logging class
Logger = logging.getLoggerClass()

# -----------------------------------------------------------------

class Conf(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This functino ...
        """

        self.level = "INFO"
        self.use_color = True
        self.log_file_format = "%(asctime)r, %(origin)r, %(levelname)r, %(message)r"
        self.show_origin = False
        self.memory = False

# -----------------------------------------------------------------

conf = Conf()

# -----------------------------------------------------------------

def find_current_module(depth=1, finddiff=False):
    """
    Determines the module/package from which this function is called.
    This function has two modes, determined by the ``finddiff`` option. it
    will either simply go the requested number of frames up the call
    stack (if ``finddiff`` is False), or it will go up the call stack until
    it reaches a module that is *not* in a specified set.
    Parameters
    ----------
    depth : int
        Specifies how far back to go in the call stack (0-indexed, so that
        passing in 0 gives back `astropy.utils.misc`).
    finddiff : bool or list
        If False, the returned ``mod`` will just be ``depth`` frames up from
        the current frame. Otherwise, the function will start at a frame
        ``depth`` up from current, and continue up the call stack to the
        first module that is *different* from those in the provided list.
        In this case, ``finddiff`` can be a list of modules or modules
        names. Alternatively, it can be True, which will use the module
        ``depth`` call stack frames up as the module the returned module
        most be different from.
    Returns
    -------
    mod : module or None
        The module object or None if the package cannot be found. The name of
        the module is available as the ``__name__`` attribute of the returned
        object (if it isn't None).
    Raises
    ------
    ValueError
        If ``finddiff`` is a list with an invalid entry.
    Examples
    --------
    The examples below assume that there are two modules in a package named
    ``pkg``. ``mod1.py``::
        def find1():
            from astropy.utils import find_current_module
            print find_current_module(1).__name__
        def find2():
            from astropy.utils import find_current_module
            cmod = find_current_module(2)
            if cmod is None:
                print 'None'
            else:
                print cmod.__name__
        def find_diff():
            from astropy.utils import find_current_module
            print find_current_module(0,True).__name__
    ``mod2.py``::
        def find():
            from .mod1 import find2
            find2()
    With these modules in place, the following occurs::
        >>> from pkg import mod1, mod2
        >>> from astropy.utils import find_current_module
        >>> mod1.find1()
        pkg.mod1
        >>> mod1.find2()
        None
        >>> mod2.find()
        pkg.mod2
        >>> find_current_module(0)
        <module 'astropy.utils.misc' from 'astropy/utils/misc.py'>
        >>> mod1.find_diff()
        pkg.mod1
    """

    frm = inspect.currentframe()
    for i in range(depth):
        frm = frm.f_back
        if frm is None:
            return None

    if finddiff:
        currmod = inspect.getmodule(frm)
        if finddiff is True:
            diffmods = [currmod]
        else:
            diffmods = []
            for fd in finddiff:
                if inspect.ismodule(fd):
                    diffmods.append(fd)
                elif isinstance(fd, str):
                    diffmods.append(importlib.import_module(fd))
                elif fd is True:
                    diffmods.append(currmod)
                else:
                    raise ValueError('invalid entry in finddiff')

        while frm:
            frmb = frm.f_back
            modb = inspect.getmodule(frmb)
            if modb not in diffmods: return modb
            frm = frmb
    else: return inspect.getmodule(frm)

# -----------------------------------------------------------------

try:
    unicode
    _unicode = True
except NameError:
    _unicode = False

# -----------------------------------------------------------------

class PTSLogRecord(logging.LogRecord):

    """
    This class ...
    """

    def getMessage(self):

        """
        Return the message for this LogRecord.
        Return the message for this LogRecord after merging any user-supplied
        arguments with the message.
        """

        if not _unicode: #if no unicode support...
            msg = str(self.msg)

        else:
            msg = self.msg
            if not isinstance(msg, basestring):
                try:
                    msg = str(self.msg)
                except UnicodeError:
                    msg = self.msg      #Defer encoding till later
        if self.args:
            msg = msg % self.args

        # Add memory
        if conf.memory and hasattr(self, "memuse"):
            msg = "{" + self.memuse + "} " + msg

        # Add origin
        if conf.show_origin:
            if hasattr(self, "origin"):
                origin = getattr(self, "origin")
                if origin is not None: #and origin != "unknown":
                    msg = "[" + origin + "] " + msg

        # Add section
        if hasattr(self, "sections"):
            sections = getattr(self, "sections")
            if sections is not None and len(sections) > 0:
                msg = " > ".join(sections) + " | " + msg

        # Return
        return msg

# -----------------------------------------------------------------

logging.LogRecord = PTSLogRecord

# -----------------------------------------------------------------

class PTSLogger(Logger):

    """
    This class is used to set up the PTS logging.
    The main functionality added by this class over the built-in
    logging.Logger class is the ability to keep track of the origin of the
    messages, the ability to enable logging of warnings.warn calls and
    exceptions, and the addition of colorized output and context managers to
    easily capture messages to a file or list.
    """

    def __init__(self, *args, **kwargs):

        """
        Initialize the logger with a name and an optional level.
        """

        # Call the constructor of the base class
        super(PTSLogger, self).__init__(*args, **kwargs)

        # Attributes
        self.sections = None

    # -----------------------------------------------------------------

    # Following code from http://stackoverflow.com/a/938800/819110:
    _proc_status = '/proc/%d/status' % os.getpid()
    _scale = {'kB': 1024.0, 'mB': 1024.0 * 1024.0,
              'KB': 1024.0, 'MB': 1024.0 * 1024.0}

    def _VmB(self, VmKey):

        """Private.
        """
        # get pseudo file  /proc/<pid>/status

        try:
            t = open(self._proc_status)
            v = t.read()
            t.close()

        except: return 0.0  # non-Linux?

        # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # whitespace
        if len(v) < 3:
            return 0.0  # invalid format?
        # convert Vm value to bytes
        return float(v[1]) * self._scale[v[2]]

    # -----------------------------------------------------------------

    def memory(self,since=0.0):

        """
        Return memory usage in bytes.
        """

        return self._VmB('VmSize:') - since

    # -----------------------------------------------------------------

    def swapsize(self,since=0.0):

        """
        Return swap size in bytes.
        """

        return self._VmB('VmSwap:') - since

    # -----------------------------------------------------------------

    def byte_to_mb(self,byte):

        """
        Return size in MB (being lazy)
        """

        return byte/(1024*1024)

    # -----------------------------------------------------------------

    def str_mem(self):

        """Return a string with the total memuse and swap size in MB
        """

        return "Total:%.0fM,Swap:%.0fM"%(self.byte_to_mb(self.memory()),self.byte_to_mb(self.swapsize()) )

    # -----------------------------------------------------------------

    def makeRecord(self, name, level, pathname, lineno, msg, args, exc_info,
                   func=None, extra=None): #, sinfo=None):

        """
        This function ...
        :param name:
        :param level:
        :param pathname:
        :param lineno:
        :param msg:
        :param args:
        :param exc_info:
        :param func:
        :param extra:
        :return:
        """

        if extra is None: extra = {}

        # Add origin
        if 'origin' not in extra:

            current_module = find_current_module(1, finddiff=[True, 'logging'])
            if current_module is not None: extra['origin'] = current_module.__name__
            else: extra['origin'] = 'unknown'

        # Add memory usage
        if conf.memory and 'memory' not in extra:
            extra["memuse"] = self.str_mem()

        # Add section
        if self.sections is not None: extra["sections"] = self.sections

        # Call the base implementation
        return Logger.makeRecord(self, name, level, pathname, lineno, msg,
                                 args, exc_info, func=func, extra=extra)
                                 #sinfo=sinfo)

    # -----------------------------------------------------------------

    def show_origins(self):

        """
        This function ...
        :return:
        """

        conf.show_origin = True

    # -----------------------------------------------------------------

    def show_memuse(self):

        """
        This function ...
        :return:
        """

        conf.memory = True

    # -----------------------------------------------------------------

    @property
    def subsection(self):

        """
        This function ...
        :return:
        """

        if self.sections is None: return None
        elif len(self.sections) == 0: return None
        else: return self.sections[-1]

    # -----------------------------------------------------------------

    @property
    def supersection(self):

        """
        This function ...
        :return:
        """

        if self.sections is None: return None
        elif len(self.sections) == 0: return None
        else: return self.sections[0]

    # -----------------------------------------------------------------

    def add_subsection(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.sections is None: self.sections = (name,)
        elif self.subsection == name: self.warning("Already a log subsection '" + name + "': not adding ...")
        else: self.sections = self.sections + (name,)

    # -----------------------------------------------------------------

    def remove_subsection(self):

        """
        This function ...
        :return:
        """

        self.sections = self.sections[:-1]

    # -----------------------------------------------------------------

    def add_supersection(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.sections is None: self.sections = (name,)
        elif self.supersection == name: self.warning("Already a log supersection '" + name + "': not adding ...")
        else: self.sections = (name,) + self.sections

    # -----------------------------------------------------------------

    def remove_supersection(self):

        """
        This function ...
        :return:
        """

        self.sections = self.sections[1:]

    # -----------------------------------------------------------------

    def set_section(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        self.sections = (name,)

    # -----------------------------------------------------------------

    def set_sections(self, *sections):

        """
        This function ...
        :param sections:
        :return:
        """

        self.sections = sections

    # -----------------------------------------------------------------

    def add_log_file(self, path, level="DEBUG"):

        """
        This function ...
        :param path:
        :param level:
        :return:
        """

        # Create file handler
        fh = logging.FileHandler(path)

        # Set the formatter
        fh.setFormatter(self.handlers[0].formatter)

        # Set the level
        conf.level = level
        fh.setLevel(level)

        # Add the handler to the log instance
        self.addHandler(fh)

        # Return the handler
        return fh

    # -----------------------------------------------------------------

    def remove_all_log_files(self):

        """
        This function ...
        :return:
        """

        # Remove all file handlers, but keep the others
        self.handlers = [h for h in self.handlers if not isinstance(h, logging.FileHandler)]

    # -----------------------------------------------------------------

    def get_filehandlers(self):

        """
        This function ...
        :return:
        """

        return [h for h in self.handlers if isinstance(h, logging.FileHandler)]

    # -----------------------------------------------------------------

    def get_handler(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        for fh in self.get_filehandlers():
            if self.get_filepath(fh) == path: return fh
        return None

    # -----------------------------------------------------------------

    def get_filepath(self, fh):

        """
        This function ...
        :param fh:
        :return:
        """

        return fh.baseFilename

    # -----------------------------------------------------------------

    def get_filepaths(self):

        """
        This function ...
        :return:
        """

        paths = []

        for fh in self.handlers:
            if not isinstance(fh, logging.FileHandler): continue
            paths.append(self.get_filepath(fh))

        return paths

    # -----------------------------------------------------------------

    def start(self, message, *args, **kwargs):

        """
        This function ...
        :param message:
        :param args:
        :param kwargs:
        :return:
        """

        if self.isEnabledFor(START): self._log(START, message, args, **kwargs)

    # -----------------------------------------------------------------

    def success(self, message, *args, **kwargs):

        """
        This function ...
        :param message:
        :param args:
        :param kwargs:
        :return:
        """

        if self.isEnabledFor(SUCCESS): self._log(SUCCESS, message, args, **kwargs)

    # -----------------------------------------------------------------

    def welcome(self):

        """
        This function ...
        :return:
        """

        # Show welcome message
        self.info("Welcome to PTS")

    # -----------------------------------------------------------------

    @property
    def is_debug(self):

        """
        This function ...
        :return:
        """

        return self.level <= 10

    # -----------------------------------------------------------------

    def _set_defaults(self):

        """
        This function ...
        :return:
        """

        # Create formatter
        # formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(message)s (%(module)s)", "%d/%m/%Y %H:%M:%S")
        formatter = PTSFormatter()

        # Create handler
        handler = ColorizingStreamHandler()
        handler.setFormatter(formatter)
        self.addHandler(handler)

        # Set default level
        self.setLevel(conf.level)

    # -----------------------------------------------------------------

    @contextmanager
    def suppress(self):

        """
        This function ...
        :return:
        """

        original_level = log.level
        log.setLevel("WARNING")
        yield
        log.setLevel(original_level)

    # -----------------------------------------------------------------

    @contextmanager
    def no_debugging(self):

        """
        This function ...
        :return:
        """

        original_level = log.level
        self.setLevel("INFO")
        yield
        log.setLevel(original_level)

    # -----------------------------------------------------------------

    @contextmanager
    def log_to_file(self, filename, filter_level=None, filter_origin=None):

        """
        Context manager to temporarily log messages to a file.
        Parameters
        ----------
        filename : str
            The file to log messages to.
        filter_level : str
            If set, any log messages less important than ``filter_level`` will
            not be output to the file. Note that this is in addition to the
            top-level filtering for the logger, so if the logger has level
            'INFO', then setting ``filter_level`` to ``INFO`` or ``DEBUG``
            will have no effect, since these messages are already filtered
            out.
        filter_origin : str
            If set, only log messages with an origin starting with
            ``filter_origin`` will be output to the file.
        Notes
        -----
        By default, the logger already outputs log messages to a file set in
        the Astropy configuration file. Using this context manager does not
        stop log messages from being output to that file, nor does it stop log
        messages from being printed to standard output.
        Examples
        --------
        The context manager is used as::
            with logger.log_to_file('myfile.log'):
                # your code here
        """

        fh = logging.FileHandler(filename)

        if filter_level is not None: fh.setLevel(filter_level)
        #fh.setLevel(level)
        if filter_origin is not None: fh.addFilter(FilterOrigin(filter_origin))

        f = logging.Formatter(conf.log_file_format)
        fh.setFormatter(f)
        self.addHandler(fh)

        yield

        fh.close()
        self.removeHandler(fh)

    # -----------------------------------------------------------------

    @contextmanager
    def log_to_list(self, filter_level=None, filter_origin=None):

        '''
        Context manager to temporarily log messages to a list.
        Parameters
        ----------
        filename : str
            The file to log messages to.
        filter_level : str
            If set, any log messages less important than ``filter_level`` will
            not be output to the file. Note that this is in addition to the
            top-level filtering for the logger, so if the logger has level
            'INFO', then setting ``filter_level`` to ``INFO`` or ``DEBUG``
            will have no effect, since these messages are already filtered
            out.
        filter_origin : str
            If set, only log messages with an origin starting with
            ``filter_origin`` will be output to the file.
        Notes
        -----
        Using this context manager does not stop log messages from being
        output to standard output.
        Examples
        --------
        The context manager is used as::
            with logger.log_to_list() as log_list:
                # your code here
        '''

        lh = ListHandler()

        if filter_level is not None: lh.setLevel(filter_level)

        if filter_origin is not None: lh.addFilter(FilterOrigin(filter_origin))

        self.addHandler(lh)

        yield lh.log_list

        self.removeHandler(lh)

# -----------------------------------------------------------------

class FilterOrigin:

    """
    A filter for the record origin
    """

    def __init__(self, origin):

        """
        This function ...
        :param origin:
        """

        self.origin = origin

    # -----------------------------------------------------------------

    def filter(self, record):

        """
        This function ...
        :param record:
        :return:
        """

        return record.origin.startswith(self.origin)

# -----------------------------------------------------------------

class ListHandler(logging.Handler):

    """
    A handler that can be used to capture the records in a list
    """

    def __init__(self, filter_level=None, filter_origin=None):

        """
        This function ...
        :param filter_level:
        :param filter_origin:
        """

        logging.Handler.__init__(self)
        self.log_list = []

    # -----------------------------------------------------------------

    def emit(self, record):

        """
        Thisn function ...
        :param record:
        :return:
        """

        self.log_list.append(record)

# -----------------------------------------------------------------

# Initialize the global logger
_init_log()

# -----------------------------------------------------------------
