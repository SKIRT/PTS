#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.filelock Multiprocessing lock using the file system.
#
# The FileLock class in this module implements a multiprocessing lock using the file system.
# This code was adjusted from:
#  https://github.com/ilastik/lazyflow/blob/master/lazyflow/utility/fileLock.py

# -----------------------------------------------------------------

import os
import sys
import time
import errno

# -----------------------------------------------------------------

## A file locking mechanism that has context-manager support so you can use it in a ``with`` statement.
class FileLock(object):

    class FileLockException(Exception):
        pass

    ## Prepare the file locker. Specify the file to lock and optionally the maximum timeout
    # and the delay between each attempt to lock.
    def __init__(self, lock_file_path, timeout=None, delay=1):
        self.is_locked = False
        self.lockfile = lock_file_path
        self.timeout = timeout
        self.delay = delay
        self._lock_file_contents = "PTS lock file"

    ## Acquire the lock, if possible. If the lock is in use, and `blocking` is False, return False.
    # Otherwise, check again every `self.delay` seconds until it either gets the lock or
    # exceeds `timeout` number of seconds, in which case it raises an exception.
    def acquire(self, blocking=True):
        start_time = time.time()
        while True:
            try:
                # Attempt to create the lockfile.
                # These flags cause os.open to raise an OSError if the file already exists.
                fd = os.open( self.lockfile, os.O_CREAT | os.O_EXCL | os.O_RDWR )
                with os.fdopen( fd, 'a' ) as f:
                    # Print some info about the current process as debug info for anyone who bothers to look.
                    f.write( self._lock_file_contents )
                break;
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if self.timeout is not None and (time.time() - start_time) >= self.timeout:
                    raise FileLock.FileLockException("Timeout occurred.")
                if not blocking:
                    return False
                time.sleep(self.delay)
        self.is_locked = True
        return True

    ## Get rid of the lock by deleting the lockfile.
    def release(self):
        self.is_locked = False
        os.unlink(self.lockfile)

    ## Activated when used in the with statement; automatically acquire a lock to be used in the with block.
    def __enter__(self):
        self.acquire()
        return self

    ## Activated at the end of the with statement; automatically release the lock.
    def __exit__(self, type, value, traceback):
        self.release()

    ## Make sure this ``FileLock`` instance doesn't leave a lock file lying around.
    def __del__(self):
        if self.is_locked:
            self.release()

# -----------------------------------------------------------------
