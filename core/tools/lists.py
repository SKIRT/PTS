#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.lists Provides functions for dealing with lists.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
from string import ascii_lowercase

# -----------------------------------------------------------------

def interleave(seqs):

    """ Interleave a sequence of sequences
    >>> list(interleave([[1, 2], [3, 4]]))
    [1, 3, 2, 4]
    >>> ''.join(interleave(('ABC', 'XY')))
    'AXBYC'
    Both the individual sequences and the sequence of sequences may be infinite
    Returns a lazy iterator
    """

    iters = itertools.cycle(map(iter, seqs))
    while True:
        try:
            for itr in iters:
                yield next(itr)
            return
        except StopIteration:
            predicate = partial(operator.is_not, itr)
            iters = itertools.cycle(itertools.takewhile(predicate, iters))

# -----------------------------------------------------------------
