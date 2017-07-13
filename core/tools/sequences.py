#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.sequences Provides functions for dealing with sequences (lists, iterables).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import operator
import itertools
from functools import partial

# -----------------------------------------------------------------

def permutate(lst, lengths=None):

    """
    This function ...
    :param lst: 
    :param lengths: 
    :return: 
    """

    if lengths is None: return itertools.permutations(lst)
    elif isinstance(lengths, int): return itertools.permutations(lst, r=lengths)
    elif isinstance(lengths, list):
        result = []
        for length in lengths: result.extend(itertools.permutations(lst, r=length))
        return result
    else: raise ValueError("Invalid value for 'lengths'")

# -----------------------------------------------------------------

def combinations(lst, lengths, repeat=False):

    """
    This function ...
    :param lst: 
    :param lengths: 
    :param repeat:
    :return: 
    """

    combination_function = itertools.combinations_with_replacement if repeat else itertools.combinations

    if isinstance(lengths, int): result = combination_function(lst, r=lengths) # tuples
    elif isinstance(lengths, list):
        result = []
        for length in lengths: result.extend(combination_function(lst, r=length)) # tuples
    else: raise ValueError("Invalid value for 'lengths")

    # Convert tuples to lists
    lists = []
    for item in result: lists.append(list(item)) # CONVERT FROM TUPLE TO LIST
    return lists

# -----------------------------------------------------------------

def iterate_lists_combinations(*lsts):

    """
    This function ...
    :param lsts: 
    :return: 
    """

    return itertools.product(*lsts)

# -----------------------------------------------------------------

def lists_combinations(*lsts):

    """
    This function ...
    :param lsts: 
    :return: 
    """

    return list(iterate_lists_combinations(*lsts))

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

def combine_unique(*args):

    """
    This function ...
    :return:
    """

    unique = set()

    # Add
    for seq in args: unique |= set(seq)

    return list(unique)

# -----------------------------------------------------------------

def find_closest_index(seq, value):

    """
    This function ...
    :param seq:
    :param value:
    :return:
    """

    closest_delta = None
    #closest_delta = float("inf")
    closest_index = None

    #column_unit = table[column_name].unit

    #value_unit = value.unit if hasattr(value, "unit") else None

    # Check units
    #if value_unit is not None:
        #if column_unit is None: raise ValueError("Value has a unit but column has not: cannot compare these values")
        #else: value = value.to(column_unit).value # for correct comparison inside loop
    #elif column_unit is not None: raise ValueError("Value has no unit but the column has: cannot compare these values")

    # Loop over all entries in the sequence
    for i in range(len(seq)):

        delta = abs(seq[i] - value)

        if closest_delta is None or delta < closest_delta:
            closest_delta = delta
            closest_index = i

    return closest_index

# -----------------------------------------------------------------

def find_exact_index(seq, value):

    """
    This function ...
    :param seq:
    :param value:
    :return:
    """

    return seq.index(value)

# -----------------------------------------------------------------

def all_equal(lst, ignore_none=False):

    """
    This function ...
    :param lst:
    :param ignore_none:
    :return:
    """

    if len(lst) == 0: raise ValueError("Cannot use empty list")

    first = lst[0]

    if first is None and ignore_none:
        try: first = find_first_not_none(lst)
        except: raise ValueError("Cannot use empty list (except for Nones)") #return True # ALL NONE, SO ALL EQUAL

    #print(first)

    for index in range(len(lst)):

        # Ignore None?
        if ignore_none and lst[index] is None: continue

        #print("comparing:")
        #print(first)
        #print(lst[index])

        #if not (lst[index] == first):
        if lst[index] != first:
            #print("HEERE")
            return False

    return True

# -----------------------------------------------------------------

def get_all_equal_value(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    if not all_equal(sequence): raise ValueError("Not all equal")
    else: return sequence[0]

# -----------------------------------------------------------------

def get_first_not_none_value(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    for item in sequence:
        if item is not None: return item
    return None

# -----------------------------------------------------------------

def all_close(lst, ignore_none=False, rtol=1.e-5, atol=1.e-8):

    """
    This fnction ...
    :param lst:
    :param ignore_none:
    :param rtol:
    :param atol:
    :return:
    """

    import numpy as np

    if len(lst) == 0: raise ValueError("Cannot use empty list")

    first = lst[0]

    if first is None and ignore_none:
        try: first = find_first_not_none(lst)
        except: raise ValueError("Cannot use empty list (except for Nones)")

    # Get first value
    if hasattr(first, "unit"): first_value = first.to(first.unit).value
    else: first_value = first

    for index in range(len(lst)):

        # Ignore None?
        if ignore_none and lst[index] is None: continue

        # If quantities, get scalar value
        if hasattr(lst[index], "unit"): value = lst[index].to(first.unit).value
        else: value = lst[index]

        #print("comparing:")
        #print(value)
        #print(first_value)

        if not np.isclose(value, first_value, rtol=rtol, atol=atol): return False

    # Return
    return True

# -----------------------------------------------------------------

def all_true(lst, ignore_none=False):

    """
    This function ...
    :param lst:
    :param ignore_none:
    :return:
    """

    if len(lst) == 0: raise ValueError("Cannot use empty list")

    indices = not_none_indices(lst) if ignore_none else range(len(lst))

    if len(lst) == 0: raise ValueError("List only contains None values")

    for index in indices:
        if not lst[index]: return False
    return True

# -----------------------------------------------------------------

def find_first_not_none(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    for item in lst:
        if item is not None: return item

    # Shouldn't get here
    raise ValueError("No not-None values")

# -----------------------------------------------------------------

def not_none_indices(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    indices = []
    for index in range(len(lst)):
        if lst[index] is not None: indices.append(index)
    return indices

# -----------------------------------------------------------------

def equal_sizes(*lsts):

    """
    THis function ...
    :param lsts: 
    :return: 
    """

    size = None
    for lst in lsts:
        if size is None: size = len(lst)
        elif size != len(lst): return False
    return True

# -----------------------------------------------------------------

def zip_into_dict(list_a, list_b):

    """
    This function ...
    :param list_a:
    :param list_b:
    :return:
    """

    dictionary = dict()
    for item_a, item_b in zip(list_a, list_b): dictionary[item_a] = item_b
    return dictionary

# -----------------------------------------------------------------

def contains_duplicates(sequence):

    """
    This function ...
    :param sequence: 
    :return: 
    """

    return len(set(sequence)) != len(sequence)

# -----------------------------------------------------------------

def contains_same_elements(sequence_a, sequence_b):

    """
    This function ...
    :param sequence_a: 
    :param sequence_b: 
    :return: 
    """

    return set(sequence_a) == set(sequence_b)

# -----------------------------------------------------------------

def elements_not_in_other(sequence_a, sequence_b):

    """
    This function ...
    :param sequence_a: 
    :param sequence_b: 
    :return: 
    """

    elements = set()

    for element in sequence_a:

        if element not in sequence_b: elements.add(element)

    return list(elements)

# -----------------------------------------------------------------

def common_elements(sequence_a, sequence_b):

    """
    This function ...
    :param sequence_a: 
    :param sequence_b: 
    :return: 
    """

    elements = set()

    for element in sequence_a:

        if element in sequence_b: elements.add(element)

    return list(elements)

# -----------------------------------------------------------------

def append_unique(lst, element):

    """
    This function ...
    :param lst: 
    :param element: 
    :return: 
    """

    if element not in lst: lst.append(element)

# -----------------------------------------------------------------

def extend_unique(lst, elements):

    """
    This function ...
    :param lst: 
    :param elements: 
    :return: 
    """

    for element in elements: append_unique(lst, element)

# -----------------------------------------------------------------

def imerge(a, b):

    """
    This function ...
    :param a: 
    :param b: 
    :return: 
    """

    for i, j in itertools.izip_longest(a,b):
        yield i
        if j is not None:
            yield j

# -----------------------------------------------------------------

def iterate_from_middle(lst):

    """
    This function ...
    :param lst: 
    :return: 
    """

    try:

        middle = len(lst)/2
        yield lst[middle]

        for shift in range(1, middle+1):

            # order is important!
            yield lst[middle - shift]
            yield lst[middle + shift]

    # occures on lst[len(lst)] or for empty list
    except IndexError: raise StopIteration

# -----------------------------------------------------------------

def rearrange_from_middle(lst):

    """
    This function ...
    :param lst: 
    :return: 
    """

    return list(iterate_from_middle(lst))

# -----------------------------------------------------------------

def multiply_all(lst):

    """
    This function ...
    :param lst: 
    :return: 
    """

    #print("lst", lst)

    result = 1.
    for element in lst: result *= element
    return result

# -----------------------------------------------------------------

def multiply_all_integers(lst):

    """
    THis function ...
    :param lst: 
    :return: 
    """

    result = 1
    for element in lst: result *= int(element)
    return result

# -----------------------------------------------------------------

def all_in(sequence, target):

    """
    Thisf unction ...
    :param sequence: 
    :param target: 
    :return: 
    """

    for element in sequence:
        if element not in target: return False
    return True

# -----------------------------------------------------------------

def any_in(sequence, target):

    """
    This function ...
    :param sequence: 
    :param target: 
    :return: 
    """

    for element in sequence:
        if element in target: return True
    return False

# -----------------------------------------------------------------

def is_empty(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return len(sequence) == 0

# -----------------------------------------------------------------

def any_empty(*sequences):

    """
    This function ...
    :param sequences:
    :return:
    """

    return any(is_empty(sequence) for sequence in sequences)

# -----------------------------------------------------------------

def find_indices(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    return [index for index, value in enumerate(sequence) if value == element]

# -----------------------------------------------------------------

def find_unique(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    indices = find_indices(sequence, element)
    if len(indices) == 0: raise ValueError("Not found")
    elif len(indices) > 1: raise ValueError("Not unique")
    else: return indices[0]

# -----------------------------------------------------------------

def pack(*sequences):

    """
    This function ...
    :param sequences:
    :return:
    """

    return [list(element) for element in zip(*sequences)]

# -----------------------------------------------------------------

def unpack(zipped, default_size=None):

    """
    This function ...
    :param zipped:
    :param default_size:
    :return:
    """

    if len(zipped) == 0:
        if default_size is None: raise ValueError("Empty input")
        else: return [[] for _ in range(default_size)]

    first_element = zipped[0]
    sequences = [[] for _ in range(len(first_element))]

    for element in zipped:
        for index in range(len(element)): sequences[index].append(element[index])

    # Return the sequences
    return sequences

# -----------------------------------------------------------------

def is_sorted(sequence, invert=False):

    """
    This function ...
    :param sequence:
    :param invert:
    :return:
    """

    if invert: return all(sequence[i] >= sequence[i+1] for i in xrange(len(sequence)-1))
    else: return all(sequence[i] <= sequence[i+1] for i in xrange(len(sequence)-1))

# -----------------------------------------------------------------

def is_minimum(sequence, value):

    """
    This function ...
    :param sequence:
    :param value:
    :return:
    """

    for item in sequence:
        if item < value: return False
    else: return True

# -----------------------------------------------------------------

def is_maximum(sequence, value):

    """
    This function ...
    :param sequence:
    :param value:
    :return:
    """

    for item in sequence:
        if item > value: return False
    else: return True

# -----------------------------------------------------------------

def same_contents(sequence_a, sequence_b):

    """
    This function ...
    :param sequence_a:
    :param sequence_b:
    :return:
    """

    if len(sequence_a) != len(sequence_b): return False

    for item in sequence_a:

        if item not in sequence_b: return False

    return True

# -----------------------------------------------------------------

def find_differences(sequence_a, sequence_b):

    """
    This function ...
    :param sequence_a:
    :param sequence_b:
    :return:
    """

    # Sort both: NO, THEN WHAT DO THE INDICES MEAN?
    #sequence_a = sorted(sequence_a)
    #sequence_b = sorted(sequence_b)

    indices = []

    for index in range(min(len(sequence_a), len(sequence_b))):
        if sequence_a[index] != sequence_b[index]: indices.append(index)

    return indices

# -----------------------------------------------------------------

def is_singleton(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return len(sequence) == 1

# -----------------------------------------------------------------

def get_singleton(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    if len(sequence) == 0: raise ValueError("Empty sequence")
    elif len(sequence) > 1: raise ValueError("Not a singleton")
    else: return sequence[0]

# -----------------------------------------------------------------

def in_all(item, sequences):

    """
    This function ...
    :param item:
    :param sequences:
    :return:
    """

    for sequence in sequences:
        if item not in sequence: return False
    return True

# -----------------------------------------------------------------

def in_some(item, sequences):

    """
    This function ...
    :param item:
    :param sequences:
    :return:
    """

    for sequence in sequences:
        if item in sequence: return True
    return False

# -----------------------------------------------------------------

def union(*sequences):
    
    """
    This function ...
    :param sequences:
    :return: 
    """

    if len(sequences) == 0: return []
    elif len(sequences) == 1: return sequences[0][:] # make copy of the list

    #else: return set().union(*sequences) # DOES PROBABLY NOT CONTAIN ORDER

    # OWN IMPLEMENTATION: quadratic complexity!
    else:

        elements = []

        # Just loop over each sequence sequentially
        for sequence in sequences:
            for item in sequence:
                if item not in elements: elements.append(item)

        # Return the list of elements
        return elements

# -----------------------------------------------------------------

def intersection(*sequences):

    """
    This function ...
    :return:
    """

    if len(sequences) == 0: return []
    elif len(sequences) == 1: return sequences[0][:] # make copy of the list
    else:
        first_sequence = sequences[0]
        other_sequences = sequences[1:]
        return [item for item in first_sequence if in_all(item, other_sequences)]

# -----------------------------------------------------------------

def prepend(sequence, item):

    """
    This function ...
    :param sequence:
    :param item:
    :return:
    """

    sequence.insert(0, item)

# -----------------------------------------------------------------
