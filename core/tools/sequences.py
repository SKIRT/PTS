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

def is_dictionary(value):

    """
    This function ...
    :param value:
    :return:
    """

    from .types import is_dictionary
    return is_dictionary(value)

# -----------------------------------------------------------------

def contains_dictionary(value):

    """
    This function ...
    :param value:
    :return:
    """

    for item in value:
        if is_dictionary(item): return True
    return False

# -----------------------------------------------------------------

def replace_dictionaries_by_tuples(value):

    """
    This function ...
    :param value:
    :return:
    """

    new = []
    for item in value:
        if is_dictionary(item): item = tuple(item.items())
        new.append(item)
    return new

# -----------------------------------------------------------------

def is_sequence(value):

    """
    This function ...
    :param value:
    :return:
    """

    from .types import is_sequence
    return is_sequence(value)

# -----------------------------------------------------------------

def contains_sequence(value):

    """
    This function ...
    :param value:
    :return:
    """

    for item in value:
        if is_sequence(item): return True
    return False

# -----------------------------------------------------------------

def replace_sequences_by_tuples(value):

    """
    This function ...
    :param value:
    :return:
    """

    new = []
    for item in value:
        if is_sequence(item): item = tuple(item)
        new.append(item)
    return new

# -----------------------------------------------------------------

def create_nested_2d(ni, nj, fill=None):

    """
    This function ...
    :param ni:
    :param nj:
    :param fill:
    :return:
    """

    return [[ fill for j in range(nj) ] for i in range(ni) ]

# -----------------------------------------------------------------

def create_nested_3d(ni, nj, nk, fill=None):

    """
    This function ...
    :param ni:
    :param nj:
    :param nk:
    :param fill:
    :return:
    """

    return [[[ fill for k in range(nk) ] for j in range(nj) ] for i in range(ni) ]

# -----------------------------------------------------------------

def iterate_2d(structure):

    """
    This function ...
    :param structure:
    :return:
    """

    for i in range(len(structure)):
        for j in range(len(structure[i])):
            yield structure[i][j]

# -----------------------------------------------------------------

def iterate_3d(structure):

    """
    This function ...
    :param structure:
    :return:
    """

    for i in range(len(structure)):
        for j in range(len(structure[i])):
            for k in range(len(structure[i][j])):
                yield structure[i][j][k]

# -----------------------------------------------------------------

def ordered(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return list(sorted(sequence))

# -----------------------------------------------------------------

def argsort(seq):

    """
    This function ...
    :param seq:
    :return:
    """

    # Determine sort key
    sortkey = seq.__getitem__

    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=sortkey)

# -----------------------------------------------------------------

def repeat(sequence, ntimes):

    """
    This function ...
    :param sequence:
    :param ntimes:
    :return:
    """

    return list(itertools.chain.from_iterable(itertools.repeat(x, ntimes) for x in sequence))

# -----------------------------------------------------------------

def before(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    elements = []
    for el in sequence:
        if el == element: break
        elements.append(el)
    else: raise ValueError("Element not in sequence")
    return elements

# -----------------------------------------------------------------

def indices_before(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    indices = []
    for index, el in enumerate(sequence):
        if el == element: break
        indices.append(index)
    else: raise ValueError("Element not in sequence")
    return indices

# -----------------------------------------------------------------

def before_and_including(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    return before(sequence, element) + [element]

# -----------------------------------------------------------------

def after(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    elements = []
    for el in reversed(sequence):
        if el == element: break
        elements.append(el)
    else: raise ValueError("Element not in sequence")
    return list(reversed(elements))

# -----------------------------------------------------------------

def after_and_including(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    return [element] + after(sequence, element)

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

def have_units(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    has_units = [hasattr(item, "unit") for item in sequence]
    if not all_equal(has_units): raise ValueError("Inconsistent sequence")
    return has_units[0]

# -----------------------------------------------------------------

def same_units(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return all_equal([item.unit for item in sequence])

# -----------------------------------------------------------------

def get_unit(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    if not same_units(sequence): raise ValueError("Not the same units")
    return sequence[0].unit

# -----------------------------------------------------------------

def without_units(sequence, unit=None, check_same=True):

    """
    This function ...
    :param sequence:
    :param unit:
    :param check_same:
    :return:
    """

    if unit is not None: return [item.to(unit).value for item in sequence]
    else: # no unit is specified
        if check_same:
            unit = get_unit(sequence)
            return [item.to(unit).value for item in sequence]
        else: return [item.value for item in sequence]

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

    # Loop over all entries in the sequence
    for i in range(len(seq)):

        delta = abs(seq[i] - value)

        if closest_delta is None or delta < closest_delta:
            closest_delta = delta
            closest_index = i

    return closest_index

# -----------------------------------------------------------------

def find_closest_value(seq, value):

    """
    This function ...
    :param seq:
    :param value:
    :return:
    """

    index = find_closest_index(seq, value)
    return seq[index]

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

def all_equal_to(lst, value, ignore=None, ignore_none=False):

    """
    This function ...
    :param lst:
    :param value:
    :param ignore:
    :param ignore_none:
    :return:
    """

    for item in lst:
        if ignore is not None and item == ignore: continue
        if ignore_none and item is None: continue
        if item != value: return False
    return True

# -----------------------------------------------------------------

def all_zero(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    return all_equal_to(lst, 0)

# -----------------------------------------------------------------

def has_equal_to(lst, value):

    """
    This function ...
    :param lst:
    :param value:
    :return:
    """

    for item in lst:
        if item == value: return True
    return False

# -----------------------------------------------------------------

def has_zero(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    return has_equal_to(lst, 0)

# -----------------------------------------------------------------

def all_identical_to(lst, value):

    """
    This function ...
    :param lst:
    :param value:
    :return:
    """

    for item in lst:
        if item is not value: return False
    return True

# -----------------------------------------------------------------

def all_none(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    return all_identical_to(lst, None)

# -----------------------------------------------------------------

def all_equal(lst, ignore_none=False, ignore=None):

    """
    This function ...
    :param lst:
    :param ignore_none:
    :param ignore:
    :return:
    """

    if len(lst) == 0: raise ValueError("Cannot use empty list")

    first = lst[0]

    if first is None and ignore_none:
        try: first = find_first_not_none(lst, ignore=ignore)
        except: raise ValueError("Cannot use empty list (except for Nones)") #return True # ALL NONE, SO ALL EQUAL

    # ELIF because ignore is also passed to find_first_not_none
    elif ignore is not None and first == ignore:
        try:
            if ignore_none: first = find_first_not_none(lst, ignore=ignore)
            else: first = find_first(lst, ignore=ignore)
        except: raise ValueError("Cannot use empty list (except for " + str(ignore) + ")")

    #print(first)

    for index in range(len(lst)):

        # Ignore None?
        if ignore_none and lst[index] is None: continue

        # Ignore other?
        if ignore is not None and lst[index] == ignore: continue

        #print("comparing:")
        #print(first)
        #print(lst[index])

        #if not (lst[index] == first):
        if lst[index] != first:
            #print("HEERE")
            return False

    # All checks passed
    return True

# -----------------------------------------------------------------

def get_all_equal_value(sequence, ignore_none=False, ignore=None, return_none=False):

    """
    This function ...
    :param sequence:
    :param ignore_none:
    :param ignore:
    :param return_none:
    :return:
    """

    if not all_equal(sequence, ignore_none=ignore_none, ignore=ignore):
        if return_none: return None
        else: raise ValueError("Not all equal: " + str(sequence))
    else:
        if ignore_none: return find_first_not_none(sequence, ignore=ignore)
        elif ignore is not None: return find_first(sequence, ignore=ignore)
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

    # Loop over the entries in the list
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

def find_first_not_none(lst, ignore=None):

    """
    This function ...
    :param lst:
    :param ignore:
    :return:
    """

    for item in lst:
        if ignore is not None and item == ignore: continue
        if item is not None: return item

    # Shouldn't get here
    raise ValueError("No not-None values")

# -----------------------------------------------------------------

def find_first(lst, ignore=None):

    """
    This function ...
    :param lst:
    :param ignore:
    :return:
    """

    for item in lst:
        if item != ignore: return item

    # Shouldn't get here
    raise ValueError("No not-" + str(ignore) + " values")

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

def has_not_none(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    for index in range(len(lst)):
        if lst[index] is not None: return True
    return False

# -----------------------------------------------------------------

def has_other(seq, test):

    """
    This function ...
    :param seq:
    :param test:
    :return:
    """

    for item in seq:
        if item not in test: return True
    return False

# -----------------------------------------------------------------

def get_other(seq, test):

    """
    This function ...
    :param seq:
    :param test:
    :return:
    """

    other = []
    for item in seq:
        if item not in test: other.append(item)
    return other

# -----------------------------------------------------------------

def equal_sequences(*sequences):

    """
    This function ...
    :param sequences:
    :return:
    """

    from . import types

    if not equal_sizes(*sequences): return False
    first_sequence = sequences[0]
    for index in range(len(first_sequence)):
        items = [sequence[index] for sequence in sequences]
        if not all_equal(items):
            if all_true([types.is_real_type(item) for item in items]):
                if not all_close(items): return False
            elif all_true([types.is_quantity(item) for item in items]):
                if not all_close(items): return False
            else: return False
    return True

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

def elements_not_in_other(sequence_a, sequence_b, check_existing=False):

    """
    This function ...
    :param sequence_a: 
    :param sequence_b:
    :param check_existing:
    :return: 
    """

    elements = set()

    # Check whether each item in sequence_b exists in sequence_a
    if check_existing:
        for element in sequence_b:
            if element not in sequence_a: raise ValueError("The element '" + str(element) + "' from the second sequence is not in the first sequence (" + str(sequence_a) + ")")

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

def unique_values(sequence, ignore=None, ignore_none=False):

    """
    This function ...
    :param sequence:
    :param ignore:
    :param ignore_none:
    :return:
    """

    # Contains other sequences or dictionaries (non-hashable)
    if contains_sequence(sequence): sequence = replace_sequences_by_tuples(sequence)
    if contains_dictionary(sequence): sequence = replace_dictionaries_by_tuples(sequence)

    result = list(set(sequence))
    if ignore_none: result = removed(result, [None])
    if ignore is not None: return removed(result, ignore)
    else: return result

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

    from . import numbers

    try:

        # Even
        if numbers.is_even(len(lst)): middle = len(lst)/2
        # Odd
        else: middle = (len(lst)-1)/2

        assert numbers.is_integer(middle)
        middle = int(middle)
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

def any_empty(*sequences):

    """
    This function ...
    :param sequences:
    :return:
    """

    return any(is_empty(sequence) for sequence in sequences)

# -----------------------------------------------------------------

def find_index(sequence, element):

    """
    This function ...
    :param sequence:
    :param element:
    :return:
    """

    index = next((i for i, t in enumerate(sequence) if element == t), None)
    return index

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

def find_unique_indices(sequence, elements):

    """
    This function ...
    :param sequence:
    :param elements:
    :return:
    """

    indices = []
    for element in elements: indices.append(find_unique(sequence, element))
    return indices

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

def is_empty(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return len(sequence) == 0

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

def in_one(item, sequences, allow_more=True):

    """
    This function ...
    :param item:
    :param sequences:
    :param allow_more:
    :return:
    """

    count = 0
    for sequence in sequences:
        if item in sequence: count += 1

    if count == 0: return False
    elif count == 1: return True
    elif allow_more: return False
    else: raise ValueError("Count > 1")

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

def pick_contains(sequences, item):

    """
    This function ...
    :param sequences:
    :param item:
    :return:
    """

    for sequence in sequences:
        if item in sequence: return sequence
    return None

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

def difference(sequence_a, sequence_b):

    """
    This function returns the elements that are in sequence a but not in sequence b
    :param sequence_a:
    :param sequence_b:
    :return:
    """

    items = []
    for item in sequence_a:
        if item not in sequence_b: items.append(item)
    return items

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

def is_subset(sequence, other_sequence):

    """
    This function ...
    :param sequence:
    :param other_sequence:
    :return:
    """

    for item in sequence:
        if item not in other_sequence: return False
    return True

# -----------------------------------------------------------------

def subset(sequence, indices):

    """
    This function ...
    :param sequence:
    :param indices:
    :return:
    """

    subset = []
    for index in indices: subset.append(sequence[index])
    return subset

# -----------------------------------------------------------------

def random_subset(sequence, nsamples, avoid_duplication=False, ignore=None):

    """
    This function ...
    :param sequence:
    :param nsamples:
    :param avoid_duplication:
    :param ignore:
    :return:
    """

    # Set probabilities
    probabilities = None
    if ignore is not None:
        indices = find_unique_indices(sequence, ignore)
        probabilities = [1] * len(sequence)
        for index in indices: probabilities[index] = 0

    # Lower the number of samples if the sequence is not long enough to take unique samples
    if avoid_duplication and nsamples > len(sequence): nsamples = len(sequence)

    import numpy as np
    return np.random.choice(sequence, nsamples, replace=not avoid_duplication, p=probabilities)

# -----------------------------------------------------------------

def get_first_values(sequence, nvalues):

    """
    This function ...
    :param sequence:
    :param nvalues:
    :return:
    """

    return sequence[:nvalues]

# -----------------------------------------------------------------

def get_last_values(sequence, nvalues):

    """
    This function ...
    :param sequence:
    :param nvalues:
    :return:
    """

    nitems = len(sequence)
    if nvalues >= nitems: return sequence[:]
    else: return sequence[nitems-nvalues:]

# -----------------------------------------------------------------

def all_except_indices(sequence, indices):

    """
    This function ...
    :param sequence:
    :param indices:
    :return:
    """

    new = []
    for index in sequence:
        if index in indices: continue
        new.append(sequence[index])
    return new

# -----------------------------------------------------------------

def all_except(sequence, ignore):

    """
    This function ...
    :param sequence:
    :param ignore:
    :return:
    """

    new = []
    for item in sequence:
        if item in ignore: continue
        new.append(item)
    return new

# -----------------------------------------------------------------

def all_except_one(sequence, ignore):

    """
    Thisf ucntion ...
    :param sequence:
    :param ignore:
    :return:
    """

    new = []
    for item in sequence:
        if item == ignore: continue
        new.append(item)
    return new

# -----------------------------------------------------------------

def make_selection(sequence, selected, not_selected, nrandom=None, all=False, none=False, indices=None, not_indices=None):

    """
    This function ...
    :param sequence:
    :param selected:
    :param not_selected:
    :param nrandom:
    :param all:
    :param none:
    :param indices:
    :param not_indices:
    :return:
    """

    if selected is not None and none: raise ValueError("Selection is made but 'none' is enabled")
    if indices is not None and none: raise ValueError("Selection is made but 'none' is enabled")
    if all and none: raise ValueError("Cannot enable 'all' and 'none' simultaneously")

    # Return empty list
    if none: return []

    # Initialize selection
    selection = None

    # Check
    if selected is not None and indices is not None: raise ValueError("Cannot specify selection and indices at the same time")
    if indices is not None and nrandom is not None: raise ValueError("Cannot specify selection and choose random")

    # Check
    if selected is not None and nrandom is not None: raise ValueError("Cannot specifiy selection and choose random")

    # Random selection
    if nrandom is not None: selection = random_subset(sequence, nrandom, ignore=not_selected, avoid_duplication=True)

    # Check
    if selected is not None and not_indices is not None: raise ValueError("Cannot specify both selection and not_indices")
    if indices is not None and not_indices is not None: raise ValueError("Cannot specify both indices and not_indices")

    # Set selections from indices
    if indices is not None: selection = subset(sequence, indices)

    # Set selection based on which indices not
    if not_indices is not None: selection = all_except_indices(sequence, not_indices)

    # Check
    if selected is not None and not_selected is not None: raise ValueError("Cannot specify both selection and not_selection")

    # Set selections
    if selected is not None: selection = selected

    # Set selections based on which not
    if not_selected is not None and nrandom is None: selection = all_except(sequence, not_selected)

    # Check
    if selected is not None and all: raise ValueError("Cannot make selection and enable 'all'")

    # Check
    if all and not_selected is not None: raise ValueError("Cannot make not_selection and enable 'all'")

    # All
    if all: selection = sequence

    # Return the selection
    return selection

# -----------------------------------------------------------------

def is_ascending(sequence):

    """
    This funciton ...
    :param sequence:
    :return:
    """

    return sorted(sequence) == sequence

# -----------------------------------------------------------------

def is_descending(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    return sorted(sequence, reverse=True) == sequence

# -----------------------------------------------------------------

# def all_true(sequence):
#
#     """
#     This function ...
#     :param sequence:
#     :return:
#     """
#
#     for item in sequence:
#         if not item: return False
#     return True

# -----------------------------------------------------------------

def all_false(sequence):

    """
    This function ...
    :param sequence:
    :return:
    """

    for item in sequence:
        if item: return False
    return True

# -----------------------------------------------------------------

def noccurences(sequence, value):

    """
    This function ...
    :param sequence:
    :param value:
    :return:
    """

    return sequence.count(value)

# -----------------------------------------------------------------

def nzeros(sequence):

    """
    Thisn function ...
    :param sequence:
    :return:
    """

    return noccurences(sequence, 0)

# -----------------------------------------------------------------

def is_unique(sequence, value):

    """
    This function ...
    :param sequence:
    :param value:
    :return:
    """

    if value not in sequence: raise ValueError("The value '" + str(value) + "' is not in the sequence '" + str(sequence) + "'")
    return noccurences(sequence, value) == 1

# -----------------------------------------------------------------

def has_any(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    return len(lst) > 0

# -----------------------------------------------------------------

def contains_any(lst, other):

    """
    This function ...
    :param lst:
    :param other:
    :return:
    """

    for item in other:
        if item in lst: return True
    return False

# -----------------------------------------------------------------

def contains_all(lst, other):

    """
    This function ...
    :param lst:
    :param other:
    :return:
    """

    for item in other:
        if item not in lst: return False
    return True

# -----------------------------------------------------------------

def removed(sequence, remove):

    """
    Thisf unction ...
    :param sequence:
    :param remove:
    :return:
    """

    #print("sequence:", sequence)
    #print("remove:", remove)

    new = []

    for item in sequence:

        if item in remove: continue
        new.append(item)

    return new

# -----------------------------------------------------------------

def removed_item(sequence, remove):

    """
    This function ...
    :param sequence:
    :param remove:
    :return:
    """

    new = []
    for item in sequence:
        if item == remove: continue
        new.append(item)
    return new

# -----------------------------------------------------------------

def remove_indices(sequence, indices):

    """
    This function ...
    :param sequence:
    :param indices:
    :return:
    """

    for index in sorted(indices, reverse=True):
        del sequence[index]

# -----------------------------------------------------------------

def find_first_in_both(seq_a, seq_b):

    """
    This function ...
    :param seq_a:
    :param seq_b:
    :return:
    """

    for item in seq_a:
        if item in seq_b: return item
    return None

# -----------------------------------------------------------------

def find_first_not_in(seq_a, seq_b):

    """
    This function ...
    :param seq_a:
    :param seq_b:
    :return:
    """

    for item in seq_a:
        if item not in seq_b: return item
    return None

# -----------------------------------------------------------------

def sort_with_first_last(sequence, first=None, last=None):

    """
    This function ...
    :param sequence:
    :param first:
    :param last:
    :return:
    """

    if first is not None and not is_sequence(first): first = [first]
    if last is not None and not is_sequence(last): last = [last]

    firsts = []
    between = []
    lasts = []

    used_indices = []

    if first is not None:
        for item in first:
            if item in sequence:
                index = sequence.index(item)
                used_indices.append(index)
                firsts.append(item)

    if last is not None:
        for item in last:
            if item in sequence:
                index = sequence.index(item)
                used_indices.append(index)
                lasts.append(item)

    for index in range(len(sequence)):
        if index in used_indices: continue
        between.append(sequence[index])

    # Return the combined new list
    return firsts + between + lasts

# -----------------------------------------------------------------

def find_endswith(sequence, endswith):

    """
    This function ...
    :param sequence:
    :param endswith:
    :return:
    """

    result = []
    for item in sequence:
        if item.endswith(endswith): result.append(item)
    return result

# -----------------------------------------------------------------

def find_unique_endswith(sequence, endswith):

    """
    This function ...
    :param sequence:
    :param endswith:
    :return:
    """

    result = find_endswith(sequence, endswith)
    if len(result) == 0: raise ValueError("None found")
    elif len(result) > 1: raise ValueError("Not unique: '" + endswith + "' in " + str(sequence))
    else: return result[0]

# -----------------------------------------------------------------

def clip_above(sequence, above):

    """
    Thisfunction ...
    :param sequence:
    :param above:
    :return:
    """

    new = []
    for item in sequence:
        if item > above: continue
        new.append(item)
    return new

# -----------------------------------------------------------------

def clip_below(sequence, below):

    """
    This function ...
    :param sequence:
    :param below:
    :return:
    """

    new = []
    for item in sequence:
        if item < below: continue
        new.append(item)
    return new

# -----------------------------------------------------------------

def all_strings(sequence, ignore=None, ignore_none=False, ignore_instance=None):

    """
    This function ...
    :param sequence:
    :param ignore:
    :param ignore_none:
    :param ignore_instance:
    :return:
    """

    from .types import is_string_type
    for item in sequence:
        if ignore is not None and item == ignore: continue
        if ignore_none and item is None: continue
        if ignore_instance is not None and isinstance(item, ignore_instance): continue
        if not is_string_type(item): return False
    return True

# -----------------------------------------------------------------
