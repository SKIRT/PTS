#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.strings Provides functions for dealing with strings.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import re
from textwrap import wrap
from string import ascii_lowercase
import itertools

# Import the relevant PTS classes and modules
from .sequences import interleave, find_indices, find_closest_index
from . import types

# -----------------------------------------------------------------

alphabet = list(ascii_lowercase)

# -----------------------------------------------------------------

superscripts = u"⁰¹²³⁴⁵⁶⁷⁸⁹"
subscripts = u"₀₁₂₃₄₅₆₇₈₉"

# -----------------------------------------------------------------

multiplication = "\u00D7".decode("utf8")

# -----------------------------------------------------------------

def split_at_last(string, pattern):

    """
    This function ...
    :param string:
    :param pattern:
    :return:
    """

    return string.rsplit(pattern, 1)

# -----------------------------------------------------------------

def split_at_first(string, pattern):

    """
    This function ...
    :param string:
    :param pattern:
    :return:
    """

    return string.split('mango', 1)

# -----------------------------------------------------------------

def noccurences(string, character):

    """
    This function ...
    :param string:
    :param character:
    :return:
    """

    return string.count(character)

# -----------------------------------------------------------------

def superscript(index):

    """
    This function ...
    :param index: 
    :return: 
    """

    return superscripts[index]

# -----------------------------------------------------------------

def subscript(index):

    """
    This function ...
    :param index: 
    :return: 
    """

    return subscripts[index]

# -----------------------------------------------------------------

def lowercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.lower()

# -----------------------------------------------------------------

def uppercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.upper()

# -----------------------------------------------------------------

def capitalize(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.lower().title()

# -----------------------------------------------------------------

def is_lowercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.islower()

# -----------------------------------------------------------------

def is_uppercase(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.isupper()

# -----------------------------------------------------------------

def is_capitalized(string):

    """
    This function ...
    :param string:
    :return:
    """

    return capitalize(string) == string

# -----------------------------------------------------------------

def case_combinations(string, also_one_letter=True):

    """
    This function ...
    :param string:
    :param also_one_letter:
    :return:
    """

    if len(string) == 1 and not also_one_letter: return [string]

    result = []

    if is_lowercase(string):

        result.append(string)
        result.append(uppercase(string))
        result.append(capitalize(string))

    elif is_uppercase(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(capitalize(string))

    elif is_capitalized(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))

    else:

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))
        result.append(capitalize(string))

    # Return the resulting list
    return result

# -----------------------------------------------------------------

def case_combinations_list(string_list, also_one_letter=True):

    """
    This function ...
    :param string_list:
    :param also_one_letter:
    :return:
    """

    result = []
    for string in string_list: result += case_combinations(string, also_one_letter)
    return result

# -----------------------------------------------------------------

def quantity_combinations(quantity):

    """
    This function ...
    :param quantity:
    :return:
    """

    result = []

    for value_string in number_representations(quantity.value):

        result.append(value_string + " " + str(quantity.unit))
        result.append(value_string + str(quantity.unit))
        result.append(value_string)

        if quantity.unit == "micron":

            result.append(value_string + "mu")
            result.append(value_string + "um")
            result.append(value_string + " mu")
            result.append(value_string + " um")

    return result

# -----------------------------------------------------------------

def iterate_alphabet():

    """
    This function ...
    :return:
    """

    for letter in alphabet: yield letter

# -----------------------------------------------------------------

def iterate_alphabet_strings():

    """
    This function ...
    :return: 
    """

    for repeat in itertools.count():
        for i in itertools.product(ascii_lowercase, repeat=repeat): yield ''.join(i)

# -----------------------------------------------------------------

class alphabet_strings_iterator(object):

    """
    This function ...
    """

    def __init__(self, dimension=3):

        """
        This function ...
        """

        self.state = [0,0,-1]

    # -----------------------------------------------------------------

    def increment_state(self):

        """
        This function ...
        :return: 
        """

        if self.state[-1] == 25:

            self.state[-1] = 0
            for index in reversed(range(len(self.state)-1)):

                if self.state[index] == 25:
                    self.state[index] = 0
                    continue
                else:
                    self.state[index] += 1
                    break

            else: raise StopIteration

        else: self.state[-1] += 1

    # -----------------------------------------------------------------

    def next(self):

        """
        This function ...
        :return: 
        """

        # Increment state
        self.increment_state()

        # Return string
        return "".join([ascii_lowercase[index] for index in self.state])

# -----------------------------------------------------------------

def split_except_within_single_quotes(text):

    """
    This function strips the whitespace from a string, except when it is within quotes
    :param text:f
    :return:
    """

    parts = []

    lst = text.split("'")

    for i, item in enumerate(lst):

        if i % 2: parts.append("'" + item + "'")
        else:

            for a in item.split(): parts.append(a)

    return parts

# -----------------------------------------------------------------

def stripwhite_around(text, around):

    """
    This function ...
    :param text:
    :param around:
    :return:
    """

    return text.replace(" " + around, around).replace(around + " ", around)

# -----------------------------------------------------------------

def stripwhite_except_quotes(text):

    """
    This function strips the whitespace from a string, except when it is within quotes
    :param text:
    :return:
    """

    lst = text.split('"')
    for i, item in enumerate(lst):
        if not i % 2:
            lst[i] = re.sub("\s+", "", item)
    return '"'.join(lst)

# -----------------------------------------------------------------

def stripwhite_except_curlybrackets(text):

    """
    This function strips the whitespace from a string, except when it is within curly brackets
    :param text:
    :return:
    """

    if "{" in text and "}" in text:
        lst = []
        for part in text.split("{"): lst += part.split("}")
        for i, item in enumerate(lst):
            if not i % 2:
                lst[i] = re.sub("\s+", "", item)
        return "".join(list(interleave([lst, ["{", "}"]])))
    else: return text.replace(" ", "")

# -----------------------------------------------------------------

def stripwhite_except_singlequotes(text):

    """
    This function strips the whitespace from a string, except when it is within single quotes
    :param text:
    :return:
    """

    lst = text.split("'")
    for i, item in enumerate(lst):
        if not i % 2:
            lst[i] = re.sub("\s+", "", item)
    return "'".join(lst)

# -----------------------------------------------------------------

def num_to_ith(num):

    """
    1 becomes 1st, 2 becomes 2nd, etc.
    """

    value             = str(num)
    last_digit        = value[-1]
    if len(value) > 1 and value[-2] == '1': return value +'th'
    if last_digit == '1': return value + 'st'
    if last_digit == '2': return value + 'nd'
    if last_digit == '3': return value + 'rd'
    return value + 'th'

# -----------------------------------------------------------------

def number_representations(value):

    """
    This function ...
    :param value:
    :return:
    """

    strings = set()

    strings.add(str(value))
    strings.add(repr(value))
    strings.add(str_from_real_or_integer(value))

    return list(strings)

# -----------------------------------------------------------------

def str_from_real_or_integer(value):

    """
    This function ...
    :param value:
    :return:
    """

    if types.is_integer_type(value): return str(value)
    elif types.is_real_type(value):

        if int(value) == value: return str(int(value))
        else: return repr(value)

    else: raise ValueError("Value must be real or integer")

# -----------------------------------------------------------------

def generate_from_two_parts(part_a, part_b, connectors=(" ", "-", ".", "_"), also_reverse=False):

    """
    This function ...
    :param part_a:
    :param part_b:
    :param connectors:
    :param also_reverse:
    :return:
    """

    for connector in connectors: yield part_a + connector + part_b

    if also_reverse:
        for connector in connectors: yield part_b + connector + part_a

# -----------------------------------------------------------------

def find_first_digit(string):

    """
    This function ...
    :param string:
    :return:
    """

    index = re.search("\d", string)
    return index

# -----------------------------------------------------------------

def find_last_digit(string):

    """
    This function ...
    :param string:
    :return:
    """

    index = re.search("\d", string[::-1])
    return len(string) - 1 - index

# -----------------------------------------------------------------

def split_in_lines(string, length=60, as_list=False):

    """
    This function ...
    :param length:
    :param as_list:
    :return:
    """

    # Split into list
    lst = wrap(string, length)

    if as_list: return lst
    else: return "\n".join(lst)

# -----------------------------------------------------------------

def startswith_any(line, patterns):

    """
    This function ...
    :param patterns:
    :return:
    """

    for pattern in patterns:
        if line.startswith(pattern): return True
    return False

# -----------------------------------------------------------------

def remove_escape_characters(string):

    """
    This function ...
    :param string:
    :return:
    """

    # A regular expression object that strips away special unicode characters, used on the remote console output
    ansi_escape = re.compile(r'\x1b[^m]*m')

    splitted = ansi_escape.sub('', string).replace('\x1b[K', '').split("\r\n")

    for i in range(len(splitted)): splitted[i] = splitted[i].replace(" \r", "").replace("\x08", "").replace("\xe2\x80\x98", "'").replace("\xe2\x80\x99", "'")

    #if splitted[-1] == "": the_output = splitted[output_start:-1]
    #else: the_output = splitted[output_start:]

    output = "".join(splitted)
    return output

# -----------------------------------------------------------------

def printed_length(string):

    """
    This function ...
    :param string:
    :return:
    """

    converted = remove_escape_characters(string)
    return len(converted)

# -----------------------------------------------------------------

def find_substring_indices(string, substring, allow_repeat=True):

    """
    This function ...
    :param string:
    :param substring:
    :param allow_repeat:
    :return:
    """

    ncharacters = len(substring)

    indices = []

    # Loop over the string
    for index in range(len(string)):

        endindex = index + ncharacters

        if string[index:endindex] != substring: continue

        if not allow_repeat:

            # Previous is also the same character
            if index > 0 and string[index - 1:endindex - 1] == substring: continue

            # Shifted one to the right is also the same substring
            if endindex < len(string) - 1 and string[index + 1:endindex + 1] == substring: continue

        # Add to count
        #count += 1

        # Add index
        indices.append(index)

    # Return the count
    #return count

    # Return the indices
    return indices

# -----------------------------------------------------------------

def get_noccurences(string, substring, allow_repeat=True):

    """
    This function ...
    :param string: 
    :param substring:
    :param allow_repeat: 
    :return: 
    """

    if allow_repeat: return string.count(substring)
    else:

        #characters = list(string)
        ncharacters = len(substring)

        count = 0

        # Loop over the string
        for index in range(len(string)):

            endindex = index + ncharacters

            if string[index:endindex] != substring: continue

            # Previous is also the same character
            if index > 0 and string[index-1:endindex-1] == substring: continue

            # Shifted one to the right is also the same substring
            if endindex < len(string) - 1 and string[index+1:endindex+1] == substring: continue

            # Add to count
            count += 1

        # Return the count
        return count

# -----------------------------------------------------------------

def split_at_indices(string, indices, ncharacters=1):

    """
    This function ...
    :param string:
    :param indices:
    :param ncharacters:
    :return:
    """

    parts = []

    last = ""
    for index in indices:

        endindex = index + ncharacters
        first = string[index:endindex]
        last = string[endindex:]

        parts.append(first)

    parts.append(last)
    return parts

# -----------------------------------------------------------------

def split(string, substring, allow_repeat=True):

    """
    This function ...
    :param string:
    :param substring:
    :param allow_repeat:
    :return:
    """

    if allow_repeat: return string.split(substring)
    else:

        indices = find_substring_indices(string, substring, allow_repeat=allow_repeat)
        ncharacters = len(substring)
        return split_at_indices(string, indices, ncharacters)

# -----------------------------------------------------------------

def find_delimiter(string, noccurences, delimiters=("-", "_", ".", "__", " ")):

    """
    This function ...
    :param string:
    :param noccurences:
    :param delimiters:
    :return:
    """

    counts = []
    for delimiter in delimiters: counts.append(get_noccurences(string, delimiter, allow_repeat=False))

    # Find the delimiters with the specified number of occurences
    indices = find_indices(counts, noccurences)

    # No delimiter
    if len(indices) == 0: raise ValueError("None of the delimiters are present " + str(noccurences) + " time(s)")

    # Just one: FINE!
    elif len(indices) == 1: return delimiters[indices[0]]

    # More than one
    else:

        # Find the delimiter which splits up the string most evenly
        possible_delimiters = [delimiters[index] for index in indices]

        ratios = []

        for delimiter in possible_delimiters:

            #parts = string.split(delimiter)
            parts = split(string, delimiter, allow_repeat=False)
            assert len(parts) == noccurences + 1

            ratio = float(len(parts[0])) / float(len(parts[1])) # TODO: now, we only compare the lengths of the first and second part

            ratios.append(ratio)

        # Find index of delimiter giving a division ratio closest to one
        the_index = find_closest_index(ratios, 1.)

        # Return the delimiter
        return possible_delimiters[the_index]

# -----------------------------------------------------------------

def find_any_case(string, sequence):

    """
    This function ...
    :param string:
    :param sequence:
    :return:
    """

    for s in sequence:
        if s.lower() == string.lower(): return s
    return None # not found

# -----------------------------------------------------------------

def is_character(string):

    """
    This function ...
    :param string:
    :return:
    """

    return len(string) == 1

# -----------------------------------------------------------------

def add_quotes_if_spaces(string):

    """
    This function ...
    :param string:
    :return:
    """

    if " " in string:
        if '"' in string:
            assert "'" not in string
            return "'" + string + "'"
        elif "'" in string: return '"' + string + '"'
        else: return "'" + string + "'"

    else: return string

# -----------------------------------------------------------------

def add_quotes_if_string_with_spaces(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not types.is_string_type(value): return value
    else: return add_quotes_if_spaces(value)

# -----------------------------------------------------------------

def add_quotes_if_string_with_spaces_else_string(value):

    """
    This function ...
    :param value:
    :return:
    """

    if not types.is_string_type(value): return str(value)
    else: return add_quotes_if_spaces(value)

# -----------------------------------------------------------------

single_quote = "'"
double_quote = '"'

# -----------------------------------------------------------------

def is_single_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """


    if "'" in string:
        if '"' in string: raise ValueError("String (" + string + ") contains both kinds of quote characters")
        else: return True
    elif '"' in string: return False
    else: return False

# -----------------------------------------------------------------

def is_double_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """

    if '"' in string:
        if "'" in string: raise ValueError("String (" + string + ") contains both kinds of quote characters")
        else: return True
    elif "'" in string: return False
    else: return False

# -----------------------------------------------------------------

def is_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """

    return is_single_quoted(string) or is_double_quoted(string)

# -----------------------------------------------------------------

def unquote(string):

    """
    This function ...
    :param string:
    :return:
    """

    if is_quoted(string): return string[1:-1]
    else: return string

# -----------------------------------------------------------------

def other_quote_character(*strings):

    """
    This function ...
    :param string:
    :return:
    """

    is_single = False
    is_double = False

    for string in strings:

        if is_single_quoted(string): is_single = True
        elif is_double_quoted(string): is_double = True
        else: pass

    if is_single:
        if is_double: raise ValueError("No consensus: both single and double quote characters are used")
        else: return double_quote
    elif is_double: return single_quote
    else: return double_quote  #doesn't matter

# -----------------------------------------------------------------

def make_single_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """

    if is_single_quoted(string): return string
    elif is_double_quoted(string): return string.replace('"', "'")
    else: return string

# -----------------------------------------------------------------

def make_double_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """

    if is_single_quoted(string): return string.replace("'", '"')
    elif is_double_quoted(string): return string
    else: return string

# -----------------------------------------------------------------

def replace_from_dict(string, replace_dict):

    """
    This function ...
    :param string:
    :param replace_dict:
    :return:
    """

    for pattern in replace_dict: string = string.replace(pattern, replace_dict[pattern])
    return string

# -----------------------------------------------------------------

def contains_which(string, patterns, case=True):

    """
    This function ...
    :param string:
    :param patterns:
    :param case: case sensitive
    :return:
    """

    for pattern in patterns:
        if case:
            if pattern in string: return pattern
        else:
            if pattern.lower() in string.lower(): return pattern

    raise ValueError("None of the patterns are contained in the string")

# -----------------------------------------------------------------
