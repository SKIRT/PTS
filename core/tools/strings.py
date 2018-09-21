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
import difflib
from collections import defaultdict
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

    return string.split(pattern, 1)

# -----------------------------------------------------------------

def replace(string, pattern, new):

    """
    This functino ...
    :param string:
    :param pattern:
    :param new:
    :return:
    """

    return string.replace(pattern, new)

# -----------------------------------------------------------------

def replace_last(string, pattern, new):

    """
    This function ...
    :param string:
    :param pattern:
    :param new:
    :return:
    """

    before, after = split_at_last(string, pattern)
    return before + new + after

# -----------------------------------------------------------------

def replace_first(string, pattern, new):

    """
    This function ...
    :param string:
    :param pattern:
    :param new:
    :return:
    """

    before, after = split_at_first(string, pattern)
    return before + new + after

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

def title(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.title()

# -----------------------------------------------------------------

def capitalize(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string.capitalize()

# -----------------------------------------------------------------

def smart_title(string):

    """
    This function ...
    :param string:
    :return:
    """

    return ' '.join(smart_capitalize(word) for word in string.split())

# -----------------------------------------------------------------

def smart_capitalize(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string[0].upper() + string[1:]

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

def is_title(string):

    """
    This function ...
    :param string:
    :return:
    """

    return title(string) == string

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
        result.append(title(string))

    elif is_uppercase(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(title(string))

    elif is_title(string):

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))

    else:

        result.append(string)
        result.append(lowercase(string))
        result.append(uppercase(string))
        result.append(title(string))

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

nletters = len(ascii_lowercase)

# -----------------------------------------------------------------

def state_to_absolute(state):

    """
    This function ...
    :param state:
    :return:
    """

    dimension = len(state)
    return sum(state[i] * nletters ** (dimension - 1 - i) for i in range(dimension))

# -----------------------------------------------------------------

def absolute_to_state(absolute, dimension):

    """
    This function ...
    :param absolute:
    :param dimension:
    :return:
    """

    # Initialize state
    new_state = []

    remainder = absolute

    for index in range(dimension):

        divider = nletters ** (dimension - 1 - index)
        quotient = remainder // divider
        remainder = remainder % divider
        new_state.append(quotient)

    # Return
    return new_state

# -----------------------------------------------------------------

class alphabet_strings_iterator(object):

    """
    This function ...
    """

    def __init__(self, dimension=3):

        """
        This function ...
        """

        # Initialize the state
        self.state = [0] * dimension

    # -----------------------------------------------------------------

    def set_state(self, state):

        """
        This function ...
        :param state:
        :return:
        """

        if len(state) != self.dimension: raise ValueError("Incorrect dimension")
        self.state = state

    # -----------------------------------------------------------------

    def set_absolute_state(self, absolute_state):

        """
        This function ...
        :param absolute_state:
        :return:
        """

        # Convert absolute state to state
        new_state = absolute_to_state(absolute_state, self.dimension)

        # Set new state
        self.set_state(new_state)

    # -----------------------------------------------------------------

    def reset(self):

        """
        This function ...
        :return:
        """

        # Set the state again
        self.set_state([0] * self.dimension)

    # -----------------------------------------------------------------

    @property
    def dimension(self):

        """
        This function ...
        :return:
        """

        return len(self.state)

    # -----------------------------------------------------------------

    @property
    def absolute_state(self):

        """
        This function ...
        :return:
        """

        return state_to_absolute(self.state)

    # -----------------------------------------------------------------

    def increment_state(self):

        """
        This function ...
        :return: 
        """

        # Loop over the digits
        for index in reversed(range(self.dimension)):

            # Last index cannot be incremented
            if self.state[index] == 25:

                self.state[index] = 0
                continue

            # Increment this digit
            else:

                self.state[index] += 1
                break

        # No digits left
        else: raise StopIteration

    # -----------------------------------------------------------------

    def decrement_state(self):

        """
        This function ...
        :return:
        """

        # Loop over the other digits
        for index in reversed(range(self.dimension-1)):

            # Last index cannot be decremented
            if self.state[index] == 0:

                self.state[index] = 25
                continue

            # Decrement the digit
            else:

                self.state[index] -= 1
                break

        # No digits left
        else: raise StopIteration

    # -----------------------------------------------------------------

    def increment_to(self, string):

        """
        This function ...
        :param string:
        :return:
        """

        # Sequence of strings
        if types.is_string_sequence(string):

            # Get absolute states
            absolute_states = []

            # Loop over the strings
            for s in string:

                # Check dimension
                if len(s) != self.dimension: raise ValueError("Invalid dimension of string: must be " + str(self.dimension))

                # Convert into state
                state = [ascii_lowercase.index(s[i]) for i in range(self.dimension)]

                # Convert into absolute state
                absolute = state_to_absolute(state)

                # Add
                absolute_states.append(absolute)

            # Take maximum absolute state, add one
            absolute_state = max(absolute_states) + 1

        # Single string
        elif types.is_string_type(string):

            # Check dimension
            if len(string) != self.dimension: raise ValueError("Invalid dimension of string: must be " + str(self.dimension))

            # Convert into state
            state = [ascii_lowercase.index(string[i]) for i in range(self.dimension)]

            # Convert into absolute state, add one
            absolute_state = state_to_absolute(state) + 1

        # Invalid
        else: raise ValueError("Invalid argument: must be string or string sequence")

        # Set the new state if it is higher
        if absolute_state > self.absolute_state: self.set_absolute_state(absolute_state)

    # -----------------------------------------------------------------

    def next(self):

        """
        This function ...
        :return: 
        """

        # Increment state
        #self.increment_state()

        # Return string
        string = "".join([ascii_lowercase[index] for index in self.state])

        # Increment
        self.increment_state()

        # Return the state
        return string

# -----------------------------------------------------------------

def split_except_within_single_quotes(text, add_quotes=True):

    """
    This function strips the whitespace from a string, except when it is within quotes
    :param text:
    :param add_quotes:
    :return:
    """

    return split_except_within(text, "'", add_pattern=add_quotes)

# -----------------------------------------------------------------

def split_except_within_double_quotes(text, add_quotes=True):

    """
    This function ...
    :param text:
    :param add_quotes:
    :return:
    """

    return split_except_within(text, '"', add_pattern=add_quotes)

# -----------------------------------------------------------------

def split_except_within(text, pattern, add_pattern=True):

    """
    This function ...
    :param text:
    :param pattern:
    :param add_pattern:
    :return:
    """

    parts = []
    lst = text.split(pattern)
    for i, item in enumerate(lst):
        if i % 2:
            if add_pattern: parts.append(pattern + item + pattern)
            else: parts.append(item)
        else:
            for a in item.split(): parts.append(a)
    return parts

# -----------------------------------------------------------------

def split_except_within_round_brackets(text, add_brackets=True):

    """
    This function ...
    :param text:
    :param add_brackets:
    :return:
    """

    return split_except_within_left_right(text, "(", ")", add_pattern=add_brackets)

# -----------------------------------------------------------------

def split_except_within_curly_brackets(text, add_brackets=True):

    """
    This function ...
    :param text:
    :param add_brackets:
    :return:
    """

    return split_except_within_left_right(text, "{", "}", add_pattern=add_brackets)

# -----------------------------------------------------------------

def split_except_within_square_brackets(text, add_brackets=True):

    """
    This function ...
    :param text:
    :param add_brackets:
    :return:
    """

    return split_except_within_left_right(text, "[", "]", add_pattern=add_brackets)

# -----------------------------------------------------------------

def split_except_within_left_right(text, left, right, add_pattern=True):

    """
    This function ...
    :param text:
    :param left:
    :param right:
    :param add_pattern:
    :return:
    """

    parts = []
    #lst = re.split(left + '|' + right, text)
    #print(lst)
    lst = text.replace(left, right).split(right)
    for i, item in enumerate(lst):
        if i % 2:
            if add_pattern: parts.append(left + item + right)
            else: parts.append(item)
        else:
            for a in item.split(): parts.append(a)
    return parts

# -----------------------------------------------------------------

def split_except_within_round_brackets_and_double_quotes(text, add_brackets=True, add_quotes=True):

    """
    This function ...
    :param text:
    :param add_brackets:
    :param add_quotes:
    :return:
    """

    return split_except_within_left_right_and_double_quotes(text, "(", ")", add_pattern=add_brackets, add_quotes=add_quotes)

# -----------------------------------------------------------------

# DOESN'T WORK WHEN ( is attached to previous word: will be splitted apart e.g. "aa(b b)" will be "aa", "(b b)"
def split_except_within_left_right_and_double_quotes_old(text, left, right, add_pattern=True, add_quotes=True):

    """
    This function ...
    :param text:
    :param left:
    :param right:
    :param add_pattern:
    :param add_quotes:
    :return:
    """

    parts = []

    opened_quotes = False

    lst = text.replace(left, right).split(right)

    for i, item in enumerate(lst):

        if i % 2:
            if add_pattern: part = left + item + right
            else: part = item
            #print(1, part)
            parts.append(part)

        else:
            for a in item.split():
                #print(2, a)
                if opened_quotes:
                    if a.endswith('"'):
                        opened_quotes = False
                        if add_quotes: parts[-1] += " " + a
                        else: parts[-1] += " " + a[:-1]
                    else: parts[-1] += " " + a
                else:
                    if a.startswith('"'):
                        opened_quotes = True
                        if add_quotes: parts.append(a)
                        else: parts.append(a[1:])
                    else: parts.append(a)
    return parts

# -----------------------------------------------------------------

def split_except_within_left_right_and_double_quotes(text, left, right, add_pattern=True, add_quotes=True):

    """
    This function ...
    :param text:
    :param left:
    :param right:
    :param add_pattern:
    :param add_quotes:
    :return:
    """

    # TEST WHETHER LEFT AND RIGHT ARE SINGLE CHARACTERS!

    parts = []

    opened_as_new = False
    new_word = True
    opened_leftright = False
    opened_quotes = False

    # Loop over the letters
    for letter in text:

        if letter == left:

            opened_leftright = True

            # Add left?
            if new_word:
                opened_as_new = True
                if add_pattern: parts.append(left)
                else: parts.append("") # empty string to start
                new_word = False
            else:
                opened_as_new = False
                parts[-1] += left

        elif letter == right:

            if not opened_leftright: raise RuntimeError("Encountered '" + right + "' without '" + left + "'")
            opened_leftright = False

            # Add right?
            if (not opened_as_new) or add_pattern: parts[-1] += right

        elif letter == double_quote:

            # Now close quotes
            if opened_quotes:

                opened_quotes = False

                # Add quote?
                if (not opened_as_new) or add_quotes: parts[-1] += '"'

            # Now open quotes
            else:

                opened_quotes = True

                # Add quote?
                if new_word:
                    opened_as_new = True
                    if add_quotes: parts.append('"')
                    else: parts.append("") # empty string to start
                    new_word = False
                else:
                    opened_as_new = False
                    parts[-1] += '"'

        # Split
        elif letter == " ":

            if opened_quotes or opened_leftright: parts[-1] += letter # add the space when inside quotes or left/right
            else: new_word = True

        # Start new word
        elif new_word:
            parts.append(letter)
            new_word = False

        #
        else: parts[-1] += letter

    # Return
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

def generate_from_two_parts(part_a, part_b, connectors=(" ", "-", ".", "_", ""), also_reverse=False):

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
    :param string:
    :param length:
    :param as_list:
    :return:
    """

    # Split into list
    lst = wrap(string, length)

    if as_list: return lst
    else: return "\n".join(lst)

# -----------------------------------------------------------------

def get_first_startswith(strings, pattern):

    """
    This function ...
    :param strings:
    :param pattern:
    :return:
    """

    for string in strings:
        if string.startswith(pattern): return string
    return None

# -----------------------------------------------------------------

def get_all_startswith(strings, pattern):

    """
    This function ...
    :param strings:
    :param pattern:
    :return:
    """

    result = []
    for string in strings:
        if string.startswith(pattern): result.append(string)
    return result

# -----------------------------------------------------------------

def get_all_not_startswith(strings, pattern):

    """
    This function ...
    :param strings:
    :param pattern:
    :return:
    """

    result = []
    for string in strings:
        if not string.startswith(pattern): result.append(string)
    return result

# -----------------------------------------------------------------

def get_startswith(line, patterns):

    """
    This function ...
    :param line:
    :param patterns:
    :return:
    """

    result = []
    for pattern in patterns:
        if line.startswith(pattern): result.append(pattern)
    return result

# -----------------------------------------------------------------

def get_unique_startswith(line, patterns, return_none=False):

    """
    This function ...
    :param line:
    :param patterns:
    :param return_none:
    :return:
    """

    result = get_startswith(line, patterns)
    if len(result) == 0:
        if return_none: return None
        else: raise ValueError("No match")
    elif len(result) > 1: raise ValueError("Line '" + line + "' starts with multiple patterns: '" + str(result) + "'")
    else: return result[0]

# -----------------------------------------------------------------

def get_endswith(line, patterns):

    """
    This function ...
    :param line:
    :param patterns:
    :return:
    """

    result = []
    for pattern in patterns:
        if line.endswith(pattern): result.append(pattern)
    return result

# -----------------------------------------------------------------

def get_unique_endswith(line, patterns, return_none=False):

    """
    This function ...
    :param line:
    :param patterns:
    :param return_none:
    :return:
    """

    result = get_endswith(line, patterns)
    if len(result) == 0:
        if return_none: return None
        else: raise ValueError("No match")
    elif len(result) > 1: raise ValueError("Line '" + line + "' ends with multiple patterns: '" + str(result) + "'")
    else: return result[0]

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

def any_startswith(lines, pattern):

    """
    This function ...
    :param lines:
    :param pattern:
    :return:
    """

    for line in lines:
        if line.startswith(pattern): return True
    return False

# -----------------------------------------------------------------

def endswith_any(line, patterns):

    """
    This function ...
    :param line:
    :param patterns:
    :return:
    """

    for pattern in patterns:
        if line.endswith(pattern): return True
    return False

# -----------------------------------------------------------------

def any_endswith(lines, pattern):

    """
    This function ...
    :param lines:
    :param pattern:
    :return:
    """

    for line in lines:
        if line.endswith(pattern): return True
    return False

# -----------------------------------------------------------------

def contains_any(line, patterns):

    """
    Thisf unction ...
    :param line:
    :param patterns:
    :return:
    """

    for pattern in patterns:
        if pattern in line: return True
    return False

# -----------------------------------------------------------------

def any_contains(lines, pattern):

    """
    This function ...
    :param lines:
    :param pattern:
    :return:
    """

    for line in lines:
        if pattern in line: return True
    return False

# -----------------------------------------------------------------

def all_contains(lines, pattern):

    """
    This function ...
    :param lines:
    :param pattern:
    :return:
    """

    for line in lines:
        if pattern not in line: return False
    return True

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

def is_quote_character(character):

    """
    Thisnfunction ...
    :param character:
    :return:
    """

    return character == single_quote or character == double_quote

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

def is_in_single_quotes(string):
    if string == "": return False
    return string[0] == "'" and string[-1] == "'"

# -----------------------------------------------------------------

def contains_single_quotes(string):

    """
    Thisf unction ...
    :param string:
    :return:
    """

    return "'" in string

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

def is_in_double_quotes(string):
    if string == "": return False
    return string[0] == '"' and string[-1] == '"'

# -----------------------------------------------------------------

def contains_double_quotes(string):

    """
    This function ...
    :param string:
    :return:
    """

    return '"' in string

# -----------------------------------------------------------------

def is_quoted(string):

    """
    This function ...
    :param string:
    :return:
    """

    return is_single_quoted(string) or is_double_quoted(string)

# -----------------------------------------------------------------

def is_in_quotes(string):
    return is_in_single_quotes(string) or is_in_double_quotes(string)

# -----------------------------------------------------------------

def unquote(string):

    """
    This function ...
    :param string:
    :return:
    """

    if is_in_quotes(string): return string[1:-1]
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

def add_other_quotes(string):

    """
    Thisf unction ...
    :param string:
    :return:
    """

    quote = other_quote_character(string)
    return quote + string + quote

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

def replace_digits_by_words(string):

    """
    This function ...
    :param string:
    :return:
    """

    from .numbers import int2text

    new_string = ""

    for char in string:
        if char.isdigit(): new_string += int2text(int(char))
        else: new_string += char

    return new_string

# -----------------------------------------------------------------

def replace_first_digit_by_word(string):

    """
    This function ...
    :param string:
    :return:
    """

    from .numbers import int2text

    if string[0].isdigit(): return int2text(int(string[0])) + string[1:]
    else: return string

# -----------------------------------------------------------------

def add_whitespace_or_ellipsis(string, length, ellipsis_position="end", whitespace_position="end"):

    """
    This function ...
    :param string:
    :param length:
    :param ellipsis_position:
    :param whitespace_position:
    :return:
    """

    if len(string) < length: return add_whitespace(string, length, position=whitespace_position)
    elif len(string) == length: return string
    else: return add_ellipsis(string, length, position=ellipsis_position)

# -----------------------------------------------------------------

def add_whitespace(string, length, position="end"):

    """
    This function ...
    :param string:
    :param length:
    :param position:
    :return:
    """

    if len(string) > length: raise ValueError("String is longer than " + str(length))
    elif len(string) == length: return string
    else:

        nspaces = length - len(string)

        if position == "end": return string + " " * nspaces
        elif position == "start": return " " * nspaces + string
        else: raise ValueError("Invalid option for 'position'")

# -----------------------------------------------------------------

def to_length(string, length, ellipsis_position="end", whitespace_position="end", fill=" "):

    """
    This function ...
    :param string:
    :param length:
    :param ellipsis_position:
    :param whitespace_position:
    :param fill:
    :return:
    """

    if fill != " ": raise NotImplementedError("Not yet implemented")
    return add_whitespace_or_ellipsis(string, length, ellipsis_position=ellipsis_position, whitespace_position=whitespace_position)

# -----------------------------------------------------------------

def add_ellipsis(string, length, position="end"):

    """
    This function ...
    :param string:
    :param length:
    :param position: 'start', 'center' or 'end'
    :return:
    """

    if len(string) < length: raise ValueError("String is shorter than " + str(length))
    elif len(string) == length: return string
    else:

        # Ellipsis at end
        if position == "end":

            ellipsis = " ..."
            return string[:length-4] + ellipsis

        # Ellipsis at start
        elif position == "start":

            ellipsis = "... "
            return ellipsis + string[length-len(string)+4]

        # Ellipsis at center
        elif position == "center":

            ellipsis = " ... "
            ncharacters = length - 5
            if ncharacters % 2 == 0: nbegin = nend = int(0.5 * ncharacters)
            else:
                nbegin = int(0.5 * ncharacters)
                nend = ncharacters - nbegin
            return string[:nbegin] + ellipsis + string[-nend:]

        # Invalid
        else: raise ValueError("Invalid option for 'position'")

# -----------------------------------------------------------------

def similarity(string_a, string_b):

    """
    This function ...
    :param string_a:
    :param string_b:
    :return:
    """

    return difflib.SequenceMatcher(None, string_a, string_b).ratio()

# -----------------------------------------------------------------

def split_cumulative(string, pattern, include_total=True):

    """
    THis function ...
    :param string:
    :param pattern:
    :param include_total:
    :return:
    """

    parts = string.split(pattern)

    if len(parts) == 1:
        if include_total: return parts
        else: return []

    result = []
    for index in range(len(parts)):

        previous_indices = range(index)
        all_indices = previous_indices + [index]

        cumulative_parts = [parts[i] for i in all_indices]
        joined = pattern.join(cumulative_parts)

        result.append(joined)

    if include_total: return result
    else: return result[:-1]

# -----------------------------------------------------------------

def hash_string(string):

    """
    This function ...
    :param string:
    :return:
    """

    import hashlib
    hash_object = hashlib.md5(string)
    return hash_object.hexdigest()

# -----------------------------------------------------------------

def alphabetize(strings, index=0, splitted_index=None, split_at=" ", include_all=False):

    """
    Thisn function ...
    :param strings:
    :param index:
    :param splitted_index:
    :param split_at:
    :param include_all:
    :return:
    """

    # Initialize dictionary
    result = defaultdict(list)

    # Include all letters?
    if include_all:
        for letter in alphabet: result[letter] = []

    # Loop over the strings
    for string in strings:

        if splitted_index is not None:

            splitted = string.split(split_at)
            original_string = string
            string = splitted[splitted_index]

        else: original_string = string

        # Get the character
        character = string[index]

        # Add string to the dictionary
        result[character].append(original_string)

    # Return the resulting dictionary
    return result

# -----------------------------------------------------------------

def remove_hash_words(string):

    """
    This function ...
    :param string:
    :return:
    """

    return " ".join(filter(lambda x: x[0] != '#', string.split()))

# -----------------------------------------------------------------

def is_integer(string):

    """
    This function ...
    :param string:
    :return:
    """

    try:
        int(string)
        return True
    except ValueError: return False

# -----------------------------------------------------------------

def integer(number, ndigits, fill="0"):

    """
    This function ...
    :param number:
    :param ndigits:
    :param fill:
    :return:
    """

    number = int(number)
    string = str(number)

    if len(string) > ndigits: raise ValueError("Cannot represent the integer number '" + str(number) + "' with only " + str(ndigits) + " digits")

    # Return
    before = fill * (ndigits - len(string))
    return before + string

# -----------------------------------------------------------------

def number(real, ndecimal, ndigits, fill="0"):

    """
    This function ...
    :param real:
    :param ndecimal:
    :param ndigits:
    :param fill:
    :return:
    """

    magnitude = len(str(float(real)).split(".")[0]) - 1
    ndig_magnitude = magnitude + 1
    if ndig_magnitude + ndecimal > ndigits: raise ValueError("Cannot represent the number '" + str(real) + "' with only " + str(ndigits) + " digits with " + str(ndecimal) + " decimal places")

    string = str(real)

    if "." in string: before, after = string.split(".")
    else: before = string, after = ""

    ndec = len(after)
    if ndec < ndecimal: after += fill * (ndecimal - ndec)
    else: after = after[:ndecimal]

    ndig = ndecimal + len(before)
    before = fill * (ndigits - ndig) + before

    # Return
    return before + "." + after

# -----------------------------------------------------------------

def get_shortest(strings):

    """
    This function ...
    :param strings:
    :return:
    """

    return min(strings, key=len)

# -----------------------------------------------------------------

def get_longest(strings):

    """
    This function ...
    :param strings:
    :return:
    """

    return max(strings, key=len)

# -----------------------------------------------------------------

def is_empty(string):

    """
    This function ...
    :param string:
    :return:
    """

    return string == ""

# -----------------------------------------------------------------

def is_all_spaces(string):

    """
    This function ...
    :param text:
    :return:
    """

    return string.isspace()

# -----------------------------------------------------------------

def is_wrapped_by(text, pattern):

    """
    This function ...
    :param text:
    :param pattern:
    :return:
    """

    if is_empty: return False

    # Pattern is a space
    if pattern == " ": return text.startswith(" ") and text.endswith(" ")

    # Other patterns
    else:

        text = text.strip()
        return text.startswith(pattern) and text.endswith(pattern)

# -----------------------------------------------------------------

def is_wrapped_by_brackets(text):

    """
    This function ...
    :param text:
    :return:
    """

    return is_wrapped_by_round_brackets(text) or is_wrapped_by_curly_brackets(text) or is_wrapped_by_squared_brackets(text)

# -----------------------------------------------------------------

def is_wrapped_by_round_brackets(text):

    """
    This function ...
    :param text:
    :return:
    """

    return text.startswith("(") and text.endswith(")")

# -----------------------------------------------------------------

def is_wrapped_by_curly_brackets(text):

    """
    This function ...
    :param text:
    :return:
    """

    return text.startswith("{") and text.endswith("}")

# -----------------------------------------------------------------

def is_wrapped_by_squared_brackets(text):

    """
    Thisn function ...
    :param text:
    :return:
    """

    return text.startswith("[") and text.endswith("]")

# -----------------------------------------------------------------

def longest_common_substring(string1, string2):

    """
    This function ...
    :param string1:
    :param string2:
    :return:
    """

    from difflib import SequenceMatcher
    match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
    return string1[match.a: match.a + match.size]

# -----------------------------------------------------------------

def common_part(*strings, **kwargs):

    """
    Thisf unction ...
    :param strings:
    :param kwargs:
    :return:
    """

    return_none = kwargs.pop("return_none", True)

    if len(strings) == 1: return strings[0]
    else:

        string = longest_common_substring(strings[0], strings[1])
        if len(strings) == 2:
            if string == "" and return_none: return None
            else: return string

        for index in range(2,len(strings)):
            #print(string, strings[index])
            string = longest_common_substring(string, strings[index])
            #print(string)
            if string == "":
                if return_none: return None
                else: return string

        # Return
        return string

# -----------------------------------------------------------------

def get_substrings(text, startswith, endswith, only_shortest=False):

    """
    This function ...
    :param text:
    :param startswith:
    :param endswith:
    :param only_shortest:
    :return:
    """

    substrings = []
    indices = find_substring_indices(text, startswith)
    for index in indices:

        #substring = text[index:].split()[0]
        substring = text[index:]

        if only_shortest: substrings.append(substring.split(endswith)[0] + endswith)
        else:
            js = find_substring_indices(substring, endswith)
            for j in js: substrings.append(substring[:j+len(endswith)])

    # Return
    return substrings

# -----------------------------------------------------------------

def get_substrings_startswith(text, pattern):

    """
    This function returns all substrings starting with a certain pattern, ending with whitespace
    :param text:
    :param pattern:
    :return:
    """

    return [s[:-1] for s in get_substrings(text, pattern, " ", only_shortest=True)] # don't include the space

# -----------------------------------------------------------------

def get_words_startswith(text, pattern):

    """
    Thisfunction returns all words (whitespace before and after) starting with a certain pattern
    :param text:
    :param pattern:
    :return:
    """

    #return re.findall(r'\b' + pattern + '\w+', text)
    return [t for t in text.split() if t.startswith(pattern)]

# -----------------------------------------------------------------
