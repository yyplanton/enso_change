# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions to check input values
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from typing import Any, Hashable

# local functions
from enso_change.tools.base_tools import plural_s
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def _types_to_string(type_or_types) -> str:
    """
    Join type names if multiple are given

    Input:
    ------
    :param type_or_types: tuple or type

    Output:
    -------
    :return: str
        Type name(s)
    """
    if isinstance(type_or_types, tuple):
        type_to_print = ", ".join(repr(k.__name__) for k in sorted(type_or_types, key=repr))
    else:
        type_to_print = type_or_types.__name__
    return type_to_print


def check_is_dim(dim: Any) -> bool:
    """
    Do a test on dimension type to decide if it is probably one

    Input:
    ------
    :param dim: Any
        Object to test, usually a string

    Output:
    -------
    :return: bool
    """
    return dim is not None and isinstance(dim, (Hashable, str)) is True


def check_list(input_value, input_name: str, list_of_defined_values, error_list: list[str]):
    """
    Check if given value is a defined value

    Input:
    ------
    :param input_value: anything
    :param input_name: str
        Name of the input value
    :param list_of_defined_values: list
        Names of defined value
    :param error_list: list
        Descriptions of errors
    """
    if input_value not in list_of_defined_values:
        error_list.append("unknown %s: %s" % (repr(input_name), repr(input_value)))
        error_list.append(str().ljust(5) + "known %s%s: %s" % (
            repr(input_name), plural_s(list_of_defined_values),
            ", ".join(repr(k) for k in sorted(list_of_defined_values, key=repr))))


def check_type(input_value, input_name: str, type_or_types, error_list: list[str]):
    """
    Check if given value has the right type

    Input:
    ------
    :param input_value: anything
    :param input_name: str
        Name of the input value
    :param type_or_types: tuple or type
    :param error_list: list
        Descriptions of errors
    """
    if not isinstance(input_value, type_or_types):
        # input is not the right type
        error_list.append("%s type error" % repr(input_name))
        error_list.append(str().ljust(5) + "%s is instance of %s, should be instance of (%s)" % (
            repr(input_name), type(input_value), _types_to_string(type_or_types)))
# ---------------------------------------------------------------------------------------------------------------------#
