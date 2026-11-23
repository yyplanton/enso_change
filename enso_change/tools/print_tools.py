# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Classes and functions for errors, warnings and logs
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from typing import Union

# local functions
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Classes
# ---------------------------------------------------------------------------------------------------------------------#
class BackgroundColors:
    blue = '\033[94m'
    green = '\033[92m'
    orange = '\033[93m'
    red = '\033[91m'
    normal = '\033[0m'
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
def plural_s(list_i: Union[list, tuple]) -> str:
    """
    Return 's' if there are multiple values in the list

    Input:
    ------
    :param list_i: list

    Output:
    -------
    :return: str
        's' if there are multiple values in the list else ''
    """
    return "s" if isinstance(list_i, (list, tuple)) is True and len(list_i) > 1 else ""


def print_fail(stack_i: list, error_i: str, fail_i: bool = True):
    """
    Print error message and, if asked, stop the code

    Input:
    ------
    :param stack_i: list
        Given by inspect.stack()
    :param error_i: str
        Encountered errors
    :param fail_i: bool
        True to stop the code
    """
    if isinstance(error_i, str) and error_i != "":
        tmp = "ERROR: file " + str(stack_i[0][1]) + " ; fct " + str(stack_i[0][3]) + " ; line " + str(stack_i[0][2])
        if fail_i:
            raise ValueError(BackgroundColors.red + str(tmp) + "\n" + str(error_i) + BackgroundColors.normal)
        else:
            print(BackgroundColors.orange + str(tmp) + "\n" + str(error_i) + BackgroundColors.normal)


def unknown_formater(name: str, unknown: str, known: Union[list, tuple]):
    """
    Format error/warning message when a given string is not among the defined (known).

    Input:
    ------
    :param name: str
        Name of the tested value
    :param unknown: str
        Unknown keyword
    :param known: list or tuple
        Listed defined (known) keywords

    Output:
    -------
    :return: str
        Formated error/warning message
    """
    str_o = str().ljust(5) + "unknown " + str(name) + " " + str(unknown)
    str_o += "\n" + str().ljust(5) + str(name) + str(plural_s(known)) + ": " + ", ".join([repr(k) for k in known])
    return str_o
# ---------------------------------------------------------------------------------------------------------------------#
