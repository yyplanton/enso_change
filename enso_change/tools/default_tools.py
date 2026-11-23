# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Setting to default value tools
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from typing import Any

# local functions
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def none_to_default(input_value: Any, default_value) -> Any:
    """
    Set input_value to given default_value if input_value is None.

    Input:
    ------
    :param input_value: Any
        Input variable, can be anything
    :param default_value: Any
        Default value for input_value, can be anything

    Output:
    -------
    :return: Any
        If input_value is None, default_value, else input_value
    """
    if input_value is None:
        input_value = copy__deepcopy(default_value)
    return input_value


def set_default_str(input_value: str, defined_values: list[str], optional_default: str, **kwargs) -> str:
    """
    Check if input_value is defined in defined_values, if not set input_value to optional_default, if optional_default
    is not in defined_values, set input_value to the first value from defined_values

    Input:
    ------
    :param input_value: str
        Input variable, instance of string; e.g., input_value = "a"
    :param defined_values: list[str]
        List of defined values, must be a list of instance of string; e.g., defined_values = ["a", "b"]
    :param optional_default: str
        Default variable, instance of string, to set to input_value if it is not a valid value
        e.g., optional_default = "b"

    Output:
    -------
    :return: str
        Valid value (as defined in defined_values)
    """
    if input_value not in defined_values and optional_default in defined_values:
        input_value = copy__deepcopy(optional_default)
    elif input_value not in defined_values and optional_default not in defined_values:
        input_value = copy__deepcopy(defined_values[0])
    return input_value


def set_instance(input_value: Any, test_type: type, test_bool: bool, default_value: Any) -> Any:
    """
    Perform an instance test on input_value and keep or modify (to default_value) depending on the result of the test

    Input:
    ------
    :param input_value: Any
        Input variable, can be anything
    :param test_type: type
        Built-in type for the instance test; e.g., test_type = str
    :param test_bool: bool
        True of False, condition of the test to modify the value of input_value; e.g., test_bool = True
    :param default_value: Any
        Default variable, to set to input_value depending on the instance test

    Output:
    -------
    :return: Any
        Input or default variable depending on the instance test
    """
    if isinstance(input_value, test_type) is test_bool:
        input_value = copy__deepcopy(default_value)
    return input_value
# ---------------------------------------------------------------------------------------------------------------------#
