# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions to handle numpy arrays
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from itertools import combinations as itertools__combinations
from math import comb as math__comb
from random import sample as random__sample
from typing import Any
# numba
from numba import njit as numba__njit
# numpy
from numpy import array as numpy__array
from numpy import column_stack as numpy__column_stack
from numpy import empty as numpy__empty
from numpy import inf as numpy__inf
from numpy import int64 as numpy__int64
from numpy import isnan as numpy__isnan
from numpy import mean as numpy__mean
from numpy import nan as numpy__nan
from numpy import ndarray as numpy__ndarray
from numpy import ones as numpy__ones
from numpy import prod as numpy__prod
from numpy import tile as numpy__tile
from numpy import vstack as numpy__vstack
from numpy import where as numpy__where
from numpy import zeros as numpy__zeros
from numpy.linalg import lstsq as numpy__linalg__lstsq
from numpy.random import default_rng as numpy__random__default_rng
from numpy.random import randint as numpy__random__randint
# scipy
from scipy.stats import gaussian_kde as scipy__stats__gaussian_kde
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def _array_flatten_but_first(arr_i: numpy__ndarray, **kwargs) -> (numpy__ndarray, numpy__ndarray):
    # get information about the shape/size of the input data
    ax1 = arr_i.shape[0]
    ax2 = arr_i.shape[1:]
    channels = numpy__prod(ax2)
    # reshape to two dimensions (ax1, channels) creating the design matrix
    arr_t = arr_i.reshape([ax1, channels])
    # find indices of values that are not missing in one row. All the rows will have missing values in the same places
    # provided the array was centered. If it wasn't then it is possible that some missing values will be missed and the
    # singular value decomposition will produce not a number for everything.
    non_nan_index = numpy__where(numpy__isnan(arr_t[0]) == False)[0]
    # remove missing values from the design matrix.
    arr_o = arr_t[:, non_nan_index]
    return arr_o, non_nan_index


def _array_reshape(arr_original, arr_computed, non_nan_index, new_axis: int = 0, **kwargs) -> numpy__ndarray:
    # get information about the shape/size of the input data
    ax2 = arr_original.shape[1:]
    channels = numpy__prod(ax2)
    # put data in flat array
    new_shape = (new_axis, channels) if isinstance(new_axis, int) is True and new_axis > 0 else (channels,)
    flat_array = numpy__ones(new_shape, dtype=arr_original.dtype) * numpy__nan
    if isinstance(new_axis, int) is True and new_axis > 0:
        flat_array[:, non_nan_index] = arr_computed
    else:
        flat_array[non_nan_index] = arr_computed
    # reshape flat array to original shape (with new the first axis)
    if isinstance(new_axis, int) is True and new_axis > 0:
        arr_o = flat_array.reshape((new_axis,) + ax2)
    else:
        arr_o = flat_array.reshape(ax2)
    return arr_o


@numba__njit
def _binomial_table(population_size: int) -> numpy__ndarray:
    table_o = numpy__zeros((population_size + 1, population_size + 1), dtype=numpy__int64)
    for ii in range(population_size + 1):
        table_o[ii, 0] = 1
        table_o[ii, ii] = 1
        for jj in range(1, ii):
            table_o[ii, jj] = table_o[ii - 1, jj - 1] + table_o[ii - 1, jj]
    return table_o


@numba__njit
def _unrank(
        rank: int,
        population_size: int,
        sample_size: int,
        binomial_table: numpy__ndarray) -> numpy__ndarray:
    combination_o = numpy__empty(sample_size, dtype=numpy__int64)
    xx = 0
    for pos in range(sample_size):
        remaining = sample_size - pos
        while True:
            count = binomial_table[population_size - xx - 1, remaining - 1]
            if rank < count:
                combination_o[pos] = xx
                xx += 1
                break
            rank -= count
            xx += 1
    return combination_o


@numba__njit
def _unrank_many(
        list_ranks: list,
        population_size: int,
        sample_size: int,
        binomial_table: numpy__ndarray) -> numpy__ndarray:
    indices_o = numpy__empty((len(list_ranks), sample_size), dtype=numpy__int64)
    for ii in range(len(list_ranks)):
        indices_o[ii] = _unrank(list_ranks[ii], population_size, sample_size, binomial_table)
    return indices_o


def combination_indices(
        population_size: int,
        nbr_combinations: int,
        sample_size: int,
        **kwargs) -> numpy__ndarray:
    """
    Compute the given statistic on a resampled array

    Inputs:
    -------
    :param population_size: int
        Number of values in the population; e.g., population_size = 100
    :param nbr_combinations: int
        Maximum number of combinations to use; e.g., nbr_combinations = 1000
    :param sample_size: int
        Number of values in each sample; e.g., sample_size = 10

    Output:
    -------
    :return: ndarray
        Indices to select from the population
    """
    # define size of array and indices of the array
    indices = list(range(population_size))
    # compute the number of combinations
    maximum_combinations = math__comb(population_size, sample_size)
    # compute and select combinations
    if maximum_combinations / nbr_combinations > 1e3:
        # there are so many possible combinations the likelihood of a duplicate combinations is small:
        # -> randomly generate combinations
        rng = numpy__random__default_rng()
        idx = numpy__array([rng.choice(population_size, sample_size, replace=False) for _ in range(nbr_combinations)])
    elif maximum_combinations / nbr_combinations > 50:
        # the number of possible combinations is large but duplicate combinations is possible:
        # -> binary-search-like walk through the implicit combinatorial tree to check for uniqueness
        ranks = random__sample(range(maximum_combinations), nbr_combinations)
        binomial_table = _binomial_table(population_size)
        idx = _unrank_many(ranks, population_size, sample_size, binomial_table)
    elif maximum_combinations < 1e7:
        # compute all possible combinations
        all_combinations = numpy__array(list(itertools__combinations(indices, sample_size)))
        if maximum_combinations <= nbr_combinations:
            # duplicate these combinations as many times as required to reach the desired number of combinations
            nn = nbr_combinations // maximum_combinations
            idx = copy__deepcopy(all_combinations) if nn == 0 else numpy__tile(all_combinations, (nn, 1))
            # if needed, randomly select some combinations to reach the desired number of combinations
            nn = idx.shape[0]
            if nn < nbr_combinations:
                # randomly generate nbr_combinations - nn numbers within [0, maximum_combinations - 1]
                to_select = numpy__array(sorted(random__sample(range(maximum_combinations), nbr_combinations - nn)))
                # select corresponding combinations
                selected_comb = all_combinations[to_select]
                # stack them on output array
                idx = numpy__vstack([selected_comb, idx])
        else:
            # randomly generate nbr_combinations numbers within [0, maximum_combinations - 1]
            to_select = numpy__array(sorted(random__sample(range(maximum_combinations), nbr_combinations)))
            # select corresponding combinations
            idx = all_combinations[to_select]
    else:
        # too many possible combinations and too many desired combinations (nbr_combinations):
        # -> random sampling with replacement
        idx = resampling_indices(population_size, nbr_combinations, sample_size)
    return idx


def do_list_to_array(input_value: Any, **kwargs) -> Any:
    """
    Make a numpy.ndarray from a list (if a list is given).

    Input:
    ------
    :param input_value: Any
        Anything

    Output:
    -------
    :return: Any
        If the input value is a list, returns a numpy.ndarray of this list, else returns input unchanged
    """
    if isinstance(input_value, list):
        input_value = numpy__array(input_value)
    return input_value


def kernel_nd(input_array: numpy__ndarray, axis: int, nbr_samples: int, **kwargs) -> numpy__ndarray:
    """
    Compute gaussian kernel on an N dimensions array.

    Input:
    ------
    :param input_array: numpy.ndarray
        numpy.ndarray to transpose
    :param axis: int
        axis number to compute the kernel on
    :param nbr_samples: int
        Number of samples to generate; e.g., nbr_samples = 1000

    Output:
    -------
    :return: numpy.ndarray
        Array, with given axis number now size nbr_samples
    """
    # compute kernel
    if len(input_array.shape) > 1:
        # set given axis to first dimension; input_array = (t, x, y), axis = 1 --> (x, t, y)
        arr_a = transpose_axis_to_first_dimension(input_array, axis)
        # flatten all but first dimensions and remove NaNs for vectorized processing --> (x, t * y) --> (x, n_non_nan)
        arr_b, non_nan_index = _array_flatten_but_first(arr_a)
        # set first axis (the one along which the kernel must be computed) as last dimension --> (n_non_nan, x)
        arr_c = transpose_axis_to_last_dimension(arr_b, 0)
        # compute kernel for all elements of the first dimension --> (n_non_nan, nbr_samples)
        arr_o = numpy__zeros((arr_c.shape[0], nbr_samples))
        for ii, arr in enumerate(arr_c):
            kernel = scipy__stats__gaussian_kde(arr)
            arr_o[ii] = kernel.resample(nbr_samples)[0]
        # set last axis (kernel axis) as first dimension --> (nbr_samples, n_non_nan)
        arr_o = transpose_axis_to_first_dimension(arr_o, 1)
        # reshape data--> (nbr_samples, t, y)
        arr_o = _array_reshape(arr_a, arr_o, non_nan_index, new_axis=nbr_samples)
    else:
        kernel = scipy__stats__gaussian_kde(input_array)
        arr_o = kernel.resample(nbr_samples)[0]
    return arr_o


def optimize_current_and_lagged_coefficients(
        arr_x: numpy__ndarray,
        arr_y: numpy__ndarray,
        max_lag: int = 100,
        **kwargs) -> ((float, float, float, int), float):
    best_error = numpy__inf
    best_params = None
    # Test all possible reference indices
    for lag in range(1, max_lag + 1):
        # Align X(i), X(i-j), and Y(i)
        x_current, x_lagged, y_valid = arr_x[lag:], arr_x[:-lag], arr_y[lag:]
        # Regression matrix
        matrix = numpy__column_stack([x_current, x_current - x_lagged, numpy__ones(len(x_current))])
        # Least squares optimization
        coeffs, _, _, _ = numpy__linalg__lstsq(matrix, y_valid, rcond=None)
        a, b, c = coeffs
        # Predictions
        y_pred = matrix @ coeffs
        # Mean squared error
        error = numpy__mean((y_valid - y_pred) ** 2)
        if error < best_error:
            best_error = float(error)
            best_params = (float(a), float(b), float(c), lag)
    return best_params, best_error


def optimize_lagged_coefficients(
        arr_x: numpy__ndarray,
        arr_y: numpy__ndarray,
        max_lag: int = 100,
        **kwargs) -> ((float, float, float, int), float):
    best_error = numpy__inf
    best_params = None
    # Test all possible reference indices
    for lag in range(0, max_lag + 1):
        # Align X(i-j) and Y(i)
        if lag == 0:
            x_lagged, y_valid = copy__deepcopy(arr_x), copy__deepcopy(arr_y)
        else:
            x_lagged, y_valid = arr_x[:-lag], arr_y[lag:]
        # Regression matrix
        matrix = numpy__column_stack([x_lagged, numpy__ones(len(x_lagged))])
        # Least squares optimization
        coeffs, _, _, _ = numpy__linalg__lstsq(matrix, y_valid, rcond=None)
        a, b = coeffs
        # Predictions
        y_pred = matrix @ coeffs
        # Mean squared error
        error = numpy__mean((y_valid - y_pred) ** 2)
        if error < best_error:
            best_error = float(error)
            best_params = (float(a), float(b), lag)
    return best_params, best_error


def resampling_indices(
        population_size: int,
        nbr_pseudo_samples: int,
        sample_size: int,
        **kwargs) -> numpy__ndarray:
    """
    Create array of indices for resampling

    Input:
    ------
    :param population_size: int
        Number of values in the population; e.g., population_size = 100
    :param nbr_pseudo_samples: int
        Number of samples to generate; e.g., nbr_pseudo_samples = 1000
        If method = 'combination', this is the maximum number of combinations to use
    :param sample_size: int
        Number of values in each sample; e.g., sample_size = 10

    Output:
    -------
    :return: ndarray
        Indices to select from the population
    """
    # create resamples
    return numpy__random__randint(0, population_size, (nbr_pseudo_samples, sample_size))


def transpose_axis_to_first_dimension(input_array: numpy__ndarray, axis: int, **kwargs) -> numpy__ndarray:
    """
    Transpose given axis to the first dimension.
    E.g., input_array.shape = (2, 3, 4) and axis = 1, output_array.shape = (3, 2, 4)

    Input:
    ------
    :param input_array: numpy.ndarray
        numpy.ndarray to transpose
    :param axis: int
        axis number to transpose at the first dimension

    Output:
    -------
    :return: numpy.ndarray
        Transposed array, with given axis number now at the first dimension
    """
    # transpose arrays (put axis 'axis' in first position)
    axis_p = copy__deepcopy(axis) if axis >= 0 else len(input_array.shape) + axis
    ordered_axes_tuple = tuple([axis_p] + [k for k in range(len(input_array.shape)) if k != axis_p])
    return input_array.transpose(ordered_axes_tuple)


def transpose_axis_to_last_dimension(input_array: numpy__ndarray, axis: int, **kwargs) -> numpy__ndarray:
    """
    Transpose given axis to the last dimension.
    E.g., input_array.shape = (2, 3, 4) and axis = 1, output_array.shape = (2, 4, 3)

    Input:
    ------
    :param input_array: numpy.ndarray
        numpy.ndarray to transpose
    :param axis: int
        axis number to transpose at the last dimension

    Output:
    -------
    :return: numpy.ndarray
        Transposed array, with given axis number now at the last dimension
    """
    # transpose arrays (put axis 'axis' in first position)
    axis_p = copy__deepcopy(axis) if axis >= 0 else len(input_array.shape) + axis
    ordered_axes_tuple = tuple([k for k in range(len(input_array.shape)) if k != axis_p] + [axis_p])
    return input_array.transpose(ordered_axes_tuple)
# ---------------------------------------------------------------------------------------------------------------------#
