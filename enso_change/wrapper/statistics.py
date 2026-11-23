# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions to compute elementary statistics
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from typing import Hashable, Literal, Union

# numpy
from numpy import array as numpy__array
from numpy import broadcast_to as numpy__broadcast_to
from numpy import clip as numpy__clip
from numpy import mean as numpy__mean
from numpy import median as numpy__median
from numpy import ndarray as numpy__ndarray
from numpy import ndindex as numpy__ndindex
from numpy import std as numpy__std
from numpy import var as numpy__var
from numpy import zeros as numpy__zeros
# scipy
from scipy.stats import anderson as scipy__stats__anderson
from scipy.stats import kurtosis as scipy__stats__kurtosis
from scipy.stats import linregress as scipy__stats__linregress
from scipy.stats import normaltest as scipy__stats__normaltest
from scipy.stats import pearsonr as scipy__stats__pearsonr
from scipy.stats import scoreatpercentile as scipy__stats__scoreatpercentile
from scipy.stats import shapiro as scipy__stats__shapiro
from scipy.stats import skew as scipy__stats__skew
from scipy.stats import t as scipy__stats__t

# local package
from enso_change.tools.default_tools import none_to_default
from enso_change.wrapper.numpy_tools import do_list_to_array, transpose_axis_to_last_dimension
from enso_change.wrapper import xarray_base
from enso_change.wrapper.xarray_base import array_wrapper
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def _normality_test(da, method="shapiro", **kwargs):
    if method == "agostino":
        _, p_value = scipy__stats__normaltest(da)
    elif method == "anderson":
        _, p_value = scipy__stats__anderson(da)
    else:
        _, p_value = scipy__stats__shapiro(da)
    return p_value


def correlation_pearson(
        arr_i1: Union[array_wrapper, numpy__ndarray],
        arr_i2: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        data_var_i1: str = None,
        data_var_i2: str = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the Pearson correlation between the two given arrays using my own code (neither numpy nor scipy allow to
    compute the Pearson correlation along a given axis)

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute Pearson correlation; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute Pearson correlation of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute Pearson correlation; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute Pearson correlation of the flattened array)
    :param data_var_i1: str, optional
        Data variable in arr_i1; e.g., data_var_i1 = "ts".
        Default is None
    :param data_var_i2: str, optional
        Data variable in arr_i2; e.g., data_var_i2 = "ts".
        Default is None
    **kwargs: dict
        Additional keyword arguments passed on to xarray.
        The following options are available (see xarray_base.correlation for details):
            'weights': array of weights.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        The Pearson correlation coefficient(s)
    """
    # compute correlation between the two nD arrays
    if isinstance(arr_i1, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
        # arguments for corr
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["skipna", "weights"])
        # compute xarray's corr
        arr_o = xarray_base.correlation(
            arr_i1, arr_i2, dim=dim, data_var_a=data_var_i1, data_var_b=data_var_i2, **xarray_kwargs)
    else:
        # input to numpy array if needed
        arr_t1, arr_t2 = do_list_to_array(arr_i1), do_list_to_array(arr_i2)
        # transpose or flatten arrays
        if axis is None:
            arr_t1, arr_t2 = arr_t1.flatten(), arr_t2.flatten()
        elif arr_t1.shape == arr_t2.shape:
            # transpose arrays (put axis 'axis' in first position)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t1.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t1.ndim) if k != axis_p])
            arr_t1, arr_t2 = arr_t1.transpose(ordered_axes_tuple), arr_t2.transpose(ordered_axes_tuple)
        else:
            # transpose arrays (put axis 'axis' in first position)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t1.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t1.ndim) if k != axis_p])
            arr_t1 = arr_t1.transpose(ordered_axes_tuple)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t2.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t2.ndim) if k != axis_p])
            arr_t2 = arr_t2.transpose(ordered_axes_tuple)
            # broadcast array if one has fewer dimensions than the other
            if arr_t1.ndim < arr_t2.ndim:
                arr_t1 = arr_t1.reshape(arr_t1.shape + (1,) * (arr_t2.ndim - arr_t1.ndim))
            if arr_t2.ndim < arr_t1.ndim:
                arr_t2 = arr_t2.reshape(arr_t2.shape + (1,) * (arr_t1.ndim - arr_t2.ndim))
        # compute r and p values
        arr_o, _ = scipy__stats__pearsonr(arr_t1, arr_t2, axis=0)
    arr_o = numpy__clip(arr_o, -1, 1)
    return arr_o


def correlation_pearson_pvalue(
        arr_i1: Union[array_wrapper, numpy__ndarray],
        arr_i2: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        data_var_i1: str = None,
        data_var_i2: str = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the 2-tailed p-value associated with Pearson’s correlation between the two given arrays using my own code
    (neither numpy nor scipy allow to compute the Pearson correlation along a given axis)

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute Pearson correlation; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute Pearson correlation of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute Pearson correlation; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute Pearson correlation of the flattened array)
    :param data_var_i1: str, optional
        Data variable in arr_i1; e.g., data_var_i1 = "ts".
        Default is None
    :param data_var_i2: str, optional
        Data variable in arr_i2; e.g., data_var_i2 = "ts".
        Default is None
    **kwargs: dict
        Additional keyword arguments passed on to xarray.
        The following options are available (see xarray_base.correlation for details):
            'skipna': if True, skip NaNs when computing function.
            'weights': array of weights.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        The 2-tailed p-value of Pearson’s correlation
    """
    # compute correlation between the two nD arrays
    if isinstance(arr_i1, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
        # arguments for corr
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["skipna", "weights"])
        # compute xarray's corr
        arr_o = xarray_base.xskillscore.pearson_r_p_value(arr_i1, arr_i2, dim=dim, **xarray_kwargs)
    else:
        # input to numpy array if needed
        arr_t1, arr_t2 = do_list_to_array(arr_i1), do_list_to_array(arr_i2)
        # transpose or flatten arrays
        if axis is None:
            arr_t1, arr_t2 = arr_t1.flatten(), arr_t2.flatten()
        elif arr_t1.shape == arr_t2.shape:
            # transpose arrays (put axis 'axis' in first position)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t1.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t1.ndim) if k != axis_p])
            arr_t1, arr_t2 = arr_t1.transpose(ordered_axes_tuple), arr_t2.transpose(ordered_axes_tuple)
        else:
            # transpose arrays (put axis 'axis' in first position)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t1.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t1.ndim) if k != axis_p])
            arr_t1 = arr_t1.transpose(ordered_axes_tuple)
            axis_p = copy__deepcopy(axis) if axis >= 0 else arr_t2.ndim + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(arr_t2.ndim) if k != axis_p])
            arr_t2 = arr_t2.transpose(ordered_axes_tuple)
            # broadcast array if one has fewer dimensions than the other
            if arr_t1.ndim < arr_t2.ndim:
                arr_t1 = arr_t1.reshape(arr_t1.shape + (1,) * (arr_t2.ndim - arr_t1.ndim))
            if arr_t2.ndim < arr_t1.ndim:
                arr_t2 = arr_t2.reshape(arr_t2.shape + (1,) * (arr_t1.ndim - arr_t2.ndim))
        # compute r and p values
        _, arr_o = scipy__stats__pearsonr(arr_t1, arr_t2, axis=0)
    return arr_o


def kurtosis(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        bias: bool = True,
        fisher: bool = True,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the kurtosis (Fisher) along the given axis using scipy.
    g2 (G2 if bias is False) as in Wright and Herrington 2011, https://doi.org/10.3758/s13428-010-0044-x

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute kurtosis; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute kurtosis of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute kurtosis; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute kurtosis of the flattened array)
    :param bias: bool, optional
        If False, then the calculations are corrected for statistical bias.
        Default is True
    :param fisher: bool, optional
        If True, Fisher’s definition is used (normal ==> 0.0).
        If False, Pearson’s definition is used (normal ==> 3.0).
        Default is True
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see
        https://docs.xarray.dev/en/stable/generated/xarray.DataArray.reduce.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the kurtosis values.
    """
    # compute kurtosis
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for reduce
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs"])
        # compute scipy.stats.kurtosis using xarray's reduce
        arr_o = arr_i.reduce(scipy__stats__kurtosis, dim=dim, bias=bias, fisher=fisher, **xarray_kwargs)
    else:
        # input to numpy array if needed
        arr_t = do_list_to_array(arr_i)
        arr_o = scipy__stats__kurtosis(arr_t, axis=axis, bias=bias, fisher=fisher)
    return arr_o


def linear_regression(
        arr_i1: Union[list, numpy__ndarray],
        arr_i2: Union[list, numpy__ndarray],
        **kwargs) -> dict[str, float]:
    """
    Compute the linear least-squares regression between the two given arrays using scipy

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like
    **kwargs - Discarded

    Output:
    -------
    :return slope: float
        Slope of the regression line
    :return intercept: float
        Intercept of the regression line
    :return correlation: float
        The Pearson correlation coefficient
    :return p_value: float
        The p-value for a hypothesis test whose null hypothesis is that the slope is zero, using Wald Test with
        t-distribution of the test statistic
    """
    # input to numpy array if needed
    arr_t1, arr_t2 = do_list_to_array(arr_i1), do_list_to_array(arr_i2)
    # compute linear regression
    slope, intercept, correlation, p_value, _ = scipy__stats__linregress(arr_t1, arr_t2)
    return {"i": intercept, "p": p_value, "r": correlation, "s": slope}


def linear_regression_slope(
        arr_i1: Union[array_wrapper, numpy__ndarray],
        arr_i2: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        data_var_i1: str = None,
        data_var_i2: str = None,
        ddof: int = 0,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the slope of the linear least-squares regression between the two given arrays using xarray

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute linear regression slope; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute linear regression slope of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute linear regression slope; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute linear regression slope of the flattened array)
    :param data_var_i1: str, optional
        Data variable in arr_i1; e.g., data_var_i1 = "ts".
        Default is None
    :param data_var_i2: str, optional
        Data variable in arr_i2; e.g., data_var_i2 = "ts".
        Default is None
    :param ddof: int, optional
        Delta Degrees of Freedom: the divisor used in the calculation is N - ddof, where N represents the number of
        elements.
        Default is 0
    **kwargs: dict
        Additional keyword arguments passed on to xarray.
        The following options are available (see xarray_base.covariance for details):
            'weights': array of weights.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        The slope of the linear least-squares regression
    """
    # compute correlation between the two nD arrays
    if isinstance(arr_i1, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
        # arguments for corr
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["weights"])
        # compute xarray's cov
        arr_t1 = xarray_base.covariance(arr_i1, arr_i2, dim=dim, data_var_a=data_var_i1, data_var_b=data_var_i2,
                                        **xarray_kwargs)
        # compute variance
        arr_t2 = arr_i1.var(ddof=ddof, dim=dim, **xarray_kwargs)
        # slope = cov(A, B) / var(A)
        arr_o = arr_t1 / arr_t2
    else:
        # input to numpy array if needed
        arr_t1, arr_t2 = do_list_to_array(arr_i1), do_list_to_array(arr_i2)
        # transpose or flatten arrays
        if axis is None:
            arr_t1, arr_t2 = arr_t1.flatten(), arr_t2.flatten()
        else:
            # transpose arrays (put axis 'axis' in first position)
            axis_p = copy__deepcopy(axis) if axis >= 0 else len(arr_i1.shape) + axis
            ordered_axes_tuple = tuple([axis_p] + [k for k in range(len(arr_t1.shape)) if k != axis_p])
            arr_t1, arr_t2 = arr_t1.transpose(ordered_axes_tuple), arr_t2.transpose(ordered_axes_tuple)
        # covariance
        xm, ym = arr_t1 - arr_t1.mean(axis=0), arr_t2 - arr_t2.mean(axis=0)
        cov_t = xm * ym
        cov_t = cov_t.sum(axis=0) / cov_t.shape[0]
        # variance of x variable
        var_x = numpy__var(arr_t1, axis=0, ddof=ddof)
        # slope of the linear regression
        arr_o = cov_t / var_x
    return arr_o


def mean(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the mean along the given axis using numpy or xarray.

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute mean; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute mean of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute mean; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute mean of the flattened array)
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see
        https://docs.xarray.dev/en/stable/generated/xarray.DataArray.mean.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the mean values.
    """
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for mean
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs", "skipna"])
        # compute xarray's mean
        arr_o = arr_i.mean(dim=dim, **xarray_kwargs)
    else:
        arr_o = numpy__mean(arr_t, axis=axis)
    return arr_o


def median(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the median along the given axis using numpy or xarray.

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute median; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute median of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute median; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute median of the flattened array)
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see
        https://docs.xarray.dev/en/latest/generated/xarray.DataArray.median.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the median values.
    """
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for median
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs"])
        # compute xarray's median
        arr_o = arr_i.median(dim=dim, **xarray_kwargs)
    else:
        arr_o = numpy__median(arr_t, axis=axis)
    return arr_o


def normality_test(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        method: Literal["agostino", "anderson", "shapiro"] = "shapiro",
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # compute xarray's median
        arr_o = xarray_base.apply_ufunc(arr_i, _normality_test, dim=dim, method=method)
    else:
        arr_o = None
    return arr_o


def percentile(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        per: Union[float, list[float], tuple[float]] = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the per th percentile(s) of the data along the specified dimension using scipy or xarray.

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute percentile; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute percentile of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute percentile; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute percentile of the flattened array)
    :param per: float, list[float], tuple[float], optional
        Percentile(s) at which to extract score. Values should be in range [0, 100].
        Default is [0, 100] (full range returned)
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see
        https://docs.xarray.dev/en/latest/generated/xarray.DataArray.quantile.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.
            'method': Interpolation method to use when the desired quantile lies between two data points.
            'skipna': True to skip missing values (as marked by NaN).

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Score at percentile(s).
    """
    per = none_to_default(per, [0, 100])
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    # compute percentile
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for quantile
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs", "method", "skipna"])
        # xarray quantile requires quantiles not percentiles
        quantile = per / 100 if isinstance(per, (float, int)) is True else [k / 100 for k in per]
        # compute xarray's quantile
        arr_o = arr_i.quantile(quantile, dim=dim, **xarray_kwargs)
    else:
        arr_o = scipy__stats__scoreatpercentile(arr_t, per, axis=axis)
    return arr_o


def skewness(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        bias: bool = True,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the skewness (Fisher-Pearson) along the given axis using scipy.
    g1 (G1 if bias is False) as in Wright and Herrington 2011, https://doi.org/10.3758/s13428-010-0044-x

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute skewness; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute skewness of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute skewness; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute skewness of the flattened array)
    :param bias: bool, optional
        If False, then the calculations are corrected for statistical bias.
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see
        https://docs.xarray.dev/en/stable/generated/xarray.DataArray.reduce.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the skewness values.
    """
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    # compute skewness
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for reduce
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs", "nan_policy"])
        # compute scipy.stats.skew using xarray's reduce
        arr_o = arr_i.reduce(scipy__stats__skew, dim=dim, bias=bias, **xarray_kwargs)
    else:
        arr_o = scipy__stats__skew(arr_t, axis=axis, bias=bias)
    return arr_o


def standard_deviation(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        ddof: int = 0,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the standard deviation along the given axis using numpy or xarray.

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute standard_deviation; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute standard_deviation of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute standard_deviation; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute standard_deviation of the flattened array)
    :param ddof: int, optional
        Delta Degrees of Freedom: the divisor used in the calculation is N - ddof, where N represents the number of
        elements.
        Default is 0
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see https://docs.xarray.dev/en/stable/generated/xarray.DataArray.std.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the standard deviation values.
    """
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    # compute standard deviation
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for std
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs", "skipna"])
        # compute xarray's std
        arr_o = arr_i.std(ddof=ddof, dim=dim, **xarray_kwargs)
    else:
        arr_o = numpy__std(arr_t, axis=axis, ddof=ddof)
    return arr_o


def student_confidence_interval(
        arr_i: Union[array_wrapper, numpy__ndarray],
        arr_v: Union[array_wrapper, numpy__ndarray, None] = None,
        axis: int = None,
        dim: Union[Hashable, str] = None,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the confidence interval of the mean along the given axis.

    :param arr_i: array_like
    :param arr_v: array_like
        Array to compute the variance (e.g., piControl used for internal variability)
    :param axis: int, optional
        Axis (numpy) along which to compute confidence interval; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute confidence interval of the flattened array)
    :param dim: Hashable or str, optional
        Name of dimension (xarray) along which to compute confidence interval; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute confidence interval of the flattened array)
    **kwargs: dict
        Additional keyword arguments for xarray and student_interval.
        The following options are available (see student_interval for details):
            'confidence_level': Percentage of the distribution that must be within the confidence interval.
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: xarray.DataArray or numpy.ndarray
        Array containing the confidence interval values.
    """
    # input to numpy array if needed
    arr_t1 = do_list_to_array(arr_i)
    arr_t2 = do_list_to_array(arr_v)
    # set axis and dimension (if axis is given find dimension name and conversely)
    if isinstance(arr_i, array_wrapper):
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
    if isinstance(arr_i, array_wrapper):
        arr_t1 = arr_t1.to_numpy()
    if isinstance(arr_v, array_wrapper):
        arr_t2 = arr_t2.to_numpy()
    # transpose array (put axis 'axis' in last position)
    arr_t1 = transpose_axis_to_last_dimension(arr_t1, axis)
    if isinstance(arr_t2, numpy__ndarray) and len(arr_t2.shape) == len(arr_t1.shape):
        arr_t2 = transpose_axis_to_last_dimension(arr_t2, axis)
    # initialize output matrix
    matrix_shape = arr_t1.shape[:-1]
    arr_o = numpy__zeros(matrix_shape + (2,))
    # iterate over all dimensions but given one
    for index in numpy__ndindex(matrix_shape):
        # 1D array
        a1 = arr_t1[index]
        # is the variance provided? (or an array to compute the variance)
        a2 = copy__deepcopy(arr_t2)
        if isinstance(arr_t2, numpy__ndarray) and len(arr_t2.shape) == len(arr_t1.shape):
            a2 = a2[index]
        arr_o[index] = student_interval(a1, arr_v=a2, **kwargs)
    # recreate array
    if isinstance(arr_i, array_wrapper):
        # numpy.ndarray to xarray.DataArray
        coordinates = {"interval": numpy__array(list(range(2)))}
        axis_added = [len(arr_o.shape)]
        data_var = kwargs["data_var"] if "data_var" in list(kwargs.keys()) else "unknown"
        arr_o = xarray_base.recreate_array(arr_o, arr_i, axis_added=axis_added, coords_added=coordinates,
                                           data_var=data_var, dim_added=["interval"], dim_removed=[dim])
    return arr_o


def student_interval(
        arr_i: Union[array_wrapper, numpy__ndarray],
        arr_v: Union[array_wrapper, numpy__ndarray, None] = None,
        axis: int = None,
        dim: Union[Hashable, str] = None,
        confidence_level: float = 95,
        **kwargs) -> numpy__ndarray:
    """
    Compute the confidence interval of the mean of the given array

    Input:
    ------
    :param arr_i: array_like
        Must be a 1D array
    :param arr_v: array_like
        1D array to compute the variance (e.g., piControl used for internal variability)
    :param axis: int, optional
        Axis (numpy) along which to compute confidence interval; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute confidence interval of the flattened array)
    :param dim: Hashable or str, optional
        Name of dimension (xarray) along which to compute confidence interval; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute confidence interval of the flattened array)
    :param confidence_level: float, optional
        Percentage of the distribution that must be within the confidence interval.
        Default is 95
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see https://docs.xarray.dev/en/stable/generated/xarray.DataArray.mean.html,
        https://docs.xarray.dev/en/stable/generated/xarray.DataArray.var.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: numpy.ndarray
        Lower and upper bounds of the confidence interval.
    """
    # input to numpy array if needed
    arr_t1 = do_list_to_array(arr_i)
    arr_t2 = do_list_to_array(arr_v)
    # set axis and dimension (if axis is given find dimension name and conversely)
    if isinstance(arr_i, array_wrapper):
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
    # array mean
    loc = mean(arr_t1, axis=axis, dim=dim, **kwargs)
    # standard deviation
    if isinstance(arr_t2, (float, int)):
        scale = copy__deepcopy(arr_t2)
    elif isinstance(arr_t2, numpy__ndarray):
        scale = variance(arr_t2, axis=axis, dim=dim, ddof=1, **kwargs)
    else:
        scale = variance(arr_t1, axis=axis, dim=dim, ddof=1, **kwargs)
    # array size
    size = arr_t1.size
    # number of standard deviations needed to obtain given significance_level
    alpha = 0.5 + confidence_level / 200
    zscore = scipy__stats__t.ppf(alpha, size - 1)
    # standard error
    se = scale ** 0.5 / (size - 1) ** 0.5
    # theoretical uncertainty of the sample mean
    return numpy__array([loc - zscore * se, loc + zscore * se])


def variance(
        arr_i: Union[array_wrapper, numpy__ndarray],
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[str], list[Hashable], tuple[Hashable], tuple[str], None] = None,
        ddof: int = 0,
        **kwargs) -> Union[float, array_wrapper, numpy__ndarray]:
    """
    Compute the variance along the given axis using numpy or xarray.

    Input:
    ------
    :param arr_i: array_like
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute variance; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute variance of the flattened array)
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute variance; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute variance of the flattened array)
    :param ddof: int, optional
        Delta Degrees of Freedom: the divisor used in the calculation is N - ddof, where N represents the number of
        elements.
        Default is 0
    **kwargs: dict
        Additional keyword arguments for xarray.
        The following options are available (see https://docs.xarray.dev/en/stable/generated/xarray.DataArray.var.html):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one.

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        Array containing the variance values.
    """
    # input to numpy array if needed
    arr_t = do_list_to_array(arr_i)
    # compute variance
    if isinstance(arr_i, array_wrapper):
        # set axis and dimension (if axis is given find dimension name and conversely)
        axis, dim = xarray_base.get_axis_dim(arr_i, axis=axis, dim=dim)
        # arguments for var
        xarray_kwargs = dict((ii, jj) for ii, jj in kwargs.items() if ii in ["keep_attrs"])
        # compute xarray's var
        arr_o = arr_i.var(ddof=ddof, dim=dim, **xarray_kwargs)
    else:
        arr_o = numpy__var(arr_t, axis=axis, ddof=ddof)
    return arr_o
# ---------------------------------------------------------------------------------------------------------------------#
