# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Basic tools built over xarray
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from typing import Hashable, Literal, Union
# numpy
from numpy import array as numpy__array
from numpy import ndarray as numpy__ndarray
from numpy.random import default_rng as numpy__random__default_rng

# regionmask
import regionmask

# local functions
from enso_change.tools import default_tools
from enso_change.wrapper import numpy_tools
from enso_change.wrapper import statistics
from enso_change.wrapper import xarray_base
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, recreate_array
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def detrend(
        da_i: array_wrapper,
        deg: int = 1,
        dim: Union[Hashable, str, None] = None,
        kwargs_polyfit: dict = None,
        **kwargs) -> array_wrapper:
    """
    Remove the least squares polynomial fit from given xarray.DataArray.

    Input:
    ------
    :param da_i: xarray.DataArray
    :param deg: int, optional
        Degree of the fitting polynomial.
        Default is 1
    :param dim: Hashable or str or None, optional
        Name of dimension(s) (xarray) along which to compute statistic; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array)
    :param kwargs_polyfit: dict, optional
        Arguments for the polyfit function.
        https://docs.xarray.dev/en/latest/generated/xarray.DataArray.polyfit.html
            - 'cov' (bool or "unscaled", default: False)
                Whether to return to the covariance matrix in addition to the coefficients. The matrix is not scaled if
                cov='unscaled'.
            - 'full' (bool, default: False)
                Whether to return the residuals, matrix rank and singular values in addition to the coefficients.
            - 'rcond' (float or None, optional)
                Relative condition number to the fit.
            - 'skipna' (bool or None, optional)
                If True, removes all invalid values before fitting each 1D slices of the array. Default is True if data
                is stored in a dask.array or if there is any invalid values, False otherwise.
            - 'w' (Hashable, array-like or None, optional)
                Weights to apply to the y-coordinate of the sample points. Can be an array-like object or the name of a
                coordinate in the dataset.
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Input object with the least squares polynomial fit removed.
    """
    kwargs_polyfit = default_tools.none_to_default(kwargs_polyfit, {})
    # compute coefficient
    p = da_i.polyfit(dim, deg, **kwargs_polyfit)
    # remove fit
    da_o = da_i.copy(deep=True) - xarray_base.polyval(da_i[dim], p["polyfit_coefficients"], degree_dim="degree")
    return da_o


def create_land_sea_mask(
        da_i: array_wrapper,
        mask_as_boolean: bool = False,
        **kwargs) -> Union[array_wrapper, None]:
    """
    Create a land-sea mask (1: land, 0: sea) for given dataset, using regionmask.
    https://regionmask.readthedocs.io/en/stable/defined_landmask.html

    Input:
    ------
    :param da_i: xarray.DataArray
    :param mask_as_boolean : bool, optional
        Define mask as (1: land, 0: sea) or (True: land, False: sea); e.g., mask_as_boolean = False.
        If True, define cells as boolean.
        If False, define cells as float.
        Default is False
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New DataArray of land-sea mask (1: land, 0: sea) or (True: land, False: sea).
    """
    # get lat-lon dimensions and arrays
    dim_x, dim_y = cf_dim_to_dim(da_i, "X"), cf_dim_to_dim(da_i, "Y")
    da_x, da_y = da_i[dim_x], da_i[dim_y]
    # create a land-sea mask using regionmask
    land_sea_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_110
    # adapt land-sea mask to input dataset or array -> land-sea mask is 0 on land and NaN on ocean
    land_sea_mask = land_sea_mask.mask(da_x, da_y)
    # convert zeros to 1 (1: land) and NaN to 0 (0: sea)
    land_sea_mask = land_sea_mask.where(land_sea_mask.isnull(), other=1)
    land_sea_mask = land_sea_mask.where(land_sea_mask == 1, other=0)
    if isinstance(mask_as_boolean, bool) and mask_as_boolean:
        # convert the land-sea mask to a boolean mask
        land_sea_mask = land_sea_mask.astype("bool")
    return land_sea_mask


def remaining_internal_variability_hs09(
        da_i: array_wrapper,
        degree: int = 1,
        dim: str = "year",
        ddof: int = 1,
        **kwargs) -> tuple[array_wrapper, array_wrapper, array_wrapper]:
    """
    Compute the remaining internal variability following Hawkins and Sutton's (2009;
    https://doi.org/10.1175/2009BAMS2607.1) computation of internal variability.
    HS09 hypothesis is that the internal variability can be defined as the departure between a simulation and a
    polynomial fit (of order N) of this simulation.
    From this I define the 'remaining internal variability' following the same method but using the ensemble mean
    instead of a single member.

    Input:
    ------
    :param da_i: xarray.DataArray
    :param degree: int
        Degree of the fitting polynomial.
        Default is 1
    :param dim: str
        Name of the temporal dimension, or the one to use for the polynomial fit.
        Default is 'year'
    :param ddof: int, optional
        Delta Degrees of Freedom to compute the variance: the divisor used in the calculation is N - ddof, where N
        represents the number of elements.
        Default is 0
    **kwargs - Discarded

    Output:
    -------
    :return polynomial_fit: xarray.DataArray
        Values of the Nth order polynomial fit derived from input array.
    :return residual: xarray.DataArray
        Departure from the Nth order polynomial fit (da_i - polynomial_fit).
    :return riv: xarray.DataArray
        Variance  (da_i - polynomial_fit).
    """
    # average across members if not already done
    if "member" in list(da_i.dims):
        da_i = da_i.mean(dim="member")
    # -- compute remaining internal variability
    # 'forced change' is a Nth order polynomial fit
    coef = da_i.polyfit(dim, degree)
    polynomial_fit = xarray_base.polyval(da_i[dim], coef["polyfit_coefficients"], degree_dim="degree")
    # residual is the departure from the polynomial fit
    residual = (da_i - polynomial_fit)
    # 'remaining internal variability' is the temporal variance of the residual
    riv = residual.var(dim=dim, ddof=ddof)
    return polynomial_fit, residual, riv


def resample(
        arr_i1: array_wrapper,
        arr_i2: Union[array_wrapper, None] = None,
        arr_i3: Union[array_wrapper, None] = None,
        axis: Union[int, None] = 0,
        data_var_i1: str = None,
        data_var_i2: str = None,
        data_var_i3: str = None,
        dim: Union[Hashable, str, None] = None,
        pseudo_sample_size: int = 10000,
        sample_method: Literal["combination", "random"] = "random",
        sample_size: int = None,
        **kwargs) -> Union[array_wrapper, tuple[array_wrapper]]:
    """
    Compute a bootstrap on input numpy.ndarray or xarray.DataArray

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like, optional
        Second array if it needs the same resampling
        Default is None
    :param arr_i3: array_like, optional
        Third array if it needs the same resampling
        Default is None
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute statistic; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array).
    :param data_var_i1: str, optional
        Data variable in arr_i1 if it is instance of xarray.DataArray; e.g., data_var_i1 = "ts".
        Default is None
    :param data_var_i2: str, optional
        Data variable in arr_i2 if it is instance of xarray.DataArray; e.g., data_var_i2 = "ts".
        Default is None
    :param data_var_i3: str, optional
        Data variable in arr_i3 if it is instance of xarray.DataArray; e.g., data_var_i3 = "ts".
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute statistic; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array)
    :param pseudo_sample_size: int, optional
        Number of resamples to randomly generate for the bootstrap; e.g., pseudo_sample_size = 1000
        Default is 10000
    :param sample_method: {"combination", "random"}, optional
        Method to create resamples:
            - "combination": creates 'pseudo_sample_size' combinations (distinct members) of 'sample_size', if possible
                             all combinations are unique but if 'pseudo_sample_size' is larger than the maximum number
                             of combinations, the combinations are repeated to fill 'pseudo_sample_size'
            - "random": creates 'pseudo_sample_size' random selection, with replacement (non distinct members), of
                        'sample_size'
        Default is "random"
    :param sample_size: int, optional
        Number of item in each sample; e.g., sample_size = 10
        Default is None (the sample size is the array size)
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray, xarray.DataArray or None
        Resampled objects, first axis of the result corresponds to the resampling and a resample dimension is added to
        the return array. The other dimensions are the dimensions of the input array.
    """
    # numpy arrays are manipulated in this function
    arr_t1 = arr_i1.to_numpy()
    arr_t2 = None
    if isinstance(arr_i2, array_wrapper):
        arr_t2 = arr_i2.to_numpy()
    arr_t3 = None
    if isinstance(arr_i3, array_wrapper):
        arr_t3 = arr_i3.to_numpy()
    # set axis and dimension (if axis is given find dimension name and conversely)
    axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
    # sample size
    if isinstance(arr_t1, numpy__ndarray) and sample_size is None:
        sample_size = arr_t1.shape[axis]
    # create random indices
    if sample_method == "random":
        idx = numpy_tools.resampling_indices(arr_t1.shape[axis], pseudo_sample_size, sample_size)
    else:
        idx = numpy_tools.combination_indices(arr_t1.shape[axis], pseudo_sample_size, sample_size)
    # transpose arrays (put axis 'axis' in first position)
    arr_t1 = numpy_tools.transpose_axis_to_first_dimension(arr_t1, axis)
    # select random indices
    sample_1 = arr_t1[idx]
    # do the same for the second and thrid arrays
    sample_2 = None
    if isinstance(arr_i2, array_wrapper):
        arr_t2 = numpy_tools.transpose_axis_to_first_dimension(arr_t2, axis)
        sample_2 = arr_t2[idx]
    sample_3 = None
    if isinstance(arr_i3, array_wrapper):
        arr_t3 = numpy_tools.transpose_axis_to_first_dimension(arr_t3, axis)
        sample_3 = arr_t3[idx]
    # recreate DataArray
    dim_removed = [dim] if isinstance(dim, (Hashable, str)) is True else copy__deepcopy(dim)
    coordinates = {
        "resample": numpy__array(list(range(pseudo_sample_size))),
        dim: numpy__array(list(range(sample_size)))}
    sample_1 = recreate_array(sample_1, arr_i1, axis_added=[0, 1], coords_added=coordinates, data_var=data_var_i1,
                              dim_added=["resample", dim], dim_removed=dim_removed)
    if isinstance(arr_i2, array_wrapper):
        sample_2 = recreate_array(sample_2, arr_i2, axis_added=[0, 1], coords_added=coordinates, data_var=data_var_i2,
                                  dim_added=["resample", dim], dim_removed=dim_removed)
    if isinstance(arr_i3, array_wrapper):
        sample_3 = recreate_array(sample_3, arr_i3, axis_added=[0, 1], coords_added=coordinates, data_var=data_var_i3,
                                  dim_added=["resample", dim], dim_removed=dim_removed)
    if isinstance(sample_1, array_wrapper) and not isinstance(sample_2, array_wrapper) and \
            not isinstance(sample_3, array_wrapper):
        output = sample_1
    else:
        output = tuple(k for k in [sample_1, sample_2, sample_3] if isinstance(k, array_wrapper))
    return output


def resample_spaced(
        arr_i1: array_wrapper,
        arr_i2: Union[array_wrapper, None] = None,
        axis: Union[int, None] = 0,
        data_var_i1: str = None,
        data_var_i2: str = None,
        dim: Union[Hashable, str, None] = None,
        pseudo_sample_size: int = 10000,
        sample_size: int = None,
        seed: int = None,
        spacing: int = 1,
        **kwargs) -> Union[array_wrapper, tuple[array_wrapper]]:
    # numpy arrays are manipulated in this function
    arr_t1 = arr_i1.to_numpy()
    arr_t2 = None
    if isinstance(arr_i2, array_wrapper):
        arr_t2 = arr_i2.to_numpy()
    # set axis and dimension (if axis is given find dimension name and conversely)
    axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
    # sample size
    if isinstance(arr_t1, numpy__ndarray) and sample_size is None:
        sample_size = arr_t1.shape[axis]
    # population size
    size = arr_t1.shape[axis]
    # minimum spacing
    spacing = copy__deepcopy(spacing) if size / spacing >= 2 * sample_size else size / (2 * sample_size)
    # random generator
    rng = numpy__random__default_rng(seed)
    # create random indices
    idx = []
    for _ in range(pseudo_sample_size):
        # Random permutation of available indices
        candidates = rng.permutation(size)
        selected = []
        for ii in candidates:
            if all(abs(ii - jj) >= spacing for jj in selected):
                selected.append(ii)
                if len(selected) == sample_size:
                    break
        if len(selected) < sample_size:
            raise ValueError(
                f"Could not select {sample_size} samples with minimum spacing {spacing} from {size} points.")
        idx.append(selected)
    idx = numpy__array(idx)
    # transpose arrays (put axis 'axis' in first position)
    arr_t1 = numpy_tools.transpose_axis_to_first_dimension(arr_t1, axis)
    # select random indices
    sample_1 = arr_t1[idx]
    # do the same for the second array
    sample_2 = None
    if isinstance(arr_i2, array_wrapper):
        arr_t2 = numpy_tools.transpose_axis_to_first_dimension(arr_t2, axis)
        sample_2 = arr_t2[idx]
    # recreate DataArray
    dim_removed = [dim] if isinstance(dim, (Hashable, str)) is True else copy__deepcopy(dim)
    coordinates = {
        "resample": numpy__array(list(range(pseudo_sample_size))),
        dim: numpy__array(list(range(sample_size)))}
    sample_1 = recreate_array(sample_1, arr_i1, axis_added=[0, 1], coords_added=coordinates, data_var=data_var_i1,
                              dim_added=["resample", dim], dim_removed=dim_removed)
    if isinstance(arr_i2, array_wrapper):
        sample_2 = recreate_array(sample_2, arr_i2, axis_added=[0, 1], coords_added=coordinates, data_var=data_var_i2,
                                  dim_added=["resample", dim], dim_removed=dim_removed)
    if isinstance(sample_1, array_wrapper) and not isinstance(sample_2, array_wrapper):
        output = sample_1
    else:
        output = tuple(k for k in [sample_1, sample_2] if isinstance(k, array_wrapper))
    return output


def select_region(
        da_i: array_wrapper,
        lats: list,
        lons: list,
        **kwargs) -> array_wrapper:
    """

    :param da_i: xarray.DataArray
    :param lats: list[float]
        List of latitudes. Either 2 points if the region is rectangular (e.g., lats = [-5, 5]) or vertices of polygonal
        region (e.g., lats = [-10, -5, 5, 10, -10])
        Note: length of lats and lons lists must be equal
    :param lons: list[float]
        List of longitudes. Either 2 points if the region is rectangular (e.g., lons = [180, 280]) or vertices of
        polygonal region (e.g., lats = [180, 280, 280, 180, 180])
        Note: length of lats and lons lists must be equal
     - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Input object with cells masked (NaN) outside the defined region.
    """
    # get lat-lon dimensions and arrays
    dim_x, dim_y = cf_dim_to_dim(da_i, "X"), cf_dim_to_dim(da_i, "Y")
    da_x, da_y = da_i[dim_x], da_i[dim_y]
    # define conditions to mask outside region
    if isinstance(lats, list) and len(lats) == 2 and isinstance(lons, list) and len(lons) == 2 and (
            max(lats) - min(lats) < 180 or max(lons) - min(lons) < 360):
        la1, la2, lo1, lo2 = min(lats), max(lats), min(lons), max(lons)
        # -- rectangular region
        if la2 - la1 < 180 and lo2 - lo1 >= 360:
            # select latitude only
            cond = (la1 <= da_y) & (da_y <= la2)
        elif la2 - la1 >= 180 and lo2 - lo1 < 360:
            # select longitude only
            cond = (lo1 <= da_x) & (da_x <= lo2)
        else:
            cond = (la1 <= da_y) & (da_y <= la2) & (lo1 <= da_x) & (da_x <= lo2)
    else:
        # -- polygonal region
        # create region using regionmask
        region = numpy__array([[lo, la] for la, lo in zip(lats, lons)])
        region = regionmask.Regions([region])
        mask = region.mask(da_x, da_y)
        cond = mask.notnull()
    # mask data outside region
    return da_i.where(cond)


def statistic_resampled(
        arr_i1: array_wrapper,
        arr_i2: Union[array_wrapper, None] = None,
        axis: Union[int, None] = 0,
        data_var_i1: str = None,
        data_var_i2: str = None,
        dim: Union[Hashable, str, None] = None,
        pseudo_sample_size: int = 10000,
        sample_method: Literal["combination", "random", "spaced"] = "random",
        sample_size: int = None,
        **kwargs) -> array_wrapper:
    """
    Compute a bootstrap for the given statistic on input numpy.ndarray or xarray.DataArray

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like, optional
        Needed for two array statistics:
            correlation (statistic = 'cor'),
        Default is None
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute statistic; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array).
    :param data_var_i1: str, optional
        Data variable in arr_i1 if it is instance of xarray.DataArray; e.g., data_var_i1 = "ts".
        Default is None
    :param data_var_i2: str, optional
        Data variable in arr_i2 if it is instance of xarray.DataArray; e.g., data_var_i2 = "ts".
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute statistic; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array)
    :param pseudo_sample_size: int, optional
        Number of resamples to randomly generate for the bootstrap; e.g., pseudo_sample_size = 1000
        Default is 10000
    :param sample_method: {"combination", "random"}, optional
        Method to create resamples:
            - "combination": creates 'pseudo_sample_size' combinations (distinct members) of 'sample_size', if possible
                             all combinations are unique but if 'pseudo_sample_size' is larger than the maximum number
                             of combinations, the combinations are repeated to fill 'pseudo_sample_size'
            - "random": creates 'pseudo_sample_size' random selection, with replacement (non distinct members), of
                        'sample_size'
        Default is "random"
    :param sample_size: int, optional
        Number of item in each sample; e.g., sample_size = 10
        Default is None (the sample size is the array size)
    **kwargs: dict
        The following options are available to select a statistic (see statistic function for details):
            'stat': Nickname of the desired statistic {"cor", "kur", "mea", "med", "per", "ske", "std", "tin", "var"}
        The following options are available for statistics using standard deviation or variance (see specific function
        for details):
            'ddof': divisor used in the calculation is N - ddof, where N represents the number of elements
        The following options are available for statistics using kurtosis or skewness (see specific function for
        details):
            'bias': False to corrected for statistical bias
        The following options are available for percentile (statistic = 'per', see stat_percentile for details):
            'per': Percentile(s) at which to extract score. Values should be in range [0, 100]
        The following options are sometimes available for xarray (see specific function for details):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one
            'method' (stat_percentile): Interpolation method to use when the desired quantile lies between two data
                points
            'weights' (stat_correlation_pearson): array of weights

    :return: xarray.DataArray
        Statistic computed on resampled data, first axis of the result corresponds to the resampling and a resample
        dimension is added to the return array. The other dimensions are the dimensions that remain after the reduction
        of the array.
    """
    # set axis and dimension (if axis is given find dimension name and conversely)
    axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
    # resample
    sample_2 = None
    if sample_method == "spaced" and isinstance(arr_i2, array_wrapper):
        sample_1, sample_2 = resample_spaced(
            arr_i1, arr_i2=arr_i2, axis=axis, data_var_i1=data_var_i1, data_var_i2=data_var_i2, dim=dim,
            pseudo_sample_size=pseudo_sample_size, sample_size=sample_size, **kwargs)
    elif sample_method == "spaced":
        sample_1 = resample_spaced(
            arr_i1, axis=axis, data_var_i1=data_var_i1, data_var_i2=data_var_i2, dim=dim,
            pseudo_sample_size=pseudo_sample_size, sample_size=sample_size, **kwargs)
    elif isinstance(arr_i2, array_wrapper):
        sample_1, sample_2 = resample(
            arr_i1, arr_i2=arr_i2, axis=axis, data_var_i1=data_var_i1, data_var_i2=data_var_i2, dim=dim,
            pseudo_sample_size=pseudo_sample_size, sample_method=sample_method, sample_size=sample_size)
    else:
        sample_1 = resample(
            arr_i1, axis=axis, data_var_i1=data_var_i1, data_var_i2=data_var_i2, dim=dim,
            pseudo_sample_size=pseudo_sample_size, sample_method=sample_method, sample_size=sample_size)
    # compute statistic
    if sample_2 is not None:
        kwargs.update({"arr_i2": sample_2})
    return statistic(sample_1, dim=dim, **kwargs)


def statistic(
        arr_i1: Union[array_wrapper, list, numpy__ndarray],
        arr_i2: Union[array_wrapper, list, numpy__ndarray] = None,
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        stat: str = "mea",
        **kwargs) -> Union[float, array_wrapper, array_wrapper, numpy__ndarray]:
    """
    Compute given statistic on input numpy.ndarray or xarray.DataArray

    Input:
    ------
    :param arr_i1: array_like
    :param arr_i2: array_like, optional
        Needed for two array statistics:
            correlation (statistic = 'cor'),
        Default is None
    :param axis: int or list[int] or tuple[int] or None, optional
        Axis(es) (numpy) along which to compute statistic; e.g., axis=0.
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array).
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] tuple[str] or None, optional
        Name of dimension(s) (xarray) along which to compute statistic; e.g., dim="time".
        Only one of the 'axis' and 'dim' arguments can be supplied.
        Default is None (compute statistic of the flattened array)
    :param stat: str, optional
        Nickname of the desired statistic (see statistics)
    **kwargs: dict
        The following options are available for statistics using standard deviation or variance (see specific function
        for details):
            'ddof': divisor used in the calculation is N - ddof, where N represents the number of elements
        The following options are available for statistics using kurtosis or skewness (see specific function for
        details):
            'bias': False to corrected for statistical bias
        The following options are available for percentile (statistic = 'per', see stat_percentile for details):
            'per': Percentile(s) at which to extract score. Values should be in range [0, 100]
        The following options are sometimes available for xarray (see specific function for details):
            'keep_attrs': True to copy the dataset's attributes (attrs) from the original object to the new one
            'method' (stat_percentile): Interpolation method to use when the desired quantile lies between two data
                points
            'weights' (stat_correlation_pearson): array of weights

    Output:
    -------
    :return: float or xarray.DataArray or numpy.ndarray
        The desired statistical value(s).
    """
    arr_i1, arr_i2 = numpy_tools.do_list_to_array(arr_i1), numpy_tools.do_list_to_array(arr_i2)
    # set axis and dimension (if axis is given find dimension name and conversely)
    if isinstance(arr_i1, array_wrapper):
        axis, dim = xarray_base.get_axis_dim(arr_i1, axis=axis, dim=dim)
    elif isinstance(arr_i2, array_wrapper):
        axis, dim = xarray_base.get_axis_dim(arr_i2, axis=axis, dim=dim)
    # known statistics
    stat_1_arr = {
        "kur": statistics.kurtosis, "mea": statistics.mean, "med": statistics.median, "per": statistics.percentile,
        "ske": statistics.skewness, "std": statistics.standard_deviation, "tin": statistics.student_confidence_interval,
        "var": statistics.variance}
    stat_2_arr = {"cor": statistics.correlation_pearson}
    # compute statistic
    if stat in list(stat_1_arr.keys()):
        arr_o = stat_1_arr[stat](arr_i1, axis=axis, dim=dim, **kwargs)
    else:
        arr_o = stat_2_arr[stat](arr_i1, arr_i2, axis=axis, dim=dim, **kwargs)
    return arr_o
# ---------------------------------------------------------------------------------------------------------------------#
