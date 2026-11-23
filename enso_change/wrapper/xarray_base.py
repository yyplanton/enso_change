# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# N-dimensional array with labeled coordinates and dimensions are needed for the package.
# In xarray:
#     - A dataset resembles an in-memory representation of a NetCDF file, and consists of variables, coordinates and
#       attributes which together form a self describing dataset.
#     - A DataArray provides a wrapper around numpy ndarrays that uses labeled dimensions and coordinates to support
#       metadata aware operations.
# This file builds on xarray functions.
# https://docs.xarray.dev/en/latest/
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from datetime import datetime as datetime__datetime
from inspect import stack as inspect__stack
from typing import Hashable, Literal, Union
# numpy
from numpy import array as numpy__array
from numpy import cos as numpy__cos
from numpy import deg2rad as numpy__deg2rad
from numpy import full as numpy__full
from numpy import nan as numpy__nan
from numpy import ndarray as numpy__ndarray
from numpy import zeros as numpy__zeros
# xarray
import xarray
import xskillscore

# local functions
from enso_change.tools.check_tools import check_is_dim
from enso_change.tools.default_tools import none_to_default
from enso_change.tools.print_tools import print_fail, unknown_formater
from enso_change.wrapper import numpy_tools
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
array_wrapper = xarray.DataArray
dataset_wrapper = xarray.Dataset


def align(
        ds_a: array_wrapper,
        ds_b: array_wrapper,
        copy: bool = True,
        fill_value: float = None,
        join: Literal["exact", "inner", "left", "outer", "override", "right"] = "inner",
        **kwargs) -> (array_wrapper, array_wrapper):
    """
    Returns new objects with aligned indexes and dimension sizes.
    Array from the aligned objects are suitable as input to mathematical operators, because along each dimension they
    have the same index and size.
    Missing values (if join != 'inner') are filled with fill_value. The default fill value is NaN.
    https://docs.xarray.dev/en/stable/generated/xarray.align.html

    Input:
    ------
    :param ds_a: xarray.DataArray
    :param ds_b: xarray.DataArray
    :param copy: bool, optional
        If copy=True, data in the return values is always copied. If copy=False and reindexing is unnecessary, or can be
        performed with only slice operations, then the output may share memory with the input. In either case, new
        xarray objects are always returned.
        Default is True
    :param fill_value: float or int or None, optional
        Value to use for newly missing values. If a dict-like, maps variable names to fill values. Use a data array’s
        name to refer to its values.
        Default is None
    :param join: {"exact", "inner", "left", "outer", "override", "right"}, optional
        Method for joining the indexes of the passed objects along each dimension:
            - “exact”: instead of aligning, raise ValueError when indexes to be aligned are not equal
            - “inner”: use the intersection of object indexes
            - “left”: use indexes from the first object with each dimension
            - “outer”: use the union of object indexes
            - “override”: if indexes are of same size, rewrite indexes to be those of the first object with that
                          dimension. Indexes for the same dimension must have the same size in all objects.
            - “right”: use indexes from the last object with each dimension
        Default is "inner"
     **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Tuple of objects with the same type as *objects with aligned coordinates.
    """
    tmp_kwargs = {} if fill_value is None else {"fill_value": fill_value}
    ds_ao, ds_bo = xarray.align(ds_a, ds_b, copy=copy, join=join, **tmp_kwargs)
    return ds_ao, ds_bo


def apply_ufunc(da: array_wrapper, func, dim: Union[Hashable, str], vectorize: bool = True, **kwargs):
    return xarray.apply_ufunc(func, da, input_core_dims=[[dim]], vectorize=vectorize, kwargs=kwargs)


def apply__kernel(da: array_wrapper, dim: Union[Hashable, str], nbr_samples: int) -> Union[array_wrapper, None]:
    da_o = None
    # fake loop to be able to break out when an error occurs
    for _ in [0]:
        # CF dimension name to DataArray dimension name
        if dim in ["T", "X", "Y"]:
            dim = cf_dim_to_dim(da, dim)
            if not check_is_dim(dim):
                break
        # find axis number corresponding to dim
        axis, dim = get_axis_dim(da, dim=dim)
        # compute kernel using numpy and scipy
        da_o = numpy_tools.kernel_nd(da.values, axis, nbr_samples)
        # numpy to DataArray
        coordinates = {"member": numpy__array(list(range(nbr_samples)))}
        da_o = recreate_array(da_o, da, axis_added=[0], coords_added=coordinates, data_var="gaussian_kde",
                              dim_added=["member"], dim_removed=["member"])
    return da_o


def cf_dim_to_dim(
        ds: array_wrapper,
        cf_dim: Literal["T", "X", "Y"],
        **kwargs) -> Union[Hashable, str, None]:
    """
    Return dimension name corresponding to CF dimension name (T is time, X is longitude, Y is latitude).

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    :param cf_dim: {"X", "Y", "T", "Z"}
        Name of a CF dimension
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str or None
        Name of given CF dimension in input xarray.DataArray or xarray.Dataset.
    """
    # list dimension
    list_dim = list(ds.coords)
    # dimension to find
    dim_to_find = {
        "T": ["time", "tim"],
        "X": ["longitude", "lon", "x"],
        "Y": ["latitude", "lat", "y"],
    }
    # find the name in the dataset
    dim_o = None
    for k1 in dim_to_find[cf_dim]:
        if k1 in list_dim:
            # k1 is a dimension name in ds
            dim_o = copy__deepcopy(k1)
            break
        else:
            # k1 is not an exact match -> check if k1 is included in a dimension name
            for k2 in list_dim:
                if k1 == "y" and k2 == "year":
                    continue
                # if k1 is a single letter it must be the first letter k2 (dimension name)
                #     e.g., k1 = 'x', if k2 = 'xt_ocean' or 'x_dim', k2 is probably the dimension we are looking for
                # otherwise, k1 must be included in k2 (dimension name)
                #     e.g., k1 = 'lon', if k2 = 'dim_lon' or 'dim_lon', k2 is probably the dimension we are looking for
                if (len(k1) == 1 and k2[:1] == k1) or (len(k1) > 1 and k1 in k2):
                    dim_o = copy__deepcopy(k2)
                    break
        if dim_o is not None and isinstance(dim_o, (Hashable, str)):
            break
    return dim_o


def check_dim(
        ds: array_wrapper,
        dim: Union[Hashable, str],
        **kwargs) -> Union[Hashable, str]:
    """
    Ensures that desired dimension exists in given object.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
    :param dim: Hashable or str
        Dimension name
    **kwargs - Discarded

    Output:
    -------
    :return: Hashable or str
        Dimension name in ds
    """
    dim_o = None
    if dim is None or not isinstance(dim, (Hashable, str)):
        pass
    elif dim in ds.dims:
        dim_o = copy__deepcopy(dim)
    elif dim in ["T", "X", "Y"]:
        dd: Literal["T", "X", "Y"] = "T" if dim == "T" else ("X" if dim == "X" else "Y")
        dim_o = cf_dim_to_dim(ds, dd)
    if dim_o is None or not isinstance(dim_o, (Hashable, str)):
        print_fail(inspect__stack(), unknown_formater("dimension", dim, ds.dims), fail_i=False)
    return dim_o


def check_multidimensional_coordinates(ds: array_wrapper, **kwargs) -> bool:
    """
    Check if input object uses multidimensional coordinates (e.g., longitude[y, x]).

    Input:
    ------
    :param ds: xarray.DataArray
    **kwargs - Discarded

    Output:
    -------
    :return: bool
        True is input object uses multidimensional coordinates
    """
    dim_lat = cf_dim_to_dim(ds, "Y")
    dim_lon = cf_dim_to_dim(ds, "X")
    bool_o = False
    if (check_is_dim(dim_lat) and len(ds[dim_lat].shape) > 1) or (check_is_dim(dim_lon) and len(ds[dim_lon].shape) > 1):
        bool_o = True
    return bool_o


def concatenate(
        obj, dim, **kwargs):
    known_kwargs = ["combine_attrs", "compat", "coords", "create_index_for_new_dim", "data_var", "fill_value", "join",
                    "positions"]
    for k in list(kwargs.keys()):
        if k not in known_kwargs:
            del kwargs[k]
    return xarray.concat(obj, dim, **kwargs)


def correlation(
        ds_a: array_wrapper,
        ds_b: array_wrapper,
        data_var_a: str = None,
        data_var_b: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> array_wrapper:
    """
    Compute the Pearson correlation coefficient between two DataArray objects along shared dimension(s).
    https://docs.xarray.dev/en/latest/generated/xarray.corr.html

    Input:
    ------
    :param ds_a: xarray.DataArray
    :param ds_b: xarray.DataArray
    :param data_var_a: str, optional
        Data variable in ds_a; e.g., data_var_a = "ts".
        Default is None
    :param data_var_b: str, optional
        Data variable in ds_b; e.g., data_var_b = "ts".
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension[s] along which to apply var; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param weights: xarray.DataArray, optional
        Array of weights.
        Default is None
     **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New object with dimensions, attributes, coordinates, name, encoding, with correlation applied to given data and
        the indicated dimension(s) removed.
    """
    # correlation value
    ds_o = xarray.corr(ds_a, ds_b, dim=dim, weights=weights)
    # attributes
    att_a, att_b = ds_a.attrs, ds_b.attrs
    for k in ["freq", "mode", "operation", "units", "weighted"]:
        if k in list(att_a.keys()):
            del att_a[k]
        if k in list(att_b.keys()):
            del att_b[k]
    attributes = copy__deepcopy(att_a)
    for k1, k2 in att_b.items():
        if k1 in list(attributes.keys()) and k1 == "long_name":
            attributes[k1] = str(att_a[k1]).split(" time series")[0] + " & " + str(k2).split(" time series")[0] + \
                             " correlation"
        elif k1 in list(attributes.keys()) and k1 == "short_name":
            attributes[k1] = "corr (" +str(att_a[k1]) + " & " + str(k2) + ")"
        elif k1 in list(attributes.keys()) and att_a[k1] != k2:
            attributes[k1] = str(data_var_a) + " -- " + str(att_a[k1]) + "\n" + str(data_var_b) + " -- " + str(k2)
        else:
            attributes[k1] = k2
    attributes["units"] = ""
    ds_o.attrs.update(**attributes)
    return ds_o


def covariance(
        ds_a: array_wrapper,
        ds_b: array_wrapper,
        data_var_a: str = None,
        data_var_b: str = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        weights: Union[array_wrapper, None] = None,
        **kwargs) -> array_wrapper:
    """
    Compute covariance between two DataArray objects along a shared dimension.
    https://docs.xarray.dev/en/stable/generated/xarray.cov.html

    Input:
    ------
    :param ds_a: xarray.DataArray
    :param ds_b: xarray.DataArray
    :param data_var_a: str, optional
        Data variable in ds_a; e.g., data_var_a = "ts".
        Default is None
    :param data_var_b: str, optional
        Data variable in ds_b; e.g., data_var_b = "ts".
        Default is None
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Name of dimension[s] along which to apply cov; e.g., dim="x" or dim=["x", "y"].
        If None, will reduce over all dimensions.
        Default is None
    :param weights: xarray.DataArray, optional
        Array of weights.
        Default is None
     **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New object with dimensions, attributes, coordinates, name, encoding, with correlation applied to given data and
        the indicated dimension(s) removed.
    """
    # covariance value
    ds_o = xarray.cov(ds_a, ds_b, dim=dim, weights=weights)
    # attributes
    att_a, att_b = ds_a.attrs, ds_b.attrs
    for k in ["freq", "mode", "operation", "weighted"]:
        if k in list(att_a.keys()):
            del att_a[k]
        if k in list(att_b.keys()):
            del att_b[k]
    attributes = copy__deepcopy(att_a)
    for k1, k2 in att_b.items():
        if k1 in list(attributes.keys()) and k1 == "long_name":
            attributes[k1] = str(att_a[k1]).split(" time series")[0] + " & " + str(k2).split(" time series")[0] + \
                             " covariance"
        elif k1 in list(attributes.keys()) and k1 == "short_name":
            attributes[k1] = "cov (" +str(att_a[k1]) + " & " + str(k2) + ")"
        elif k1 in list(attributes.keys()) and att_a[k1] != k2:
            attributes[k1] = str(data_var_a) + " -- " + str(att_a[k1]) + "\n" + str(data_var_b) + " -- " + str(k2)
        else:
            attributes[k1] = k2
    uni_a = att_a["units"] if "units" in list(att_a.keys()) else ""
    uni_b = att_b["units"] if "units" in list(att_b.keys()) else ""
    attributes["units"] = str(uni_a) + "." + str(uni_b) if len(uni_a) > 0 and len(uni_b) > 0 else (
        str(uni_a) if len(uni_a) > 0 else str(uni_b))
    ds_o.attrs.update(**attributes)
    return ds_o


def create_array(
        ds: array_wrapper,
        data_var: str = "",
        value: Union[float, int, None] = 0,
        **kwargs) -> array_wrapper:
    """
    Return a new DataArray of value with the same shape, axes, coordinates, attributes,... as input DataArray.
    https://docs.xarray.dev/en/latest/generated/xarray.zeros_like.html

    Input:
    ------
    :param ds: xarray.DataArray
    :param data_var: str, optional
        Name of the output data variable.
        Default is ""
    :param value: float or int or None, optional
        Value to fill the array with; e.g., value = 1.
        If None, array filled with NaN.
        Default is 0
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New DataArray of ones with the same shape and type as ds.
    """
    # create an array of zeros with the same shape, axes, coordinates, attributes,... as input DataArray
    da_o = xarray.zeros_like(ds)
    # name this new array
    da_o.name = data_var
    if value is None:
        # fill with NaN
        da_o = da_o.where(da_o == 1)
    elif isinstance(value, (float, int)) and value != 0:
        # fill with value
        da_o = da_o.where(da_o == value, other=value)
    return da_o


def create_array_time_series(
        dict_i: dict[str, dict[int, list[Union[float, int]]]],
        attrs: Union[dict[str, str], None] = None,
        data_var: Union[str, None] = None,
        **kwargs) -> array_wrapper:
    """
    Create a DataArray from data extracted from json files

    :param dict_i: dict[str, dict[str, list[Union[float, int]]]]
        Dictionary with member and epoch (as int) as keys
    :param attrs: dict[str, str], optional
        Variable attributes for DataArray; e.g, attrs = {"attr_name": "this is an attribute"}
        Default is None (no attribute)
    :param data_var: str, optional
        Name of the variable (i.e., diagnostic name); data_var = "ave_ts_val_n30e"
        Default is None (no name)
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Array [time, member] for a dataset
    """
    attrs = none_to_default(attrs, {})
    data_var = none_to_default(data_var, "no_name")
    # list keys lev 1
    list_k1 = list(dict_i.keys())
    # create array
    if isinstance(dict_i[list_k1[0]], dict):
        # -> two nested levels, they are members then years
        # list years
        list_k2 = sorted(list(dict_i[list_k1[0]].keys()))
        # if exp is piControl, the number of years may differ from one member to another
        max_length = max([len(list(dict_i[mem].keys())) for mem in list_k1])
        # output array
        arr = numpy__full((len(list_k1), max_length), numpy__nan)
        for ii, mem in enumerate(list_k1):
            for jj, yea in enumerate(sorted(list(dict_i[mem].keys()))):
                arr[ii, jj] = dict_i[mem][yea]
        # member dimension is the 'r' number
        list_k1 = [int(mem[1:].split("i")[0]) for mem in list_k1]
        if len(list_k2) != max_length:
            # the time axis of the first member is used and extended if too short
            yn = list_k2[-1]
            dy = list_k2[1] - list_k2[0]
            list_k2 += [yn + (ii + 1) * dy for ii in range(max_length - len(list_k2))]
        # list dimensions
        dimensions = ["member", "year"]
        # get coordinates corresponding to dimensions
        coordinates = {"member": numpy__array(list_k1), "year": numpy__array(list_k2)}
    else:
        # -> one level, it is years
        list_k1 = sorted(list_k1)
        # output array
        arr = numpy__array([dict_i[yea] for yea in list_k1])
        # list dimensions
        dimensions = ["year"]
        # get coordinates corresponding to dimensions
        coordinates = {"year": numpy__array(list_k1)}
    # numpy.ndarray to xarray.DataArray
    return array_wrapper(attrs=attrs, coords=coordinates, data=arr, dims=dimensions, name=data_var)


def create_time_plot(ds: array_wrapper, year_start: int = None, **kwargs) -> Union[numpy__ndarray, None]:
    """
    Create a time axis based on 'year fraction', e.g., 2025-01-16 0:0:0 ~ 2025.0.04, 2025-02-14 0:0:0 ~ 2025.12
    Always uses proleptic Gregorian calendar

    Input:
    ------
    :param ds: xarray.DataArray
    :param year_start: int
        Start year of the time axis
    **kwargs - Discarded

    Output:
    -------
    :return: numpy__ndarray
        Time axis based on 'year fraction'
    """
    time_o = None
    # time dimension
    dim_tim = cf_dim_to_dim(ds, "T")
    if check_is_dim(dim_tim):
        # get time array
        time = ds[dim_tim].to_index()
        # create output array
        time_o = numpy__zeros(len(time))
        # time as year fraction
        dy = year_start - time.year[0]
        for k, (y, m, d, h, mi, se) in enumerate(
                zip(time.year, time.month, time.day, time.hour, time.minute, time.second)):
            # create datetime for this year and the current time
            year_start = datetime__datetime(y + dy, 1, 1, 0, 0, 0)
            year_end = datetime__datetime(y + dy + 1, 1, 1, 0, 0, 0)
            year_now = datetime__datetime(y + dy, m, d, h, mi, se)
            # number of seconds from the beginning of the year and total number of seconds during this year
            sec_now = (year_now - year_start).total_seconds()
            sec_tot = (year_end - year_start).total_seconds()
            # time is year + second fraction
            time_o[k] = y + dy + sec_now / sec_tot
    return time_o


def flatten_to_list(
        arr_i: Union[array_wrapper, float, int, list, numpy__ndarray, list],
        list_values: list = None,
        **kwargs) -> list:
    """
    Flatten to a list of values

    Input:
    ------
    :param arr_i: array_like
    :param list_values: list, optional
    **kwargs - Discarded

    Output:
    -------
    :return list_values: list
        Flatten list
    """
    # set default values if not given
    list_values = none_to_default(list_values, [])
    # if list contains list, flatten
    if isinstance(arr_i, (array_wrapper, numpy__ndarray)):
        list_values += list(arr_i.flatten())
    elif isinstance(arr_i, list):
        if all(isinstance(k, (float, int)) for k in arr_i):
            list_values += arr_i
        else:
            for k in arr_i:
                flatten_to_list(k, list_values=list_values)
    else:
        list_values.append(arr_i)
    return list_values


def get_axis_dim(
        ds: array_wrapper,
        axis: Union[int, list[int], tuple[int], None] = None,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        **kwargs) -> Union[tuple[int, str], tuple[int, Hashable], tuple[list[int], list[str]]]:
    """
    Return axis(es) corresponding to given dimension name(s) and conversely.

    Input:
    ------
    :param ds: xarray.DataArray
    :param axis: int or list[int] or tuple[int], optional
        Axis(es) for which the dimension name(s) must be found; e.g., axis = 0
    :param dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str], optional
        Dimension(s) for which the axis number(s) must be found; e.g., dim = "time"
        Only one of the 'axis' and 'dim' arguments can be supplied
    **kwargs - Discarded

    Outputs:
    --------
    :return axis: int or list[int] or tuple[int]
        Axis(es) corresponding to given dimension name(s) or input axis(es)
    :return dim: Hashable or str or list[Hashable] or list[str] or tuple[Hashable] or tuple[str]
        Dimension(s) corresponding to given axis(es) or input dimension(s)
    """
    # get dimension listed as in xarray.DataArray
    list_dim = list(ds.dims)
    if dim is not None:
        # get axis(es) corresponding to dimension name(s)
        if isinstance(dim, str):
            axis = list_dim.index(dim)
        else:
            axis = [list_dim.index(k) for k in dim]
    elif axis is not None:
        # get dimension name(s) corresponding to axis(es)
        if isinstance(axis, int):
            dim = list_dim[axis]
        else:
            dim = [list_dim[k] for k in axis]
    return axis, dim


def mean_weighted(
        da: array_wrapper,
        dim: Union[Hashable, str, list[Hashable], list[str], tuple[Hashable], tuple[str]] = "T",
        keep_attrs: bool = True,
        skipna: bool = False,
        weights: Union[array_wrapper, numpy__ndarray, None] = None,
        **kwargs):
    """
    Reduce this DataArray’s data by a weighted mean along given dimension(s).
    https://docs.xarray.dev/en/stable/generated/xarray.DataArray.weighted.html
    https://docs.xarray.dev/en/stable/generated/xarray.computation.weighted.DataArrayWeighted.mean.html
    https://docs.xarray.dev/en/stable/generated/xarray.DataArray.mean.html

    Input:
    ------
    :param da: xarray.DataArray
    :param dim: Hashable or Iterable of Hashable, optional
        Dimension(s) over which to apply the weighted mean; e.g., dim = "X" or dim = ("X", "Y")
        Default is "T"
    :param keep_attrs: bool or None, optional
        If True, attrs will be copied from the original object to the new one.
        If False, the new object will be returned without attributes.
        Default is True
    :param skipna: bool or None, optional
        If True, skip missing values (as marked by NaN).
        Default is None (i.e., only skips missing values for float dtypes)
    :param weights: xarray.DataArray or numpy.ndarray or None, optional
        An array of weights associated with the values in this Dataset. Each value in the data contributes to the
        reduction operation according to its associated weight.
        If not given and latitude in dim, weights are computed using cos(latitude).
        If not given and time in dim, weights are computed using time_weights.
        Then da.weighted(weights).mean() is computed if weights are instance of xarray.DataArray or numpy.ndarray
        Else, ds.mean() is computed.
        Default is None
    **kwargs – Additional keyword arguments passed on to the appropriate array function for calculating mean on this
               object's data. These could include dask-specific kwargs like split_every.

    Output:
    -------
    :return: xarray.DataArray
         New DataArray object with mean or weighted mean applied to its data and the indicated dimension(s) removed.
    """
    da_o = None
    # fake loop to be able to break out when an error occurs
    for _ in [0]:
        # get dimension name in da
        if dim is not None and isinstance(dim, (list, tuple)):
            dim = [check_dim(da, k) for k in dim]
        else:
            dim = check_dim(da, dim)
        # create an iterable of dimensions to facilitate tests
        dim_list = copy__deepcopy(dim) if isinstance(dim, (list, tuple)) else [dim]
        if not all([check_is_dim(k) for k in dim_list]):
            break
        # get 'theoretical' dimension names of time, longitude and latitude (None is returned if they don't exist)
        time_dim = cf_dim_to_dim(da, "T")
        lon_dim = cf_dim_to_dim(da, "X")
        lat_dim = cf_dim_to_dim(da, "Y")
        # compute weights if not provided and if possible
        if isinstance(weights, bool) and weights:
            # list longitude and latitude 'theoretical' dimension names if they exist
            hor_dim = [k for k in [lon_dim, lat_dim] if check_is_dim(k)]
            # check if weights should be computed
            if check_is_dim(lat_dim) and lat_dim in dim_list and len(list(set(dim_list) - set(hor_dim))) == 0:
                # compute weights using cos(latitude) only if dim_list = ['Y'] or ['X', 'Y']
                weights = numpy__cos(numpy__deg2rad(da[lat_dim]))
            # elif check_is_dim(time_dim) and time_dim in dim_list and len(dim_list) == 1:
            #     # compute time weights using days per month only if dim_list = ['T']
            #     weights = xt.weights_time(da)
        # test if multidimensional coordinates
        if (check_is_dim(lat_dim) and lat_dim in dim_list and len(da[lat_dim].shape) > 1) or (
                check_is_dim(lon_dim) and lon_dim in dim_list and len(da[lon_dim].shape) > 1):
            # is multidimensional
            dim_old = copy__deepcopy(dim) if isinstance(dim, (list, tuple)) else [dim]
            # convert lat/lon to y/x
            dim_new = []
            for k in dim_old:
                if check_is_dim(lat_dim) and k == lat_dim and len(da[lat_dim].shape) > 1:
                    dim_new.append(list(da[lat_dim].dims)[0])
                elif check_is_dim(lon_dim) and k == lon_dim and len(da[lon_dim].shape) > 1:
                    dim_new.append(list(da[lat_dim].dims)[1])
                else:
                    dim_new.append(k)
            dim = copy__deepcopy(dim_new) if isinstance(dim, (list, tuple)) else dim_new[0]
        # compute mean
        if isinstance(da, array_wrapper) and isinstance(weights, (array_wrapper, numpy__ndarray)):
            da_o = da.weighted(weights).mean(dim=dim, keep_attrs=keep_attrs, skipna=skipna)
        else:
            da_o = da.mean(dim=dim, keep_attrs=keep_attrs, skipna=skipna, **kwargs)
    return da_o


def netcdf_open(
        filename: Union[str, list[str]],
        variable: str,
        bounds: dict[str, tuple],
        **kwargs) -> Union[dict[str, array_wrapper], None]:
    """
    Open file(s) and read given variable
    https://docs.xarray.dev/en/latest/generated/xarray.open_mfdataset.html

    Input:
    ------
    :param filename: str or list[str]
        Paths to dataset files:
            - Directory path (e.g., "path/to/files"), which is converted to a string glob of *.nc files
            - String glob (e.g., "path/to/files/*.nc"), which is expanded to a 1-dimensional list of file paths
            - File path to dataset (e.g., "path/to/files/file1.nc")
            - List of file paths (e.g., ["path/to/files/file1.nc", ...])
        If multiple files, concatenation along the time dimension recommended
    :param variable: str
        Name of a variable to read in filename; e.g., variable = 'sst'
    :param bounds: dict[str, tuple]
        Dictionary of desired bounds to select along the specified dimension(s);
        e.g., bounds = {"X": (160., 280.)}
    **kwargs - Discarded

    Output:
    -------
    :return: dict[str, xarray.DataArray]
        Dictionary with variables as keys and xarray.DataArray as values
    """
    # create dataset
    ds = xarray.open_mfdataset(filename, coords="different")
    # read dataarray
    dict_o = None
    if variable is not None and isinstance(variable, str) and variable in list(ds.keys()):
        dict_o = {variable: ds[variable]}
    elif variable is not None and isinstance(variable, str) and variable == "*":
        dict_o = dict((k, ds[k]) for k in list(ds.keys()) if "bound" not in k and "bnd" not in k)
    # select region (if given)
    if isinstance(dict_o, dict) and len(list(dict_o.keys())) > 0 and isinstance(bounds, dict) and \
            len(list(bounds.keys())) > 0:
        for va, da in dict_o.items():
            # test if input is a map or a time series
            if cf_dim_to_dim(da, "X") is None and cf_dim_to_dim(da, "Y") is None:
                dict_o[va] = da
                continue
            # roll longitude if necessary
            if "X" in list(bounds.keys()) and (min(bounds["X"]) < 0 or max(bounds["X"]) > 360):
                da = roll_longitude(da, new_lon_min=min(bounds["X"]))
            # select
            coordinates = {}
            for k1, k2 in bounds.items():
                dim = cf_dim_to_dim(da, k1)
                if not check_is_dim(dim):
                    break
                coordinates[dim] = slice(*k2)
            if len(list(bounds.keys())) != len(list(coordinates.keys())):
                dict_o = None
                break
            if not check_multidimensional_coordinates(da):
                # normal selection method
                da = da.sel(indexers=coordinates)
            else:
                # for multidimensional coordinates (e.g., curvilinear grids)
                dim_lon = cf_dim_to_dim(da, "X")
                dim_lat = cf_dim_to_dim(da, "Y")
                # create condition
                if "Y" in list(bounds.keys()) and "X" not in list(bounds.keys()):
                    # select Y
                    arr_lat = da[dim_lat]
                    cond = (bounds["Y"][0] < arr_lat) & (arr_lat < bounds["Y"][1])
                elif "Y" not in list(bounds.keys()) and "X" in list(bounds.keys()):
                    # select X
                    arr_lon = da[dim_lon]
                    cond = (bounds["X"][0] < arr_lon) & (arr_lon < bounds["X"][1])
                else:
                    # select X / Y
                    arr_lat, arr_lon = da[dim_lat], da[dim_lon]
                    cond = (bounds["Y"][0] <= arr_lat) & (arr_lat <= bounds["Y"][1]) & \
                           (bounds["X"][0] <= arr_lon) & (arr_lon <= bounds["X"][1])
                da = da.where(cond)
            dict_o[va] = da
    return dict_o


def polyval(
        da_coords: array_wrapper,
        da_coeffs: array_wrapper,
        degree_dim: str = "degree",
        **kwargs) -> array_wrapper:
    """
    Evaluate a polynomial at specific values.
    https://docs.xarray.dev/en/stable/generated/xarray.polyval.html

    Input:
    ------
    :param da_coords: xarray.DataArray
        Values at which to evaluate the polynomial.
    :param da_coeffs: xarray.DataArray
        Coefficients of the polynomial.
    :param degree_dim: Hashable or str, optional
        Name of the polynomial degree dimension in da_coeffs.
        Default is 'degree'
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Evaluated polynomial.
    """
    return xarray.polyval(da_coords, da_coeffs, degree_dim=degree_dim)


def recreate_array(
        arr: numpy__ndarray,
        ds: array_wrapper,
        attrs_added: dict[str, str] = None,
        attrs_removed: list[str] = None,
        axis_added: Union[list[int], tuple[int], None] = None,
        coords_added: dict[str, Union[numpy__ndarray, array_wrapper]] = None,
        data_var: str = None,
        dim_added: Union[list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        dim_removed: Union[list[Hashable], list[str], tuple[Hashable], tuple[str], None] = None,
        **kwargs) -> array_wrapper:
    """
    Recreate output xarray.DataArray from input xarray.DataArray.

    Input:
    ------
    :param arr: numpy.ndarray
        Array derived from ds (e.g., a statistic was computed) that was transformed into a numpy.ndarray in the
        process
    :param ds: xarray.DataArray
        Original xarray.DataArray from which arr is derived
    :param attrs_added: dict[str, str] or None, optional
        Variable attributes to add to the DataArray; e.g, attrs_added = {"attr_name": "new attribute"}.
        Default is None (no attribute added)
    :param attrs_removed: list[str] or None, optional
        Variable attributes to remove from the DataArray; e.g, attrs_removed = ["attr_name"].
        Default is None (no attribute removed)
    :param axis_added: list[int] or tuple[int] or None
        Position(s) of added dimension(s) (if any); e.g, axis_added = [0] or axis_added = [0, 1].
        If given, both axis_added and dim_added must be provided.
        Default is None (no dimension has been added)
    :param coords_added: dict[str, numpy.ndarray or xarray.DataArray] or None, optional
        Coordinates (tick labels) to use for indexing along each dimension.
        If dim_added is not in coords_added, coordinates will be a sequence of numbers.
        Default is None
    :param data_var: str, optional
        Name of the output data variable.
        Default is ""
    :param dim_added: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been added from ds to arr (if any);
        e.g, dim_added = ["x"] or dim_added = ["x", "y"].
        If given, both axis_added and dim_added must be provided.
        Default is None (no dimension has been added)
    :param dim_removed: list[Hashable] or list[str] or tuple[Hashable] or tuple[str] or None, optional
        Dimension name(s) that has been removed from ds to arr (e.g., to compute a statistic);
        e.g., dim_removed = ["x"] or dim_removed = ["x", "y"].
        Default is None (no dimension was removed)
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Input arr wrapped in a xarray.DataArray.
    """
    dim_removed = none_to_default(dim_removed, [])
    # list dimensions
    dimensions = list(ds.dims)
    # delete removed dimension(s)
    for k in dim_removed:
        if isinstance(k, str) and k in dimensions:
            dimensions.remove(k)
    # get coordinates corresponding to dimensions
    coordinates = dict((k, ds[k]) for k in dimensions)
    # add given dimension(s)
    if isinstance(axis_added, (list, tuple)) and isinstance(dim_added, (list, tuple)) and \
            len(axis_added) == len(dim_added):
        for k1, k2 in zip(axis_added, dim_added):
            # add given dimension at given position
            dimensions.insert(k1, k2)
            # add coordinates in dictionary
            if isinstance(coords_added, dict) and k2 in coords_added.keys():
                coordinates[k2] = coords_added[str(k2)]
            else:
                coordinates[k2] = array_wrapper(numpy__array(list(range(arr.shape[k1]))))
    # get input attributes
    attributes = ds.attrs
    # remove attributes
    if isinstance(attrs_removed, list):
        for k in attrs_removed:
            if k in list(attributes.keys()):
                del attributes[k]
    # add given attribute(s)
    if isinstance(attrs_added, dict):
        attributes.update(attrs_added)
    attributes = dict((k, attributes[k]) for k in sorted(list(attributes.keys()), key=lambda v: v.lower()))
    # numpy.ndarray to xarray.DataArray
    ds_o = array_wrapper(attrs=attributes, coords=coordinates, data=arr, dims=dimensions, name=data_var)
    return ds_o


def roll_longitude(
        ds: array_wrapper,
        new_lon_min: Union[float, int, None] = None,
        **kwargs) -> array_wrapper:
    """
    Roll longitude dimension, either to ensure that longitude ranges from 0 to 360E or to start from new_lon_min (e.g.,
    new_lon_min = -70, longitude = [-70; 290])

    Input:
    ------
    :param ds: xarray.DataArray
    :param new_lon_min: float or int or None
        New minimum longitude
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New object (as input) with rolled longitude
    """
    dim_lon = cf_dim_to_dim(ds, "X")
    if check_is_dim(dim_lon):
        # ensure that longitude ranges from 0 to 360E
        ds = ds.assign_coords({dim_lon: (360 + (ds[dim_lon] % 360)) % 360})
        # update longitude
        if isinstance(new_lon_min, (float, int)) and new_lon_min != 0:
            # longitude DataArray
            da_lon = ds[dim_lon]
            # change longitude DataArray according to new_lon_min
            if new_lon_min < 0:
                # e.g., initial longitude = [0; 360], new_lon_min = -70, new longitude = [-70; 290]
                da_lon = da_lon.where(da_lon < 360 + new_lon_min, other=da_lon - 360)
            else:
                # e.g., initial longitude = [0; 360], new_lon_min = 20, new longitude = [20; 380]
                da_lon = da_lon.where(da_lon >= new_lon_min, other=da_lon + 360)
            # set new longitude array in input object
            ds = ds.assign_coords({dim_lon: da_lon})
        # roll so that the first longitude of the dimension is the minimum longitude
        if not check_multidimensional_coordinates(ds):
            # normal roll method
            ds = ds.roll({dim_lon: -ds[dim_lon].argmin().values}, roll_coords=True)
        else:
            # for multidimensional coordinates (e.g., curvilinear grids)
            # average lon along Y
            arr_lon = ds[dim_lon]
            lon_x = arr_lon.to_numpy().mean(axis=0)
            # find minimum value
            min_x = lon_x.argmin()
            # roll
            ds = ds.roll({arr_lon.dims[-1]: -min_x}, roll_coords=True)
    return ds


def set_dim_as_first_or_last(
        ds: array_wrapper,
        dim: Union[Hashable, str],
        as_first: bool = True,
        **kwargs) -> Union[array_wrapper, None]:
    """
    Return a new object (as input) with all array dimensions transposed.
    https://docs.xarray.dev/en/stable/generated/xarray.DataArray.transpose.html

    Input:
    ------
    :param ds: xarray.DataArray
    :param dim: Hashable or str
        Dimension name to transpose
    :param as_first: bool
        True to set given dimension as first, else set given dimension as last
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        Object (as input) transposed dimensions.
    """
    # get dimension name
    dim_name = check_dim(ds, dim)
    # list dimensions in ds
    list_dim = list(ds.dims)
    if check_is_dim(dim_name):
        if len(list_dim) == 1 or (as_first and list_dim[0] == dim_name) or \
                (not as_first and list_dim[-1] == dim_name):
            return ds
        elif as_first and list_dim[0] != dim_name:
            # set dim as the first dimension of ds
            return ds.transpose(dim_name, "...")
        elif not as_first and list_dim[-1] != dim_name:
            # set dim as the last dimension of ds
            order = tuple(k for k in list_dim if k != dim_name) + (dim_name,)
            return ds.transpose(*order)
    return None


def to_netcdf(
        ds: Union[array_wrapper, dataset_wrapper],
        filename: str,
        file_format: Literal["NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", "NETCDF3_CLASSIC"] = "NETCDF4",
        list_of_variables: list[str] = None,
        mode: Literal["a", "w"] = "a",
        **kwargs):
    """
    Write dataset contents to a netCDF file.
    https://docs.xarray.dev/en/latest/generated/xarray.Dataset.to_netcdf.html

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param filename: str
        Path and filename to save this dataset
    :param file_format: {"NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", "NETCDF3_CLASSIC"}, optional
        File format for the resulting netCDF file:
            - NETCDF4: Data is stored in an HDF5 file, using netCDF4 API features.
            - NETCDF4_CLASSIC: Data is stored in an HDF5 file, using only netCDF 3 compatible API features.
            - NETCDF3_64BIT: 64-bit offset version of the netCDF 3 file format, which fully supports 2+ GB files, but is
                             only compatible with clients linked against netCDF version 3.6.0 or later.
            - NETCDF3_CLASSIC: The classic netCDF 3 file format. It does not handle 2+ GB files very well.
        All formats are supported by the netCDF4-python library. scipy.io.netcdf only supports the last two formats.
        Default is "NETCDF4"
    :param list_of_variables: list[str]
        List of data variable(s) in ds to save
    :param mode: {"w", "a"}, optional
        Write ("w") or append ("a") mode.
        If mode="w", any existing file at this location will be overwritten.
        If mode="a", existing variables will be overwritten.
        Default is "a"
    **kwargs - Discarded
    """
    # Write dataset contents to a netCDF file
    if isinstance(ds, array_wrapper) and isinstance(list_of_variables, list) and len(list_of_variables) > 0:
        ds[list_of_variables].to_netcdf(filename, format=file_format, mode=mode)
    else:
        ds.to_netcdf(filename, format=file_format, mode=mode)
# ---------------------------------------------------------------------------------------------------------------------#
