# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Functions built over xeofs
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from typing import Union, Hashable

# xeofs
import xeofs

# local functions
from enso_change.wrapper.xarray_base import cf_dim_to_dim
from enso_change.wrapper.xarray_tools import array_wrapper
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def eofs(
        ds: array_wrapper,
        dim: Union[Hashable, str] = None,
        **kwargs) -> Union[array_wrapper, None]:
    """
    Compute EOF eigenvectors along given dimension.

    Input:
    ------
    :param ds: xarray.DataArray or xarray.Dataset
        DataArray or Dataset
    :param dim: Hashable or str
        Name of dimension along which to apply EOF; e.g., dim = "time"
    **kwargs - Key arguments passed to eofs_compute

    Output:
    -------
    :return: xarray.DataArray or xarray.Dataset or None
        New object as input with EOF eigenvectors computed along the specified dimension
        None is returned if input dim is CF dim and wasn't found
    """
    # find dimension name if CF dim
    dim = cf_dim_to_dim(ds, dim) if dim in ["T", "X", "Y"] else copy__deepcopy(dim)
    # compute EOF eigenvectors along given dimension
    ds_o = None
    # compute eofs
    if isinstance(dim, (Hashable, str)) and isinstance(ds, array_wrapper):
        ds_o = eofs_compute(ds, dim, **kwargs)
    return ds_o


def eofs_compute(
        da: array_wrapper,
        dim: Union[Hashable, str] = None,
        kwargs_fit: dict = None,
        kwargs_std: dict = None,
        kwargs_xeofs: dict = None,
        **kwargs):
    """
    Compute EOF eigenvectors along given dimension.
    https://xeofs.readthedocs.io/en/develop/content/api_reference/_autosummary/xeofs.single.EOF.html
    https://docs.xarray.dev/en/latest/generated/xarray.DataArray.std.html

    Input:
    ------
    :param da: xarray.DataArray
    :param dim: Hashable or str
        Name of dimension along which to apply EOF; e.g., dim = "time"
        This must be a dimension in da (i.e., cannot be a CD dimension)
    :param kwargs_fit: dict, optional
        Key arguments passed to the function xeofs.single.EOF.fit (see url above)
        Default is None
    :param kwargs_std: dict, optional
        Key arguments passed to the function xarray std (see url above)
        Default is None
    :param kwargs_xeofs: dict, optional
        Key arguments passed to the function xeofs.single.EOF (see url above)
        Default is None
    **kwargs - Discarded

    Output:
    -------
    :return: xarray.DataArray
        New object as input with EOF eigenvectors computed along the specified dimension
    """
    # set key arguments dictionary to instance of dictionary if the wrong instance was given
    if not isinstance(kwargs_xeofs, dict):
        kwargs_xeofs = {}
    if not isinstance(kwargs_fit, dict):
        kwargs_fit = {}
    if not isinstance(kwargs_std, dict):
        kwargs_std = {}
    # set some default values to kwargs_xeofs
    default = {"n_modes": 1, "use_coslat": True}
    for k1, k2 in default.items():
        if k1 not in list(kwargs_xeofs.keys()):
            kwargs_xeofs[k1] = copy__deepcopy(k2)
    # generate EOF object
    eof_object = xeofs.single.EOF(**kwargs_xeofs)
    # fit EOF object to given DataArray
    fit = eof_object.fit(da, dim, **kwargs_fit)
    # compute eigenvectors (components), PCs (scores) and explained variance
    components = fit.components()
    scores = fit.scores()
    explained_variance = fit.explained_variance_ratio().values * 100
    explained_variance = [float(k) for k in explained_variance]
    # compute scores' std along given dimension
    scores_std = scores.std(dim=dim, **kwargs_std)
    # multiply components by scores std (to have components as input units)
    da_o = components * scores_std
    da_o.attrs.update({"explained_variance": explained_variance})
    # input dimension attrs on new array
    for dim in da.dims:
        if dim in da_o.dims:
            da_o[dim].attrs.update(da[dim].attrs)
    return da_o
# ---------------------------------------------------------------------------------------------------------------------#
