# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Processors built over xarray
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from glob import iglob as glob__iglob
from inspect import stack as inspect__stack
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
from typing import Union

# local package
from enso_change.tools import base_tools as bt
from enso_change.tools import check_list, check_type, flatten_dict, none_to_default, print_fail, put_in_dict, \
    sort_dict, tuple_for_dict
from enso_change.wrapper.xarray_base import array_wrapper, create_array_time_series, netcdf_open, numpy__ndarray
from enso_change.wrapper import xarray_tools
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def compute_1var(
        dict_i: Union[dict, array_wrapper, numpy__ndarray],
        function_to_compute: str,
        function_arguments: dict = None,
        dict_o: dict = None,
        tuple_k: tuple = None,
        tuple_k_last: tuple = None) -> (dict, tuple, tuple):
    """
    Compute the given 'function_to_compute' using given 'function_arguments' on the array at the bottom of the given
    'dict_i'

    Input:
    ------
    :param dict_i: dict or array_wrapper or numpy__ndarray
        Dictionary with nested levels, filled with an array
    :param function_to_compute: str
        Name of a function to compute; e.g., function_to_compute = 'statistic'
        The following options are available: see enso_and_trend.wrapper.xarray_tools
    :param function_arguments: dict, optional
        Dictionary of arguments for the function; e.g., function_arguments = {'axis': 0, 'statistic': 'mea'}
        The following options are available: see enso_and_trend.wrapper.xarray_tools
        Default is None (no argument passed to the function)
    :param dict_o: dict or None, optional
        Dictionary in which output values will be stored
    :param tuple_k: tuple or None, optional
        List of keys, in order, for the output nested dictionary
    :param tuple_k_last: tuple or None, optional
        List of the last key of each nested level, in order, to keep track of the position within the input nested
        dictionary

    Output:
    -------
    :return dict_o: dict
        Dictionary with the same nested levels as input dictionary filled with the output of function_to_compute
    :return tuple_k: tuple
    :return tuple_k_last: tuple
    """
    # list statistics functions
    list_fct = [k for k in dir(xarray_tools) if "__" not in k]
    # set None input to its default value
    function_arguments = none_to_default(function_arguments, {})
    dict_o = none_to_default(dict_o, {})
    tuple_k = none_to_default(tuple_k, ())
    tuple_k_last = none_to_default(tuple_k_last, ())
    # check input
    error = []
    check_type(dict_i, "dict_i", (dict, array_wrapper, numpy__ndarray), error)
    check_list(function_to_compute, "function_to_compute", list_fct, error)
    check_type(function_arguments, "function_arguments", dict, error)
    check_type(dict_o, "dict_o", dict, error)
    check_type(tuple_k, "tuple_k", tuple, error)
    check_type(tuple_k_last, "tuple_k_last", tuple, error)
    print_fail(inspect__stack(), "\n".join(k for k in error))
    # loop through nested levels
    if isinstance(dict_i, dict):
        list_keys = list(dict_i.keys())
        list_keys = sorted(list_keys, key=lambda s: s.lower()) if isinstance(list_keys[0], str) else sorted(list_keys)
        for k in list_keys:
            dict_o, tuple_k, tuple_k_last = compute_1var(
                dict_i[k], function_to_compute, function_arguments=function_arguments, dict_o=dict_o,
                tuple_k=tuple_k + (k,), tuple_k_last=tuple_k_last + (list_keys[-1],))
    else:
        # compute given function
        dict_i = dict_i.load()
        arr_o = getattr(xarray_tools, function_to_compute)(dict_i, **function_arguments)
        # save values
        put_in_dict(dict_o, arr_o, *tuple_k)
        # remove relevant keys from the tuples of keys
        tuple_k, tuple_k_last = tuple_for_dict(tuple_k, tuple_k_last)
    return dict_o, tuple_k, tuple_k_last


def reader_json(
        data_diagnostics: list[str] = None,
        data_epoch_lengths: list[str] = None,
        data_experiments: list[str] = None,
        data_observational_datasets: dict = None,
        data_projects: list[str] = None,
        smile_array: bool = True,
        smile_minimum_size: int = 1,
        smile_rejected: list[str] = None,
        smile_require_picontrol: bool = False,
        **kwargs) -> (dict, dict):
    """
    Read json dictionary and select values

    Input:
    ------
    :param data_diagnostics: list[str], optional
        Diagnostic names; e.g., data_diagnostics = ['var_pr_ano_n30e', 'var_ts_ano_n30e']
        Default is None (all available diagnostics read)
    :param data_epoch_lengths: list[str], optional
        Epoch length names; e.g., data_epoch_lengths = ['030_year_epoch', '150_year_epoch']
        Default is None (all available epoch lengths read)
    :param data_experiments: list[str], optional
        Experiment names; e.g., data_experiments = ['historical', 'piControl']
        Default is None (all available experiments read)
    :param data_observational_datasets: dict[str, list[str]], optional
        Read only desired observational datasets for each diagnostic;
        e.g., data_observational_datasets = {'var_pr_ano_n30e': ['CMAP'], 'var_ts_ano_n30e': ['OISSTv2']}
        Default is None (all available observational datasets read)
    :param data_projects: list[str], optional
        Project names; e.g., data_projects = ['cmip6', 'observations']
        Default is None (all available projects read)
    :param smile_array: bool, optional
        True to create a xarray.DataArray with dimensions [time, member]; e.g., smile_array = True
        Default is True (xarray.DataArray created)
    :param smile_minimum_size: int, optional
        Minimum number of members (epochs for piControl) required to keep a SMILE; e.g., smile_minimum_size = 10
        Default is 1 (all SMILEs are read)
    :param smile_rejected: list[str], optional
        Reject given SMILEs; e.g., smile_rejected = ['KACE-1-0-G']
        Default is None (all SMILEs are read)
    :param smile_require_picontrol: bool, optional
        True to keep SMILEs for which piControl is available (smile_minimum_size threshold must be met);
        e.g., smile_require_picontrol = True
        Default is False (all SMILEs are read)
    **kwargs - Discarded

    Output:
    -------
    :return dict_diagnostics: dict
        Dictionary with five nested levels [diagnostic, epoch_length, project, experiment, dataset], filled with an
        2D array [time, member] of diagnostic values
    :return dict_metadata: dict
        Dictionary with two nested levels [diagnostic, metadata], filled with a string, metadata keys are: 'method',
        'name_long', 'name_short' and 'units'
    """
    # set None input to its default value
    data_diagnostics = none_to_default(data_diagnostics, [])
    data_epoch_lengths = none_to_default(data_epoch_lengths, [])
    data_experiments = none_to_default(data_experiments, [])
    data_observational_datasets = none_to_default(data_observational_datasets, {})
    data_projects = none_to_default(data_projects, [])
    smile_rejected = none_to_default(smile_rejected, [])
    # check input
    error = []
    check_type(data_diagnostics, "data_diagnostics", list, error)
    check_type(data_epoch_lengths, "data_epoch_lengths", list, error)
    check_type(data_experiments, "data_experiments", list, error)
    check_type(data_observational_datasets, "data_observational_datasets", dict, error)
    check_type(data_projects, "data_projects", list, error)
    check_type(smile_rejected, "smile_rejected", list, error)
    print_fail(inspect__stack(), "\n".join(k for k in error))
    # read input json file
    dict_i = bt.read_json()
    # output metadata and value dictionaries
    dict_diagnostics, dict_metadata = {}, {}
    # list desired diagnostics that are available
    list_dia = list(dict_i.keys())
    if len(data_diagnostics) > 0:
        list_dia = [k for k in data_diagnostics if k in list_dia]
    for dia in list_dia:  # loop on diagnostics
        # read metadata
        method = dict_i[dia]["metadata"]["method"]
        if method.split(";; ")[-1].split(") ")[-1] == "average computed":
            method = method.replace(") average computed", ") temporal average computed")
        name_long = dict_i[dia]["metadata"]["diagnostic_long_name"]
        if name_long[:11] == "average of ":
            name_long = name_long.replace("average of ", "temporal average of ")
        name_short = bt.short_name(dict_i[dia]["metadata"]["diagnostic_short_name"])
        name_units = bt.units(dict_i[dia]["metadata"]["units"])
        region = flatten_dict({"region": dict_i[dia]["metadata"]["region"]})
        variable = dict_i[dia]["metadata"]["variable"]
        for var in list(variable.keys()):
            if var.lower() == "pr":
                variable[var]["units"] = "mm.day**-1"
            elif var.lower() in ["lhf", "lwr", "nhf", "nhfo", "shf", "swr", "thf"]:
                variable[var]["units"] = "W.m**-2"
            elif var.lower() in ["slp", "taum", "taux", "tauy"]:
                variable[var]["units"] = "Pa"
            elif var.lower() == "ssh":
                variable[var]["units"] = "m"
            elif var.lower() in ["sst", "tos", "ts"]:
                variable[var]["units"] = "degC"
        variable = flatten_dict({"variable": variable})
        dt = {"method": method, "name_long": name_long, "name_short": name_short, "units": name_units}
        dt, _, _ = sort_dict({**dt, **region, **variable})
        dict_metadata[dia] = dt
        # list projects
        list_pro = list(dict_i[dia]["diagnostic"]["value"].keys())
        if len(data_projects) > 0:
            list_pro = [k for k in data_projects if k in list_pro]
        for pro in list_pro:  # loop on projects
            # list datasets
            list_dat = list(dict_i[dia]["diagnostic"]["value"][pro].keys())
            # select datasets
            if pro != "observations" and len(smile_rejected) > 0:
                list_dat = list(set(list_dat) - set(smile_rejected))
            elif pro == "observations" and dia in list(data_observational_datasets.keys()):
                # list desired observational datasets that are available
                list_dat = [k for k in list_dat if k in data_observational_datasets[dia]]
            list_dat = sorted(list_dat, key=lambda s: s.lower())
            for dat in list_dat:  # loop on datasets
                # temporary dictionary
                dict_t = dict_i[dia]["diagnostic"]["value"][pro][dat]
                # list desired experiments that are available
                list_ex1 = list(dict_t.keys())
                if len(data_experiments) > 0:
                    list_ex1 = [k for k in data_experiments if k in list_ex1]
                for ex1 in list_ex1:  # loop on experiments (main)
                    # 'rcp' and 'ssp' experiments must be added to historical time series
                    # list_ex2 = ["historical"] + [ex1] if ex1[:3] in ["rcp", "ssp"] else [ex1]
                    list_ex2 = [ex1]
                    # list members
                    list_mem = list(dict_t[list_ex2[0]].keys())
                    for ex2 in list_ex2[1:]:
                        list_mem = list(set(list_mem) & set(list(dict_t[ex2].keys())))
                    list_mem = bt.sort_members(dat, list_mem)
                    dict_v = {}
                    for ii, mem in enumerate(list_mem):  # loop on members
                        # list desired epoch lengths that are available
                        list_dur = list(dict_t[ex1][mem].keys())
                        if len(data_epoch_lengths) > 0:
                            list_dur = [k for k in data_epoch_lengths if k in list_dur]
                        for dur in list_dur:  # loop on epoch lengths
                            for ex2 in list_ex2:  # loop on experiments (secondary)
                                # list epochs
                                list_epo = sorted(list(dict_t[ex2][mem][dur].keys()), key=lambda s: s.lower())
                                # The HadGEM3-GC31-LL's piControl simulation has an increase of GMST between year 400 &
                                # 600
                                # All historical members start during the first 500 years of simulations, so only the
                                # first 500 years are kept
                                if ex2 == "piControl" and dat == "HadGEM3-GC31-LL":
                                    list_epo = [k for k in list_epo
                                                if int(k[1:]) <= int(list_epo[0][1:]) + 500 - int(dur.split("_")[0])]
                                for epo in list_epo:  # loop on epochs
                                    epo_new = int(epo[1:]) + int(dur.split("_")[0])
                                    mem_new = copy__deepcopy(mem)
                                    if isinstance(smile_array, bool) and smile_array and ex2 == "piControl":
                                        # piControl members must be created
                                        mem_new = "r" + str(ii * 10000 + epo_new) + "i1p1f1"
                                        # for piControl experiments, epochs are used as if they were members. so there
                                        # is only one epoch
                                        epo_new = int(dur.split("_")[0])
                                        if ii > 0:
                                            continue
                                    put_in_dict(dict_v, dict_t[ex2][mem][dur][epo], dur, mem_new, epo_new)
                    # create a 2D array [time, member] with diagnostic values
                    list_dur = sorted(list(dict_v.keys()), key=lambda s: s.lower())
                    for dur in list_dur:
                        # check if enough members are available
                        ld = dict_v[dur]
                        if pro not in ["cmip", "observations"] and (
                                (ex1 != "piControl" and len(list(ld.keys())) < smile_minimum_size) or
                                (ex1 == "piControl" and (
                                        (len(list(ld.keys())) == 1 and
                                         len(list(ld[list(ld.keys())[0]].keys())) < smile_minimum_size) or
                                        (1 < len(list(ld.keys())) < smile_minimum_size and
                                         len(list(ld[list(ld.keys())[0]].keys())) < 2)))):
                            continue
                        if isinstance(smile_array, bool) and smile_array:
                            # create array
                            arr = create_array_time_series(dict_v[dur], attrs=dict_metadata[dia], data_var=dia)
                            # save value
                            put_in_dict(dict_diagnostics, arr, dia, dur, pro, ex1, dat)
                        else:
                            for mem, list_yv in dict_v[dur].items():
                                # create array
                                arr = create_array_time_series(list_yv, attrs=dict_metadata[dia], data_var=dia)
                                # save value
                                put_in_dict(dict_diagnostics, arr, dia, dur, pro, ex1, dat, mem)
            # delete SMILEs with few members and/or SMILEs for which piControl is no available
            if pro != "observations" and dia in list(dict_diagnostics.keys()) and smile_require_picontrol is True:
                for dur in list(dict_diagnostics[dia].keys()):
                    if pro not in list(dict_diagnostics[dia][dur].keys()):
                        continue
                    # list datasets
                    list_datasets = list(set([dat for exp in list(dict_diagnostics[dia][dur][pro].keys())
                                              for dat in list(dict_diagnostics[dia][dur][pro][exp].keys())]))
                    for dat in list_datasets:
                        # list experiments
                        list_experiments = [exp for exp in list(dict_diagnostics[dia][dur][pro].keys())
                                            if dat in list(dict_diagnostics[dia][dur][pro][exp].keys())]
                        # delete dataset
                        if "piControl" not in list_experiments:
                            for exp in list_experiments:
                                del dict_diagnostics[dia][dur][pro][exp][dat]
                    # delete experiment if there is no dataset available
                    for exp in list(dict_diagnostics[dia][dur][pro].keys()):
                        if len(list(dict_diagnostics[dia][dur][pro][exp].keys())) == 0:
                            del dict_diagnostics[dia][dur][pro][exp]
                    # delete project if there is no experiment available
                    if len(dict_diagnostics[dia][dur][pro].keys()) == 0:
                        del dict_diagnostics[dia][dur][pro]
        # delete diagnostic or epoch length if needed
        if dia in list(dict_diagnostics.keys()):
            for dur in list(dict_diagnostics[dia].keys()):
                if len(dict_diagnostics[dia][dur].keys()) == 0:
                    del dict_diagnostics[dia][dur]
            # delete diagnostic if there is no epoch length available
            if len(dict_diagnostics[dia].keys()) == 0:
                del dict_diagnostics[dia]
    return dict_diagnostics, dict_metadata


def reader_netcdf(
        data_datasets: list[str] = None,
        data_diagnostics: list[str] = None,
        data_epoch_lengths: list[str] = None,
        data_experiments: list[str] = None,
        data_observational_datasets: dict = None,
        data_projects: list[str] = None,
        data_region: dict[str, tuple] = None,
        smile_minimum_size: int = 1,
        smile_rejected: list[str] = None,
        smile_require_picontrol: bool = False,
        **kwargs) -> (dict, dict):
    # set None input to its default value
    data_datasets = none_to_default(data_datasets, ["*"])
    data_diagnostics = none_to_default(data_diagnostics, ["*"])
    data_epoch_lengths = none_to_default(data_epoch_lengths, ["*"])
    data_experiments = none_to_default(data_experiments, ["*"])
    data_observational_datasets = none_to_default(data_observational_datasets, {})
    data_projects = none_to_default(data_projects, ["*"])
    smile_rejected = none_to_default(smile_rejected, [])
    # check input
    error = []
    check_type(data_datasets, "data_datasets", list, error)
    check_type(data_diagnostics, "data_diagnostics", list, error)
    check_type(data_epoch_lengths, "data_epoch_lengths", list, error)
    check_type(data_experiments, "data_experiments", list, error)
    check_type(data_observational_datasets, "data_observational_datasets", dict, error)
    check_type(data_projects, "data_projects", list, error)
    check_type(smile_minimum_size, "smile_minimum_size", int, error)
    check_type(smile_rejected, "smile_rejected", list, error)
    check_type(smile_require_picontrol, "smile_require_picontrol", bool, error)
    print_fail(inspect__stack(), "\n".join(k for k in error))
    # plot directory (relative to current file directory)
    data_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/data"
    # get data
    dict_diagnostics, dict_metadata = {}, {}
    for dia in data_diagnostics:
        for pro in data_projects:  # loop on projects
            # list datasets
            list_dat = copy__deepcopy(data_datasets)
            # select datasets
            if pro != "observations" and list_dat != ["*"] and len(smile_rejected) > 0:
                list_dat = list(set(list_dat) - set(smile_rejected))
            elif pro == "observations" and dia in list(data_observational_datasets.keys()):
                # list desired observational datasets that are available
                list_dat = data_observational_datasets[dia]
            list_dat = sorted(list_dat, key=lambda s: s.lower())
            for dat in list_dat:  # loop on datasets
                for exp in data_experiments:  # loop on experiments
                    for dur in data_epoch_lengths:  # loop on epoch lengths
                        # construct filename
                        # filenames = "_".join([pro.lower(), dat, exp, dur, "extremes_full"]) + ".nc"
                        filenames = "_".join([pro.lower(), dat, exp, dur, dia]) + ".nc"
                        filenames = os__path__join(data_directory, filenames)
                        filenames = sorted(list(glob__iglob(filenames)), key=lambda s: s.lower())
                        for fil in filenames:  # loop on files
                            # get names
                            project, dataset, experiment = fil.split("/")[-1].split("_")[:3]
                            list_ens = ["i1", "i1p1", "i2", "i2p2", "p1", "p1f1", "p1f2", "p1f3", "p2", "p3", "p3f1",
                                        "p4", "p5f1", "f1", "f2"]
                            if experiment in list_ens:
                                dataset += "_" + str(experiment)
                                experiment = fil.split(str(project) + "_" + str(dataset) + "_")[-1].split("_")[0]
                            epoch_length = fil.split("_" + str(dataset) + "_" + str(experiment) + "_")[-1]
                            epoch_length = epoch_length.split("_" + str(dia))[0]
                            if dataset in smile_rejected:
                                continue
                            # read dataset
                            arr = netcdf_open(fil, dia, data_region)
                            if arr is None:
                                continue
                            # save value
                            for kk, da in arr.items():
                                if da is None or (pro not in ["cmip", "observations"] and isinstance(da, array_wrapper)
                                                  and da.shape[0] < smile_minimum_size):
                                    continue
                                put_in_dict(dict_diagnostics, da, kk, epoch_length, project, experiment, dataset)
    # delete SMILEs with few members and/or SMILEs for which piControl is no available
    if smile_require_picontrol:
        for dia in list(dict_diagnostics.keys()):
            for dur in list(dict_diagnostics[dia].keys()):
                for pro in list(dict_diagnostics[dia][dur].keys()):
                    if pro == "observations":
                        continue
                    # list all datasets
                    list_dat = [dat for exp in list(dict_diagnostics[dia][dur][pro].keys())
                                for dat in list(dict_diagnostics[dia][dur][pro][exp].keys())]
                    list_dat = sorted(list(set(list_dat)), key=lambda s: s.lower())
                    for dat in list_dat:
                        # list experiments for given dataset
                        list_exp = [exp for exp in list(dict_diagnostics[dia][dur][pro].keys())
                                    if dat in list(dict_diagnostics[dia][dur][pro][exp].keys())]
                        # check if piControl is available
                        if "piControl" not in list_exp:
                            # delete given dataset
                            for exp in list_exp:
                                del dict_diagnostics[dia][dur][pro][exp][dat]
                    # check if all datasets have been deleted
                    for exp in list(dict_diagnostics[dia][dur][pro].keys()):
                        if len(list(dict_diagnostics[dia][dur][pro][exp].keys())) == 0:
                            del dict_diagnostics[dia][dur][pro][exp]
                # check if all experiments have been deleted
                for pro in list(dict_diagnostics[dia][dur].keys()):
                    if len(list(dict_diagnostics[dia][dur][pro].keys())) == 0:
                        del dict_diagnostics[dia][dur][pro]
            # check if all epoch lengths have been deleted
            for dur in list(dict_diagnostics[dia].keys()):
                if len(list(dict_diagnostics[dia][dur].keys())) == 0:
                    del dict_diagnostics[dia][dur]
        # check if all projects have been deleted
        for dia in list(dict_diagnostics.keys()):
            if len(list(dict_diagnostics[dia].keys())) == 0:
                del dict_diagnostics[dia]
    # read metadata
    for dia in list(dict_diagnostics.keys()):
        if dia in list(dict_metadata.keys()):
            continue
        for dur in list(dict_diagnostics[dia].keys()):
            if dia in list(dict_metadata.keys()):
                continue
            for pro in list(dict_diagnostics[dia][dur].keys()):
                if dia in list(dict_metadata.keys()):
                    break
                for exp in list(dict_diagnostics[dia][dur][pro].keys()):
                    if dia in list(dict_metadata.keys()):
                        break
                    for dat in list(dict_diagnostics[dia][dur][pro][exp].keys()):
                        if dia in list(dict_metadata.keys()):
                            break
                        # get array
                        da = dict_diagnostics[dia][dur][pro][exp][dat]
                        # read metadata
                        method = da.attrs["method"] if "method" in list(da.attrs.keys()) else ""
                        if " read during " in method:
                            method = method.split(" read during ")[0][3:] + ": 1) read in Global [90S-90N ; 0-360E]" + \
                                ";; ".join([""] + method.split(";; ")[1:])
                        name_long = da.attrs["long_name"] if "long_name" in list(da.attrs.keys()) else ""
                        if "name_long" in list(da.attrs.keys()) and name_long == "":
                            name_long = da.attrs["name_long"]
                        if name_long[:11] == "average of ":
                            name_long = name_long.replace("average of ", "temporal average of ")
                        name_short = ""
                        if "short_name" in list(da.attrs.keys()):
                            name_short = bt.short_name(da.attrs["short_name"]).replace("pr ", "PR ")
                            name_short = name_short.replace("sst ", "SST ")
                        if "name_short" in list(da.attrs.keys()):
                            name_short = da.attrs["name_short"]
                        name_units = bt.units(da.attrs["units"]) if "units" in list(da.attrs.keys()) else ""
                        # time_period = da.attrs["time_period"]
                        # dict_metadata[dia] = {"method": method, "name_long": name_long, "name_short": name_short,
                        #                       "time_period": time_period, "units": name_units}
                        dict_metadata[dia] = {"method": method, "name_long": name_long, "name_short": name_short,
                                              "units": name_units}
                        for k1, v1 in da.attrs.items():
                            if "__" in k1 and k1.split("__")[0] in ["region", "variable"] and \
                                    k1 not in ["variable__pr__standard_name", "variable__sst__standard_name"]:
                                # modify key
                                k2 = "__".join(["name_long" if i2 == "long_name" else (
                                    "name_short" if i2 == "short_name" else
                                    ("cf_name_standard" if i2 == "cf_standard_name" else i2)) for i2 in k1.split("__")])
                                if k2.split("__")[0] in ["region"] and len(k2.split("__")) == 3:
                                    k2 = "__".join([i2.upper() if i1 == 1 else i2
                                                    for i1, i2 in enumerate(k2.split("__"))])
                                # modify value
                                v2 = copy__deepcopy(v1)
                                if "__latitude" in k2 or "__longitude" in k2:
                                    v2 = v2.replace("(", "[").replace(")", "]")
                                elif k2 == "variable__sst__name_short":
                                    v2 = v2.upper()
                                dict_metadata[dia][k2] = v2
                        if "region__GLOB__name_long" in list(dict_metadata[dia].keys()) and \
                                "region__GLOB__name_short" not in list(dict_metadata[dia].keys()):
                            dict_metadata[dia]["region__GLOB__name_short"] = "GLOB"
                        for kk in ["pr", "lhf", "lwr", "nhf", "nhfo", "shf", "slp", "ssh", "sst", "swr", "taum", "taux",
                                   "tauy", "thf", "tos", "ts"]:
                            k1, k2 = "variable__" + str(kk) + "__name_short", "variable__" + str(kk) + "__units"
                            if k1 in list(dict_metadata[dia].keys()) and k2 not in list(dict_metadata[dia].keys()):
                                if kk == "pr":
                                    dict_metadata[dia][k2] = "mm.day**-1"
                                elif kk in ["lhf", "lwr", "nhf", "nhfo", "shf", "swr", "thf"]:
                                    dict_metadata[dia][k2] = "W.m**-2"
                                elif kk in ["slp", "taum", "taux", "tauy"]:
                                    dict_metadata[dia][k2] = "Pa"
                                elif kk == "ssh":
                                    dict_metadata[dia][k2] = "m"
                                elif kk in ["sst", "tos", "ts"]:
                                    dict_metadata[dia][k2] = "degC"
    return dict_diagnostics, dict_metadata
# ---------------------------------------------------------------------------------------------------------------------#
