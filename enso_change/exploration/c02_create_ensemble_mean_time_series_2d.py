# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# NetCDF: takes json data and create ensemble mean time series
# 2D netCDF files are used as input
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from inspect import stack as inspect__stack
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
import sys
from typing import Literal

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, print_fail, sort_dict
from enso_change.wrapper.processors import compute_1var, reader_json, reader_netcdf
from enso_change.wrapper.xarray_base import numpy__array, concatenate, to_netcdf
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["ave_pr_val_n34e", "ave_tx_val_n34e", "ave_ts_val_n34e"],
    "data_diagnostics2": ["ave_taux_val"],
    # list of epoch lengths: list[str]
    # "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # "data_experiments": ["piControl", "historical"],
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    # "data_projects": ["cmip5", "observations"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    "smile_rejected": ["EC-EARTH"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 100.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 40,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {},
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {},
    # titles
    "fig_titles": {},
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {},
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def c02_create_ensemble_mean_time_series_2d_plot(
        data_diagnostics1: list = None,
        data_diagnostics2: list = None,
        data_epoch_lengths: list = None,
        data_experiments: list = None,
        data_observational_datasets: dict = None,
        data_projects: list = None,
        smile_minimum_size: int = None,
        smile_rejected: list = None,
        smile_require_picontrol: bool = None,
        proc_ci: float = None,
        proc_co: dict[str: dict[str, list[float]]] = None,
        proc_md: Literal["members", "percent"] = None,
        proc_ps: int = None,
        proc_th: int = None,
        fig_datasets: dict[str, dict[str, dict[str, dict[str, str]]]] = None,
        fig_format: Literal["eps", "pdf", "png", "svg"] = None,
        fig_colors: dict = None,
        fig_linestyles: dict = None,
        fig_markers: dict = None,
        fig_panel_size: dict[str, dict[str, float | int]] = None,
        fig_ticks: dict[
            str, list[float | int] | str | None | dict[
                str, list[float | int] | str | None | dict[
                    str, list[float | int] | str | None | dict[
                        str, list[float | int] | str | None | dict[
                            str, list[float | int] | str | None]]]]] = None,
        fig_titles: dict[str, dict[str, str] | str] = None,
        panel_param: dict = None,
        **kwargs):
    #
    # -- Set to default
    #
    data_diagnostics1 = none_to_default(data_diagnostics1, default["data_diagnostics1"])
    data_diagnostics2 = none_to_default(data_diagnostics2, default["data_diagnostics2"])
    data_epoch_lengths = none_to_default(data_epoch_lengths, default["data_epoch_lengths"])
    data_experiments = none_to_default(data_experiments, default["data_experiments"])
    data_observational_datasets = none_to_default(data_observational_datasets, default["data_observational_datasets"])
    data_projects = none_to_default(data_projects, default["data_projects"])
    smile_minimum_size = none_to_default(smile_minimum_size, default["smile_minimum_size"])
    smile_rejected = none_to_default(smile_rejected, default["smile_rejected"])
    smile_require_picontrol = none_to_default(smile_require_picontrol, default["smile_require_picontrol"])
    proc_ci = none_to_default(proc_ci, default["proc_ci"])
    proc_co = none_to_default(proc_co, default["proc_co"])
    proc_md = none_to_default(proc_md, default["proc_md"])
    proc_ps = none_to_default(proc_ps, default["proc_ps"])
    proc_th = none_to_default(proc_th, default["proc_th"])
    fig_colors = none_to_default(fig_colors, default["fig_colors"])
    fig_datasets = none_to_default(fig_datasets, default["fig_datasets"])
    fig_format = none_to_default(fig_format, default["fig_format"])
    fig_linestyles = none_to_default(fig_linestyles, default["fig_linestyles"])
    fig_markers = none_to_default(fig_markers, default["fig_markers"])
    fig_panel_size = none_to_default(fig_panel_size, default["fig_panel_size"])
    fig_ticks = none_to_default(fig_ticks, default["fig_ticks"])
    fig_titles = none_to_default(fig_titles, default["fig_titles"])
    panel_param = none_to_default(panel_param, default["panel_param"])
    #
    # -- Read data
    #
    # values_json = {"diagnostic": {"epoch_length": {"project": {"experiment": {"dataset": DataArray}}}}}
    values_json, metadata_json = reader_json(
        data_diagnostics=data_diagnostics1, data_epoch_lengths=data_epoch_lengths, data_experiments=data_experiments,
        data_observational_datasets=data_observational_datasets, data_projects=data_projects,
        smile_minimum_size=smile_minimum_size, smile_rejected=smile_rejected,
        smile_require_picontrol=smile_require_picontrol)
    print("reader_json", len(list(values_json.keys())), list(values_json.keys()))
    # values_netcdf = {"diagnostic": {"epoch_length": {"project": {"experiment": {"dataset": DataArray}}}}}
    values_netcdf, metadata_netcdf = reader_netcdf(
        data_diagnostics=data_diagnostics2, data_epoch_lengths=data_epoch_lengths, data_experiments=data_experiments,
        data_observational_datasets=data_observational_datasets, data_projects=data_projects,
        smile_minimum_size=smile_minimum_size, smile_rejected=smile_rejected,
        smile_require_picontrol=smile_require_picontrol)
    print("reader_netcdf", len(list(values_netcdf.keys())), list(values_netcdf.keys()))
    # concatenate historical and SSPs
    print("concatenate historical and SSPs")
    for dd in [values_json, values_netcdf]:
        for dia in list(dd.keys()):
            for dur in list(dd[dia].keys()):
                for pro in list(dd[dia][dur].keys()):
                    if pro == "observations":
                        continue
                    # dictionary
                    d1 = dd[dia][dur][pro]
                    # list datasets
                    list_dataset = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                    for dat in list_dataset:
                        # list experiments
                        list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                        list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                        # delete datasets
                        if "historical" not in list_exp:
                            for exp in list_exp:
                                del dd[dia][dur][pro][exp][dat]
                            continue
                        # concatenate historical and SSPs
                        for exp in list_ssp:
                            # arrays
                            da_h = copy__deepcopy(dd[dia][dur][pro]["historical"][dat])
                            da_s = copy__deepcopy(dd[dia][dur][pro][exp][dat])
                            del dd[dia][dur][pro][exp][dat]
                            # list members in both historical and ssp
                            lm = sorted(list(set(da_h["member"].values.tolist()) & set(da_s["member"].values.tolist())))
                            # select members
                            da_h, da_s = da_h.sel(indexers={"member": lm}), da_s.sel(indexers={"member": lm})
                            # concatenate along year axis
                            da_o = concatenate((da_h, da_s), "year")
                            # resave
                            put_in_dict(dd, copy__deepcopy(da_o), dia, dur, pro, exp, dat)
                    # delete experiment if no dataset available
                    for exp in list(dd[dia][dur][pro].keys()):
                        if len(list(dd[dia][dur][pro][exp].keys())) == 0:
                            del dd[dia][dur][pro][exp]
    n_ensemble = {}
    for dia in list(values_netcdf.keys()):
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                if pro == "observations":
                    continue
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_netcdf[dia][dur][pro][exp].keys()), pro, exp)
    for pro in data_projects:
        for exp in data_experiments:
            if not check_dict_levels(n_ensemble, pro, exp):
                continue
            list_ensembles = copy__deepcopy(n_ensemble[pro][exp])
            del n_ensemble[pro][exp]
            list_ensembles = sorted(list(set(list_ensembles)), key=lambda s: s.lower())
            put_in_dict(n_ensemble, copy__deepcopy(list_ensembles), pro, exp)
            print(pro, exp.ljust(10), "number of ensembles: " + str(len(list_ensembles)).rjust(2))
            # if exp == "historical":
            #     print(list_ensembles)
    list_datasets, list_numbers = [], []
    for pro in data_projects:
        if not check_dict_levels(n_ensemble, pro):
            continue
        # list datasets for given project
        list_dat = sorted(list(set([dat for exp in list(n_ensemble[pro].keys()) for dat in n_ensemble[pro][exp]])),
                          key=lambda s: s.lower())
        list_datasets += copy__deepcopy(list_dat)
        list_numbers += [k + 1 + (100 if pro == "cmip6" else 0) for k in range(len(list_dat))]
    dict_datasets = {k: v for k, v in zip(list_datasets, list_numbers)}
    # select pr/ts members for some models (members where not selected for the availability of both for 2D data)
    # dia2 = data_diagnostics1[0]
    # for dia1 in list(values_netcdf.keys()):
    #     for dur in list(values_netcdf[dia1].keys()):
    #         for pro in list(values_netcdf[dia1][dur].keys()):
    #             for exp in list(values_netcdf[dia1][dur][pro].keys()):
    #                 if exp not in ["historical", "rcp26", "rcp45", "rcp85", "ssp245"]:
    #                     continue
    #                 for dat in list(values_netcdf[dia1][dur][pro][exp].keys()):
    #                     # datasets manually checked
    #                     if not ((pro == "cmip5" and ((exp in ["historical", "rcp45"] and dat == "MPI-ESM") or
    #                                             (exp in ["rcp26", "rcp45", "rcp85"] and dat == "CCSM4") or
    #                                             (exp in ["rcp45"] and dat == "GFDL-CM3"))) or
    #                             (pro == "cmip6" and exp == "ssp245" and dat == "UKESM1-0-LL")) or \
    #                             not check_dict_levels(values_netcdf, dia1, dur, pro, exp, dat) or \
    #                             not check_dict_levels(values_json, dia2, dur, pro, exp, dat):
    #                         continue
    #                     # arrays
    #                     da1, da2 = values_netcdf[dia1][dur][pro][exp][dat], values_json[dia2][dur][pro][exp][dat]
    #                     print(dia1, dur, pro, exp, dat)
    #                     del values_netcdf[dia1][dur][pro][exp][dat]
    #                     # select common members
    #                     list_mem = list(set(da1["member"].values.tolist()) & set(da2["member"].values.tolist()))
    #                     da1 = da1.sel(indexers={"member": list_mem})
    #                     # resave
    #                     put_in_dict(values_netcdf, da1, dia1, dur, pro, exp, dat)
    # check ensemble size of datasets
    for dur in data_epoch_lengths:
        for pro in data_projects:
            for exp in data_experiments:
                # list datasets
                list_dat = []
                for dia1 in data_diagnostics1:
                    if check_dict_levels(values_json, dia1, dur, pro, exp):
                        list_dat += list(values_json[dia1][dur][pro][exp].keys())
                for dia1 in data_diagnostics2:
                    if check_dict_levels(values_netcdf, dia1, dur, pro, exp):
                        list_dat += list(values_netcdf[dia1][dur][pro][exp].keys())
                list_dat = sorted(list(set(list_dat)), key=lambda s: s.lower())
                for dat in list_dat:
                    # list diagnostics available for given dataset
                    list_diagnostics = []
                    for dia1 in data_diagnostics2:
                        if check_dict_levels(values_netcdf, dia1, dur, pro, exp, dat):
                            list_diagnostics.append(dia1)
                    missing = sorted(list(set(data_diagnostics2) - set(list_diagnostics)), key=lambda s: s.lower())
                    if len(missing) != 0:
                        error = "missing diagnostics\n" + ", ".join([dur, pro, exp.ljust(10), dat]) + "\nmissing: " + \
                                ", ".join(missing)
                        print_fail(inspect__stack(), error)
                    print(dur, pro, exp, dat)
                    # ensemble size
                    dt1 = {dia1.replace("_tx_", "_taux_").replace("_n34e", ""):
                               values_json[dia1][dur][pro][exp][dat]["member"].size for dia1 in data_diagnostics1}
                    dt2 = {dia1: values_netcdf[dia1][dur][pro][exp][dat]["member"].size for dia1 in data_diagnostics2}
                    maximum_ensemble_size = max(list(dt1.values()))
                    missing = [str(dia1) + " (" + str(dt2[dia1]) + " vs. " + str(dt1[dia1]) + ")"
                               for dia1 in data_diagnostics2 if (exp != "piControl" and dt2[dia1] != dt1[dia1]) or
                               (exp == "piControl" and dt2[dia1] != maximum_ensemble_size and dt2[dia1] < 40)]
                    if len(missing) != 0:
                        error = "missing members\n" + ", ".join([dur, pro, exp, dat]) + "\nmissing: " + \
                                ", ".join(missing)
                        print_fail(inspect__stack(), error)
                    # check members
                    for dia1 in data_diagnostics2:
                        dia2 = dia1.replace("_taux_", "_tx_") + "_n34e"
                        if not check_dict_levels(values_netcdf, dia1, dur, pro, exp, dat) or \
                                not check_dict_levels(values_json, dia2, dur, pro, exp, dat):
                            continue
                        # arrays
                        da1, da2 = values_netcdf[dia1][dur][pro][exp][dat], values_json[dia2][dur][pro][exp][dat]
                        # select common members
                        l1, l2 = da1["member"].values.tolist(), da2["member"].values.tolist()
                        list_mem1, list_mem2 = list(set(l1) - set(l2)), list(set(l2) - set(l1))
                        if len(list_mem1) != 0 or len(list_mem2) != 0:
                            error = "members non matching\n" + ", ".join([dia1, dur, pro, exp, dat])
                            if len(list_mem1) != 0:
                                error += "\nin json not in netCDF: " + ", ".join(list_mem1)
                            if len(list_mem2) != 0:
                                error += "\nin netCDF not in json: " + ", ".join(list_mem2)
                            print_fail(inspect__stack(), error)
    #
    # -- Processing
    #
    # mask low temperature for sst (as I use the ts variable sometimes cold points, i.e., land or ice points are
    # included)
    for dia in list(values_netcdf.keys()):
        if dia != "ave_sst_val":
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        del values_netcdf[dia][dur][pro][exp][dat]
                        # mask cold values
                        da = da.where(da >= -2)
                        # resave
                        put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # ensemble mean
    mean_netcdf, _, _ = compute_1var(values_netcdf, "statistic",
                             function_arguments={"dim": "member", "stat": "mea"})
    # regroup experiments
    print("regroup experiments")
    mean_netcdf_exp = {}
    for dia in list(mean_netcdf.keys()):
        for dur in list(mean_netcdf[dia].keys()):
            for pro in list(mean_netcdf[dia][dur].keys()):
                for exp in list(mean_netcdf[dia][dur][pro].keys()):
                    for dat in list(mean_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da1 = values_netcdf[dia][dur][pro][exp][dat]
                        da2 = mean_netcdf[dia][dur][pro][exp][dat]
                        if "statistics" in list(da2.dims):
                            da2 = da2.isel(indexers={"statistics": 0})
                        ensemble_size = da1["member"].size
                        da2.attrs.update({"ensemble_size": str(dat) + "=" + str(ensemble_size)})
                        # save
                        exp_o = "scenario" + str(exp[-2:]) if exp[:3] in ["rcp", "ssp"] else copy__deepcopy(exp)
                        pro_o = copy__deepcopy(pro[:4])
                        # for size in [1, 2, 3, 4, 5, 10]:
                        # for size in [1, 2, 3, 4, 5]:
                        for size in [1]:
                            if ensemble_size >= size:
                                size_o = str(size).zfill(3)
                                put_in_dict(mean_netcdf_exp, da2, dia, dur, size_o, pro_o, exp_o, dat)
    print("regroup experiments: done", len(list(mean_netcdf_exp.keys())))
    # select datasets for piControl and another experiment
    # output plot directory (relative to current file directory)
    data_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/data"
    print("select datasets for piControl and another experiment")
    for dia in list(mean_netcdf_exp.keys()):
        for dur in list(mean_netcdf_exp[dia].keys()):
            for siz in list(mean_netcdf_exp[dia][dur].keys()):
                for pro in list(mean_netcdf_exp[dia][dur][siz].keys()):
                    print(dia, dur, siz, pro)
                    # -- create multiple piControl experiments corresponding to models in each experiment
                    # list experiments
                    list_experiments = [k for k in list(mean_netcdf_exp[dia][dur][siz][pro].keys()) if k != "piControl"]
                    dt = {}
                    for exp in list_experiments:
                        for dat, da in mean_netcdf_exp[dia][dur][siz][pro][exp].items():
                            put_in_dict(dt, da, exp, dat)
                    # -- regroup datasets in a single array
                    for exp in list(dt.keys()):
                        # concatenate
                        list_dat = [dat for dat in list(dict_datasets.keys()) if dat in list(dt[exp].keys())]
                        da = concatenate(tuple(dt[exp][dat] for dat in list_dat), "dataset")
                        da = da.assign_coords({"dataset": numpy__array([dict_datasets[dat] for dat in list_dat])})
                        # attributes
                        dat = list_dat[0]
                        att = str(dict_datasets[dat]).zfill(3) + " " + str(dt[exp][dat].attrs["ensemble_size"])
                        for dat in list_dat[1:]:
                            att += "; " + str(dict_datasets[dat]).zfill(3) + " " + str(
                                dt[exp][dat].attrs["ensemble_size"])
                        att, _, _ = sort_dict({"dataset": att, **metadata_netcdf[dia]})
                        if "ensemble_size" in list(att.keys()):
                            del att["ensemble_size"]
                        # from json import dumps as json__dumps
                        # print(json__dumps(att, indent=4))
                        # stop
                        da.attrs.update(att)
                        # save
                        dia2 = copy__deepcopy(dia)
                        if dia in ["ave_nhf_val", "ave_taux_val", "ave_tauy_val", "tre_sst_val"]:
                            dia2 = dia2.replace("_nhf_", "_hf_").replace("_sst_", "_ts_")
                            dia2 = dia2.replace("_taux_", "_tx_").replace("_tauy_", "_ty_")
                            da = da.rename(dia2)
                        dat_o = "multi" + str("model" if pro == "cmip" else "obs") + str(siz)
                        print(pro, dat_o, str(exp).ljust(22), dur, dia2, str(len(list_dat)).rjust(3))
                        file_o = "_".join([pro, dat_o, exp, dur, dia2]) + ".nc"
                        to_netcdf(da, os__path__join(data_directory, file_o))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
c02_create_ensemble_mean_time_series_2d_plot()
