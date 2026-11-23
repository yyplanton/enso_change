# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Uses PR bias and zonal gradient change to explain ENSO change
# Figure:
#  - ???
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
import sys
from typing import Literal

# numpy
from numpy.random import randint as numpy__random__randint
# scipy
from scipy.stats import percentileofscore as scipy__stats__percentileofscore

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import pandas, default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, correlation_pearson, median, linear_regression, percentile, \
    linear_regression_slope
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, mean_weighted, polyval
from enso_change.wrapper.xarray_tools import create_land_sea_mask, select_region
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": [],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_glob", "ave_ts_val"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "scenario26", "scenario45", "scenario85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    # "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    "data_projects": ["cmip"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    # "smile_rejected": ["CMCC-CM", "EC-EARTH", "KACE-1-0-G"],
    # "smile_rejected": ["EC-EARTH", "KACE-1-0-G"],
    "smile_rejected": ["EC-EARTH", "MPI-ESM-LR"],
    # "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 60.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co":{
        "std_pr_ano_n34e": {"hs09--1": 85, "hs09--2": 0.61, "hs09--3": 68},
        "std_ts_ano_n34e": {"hs09--1": 97, "hs09--2": 0.66, "hs09--3": 69}
    },
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": {"ave": 4, "std": 2, "ske": 2},
    # threshold for selection
    "proc_se": 50,  # 100 / 3,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.005, "y": 0.005},
        "panel_01": {"x_delt": 300, "x_size": 1440, "y_delt": 150, "y_size": 480},
        "panel_11": {"x_delt": 940, "x_size": 800, "y_delt": 300, "y_size": 800},
        "panel_21": {"x_delt": 320, "x_size": 800, "y_delt": 300, "y_size": 800},
        "panel_31": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_32": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        # "panel_01": {"x_delt": 300, "x_size": 1440, "y_delt": 150, "y_size": 480},
        # "panel_11": {"x_delt": 300, "x_size": 800, "y_delt": 300, "y_size": 800},
        # "panel_21": {"x_delt": 700, "x_size": 800, "y_delt": 300, "y_size": 800},
        # "panel_31": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        # "panel_32": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "cmip5": "green",
        "cmip6": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "8",
        "cmip6": "D",
        "scenario26": "8",
        "scenario45": "D",
        "scenario85": "s",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {},
    # titles
    "fig_titles": {
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 30,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e73_multi_linear_regression_plot(
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
        proc_se: float = None,
        proc_th: dict[str: int] = None,
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
    proc_se = none_to_default(proc_se, default["proc_se"])
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
    print("reader_json", str(len(list(values_json.keys()))).rjust(3))
    # values_netcdf = {"diagnostic": {"epoch_length": {"project": {"experiment": {"dataset": DataArray}}}}}
    values_netcdf, metadata_netcdf = reader_netcdf(
        data_diagnostics=data_diagnostics2, data_epoch_lengths=data_epoch_lengths, data_experiments=data_experiments,
        data_observational_datasets=data_observational_datasets, data_projects=data_projects,
        data_region={"X": (20, 380), "Y": (-90, 90)}, smile_minimum_size=smile_minimum_size,
        smile_rejected=smile_rejected, smile_require_picontrol=smile_require_picontrol)
    print("reader_netcdf", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # remove datasets
    dia1, dur1, pro1 = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat1, exp1 = "multimodel001", "historical"
    da = values_netcdf[dia1][dur1][pro1][exp1][dat1]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    to_remove = [dict_datasets[k] for k in smile_rejected if k in list(dict_datasets.keys())]
    if len(to_remove) > 0:
        for dia in list(values_netcdf.keys()):
            for dur in list(values_netcdf[dia].keys()):
                for pro in list(values_netcdf[dia][dur].keys()):
                    if pro != "cmip":
                        continue
                    for exp in list(values_netcdf[dia][dur][pro].keys()):
                        for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                            # array
                            da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                            del values_netcdf[dia][dur][pro][exp][dat]
                            # list models by number
                            list_mod = da["dataset"].values.tolist()
                            # remove model numbers
                            list_mod = sorted(list(set(list_mod) - set(to_remove)))
                            # select models
                            da = da.sel(indexers={"dataset": list_mod})
                            # save
                            put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # number of datasets
    n_ensemble = {}
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 4:
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                if pro != "cmip":
                    continue
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat, da in values_netcdf[dia][dur][pro][exp].items():
                        n = da.isel(year=0).notnull().sum().values if exp == "historical" else da["dataset"].size
                        put_in_dict(n_ensemble, n, exp, dat)
    for exp in data_experiments:
        if not check_dict_levels(n_ensemble, exp):
            continue
        for dat in list(n_ensemble[exp].keys()):
            print(exp.ljust(10), dat, "number of datasets: " + str(n_ensemble[exp][dat]).rjust(3))
    # print("I need to find the best regions for mean state changes across CMIP")
    # stop
    #
    # -- Processing
    #
    # -- compute rSST
    dia2, dia3 = "ave_ts_val_tron", data_diagnostics2[-1]
    if dia2 in list(values_netcdf.keys()) and dia3 in list(values_netcdf.keys()):
        dia_o = dia3.replace("_sst_", "_rsst_").replace("_ts_", "_rts_")
        for dur in list(values_netcdf[dia3].keys()):
            for pro in list(values_netcdf[dia3][dur].keys()):
                for exp in list(values_netcdf[dia3][dur][pro].keys()):
                    for dat in list(values_netcdf[dia3][dur][pro][exp].keys()):
                        if not check_dict_levels(values_netcdf, dia2, dur, pro, exp, dat):
                            continue
                        # arrays
                        da1 = copy__deepcopy(values_netcdf[dia3][dur][pro][exp][dat])
                        da2 = copy__deepcopy(values_netcdf[dia2][dur][pro][exp][dat])
                        # select common years
                        list_y = sorted(list(set(da1["year"].values.tolist()) & set(da2["year"].values.tolist())))
                        da1, da2 = da1.sel(indexers={"year": list_y}), da2.sel(indexers={"year": list_y})
                        # remove tropical mean SST
                        da1 -= da2
                        # save
                        put_in_dict(values_netcdf, da1, dia_o, dur, pro, exp, dat)
        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia3])
        val = metadata_netcdf[dia3]["name_short"].replace("SST", "rSST").replace("ts", "rTS")
        metadata_netcdf[dia_o]["name_short"] = val
    for dia1 in [dia2]:  # [dia2, dia3]:
        if dia1 in list(values_netcdf.keys()):
            del values_netcdf[dia1]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- create land-sea mask
    dia2 = [k for k in data_diagnostics2 if len(k.split("_")) == 3][0]
    dur2 = data_epoch_lengths[0]
    pro2 = data_projects[0]
    exp2 = data_experiments[0]
    dat2 = list(values_netcdf[dia2][dur2][pro2][exp2].keys())[0]
    mask = create_land_sea_mask(values_netcdf[dia2][dur2][pro2][exp2][dat2])
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 3:
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        del values_netcdf[dia][dur][pro][exp][dat]
                        # mask land
                        da = da.where(mask == 0)
                        # save
                        put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # -- select bias year
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_bias = {}
    for dia1 in list(values_netcdf.keys()):
        if len(dia1.split("_")) != 3:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp in list(values_netcdf[dia1][dur1][pro].keys()):
                    for dat in list(values_netcdf[dia1][dur1][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia1][dur1][pro][exp][dat])
                        # select year
                        if pro == "cmip" or (pro == "obse" and dur1 == "050_year_epoch"):
                            yy = int(dur1.split("_")[0]) if exp == "piControl" else 2006
                        else:
                            yy = 2025
                        da = da.sel(indexers={"year": yy})
                        # save
                        put_in_dict(mean_bias, da, dur1, dia1, pro, exp, dat)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical"]:
                        continue
                    for dat in list(values_netcdf[dia1][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_netcdf, dia1, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia1, dur2, pro, exp2) or \
                                not check_dict_levels(values_netcdf, dia2, dur1, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_netcdf[dia1][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp2][dat])
                        da3 = copy__deepcopy(values_netcdf[dia2][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp1][dat])
                        else:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp2][dat])
                        # select experiment epoch
                        yy = 1850 + int(dur1.split("_")[0]) if exp1 == "historical" else 2100
                        da1, da3 = da1.sel(indexers={"year": yy}), da3.sel(indexers={"year": yy})
                        # select reference epoch
                        yy = 2006
                        da2, da4 = da2.sel(indexers={"year": yy}), da4.sel(indexers={"year": yy})
                        # change
                        da1 = da1 - da2 if exp1 in ["scenario26", "scenario45", "scenario85"] else da2 - da1
                        if dia1.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # degree of GMST
                        da3 = da3 - da4 if exp1 in ["scenario26", "scenario45", "scenario85"] else da4 - da3
                        # change per CC or not
                        l1 = ["raw", "pcc"] if exp1 in ["scenario26", "scenario45", "scenario85"] else ["raw"]
                        for war in l1:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if war == "pcc":
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- compute regional mean
    default_region = {
        "east_pla2": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EaEqP$_{Pla27}$"},
        "west_pla2": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "EaEqP$_{Pla27}$"},
        "peqe": {"latitude": [-6., 6.], "longitude": [215., 285.], "short_name": "EEP"},
        "peqw": {"latitude": [-6., 6.], "longitude": [145., 215.], "short_name": "WEP"},
        "pepr": {"latitude": [-8., 3.], "longitude": [190., 280.], "short_name": "PEQS"},
        # "pxts": {"regions": ["peqw", "peqe"], "short_name": "$\\Delta_{x}$EP"},
        "pxts": {"regions": ["west_pla2", "east_pla2"], "short_name": "Pla27"},
    }
    list_regions = ["pepr"]
    mean_bias_ave = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            if dia not in ["ave_pr_val"]:
                continue
            for pro in list(mean_bias[dur][dia].keys()):
                for exp in list(mean_bias[dur][dia][pro].keys()):
                    for dat in list(mean_bias[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias[dur][dia][pro][exp][dat]).compute()
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            if not (reg[-2:] == "pr" and "_pr_" in dia) and not (reg[-2:] == "ts" and "_ts_" in dia) \
                                    and not (reg[-2:] == "tx" and "_taux_" in dia):
                                continue
                            list_r2 = copy__deepcopy(default_region[reg]["regions"]) \
                                if reg in list(default_region.keys()) and \
                                   "regions" in list(default_region[reg].keys()) else [reg]
                            da_o, da_t = None, None
                            for r2 in list_r2:
                                # mask outside region
                                da_t = select_region(da, default_region[r2]["latitude"],
                                                     default_region[r2]["longitude"])
                                # horizontal average
                                da_t = mean_weighted(da_t, dim=[dim_y, dim_x], skipna=True)
                                # sum
                                da_o = copy__deepcopy(da_t) if da_o is None else da_o - da_t
                            if isinstance(da_o, array_wrapper) and isinstance(da_t, array_wrapper):
                                da_o.attrs.update(da_t.attrs)
                            # save
                            put_in_dict(mean_bias_ave, da_o, dur, str(dia) + "_" + str(reg), pro, exp, dat)
            if check_dict_levels(mean_bias_ave, dur):
                for dia_o in list(mean_bias_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_bias_ave", list(mean_bias_ave.keys()))
    list_regions = ["pxts"]
    mean_change_ave = {}
    war2 = "pcc"
    for dur1 in list(mean_change.keys()):
        for dia1 in list(mean_change[dur1].keys()):
            print("mean_change_ave", dur1, dia1)
            if dia1 not in ["ave_ts_val"]:
                continue
            for pro in list(mean_change[dur1][dia1].keys()):
                for exp1 in list(mean_change[dur1][dia1][pro].keys()):
                    for dat in list(mean_change[dur1][dia1][pro][exp1].keys()):
                        print("mean_change_ave", dur1, dia1, pro, exp1)
                        if not check_dict_levels(values_netcdf, dia1, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia1, dur2, pro, exp2) or \
                                not check_dict_levels(values_netcdf, dia2, dur1, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp2):
                            continue
                        # array
                        da = copy__deepcopy(mean_change[dur1][dia1][pro][exp1][dat][war2]).compute()
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            if not (reg[-2:] == "pr" and "_pr_" in dia1) and not (reg[-2:] == "ts" and "_ts_" in dia1) \
                                    and not (reg[-2:] == "tx" and "_tx_" in dia1):
                                continue
                            list_r2 = copy__deepcopy(default_region[reg]["regions"]) \
                                if reg in list(default_region.keys()) and \
                                   "regions" in list(default_region[reg].keys()) else [reg]
                            print("mean_change_ave", dur1, dia1, pro, exp1, reg, list_r2)
                            da_o, da_t = None, None
                            for r2 in list_r2:
                                # mask outside region
                                da_t = select_region(da, default_region[r2]["latitude"],
                                                     default_region[r2]["longitude"])
                                # horizontal average
                                da_t = mean_weighted(da_t, dim=[dim_y, dim_x], skipna=True)
                                print("mean_change_ave", dur1, dia1, pro, exp1, reg, r2, float(da_t.min().values),
                                      float(da_t.max().values))
                                # sum
                                da_o = copy__deepcopy(da_t) if da_o is None else da_o - da_t
                            if isinstance(da_o, array_wrapper) and isinstance(da_t, array_wrapper):
                                da_o.attrs.update(da_t.attrs)
                            print("mean_change_ave", dur1, dia1, pro, exp1, reg, float(da_o.min().values),
                                  float(da_o.max().values))
                            # save
                            put_in_dict(mean_change_ave, da_o, dur1, str(dia1) + "_" + str(reg), pro, exp1, dat,
                                        war2)
            if check_dict_levels(mean_bias_ave, dur1):
                for dia_o in list(mean_bias_ave[dur1].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia1])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_change_ave", list(mean_change_ave.keys()))
    # -- test
    dia_y = "std_ts_ano_n34e"
    war1 = "pcc"
    for dur in list(mean_change_ave.keys()):
        for dia_x in list(mean_change_ave[dur].keys()):
            for pro in list(mean_change_ave[dur][dia_x].keys()):
                for exp in list(mean_change_ave[dur][dia_x][pro].keys()):
                    for dat in list(mean_change_ave[dur][dia_x][pro][exp].keys()):
                        if not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat, war1):
                            continue
                        da_x = copy__deepcopy(mean_change_ave[dur][dia_x][pro][exp][dat][war1])
                        da_y = copy__deepcopy(mean_change[dur][dia_y][pro][exp][dat][war1])
                        corr = correlation_pearson(da_x, da_y, dim="dataset")
                        slop = linear_regression_slope(da_x, da_y, dim="dataset")
                        print(dur, dia_x, dia_y, pro, exp, dat,
                              "r=" + "{0:.3f}".format(round(float(corr.values), 3)).rjust(6),
                              "s=" + "{0:.3f}".format(round(float(slop.values), 3)).rjust(10))
    # stop
    # -- data for model
    # list dataset numbers in all scenarios
    list_dat_nbr = []
    for dur in list(mean_bias_ave.keys()):
        for dia in list(mean_bias_ave[dur].keys()):
            for pro in list(mean_bias_ave[dur][dia].keys()):
                for exp in list(mean_bias_ave[dur][dia][pro].keys()):
                    if exp[:3] not in ["sce"]:
                        continue
                    for da in mean_bias_ave[dur][dia][pro][exp].values():
                        list_dat_nbr += da["dataset"].values.tolist()
    list_dat_nbr = sorted(list(set(list_dat_nbr)))
    # format panda dataframe
    dat1, dur1, pro1, war1 = "multimodel001", data_epoch_lengths[0], data_projects[0], "pcc"
    dfs = []
    dict_var = {"ave_pr_val_pepr": mean_bias_ave[dur1],
                "ave_ts_val_pxts": mean_change_ave[dur1],
                "std_ts_ano_n34e": mean_change[dur1]}
    var1 = list(dict_var.keys())[0]
    list_exp = sorted([k for k in list(dict_var[var1][var1][pro1].keys()) if k[:3] in ["sce"]],
                      key=lambda s: s.lower())
    for exp in list_exp:
        df_scenario = None
        for var, dic in dict_var.items():
            # get xarray data array
            if var == "ave_pr_val_pepr":
                da = dic[var][pro1]["historical"][dat1]
                da = da.sel(indexers={"dataset": list_dat_nbr})
            else:
                da = dic[var][pro1][exp][dat1][war1]
            # convert to pandas Series
            series = da.to_series()
            series.name = "X1" if var == "ave_pr_val_pepr" else ("X2" if var == "ave_ts_val_pxts" else "Y")
            if df_scenario is None:
                df_scenario = series.to_frame()
            else:
                df_scenario = df_scenario.join(series, how="outer")
        # reset model index
        df_scenario = df_scenario.reset_index()
        # rename model identifier
        df_scenario = df_scenario.rename(columns={"dataset": "model"})
        df_scenario["scenario"] = exp
        dfs.append(df_scenario)
    # combine all scenarios
    df = pandas.concat(dfs, ignore_index=True)
    # remove incomplete observations
    df = df.dropna(subset=["model", "scenario", "X1", "X2", "Y"]).copy()
    # normalize
    for col in ["X1", "X2", "Y"]:
        df[col + "_z"] = (df[col] - df[col].mean()) / df[col].std()
    # calculate the within - between decomposition
    # model-specific average X2
    df["X2_between"] = df.groupby("model")["X2_z"].transform("mean")
    # deviation from model-specific average
    df["X2_within"] = df["X2_z"] - df["X2_between"]
    # fit regression models
    import statsmodels.formula.api as smf
    # model 1: X1 only
    m1 = smf.mixedlm("Y_z ~ X1_z + C(scenario)", data=df, groups=df["model"]).fit(reml=False)
    # model 2: X2 only
    m2 = smf.mixedlm("Y_z ~ X2_z + C(scenario)", data=df, groups=df["model"]).fit(reml=False)
    # model 3: X1 + X2
    m3 = smf.mixedlm("Y_z ~ X1_z + X2_z + C(scenario)", data=df, groups=df["model"]).fit(reml=False)
    # model 4: Within - between decomposition
    m4 = smf.mixedlm(
        "Y_z ~ X1_z + X2_between + X2_within + C(scenario)", data=df, groups=df["model"]).fit(reml=False)
    # model 5: Include interactions
    m5 = smf.mixedlm(
        "Y_z ~ X1_z * X2_between + X1_z * X2_within + C(scenario)", data=df, groups=df["model"]).fit(reml=False)
    # model 6: Scenario - dependent coefficients
    m6 = smf.mixedlm("Y_z ~ X1_z * C(scenario) + X2_between + X2_within * C(scenario)", data=df,
                     groups=df["model"]).fit(reml=False)
    # compare model performance
    models = {
        "M1: X1": m1,
        "M2: X2": m2,
        "M3: X1 + X2": m3,
        "M4: Within-between": m4,
        "M5: Interactions": m5,
        "M6: Scenario interactions": m6
    }
    comparison = pandas.DataFrame(
        {name: {"AIC": model.aic, "BIC": model.bic, "Log likelihood": model.llf} for name, model in models.items()}).T
    print(comparison.sort_values("AIC"))
    # Extract coefficients
    coef_table = pandas.DataFrame({name: model.fe_params for name, model in models.items()})
    print("coef_table")
    print(coef_table)
    print("conf_int")
    for k, d in models.items():
        print(k)
        print(d.conf_int())
    stop
    # -- compute mean maps
    mean_map = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            for pro1 in list(mean_bias[dur][dia].keys()):
                for exp in list(mean_bias[dur][dia][pro1].keys()):
                    if exp != "historical":
                        continue
                    for dat in list(mean_bias[dur][dia][pro1][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias[dur][dia][pro1][exp][dat])
                        # list models by number
                        list_mod = da["dataset"].values.tolist()
                        # save
                        pro_o = "obse" if pro1 == "obse" else "CMIP6 minus CMIP5"
                        tuple_k = (dur, exp, dat, dia, pro_o)
                        put_in_dict(mean_map, copy__deepcopy(da.mean(dim="dataset")), *tuple_k + ("l",))
                        if pro1 != "cmip":
                            continue
                        # separate cmip5/6
                        dt = {}
                        l2 = ["cmip5", "cmip6"] if pro1 == "cmip" else [pro1]
                        for pro2 in l2:
                            # model numbers corresponding to cmip5/6
                            l3 = [k for k in list_mod if (pro2 != "cmip6" and k < 100) or (pro2 == "cmip6" and k > 100)]
                            # select historical bias of models with projections
                            dt[pro2] = copy__deepcopy(da).sel(indexers={"dataset": l3})
                        # compare cmip5/6
                        da5, da6 = dt["cmip5"], dt["cmip6"]
                        # difference of multimodel mean
                        da_diff = da6.mean(dim="dataset") - da5.mean(dim="dataset")
                        # save
                        put_in_dict(mean_map, copy__deepcopy(da_diff), *tuple_k + ("s",))
    # -- region distributions
    mean_box = {}
    for dur in list(mean_bias_ave.keys()):
        for dia in list(mean_bias_ave[dur].keys()):
            for pro1 in list(mean_bias_ave[dur][dia].keys()):
                for exp in list(mean_bias_ave[dur][dia][pro1].keys()):
                    if exp != "historical":
                        continue
                    for dat in list(mean_bias_ave[dur][dia][pro1][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias_ave[dur][dia][pro1][exp][dat])
                        # list models by number
                        list_mod = da["dataset"].values.tolist()
                        # separate cmip5/6
                        l2 = ["cmip5", "cmip6"] if pro1 == "cmip" else [pro1]
                        for pro2 in l2:
                            # model numbers corresponding to cmip5/6
                            l3 = [k for k in list_mod if (pro2 != "cmip6" and k < 100) or (pro2 == "cmip6" and k > 100)]
                            # select historical bias of models with projections
                            da_o = copy__deepcopy(da).sel(indexers={"dataset": l3})
                            # save
                            put_in_dict(mean_box, copy__deepcopy(da_o), dur, exp, dat, dia, pro2)
    print("mean_box", list(mean_box.keys()))
    # -- select top / bottom models
    change_m, change_r = {}, {}
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    dia_x = None
    if check_dict_levels(mean_bias_ave, dur2):
        dia_x = list(mean_bias_ave[dur2].keys())[0]
    for dia_y in list(values_netcdf.keys()):
        if len(dia_y.split("_")) != 4 or dia_y.split("_")[0] not in ["std", "var"]:
            continue
        for dur1 in list(values_netcdf[dia_y].keys()):
            if dur1 != data_epoch_lengths[0]:
                continue
            for pro in list(values_netcdf[dia_y][dur1].keys()):
                if pro != "cmip":
                    continue
                for exp1 in list(values_netcdf[dia_y][dur1][pro].keys()):
                    if exp1 == "piControl":
                        continue
                    for dat in list(values_netcdf[dia_y][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_netcdf, dia_y, dur2, pro, exp1, dat) or \
                                not check_dict_levels(values_netcdf, dia_y, dur2, pro, exp2, dat) or \
                                not check_dict_levels(mean_bias_ave, dur2, dia_x, pro, exp1, dat):
                            continue
                        # -- compute change time series
                        # arrays
                        da1 = copy__deepcopy(values_netcdf[dia_y][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da2 = copy__deepcopy(values_netcdf[dia_y][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_netcdf[dia_y][dur2][pro][exp2][dat])
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 2006})
                        # change
                        da1 -= da2
                        if dia_y.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # -- select bottom / top
                        # array
                        da3 = copy__deepcopy(mean_bias_ave[dur2][dia_x][pro][exp1][dat])
                        print(dur2, dia_x, pro, exp1, dat, "da3", da3.shape, da3.dims)
                        ens_size = {int(k.split(" ")[0]): int(k.split("=")[1])
                                    for k in da3.attrs["dataset"].split("; ")}
                        # model nbr with dia_x value
                        dict_t = {int(k): float(v) for k, v in zip(da3["dataset"].values, da3.values)}
                        dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
                        list_datasets = list(dict_t.keys())
                        list_datasets = [k for k in list_datasets if k > 100]
                        nn = int(round(len(list_datasets) * proc_se / 100, 0))
                        list_d_bot, list_d_top = list_datasets[:nn], list_datasets[-nn:]
                        l1 = ["strongest", "weakest"] if "ave_ts_val" in dia_x else ["driest", "wettest"]
                        l2 = [list_d_top, list_d_bot] if "ave_ts_val" in dia_x else [list_d_bot, list_d_top]
                        for grp, l3 in zip(l1, l2):
                            # select models
                            da_t = da1.sel(indexers={"dataset": l3})
                            size = numpy__array([ens_size[k] for k in l3])
                            # -- remaining internal variability
                            # order of the polynomial
                            nn = proc_th[dia_y.split("_")[0]]
                            # 'forced change' is a Nth order polynomial fit computed on the ensemble mean
                            coef = da_t.polyfit("year", nn)
                            da_p = polyval(da_t["year"], coef["polyfit_coefficients"], degree_dim="degree")
                            # 'residual internal variability' is the temporal variance of the residual (da - da_p)
                            riv = (da_t - da_p).var(dim="year", ddof=1)
                            # correction: riv is underestimated as IV is underestimated for small ensembles
                            riv *= (1 / proc_co[dia_y]["hs09--" + str(nn)]) * size ** -0.2
                            # -- scaled distribution
                            # initial multimodel mean using ensemble mean -> will not change
                            da_m = da_t.mean(dim="dataset", skipna=True)
                            # initial multimodel spread using ensemble mean -> that what will be scaled
                            si = da_t.var(dim="dataset", skipna=True)
                            # initial multimodel mean residual internal variability -> will be used for scaling
                            ri = copy__deepcopy(riv).mean(dim="dataset", skipna=True)
                            # final spread (initial spread minus residual internal variability) -> must be set to 0 if
                            # negative
                            sf = (si - ri) / si
                            sf = sf.where(sf > 0, other=0)
                            # scaled distribution
                            da_r = da_m + (da_t - da_m) * sf
                            # confidence interval
                            da_r = percentile(da_r, dim="dataset", per=percentiles, skipna=True)
                            put_in_dict(change_m, da_m, dur1, dia_y, exp1, pro, grp)
                            put_in_dict(change_r, da_r, dur1, dia_y, exp1, pro, grp)
                        # da = da3.sel(indexers={"dataset": list_d_bot})
                        # print(proc_se, "min-max", "{0:.2f}".format(round(float(da.min().values), 2)).rjust(5),
                        #       "{0:.2f}".format(round(float(da.max().values), 2)).rjust(5), "mean",
                        #       "{0:.2f}".format(round(float(da.mean().values), 2)).rjust(5))
                        # da4 = copy__deepcopy(mean_bias_ave["045_year_epoch"][dia_x]["obse"]["historical"]["multiobs001"])
                        # print("045_year_epoch", dia_x, "obse", "historical", "multiobs001", da4.shape, da4.dims)
                        # print("min-max", "{0:.2f}".format(round(float(da4.min().values), 2)).rjust(5),
                        #       "{0:.2f}".format(round(float(da4.max().values), 2)).rjust(5), "mean",
                        #       "{0:.2f}".format(round(float(da4.mean().values), 2)).rjust(5))
                        # print("m+/-50%", "{0:.2f}".format(round(float((da4.mean() * 0.5).values), 2)).rjust(5),
                        #       "{0:.2f}".format(round(float((da4.mean() * 1.5).values), 2)).rjust(5))
                        # stop
    # best_param = {}
    # for dia_x in list(mean_bias_ave[dur2].keys()):
    #     if "ave_ts_val" in dia_x:
    #         da_mod = copy__deepcopy(mean_bias_ave["050_year_epoch"][dia_x]["cmip"]["historical"]["multimodel001"])
    #         da_obs = copy__deepcopy(mean_bias_ave["045_year_epoch"][dia_x]["obse"]["historical"]["multiobs001"])
    #         rank = numpy__array([scipy__stats__percentileofscore(da_mod.values, k) for k in da_obs.values])
    #         if median(rank) < 60:
    #             continue
    #         print(dia_x, rank)
    #     for dia_y in list(values_netcdf.keys()):
    #         if len(dia_y.split("_")) != 4 or dia_y.split("_")[0] not in ["std", "var"]:
    #             continue
    #         for dur1 in list(values_netcdf[dia_y].keys()):
    #             if dur1 != data_epoch_lengths[0]:
    #                 continue
    #             for pro in list(values_netcdf[dia_y][dur1].keys()):
    #                 for exp1 in list(values_netcdf[dia_y][dur1][pro].keys()):
    #                     if pro == "cmip" and exp1 in ["historical", "piControl"]:
    #                         continue
    #                     for dat in list(values_netcdf[dia_y][dur1][pro][exp1].keys()):
    #                         if not check_dict_levels(values_netcdf, dia_y, dur2, pro, exp1, dat) or \
    #                                 not check_dict_levels(values_netcdf, dia_y, dur2, pro, exp2, dat) or \
    #                                 not check_dict_levels(mean_bias_ave, dur2, dia_x, pro, exp1, dat):
    #                             continue
    #                         # -- compute change time series
    #                         # arrays
    #                         da1 = copy__deepcopy(values_netcdf[dia_y][dur1][pro][exp1][dat])
    #                         if exp1 in ["scenario26", "scenario45", "scenario85"]:
    #                             da2 = copy__deepcopy(values_netcdf[dia_y][dur2][pro][exp1][dat])
    #                         else:
    #                             da2 = copy__deepcopy(values_netcdf[dia_y][dur2][pro][exp2][dat])
    #                         # select reference epoch
    #                         da1 = da1.sel(indexers={"year": slice(2071, 2100)})
    #                         da2 = da2.sel(indexers={"year": 2006})
    #                         # change
    #                         da1 -= da2
    #                         if dia_y.split("_")[0] in ["std", "var"]:
    #                             da1 *= 100 / da2
    #                         # -- select bottom / top
    #                         # array
    #                         da3 = copy__deepcopy(mean_bias_ave[dur2][dia_x][pro][exp1][dat])
    #                         # model nbr with dia_x value
    #                         dict_t = {int(k): float(v) for k, v in zip(da3["dataset"].values, da3.values)}
    #                         dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
    #                         list_datasets = list(dict_t.keys())
    #                         list_datasets = [k for k in list_datasets if k > 100]
    #                         nn = int(round(len(list_datasets) * proc_se / 100, 0))
    #                         list_d_bot, list_d_top = list_datasets[:nn], list_datasets[-nn:]
    #                         dict_t = {}
    #                         l1 = ["strongest", "weakest"] if "ave_ts_val" in dia_x else ["driest", "wettest"]
    #                         l2 = [list_d_top, list_d_bot] if "ave_ts_val" in dia_x else [list_d_bot, list_d_top]
    #                         for grp, l3 in zip(l1, l2):
    #                             # select models
    #                             da_t = da1.sel(indexers={"dataset": l3})
    #                             # initial multimodel mean using ensemble mean -> will not change
    #                             dict_t[grp] = da_t.mean(dim="dataset", skipna=True)
    #                         # difference
    #                         diff = dict_t[l1[0]] - dict_t[l1[-1]]
    #                         # save
    #                         put_in_dict(best_param, [float(diff.mean().values)], dia_x)
    # best_param = {k: float(numpy__array(v).mean()) for k, v in best_param.items()}
    # best_param = dict(sorted(best_param.items(), key=lambda item: item[1]))
    # for dia_x in list(best_param.keys())[-20:]:
    #     da_mod = copy__deepcopy(mean_bias_ave["050_year_epoch"][dia_x]["cmip"]["historical"]["multimodel001"])
    #     da_obs = copy__deepcopy(mean_bias_ave["045_year_epoch"][dia_x]["obse"]["historical"]["multiobs001"])
    #     rank = sorted([scipy__stats__percentileofscore(da_mod.values, k) for k in da_obs.values])
    #     rank = float(numpy__array(rank[1:] if "ave_ts_val" in dia_x else copy__deepcopy(rank)).mean())
    #     print(dia_x, "{0:.3f}".format(round(best_param[dia_x], 3)).rjust(6), "obs rank",
    #           "{0:.1f}".format(round(rank, 1)).rjust(5))
    # stop
    print("change ts min-max")
    for dur in list(change_m.keys()):
        for dia in list(change_m[dur].keys()):
            tmp = []
            for exp in list(change_m[dur][dia].keys()):
                for pro in list(change_m[dur][dia][exp].keys()):
                    for da in change_m[dur][dia][exp][pro].values():
                        tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group, panel = "group_1", "panel_01"
    dur2, exp2 = "050_year_epoch", "historical"
    dat2 = list(mean_map[dur2][exp2].keys())[0]
    dur3, pro3 = "045_year_epoch", "obse"
    dat3 = list(mean_map[dur3][exp2].keys())[0]
    ny2 = int(dur2.split("_")[0])
    for jj, dia in enumerate(list(mean_map[dur2][exp2][dat2].keys())):
        for ii, pro in enumerate(list(mean_map[dur2][exp2][dat2][dia].keys())):
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
            tuple_k = (grp, pan)
            print(dur2, exp2, dat2, dia.ljust(12), pro, tuple_k)
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, k2, *tk + (k1,))
            # panel axes definition
            if check_dict_levels(fig_ticks, dia, panel):
                for k1, k2 in fig_ticks[dia][panel].items():
                    if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                        continue
                    tmp = copy__deepcopy(k2)
                    if isinstance(tmp, str) and k1 == "map_c_nam":
                        tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                        tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                        tmp = tmp.replace(" ()", "")
                        tmp = tmp.replace("TS", "SST")
                    put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                    if k1 == "y_tic" and jj != 0 and isinstance(tmp, list):
                        put_in_dict(figure_axes, ["" for _ in tmp], *tuple_k + ("y_lab",))
            dt = {}
            if check_dict_levels(figure_axes, *tuple_k):
                dt = read_from_dict(figure_axes, *tuple_k)
            # test possibility to plot text / legend
            do_legend = False
            if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                    dt["y_lim"] is not None:
                do_legend = True
            # -- title
            if ii == 0:
                val = copy__deepcopy(metadata_netcdf[dia]["name_short"]) + ": " + str(2006 - ny2) + "-2005"
                val += " (" + str(metadata_netcdf[dia]["units"]) + ")"
                val = val.replace(" ()", "")
                val = val.replace("TS", "SST")
                put_in_dict(plot_data, val, *tuple_k + ("title_col",))
            if jj == 0:
                put_in_dict(plot_data, pro, *tuple_k + ("title_row",))
            # -- map
            d1 = mean_map[dur2][exp2][dat2][dia][pro]
            da_s = d1["s"]
            plot_type = "map"
            # shading
            put_in_dict(plot_data, da_s, *tuple_k + (str(plot_type) + "_c",))
            put_in_dict(plot_data, "right", *tuple_k + ("legend_position",))
            # contours
            da_l = None
            contour_type = "l2"
            if check_dict_levels(mean_map, dur3, exp2, dat3, dia, pro3, "l"):
                da_l = mean_map[dur3][exp2][dat3][dia][pro3]["l"]
            elif check_dict_levels(mean_map, dur2, exp2, dat3, dia, pro3, "l"):
                da_l = mean_map[dur2][exp2][dat3][dia][pro3]["l"]
            put_in_dict(plot_data, da_l, *tuple_k + (str(plot_type) + "_" + str(contour_type),))
            dd = {
                "colors": ["k"] * len(dt["map_" + str(contour_type) + "_tic"]),
                "linestyles": ["-" if k >= 0 else "--" for k in dt["map_" + str(contour_type) + "_tic"]],
                "linewidths": [1 if k != 0 else 2 for k in dt["map_" + str(contour_type) + "_tic"]],
            }
            put_in_dict(plot_data, copy__deepcopy(dd), *tuple_k + ("kwargs_" + str(contour_type),))
            # -- text
            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
            plot_type = "text"
            if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # panel 11: distribution
    print("distribution")
    list_pro = list(set([pro for dia in list(mean_map[dur2][exp2][dat2].keys())
                         for pro in list(mean_map[dur2][exp2][dat2][dia].keys())]))
    group = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + len(list_pro))
    panel = "panel_11"
    if check_dict_levels(mean_box, dur2, exp2, dat2):
        # list_dia = [k for k in list(mean_box[dur2][exp2][dat2].keys())
        #             if k in ["ave_pr_val_n30e", "ave_pr_val_n30s", "ave_pr_val_peqs"]]
        list_dia = [k for k in list(mean_box[dur2][exp2][dat2].keys())
                    if k in ["ave_pr_val_peqspr", "ave_pr_val_p063pr", "ave_pr_val_p074pr", "ave_pr_val_p117pr",
                             "ave_pr_val_p120pr", "ave_ts_val_gycaits", "ave_ts_val_p043ts", "ave_ts_val_p068ts",
                             "ave_taux_val_peqntx", "ave_taux_val_peqstx"]]
        for ii, dia2 in enumerate(list_dia):
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + ii).zfill(2)
            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1])).zfill(2)
            tuple_k = (grp, pan)
            print(dur2, exp2, dat2, dia2.ljust(12), tuple_k)
            # dictionary
            d1 = mean_box[dur2][exp2][dat2][dia2]
            list_d = list(d1.keys())
            nn = len(list_d)
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, k2, *tk + (k1,))
            # panel axes definition
            if check_dict_levels(fig_ticks, dia2, panel):
                for k1, k2 in fig_ticks[dia2][panel].items():
                    if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                        continue
                    tmp = copy__deepcopy(k2)
                    if isinstance(tmp, str) and k1 == "y_nam":
                        tmp = tmp.replace("VALUE", metadata_netcdf[dia2]["name_short"])
                        tmp = tmp.replace("UNITS", metadata_netcdf[dia2]["units"])
                        tmp = tmp.replace(" ()", "")
                    elif tmp is None and k1 == "x_lim":
                        tmp = [-0.5, nn - 0.5]
                    elif tmp is None and k1 == "x_tic":
                        tmp = list(range(nn))
                    put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                    if isinstance(tmp, list) and k1 == "x_tic":
                        put_in_dict(figure_axes, ["" for _ in list_d], *tuple_k + ("x_lab",))
            dt = {}
            if check_dict_levels(figure_axes, *tuple_k):
                dt = read_from_dict(figure_axes, *tuple_k)
            # test possibility to plot text / legend
            do_legend = False
            if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                    dt["y_lim"] is not None:
                do_legend = True
            # -- title
            # -- marker
            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
            pro3 = "cmip"
            dict_t1 = {}
            if check_dict_levels(mean_bias_ave, dur2, dia2, pro3, exp2, dat2):
                da3 = copy__deepcopy(mean_bias_ave[dur2][dia2][pro3][exp2][dat2])
                # model nbr with dia_x value
                dict_t2 = {int(k): float(v) for k, v in zip(da3["dataset"].values, da3.values)}
                dict_t2 = {k: v for k, v in sorted(dict_t2.items(), key=lambda item: item[1])}
                list_datasets = list(dict_t2.keys())
                list_datasets = [k for k in list_datasets if k > 100]
                nn = int(round(len(list_datasets) * proc_se / 100, 0))
                list_d_bot, list_d_top = list_datasets[:nn], list_datasets[-nn:]
                l1 = ["strongest", "weakest"] if "ave_ts_val" in dia_x else ["driest", "wettest"]
                l2 = [list_d_top, list_d_bot] if "ave_ts_val" in dia_x else [list_d_bot, list_d_top]
                dict_t1 = {grp: da3.sel(indexers={"dataset": l3}) for grp, l3 in zip(l1, l2)}
            list_grp = list(dict_t1.keys())
            if isinstance(do_legend, bool) and do_legend and len(list_grp) > 0:
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                for grp in list_grp:
                    # lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                    # lm += [fig_markers[grp] if grp in list(fig_markers.keys()) else "o"]
                    # ls += [150]
                    # lw += [2] * 2
                    # lw += [0.5]
                    # lx += [(x1 + x2) / 2]
                    # ly += [dict_t1[grp].mean().values]
                    # lz += [3]
                    lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                    lm += [fig_markers[grp] if grp in list(fig_markers.keys()) else "o"]
                    ls += [150]
                    lw += [2] * 2
                    lw += [0.5]
                    lx += [1]
                    ly += [dict_t1[grp].mean().values]
                    lz += [3]
            plot_type = "mar"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend and len(list_grp) > 0:
                for grp in list_grp:
                    # lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"] * 2
                    # ls += [fig_linestyles[grp] if grp in list(fig_linestyles.keys()) else ":"] * 2
                    # lw += [2] * 2
                    # lx += [dt["x_lim"]] * 2
                    # ly += [[dict_t1[grp].min().values] * 2, [dict_t1[grp].max().values] * 2]
                    # lz += [2]
                    lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"] * 2
                    ls += [fig_linestyles[grp] if grp in list(fig_linestyles.keys()) else ":"] * 2
                    lw += [2] * 2
                    lx += [[0.5, max(dt["x_lim"])]] * 2
                    ly += [[dict_t1[grp].min().values] * 2, [dict_t1[grp].max().values] * 2]
                    lz += [2]
            dur3, exp3 = "045_year_epoch", "historical"
            dat3 = list(mean_box[dur3][exp3].keys())[0]
            pro3 = "obse"
            if isinstance(do_legend, bool) and do_legend and \
                    check_dict_levels(mean_box, dur3, exp3, dat3, dia2, pro3):
                da1 = copy__deepcopy(mean_box[dur3][exp3][dat3][dia2][pro3])
                for da2 in da1:
                    lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                    ls += [fig_linestyles[pro3] if pro3 in list(fig_linestyles.keys()) else ":"]
                    lw += [2]
                    lx += [dt["x_lim"]]
                    ly += [[da2.values] * 2]
                    lz += [5]
            plot_type = "cur"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- violin
            la, lbod_fc, lc, lf, lls, llw, lmea, lmea_lw = [], [], [], [], [], [], [], []
            lmed, lmed_lw, lq, lq_ls, lq_lw, lw, lx, ly, lz = [], [], [], [], [], [], [], [], []
            for kk, pro in enumerate(list_d):
                la += [0.4]
                lbod_fc += ["k"]
                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                lf += [True]
                lls += ["-"]
                llw += [2]
                lmea += [False]
                lmea_lw += [0]
                lmed += [True]
                lmed_lw += [4]
                lq += [[0.25, 0.75]]
                lq_ls += ["--"]
                lq_lw += [1]
                lw += [0.9]
                ly += [d1[pro].values]
                lx += [kk]
                lz += [3]
            plot_type = "vio"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
                put_in_dict(plot_data, copy__deepcopy(lbod_fc), *tuple_k + (str(plot_type) + "_bod_fc",))
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_flie",))
                put_in_dict(plot_data, copy__deepcopy(lls), *tuple_k + (str(plot_type) + "_ls",))
                put_in_dict(plot_data, copy__deepcopy(llw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lmea), *tuple_k + (str(plot_type) + "_mea",))
                put_in_dict(plot_data, copy__deepcopy(lmea_lw), *tuple_k + (str(plot_type) + "_mea_lw",))
                put_in_dict(plot_data, copy__deepcopy(lmed), *tuple_k + (str(plot_type) + "_med",))
                put_in_dict(plot_data, copy__deepcopy(lmed_lw), *tuple_k + (str(plot_type) + "_med_lw",))
                put_in_dict(plot_data, copy__deepcopy(lq), *tuple_k + (str(plot_type) + "_q",))
                put_in_dict(plot_data, copy__deepcopy(lq_ls), *tuple_k + (str(plot_type) + "_qua_ls",))
                put_in_dict(plot_data, copy__deepcopy(lq_lw), *tuple_k + (str(plot_type) + "_qua_lw",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- text
            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # x lab
                for kk, pro in enumerate(list_d):
                    lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                    lf += [15]
                    lh += ["center"]
                    lr += [0]
                    lt += [fig_titles[pro]]
                    lv += ["top"]
                    lw += ["regular"]
                    lx += [kk]
                    ly += [y1 - 5 * dy]
                # legend
                for grp in list_grp:
                    lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                    lf += [15]
                    lh += ["left"]
                    lr += [0]
                    val = str(int(proc_se)) + "%" if proc_se == int(proc_se) else (
                        "1/3" if proc_se == 100 / 3 else str(proc_se) + "%")
                    lt += [str(fig_titles[grp].replace("THR%", val)) + "\nof models"]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x2 + 10 * dx]
                    ly += [(dict_t1[grp].min().values + dict_t1[grp].max().values) / 2]
            plot_type = "text"
            if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # panel 21: scatter
    print("scatter")
    panel = "panel_21"
    dia_y, war2 = "std_ts_ano_n34e", "pcc"
    if check_dict_levels(mean_bias_ave, dur2) and check_dict_levels(mean_change, dur2, dia_y):
        # list_dia = [k for k in list(mean_bias_ave[dur2].keys())
        #             if k in ["ave_pr_val_n30e", "ave_pr_val_n30s", "ave_pr_val_peqs"]]
        list_dia = [k for k in list(mean_bias_ave[dur2].keys())
                    if k in ["ave_pr_val_peqspr", "ave_pr_val_p063pr", "ave_pr_val_p074pr", "ave_pr_val_p117pr",
                             "ave_pr_val_p120pr", "ave_ts_val_gycaits", "ave_ts_val_p043ts", "ave_ts_val_p068ts",
                             "ave_taux_val_peqntx", "ave_taux_val_peqstx"]]
        for ii, dia_x in enumerate(list_dia):
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + ii).zfill(2)
            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1])).zfill(2)
            tuple_k = (grp, pan)
            print(dur2, exp2, dat2, dia_x.ljust(12), tuple_k)
            # dictionary
            dix = mean_bias_ave[dur2][dia_x]
            diy = mean_change[dur2][dia_y]
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, k2, *tk + (k1,))
            # panel axes definition
            dia = "--".join([dia_x, dia_y, war2])
            if check_dict_levels(fig_ticks, dia, panel):
                for k1, k2 in fig_ticks[dia][panel].items():
                    if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                        continue
                    tmp = copy__deepcopy(k2)
                    if isinstance(tmp, str) and k1 in ["x_nam"]:
                        tmp = tmp.replace("VALUE", metadata_netcdf[dia_x]["name_short"])
                        tmp = tmp.replace("UNITS", metadata_netcdf[dia_x]["units"])
                        tmp = tmp.replace(" ()", "")
                    elif isinstance(tmp, str) and k1 in ["y_nam"]:
                        tmp = tmp.replace("VALUE", metadata_netcdf[dia_y]["name_short"])
                        tmp = tmp.replace("UNITS", metadata_netcdf[dia_y]["units"])
                        tmp = tmp.replace(" ()", "")
                    if ii > 0 and k1 == "y_nam":
                        tmp = ""
                    put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
            dt = {}
            if check_dict_levels(figure_axes, *tuple_k):
                dt = read_from_dict(figure_axes, *tuple_k)
            # test possibility to plot text / legend
            do_legend = False
            if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                    dt["y_lim"] is not None:
                do_legend = True
            # -- title
            # -- marker
            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
            for pro in list(dix.keys()):
                for exp in list(dix[pro].keys()):
                    if exp[:3] not in ["rcp", "sce", "ssp"]:
                        continue
                    for dat in list(dix[pro][exp].keys()):
                        if not check_dict_levels(diy, pro, exp, dat, war2):
                            continue
                        da_x, da_y = copy__deepcopy(dix[pro][exp][dat]), copy__deepcopy(diy[pro][exp][dat][war2])
                        # model numbers corresponding to cmip5/6
                        list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                             set(da_y["dataset"].values.tolist())))
                        # list_m = [k for k in list_m if k > 100]
                        da_x, da_y = da_x.sel(indexers={"dataset": list_m}), da_y.sel(indexers={"dataset": list_m})
                        for dd, xx, yy in zip(list_m, da_x.values, da_y.values):
                            pro2 = "cmip5" if dd < 100 else "cmip6"
                            lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                            lm += [fig_markers[exp] if exp in list(fig_markers.keys()) else "o"]
                            ls += [80]
                            lw += [0.5]
                            lx += [xx]
                            ly += [yy]
                            lz += [numpy__random__randint(2, 9)]
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                k1 = 0
                for pro in list(dix.keys()):
                    for exp in list(dix[pro].keys()):
                        if exp[:3] not in ["rcp", "sce", "ssp"]:
                            continue
                        for dat in list(dix[pro][exp].keys()):
                            if not check_dict_levels(diy, pro, exp, dat, war2):
                                continue
                            lc += ["k"]
                            lm += [fig_markers[exp] if exp in list(fig_markers.keys()) else "o"]
                            ls += [80]
                            lw += [0.5]
                            lx += [x2 + 5 * dx]
                            ly += [y2 - (k1 * 8 + 6) * dy]
                            k1 += 1
            plot_type = "mar"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                for pro in list(dix.keys()):
                    for exp in list(dix[pro].keys()):
                        if exp[:3] not in ["rcp", "sce", "ssp"]:
                            continue
                        for dat in list(dix[pro][exp].keys()):
                            if not check_dict_levels(diy, pro, exp, dat, war2):
                                continue
                            da_x, da_y = copy__deepcopy(dix[pro][exp][dat]), copy__deepcopy(diy[pro][exp][dat][war2])
                            list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                 set(da_y["dataset"].values.tolist())))
                            # list_m = [k for k in list_m if k > 100]
                            da_x, da_y = da_x.sel(indexers={"dataset": list_m}), da_y.sel(indexers={"dataset": list_m})
                            ln = linear_regression(da_x.values, da_y.values)
                            xx = [x1 + 5 * dx, x2 - 5 * dx]
                            lc += ["k"]
                            ls += ["-"]
                            lw += [2]
                            lx += [xx]
                            ly += [[ln["s"] * k + ln["i"] for k in xx]]
                            lz += [11]
            dur3, exp3 = "045_year_epoch", "historical"
            dat3 = list(mean_box[dur3][exp3].keys())[0]
            pro3 = "obse"
            if isinstance(do_legend, bool) and do_legend and \
                    check_dict_levels(mean_box, dur3, exp3, dat3, dia_x, pro3):
                da1 = copy__deepcopy(mean_box[dur3][exp3][dat3][dia_x][pro3])
                for da2 in da1:
                    lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                    ls += [fig_linestyles[pro3] if pro3 in list(fig_linestyles.keys()) else ":"]
                    lw += [2]
                    lx += [[da2.values] * 2]
                    ly += [dt["y_lim"]]
                    lz += [12]
            if isinstance(do_legend, bool) and do_legend:
                # x = 0 and y = 0
                for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                    if k1[0] < 0 < k1[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                        lz += [13]
            plot_type = "cur"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            if isinstance(do_legend, bool) and do_legend:
                # x = 0 and y = 0
                for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                    if k1[0] < 0 < k1[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                        lz += [13]
            # -- text
            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # linear regression
                k1 = 0
                for pro in list(dix.keys()):
                    for exp in list(dix[pro].keys()):
                        if exp[:3] not in ["rcp", "sce", "ssp"]:
                            continue
                        for dat in list(dix[pro][exp].keys()):
                            if not check_dict_levels(diy, pro, exp, dat, war2):
                                continue
                            da_x, da_y = copy__deepcopy(dix[pro][exp][dat]), copy__deepcopy(diy[pro][exp][dat][war2])
                            list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                 set(da_y["dataset"].values.tolist())))
                            # list_m = [k for k in list_m if k > 100]
                            da_x, da_y = da_x.sel(indexers={"dataset": list_m}), da_y.sel(indexers={"dataset": list_m})
                            ln = linear_regression(da_x.values, da_y.values)
                            lc += ["k"] * 4
                            lf += [15] * 4
                            # lh += ["right"] * 3 + ["left"]
                            # lr += [0] * 4
                            # val = "***" if ln["p"] < 1e-3 else ("**" if ln["p"] < 1e-2 else (
                            #     "*" if ln["p"] < 0.05 else ""))
                            # lt += [fig_titles[exp], "r=", "{0:.2f}".format(round(ln["r"], 2)), str(val)]
                            # lv += ["center"] * 4
                            # lw += ["regular"] * 4
                            # lx += [x2 - 42 * dx, x2 - 30 * dx, x2 - 12.5 * dx, x2 - 12 * dx]
                            # ly += [y2 - (k1 * 8 + 6) * dy] * 4
                            lh += ["left", "left", "right", "left"]
                            lr += [0] * 4
                            val = "***" if ln["p"] < 1e-3 else ("**" if ln["p"] < 1e-2 else (
                                "*" if ln["p"] < 0.05 else ""))
                            lt += [fig_titles[exp], "r=", "{0:.2f}".format(round(ln["r"], 2)), str(val)]
                            lv += ["center"] * 4
                            lw += ["regular"] * 4
                            lx += [x2 + 10 * dx, x2 + 36 * dx, x2 + 64 * dx, x2 + 64.5 * dx]
                            ly += [y2 - (k1 * 8 + 6) * dy] * 4
                            k1 += 1
                for k1, pro in enumerate(["cmip5", "cmip6"]):
                    lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                    lf += [15]
                    lh += ["right"]
                    lr += [0]
                    lt += [fig_titles[pro]]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x2 - 2 * dx]
                    ly += [y2 - (k1 * 8 + 6) * dy]
            plot_type = "text"
            if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # panel 31-32: time series
    print("time series")
    group = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + 1)
    dur2, pro2, war2 = "030_year_epoch", "cmip", "raw"
    if check_dict_levels(change_m, dur2):
        list_dia = [k for k in data_diagnostics2 if k in list(change_m[dur2].keys())]
        for ii, dia in enumerate(list_dia):
            list_exp = [k for k in data_experiments if k in list(change_m[dur2][dia].keys())]
            for jj, exp in enumerate(list_exp):
                tp = "panel_31" if exp == "historical" else "panel_32"
                pan = str(tp.split("_")[0]) + "_" + str(int(tp.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (grp, pan)
                print(dur2, dia, exp.ljust(10), tuple_k)
                # panel size
                tk = tuple_k[:-2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    for k1, d1 in fig_panel_size.items():
                        if k1 == tp and not check_dict_levels(figure_psize, *tk + (pan,)):
                            for k2, d2 in d1.items():
                                if check_dict_levels(figure_psize, *tk + (pan, k2,)):
                                    continue
                                tmp = copy__deepcopy(d2)
                                if k2 == "x_size" and check_dict_levels(fig_ticks, "time", tp, "lim"):
                                    k3 = read_from_dict(fig_ticks, "time", tp, "lim")
                                    tmp *= (k3[-1] - k3[0])
                                put_in_dict(figure_psize, tmp, *tk + (pan, k2,))
                        elif k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, copy__deepcopy(d1), *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time", dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dd, tp):
                        for k1, k2 in fig_ticks[dd][tp].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
                                tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and dia != list_dia[-1]):
                                tmp = ""
                            put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
                dt = {}
                if check_dict_levels(figure_axes, *tuple_k):
                    dt = read_from_dict(figure_axes, *tuple_k)
                # test possibility to plot text / legend
                do_legend = False
                if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                        dt["y_lim"] is not None:
                    do_legend = True
                # -- title
                if ii == 0:
                    put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                for pro in list(change_m[dur2][dia][exp].keys()):
                    for grp in list(change_m[dur2][dia][exp][pro].keys()):
                        da = copy__deepcopy(change_m[dur2][dia][exp][pro][grp])
                        lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                        ls += [fig_linestyles[grp] if grp in list(fig_linestyles.keys()) else ":"]
                        lw += [2]
                        lx += [da["year"].values]
                        ly += [da.values]
                        lz += [4]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # x = 0 and y = 0
                    for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                        if k1[0] < 0 < k1[1]:
                            lc += ["k"]
                            ls += [":"]
                            lw += [1]
                            lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                            ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                            lz += [10]
                    # legend
                    if jj == 0:
                        lc += ["k"]
                        ls += ["-"]
                        lw += [2]
                        lx += [[x1 + 10 * dx, x1 + 18 * dx]]
                        ly += [[y2 - (5) * dy] * 2]
                        lz += [1]
                plot_type = "cur"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- shading
                la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
                if check_dict_levels(change_r, dur2, dia, exp):
                    for pro in list(change_r[dur2][dia][exp].keys()):
                        for grp in list(change_r[dur2][dia][exp][pro].keys()):
                            da = copy__deepcopy(change_r[dur2][dia][exp][pro][grp])
                            la += [0.4]
                            lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                            lx += [da["year"].values]
                            ly1 += [da.min(dim="quantile").values]
                            ly2 += [da.max(dim="quantile").values]
                            lz += [3]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    if jj == 0:
                        la += [0.4]
                        lc += ["k"]
                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                        ly1 += [[y2 - (8 + 5 + 3) * dy] * 2]
                        ly2 += [[y2 - (8 + 5 - 3) * dy] * 2]
                        lz += [1]
                plot_type = "sha"
                if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
                    put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(plot_type) + "_y1",))
                    put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(plot_type) + "_y2",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- text
                lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    if jj == 0:
                        # legend
                        val = [str(int(k3)) if k3 == int(k3) else str(k3) for k3 in percentiles]
                        l1 = ["mean", "range"]
                        for k1, txt in enumerate(l1):
                            lc += ["k"]
                            lf += [15]
                            lh += ["left"]
                            lr += [0]
                            lt += [fig_titles[txt].replace("THR", "-".join(val))]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x1 + (20 if k1 == 0 else 12) * dx]
                            ly += [y2 - (k1 * 8 + 5) * dy]
                        # groups
                        list_grp = list(reversed(sorted(list(set(
                            [grp for pro in list(change_m[dur2][dia][exp].keys())
                             for grp in list(change_m[dur2][dia][exp][pro].keys())])), key=lambda s: s.lower())))
                        for k1, txt in enumerate(list_grp):
                            lc += [fig_colors[txt]]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            val = str(int(proc_se)) + "%" if proc_se == int(proc_se) else (
                                "1/3" if proc_se == 100 / 3 else str(proc_se) + "%")
                            lt += [fig_titles[txt].replace("THR%", val)]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x2 - 2 * dx]
                            ly += [y1 + (k1 * 7 + 4) * dy if dia[:10] == "std_ts_ano" else y2 - ((k1 + 2) * 8 + 5) * dy]
                plot_type = "text"
                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                    put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                    put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1, a4 = str(smile_minimum_size).zfill(2), str(len(smile_rejected)).zfill(2)
    a1 = int(proc_ci) if proc_ci == int(proc_ci) else copy__deepcopy(proc_ci)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a4) + "rem", dia_x])
    # figure loop
    if ' ':
        # list plot names
        list_names = list(plot_data.keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name)
        # plot
        fig_basic(plot_data, list_names, 1, figure_axes, fig_format, figure_name_and_path, figure_psize,
                  fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e73_multi_linear_regression_plot()
