# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Relationship between selected mean state change and ENSO SSTA change using only models with weak PR bias
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

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, scipy__stats__t, linear_regression
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, mean_weighted
from enso_change.wrapper.xarray_tools import select_region, create_land_sea_mask
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_tx_val", "ave_ts_val", "ave_ts_val_glob"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "045_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "scenario26", "scenario45", "scenario85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    # "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    "data_projects": ["cmip", "obse"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
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
        # no resample
        "ske_pr_ano_n30e": {"hs09--1": 0.77, "hs09--2": 0.63, "hs09--3": 54},
        "ske_pr_ano_n34e": {"hs09--1": 0.78, "hs09--2": 0.64, "hs09--3": 55},
        "ske_pr_ano_n40e": {"hs09--1": 0.76, "hs09--2": 0.63, "hs09--3": 53},
        "ske_ts_ano_n30e": {"hs09--1": 0.74, "hs09--2": 0.62, "hs09--3": 53},
        "ske_ts_ano_n34e": {"hs09--1": 0.74, "hs09--2": 0.62, "hs09--3": 53},
        "ske_ts_ano_n40e": {"hs09--1": 0.75, "hs09--2": 0.63, "hs09--3": 54},
        "std_pr_ano_n30e": {"hs09--1": 1.01, "hs09--2": 0.64, "hs09--3": 51},
        "std_pr_ano_n34e": {"hs09--1": 1.07, "hs09--2": 0.66, "hs09--3": 53},
        "std_pr_ano_n40e": {"hs09--1": 0.93, "hs09--2": 0.63, "hs09--3": 51},
        "std_ts_ano_n30e": {"hs09--1": 75, "hs09--2": 0.61, "hs09--3": 49},
        "std_ts_ano_n34e": {"hs09--1": 76, "hs09--2": 0.61, "hs09--3": 49},
        "std_ts_ano_n40e": {"hs09--1": 76, "hs09--2": 0.61, "hs09--3": 50},
        # resample 8000
        # "std_pr_ano_n34e": {"hs09--1": 1.11, "hs09--2": 0.68, "hs09--3": 54},
        # "std_ts_ano_n34e": {"hs09--1": 78, "hs09--2": 0.63, "hs09--3": 51}
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
        "panel_41": {"x_delt": 600, "x_size": 700, "y_delt": 200, "y_size": 700},
    },
    # colors
    "fig_colors": {
        "cmip": "k",  # "peru",  #
        "cmip5": "green",
        "cmip6": "r",
        "obse": "k",
        "driest": "sienna",
        "largest": "orange",
        "smallest": "darkviolet",
        "wettest": "dodgerblue",
        "scenario26": "grey",
        "scenario45": "grey",
        "scenario85": "grey",
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        "obse": "--",
        "driest": "-",
        "largest": "-",
        "smallest": "-",
        "wettest": "-",
    },
    # markers
    "fig_markers": {
        "cmip5": "H",
        "cmip6": "D",
        "scenario26": "H",
        "scenario45": "D",
        "scenario85": "s",
        "driest": "D",
        "largest": "D",
        "smallest": "D",
        "wettest": "D",
    },
    # ticks
    "fig_ticks": {
        "ave_pr_val_pepr--bias--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [0, 6],
                "x_nam": "PR bias (UNITS)",
                "x_nbr_minor": 1,
                "x_tic": list(range(0, 7, 1)),
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_pr_val_p2pr--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [0.0, 0.8],
                "x_nam": "PR change (UNITS)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(0, 9, 2)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_ts_val_pgts--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.4, 0.4],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_ts_val_p2ts--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.6, 0.6],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_ts_val_p3ts--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.6, 0.6],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_ts_val_p4ts--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.6, 0.6],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_ts_val_p5ts--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.6, 0.6],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_tx_val_p2tx--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.006, 0.006],
                "x_nam": "Taux change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 1000, 3) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_tx_val_p3tx--chan--std_ts_ano_n34e--chan": {
            "panel_41": {
                "x_lim": [-0.006, 0.006],
                "x_nam": "Taux change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 1000, 3) for k in range(-6, 7, 3)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
    },
    # titles
    "fig_titles": {
        "050": "models PR within 50% obs",
        "075": "models PR within 75% obs",
        "100": "models PR within 100% obs",
        "125": "models PR within 125% obs",
        "150": "models PR within 150% obs",
        "200": "models PR within 200% obs",
        "500": "all models",
        "ave_pr_val_pepr": "PR bias",
        "ave_pr_val_p2pr": "PR change",
        "ave_ts_val_pgts": "1) $\\Delta_{x}$SST change",
        "ave_ts_val_p2ts": "2) $\\Delta_{x}$SST change",
        "ave_ts_val_p3ts": "3) $\\Delta_{x}$SST change",
        "ave_ts_val_p4ts": "4) $\\Delta_{x}$SST change",
        "ave_ts_val_p5ts": "5) $\\Delta_{x}$SST change",
        "ave_tx_val_p2tx": "1) Taux change",
        "ave_tx_val_p3tx": "2) Taux change",
        "std_ts_ano_n34e": "$\\Delta_{t}$E relative to 1956-2005 (%)",
        "cmip": "CMIP 5&6",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        # figure title
        "figure_title": "Warming pattern linked to ENSO change when small PR bias",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 3,
        "title_row_x": 30,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e81_change_vs_change_bias_control_plot(
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
    for dia in list(metadata_netcdf.keys()):
        if "_ts_" in dia and "name_short" in list(metadata_netcdf[dia].keys()):
            metadata_netcdf[dia]["name_short"] = metadata_netcdf[dia]["name_short"].replace("TS", "SST")
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
    #
    # -- Processing
    #
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
                        put_in_dict(mean_bias, copy__deepcopy(da).compute(), dur1, dia1, pro, exp, dat)
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
                        # change per CC
                        da1 /= da3
                        # save
                        put_in_dict(mean_change, copy__deepcopy(da1).compute(), dur1, dia1, pro, exp1, dat)
    print("mean_change", list(mean_change.keys()))
    default_region = {
        # "p068n": {"latitude": [2., 8.], "longitude": [229., 279.], "short_name": "EEP$_{N}$"},
        # "p068s": {"latitude": [-4., 2.], "longitude": [229., 279.], "short_name": "EEP$_{S}$"},
        "pl2e": {"latitude": [-6., 6.], "longitude": [230., 270.], "short_name": "EEP"},
        "pl2w": {"latitude": [-5., 5.], "longitude": [135., 185.], "short_name": "WEP"},
        "pl3e": {"latitude": [-6., 6.], "longitude": [230., 270.], "short_name": "EEP"},
        "pl3w": {"latitude": [-4., 4.], "longitude": [130., 185.], "short_name": "WEP"},
        "pl4e": {"latitude": [-8., 8.], "longitude": [230., 270.], "short_name": "EEP"},
        "pl4w": {"latitude": [-5., 5.], "longitude": [130., 185.], "short_name": "WEP"},
        "pl5e": {"latitude": [-8., 8.], "longitude": [230., 270.], "short_name": "EEP"},
        "pl5w": {"latitude": [-6., 6.], "longitude": [135., 185.], "short_name": "WEP"},
        "plae": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EEP"},
        "plaw": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "WEP"},
        "pgxe": {"latitude": [-5., 5.], "longitude": [220., 280.], "short_name": "EEP"},
        "pgxw": {"latitude": [-5., 5.], "longitude": [140., 190.], "short_name": "WEP"},
        "pepr": {"latitude": [-8., 3.], "longitude": [190., 280.], "short_name": ""},
        "p2pr": {"latitude": [-14., 14.], "longitude": [235., 275.], "short_name": ""},
        "p2tx": {"latitude": [-8., 8.], "longitude": [180., 270.], "short_name": ""},
        "p3tx": {"latitude": [-8., 8.], "longitude": [175., 270.], "short_name": ""},
        # "pgts": {"regions": ["pgxw", "pgxe"], "short_name": "$\\Delta_{x}$EP"},
        "pgts": {"regions": ["plaw", "plae"], "short_name": "W - E"},
        "p2ts": {"regions": ["pl2w", "pl2e"], "short_name": "W - E"},
        "p3ts": {"regions": ["pl3w", "pl3e"], "short_name": "W - E"},
        "p4ts": {"regions": ["pl4w", "pl4e"], "short_name": "W - E"},
        "p5ts": {"regions": ["pl5w", "pl5e"], "short_name": "W - E"},
        # "p012pr": {"latitude": [-7., 3.], "longitude": [190., 280.], "short_name": "PEQS"},
        # "p120pr": {"latitude": [-8., 3.], "longitude": [190., 280.], "short_name": "PEQS"},
        # "pets": {"regions": ["p068n", "p068s"], "short_name": "$\\Delta_{y}$EEP"},
    }
    list_regions = ["pepr"]  # list(default_region.keys())
    mean_bias_ave = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            if dia not in ["ave_pr_val", "ave_ts_val"]:
                continue
            for pro in list(mean_bias[dur][dia].keys()):
                for exp in list(mean_bias[dur][dia][pro].keys()):
                    for dat in list(mean_bias[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias[dur][dia][pro][exp][dat]).compute()
                        # land-sea mask
                        da = da.where(create_land_sea_mask(da) == 0)
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            if not (reg[-2:] == "pr" and "_pr_" in dia) and not (reg[-2:] == "ts" and "_ts_" in dia):
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
                            put_in_dict(mean_bias_ave, da_o.compute(),
                                        dur, str(dia) + "_" + str(reg), pro, exp, dat)
            if check_dict_levels(mean_bias_ave, dur):
                for dia_o in list(mean_bias_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_bias_ave", list(mean_bias_ave.keys()))
    list_regions = ["pgts", "p2ts", "p3ts", "p4ts", "p5ts", "p2pr", "p2tx", "p3tx"]  # list(default_region.keys())
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if dia not in ["ave_pr_val", "ave_ts_val", "ave_tx_val"]:
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat]).compute()
                        # land-sea mask
                        da = da.where(create_land_sea_mask(da) == 0)
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            if not (reg[-2:] == "pr" and "_pr_" in dia) and not (reg[-2:] == "ts" and "_ts_" in dia) \
                                    and not (reg[-2:] == "tx" and "_tx_" in dia):
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
                            put_in_dict(mean_change_ave, da_o.compute(),
                                        dur, str(dia) + "_" + str(reg), pro, exp, dat)
            if check_dict_levels(mean_change_ave, dur):
                for dia_o in list(mean_change_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_change_ave", list(mean_change_ave.keys()))
    # -- model selection
    for dur in list(mean_bias_ave.keys()):
        if dur != data_epoch_lengths[-1]:
            continue
        # select in reference variable
        dia_r = "ave_pr_val_pepr"
        for dia in list(mean_bias_ave[dur].keys()):
            if dia != dia_r:
                continue
            for pro in list(mean_bias_ave[dur][dia].keys()):
                for exp in list(mean_bias_ave[dur][dia][pro].keys()):
                    for dat in list(mean_bias_ave[dur][dia][pro][exp].keys()):
                        dur3 = copy__deepcopy(dur) if "_ts_" in dia else "045_year_epoch"
                        dat3, exp3, pro3 = "multiobs001", "historical", "obse"
                        if not check_dict_levels(mean_bias_ave, dur3, dia, pro3, exp3, dat3):
                            continue
                        # array
                        da_m = copy__deepcopy(mean_bias_ave[dur][dia][pro][exp][dat])
                        da_r = copy__deepcopy(mean_bias_ave[dur3][dia][pro3][exp3][dat3])
                        del mean_bias_ave[dur][dia][pro][exp][dat]
                        # obs mean
                        da_r = float(da_r.mean().values)
                        # list models
                        list_models = da_m["dataset"].values.tolist()
                        # select
                        for thr in list(range(100, 151, 25)) + [500]:
                            da_t = abs((da_m - da_r) * 100 / da_r)
                            list_m = [k for k, v in zip(list_models, da_t.values) if v <= thr]
                            da_o = da_m.sel(indexers={"dataset": list_m})
                            put_in_dict(mean_bias_ave, da_o.compute(), dur, dia, pro, exp, dat, str(thr).zfill(3))
        # select for all variables
        for dia in list(mean_bias_ave[dur].keys()):
            if dia == dia_r:
                continue
            for pro in list(mean_bias_ave[dur][dia].keys()):
                for exp in list(mean_bias_ave[dur][dia][pro].keys()):
                    for dat in list(mean_bias_ave[dur][dia][pro][exp].keys()):
                        if not check_dict_levels(mean_bias_ave, dia_r, dia, pro, exp, dat):
                            continue
                        # array
                        da_m = copy__deepcopy(mean_bias_ave[dur][dia][pro][exp][dat])
                        del mean_bias_ave[dur][dia][pro][exp][dat]
                        for thr, da_r in mean_bias_ave[dur][dia_r][pro][exp][dat].items():
                            # select
                            list_models = sorted(list(set(da_m["dataset"].values.tolist()) &
                                                      set(da_r["dataset"].values.tolist())))
                            # save
                            da_o = da_m.sel(indexers={"dataset": list_models})
                            put_in_dict(mean_bias_ave, da_o.compute(), dur, dia, pro, exp, dat, thr)
        for dia in list(mean_change_ave[dur].keys()):
            for pro in list(mean_change_ave[dur][dia].keys()):
                for exp in list(mean_change_ave[dur][dia][pro].keys()):
                    for dat in list(mean_change_ave[dur][dia][pro][exp].keys()):
                        if not check_dict_levels(mean_bias_ave, dur, dia_r, pro, exp, dat):
                            continue
                        # array
                        da_m = copy__deepcopy(mean_change_ave[dur][dia][pro][exp][dat])
                        del mean_change_ave[dur][dia][pro][exp][dat]
                        for thr, da_r in mean_bias_ave[dur][dia_r][pro][exp][dat].items():
                            # select
                            list_models = sorted(list(set(da_m["dataset"].values.tolist()) &
                                                      set(da_r["dataset"].values.tolist())))
                            # save
                            da_o = da_m.sel(indexers={"dataset": list_models})
                            put_in_dict(mean_change_ave, da_o.compute(), dur, dia, pro, exp, dat, thr)
        for dia in list(mean_change[dur].keys()):
            if len(dia.split("_")) != 4 or "_".join(dia.split("_")[:-1]) != "std_ts_ano":
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        if not check_dict_levels(mean_bias_ave, dur, dia_r, pro, exp, dat):
                            continue
                        # array
                        da_m = copy__deepcopy(mean_change[dur][dia][pro][exp][dat])
                        del mean_change[dur][dia][pro][exp][dat]
                        for thr, da_r in mean_bias_ave[dur][dia_r][pro][exp][dat].items():
                            # select
                            list_models = sorted(list(set(da_m["dataset"].values.tolist()) &
                                                      set(da_r["dataset"].values.tolist())))
                            # save
                            da_o = da_m.sel(indexers={"dataset": list_models})
                            put_in_dict(mean_change_ave, da_o.compute(), dur, dia, pro, exp, dat, thr)
    # -- correlation
    print("correlations")
    p_value, r_value, model_nbr = {}, {}, {}
    for dur in list(mean_change_ave.keys()):
        dx = [mean_bias_ave] * len(list(mean_change_ave[dur].keys())) + \
             [mean_change_ave] * (len(list(mean_change_ave[dur].keys())) - 1)
        tx = ["bias"] * len(list(mean_change_ave[dur].keys())) + ["chan"] * (len(list(mean_change_ave[dur].keys())) - 1)
        lx = ["ave_pr_val_pepr"] * len(list(mean_change_ave[dur].keys())) + list(mean_change_ave[dur].keys())[:-1]
        dy = [mean_change_ave] * len(list(mean_change_ave[dur].keys())) + \
             [mean_change_ave] * (len(list(mean_change_ave[dur].keys())) - 1)
        ty = ["chan"] * len(list(mean_change_ave[dur].keys())) + ["chan"] * (len(list(mean_change_ave[dur].keys())) - 1)
        ly = list(mean_change_ave[dur].keys()) +  ["std_ts_ano_n34e"] * (len(list(mean_change_ave[dur].keys())) - 1)
        for nx, dx, dia_x, ny, dy, dia_y in zip(tx, dx, lx, ty, dy, ly):
            if not check_dict_levels(dx, dur, dia_x) or not check_dict_levels(dy, dur, dia_y):
                continue
            for pro1 in list(dy[dur][dia_y].keys()):
                for exp in list(dy[dur][dia_y][pro1].keys()):
                    for dat in list(dy[dur][dia_y][pro1][exp].keys()):
                        for thr in list(dy[dur][dia_y][pro1][exp][dat].keys()):
                            if not check_dict_levels(dx, dur, dia_x, pro1, exp, dat, thr):
                                continue
                            da_x = copy__deepcopy(dx[dur][dia_x][pro1][exp][dat][thr])
                            da_y = copy__deepcopy(dy[dur][dia_y][pro1][exp][dat][thr])
                            # model numbers
                            list_m = sorted(
                                list(set(da_x["dataset"].values.tolist()) & set(
                                    da_y["dataset"].values.tolist())))
                            for pro2 in ["cmip", "cmip5", "cmip6"]:
                                # select models
                                tmp = [k for k in list_m if pro2 == "cmip" or (pro2 == "cmip5" and k < 100) or
                                       (pro2 == "cmip6" and k > 100)]
                                da_xt = da_x.sel(indexers={"dataset": tmp})
                                da_yt = da_y.sel(indexers={"dataset": tmp})
                                # regression
                                ln = linear_regression(da_xt.values, da_yt.values)
                                lk = (dur, thr, pro2, exp, dia_x, dia_y)
                                put_in_dict(p_value, ln["p"], *lk)
                                put_in_dict(r_value, ln["r"], *lk)
                                put_in_dict(model_nbr, len(tmp), *lk)
    p_partial, r_partial = {}, {}
    for dur in list(r_value.keys()):
        for thr in list(r_value[dur].keys()):
            for pro in list(r_value[dur][thr].keys()):
                for exp in list(r_value[dur][thr][pro].keys()):
                    for dia_x in list(r_value[dur][thr][pro][exp].keys()):
                        if dia_x == "ave_pr_val_pepr":
                            continue
                        nn = model_nbr[dur][thr][pro][exp]["ave_pr_val_pepr"]["std_ts_ano_n34e"]
                        rpe = r_value[dur][thr][pro][exp]["ave_pr_val_pepr"]["std_ts_ano_n34e"]
                        rge = r_value[dur][thr][pro][exp][dia_x]["std_ts_ano_n34e"]
                        rpg = r_value[dur][thr][pro][exp]["ave_pr_val_pepr"][dia_x]
                        # partial correlation PR-ENSO accounting for Grad-ENSO
                        rpeg = (rpe - rge * rpg) / ((1 - rge ** 2) * (1 - rpg ** 2)) ** 0.5
                        tpeg = rpeg * ((nn - 3) / (1 - rpeg ** 2)) ** 0.5
                        ppeg = 2 * scipy__stats__t.sf(abs(tpeg), df=nn - 3)
                        # put_in_dict(p_partial, ppeg, dur, thr, pro, exp, "ave_pr_val_pepr", "std_ts_ano_n34e")
                        # put_in_dict(r_partial, rpeg, dur, thr, pro, exp, "ave_pr_val_pepr", "std_ts_ano_n34e")
                        # partial correlation Grad-ENSO accounting for PR-ENSO
                        rgep = (rge - rpe * rpg) / ((1 - rpe ** 2) * (1 - rpg ** 2)) ** 0.5
                        tgep = rgep * ((nn - 3) / (1 - rgep ** 2)) ** 0.5
                        pgep = 2 * scipy__stats__t.sf(abs(tgep), df=nn - 3)
                        put_in_dict(p_partial, pgep, dur, thr, pro, exp, dia_x, "std_ts_ano_n34e")
                        put_in_dict(r_partial, rgep, dur, thr, pro, exp, dia_x, "std_ts_ano_n34e")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 41: scatter
    print("scatter")
    group, panel = "scatter_01", "panel_41"
    if check_dict_levels(mean_bias_ave, dur2) and check_dict_levels(mean_change_ave, dur2):
        ldic_x = [mean_bias_ave[dur2]] + [mean_change_ave[dur2]] * (len(list(mean_change_ave[dur2].keys())) - 1)
        ltyp_x = ["bias"]  + ["chan"] * (len(list(mean_change_ave[dur2].keys())) - 1)
        ldia_x = ["ave_pr_val_pepr"] + list(mean_change_ave[dur2].keys())[:-1]
        ldic_y = [mean_change_ave[dur2]] * len(ldic_x)
        ltyp_y = ["chan"] * len(ldic_x)
        ldia_y = ["std_ts_ano_n34e"] * len(ldic_x)
        for ii, (dia_x, dia_y, dic_x, dic_y, typ_x, typ_y) in enumerate(
                zip(ldia_x, ldia_y, ldic_x, ldic_y, ltyp_x, ltyp_y)):
            for jj, thr in enumerate([str(k).zfill(3) for k in list(range(100, 151, 25)) + [500]]):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (grp, pan)
                print(dur2, dia_x.ljust(12), typ_x, dia_y.ljust(12), typ_y, thr, tuple_k)
                # panel size
                tk = tuple_k[:-2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, k2, *tk + (k1,))
                # panel axes definition
                dia = "--".join([dia_x, typ_x, dia_y, typ_y])
                if not check_dict_levels(fig_ticks, dia, panel):
                    print(dia, panel)
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
                            if len(ldia_x) == 3 and ii == 2:
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
                if ii == 0:
                    put_in_dict(plot_data, fig_titles[thr], *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, fig_titles[dia_x], *tuple_k + ("title_row",))
                # -- marker
                leg_col, leg_mar = [], []
                lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                for pro in list(dic_x[dia_x].keys()):
                    if pro != "cmip":
                        continue
                    if not check_dict_levels(dic_y, dia_y, pro):
                        continue
                    for exp_x in list(dic_x[dia_x][pro].keys()):
                        if (typ_x == "bia" and typ_y == "bia" and exp_x != "historical") or (
                                (typ_x == "cha" or typ_y == "cha") and exp_x == "historical"):
                            continue
                        for exp_y in list(dic_y[dia_y][pro].keys()):
                            if (typ_x == "bia" and typ_y == "bia" and exp_y != "historical") or (
                                    (typ_x == "cha" or typ_y == "cha") and exp_x != exp_y):
                                continue
                            for dat in list(dic_x[dia_x][pro][exp_x].keys()):
                                if not check_dict_levels(dic_x, dia_x, pro, exp_x, dat, thr) or \
                                        not check_dict_levels(dic_y, dia_y, pro, exp_y, dat, thr):
                                    continue
                                da_x = copy__deepcopy(dic_x[dia_x][pro][exp_x][dat][thr])
                                da_y = copy__deepcopy(dic_y[dia_y][pro][exp_y][dat][thr])
                                # model numbers corresponding to cmip5/6
                                list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                     set(da_y["dataset"].values.tolist())))
                                # list_m = [k for k in list_m if k > 100]
                                da_x = da_x.sel(indexers={"dataset": list_m})
                                da_y = da_y.sel(indexers={"dataset": list_m})
                                for dd, xx, yy in zip(list_m, da_x.values.tolist(), da_y.values.tolist()):
                                    pro2 = "cmip5" if dd < 100 else "cmip6"
                                    leg_col.append(pro2)
                                    leg_mar.append(exp_y)
                                    lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                                    lm += [fig_markers[exp_y] if exp_y in list(fig_markers.keys()) else "o"]
                                    ls += [80]
                                    lw += [0.5]
                                    lx += [xx]
                                    ly += [yy]
                                    lz += [numpy__random__randint(2, 20)]
                leg_col = sorted(list(set(leg_col)), key=lambda s: s.lower())
                leg_mar = sorted(list(set(leg_mar)), key=lambda s: s.lower())
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    # if ii == 0 and jj == 0:
                    #     for kk, mm in enumerate(leg_mar):
                    #         lc += [fig_colors[mm] if mm in list(fig_colors.keys()) else "cyan"]
                    #         lm += [fig_markers[mm] if mm in list(fig_markers.keys()) else "o"]
                    #         ls += [80]
                    #         lw += [0.5]
                    #         lx += [x2 + 10 * dx]
                    #         ly += [(y1 + y2) / 2 + (8 - kk * 8) * dy]
                    #         lz += [1]
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
                lc, lo, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    # x = 0 and y = 0
                    for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                        if k1[0] < 0 < k1[1]:
                            lc += ["k"]
                            lo += [True]
                            ls += [":"]
                            lw += [1]
                            lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                            ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                            lz += [23]
                    # observations
                    for k1, k2 in zip(["x", "y"], [dia_x, dia_y]):
                        dur3 = copy__deepcopy(dur2) if "_ts_" in k2 else "045_year_epoch"
                        dat3, exp3, pro3 = "multiobs001", "historical", "obse"
                        if not check_dict_levels(mean_bias_ave, dur3, k2, pro3, exp3, dat3):
                            continue
                        da1 = copy__deepcopy(mean_bias_ave[dur3][k2][pro3][exp3][dat3])
                        for da2 in da1:
                            lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                            lo += [True]
                            ls += [fig_linestyles[pro3] if pro3 in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [[da2.values] * 2 if k1 == "x" else dt["x_lim"]]
                            ly += [dt["y_lim"] if k1 == "x" else [da2.values] * 2]
                            lz += [22]
                plot_type = "cur"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lo), *tuple_k + (str(plot_type) + "_on",))
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
                    if check_dict_levels(p_value, dur2, thr) and check_dict_levels(r_value, dur2, thr):
                        list_pro = sorted(list(r_value[dur2][thr].keys()), key=lambda s: s.lower())
                        for kk, pro in enumerate(list_pro):
                            if not check_dict_levels(p_value, dur2, thr, pro):
                                continue
                            pp, rr = p_value[dur2][thr][pro], r_value[dur2][thr][pro]
                            tp , tr = [], []
                            for exp in list(pp.keys()):
                                if check_dict_levels(pp, exp, dia_x, dia_y):
                                    tp += [pp[exp][dia_x][dia_y]]
                                    tr += [rr[exp][dia_x][dia_y]]
                            if len(tp) == 0:
                                continue
                            tr = float(numpy__array(tr).mean())
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "k"] * 3
                            lf += [15] * 3
                            lh += ["right", "right", "left"]
                            lr += [0] * 3
                            lt += ["r=", "{0:.2f}".format(round(tr, 2)),
                                   "".join(["*" if vv < 1e-2 else "" for vv in tp])]
                            lv += ["center"] * 3
                            lw += ["regular"] * 3
                            lx += [x2 + 15. * dx, x2 + 32.5 * dx, x2 + 33 * dx]
                            ly += [y2 - (kk * 8 + 10) * dy] * 3
                        lc += ["grey"]
                        lf += [15]
                        lh += ["left"]
                        lr += [0]
                        lt += ["correlation"]
                        lv += ["bottom"]
                        lw += ["regular"]
                        lx += [x2 + 4 * dx]
                        ly += [y2 - (10 - 6) * dy]
                    if check_dict_levels(p_partial, dur2, thr) and check_dict_levels(r_partial, dur2, thr):
                        list_pro = sorted(list(r_partial[dur2][thr].keys()), key=lambda s: s.lower())
                        for kk, pro in enumerate(list_pro):
                            if not check_dict_levels(p_partial, dur2, thr, pro):
                                continue
                            pp, rr = p_partial[dur2][thr][pro], r_partial[dur2][thr][pro]
                            tp, tr = [], []
                            for exp in list(pp.keys()):
                                if check_dict_levels(pp, exp, dia_x, dia_y):
                                    tp += [pp[exp][dia_x][dia_y]]
                                    tr += [rr[exp][dia_x][dia_y]]
                            if len(tp) == 0:
                                continue
                            tr = float(numpy__array(tr).mean())
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "k"] * 3
                            lf += [15] * 3
                            lh += ["right", "right", "left"]
                            lr += [0] * 3
                            lt += ["r$_p$=", "{0:.2f}".format(round(tr, 2)),
                                   "".join(["*" if vv < 1e-2 else "" for vv in tp])]
                            lv += ["center"] * 3
                            lw += ["regular"] * 3
                            lx += [x2 + 15. * dx, x2 + 32.5 * dx, x2 + 33 * dx]
                            ly += [y2 - ((kk + 6) * 8 + 10) * dy] * 3
                        lc += ["grey"]
                        lf += [15]
                        lh += ["left"]
                        lr += [0]
                        lt += ["partial corr."]
                        lv += ["bottom"]
                        lw += ["regular"]
                        lx += [x2 + 4 * dx]
                        ly += [y2 - (6 * 8 + 10 - 6) * dy]
                    # legend
                    # if ii == 0 and jj == 0:
                    #     for kk, mm in enumerate(["cmip"] + leg_col + leg_mar):
                    #         lc += [fig_colors[mm] if mm in list(fig_colors.keys()) else "cyan"]
                    #         lf += [15]
                    #         lh += ["right" if mm not in leg_mar else "left"]
                    #         lr += [0]
                    #         lt += [fig_titles[mm]]
                    #         lv += ["center"]
                    #         lw += ["regular"]
                    #         lx += [x2 - 2 * dx if mm not in leg_mar else x2 + 14 * dx]
                    #         ly += [y2 - (kk * 7 + 4.5) * dy if mm not in leg_mar else
                    #                (y1 + y2) / 2 + (8 - (kk - 3) * 8) * dy]
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
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a4) + "rem"])
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
e81_change_vs_change_bias_control_plot()
