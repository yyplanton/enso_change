# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: Correlations between bias or change and ENSO
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
from numpy import sign as numpy__sign

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import correlation_pearson
from enso_change.wrapper.xarray_base import concatenate, xarray, array_wrapper, numpy__nan, cf_dim_to_dim, mean_weighted
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_tx_val", "ave_ts_val",
                          "ave_ts_val_glob", "ave_ts_val_tron"],
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
        "frac": {"x": 0.25, "y": 0.25},
        "panel_01": {"x_delt": 5, "x_size": 20, "y_delt": 3, "y_size": 6},
    },
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": {},
    # ticks
    "fig_ticks": {
        "ave_pr_val--bias": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": list(range(4, 9, 4)),
                "map_l2_tic": list(range(4, 9, 4)),
                "region": None,
            },
        },
        "ave_pr_val--chan": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": [round(k / 10, 1) for k in range(0, 9, 4)],
                "map_l2_tic": [round(k / 10, 1) for k in range(0, 9, 4)],
                "region": None,
            },
        },
        "ave_rts_val--bias": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": list(range(-2, 3, 2)),
                "map_l2_tic": list(range(-2, 3, 2)),
                "region": None,
            },
        },
        "ave_rts_val--chan": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": [round(k / 100, 2) for k in range(-15, 16, 15)],
                "map_l2_tic": [round(k / 100, 2) for k in range(-15, 16, 15)],
                "region": None,
            },
        },
        "ave_ts_val--bias": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": list(range(24, 29, 2)),
                "map_l2_tic": list(range(24, 29, 2)),
                "region": None,
            },
        },
        "ave_ts_val--chan": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": [round(k / 10, 1) for k in range(7, 10, 1)],
                "map_l2_tic": [round(k / 10, 1) for k in range(7, 10, 1)],
                "region": None,
            },
        },
        "ave_tx_val--bias": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": list(range(-6, 1, 2)),
                "map_l2_tic": list(range(-6, 1, 2)),
                "region": None,
            },
        },
        "ave_tx_val--chan": {
            "panel_01": {
                "x_lim": [100, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": [round(k / 10, 1) for k in range(-2, 3, 2)],
                "map_l2_tic": [round(k / 10, 1) for k in range(-2, 3, 2)],
                "region": None,
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_gengts": "Gen24",
        "ave_ts_val_heedts": "HF21",
        "ave_ts_val_seagts": "Sea22",
        "ave_ts_val_wangts": "Wan24",
        "ave_ts_val_watats": "Wat21",
        "cmip": "CMIP 5&6",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        # figure title
        "figure_title": "Gradient change and ENSO amplitude change",
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
def e77_maps_bias_change_vs_change_plot(
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
        data_region={"X": (80, 320), "Y": (-50, 50)}, smile_minimum_size=smile_minimum_size,
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
    # -- compute rSST
    dia2, dia3 = "ave_ts_val_tron", [k for k in data_diagnostics2 if k in ["ave_sst_val", "ave_ts_val"]][0]
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
        if check_dict_levels(values_netcdf, dia_o):
            if not check_dict_levels(metadata_netcdf, dia_o):
                metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia3])
                val = metadata_netcdf[dia_o]["name_short"].replace("TS", "SST").replace("SST", "rSST")
                metadata_netcdf[dia_o]["name_short"] = val
    # -- select bias year
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_bias = {}
    for dia1 in list(values_netcdf.keys()):
        if len(dia1.split("_")) != 3:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp in ["piControl", "historical"]:
                        continue
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
                        put_in_dict(mean_bias, copy__deepcopy(da.compute()), dur1, dia1, pro, exp, dat)
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
                        put_in_dict(mean_change, copy__deepcopy(da1.compute()), dur1, dia1, pro, exp1, dat)
    print("mean_change", list(mean_change.keys()))
    default_region = {
        "east_pla2": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EaEqP$_{Pla27}$"},
        "west_pla2": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "EaEqP$_{Pla27}$"},
        "pla2ts": {"regions": ["west_pla2", "east_pla2"], "short_name": "Pla27$_{10-10}$"},
    }
    list_regions = ["pla2ts"]
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if dia not in ["ave_pr_val", "ave_ts_val"]:
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat])
                        mask = create_land_sea_mask(da)
                        da = da.where(mask == 0)
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
                                if da_o is None:
                                    da_o = copy__deepcopy(da_t)
                                elif reg == "geng_ts" and r2 in ["west_gen", "nort_gen", "sout_gen"]:
                                    da_o += da_t
                                    if r2 == "sout_gen":
                                        da_o /= 3
                                else:
                                    da_o -= da_t
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
    # -- correlations
    correlations = {}
    dia_y = "ave_ts_val_pla2ts"  # "std_ts_ano_n34e"  #
    for tt, dd in zip(["bias", "chan"], [mean_bias, mean_change]):
        for dur in list(dd.keys()):
            for dia_x in list(dd[dur].keys()):
                if len(dia_x.split("_")) != 3:
                    continue
                for pro1 in list(dd[dur][dia_x].keys()):
                    for exp in list(dd[dur][dia_x][pro1].keys()):
                        for dat in list(dd[dur][dia_x][pro1][exp].keys()):
                            if not check_dict_levels(mean_change, dur, dia_y, pro1, exp1, dat) and \
                                    not check_dict_levels(mean_change_ave, dur, dia_y, pro1, exp1, dat):
                                continue
                            # array
                            da_x = copy__deepcopy(dd[dur][dia_x][pro1][exp][dat])
                            if check_dict_levels(mean_change, dur, dia_y, pro1, exp1, dat):
                                da_y = copy__deepcopy(mean_change[dur][dia_y][pro1][exp][dat])
                            else:
                                da_y = copy__deepcopy(mean_change_ave[dur][dia_y][pro1][exp][dat])
                            # correlation CMIP5 & CMIP6
                            rr = correlation_pearson(da_x, da_y, dim="dataset")
                            # save
                            lk = (tt, dur, dia_y, dia_x)
                            put_in_dict(correlations, copy__deepcopy(rr.compute()),*lk + (pro1, exp,))
                            # separate CMIP5 and CMIP6 datasets
                            list_mod = da_x["dataset"].values.tolist()
                            for pro2 in ["cmip5", "cmip6"]:
                                list_t = sorted([k for k in list_mod
                                                 if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)])
                                # select models
                                da_xt = da_x.sel(indexers={"dataset": list_t})
                                da_yt = da_y.sel(indexers={"dataset": list_t})
                                # correlation CMIP5 or CMIP6
                                rr = correlation_pearson(da_xt, da_yt, dim="dataset")
                                # save
                                put_in_dict(correlations, copy__deepcopy(rr.compute()), *lk + (pro2, exp,))
    print("correlations", len(list(correlations.keys())))
    # -- average correlations across experiments
    to_plot = {}
    for typ in list(correlations.keys()):
        for dur in list(correlations[typ].keys()):
            for dia_y in list(correlations[typ][dur].keys()):
                for dia_x in list(correlations[typ][dur][dia_y].keys()):
                    for pro in list(correlations[typ][dur][dia_y][dia_x].keys()):
                        # array
                        da = tuple(da for da in correlations[typ][dur][dia_y][dia_x][pro].values())
                        da = concatenate(da, dim="experiment")
                        # median correlation across maps
                        mean = da.mean(dim="experiment", skipna=True)
                        # sign of the median
                        sign = xarray.where(mean > 0, 1, xarray.where(mean < 0, -1, numpy__nan))
                        # does each map have the same sign as the median?
                        same_sign = xarray.where(numpy__sign(da) == sign, 1.0, 0.0)
                        # fraction of maps agreeing with the median sign
                        agre = same_sign.where(da.notnull()).mean(dim="experiment")
                        # save
                        put_in_dict(to_plot, agre.compute(), typ, dur, dia_y, dia_x, pro, "p")
                        put_in_dict(to_plot, mean.compute(), typ, dur, dia_y, dia_x, pro, "r")
    # -- average bias / changes across experiments
    for tt, dd in zip(["bias", "chan"], [mean_bias, mean_change]):
        for dur in list(dd.keys()):
            for dia_x in list(dd[dur].keys()):
                if len(dia_x.split("_")) != 3:
                    continue
                for pro1 in list(dd[dur][dia_x].keys()):
                    list_dat = sorted(list(set([k2 for k1 in list(dd[dur][dia_x][pro1].keys())
                                                for k2 in list(dd[dur][dia_x][pro1][k1].keys())])),
                                      key=lambda s: s.lower())
                    for dat in list_dat:
                        list_exp = sorted(list(set([k1 for k1 in list(dd[dur][dia_x][pro1].keys())
                                                    if dat in list(dd[dur][dia_x][pro1][k1].keys())])),
                                          key=lambda s: s.lower())
                        for pro2 in ["cmip", "cmip5", "cmip6"]:
                            # get data for each experiment
                            da_o = ()
                            for exp in list_exp:
                                if not check_dict_levels(mean_change, dur, dia_y, pro1, exp, dat) and \
                                        not check_dict_levels(mean_change_ave, dur, dia_y, pro1, exp, dat):
                                    continue
                                # array
                                da_x = copy__deepcopy(dd[dur][dia_x][pro1][exp][dat])
                                # separate CMIP5 and CMIP6 datasets
                                list_mod = da_x["dataset"].values.tolist()
                                list_mod = sorted([k for k in list_mod if pro2 == "cmip" or
                                                   (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)])
                                # select models
                                da_o += (da_x.sel(indexers={"dataset": list_mod}).mean(dim="dataset"),)
                            da_o = concatenate(da_o, dim="experiment")
                            # mean
                            da_o = da_o.mean(dim="experiment", skipna=True)
                            if "_tx_" in dia_x:
                                da_o *= 100
                            print(tt, dur, dia_y, dia_x.ljust(11), pro2.ljust(5),
                                  "{0:.2f}".format(round(float(da_o.min().values), 2)).rjust(10),
                                  "{0:.2f}".format(round(float(da_o.max().values), 2)).rjust(10))
                            # save
                            put_in_dict(to_plot, da_o.compute(), tt, dur, dia_y, dia_x, pro2, "m")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    group, panel = "goup_01", "panel_01"
    print("map", group, panel)
    for typ in list(to_plot.keys()):
        for dur in list(to_plot[typ].keys()):
            for dia_y in list(to_plot[typ][dur].keys()):
                list_dia_x = list(to_plot[typ][dur][dia_y].keys())
                for ii, dia_x in enumerate(list_dia_x):
                    list_pro = list(to_plot[typ][dur][dia_y][dia_x].keys())
                    for jj, pro in enumerate(list_pro):
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = ("_".join([typ, dur, dia_y]), grp, pan)
                        print(typ, dur, dia_y, dia_x.ljust(11), pro.ljust(5), tuple_k)
                        # panel size
                        tk = tuple_k[:-2]
                        if not check_dict_levels(figure_psize, *tk + (pan,)):
                            put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                            for k1, k2 in fig_panel_size.items():
                                if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                                    put_in_dict(figure_psize, k2, *tk + (k1,))
                        # panel axes definition
                        dia_t = "--".join([dia_x, typ])
                        if check_dict_levels(fig_ticks, dia_t, panel):
                            for k1, k2 in fig_ticks[dia_t][panel].items():
                                if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                                    continue
                                tmp = copy__deepcopy(k2)
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
                            put_in_dict(plot_data, fig_titles[pro], *tuple_k + ("title_col",))
                        if jj == 0:
                            val = copy__deepcopy(metadata_netcdf[dia_x]["name_short"])
                            put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                        # -- map
                        da_m = to_plot[typ][dur][dia_y][dia_x][pro]["m"]
                        da_p = to_plot[typ][dur][dia_y][dia_x][pro]["p"]
                        da_r = to_plot[typ][dur][dia_y][dia_x][pro]["r"]
                        cond = (da_r > -0.2) & (da_r < 0.2)
                        da_r = da_r.where(~cond)
                        plot_type = "map"
                        # shading
                        put_in_dict(plot_data, copy__deepcopy(da_r), *tuple_k + (str(plot_type) + "_c",))
                        val = "bottom" if dia_x == list_dia_x[-1] else None
                        put_in_dict(plot_data, val, *tuple_k + ("legend_position",))
                        # hatching
                        put_in_dict(plot_data, copy__deepcopy(da_p), *tuple_k + (str(plot_type) + "_h",))
                        density = 3
                        kwargs_h = {"hatches": [density * "/"], "levels": [0.0, 0.9]}
                        put_in_dict(plot_data, kwargs_h, *tuple_k + ("kwargs_h",))
                        # contour
                        ct = "l2"
                        put_in_dict(plot_data, copy__deepcopy(da_m), *tuple_k + (str(plot_type) + "_" + str(ct),))
                        dd = {
                            "colors": ["k"] * len(dt["map_" + str(ct) + "_tic"]),
                            "linestyles": ["-" if k >= 0 else "--" for k in dt["map_" + str(ct) + "_tic"]],
                            "linewidths": [1 if k != 0 else 2 for k in dt["map_" + str(ct) + "_tic"]],
                        }
                        put_in_dict(plot_data, copy__deepcopy(dd), *tuple_k + ("kwargs_" + str(ct),))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1, a4 = str(smile_minimum_size).zfill(2), str(len(smile_rejected)).zfill(2)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a4) + "rem"])
    # figure loop
    for k in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[k].keys())
        numb = False if len(list_names) < 2 else True
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", k)))
        # plot
        fig_basic(plot_data[k], list_names, 1, figure_axes[k], fig_format, figure_name_and_path,
                  figure_psize[k], fig_title_bool=numb, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e77_maps_bias_change_vs_change_plot()
