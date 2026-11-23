# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Maps of biases (values at the end of the 20th century) and changes (departure from the end of the 20th century)
#  - Distributions of selected bias and change
#  - Relationship between selected bias and mean stat change and ENSO SSTA change
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
from numpy.random import default_rng as numpy__random__default_rng
from numpy.random import randint as numpy__random__randint
# scipy
from scipy.stats import percentileofscore as scipy__stats__percentileofscore

# scipy
from scipy.stats import scoreatpercentile as scipy__stats__scoreatpercentile

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, numpy__ndarray, scipy__stats__t, \
    linear_regression, percentile
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, mean_weighted, polyval
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val", "ave_ts_val_glob", "ave_ts_val_tron"],
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
        "panel_01": {"x_delt": 250, "x_size": 1440, "y_delt": 220, "y_size": 480},
        "panel_11": {"x_delt": 250, "x_size": 1440, "y_delt": 250, "y_size": 480},
        "panel_21": {"x_delt": 250, "x_size": 630, "y_delt": 220, "y_size": 600},
        "panel_31": {"x_delt": 250, "x_size": 630, "y_delt": 250, "y_size": 600},
        "panel_41": {"x_delt": 730, "x_size": 700, "y_delt": 300, "y_size": 700},
        "panel_51": {"x_delt": 180, "x_size": 4, "y_delt": 100, "y_size": 700},
        "panel_52": {"x_delt": 180, "x_size": 4, "y_delt": 100, "y_size": 700},
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
        "ave_pr_val": {
            "panel_01": {
                "x_lim": [115, 295],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "BrBG",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": [3, 6],  # list(range(3, 10, 2)),
                "map_l2_tic": [4, 8],  # list(range(3, 10, 2)),
                "region": ["p120pr"],
            },
        },
        "ave_pr_val_pepr": {
            "panel_21": {
                "x_lim": None,
                "x_nam": "",
                "x_nbr_minor": 0,
                "x_tic": None,
                "y_lim": [0.5, 6.5],
                "y_nam": "PR bias (UNITS)",
                "y_nbr_minor": 2,
                "y_tic": [round(k / 10, 1) for k in range(5, 66, 15)],
                # "y_lim": [0.0, 6.0],
                # "y_nam": "VALUE (UNITS)",
                # "y_tic": [round(k / 10, 1) for k in range(0, 61, 15)],
            },
        },
        "ave_pr_val_pepr--bia--std_ts_ano_n34e--cha": {
            "panel_41": {
                # "x_lim": [0.2, 5.4],
                # "x_nam": "VALUE (UNITS)",
                # "x_tic": [round(k / 10, 1) for k in range(2, 55, 13)],
                "x_lim": [0., 6.0],
                "x_nam": "PR bias (UNITS)",
                "x_nbr_minor": 2,
                "x_tic": [round(k / 10, 1) for k in range(0, 61, 15)],
                # "x_lim": [0.0, 6.0],
                # "x_nam": "VALUE (UNITS)",
                # "x_tic": [round(k / 10, 1) for k in range(0, 61, 15)],
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "ave_rts_val": {
            "panel_01": {
                "x_lim": [115, 295],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-2, 3, 1)],
                # "map_l1_tic": [round(k / 100, 2) for k in range(-10, 11, 5)],
                "map_l1_tic": [round(k / 10, 1) for k in range(-1, 2, 1)],
                "map_l2_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "region": ["pl2w", "pl2e"], #  ["plaw", "plae"],
            },
        },
        # "ave_ts_val_pgts": {
        #     "panel_21": {
        #         "x_lim": None,
        #         "x_nam": "",
        #         "x_nbr_minor": 0,
        #         "x_tic": None,
        #         "y_lim": [-0.4, 0.4],
        #         "y_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
        #         "y_tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
        #     },
        # },
        # "ave_ts_val_pgts--cha--std_ts_ano_n34e--cha": {
        #     "panel_41": {
        #         "x_lim": [-0.4, 0.4],
        #         "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
        #         "x_nbr_minor": 1,
        #         "x_tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
        #         "y_axe": "right",
        #         "y_lim": [-45, 30],
        #         "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
        #         "y_nbr_minor": 2,
        #         "y_tic": list(range(-45, 31, 15)),
        #     },
        # },
        "ave_ts_val_pgts": {
            "panel_21": {
                "x_lim": None,
                "x_nam": "",
                "x_nbr_minor": 0,
                "x_tic": None,
                "y_lim": [-0.6, 0.3],
                "y_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": [round(k / 10, 1) for k in range(-6, 4, 3)],
            },
        },
        "ave_ts_val_pgts--cha--std_ts_ano_n34e--cha": {
            "panel_41": {
                "x_lim": [-0.8, 0.4],
                "x_nam": "$\\Delta_{x}$SST change (UNITS.$^\\circ$C$^{-1}$)",
                "x_nbr_minor": 1,
                "x_tic": [round(k / 10, 1) for k in range(-8, 5, 4)],
                "y_axe": "right",
                "y_lim": [-45, 30],
                "y_nam": "$\\Delta_{t}$E(t=2050-2099) (%.$^\\circ$C$^{-1}$)",
                "y_nbr_minor": 2,
                "y_tic": list(range(-45, 31, 15)),
            },
        },
        "std_ts_ano_n34e": {
            "panel_51": {
                "lim": [-20, 40],
                "nam": "$\\Delta_{t}$E relative to 1956-2005 (%)",
                "tic": list(range(-40, 41, 20)),
            },
            "panel_52": {
                "lim": [-20, 40],
                "nam": "",
                "tic": list(range(-40, 41, 20)),
            },
        },
        "time": {
            "panel_51": {
                "lim": [1860, 2005],
                "nam": "Year",
                "nbr_minor": 4,
                "tic": list(range(1850, 2006, 50)),
            },
            "panel_52": {
                "lim": [1990, 2090],
                "nam": "Year",
                "nbr_minor": 4,
                "tic": list(range(2000, 2101, 50)),
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_pr_val": "VALUE (UNITS)",
        "ave_pr_val_pepr": "VALUE (UNITS)",
        "ave_rts_val": "VALUE (UNITS.$^\\circ$C$^{-1}$)",
        "ave_ts_val": "VALUE (UNITS)",
        "ave_ts_val_pets": "VALUE (UNITS)",
        "ave_ts_val_pgxd": "VALUE (UNITS)",
        "std_ts_ano_n34e": "$\\Delta_{t}$E relative to 1956-2005 (%)",
        "cmip": "CMIP 5&6",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        "diffe": "CMIP6 minus CMIP5",
        "scenario": "RCP-SSP",
        # "driest": "driest THR%",
        # "wettest": "wettest THR%",
        "driest": "driest",
        "wettest": "wettest",
        # "strongest": "strongest $\\Delta_{y}$ THR%",
        # "weakest": "weakest $\\Delta_{y}$ THR%",
        "largest": "largest",
        "smallest": "smallest",
        "mean": "mean",
        "range": "Forced-response\nspread",
        "bia": "mean-state bias\n(1956-2005)",
        "cha": "warming pattern\n(2050-2099 minus 1956-2005)",
        # figure title
        # "figure_title": "Mean state alter ENSO amplitude change",
        "figure_title": "Pacific mean-state bias and warming pattern linked to ENSO change",
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
def f03c_bias_vs_change_plot(
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
        # "plae": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EEP"},
        # "plaw": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "WEP"},
        "plae": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EEP"},
        "plaw": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "WEP"},
        "pl2e": {"latitude": [-6., 6.], "longitude": [230., 270.], "short_name": "EEP"},
        "pl2w": {"latitude": [-5., 5.], "longitude": [135., 185.], "short_name": "WEP"},
        "pgxe": {"latitude": [-5., 5.], "longitude": [220., 280.], "short_name": "EEP"},
        "pgxw": {"latitude": [-5., 5.], "longitude": [140., 190.], "short_name": "WEP"},
        "pepr": {"latitude": [-8., 3.], "longitude": [190., 280.], "short_name": ""},
        # "pgts": {"regions": ["pgxw", "pgxe"], "short_name": "$\\Delta_{x}$EP"},
        # "pgts": {"regions": ["plaw", "plae"], "short_name": "W - E"},
        "pgts": {"regions": ["pl2w", "pl2e"], "short_name": "W - E"},
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
    list_regions = ["pgts"]  # list(default_region.keys())
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if dia not in ["ave_pr_val", "ave_ts_val"]:
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
    # -- compute mean maps
    mean_map = {}
    for nn, dd in zip(["bia", "cha"], [mean_bias, mean_change]):
        for dur in list(dd.keys()):
            for dia in list(dd[dur].keys()):
                if len(dia.split("_")) != 3:
                    continue
                for pro1 in list(dd[dur][dia].keys()):
                    for exp in list(dd[dur][dia][pro1].keys()):
                        if (nn == "bia" and exp != "historical") or (nn == "cha" and exp == "historical"):
                            continue
                        for dat in list(dd[dur][dia][pro1][exp].keys()):
                            # array
                            da = copy__deepcopy(dd[dur][dia][pro1][exp][dat])
                            # list models by number
                            list_mod = da["dataset"].values.tolist()
                            # mean
                            da_m = da.mean(dim="dataset")
                            da_m = da_m.where(create_land_sea_mask(da_m) == 0)
                            # save
                            pro_o = "obse" if pro1 == "obse" else "CMIP6 minus CMIP5"
                            tuple_k = (dur, nn, exp, dat, dia, pro_o)
                            put_in_dict(mean_map, copy__deepcopy(da_m.compute()), *tuple_k + ("l",))
                            if pro1 != "cmip":
                                continue
                            # separate cmip5/6
                            dt = {}
                            l2 = ["cmip5", "cmip6"] if pro1 == "cmip" else [pro1]
                            for pro2 in l2:
                                # model numbers corresponding to cmip5/6
                                l3 = [k for k in list_mod
                                      if (pro2 != "cmip6" and k < 100) or (pro2 == "cmip6" and k > 100)]
                                # select historical bias of models with projections
                                dt[pro2] = da.sel(indexers={"dataset": l3}).mean(dim="dataset")
                                # save
                                put_in_dict(mean_map, dt[pro2], dur, nn, exp, dat, dia, pro2, "s")
                            # compare cmip5/6
                            da5, da6 = dt["cmip5"], dt["cmip6"]
                            # difference of multimodel mean
                            da_diff = da6 - da5
                            # save
                            put_in_dict(mean_map, copy__deepcopy(da_diff.compute()), *tuple_k + ("s",))
    # -- region distributions
    mean_box = {}
    for nn, dd in zip(["bia", "cha"], [mean_bias_ave, mean_change_ave]):
        for dur in list(dd.keys()):
            for dia in list(dd[dur].keys()):
                for pro1 in list(dd[dur][dia].keys()):
                    for exp in list(dd[dur][dia][pro1].keys()):
                        if (nn == "bia" and exp != "historical") or (nn == "cha" and exp == "historical"):
                            continue
                        for dat in list(dd[dur][dia][pro1][exp].keys()):
                            # array
                            da = copy__deepcopy(dd[dur][dia][pro1][exp][dat])
                            # list models by number
                            list_mod = da["dataset"].values.tolist()
                            # separate cmip5/6
                            l2 = ["cmip5", "cmip6"] if pro1 == "cmip" else [pro1]
                            for pro2 in l2:
                                # model numbers corresponding to cmip5/6
                                l3 = [k for k in list_mod
                                      if (pro2 != "cmip6" and k < 100) or (pro2 == "cmip6" and k > 100)]
                                # select historical bias of models with projections
                                da_o = copy__deepcopy(da).sel(indexers={"dataset": l3})
                                # save
                                lk = (dur, nn, exp, dat, dia, pro2)
                                put_in_dict(mean_box, copy__deepcopy(da_o.compute()), *lk)
    print("box ave min-max")
    dt = {}
    for dur in list(mean_box.keys()):
        for typ in list(mean_box[dur].keys()):
            for exp in list(mean_box[dur][typ].keys()):
                for dat in list(mean_box[dur][typ][exp].keys()):
                    for dia in list(mean_box[dur][typ][exp][dat].keys()):
                        for da in mean_box[dur][typ][exp][dat][dia].values():
                            put_in_dict(dt, da.values.flatten().tolist(), typ, dia)
    for typ in sorted(list(dt.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(dt[typ].keys()), key=lambda s: s.lower()):
            if len(dt[typ][dia]) > 0:
                print(dia.ljust(12), "{0:.2f}".format(round(min(dt[typ][dia]), 2)).rjust(8),
                      "{0:.2f}".format(round(max(dt[typ][dia]), 2)).rjust(8))
    # -- select top / bottom models
    change_m, change_r = {}, {}
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    for dia_x in list(mean_bias_ave[dur2].keys()):
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
                            ens_size = {int(k.split(" ")[0]): int(k.split("=")[1])
                                        for k in da3.attrs["dataset"].split("; ")}
                            # model nbr with dia_x value
                            dict_t = {int(k): float(v) for k, v in zip(da3["dataset"].values, da3.values)}
                            dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
                            list_datasets = list(dict_t.keys())
                            list_datasets = [k for k in list_datasets if k > 100]
                            pro_o = "cmip6"
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
                                # final spread (initial spread minus residual internal variability) -> must be set to 0
                                # if negative
                                sf = (si - ri) / si
                                sf = sf.where(sf > 0, other=0)
                                # scaled distribution
                                da_r = da_m + (da_t - da_m) * sf
                                # confidence interval
                                da_r = percentile(da_r, dim="dataset", per=percentiles, skipna=True)
                                # centered time
                                da_m = da_m.assign_coords({"year": da_m["year"] - int(dur1.split("_")[0]) / 2})
                                da_r = da_r.assign_coords({"year": da_r["year"] - int(dur1.split("_")[0]) / 2})
                                lk = (dur1, dia_y, dia_x, exp1, pro_o, grp)
                                put_in_dict(change_m, copy__deepcopy(da_m.compute()), *lk)
                                put_in_dict(change_r, copy__deepcopy(da_r.compute()), *lk)
    print("change ts min-max")
    for dur in list(change_r.keys()):
        dt = {}
        for dia_y in list(change_r[dur].keys()):
            for dia_x in list(change_r[dur][dia_y].keys()):
                for exp in list(change_r[dur][dia_y][dia_x].keys()):
                    for pro in list(change_r[dur][dia_y][dia_x][exp].keys()):
                        for da in change_r[dur][dia_y][dia_x][exp][pro].values():
                            put_in_dict(change_r, da.values.flatten().tolist(), dia_y)
        for dia in sorted(list(dt.keys()), key=lambda s: s.lower()):
            print(dia, "{0:.2f}".format(round(min(dt[dia]), 2)).rjust(8),
                  "{0:.2f}".format(round(max(dt[dia]), 2)).rjust(8))
    # -- correlation
    print("correlations")
    p_value, r_value, model_nbr = {}, {}, {}
    for dur in list(mean_change_ave.keys()):
        for nx, dx, ny, dy in zip(["bias", "bias", "chan"], [mean_bias_ave, mean_bias_ave, mean_change_ave],
                                  ["chan", "chan", "chan"], [mean_change_ave, mean_change, mean_change]):
            if not check_dict_levels(dx, dur) or not check_dict_levels(dy, dur):
                continue
            for dia_x in list(dx[dur].keys()):
                if "ave_pr_val" not in dia_x and "ave_ts_val_pgts" not in dia_x:
                    continue
                for dia_y in list(dy[dur].keys()):
                    if len(dia_y.split("_")) != 4 or ("ave_ts_val_pgts" not in dia_y and "std_ts_ano" not in dia_y):
                        continue
                    for pro1 in list(dy[dur][dia_y].keys()):
                        for exp in list(dy[dur][dia_y][pro1].keys()):
                            for dat in list(dy[dur][dia_y][pro1][exp].keys()):
                                if not check_dict_levels(dx, dur, dia_x, pro1, exp, dat):
                                    continue
                                da_x = copy__deepcopy(dx[dur][dia_x][pro1][exp][dat])
                                da_y = copy__deepcopy(dy[dur][dia_y][pro1][exp][dat])
                                # model numbers
                                list_m = sorted(
                                    list(set(da_x["dataset"].values.tolist()) & set(da_y["dataset"].values.tolist())))
                                for pro2 in ["cmip", "cmip5", "cmip6"]:
                                    # select models
                                    tmp = [k for k in list_m
                                           if pro2 == "cmip" or (pro2 == "cmip5" and k < 100) or (
                                                       pro2 == "cmip6" and k > 100)]
                                    da_xt = da_x.sel(indexers={"dataset": tmp})
                                    da_yt = da_y.sel(indexers={"dataset": tmp})
                                    # regression
                                    ln = linear_regression(da_xt.values, da_yt.values)
                                    put_in_dict(p_value, ln["p"], dur, pro2, exp, "--".join([dia_x, dia_y]))
                                    put_in_dict(r_value, ln["r"], dur, pro2, exp, "--".join([dia_x, dia_y]))
                                    put_in_dict(model_nbr, len(tmp), dur, pro2, exp, "--".join([dia_x, dia_y]))
    p_partial, r_partial = {}, {}
    for dur in list(r_value.keys()):
        for pro in list(r_value[dur].keys()):
            for exp in list(r_value[dur][pro].keys()):
                nn = model_nbr[dur][pro][exp]["ave_pr_val_pepr--std_ts_ano_n34e"]
                rpe = r_value[dur][pro][exp]["ave_pr_val_pepr--std_ts_ano_n34e"]
                rge = r_value[dur][pro][exp]["ave_ts_val_pgts--std_ts_ano_n34e"]
                rpg = r_value[dur][pro][exp]["ave_pr_val_pepr--ave_ts_val_pgts"]
                # partial correlation PR-ENSO accounting for Grad-ENSO
                rpeg = (rpe - rge * rpg) / ((1 - rge ** 2) * (1 - rpg ** 2)) ** 0.5
                tpeg = rpeg * ((nn - 3) / (1 - rpeg**2))**0.5
                ppeg = 2 * scipy__stats__t.sf(abs(tpeg), df=nn - 3)
                # print(dur, pro.ljust(5), exp, "r=", "{0:.2f}".format(round(rpeg, 2)).rjust(6), "t=",
                #       "{0:.3f}".format(round(ppeg, 3)).rjust(6))
                put_in_dict(p_partial, ppeg, dur, pro, exp, "ave_pr_val_pepr--std_ts_ano_n34e")
                put_in_dict(r_partial, rpeg, dur, pro, exp, "ave_pr_val_pepr--std_ts_ano_n34e")
                # partial correlation Grad-ENSO accounting for PR-ENSO
                rgep = (rge - rpe * rpg) / ((1 - rpe**2) * (1 - rpg**2))**0.5
                tgep = rgep * ((nn - 3) / (1 - rgep ** 2)) ** 0.5
                pgep = 2 * scipy__stats__t.sf(abs(tgep), df=nn - 3)
                # print(dur, pro.ljust(5), exp, "r=", "{0:.2f}".format(round(rgep, 2)).rjust(6), "t=",
                #       "{0:.3f}".format(round(pgep, 3)).rjust(6))
                put_in_dict(p_partial, pgep, dur, pro, exp, "ave_ts_val_pgts--std_ts_ano_n34e")
                put_in_dict(r_partial, rgep, dur, pro, exp, "ave_ts_val_pgts--std_ts_ano_n34e")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    dat2, dur2 = "multimodel001", "050_year_epoch"
    dat3, dur3, pro3 = "multiobs001", "045_year_epoch", "obse"
    list_typ = ["bia", "cha"]
    for typ in list_typ:
        exp2 = "historical" if typ == "bia" else "scenario85"
        list_dia = list(reversed(sorted([k for k in list(mean_map[dur2][typ][exp2][dat2].keys())
                                         if len(k.split("_")) == 3], key=lambda s: s.lower())))
        list_dia = [k for k in list_dia if (typ == "bia" and "_pr_" in k) or (typ == "cha" and "_rts_" in k)]
        for ii, dia in enumerate(list_dia):
            list_pro = list(mean_map[dur2][typ][exp2][dat2][dia].keys())
            list_pro = [k for k in list_pro if k not in ["cmip5", "cmip5", "cmip6"]]
            for jj, pro in enumerate(list_pro):
                panel = "panel_11" if typ == list_typ[-1] and dia == list_dia[-1] else "panel_01"
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = copy__deepcopy(dia).replace("_rts_", "_ts_")
                tuple_k = (grp, pan)
                print(typ, dur2, exp2, dat2, dia.ljust(12), pro, tuple_k)
                # panel size
                tk = tuple_k[:-2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, k2, *tk + (k1,))
                # panel axes definition
                if check_dict_levels(fig_ticks, dia, "panel_01"):
                    for k1, k2 in fig_ticks[dia]["panel_01"].items():
                        if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "map_c_nam":
                            tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                            tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                            tmp = tmp.replace(" ()", "")
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
                    val = copy__deepcopy(fig_titles[dia]) if dia in list(fig_titles.keys()) else copy__deepcopy(dia)
                    val = val.replace("VALUE", metadata_netcdf[dia]["name_short"])
                    val = val.replace("UNITS", metadata_netcdf[dia]["units"]).replace(" ()", "")
                    val += ": " + str(pro) + " (shading)"
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    val = fig_titles[typ] if typ in list(fig_titles.keys()) else copy__deepcopy(typ)
                    put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                # -- map
                d1 = mean_map[dur2][typ][exp2][dat2][dia][pro]
                da_s = d1["s"]
                plot_type = "map"
                # shading
                put_in_dict(plot_data, da_s, *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, "bottom", *tuple_k + ("legend_position",))
                # contours
                # da_l = None
                # contour_type = "l2"
                # if check_dict_levels(mean_map, dur3, typ, exp2, dat3, dia, pro3, "l"):
                #     da_l = mean_map[dur3][typ][exp2][dat3][dia][pro3]["l"]
                # elif check_dict_levels(mean_map, dur2, typ, exp2, dat3, dia, pro3, "l"):
                #     da_l = mean_map[dur2][typ][exp2][dat3][dia][pro3]["l"]
                # put_in_dict(plot_data, da_l, *tuple_k + (str(plot_type) + "_" + str(contour_type),))
                # dd = {
                #     "colors": ["k"] * len(dt["map_" + str(contour_type) + "_tic"]),
                #     "linestyles": ["-" if k >= 0 else "--" for k in dt["map_" + str(contour_type) + "_tic"]],
                #     "linewidths": [1 if k != 0 else 2 for k in dt["map_" + str(contour_type) + "_tic"]],
                # }
                # put_in_dict(plot_data, copy__deepcopy(dd), *tuple_k + ("kwargs_" + str(contour_type),))
                da_l = None
                contour_type = "l1"
                if check_dict_levels(mean_map, dur2, typ, exp2, dat2, dia, pro, "l"):
                    da_l = mean_map[dur2][typ][exp2][dat2][dia][pro]["l"]
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
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # plot title
                    # if ii == 0 and jj == 0 and typ == "bia":
                    #     # add to text
                    #     lc += ["k"]
                    #     lf += [19]
                    #     lh += ["center"]
                    #     lr += [0]
                    #     lt += [fig_titles["figure_title"]]
                    #     lv += ["center"]
                    #     lw += ["bold"]
                    #     # maps x = 0 is at the center of the map
                    #     lx += [0]
                    #     ly += [y2 + 20 * dy]
                    # if typ == list_typ[0] and ii == 0 and jj == 0:
                    #     # plot size
                    #     d1, s1, s2 = df[pan]["x_delt"], df[pan]["x_size"], fig_panel_size["panel_21"]["x_size"]
                    #     dx = (x2 - x1) / (100 * s1)
                    #     # add to text
                    #     lc += ["k"]
                    #     lf += [19]
                    #     lh += ["center"]
                    #     lr += [0]
                    #     lt += [fig_titles["figure_title"]]
                    #     lv += ["center"]
                    #     lw += ["bold"]
                    #     # maps x = 0 is at the center of the map, to start from the left of the map, 1/2 map size must
                    #     # be removed to the x coordinate
                    #     lx += [(x1 - x2) / 2 + ((d1 + s1) * len(list_pro) + s2) * 50 * dx]
                    #     ly += [y2 + 20 * dy]
                    if typ == list_typ[0] and ii == 0 and jj == 0:
                        # plot size
                        d1, s1, s2 = df[pan]["x_delt"], df[pan]["x_size"], fig_panel_size["panel_21"]["x_size"]
                        dx = (x2 - x1) / (100 * s1)
                        # add to text
                        lc += ["k"]
                        lf += [18.5]
                        lh += ["right"]
                        lr += [0]
                        lt += [fig_titles["figure_title"]]
                        lv += ["center"]
                        lw += ["bold"]
                        # maps x = 0 is at the center of the map, to start from the left of the map, 1/2 map size must
                        # be removed to the x coordinate
                        lx += [(x1 - x2) / 2 + ((d1 + s1) * len(list_pro) + s2) * 100 * dx]
                        ly += [y2 + 20 * dy]
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
    # panel 21: distribution
    print("distribution")
    list_typ = ["bia", "cha"]
    for typ in list_typ:
        exp2 = "historical" if typ == "bia" else "scenario85"
        if check_dict_levels(mean_box, dur2, typ, exp2, dat2):
            list_dia1 = list(reversed(sorted(list(mean_map[dur2][typ][exp2][dat2].keys()), key=lambda s: s.lower())))
            list_dia1 = [k for k in list_dia1 if (typ == "bia" and "_pr_" in k) or (typ == "cha" and "_ts_" in k)]
            ii = 0
            for dia1 in list_dia1:
                list_dia2 = [k for k in list(mean_box[dur2][typ][exp2][dat2].keys()) if dia1 in k]
                for dia2 in list_dia2:
                    ii += 1
                    panel = "panel_31" if typ == list_typ[-1] and dia1 == list_dia1[-1] and dia2 == list_dia2[-1] else \
                        "panel_21"
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + ii).zfill(2)
                    grp = copy__deepcopy(dia1).replace("_rts_", "_ts_")
                    tuple_k = (grp, pan)
                    print(typ, dur2, exp2, dat2, dia1.ljust(12), dia2.ljust(12), tuple_k)
                    # dictionary
                    d1 = mean_box[dur2][typ][exp2][dat2][dia2]
                    list_d = list(d1.keys())
                    nn = len(list_d) + 1
                    # nn = len(list_d)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_psize, *tk + (pan,)):
                        put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                                put_in_dict(figure_psize, k2, *tk + (k1,))
                    # panel axes definition
                    if check_dict_levels(fig_ticks, dia2, "panel_21"):
                        for k1, k2 in fig_ticks[dia2]["panel_21"].items():
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
                            if isinstance(tmp, list) and k1 in ["x_tic"]:
                                put_in_dict(figure_axes, ["" for _ in range(nn)], *tuple_k + ("x_lab",))
                    dt = {}
                    if check_dict_levels(figure_axes, *tuple_k):
                        dt = read_from_dict(figure_axes, *tuple_k)
                    # test possibility to plot text / legend
                    do_legend = False
                    if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                            dt["y_lim"] is not None:
                        do_legend = True
                    # -- title
                    put_in_dict(plot_data, "CMIP5 vs. CMIP6", *tuple_k + ("title_col",))
                    # -- marker
                    lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                    pro3 = "cmip6"
                    dict_t1 = {}
                    if check_dict_levels(d1, pro3) and dia1 == "ave_pr_val":
                        da1 = copy__deepcopy(d1[pro3])
                        # model nbr with dia_x value
                        dict_t2 = {int(k): float(v) for k, v in zip(da1["dataset"].values, da1.values)}
                        dict_t2 = {k: v for k, v in sorted(dict_t2.items(), key=lambda item: item[1])}
                        list_datasets = list(dict_t2.keys())
                        list_datasets = [k for k in list_datasets if k > 100]
                        nbr = int(round(len(list_datasets) * proc_se / 100, 0))
                        list_d_bot, list_d_top = list_datasets[:nbr], list_datasets[-nbr:]
                        l1 = ["largest", "smallest"] if dia1 == "ave_ts_val" else ["driest", "wettest"]
                        l2 = [list_d_bot, list_d_top]
                        dict_t1 = {grp: da1.sel(indexers={"dataset": l3}) for grp, l3 in zip(l1, l2)}
                    list_grp = list(dict_t1.keys())
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
                    dat3, dur3, exp3 = "multiobs001", "045_year_epoch", "historical"
                    pro3 = "obse"
                    if isinstance(do_legend, bool) and do_legend:
                        # x = 0 and y = 0
                        for n1, k1 in zip(["y"], [dt["y_lim"]]):
                            if k1[0] < 0 < k1[1]:
                                lc += ["k"]
                                ls += [":"]
                                lw += [1]
                                lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                lz += [23]
                        # observations
                        if check_dict_levels(mean_box, dur3, typ, exp3, dat3, dia2, pro3):
                            da1 = copy__deepcopy(mean_box[dur3][typ][exp3][dat3][dia2][pro3])
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
                    lmed, lmed_lw, lo, lq, lq_ls, lq_lw, lw, lx, ly, lz = [], [], [], [], [], [], [], [], [], []
                    for kk, pro in enumerate(list_d):
                        rank = scipy__stats__percentileofscore(d1[pro].values, 0)
                        la += [0.3]
                        lbod_fc += ["k"]
                        lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                        lf += [True]
                        lls += ["-"]
                        llw += [2]
                        lmea += [False]
                        lmea_lw += [0]
                        lmed += [True]
                        lmed_lw += [4]
                        lo += ["vertical"]
                        lq += [[0.25, 0.75]]
                        lq_ls += [":"]
                        lq_lw += [2]
                        lw += [0.8]
                        ly += [d1[pro].values]
                        lx += [kk]
                        lz += [3]
                    for kk, dd in dict_t1.items():
                        la += [0.3]
                        lbod_fc += ["k"]
                        lc += [fig_colors[kk] if kk in list(fig_colors.keys()) else "cyan"]
                        lf += [True]
                        lls += ["-"]
                        llw += [2]
                        lmea += [False]
                        lmea_lw += [0]
                        lmed += [True]
                        lmed_lw += [4]
                        lo += ["vertical"]
                        lq += [[0.25, 0.75]]
                        lq_ls += [":"]
                        lq_lw += [2]
                        lw += [0.8]
                        lx += [len(list_d)]
                        ly += [dd.values]
                        lz += [3]
                    normal = None
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_psize, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        loc, scale = (x1 + x2) / 2, abs(11 * dx)
                        normal = numpy__random__default_rng().normal(loc=loc, scale=scale, size=10000)
                        tmin = scipy__stats__scoreatpercentile(normal, 1)
                        tmax = scipy__stats__scoreatpercentile(normal, 99)
                        normal = numpy__array([k for k in normal if tmin <= k <= tmax])
                        la += [0.35]
                        lbod_fc += ["grey"]
                        lc += ["grey"]
                        lf += [True]
                        lls += ["-"]
                        llw += [1]
                        lmea += [False]
                        lmea_lw += [0]
                        lmed += [True]
                        lmed_lw += [2]
                        lo += ["horizontal"]
                        lq += [[0.25, 0.75]]
                        lq_ls += [":"]
                        lq_lw += [1.5]
                        lw += [abs(6 * dy)]
                        lx += [y1 + 0.05 * (y2 - y1) if dia2 == "ave_ts_val_pgts" else y2 - 0.05 * (y2 - y1)]
                        ly += [normal]
                        lz += [2]
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
                        put_in_dict(plot_data, copy__deepcopy(lo), *tuple_k + (str(plot_type) + "_o",))
                        put_in_dict(plot_data, copy__deepcopy(lq), *tuple_k + (str(plot_type) + "_q",))
                        put_in_dict(plot_data, copy__deepcopy(lq_ls), *tuple_k + (str(plot_type) + "_qua_ls",))
                        put_in_dict(plot_data, copy__deepcopy(lq_lw), *tuple_k + (str(plot_type) + "_qua_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                    # -- annotate
                    if isinstance(do_legend, bool) and do_legend:
                        kk = list_typ.index(typ)
                        df = read_from_dict(figure_psize, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        val = [dict(edgecolor="k", lw=1, facecolor="magenta", shrink=0, headlength=20, headwidth=26,
                                    width=10)]
                        put_in_dict(plot_data, copy__deepcopy(val), *tuple_k + ("ann_arr",))
                        put_in_dict(plot_data, ["magenta"], *tuple_k + ("ann_t_c",))
                        val = [(x1 - (kk * 4 + 24) * dx, (y1 + y2) / 2 + (1.1 - kk * -2) * dy)]
                        put_in_dict(plot_data, copy__deepcopy(val), *tuple_k + ("ann_xy1",))
                        val = [(x1 - (kk * 10 + 46) * dx, (y1 + y2) / 2 + (1.1 - kk * -2) * dy)]
                        put_in_dict(plot_data, copy__deepcopy(val), *tuple_k + ("ann_xy2",))
                        put_in_dict(plot_data, [20], *tuple_k + ("ann_z",))
                    # -- text
                    lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_psize, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # x labels
                        list_leg = list_d + list_grp[:1]
                        for kk, pro in enumerate(list_leg):
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["center"]
                            lr += [0]
                            val = str(int(proc_se)) + "%" if proc_se == int(proc_se) else (
                                "1/3" if proc_se == 100 / 3 else str(proc_se) + "%")
                            lt += [fig_titles[pro].replace("THR%", val)]
                            lv += ["top"]
                            lw += ["regular"]
                            lx += [kk]
                            ly += [y1 - 4 * dy]
                        if len(list_grp) > 1:
                            for kk, pro in enumerate(list_grp[1:]):
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["right"]
                                lr += [0]
                                val = str(int(proc_se)) + "%" if proc_se == int(proc_se) else (
                                    "1/3" if proc_se == 100 / 3 else str(proc_se) + "%")
                                lt += [fig_titles[pro].replace("THR%", val)]
                                lv += ["bottom"]
                                lw += ["regular"]
                                lx += [x2 - 2 * dx]
                                ly += [float(dict_t1[pro].values.max()) + 2 * dy]
                        # observations
                        if check_dict_levels(mean_box, dur3, typ, exp3, dat3, dia2, pro3):
                            da1 = copy__deepcopy(mean_box[dur3][typ][exp3][dat3][dia2][pro3])
                            lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["center"]
                            lr += [0]
                            lt += ["obs."]
                            lv += ["top"]
                            lw += ["regular"]
                            lx += [0.5]
                            ly += [float(da1.min().values) - 2 * dy]
                            for kk, pro in enumerate(list_d):
                                rank = 100 - scipy__stats__percentileofscore(d1[pro].values, float(da1.max().values))
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left" if kk == 0 else "center"]
                                lr += [0]
                                lt += [str(int(round(rank, 0))) + "%"]
                                lv += ["bottom"]
                                lw += ["regular"]
                                lx += [kk + 2 * dx if kk == 0 else kk]
                                ly += [4]
                            lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["left"]
                            lr += [0]
                            lt += ["values > obs"]
                            lv += ["bottom"]
                            lw += ["regular"]
                            lx += [0.0 + 2 * dx]
                            ly += [4.6]
                        if typ == list_typ[-1]:
                            for kk, pro in enumerate(list_d):
                                rank = scipy__stats__percentileofscore(d1[pro].values, 0)
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left" if kk == 0 else "center"]
                                lr += [0]
                                lt += [str(int(round(100 - rank, 0))) + "%"]
                                lv += ["bottom"]
                                lw += ["regular"]
                                lx += [kk + 2 * dx if kk == 0 else kk]
                                ly += [0.07]
                            lc += ["k"]
                            lf += [15]
                            lh += ["left"]
                            lr += [0]
                            lt += ["values > 0"]
                            lv += ["bottom"]
                            lw += ["regular"]
                            lx += [0.0 + 2 * dx]
                            ly += [0.15]
                        # legend distribution
                        if isinstance(normal, numpy__ndarray):
                            lc += ["grey"] * 5
                            lf += [15] * 5
                            lh += ["center"] + ["right"] + ["center"] + ["left"] + ["center"]
                            lr += [0] * 5
                            lt += [str(k) + ("" if k > 1 else "%") for k in [0, 25, 50, 75, 100]]
                            lv += ["bottom" if dia2 == "ave_ts_val_pgts" else "top"] * 5
                            lw += ["regular"] * 5
                            lx += [scipy__stats__scoreatpercentile(normal, k) if 1 < k < 99 else
                                   (float(normal.min()) if k <= 1 else float(normal.max())) for k in
                                   [0, 25, 50, 75, 100]]
                            ly += [y1 + 0.05 * (y2 - y1) + abs(2.2 * dy) if dia2 == "ave_ts_val_pgts" else
                                   y2 - 0.05 * (y2 - y1) - abs(4. * dy)] * 5
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
    # panel 41: scatter
    print("scatter")
    group, panel = "scatter_01", "panel_41"
    if check_dict_levels(mean_bias_ave, dur2) and check_dict_levels(mean_change, dur2):
        l1 = list(reversed(sorted(list(mean_bias_ave[dur2].keys()), key=lambda s: s.lower())))
        l2 = sorted(list(mean_change_ave[dur2].keys()), key=lambda s: s.lower())
        l3 = sorted(list(mean_change[dur2].keys()), key=lambda s: s.lower())
        l3 = [k2 for k1 in ["std_ts_", "std_pr_"] for k2 in l3 if k1 in k2]
        ldia_x = [k1 for k1 in l1 for _ in l3] + [k1 for k1 in l2 for _ in l3]
        ldic_x = [mean_bias_ave[dur2] for _ in l1 for _ in l3] + [mean_change_ave[dur2] for _ in l2 for _ in l3]
        ltyp_x = ["bia" for _ in l1 for _ in l3] + ["cha" for _ in l2 for _ in l3]
        ldia_y = [k2 for _ in l1 for k2 in l3] + [k2 for _ in l2 for k2 in l3]
        ldic_y = [mean_change[dur2] for _ in l1 for _ in l3] + [mean_change[dur2] for _ in l2 for _ in l3]
        ltyp_y = ["cha" for _ in l1 for _ in l3] + ["cha" for _ in l2 for _ in l3]
        for ii, (dia_x, dia_y, dic_x, dic_y, typ_x, typ_y) in enumerate(
                zip(ldia_x, ldia_y, ldic_x, ldic_y, ltyp_x, ltyp_y)):
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + ii % 2).zfill(2)
            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii // 2).zfill(2)
            tuple_k = (grp, pan)
            print(dur2, dia_x.ljust(12), typ_x, dia_y.ljust(12), typ_y, tuple_k)
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, k2, *tk + (k1,))
            # panel axes definition
            dia = "--".join([dia_x, typ_x, dia_y, typ_y])
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
            # -- marker
            leg_col, leg_mar = [], []
            dic_xt, dic_yt = {}, {}
            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
            for pro in list(dic_x[dia_x].keys()):
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
                            if not check_dict_levels(dic_y, dia_y, pro, exp_y, dat):
                                continue
                            da_x = copy__deepcopy(dic_x[dia_x][pro][exp_x][dat])
                            da_y = copy__deepcopy(dic_y[dia_y][pro][exp_y][dat])
                            # model numbers corresponding to cmip5/6
                            list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                 set(da_y["dataset"].values.tolist())))
                            # list_m = [k for k in list_m if k > 100]
                            da_x = da_x.sel(indexers={"dataset": list_m})
                            da_y = da_y.sel(indexers={"dataset": list_m})
                            put_in_dict(dic_xt, da_x, exp_y)
                            put_in_dict(dic_yt, da_y, exp_y)
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
            if ii == 1:
                leg_mar = list(reversed(leg_mar))
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # legend
                if ii % 2 == 0:
                    for kk, mm in enumerate(leg_mar):
                        lc += [fig_colors[mm] if mm in list(fig_colors.keys()) else "cyan"]
                        lm += [fig_markers[mm] if mm in list(fig_markers.keys()) else "o"]
                        ls += [80]
                        lw += [0.5]
                        lx += [x2 + 10 * dx]
                        ly += [(y1 + y2) / 2 + (8 - kk * 8) * dy]
                        lz += [1]
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
            dat3, dur3, exp3, pro3 = "multiobs001", "045_year_epoch", "historical", "obse"
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
                for k1, k2, k3 in zip(["x", "y"], [dia_x, dia_y], [typ_x, typ_y]):
                    if not check_dict_levels(mean_box, dur3, k3, exp3, dat3, k2, pro3) or \
                            (k1 == "x" and typ_x == "chan") or (k1 == "y" and typ_y == "chan"):
                        continue
                    da1 = copy__deepcopy(mean_box[dur3][k3][exp3][dat3][k2][pro3])
                    for da2 in da1:
                        lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                        lo += [True]
                        ls += [fig_linestyles[pro3] if pro3 in list(fig_linestyles.keys()) else ":"]
                        lw += [2]
                        lx += [[da2.values] * 2 if k1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if k1 == "x" else [da2.values] * 2]
                        lz += [22]
                # box
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                lc += ["k"] * 3
                lo += [False] * 3
                ls += ["-"] * 3
                lw += [1] * 3
                if ii == 0:
                    lx += [[x2, x2 + 84 * dx], [x2 + 84 * dx] * 2, [x2, x2 + 84 * dx]]
                    ly += [[y2, y2], [y2, y2 - 22 * dy], [y2 - 22 * dy] * 2]
                else:
                    lx += [[x1, x1 - 84 * dx], [x1 - 84 * dx] * 2, [x1 - 84 * dx, x1]]
                    ly += [[y1, y1], [y1, y1 + 22 * dy], [y1 + 22 * dy] * 2]
                lz += [1] * 3
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
                # linear regression
                # ln = {}
                # for exp1 in list(dic_xt.keys()):
                #     da_x, da_y = dic_xt[exp1], dic_yt[exp1]
                #     # model numbers
                #     list_m = sorted(list(set(da_x["dataset"].values.tolist()) & set(da_y["dataset"].values.tolist())))
                #     for pro1 in ["cmip", "cmip5", "cmip6"]:
                #         # select models
                #         tmp = [k for k in list_m
                #               if pro1 == "cmip" or (pro1 == "cmip5" and k < 100) or (pro1 == "cmip6" and k > 100)]
                #         da_xt, da_yt = da_x.sel(indexers={"dataset": tmp}), da_y.sel(indexers={"dataset": tmp})
                #         # regression
                #         tmp = linear_regression(da_xt.values, da_yt.values)
                #         put_in_dict(ln, tmp, pro1, exp1)
                if check_dict_levels(p_value, dur2) and check_dict_levels(r_value, dur2):
                    list_pro = sorted(list(r_value[dur2].keys()), key=lambda s: s.lower())
                    if ii != 0:
                        list_pro = list(reversed(list_pro))
                    for kk, pro in enumerate(list_pro):
                        if not check_dict_levels(p_value, dur2, pro):
                            continue
                        pp, rr = p_value[dur2][pro], r_value[dur2][pro]
                        tp = [pp[exp]["--".join([dia_x, dia_y])] for exp in list(pp.keys())]
                        tr = float(numpy__array([rr[exp]["--".join([dia_x, dia_y])] for exp in list(rr.keys())]).mean())
                        lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "k"] * 3
                        lf += [15] * 3
                        lh += ["right", "right", "left"]
                        lr += [0] * 3
                        lt += ["r=", "{0:.2f}".format(round(tr, 2)), "".join(["*" if vv < 1e-2 else "" for vv in tp])]
                        lv += ["center"] * 3
                        lw += ["regular"] * 3
                        if ii == 0:
                            lx += [x2 + 12. * dx, x2 + 29.5 * dx, x2 + 30 * dx]
                            ly += [y2 - (kk * 7 + 4.5) * dy] * 3
                        else:
                            lx += [x1 - 75. * dx, x1 - 57.5 * dx, x1 - 57. * dx]
                            ly += [y1 + (kk * 7 + 3.5) * dy] * 3
                        if ii == 0:
                            lc += ["grey"]
                            lf += [15]
                            lh += ["left"]
                            lr += [0]
                            lt += ["correlation"]
                            lv += ["bottom"]
                            lw += ["regular"]
                            lx += [x2 + 3. * dx]
                            ly += [y2 + 1 * dy]
                if check_dict_levels(p_partial, dur2) and check_dict_levels(r_partial, dur2):
                    list_pro = sorted(list(r_partial[dur2].keys()), key=lambda s: s.lower())
                    if ii != 0:
                        list_pro = list(reversed(list_pro))
                    for kk, pro in enumerate(list_pro):
                        if not check_dict_levels(p_partial, dur2, pro):
                            continue
                        pp, rr = p_partial[dur2][pro], r_partial[dur2][pro]
                        tp = [pp[exp]["--".join([dia_x, dia_y])] for exp in list(pp.keys())]
                        tr = float(numpy__array([rr[exp]["--".join([dia_x, dia_y])] for exp in list(rr.keys())]).mean())
                        lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "k"] * 3
                        lf += [15] * 3
                        lh += ["right", "right", "left"]
                        lr += [0] * 3
                        lt += ["r$_p$=", "{0:.2f}".format(round(tr, 2)),
                               "".join(["*" if vv < 1e-2 else "" for vv in tp])]
                        lv += ["center"] * 3
                        lw += ["regular"] * 3
                        if ii == 0:
                            lx += [x2 + 55. * dx, x2 + 72.5 * dx, x2 + 73 * dx]
                            ly += [y2 - (kk * 7 + 4.5) * dy] * 3
                        else:
                            lx += [x1 - 32. * dx, x1 - 14.5 * dx, x1 - 14. * dx]
                            ly += [y1 + (kk * 7 + 3.5) * dy] * 3
                        if ii == 0:
                            lc += ["grey"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            lt += ["partial corr."]
                            lv += ["bottom"]
                            lw += ["regular"]
                            lx += [x2 + 83. * dx]
                            ly += [y2 + 1 * dy]
                # legend
                if ii % 2 == 0:
                    for kk, mm in enumerate(["cmip"] + leg_col + leg_mar):
                        lc += [fig_colors[mm] if mm in list(fig_colors.keys()) else "cyan"]
                        lf += [15]
                        lh += ["right" if mm not in leg_mar else "left"]
                        lr += [0]
                        lt += [fig_titles[mm]]
                        lv += ["center"]
                        lw += ["regular"]
                        lx += [x2 - 2 * dx if mm not in leg_mar else x2 + 14 * dx]
                        ly += [y2 - (kk * 7 + 4.5) * dy if mm not in leg_mar else
                               (y1 + y2) / 2 + (8 - (kk - 3) * 8) * dy]
                # plot title
                if typ_x == "bia":
                    lc += ["k"]
                    lf += [15]
                    lh += ["right"]
                    lr += [0]
                    lt += ["PR bias vs. ENSO change"]
                    lv += ["bottom"]
                    lw += ["bold"]
                    lx += [x2]
                    ly += [y2 + 4 * dy]
                else:
                    lc += ["k"] * 3
                    lf += [15] * 3
                    lh += ["center", "center", "right"]
                    lr += [0] * 3
                    lt += ["warming", "pattern", "vs. ENSO change"]
                    lv += ["bottom"] * 3
                    lw += ["bold"] * 3
                    lx += [x1 + 15 * dx, x1 + 15 * dx, x2 + 10 * dx]
                    ly += [y2 + 11 * dy, y2 + 4 * dy, y2 + 4 * dy]
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
    # panel 51-52: time series
    group = "time_01"
    print("time series")
    dur2, pro2, war2 = "030_year_epoch", "cmip", "raw"
    if check_dict_levels(change_m, dur2):
        ldia_y = list(reversed(sorted(list(change_m[dur2].keys()), key=lambda s: s.lower())))
        ii = -1
        for dia_y in ldia_y:
            ldia_x = list(reversed(sorted(list(change_m[dur2][dia_y].keys()), key=lambda s: s.lower())))
            for dia_x in ldia_x:
                ii += 1
                list_exp = [k for k in data_experiments if k in list(change_m[dur2][dia_y][dia_x].keys())]
                for jj, exp in enumerate(list_exp):
                    tp = "panel_51" if exp == "historical" else "panel_52"
                    pan = str(tp.split("_")[0]) + "_" + str(int(tp.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (grp, pan)
                    print(dur2, dia_y.ljust(12), dia_x.ljust(12), exp.ljust(10), tuple_k)
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
                    l1, l2 = ["x_", "y_"], ["time", dia_y]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dd, tp):
                            for k1, k2 in fig_ticks[dd][tp].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if (k3 == "y_nam" and jj > 0) or (
                                        k3 == "x_nam" and not (dia_y == ldia_y[-1] and dia_x == ldia_x[-1])):
                                    tmp = ""
                                elif k3 == "y_nam":
                                    tmp = tmp.replace("VALUE", metadata_json[dia_y]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dia_y]["units"]).replace(" ()", "")
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
                    pro3 = "unknown"
                    if check_dict_levels(change_m, dur2, dia_y, dia_x, exp):
                        pro3 = list(change_m[dur2][dia_y][dia_x][exp].keys())[0]
                        for grp in list(change_m[dur2][dia_y][dia_x][exp][pro3].keys()):
                            da = copy__deepcopy(change_m[dur2][dia_y][dia_x][exp][pro3][grp])
                            lc += [fig_colors[grp] if grp in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[grp] if grp in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [4]
                    legend = []
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
                            ly += [[y2 - 5 * dy] * 2]
                            lz += [1]
                            legend.append("mean")
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
                    if check_dict_levels(change_r, dur2, dia_y, dia_x, exp, pro3):
                        for grp in list(change_r[dur2][dia_y][dia_x][exp][pro3].keys()):
                            da = copy__deepcopy(change_r[dur2][dia_y][dia_x][exp][pro3][grp])
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
                            la += [0.3]
                            lc += ["k"]
                            lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                            ly1 += [[y2 - (10 + 5 + 3) * dy] * 2]
                            ly2 += [[y2 - (10 + 5 - 3) * dy] * 2]
                            lz += [1]
                            legend.append("range")
                        # reference epoch
                        la += [0.3]
                        lc += ["grey"]
                        lx += [[1956, 2005]]
                        ly1 += [[y1] * 2]
                        ly2 += [[y1 + (y2 - y1) / 6] * 2]
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
                        # reference epoch
                        if jj == 0:
                            lc += ["k"] * 2
                            lf += [14] * 2
                            lh += ["right"] * 2
                            lr += [0] * 2
                            lt += ["reference", "epoch"]
                            lv += ["bottom"] * 2
                            lw += ["regular"] * 2
                            lx += [x2 - dx] * 2
                            ly += [y1 + 0.1 * dy, y1 + 5.5 * dy]
                        # legend
                        if jj == 0:
                            # legend
                            for kk, tt in enumerate(legend):
                                tmp = fig_titles[tt].split("\n")
                                lc += ["k"] * len(tmp)
                                lf += [15] * len(tmp)
                                lh += ["left"] * len(tmp)
                                lr += [0] * len(tmp)
                                lt += copy__deepcopy(tmp)
                                lv += ["center"] if len(tmp) == 1 else ["bottom", "top"]
                                lw += ["regular"] * len(tmp)
                                lx += [x1 + (20 if kk == 0 else 12) * dx] * len(tmp)
                                ly += [y2 - (kk * 10 + 5 + 0.4 * (1 - 2 * k4 if len(tmp) == 2 else 0)) * dy
                                       for k4 in range(len(tmp))]
                            # groups
                            list_grp = list(reversed(sorted(list(set(
                                [grp for grp in list(change_m[dur2][dia_y][dia_x][exp][pro3].keys())])),
                                key=lambda s: s.lower())))
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
                                ly += [y2 - ((k1 + 2.5) * 8 + 5) * dy]
                            if jj == 0:
                                lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left"]
                                lr += [0]
                                lt += [fig_titles[pro3] if pro3 in list(fig_colors.keys()) else pro3]
                                lv += ["bottom"]
                                lw += ["regular"]
                                lx += [x1 + 2 * dx]
                                ly += [0 + 1 * dy]
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
