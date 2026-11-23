# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Maps of biases: PR, TS, rTS, Tx (values at the end of the 20th century)
#  - Distributions of selected bias
#  - Relationship between selected bias and ENSO SSTA change (departure from the end of the 20th century)
#
# Layout
#
# bias difference
#        ------------ ------------ ------------
#        |    PR    | |    TS    | |          |
#        ------------ ------------ ------------
#          boxplot      scatter
#        ------------ ------------ ------------
#        |          | |          | |          |
#        ------------ ------------ ------------
# time series bottom / top (like Fig. 2)
#        ------------ ------------ ------------
#        |          | |          | |          |
#        ------------ ------------ ------------
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
from enso_change.wrapper.statistics import linear_regression
from enso_change.wrapper.xarray_base import cf_dim_to_dim, concatenate, mean_weighted
from enso_change.wrapper.xarray_tools import select_region
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_glob"],
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
    "smile_rejected": ["EC-EARTH"],
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
    "proc_se": 100 / 3,
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
        "panel_01": {"x_delt": 300, "x_size": 1440, "y_delt": 250, "y_size": 480},
        "panel_11": {"x_delt": 300, "x_size": 1440, "y_delt": 250, "y_size": 480},
        "panel_21": {"x_delt": 320, "x_size": 800, "y_delt": 120, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "cmip5": "green",
        "cmip6": "r",
        "obse": "k",
        "driest": "sienna",
        "wettest": "dodgerblue",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        "obse": "--",
        "driest": "-",
        "wettest": "-",
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
        "driest": "D",
        "wettest": "D",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "ave_pr_val": {
            "panel_01": {
                "x_lim": [120, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "BrBG",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": list(range(-4, 5, 2)),
                "map_l1_tic": list(range(3, 10, 2)),
                "map_l2_tic": list(range(3, 10, 2)),
                "region":["peqspr"],
            },
            "panel_11": {
                "x_lim": [120, 300],
                "x_nam": "",
                "x_tic": [130, 180, 230, 280],
                "y_lim": [-30, 30],
                "y_nam": "",
                "y_tic": list(range(-30, 31, 15)),
                "map_c_cbar": "BrBG_r",
                "map_c_land": "dimgrey",
                "map_c_nam": "",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "map_l1_tic": list(range(3, 10, 2)),
                "map_l2_tic": list(range(3, 10, 2)),
                "region":["peqspr"],
            },
        },
        "ave_pr_val_peqspr--pcc": {
            "panel_21": {
                "x_lim": [0, 6],
                "x_nam": "VALUE (UNITS)",
                "x_tic": [round(k / 10, 1) for k in range(0, 61, 15)],
                "y_lim": [-1, 3],
                "y_nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "y_tic": list(range(-1, 4, 1)),
            },
        },
        "ave_pr_val_peqspr--raw": {
            "panel_21": {
                "x_lim": [0, 6],
                "x_nam": "VALUE (UNITS)",
                "x_tic": [round(k / 10, 1) for k in range(0, 61, 15)],
                "y_lim": [-1.5, 4.5],
                "y_nam": "change VALUE\n(UNITS)",
                "y_tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
            },
        },
    },
    # titles
    "fig_titles": {
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        "pcc": "per degree of CC",
        "raw": "values",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 45,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e65_change_pr_plot(
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
                    if exp1 == "piControl":
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
    # for dur in list(mean_bias.keys()):
    #     if dur != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_bias[dur].keys()):
    #         if dia_x != "ave_taux_val":
    #             continue
    #         for pro in list(mean_bias[dur][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             list_exp = [k for k in list(mean_bias[dur][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat = list(mean_bias[dur][dia_x][pro][list_exp[0]].keys())[0]
    #             dia_y = "std_ts_ano_n34e"
    #             print(dur, dia_x, dia_y, pro, dat, list_exp)
    #             dt = {}
    #             for exp in list_exp:
    #                 if not check_dict_levels(mean_bias, dur, dia_x, pro, exp, dat) or \
    #                         not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat, "pcc"):
    #                     continue
    #                 print(dur, dia_x, dia_y, pro, dat, exp)
    #                 # array
    #                 da = copy__deepcopy(mean_bias[dur][dia_x][pro][exp][dat])
    #                 da_y = copy__deepcopy(mean_change[dur][dia_y][pro][exp][dat]["pcc"])
    #                 # get lat-lon dimensions
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # region average
    #                 for la1 in range(-15, 6, 1):
    #                     for lo1 in range(120, 241, 5):
    #                         for lo2 in range(lo1 + 50, min(lo1 + 150, 280) + 1, 5):
    #                             la2 = la1 + 10
    #                             # mask outside region
    #                             da_x = select_region(da, [la1, la2], [lo1, lo2])
    #                             # horizontal average
    #                             da_x = mean_weighted(da_x, dim=[dim_y, dim_x], skipna=True)
    #                             # correlation with ENSO change
    #                             corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
    #                             # save correlation
    #                             reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2]])
    #                             put_in_dict(dt, [corr], reg)
    #             # rank correlations
    #             dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    #             dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    #             for k in list(dt.keys()):
    #                 la1, la2, lo1, lo2 = int(k[:4]), int(k[4:8]), int(k[8:12]), int(k[12:])
    #                 print(dur, dia_x, dia_y, pro, str(la1).rjust(10) + "N", str(la2).rjust(3) + "N",
    #                       str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E", "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    default_region = {
        # "n30e": {"latitude": [-5., 5.], "longitude": [210., 270.], "short_name": "N3"},
        # "n30s": {"latitude": [-7., 3.], "longitude": [210., 270.], "short_name": "N3[7S-3N]"},
        # "peqntx": {"latitude": [5., 15.], "longitude": [215., 275.], "short_name": "PEQN$_{tx}$"},
        # "peqstx": {"latitude": [-14., -4.], "longitude": [240., 290.], "short_name": "PEQS$_{tx}$"},
        "peqspr": {"latitude": [-7., 3.], "longitude": [190., 280.], "short_name": "PEQS$_{pr}$"},
    }
    list_regions = list(default_region.keys())
    mean_bias_ave = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            if dia not in ["ave_pr_val"]:
                continue
            for pro in list(mean_bias[dur][dia].keys()):
                for exp in list(mean_bias[dur][dia][pro].keys()):
                    for dat in list(mean_bias[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias[dur][dia][pro][exp][dat])
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            da_o = copy__deepcopy(da)
                            # mask outside region
                            da_o = select_region(da_o, default_region[reg]["latitude"],
                                                 default_region[reg]["longitude"])
                            # horizontal average
                            da_o = mean_weighted(da_o, dim=[dim_y, dim_x], skipna=True)
                            # save
                            put_in_dict(mean_bias_ave, da_o, dur, str(dia) + "_" + str(reg), pro, exp)
            if check_dict_levels(mean_bias_ave, dur):
                for dia_o in list(mean_bias_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_bias_ave", list(mean_bias_ave.keys()))
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if dia not in ["ave_pr_val"]:
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        for war in list(mean_change[dur][dia][pro][exp][dat].keys()):
                            # array
                            da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat][war])
                            # get lat-lon dimensions
                            dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                            # loop on regions
                            for reg in list_regions:
                                da_o = copy__deepcopy(da)
                                # mask outside region
                                da_o = select_region(da_o, default_region[reg]["latitude"],
                                                     default_region[reg]["longitude"])
                                # horizontal average
                                da_o = mean_weighted(da_o, dim=[dim_y, dim_x], skipna=True)
                                # save
                                dia_o = str(dia) + "_" + str(reg)
                                put_in_dict(mean_change_ave, da_o, dur, dia_o, pro, exp, war)
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
    dat2, dur2, exp2, pro2 = "multiobs001", "045_year_epoch", "historical", "obse"
    for dur1 in list(mean_bias.keys()):
        for dia in list(mean_bias[dur1].keys()):
            for pro1 in list(mean_bias[dur1][dia].keys()):
                for exp1 in list(mean_bias[dur1][dia][pro1].keys()):
                    if exp1 != "historical":
                        continue
                    for dat1 in list(mean_bias[dur1][dia][pro1][exp1].keys()):
                        if not check_dict_levels(mean_bias, dur2, dia, pro2, exp2, dat2):
                            continue
                        # array
                        da = copy__deepcopy(mean_bias[dur1][dia][pro1][exp1][dat1]).mean(dim="dataset")
                        da_r = copy__deepcopy(mean_bias[dur2][dia][pro2][exp2][dat2])
                        # compute bias
                        if pro1 != pro2:
                            da -= da_r.mean(dim="dataset")
                        # save
                        put_in_dict(mean_map, copy__deepcopy(da), dur1, dia, pro1, "bias")
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                list_dat = sorted(list(set([k2 for k1 in list(mean_change[dur][dia][pro].keys())
                                            for k2 in list(mean_change[dur][dia][pro][k1].keys())])),
                                  key=lambda s: s.lower())
                dat1 = list_dat[0]
                war1 = "pcc"
                list_exp = [k for k in list(mean_change[dur][dia][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
                # concatenate experiments
                da_o = ()
                for exp in list_exp:
                    if not check_dict_levels(mean_change, dur, dia, pro, exp, dat1, war1):
                        continue
                    # array
                    da_o += (mean_change[dur][dia][pro][exp][dat1][war1].mean(dim="dataset"),)
                if len(da_o) == 0:
                    continue
                da_o = concatenate(da_o, dim="experiment").mean(dim="experiment")
                # save
                put_in_dict(mean_map, copy__deepcopy(da_o), dur, dia, pro, "change")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group = "group_1"
    dur2 = "050_year_epoch"
    dur3, pro3 = "045_year_epoch", "obse"
    ii = 0
    for dia in list(mean_map[dur2].keys()):
        if "_pr_" not in dia:
            continue
        for pro in list(mean_map[dur2][dia].keys()):
            if pro != "cmip":
                continue
            for typ in list(mean_map[dur2][dia][pro].keys()):
                panel = "panel_01" if typ == "bias" else "panel_11"
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + ii % 3).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii // 3).zfill(2)
                tuple_k = (grp, pan)
                print(dur2, dia.ljust(12), pro, typ, tuple_k)
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
                        if k1 == "y_tic" and ii % 3 != 0 and isinstance(tmp, list):
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
                val = copy__deepcopy(metadata_netcdf[dia]["name_short"]) + " " + str(typ)
                val += (" per CC" if typ == "change" else "") + " (" + str(metadata_netcdf[dia]["units"])
                val += ("$^\circ$C$^{-1}$" if typ == "change" else "") + ")"
                val = val.replace(" ()", "")
                val = val.replace("TS", "SST")
                put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                # -- map
                da = mean_map[dur2][dia][pro][typ]
                print(da.shape, da.dims)
                plot_type = "map"
                # shading
                put_in_dict(plot_data, copy__deepcopy(da), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, "right", *tuple_k + ("legend_position",))
                # contours
                da = None
                contour_type = "l2"
                if check_dict_levels(mean_map, dur3, dia, pro3, "bias"):
                    da = mean_map[dur3][dia][pro3]["bias"]
                    print(da.shape, da.dims)
                put_in_dict(plot_data, da, *tuple_k + (str(plot_type) + "_" + str(contour_type),))
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
                ii += 1
    # panel 21: scatter
    print("scatter")
    group = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + (ii // 3) + 1)
    panel = "panel_21"
    pro2 = "cmip"
    print(check_dict_levels(mean_bias_ave, dur2), check_dict_levels(mean_change_ave, dur2))
    if check_dict_levels(mean_bias_ave, dur2) and check_dict_levels(mean_change_ave, dur2):
        list_dia = [k for k in list(mean_bias_ave[dur2].keys()) if k in ["ave_pr_val_peqspr"]]
        ii = -1
        for dia in list_dia:
            if not check_dict_levels(mean_bias_ave, dur2, dia, pro2) or \
                    not check_dict_levels(mean_change_ave, dur2, dia, pro2):
                continue
            list_war = list(reversed(sorted(list(set([k2 for k1 in list(mean_change_ave[dur2][dia][pro2].keys())
                                                      for k2 in list(mean_change_ave[dur2][dia][pro2][k1].keys())])),
                                            key=lambda s: s.lower())))
            for war in list_war:
                ii += 1
                list_exp = sorted([k for k in list(mean_change_ave[dur2][dia][pro2].keys())
                                   if k[:3] in ["rcp", "sce", "ssp"]],  key=lambda s: s.lower())
                for jj, exp in enumerate(list_exp):
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (grp, pan)
                    print(dur2, dia.ljust(12), war, exp, tuple_k)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_psize, *tk + (pan,)):
                        put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                                put_in_dict(figure_psize, k2, *tk + (k1,))
                    # panel axes definition
                    dia2 = "--".join([dia, war])
                    if check_dict_levels(fig_ticks, dia2, panel):
                        for k1, k2 in fig_ticks[dia2][panel].items():
                            if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 in ["x_nam", "y_nam"]:
                                tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                                tmp = tmp.replace(" ()", "")
                            if (jj > 0 and k1 == "y_nam") or (k1 == "x_nam" and ii == 0):
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
                        put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, fig_titles[war], *tuple_k + ("title_row",))
                    # -- marker
                    lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                    if check_dict_levels(mean_bias_ave, dur2, dia, pro2, exp) and \
                            check_dict_levels(mean_change_ave, dur2, dia, pro2, exp, war):
                        da_x = copy__deepcopy(mean_bias_ave[dur2][dia][pro2][exp])
                        da_y = copy__deepcopy(mean_change_ave[dur2][dia][pro2][exp][war])
                        # model numbers corresponding to cmip5/6
                        list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                             set(da_y["dataset"].values.tolist())))
                        da_x, da_y = da_x.sel(indexers={"dataset": list_m}), da_y.sel(indexers={"dataset": list_m})
                        for dd, xx, yy in zip(list_m, da_x.values, da_y.values):
                            pro3 = "cmip5" if dd < 100 else "cmip6"
                            lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"]
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
                        if check_dict_levels(mean_bias_ave, dur2, dia, pro2, exp) and \
                                check_dict_levels(mean_change_ave, dur2, dia, pro2, exp, war):
                            da_x = copy__deepcopy(mean_bias_ave[dur2][dia][pro2][exp])
                            da_y = copy__deepcopy(mean_change_ave[dur2][dia][pro2][exp][war])
                            # model numbers corresponding to cmip5/6
                            list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                 set(da_y["dataset"].values.tolist())))
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
                    pro3 = "obse"
                    if isinstance(do_legend, bool) and do_legend and \
                            check_dict_levels(mean_bias_ave, dur3, dia, pro3, exp3):
                        da_x = copy__deepcopy(mean_bias_ave[dur3][dia][pro3][exp3])
                        for da2 in da_x:
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
                    # -- text
                    lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_psize, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # linear regression
                        if check_dict_levels(mean_bias_ave, dur2, dia, pro2, exp) and \
                                check_dict_levels(mean_change_ave, dur2, dia, pro2, exp, war):
                            da_x = copy__deepcopy(mean_bias_ave[dur2][dia][pro2][exp])
                            da_y = copy__deepcopy(mean_change_ave[dur2][dia][pro2][exp][war])
                            # model numbers corresponding to cmip5/6
                            list_m = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                 set(da_y["dataset"].values.tolist())))
                            da_x, da_y = da_x.sel(indexers={"dataset": list_m}), da_y.sel(indexers={"dataset": list_m})
                            ln = linear_regression(da_x.values, da_y.values)
                            lc += ["k"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            val = "***" if ln["p"] < 1e-3 else ("**" if ln["p"] < 1e-2 else (
                                "*" if ln["p"] < 0.05 else ""))
                            lt += ["r=" + "{0:.2f}".format(round(ln["r"], 2)) + str(val)]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x2 - 2 * dx]
                            ly += [y2 - 5 * dy]
                        for k1, pro in enumerate(["cmip6", "cmip5"]):
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            lt += [fig_titles[pro]]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x2 - 2 * dx]
                            ly += [y1 + (k1 * 8 + 5) * dy]
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
e65_change_pr_plot()
