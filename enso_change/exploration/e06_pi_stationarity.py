# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Test if piControl's values are stationary
#    * time series
#    * trends in each model
#    * trends in piControl vs. trend in historical
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

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.statistics import numpy__array, linear_regression, linear_regression_slope, percentile
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["ave_ts_val_n34e", "std_ts_ano_n34e", "ske_ts_ano_n34e"],
    "data_diagnostics1": ["ave_ts_val_glob", "ave_ts_val_n34e", "std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["piControl", "historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    "smile_rejected": [],  # default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": True,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 90.,
    # method: {"members", "percent"}
    "proc_md": "percent",
    # number of pseudo samples: int
    "proc_ps": 100000,
    # threshold for classification
    "proc_th": 10,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.5, "y": 0.5},
        "panel_01": {"x_delt": 4, "x_size": 16, "y_delt": 1, "y_size": 8},
        "panel_11": {"x_delt": 4, "x_size": 18, "y_delt": 1, "y_size": 8},
        "panel_21": {"x_delt": 4, "x_size": 8, "y_delt": 1, "y_size": 8},
    },
    # colors
    "fig_colors": {
        "internal": "limegreen",
        "model": "blue",
        "scenario": "b",
        "CESM2-WACCM-FV2": "r",
        "CMCC-CM2-SR5": "m",
        "CNRM-CM6-1-HR": "y",
        "GFDL-CM4": "limegreen",
        "HadGEM3-GC31-MM": "mediumorchid",
        "CMCC-ESM2": "c",
        # "GFDL-CM4": "limegreen",
        "GFDL-ESM4": "orange",
        "AWI-CM-1-1-MR": "g",
        "E3SM-2-1": "cornflowerblue",
        "NorESM2-LM_p4": "peru",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "internal": "-",
        "model": "-",
        "scenario": "-",
        "CESM2-WACCM-FV2": "-",
        "CMCC-CM2-SR5": "--",
        "CNRM-CM6-1-HR": ":",
        "GFDL-CM4": "-.",
        "HadGEM3-GC31-MM": "-",
        "CMCC-ESM2": "-",
        # "GFDL-CM4": "-.",
        "GFDL-ESM4": "--",
        "AWI-CM-1-1-MR": "-",
        "E3SM-2-1": "--",
        "NorESM2-LM_p4": ":",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "CESM2-WACCM-FV2": "D",
        "CMCC-CM2-SR5": "s",
        "CNRM-CM6-1-HR": "P",
        "GFDL-CM4": ">",
        "HadGEM3-GC31-MM": "^",
        "CMCC-ESM2": "D",
        # "GFDL-CM4": ">",
        "GFDL-ESM4": "s",
        "AWI-CM-1-1-MR": "D",
        "E3SM-2-1": "s",
        "NorESM2-LM_p4": "P",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_glob_change": {
                "panel_01": {
                    "lim": [-1.4, 1.4],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-14, 15, 7)],
                },
            },
            "ave_ts_val_glob_trend": {
                "panel_11": {
                    "lim": [-0.4, 0.4],
                    "nam": "trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            "ave_ts_val_glob_trend comparison": {
                "panel_21": {
                    "x_lim": [-0.4, 0.4],
                    "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                    "y_lim": [-0.4, 1.2],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-4, 13, 4)],
                },
            },
            "ave_ts_val_n34e_change": {
                "panel_01": {
                    "lim": [-1.0, 1.0],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                },
            },
            "ave_ts_val_n34e_trend": {
                "panel_11": {
                    "lim": [-0.26, 0.26],
                    "nam": "trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-26, 27, 13)],
                },
            },
            "ave_ts_val_n34e_trend comparison": {
                "panel_21": {
                    "x_lim": [-0.26, 0.26],
                    "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_tic": [round(k / 100, 2) for k in range(-26, 27, 13)],
                    "y_lim": [-0.3, 0.9],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-3, 10, 3)],
                },
            },
            "ske_ts_ano_n34e_change": {
                "panel_01": {
                    "lim": [-0.8, 0.8],
                    "nam": "change",
                    "tic": [round(k / 10, 1) for k in range(-8, 9, 4)],
                },
            },
            "ske_ts_ano_n34e_trend": {
                "panel_11": {
                    "lim": [-0.2, 0.2],
                    "nam": "trend (century$^{-1}$)",
                    "tic": [round(k / 10, 1) for k in range(-2, 3, 1)],
                },
            },
            "ske_ts_ano_n34e_trend comparison": {
                "panel_21": {
                    "x_lim": [-0.2, 0.2],
                    "x_nam": "EXP trend (century$^{-1}$)",
                    "x_tic": [round(k / 10, 1) for k in range(-2, 3, 1)],
                    "y_lim": [-0.4, 0.4],
                    "y_nam": "EXP trend\n(century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            "std_ts_ano_n34e_change": {
                "panel_01": {
                    "lim": [-100, 100],
                    "nam": "relative change (%)",
                    "tic": list(range(-100, 101, 50)),
                },
            },
            "std_ts_ano_n34e_trend": {
                "panel_11": {
                    "lim": [-10, 10],
                    "nam": "relative trend (%.century$^{-1}$)",
                    "tic": list(range(-10, 11, 5)),
                },
            },
            "std_ts_ano_n34e_trend comparison": {
                "panel_21": {
                    "x_lim": [-10, 10],
                    "x_nam": "EXP relative trend (%.century$^{-1}$)",
                    "x_tic": list(range(-10, 11, 5)),
                    "y_lim": [-50, 50],
                    "y_nam": "EXP relative trend\n(%.century$^{-1}$)",
                    "y_tic": list(range(-50, 51, 25)),
                },
            },
            "dataset": {
                "panel_11": {
                    "lab": None,
                    "lab_rot": 45,
                    "lim": None,
                    "nam": "",
                    "nbr_minor": 0,
                    "tic": None,
                },
            },
            "time": {
                "panel_01": {
                    "lim": [0, 2000],
                    "nam": "epochs",
                    "tic": list(range(0, 2001, 250)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 4,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e06_pi_stationarity_plot(
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
    n_ensemble = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_experiments = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(12), "number of ensembles: " + str(len(list_experiments)).rjust(2))
    for dia in list(values_json.keys()):
        if dia != data_diagnostics1[0]:
            continue
        for dur in list(values_json[dia].keys()):
            if dur != data_epoch_lengths[0]:
                continue
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    list_dat = sorted(list(values_json[dia][dur][pro][exp].keys()), key=lambda s: s.lower())
                    for dat in list_dat:
                        da = values_json[dia][dur][pro][exp][dat]
                        print(dia, dur, pro, exp.ljust(12), dat.ljust(17), da.shape, da.dims)
    #
    # -- Processing
    #
    # -- rolling mean
    # compute the average over a moving window
    change_mea = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    if exp != "piControl":
                        continue
                    for dat in list(values_json[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        if exp == "piControl":
                            da = da.isel(indexers={"year": 0})
                        # reference
                        da_r = da.isel(indexers={"member": list(range(5))}).mean(dim="member")
                        # change
                        da -= da_r
                        if dia.split("_")[0] in ["std", "var"]:
                            da *= 100 / da_r
                        # rolling average
                        da = da.rolling(dim={"member": 5}, min_periods=5, center=True).mean()
                        # save
                        put_in_dict(change_mea, copy__deepcopy(da), dur, pro, exp, dia, dat)
    # min-max
    print("change min-max")
    for dur in list(change_mea.keys()):
        dt = {}
        for pro in list(change_mea[dur].keys()):
            for exp in list(change_mea[dur][pro].keys()):
                for dia in list(change_mea[dur][pro][exp].keys()):
                    tmp = []
                    for da in change_mea[dur][pro][exp][dia].values():
                        tmp += [float(da.min(skipna=True).values), float(da.max(skipna=True).values)]
                    if len(tmp) > 0:
                        put_in_dict(dt, tmp,dia, exp)
        for dia in list(dt.keys()):
            for exp in list(dt[dia].keys()):
                print(dur, dia, exp.ljust(12), "{0:.2f}".format(round(min(dt[dia][exp]), 2)).rjust(8),
                      "{0:.2f}".format(round(max(dt[dia][exp]), 2)).rjust(8))
    # -- trend
    # compute the trend
    change_tre = {}
    exp2 = "piControl"
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp1 in list(values_json[dia][dur][pro].keys()):
                    for dat in list(values_json[dia][dur][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur, pro, exp2, dat):
                            continue
                        dy = int(dur.split("_")[0])
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur][pro][exp1][dat])
                        da2 = copy__deepcopy(values_json[dia][dur][pro][exp2][dat])
                        # average or select
                        da1 = da1.isel(indexers={"year": 0}) if exp1 == "piControl" else da1.mean(dim="member")
                        da2 = da2.isel(indexers={"year": 0})
                        da2 = float(da2.isel(indexers={"member": list(range(5))}).mean(dim="member"))
                        # trend
                        # if exp1 == "piControl":
                        #     da1 = float(linear_regression_slope(da1["member"], da1).values) * 100 / dy
                        # else:
                        #     da1 = float(linear_regression_slope(da1["year"], da1).values) * 100
                        if exp1 == "piControl":
                            da1 = float(linear_regression_slope(da1["member"], da1).values) * 100
                        else:
                            da1 = float(linear_regression_slope(da1["year"], da1).values) * 100
                        # relative trend
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(change_tre, copy__deepcopy(da1), dur, pro, exp1, dia, dat)
    # min-max
    print("trend min-max")
    percentiles_to_compute = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    outliers = {}
    for dur in list(change_tre.keys()):
        dt = {}
        for pro in list(change_tre[dur].keys()):
            for exp in list(change_tre[dur][pro].keys()):
                for dia in list(change_tre[dur][pro][exp].keys()):
                    t1 = {dat: da for dat, da in change_tre[dur][pro][exp][dia].items()}
                    t2 = list(t1.values())
                    if len(t2) > 0:
                        put_in_dict(dt, t2, dia, exp)
                        t3 = percentile(numpy__array(t2), stat="per", per=percentiles_to_compute)
                        a0 = [dat for dat, da in t1.items() if da < t3[0]]
                        a1 = [str(dat) + "(" + "{0:.3f}".format(round(da, 3)).rjust(7) + ")"
                              for dat, da in t1.items() if da < t3[0]]
                        b0 = [dat for dat, da in t1.items() if da > t3[1]]
                        b1 = [str(dat) + "(" + "{0:.3f}".format(round(da, 3)).rjust(7) + ")"
                              for dat, da in t1.items() if da > t3[1]]
                        "{0:.2f}".format(round(min(dt[dia][exp]), 2)).rjust(8)
                        if len(a0) > 0:
                            put_in_dict(outliers, copy__deepcopy(a0), dur, pro, exp, dia, "bottom")
                            print(dur, pro, exp.ljust(10), dia, "outliers", "bot", ", ".join(a1))
                        if len(b0) > 0:
                            put_in_dict(outliers, copy__deepcopy(b0), dur, pro, exp, dia, "top")
                            print(dur, pro, exp.ljust(10), dia, "outliers", "top", ", ".join(b1))
        for dia in list(dt.keys()):
            for exp in list(dt[dia].keys()):
                print(dur, dia, exp.ljust(12), "{0:.2f}".format(round(min(dt[dia][exp]), 2)).rjust(8),
                      "{0:.2f}".format(round(max(dt[dia][exp]), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group = "group_1"
    for dur in data_epoch_lengths:
        for pro in data_projects:
            for exp in data_experiments:
                if exp != "piControl":
                    continue
                for ii, dia in enumerate(data_diagnostics1):
                    l1, l2 = ["change", "trend", "trend comparison"], [change_mea, change_tre, change_tre]
                    for jj, (nam, dict_i) in enumerate(zip(l1, l2)):
                        panel = "panel_" + str(jj * 10 + 1).zfill(2)
                        pan = copy__deepcopy(panel)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = (dur, pro, exp, grp, pan)
                        print(dur, pro, exp.ljust(12), nam.ljust(6), dia, tuple_k)
                        # dictionaries
                        dict_t1 = {}
                        if check_dict_levels(dict_i, dur, pro, exp, dia):
                            dict_t1 = read_from_dict(dict_i, dur, pro, exp, dia)
                        list_dat = sorted(list(dict_t1.keys()), key=lambda s: s.lower())
                        dict_t2 = {}
                        if check_dict_levels(outliers, dur, pro, exp, dia):
                            dict_t2 = read_from_dict(outliers, dur, pro, exp, dia)
                        # panel size
                        tk = tuple_k[:3]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                            for k1, k2 in fig_panel_size.items():
                                if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                    put_in_dict(figure_panel_size, k2, *tk + (k1,))
                        # panel axes definition
                        l1 = [""] if nam == "trend comparison" else ["x_", "y_"]
                        l2 = [str(dia) + "_" + str(nam)] if nam == "trend comparison" else [
                            "time" if nam == "change" else "dataset", str(dia) + "_" + str(nam)]
                        for vv, dd in zip(l1, l2):
                            d2 = copy__deepcopy(dd).split("_" + str(nam))[0]
                            if check_dict_levels(fig_ticks, dur, dd, panel):
                                for k1, k2 in fig_ticks[dur][dd][panel].items():
                                    k3 = str(vv) + str(k1)
                                    if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                        continue
                                    tmp = copy__deepcopy(k2)
                                    if isinstance(tmp, str) and k1.split("_")[-1] == "nam" and \
                                            dd not in ["dataset", "time"]:
                                        tmp = tmp.replace("VALUE", metadata_json[d2]["name_short"])
                                        tmp = tmp.replace("UNITS", metadata_json[d2]["units"])
                                        tmp = tmp.replace(" ()", "")
                                    if k3 == "x_nam" and dia != data_diagnostics1[-1]:
                                        tmp = ""
                                    elif k3 == "x_nam" and nam == "trend comparison":
                                        tmp = tmp.replace("EXP", exp)
                                    elif k3 == "y_nam" and nam == "trend comparison":
                                        tmp = tmp.replace("EXP", data_experiments[-1])
                                    elif dd == "dataset" and k3 in ["x_lab"]:
                                        tmp = ["" for _ in list_dat]
                                    elif dd == "dataset" and k3 in ["x_lim"]:
                                        tmp = [-0.5, len(list_dat) - 0.5]
                                    elif dd == "dataset" and k3 in ["x_tic"]:
                                        tmp = list(range(len(list_dat)))
                                    put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
                        dt = {}
                        if check_dict_levels(figure_axes, *tuple_k):
                            dt = read_from_dict(figure_axes, *tuple_k)
                        # test possibility to plot text / legend
                        do_legend = False
                        if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None \
                                and dt["y_lim"] is not None:
                            do_legend = True
                        # -- title
                        if ii == 0:
                            if nam == "trend comparison":
                                val = copy__deepcopy(nam)
                            else:
                                val = ("evolution" if nam == "change" else "trend") + " in " + str(exp)
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0:
                            put_in_dict(plot_data, metadata_json[dia]["name_short"], *tuple_k + ("title_row",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        if nam == "change":
                            for dat in list_dat:
                                lc += ["dodgerblue" if dat in dat in dict_t2["bottom"] else (
                                    "r" if dat in dat in dict_t2["top"] else "grey")]
                                ls += ["-" if dat in dat in dict_t2["bottom"] + dict_t2["top"] else "--"]
                                lw += [2 if dat in dat in dict_t2["bottom"] + dict_t2["top"] else 1]
                                lx += [dict_t1[dat]["member"].values - dict_t1[dat]["member"].values[0]]
                                ly += [dict_t1[dat].values]
                                lz += [3 if dat in dat in dict_t2["bottom"] + dict_t2["top"] else 2]
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # x = 0 and y = 0
                            for k1, k2 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                                if k2[0] < 0 < k2[1] and (k1 == "y" or (k1 == "x" and nam == "trend comparison")):
                                    lc += ["k"]
                                    ls += [":"]
                                    lw += [1]
                                    lx += [[0] * 2 if k1 == "x" else dt["x_lim"]]
                                    ly += [dt["y_lim"] if k1 == "x" else [0] * 2]
                                    lz += [10]
                            # legend
                            if nam == "change":
                                lc += ["dodgerblue", "r"]
                                ls += ["-", "-"]
                                lw += [2, 2]
                                lx += [[x1 + 2 * dx, x1 + 10 * dx], [x1 + 42 * dx, x1 + 50 * dx]]
                                ly += [[y1 + 5 * dy] * 2] * 2
                                lz += [1] * 2
                            elif nam == "trend comparison":
                                xx, yy = [], []
                                for kk, dat in enumerate(list_dat):
                                    if not check_dict_levels(dict_i, dur, pro, data_experiments[-1], dia, dat):
                                        continue
                                    xx.append(dict_t1[dat])
                                    yy.append(dict_i[dur][pro][data_experiments[-1]][dia][dat])
                                ln = linear_regression(xx, yy)
                                lc += ["k"]
                                ls += ["-"]
                                lw += [1]
                                lx += [[min(xx), max(yy)]]
                                ly += [[ln["s"] * k + ln["i"] for k in [min(xx), max(yy)]]]
                                lz += [20]
                        plot_type = "cur"
                        if len(lx) > 0 and len(lx) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                            put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                            put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                        # -- marker
                        lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                        if nam == "trend":
                            for kk, dat in enumerate(list_dat):
                                lc += ["dodgerblue" if dat in dat in dict_t2["bottom"] else (
                                    "r" if dat in dat in dict_t2["top"] else "grey")]
                                lm += ["D" if dat in dict_t2["bottom"] + dict_t2["top"] else "o"]
                                ls += [60]
                                lw += [0.5]
                                lx += [kk]
                                ly += [dict_t1[dat]]
                                lz += [2]
                        elif nam == "trend comparison":
                            for kk, dat in enumerate(list_dat):
                                if not check_dict_levels(dict_i, dur, pro, data_experiments[-1], dia, dat):
                                    continue
                                lc += ["dodgerblue" if dat in dat in dict_t2["bottom"] else (
                                    "r" if dat in dat in dict_t2["top"] else "grey")]
                                lm += ["D" if dat in dict_t2["bottom"] + dict_t2["top"] else "o"]
                                ls += [60]
                                lw += [0.5]
                                lx += [dict_t1[dat]]
                                ly += [dict_i[dur][pro][data_experiments[-1]][dia][dat]]
                                lz += [2]
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            if nam == "trend":
                                lc += ["dodgerblue", "r"]
                                lm += ["D"] * 2
                                ls += [60] * 2
                                lw += [0.5] * 2
                                lx += [x1 + 6 * dx, x1 + 46 * dx]
                                ly += [y1 + 5 * dy] * 2
                                lz += [1] * 2
                        plot_type = "mar"
                        if len(lx) > 0 and len(lx) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                            put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                            put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                        # -- text
                        lc, lf, lh, lr, lt, lv, lx, ly = [], [], [], [], [], [], [], []
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            if nam in ["change", "trend"]:
                                lc += ["dodgerblue", "r"]
                                lf += [15] * 2
                                lh += ["left"] * 2
                                lr += [0] * 2
                                val = 50 - proc_ci / 2
                                lt += ["bottom", "top " + str(int(val) if int(val) == val else val) + "%"]
                                lv += ["center"] * 2
                                lx += [x1 + 12 * dx, x1 + 52 * dx]
                                ly += [y1 + 5 * dy] * 2
                            elif nam == "trend comparison":
                                xx, yy = [], []
                                for kk, dat in enumerate(list_dat):
                                    if not check_dict_levels(dict_i, dur, pro, data_experiments[-1], dia, dat):
                                        continue
                                    xx.append(dict_t1[dat])
                                    yy.append(dict_i[dur][pro][data_experiments[-1]][dia][dat])
                                ln = linear_regression(xx, yy)
                                l3 = ["r", "s", "p"]
                                for kk, ll in enumerate(l3):
                                    lc += ["k"] * 2
                                    lf += [15] * 2
                                    lh += ["left" if dia.split("_")[0] in ["ave", "ske"] else "right"] * 2
                                    lr += [0] * 2
                                    lt += [str(ll) + " = ", "{0:+.2f}".format(round(ln[ll], 2))]
                                    lv += ["bottom" if dia.split("_")[0] in ["ave", "ske"] else "top"] * 2
                                    lx += [x1 + 2 * dx, x1 + 15 * dx] if dia.split("_")[0] in ["ave", "ske"] else \
                                        [x2 - 22 * dx, x2 - 2 * dx]
                                    ly += [y1 + ((len(l3) - 1 - kk) * 7.5 + 1) * dy if dia.split("_")[0] in ["ave",
                                                                                                             "ske"]
                                           else y2 - (kk * 7.5 + 2) * dy] * 2
                        plot_type = "text"
                        if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                            put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                            put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                            put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(smile_minimum_size).zfill(2) + "mem"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for pro in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for exp in sorted(list(plot_data[dur][pro].keys()), key=lambda s: s.lower()):
                # list plot names
                list_names = list(plot_data[dur][pro][exp].keys())
                # path to output figure
                figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, pro, exp)))
                # plot
                fig_basic(plot_data[dur][pro][exp], list_names, 1, figure_axes[dur][pro][exp], fig_format,
                          figure_name_and_path, figure_panel_size[dur][pro][exp], fig_title_bool=True,
                          panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e06_pi_stationarity_plot()
