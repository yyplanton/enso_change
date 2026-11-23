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

# numpy
from numpy.random import randint as numpy__random__randint

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.statistics import numpy__array, linear_regression, linear_regression_slope, mean, median, \
    percentile, standard_deviation
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["ave_ts_val_glob", "ave_ts_val_n34e", "std_ts_ano_n34e"],
    # "data_diagnostics1": ["ave_pr_val_glob", "ave_pr_val_n34e", "std_pr_ano_n34e"],
    # "data_diagnostics1": ["ave_pr_val_glob", "ave_pr_val_n34e", "std_pr_ano_n34e", "ave_ts_val_glob",
    #                       "ave_ts_val_n34e", "std_ts_ano_n34e"],
    "data_diagnostics1": ["ave_ts_val_glob", "ave_pr_val_glob", "ave_ts_val_n34e", "ave_pr_val_n34e",
                          "std_ts_ano_n34e", "std_pr_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": default_arguments["data_epoch_lengths"],
    # list of experiments: list[str]
    "data_experiments": ["piControl", "historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    "smile_rejected": ["EC-EARTH"],
    # "smile_rejected": [],
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
        "panel_11": {"x_delt": 4, "x_size": 8, "y_delt": 1, "y_size": 8},
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
        "cmip5": "--",
        "cmip6": "-.",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "s",
        "cmip6": "D",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_pr_val_glob_change": {
                "panel_01": {
                    # "lim": [-0.05, 0.05],
                    # "nam": "change (UNITS)",
                    # "tic": [round(k / 1000, 3) for k in range(-50, 51, 25)],
                    "lim": [-0.03, 0.09],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 100, 2) for k in range(-3, 10, 3)],
                },
            },
            "ave_pr_val_glob_trend comparison": {
                "panel_11": {
                    "x_lim": [-0.01, 0.01],
                    # "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": [round(k / 1000, 3) for k in range(-10, 11, 5)],
                    "y_lim": [-0.06, 0.06],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 100, 2) for k in range(-6, 7, 3)],
                },
            },
            "ave_ts_val_glob_change": {
                "panel_01": {
                    # "lim": [-0.8, 0.8],
                    # "nam": "change (UNITS)",
                    # "tic": [round(k / 10, 1) for k in range(-8, 9, 4)],
                    "lim": [-0.5, 1.5],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-5, 16, 5)],
                },
            },
            "ave_ts_val_glob_trend comparison": {
                "panel_11": {
                    "x_lim": [-0.16, 0.16],
                    # "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": [round(k / 100, 2) for k in range(-16, 17, 8)],
                    "y_lim": [-0.4, 1.2],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-4, 13, 4)],
                },
            },
            "ave_pr_val_n34e_change": {
                "panel_01": {
                    # "lim": [-0.5, 0.5],
                    # "nam": "change (UNITS)",
                    # "tic": [round(k / 100, 2) for k in range(-50, 51, 25)],
                    "lim": [-0.8, 0.8],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-8, 9, 4)],
                },
            },
            "ave_pr_val_n34e_trend comparison": {
                "panel_11": {
                    "x_lim": [-0.1, 0.1],
                    # "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": [round(k / 100, 2) for k in range(-10, 11, 5)],
                    "y_lim": [-1.0, 1.0],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                },
            },
            "ave_ts_val_n34e_change": {
                "panel_01": {
                    # "lim": [-0.6, 0.6],
                    # "nam": "change (UNITS)",
                    # "tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                    "lim": [-0.8, 0.8],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-8, 9, 4)],
                },
            },
            "ave_ts_val_n34e_trend comparison": {
                "panel_11": {
                    "x_lim": [-0.10, 0.10],
                    # "x_nam": "EXP trend (UNITS.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": [round(k / 100, 2) for k in range(-10, 11, 5)],
                    "y_lim": [-0.3, 0.9],
                    "y_nam": "EXP trend\n(UNITS.century$^{-1}$)",
                    "y_tic": [round(k / 10, 1) for k in range(-3, 10, 3)],
                },
            },
            "std_pr_ano_n34e_change": {
                "panel_01": {
                    "lim": [-70, 70],
                    "nam": "relative change (%)",
                    "tic": list(range(-70, 71, 35)),
                },
            },
            "std_pr_ano_n34e_trend comparison": {
                "panel_11": {
                    "x_lim": [-22, 22],
                    # "x_nam": "EXP trend (%.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": list(range(-22, 23, 11)),
                    "y_lim": [-70, 70],
                    "y_nam": "EXP trend\n(%.century$^{-1}$)",
                    "y_tic": list(range(-70, 71, 35)),
                },
            },
            "std_ts_ano_n34e_change": {
                "panel_01": {
                    "lim": [-40, 40],
                    "nam": "relative change (%)",
                    "tic": list(range(-40, 41, 20)),
                },
            },
            "std_ts_ano_n34e_trend comparison": {
                "panel_11": {
                    "x_lim": [-8, 8],
                    # "x_nam": "EXP trend (%.century$^{-1}$)",
                    "x_nam": "EXP trend",
                    "x_tic": list(range(-8, 9, 4)),
                    "y_lim": [-40, 40],
                    "y_nam": "EXP trend\n(%.century$^{-1}$)",
                    "y_tic": list(range(-40, 41, 20)),
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
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 4,
        "title_row_x": 40,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e57_pi_stationarity_plot(
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
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
    #
    # -- Processing
    #
    # -- rolling mean
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # compute the average over a moving window
    change_mea = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 != data_epoch_lengths[0]:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if exp1 != "piControl":
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        if exp1 == "piControl":
                            da1 = da1.isel(indexers={"year": 0})
                        if exp2 == "piControl":
                            da2 = da2.isel(indexers={"year": 0})
                        elif exp2 == "historical":
                            da2 = da2.sel(indexers={"year": 2006}).mean(dim="member")
                        # reference
                        if exp1 == "piControl":
                            dt = da1["member"].values[1] - da1["member"].values[0]
                        else:
                            dt = da1["year"].values[1] - da1["year"].values[0]
                        dt = int(round(100 / dt, 0))
                        # change
                        da1 -= da1.isel(indexers={"member": list(range(dt))}).mean(dim="member")
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # rolling average
                        da1 = da1.rolling(dim={"member": dt}, min_periods=dt, center=True).mean()
                        # save
                        put_in_dict(change_mea, copy__deepcopy(da1), dur1, exp1, dia, pro, dat)
    # min-max
    print("change min-max")
    for dur in list(change_mea.keys()):
        dt = {}
        for exp in list(change_mea[dur].keys()):
            for dia in list(change_mea[dur][exp].keys()):
                tmp = []
                for pro in list(change_mea[dur][exp][dia].keys()):
                    for da in change_mea[dur][exp][dia][pro].values():
                        tmp += [float(da.min(skipna=True).values), float(da.max(skipna=True).values)]
                if len(tmp) > 0:
                    put_in_dict(dt, tmp, dia, exp)
        for dia in list(dt.keys()):
            for exp in list(dt[dia].keys()):
                print(dur, dia, exp.ljust(12), "{0:.3f}".format(round(min(dt[dia][exp]), 3)).rjust(9),
                      "{0:.3f}".format(round(max(dt[dia][exp]), 3)).rjust(9))
    # -- trend
    # compute the trend
    change_tre = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 != data_epoch_lengths[0]:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # average or select
                        da1 = da1.isel(indexers={"year": 0}) if exp1 == "piControl" else da1.mean(dim="member")
                        if exp2 == "piControl":
                            da2 = da2.isel(indexers={"year": 0})
                        elif exp2 == "historical":
                            da2 = da2.sel(indexers={"year": 2006}).mean(dim="member")
                        # trend
                        if exp1 == "piControl":
                            da1 = float(linear_regression_slope(da1["member"], da1).values) * 100
                        else:
                            da1 = float(linear_regression_slope(da1["year"], da1).values) * 100
                        # relative trend
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / float(da2.values)
                        # save
                        put_in_dict(change_tre, copy__deepcopy(da1), dur1, exp1, dia, pro, dat)
    # outliers
    print("outliers")
    outliers = {}
    for dur in list(change_tre.keys()):
        for exp in list(change_tre[dur].keys()):
            if exp != "piControl":
                continue
            for dia in list(change_tre[dur][exp].keys()):
                # get data
                da = numpy__array([da for pro in list(change_tre[dur][exp][dia].keys())
                                   for dat, da in change_tre[dur][exp][dia][pro].items()])
                # compute mean and IQR (less influenced by outliers)
                mea, med = mean(da), median(da)
                iqr, std = abs(percentile(da, per=75) - percentile(da, per=25)), standard_deviation(da)
                # find outliers
                for pro in list(change_tre[dur][exp][dia].keys()):
                    for dat, da in change_tre[dur][exp][dia][pro].items():
                        tt = abs(da - mea) / std
                        # tt = abs(da - med) / iqr
                        if tt > 4:
                            put_in_dict(outliers, tt, dur, dia, pro, dat)
                            print(dur, dia, pro, dat.ljust(17), "{0:.1f}".format(round(tt, 1)).rjust(4))
    # min-max
    print("trend min-max")
    for dur in list(change_tre.keys()):
        list_outliers = []
        if check_dict_levels(outliers, dur):
            list_outliers = sorted(list(set([
                k3 for k1 in list(outliers[dur].keys()) for k2 in list(outliers[dur][k1].keys())
                for k3 in list(outliers[dur][k1][k2].keys())])), key=lambda s: s.lower())
        dt = {}
        for exp in list(change_tre[dur].keys()):
            for dia in list(change_tre[dur][exp].keys()):
                for pro in list(change_tre[dur][exp][dia].keys()):
                    t1 = []
                    for dat, da in change_tre[dur][exp][dia][pro].items():
                        if dat in list_outliers:
                            continue
                        # if (dia == "ave_pr_val_n34e" and exp == "piControl" and abs(da) > 0.15) or (
                        #         dia == "std_pr_ano_n34e" and
                        #         ((exp == "piControl" and abs(da) > 20) or (exp == "historical" and abs(da) > 70))):
                        #     print(dur, exp.ljust(10), dia, pro, dat.ljust(17),
                        #           "{0:.3f}".format(round(da, 3)).rjust(9))
                        # if (dia == "ave_ts_val_n34e" and exp == "piControl" and abs(da) > 0.1) or (
                        #         dia == "std_pr_ano_n34e" and exp == "piControl" and abs(da) > 11) or (
                        #         dia == "std_ts_ano_n34e" and exp == "piControl" and abs(da) > 7):
                        #     print(dur, exp.ljust(10), dia, pro, dat.ljust(17),
                        #           "{0:.3f}".format(round(da, 3)).rjust(9))
                        t1 += [da]
                    if len(t1) > 0:
                        put_in_dict(dt, t1, dia, exp)
        for dia in list(dt.keys()):
            for exp in list(dt[dia].keys()):
                print(dur, dia, exp.ljust(12), "{0:.3f}".format(round(min(dt[dia][exp]), 3)).rjust(9),
                      "{0:.3f}".format(round(max(dt[dia][exp]), 3)).rjust(9))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group = "group_1"
    for dur in data_epoch_lengths:
        for exp in data_experiments:
            if exp != "piControl":
                continue
            for ii, dia in enumerate(data_diagnostics1):
                l1, l2 = ["change", "trend comparison"], [change_mea, change_tre]
                for jj, (nam, dict_i) in enumerate(zip(l1, l2)):
                    if not check_dict_levels(dict_i, dur, exp, dia):
                        continue
                    panel = "panel_" + str(jj * 10 + 1).zfill(2)
                    pan = copy__deepcopy(panel)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, exp, grp, pan)
                    print(dur, exp.ljust(12), nam.ljust(6), dia, tuple_k)
                    # dictionaries
                    dict_t1 = {}
                    if check_dict_levels(dict_i, dur, exp, dia):
                        dict_t1 = read_from_dict(dict_i, dur, exp, dia)
                    list_pro = sorted(list(dict_t1.keys()), key=lambda s: s.lower())
                    list_outliers = []
                    if check_dict_levels(outliers, dur):
                        list_outliers = sorted(list(set([
                            k3 for k1 in list(outliers[dur].keys()) for k2 in list(outliers[dur][k1].keys())
                            for k3 in list(outliers[dur][k1][k2].keys())])), key=lambda s: s.lower())
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                        put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                put_in_dict(figure_panel_size, k2, *tk + (k1,))
                    # panel axes definition
                    l1 = [""] if nam == "trend comparison" else ["x_", "y_"]
                    l2 = [str(dia) + "_" + str(nam)] if nam == "trend comparison" else \
                        ["time", str(dia) + "_" + str(nam)]
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
                        val = copy__deepcopy(nam) if nam == "trend comparison" else "evolution" + " in " + str(exp)
                        put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, metadata_json[dia]["name_short"], *tuple_k + ("title_row",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    if nam == "change":
                        for pro in list(dict_t1.keys()):
                            for dat, da in dict_t1[pro].items():
                                val = "grey"
                                if check_dict_levels(outliers, dur, dia, pro, dat):
                                    val = fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"
                                lc += [val]
                                val = "-"
                                if check_dict_levels(outliers, dur, dia, pro, dat):
                                    val = fig_linestyles[pro] if pro in list(fig_linestyles.keys()) else "-."
                                ls += [val]
                                lw += [2 if check_dict_levels(outliers, dur, dia, pro, dat) else 1]
                                lx += [da["member"].values - da["member"].values[0]]
                                ly += [da.values]
                                val = numpy__random__randint(2, 9)
                                if check_dict_levels(outliers, dur, dia, pro, dat):
                                    val = 10
                                lz += [val]
                    elif nam == "trend comparison":
                        xx, yy = [], []
                        for pro in list(dict_t1.keys()):
                            for dat, da in dict_t1[pro].items():
                                if not check_dict_levels(dict_i, dur, data_experiments[-1], dia, pro, dat) or \
                                        dat in list_outliers:
                                    continue
                                xx.append(da)
                                yy.append(dict_i[dur][data_experiments[-1]][dia][pro][dat])
                        ln = linear_regression(xx, yy)
                        lc += ["k"]
                        ls += ["-"]
                        lw += [1]
                        lx += [[min(xx), max(xx)]]
                        ly += [[ln["s"] * k + ln["i"] for k in [min(xx), max(xx)]]]
                        lz += [15]
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
                                lz += [11]
                        if nam == "change":
                            t1 = dia.split("_")[-1] == "glob"
                            l3 = copy__deepcopy(list_pro) if t1 else list(reversed(list_pro))
                            for kk, pro in enumerate(l3):
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[pro] if pro in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [[x2 - 10 * dx, x2 - 2 * dx]]
                                ly += [[y2 - (kk * 8 + 5) * dy if t1 else y1 + (kk * 8 + 5) * dy] * 2]
                                lz += [1]
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
                    if nam == "trend comparison":
                        for pro in list(dict_t1.keys()):
                            for dat, da in dict_t1[pro].items():
                                if not check_dict_levels(dict_i, dur, data_experiments[-1], dia, pro, dat) or \
                                        dat in list_outliers:
                                    continue
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                lm += [fig_markers[pro] if pro in list(fig_markers.keys()) else "*"]
                                ls += [60]
                                lw += [0.5]
                                lx += [da]
                                ly += [dict_i[dur][data_experiments[-1]][dia][pro][dat]]
                                lz += [numpy__random__randint(2, 9)]
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
                        # project
                        t1 = nam == "change" and dia.split("_")[-1] == "glob"
                        l3 = copy__deepcopy(list_pro) if t1 else list(reversed(list_pro))
                        for kk, pro in enumerate(l3):
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            lt += [fig_titles[pro]]
                            lv += ["center"]
                            lx += [x2 - 12 * dx if nam == "change" else x2 - 2 * dx]
                            ly += [y2 - (kk * 8 + 5) * dy if t1 else y1 + (kk * 8 + 5) * dy]
                        if nam == "change":
                            lc += ["k"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            lt += ["trend outliers: " + r'$(\frac{x - \bar{x}}{\sigma})$>3']
                            lv += ["center"]
                            lx += [x2 - 40 * dx]
                            ly += [y2 - 9 * dy if t1 else y1 + 9 * dy]
                        # linear regression
                        if nam == "trend comparison":
                            xx, yy = [], []
                            for pro in list(dict_t1.keys()):
                                for dat, da in dict_t1[pro].items():
                                    if not check_dict_levels(dict_i, dur, data_experiments[-1], dia, pro, dat) \
                                            or dat in list_outliers:
                                        continue
                                    xx.append(da)
                                    yy.append(dict_i[dur][data_experiments[-1]][dia][pro][dat])
                            ln = linear_regression(xx, yy)
                            l3 = ["r", "s", "p"]
                            for kk, ll in enumerate(list(reversed(l3))):
                                lc += ["k"] * 2
                                lf += [15] * 2
                                # lh += ["left" if dia.split("_")[0] in ["ave", "ske"] else "right"] * 2
                                lh += ["left"] * 2
                                lr += [0] * 2
                                lt += [str(ll) + " = ", "{0:+.2f}".format(round(ln[ll], 2))]
                                # lv += ["bottom" if dia.split("_")[0] in ["ave", "ske"] else "top"] * 2
                                lv += ["center"] * 2
                                # lx += [x1 + 2 * dx, x1 + 15 * dx] if dia.split("_")[0] in ["ave", "ske"] else \
                                #     [x2 - 22 * dx, x2 - 2 * dx]
                                lx += [x1 + 2 * dx, x1 + 15 * dx]
                                # ly += [y1 + ((len(l3) - 1 - kk) * 7.5 + 1) * dy if dia.split("_")[0] in ["ave", "ske"]
                                #        else y2 - (kk * 7.5 + 2) * dy] * 2
                                ly += [y1 + (kk * 7.5 + 5) * dy] * 2
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
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(smile_minimum_size).zfill(2) + "mem", str(data_diagnostics1[0].split("_")[1])])
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(smile_minimum_size).zfill(2) + "mem"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for exp in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][exp].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, exp)))
            # plot
            fig_basic(plot_data[dur][exp], list_names, 1, figure_axes[dur][exp], fig_format,
                      figure_name_and_path, figure_panel_size[dur][exp], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e57_pi_stationarity_plot()
