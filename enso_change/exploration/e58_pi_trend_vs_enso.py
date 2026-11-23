# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Check if piControl trends influence historical and scenario trends
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
from enso_change.wrapper.statistics import linear_regression, linear_regression_slope
from enso_change.wrapper.xarray_base import concatenate
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
    "data_diagnostics1": ["ave_pr_val_glob", "ave_ts_val_glob", "ave_pr_val_n34e", "ave_ts_val_n34e",
                          "std_pr_ano_n34e", "std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": default_arguments["data_epoch_lengths"],
    # list of experiments: list[str]
    "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
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
        "frac": {"x": 0.2, "y": 0.2},
        "panel_1": {"x_delt": 6, "x_size": 20, "y_delt": 5, "y_size": 20},
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
        "cmip6": "-",
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
            "ave_pr_val_glob--piControl": {
                "panel_1": {
                    "lim": [-0.01, 0.01],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 1000, 3) for k in range(-10, 11, 5)],
                },
            },
            "ave_pr_val_glob--future": {
                "panel_1": {
                    "lim": [-0.12, 0.36],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-12, 37, 12)],
                },
            },
            "ave_ts_val_glob--piControl": {
                "panel_1": {
                    "lim": [-0.16, 0.16],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-16, 17, 8)],
                },
            },
            "ave_ts_val_glob--future": {
                "panel_1": {
                    "lim": [-2.5, 7.5],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 10, 1) for k in range(-25, 76, 25)],
                },
            },
            "ave_pr_val_n34e--piControl": {
                "panel_1": {
                    "lim": [-0.1, 0.1],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-10, 11, 5)],
                },
            },
            "ave_pr_val_n34e--future": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_n34e--piControl": {
                "panel_1": {
                    "lim": [-0.1, 0.1],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": [round(k / 100, 2) for k in range(-10, 11, 5)],
                },
            },
            "ave_ts_val_n34e--future": {
                "panel_1": {
                    "lim": [-2, 6],
                    "nam": "VALUE trend (UNITS.century$^{-1}$)",
                    "tic": list(range(-2, 7, 2)),
                },
            },
            "std_pr_ano_n34e--piControl": {
                "panel_1": {
                    "lim": [-22, 22],
                    "nam": "VALUE trend (%.century$^{-1}$)",
                    "tic": list(range(-22, 23, 11)),
                },
            },
            "std_pr_ano_n34e--future": {
                "panel_1": {
                    "lim": [-200, 600],
                    "nam": "VALUE trend (%.century$^{-1}$)",
                    "tic": list(range(-200, 601, 200)),
                },
            },
            "std_ts_ano_n34e--piControl": {
                "panel_1": {
                    "lim": [-8, 8],
                    "nam": "VALUE trend (%.century$^{-1}$)",
                    "tic": list(range(-8, 9, 4)),
                },
            },
            "std_ts_ano_n34e--future": {
                "panel_1": {
                    "lim": [-90, 90],
                    "nam": "VALUE trend (%.century$^{-1}$)",
                    "tic": list(range(-90, 91, 45)),
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
def e58_pi_trend_vs_enso_plot(
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
    # concatenate historical and SSPs
    print("concatenate historical and SSPs")
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                if pro == "observations":
                    continue
                # dictionary
                d1 = values_json[dia][dur][pro]
                # list datasets
                list_dataset = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                for dat in list_dataset:
                    # list experiments
                    list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                    list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                    # delete datasets
                    if "historical" not in list_exp:
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                    # concatenate historical and SSPs
                    for exp in list_ssp:
                        # arrays
                        da_h = copy__deepcopy(values_json[dia][dur][pro]["historical"][dat])
                        da_s = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        del values_json[dia][dur][pro][exp][dat]
                        # list members in both historical and ssp
                        lm = sorted(list(set(da_h["member"].values.tolist()) & set(da_s["member"].values.tolist())))
                        # select members
                        da_h, da_s = da_h.sel(indexers={"member": lm}), da_s.sel(indexers={"member": lm})
                        # concatenate along year axis
                        da_o = concatenate((da_h, da_s), "year")
                        # resave
                        put_in_dict(values_json, copy__deepcopy(da_o), dia, dur, pro, exp, dat)
                # delete experiment if no dataset available
                for exp in list(values_json[dia][dur][pro].keys()):
                    if len(list(values_json[dia][dur][pro][exp].keys())) == 0:
                        del values_json[dia][dur][pro][exp]
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
    dur2, exp2 = data_epoch_lengths[-1], "historical"
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
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # average or select
                        da1 = da1.isel(indexers={"year": 0}) if exp1 == "piControl" else da1.mean(dim="member")
                        if exp2 == "piControl":
                            da2 = da2.isel(indexers={"year": 0})
                        elif exp2 == "historical":
                            da2 = da2.sel(indexers={"year": 2006}).mean(dim="member")
                        if exp1[:3] in ["rcp", "ssp"]:
                            da1 = da1.sel(indexers={"year": [k for k in da1["year"].values if k >= 2000]})
                        # trend
                        if exp1 == "piControl":
                            da1 = float(linear_regression_slope(da1["member"], da1).values) * 100
                        else:
                            da1 = float(linear_regression_slope(da1["year"], da1).values) * 100
                        # relative trend
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / float(da2.values)
                        # save
                        exp_o = "scenario" + str(exp1[-2:]) if exp1[:3] in ["rcp", "ssp"] else copy__deepcopy(exp1)
                        put_in_dict(change_tre, copy__deepcopy(da1), dur1, dia, exp_o, pro, dat)
    print("change_tre", list(change_tre.keys()))
    # min-max
    print("trend min-max")
    for dur in list(change_tre.keys()):
        dt = {}
        for dia in list(change_tre[dur].keys()):
            for exp in list(change_tre[dur][dia].keys()):
                for pro in list(change_tre[dur][dia][exp].keys()):
                    t1 = [da for da in change_tre[dur][dia][exp][pro].values()]
                    if len(t1) > 0:
                        exp_o = "future" if exp != "piControl" else copy__deepcopy(exp)
                        put_in_dict(dt, t1, dia, exp_o)
        for dia in list(dt.keys()):
            for exp in list(dt[dia].keys()):
                print(dur, dia, exp.ljust(9), "{0:.3f}".format(round(min(dt[dia][exp]), 3)).rjust(9),
                      "{0:.3f}".format(round(max(dt[dia][exp]), 3)).rjust(9))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_1"
    for dur in list(change_tre.keys()):
        list_dia_x = [k for k in data_diagnostics1 if k in list(change_tre[dur].keys())]
        list_dia_y = [k for k in list_dia_x if k.split("_")[0] in ["std", "var"]]
        list_exp_y = sorted(list(set([k2 for k1 in list_dia_x for k2 in list(change_tre[dur][k1].keys())
                                      if k2 != "piControl"])), key=lambda s: s.lower())
        for dia_y in list_dia_y:
            for ii, dia_x in enumerate(list_dia_x):
                for jj, exp_y in enumerate(list_exp_y):
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, dia_y, grp, pan)
                    print(dur, dia_y, dia_x, exp_y.ljust(10), tuple_k)
                    # dictionaries
                    dix, diy = {}, {}
                    if check_dict_levels(change_tre, dur, dia_x, "piControl"):
                        dix = read_from_dict(change_tre, dur, dia_x, "piControl")
                    if check_dict_levels(change_tre, dur, dia_y, exp_y):
                        diy = read_from_dict(change_tre, dur, dia_y, exp_y)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                        put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                put_in_dict(figure_panel_size, k2, *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], [str(dia_x) + "--piControl", str(dia_y) + "--future"]
                    for vv, dd1 in zip(l1, l2):
                        dd2 = copy__deepcopy(dd1).split("--")[0]
                        if check_dict_levels(fig_ticks, dur, dd1, panel):
                            for k1, k2 in fig_ticks[dur][dd1][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k1 == "nam":
                                    tmp = tmp.replace("VALUE", metadata_json[dd2]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dd2]["units"])
                                    tmp = tmp.replace(" ()", "")
                                if k3 == "y_nam" and jj != 0:
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
                        put_in_dict(plot_data, fig_titles[exp_y], *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, metadata_json[dia_x]["name_short"], *tuple_k + ("title_row",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    xx, yy = [], []
                    for pro in list(dix.keys()):
                        for dat, da in dix[pro].items():
                            if check_dict_levels(diy, pro, dat):
                                xx.append(da)
                                yy.append(diy[pro][dat])
                    ln = linear_regression(xx, yy)
                    lc += ["k"]
                    ls += ["-"]
                    lw += [1]
                    lx += [[min(xx), max(xx)]]
                    ly += [[ln["s"] * k + ln["i"] for k in [min(xx), max(xx)]]]
                    lz += [11]
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # x = 0 and y = 0
                        for k1, k2 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                            if k2[0] < 0 < k2[1]:
                                lc += ["k"]
                                ls += [":"]
                                lw += [1]
                                lx += [[0] * 2 if k1 == "x" else dt["x_lim"]]
                                ly += [dt["y_lim"] if k1 == "x" else [0] * 2]
                                lz += [10]
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
                    for pro in list(dix.keys()):
                        for dat, da in dix[pro].items():
                            if check_dict_levels(diy, pro, dat):
                                lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                                lm += [fig_markers[pro] if pro in list(fig_markers.keys()) else "*"]
                                ls += [60]
                                lw += [0.5]
                                lx += [da]
                                ly += [diy[pro][dat]]
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
                        for kk, pro in enumerate(list(reversed(list(dix.keys())))):
                            lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["left"]
                            lr += [0]
                            lt += [fig_titles[pro]]
                            lv += ["bottom"]
                            lx += [x1 + 2 * dx]
                            ly += [y1 + (kk * 7.5 + 1) * dy]
                        # linear regression
                        l3 = ["r", "s", "p"]
                        for kk, ll in enumerate(l3):
                            lc += ["k"] * 2
                            lf += [15] * 2
                            lh += ["right"] * 2
                            lr += [0] * 2
                            val = "{0:+.2f}".format(round(ln[ll], 2)) if abs(ln[ll]) < 2 else (
                                "{0:+.1f}".format(round(ln[ll], 1)) if abs(ln[ll]) < 10 else (
                                    str(int(round(ln[ll], 0))) if abs(ln[ll]) < 10000 else
                                    "{0:+.2e}".format(round(ln[ll], 2))))
                            lt += [str(ll) + " = ", val]
                            lv += ["top"] * 2
                            lx += [x2 - 23 * dx, x2 - 2 * dx]
                            ly += [y2 - (kk * 7.5 + 1) * dy] * 2
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
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dia].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia)))
            # plot
            fig_basic(plot_data[dur][dia], list_names, 1, figure_axes[dur][dia], fig_format,
                      figure_name_and_path, figure_panel_size[dur][dia], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e58_pi_trend_vs_enso_plot()
