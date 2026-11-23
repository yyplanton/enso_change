# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: Scatterplots of WEP & EEP vs. zonal gradient (i.e., WEP minus EEP), present, near future, far future
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
from enso_change.wrapper.statistics import numpy__array, linear_regression
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
    "data_diagnostics2": ["ave_ts_val_gxpl", "ave_ts_val_plee", "ave_ts_val_plew"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 60.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co":{},
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
        "frac": {"x": 0.25, "y": 0.25},
        "panel_1": {"x_delt": 4, "x_size": 16, "y_delt": 2, "y_size": 16},
    },
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": {},
    # ticks
    "fig_ticks": {
        "ave_ts_val_gxpl": {
            "panel_1": {
                "lim": [0, 4],
                "nam": "reference epoch\nVALUE (UNITS)",
                "tic": list(range(0, 5, 1)),
            },
        },
        "ave_ts_val_gxpl--present": {
            "panel_1": {
                "lim": [0, 4],
                "nam": "VALUE (UNITS)",
                "tic": list(range(0, 5, 1)),
            },
        },
        "ave_ts_val_gxpl--near": {
            "panel_1": {
                # "lim": [-0.6, 0.6],
                # "nam": "VALUE change (UNITS)",
                # "tic": [round(k / 10, 1) for k in range(-6, 7, 3)],
                "lim": [-0.4, 0.4],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
            },
        },
        "ave_ts_val_gxpl--far": {
            "panel_1": {
                "lim": [-1.5, 0.5],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
            },
        },
        "ave_ts_val_plee": {
            "panel_1": {
                "lim": [24, 28],
                "nam": "reference epoch\nVALUE (UNITS)",
                "tic": list(range(24, 29, 1)),
            },
        },
        "ave_ts_val_plee--present": {
            "panel_1": {
                "lim": [24, 28],
                "nam": "VALUE (UNITS)",
                "tic": list(range(24, 29, 1)),
            },
        },
        "ave_ts_val_plee--near": {
            "panel_1": {
                "lim": [0.0, 1.6],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(0, 17, 4)],
            },
        },
        "ave_ts_val_plee--far": {
            "panel_1": {
                "lim": [1.0, 7.0],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(10, 71, 15)],
            },
        },
        "ave_ts_val_plew": {
            "panel_1": {
                "lim": [27, 31],
                "nam": "reference epoch\nVALUE (UNITS)",
                "tic": list(range(27, 32, 1)),
            },
        },
        "ave_ts_val_plew--present": {
            "panel_1": {
                "lim": [27, 31],
                "nam": "VALUE (UNITS)",
                "tic": list(range(27, 32, 1)),
            },
        },
        "ave_ts_val_plew--near": {
            "panel_1": {
                "lim": [0.2, 1.4],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(2, 15, 3)],
            },
        },
        "ave_ts_val_plew--far": {
            "panel_1": {
                "lim": [1.5, 5.5],
                "nam": "VALUE change (UNITS)",
                "tic": [round(k / 10, 1) for k in range(15, 56, 10)],
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_gxpl": "$\Delta_{x}$ SST (WEP minus EEP)",
        "ave_ts_val_plee": "EEP SST",
        "ave_ts_val_plew": "WEP SST",
        "present": "reference epoch",
        "near": "near future",
        "far": "far future",
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
def z01_alessandra_gradx_plot(
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
    #
    # -- Processing
    #
    # -- select bias year
    mean_bias = {}
    for dia in list(values_netcdf.keys()):
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        # select year
                        da = da.sel(indexers={"year": 1980})
                        # ensemble mean
                        da = float(da.mean(dim="member").values)
                        # save
                        put_in_dict(mean_bias, da, dia, dat)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia in list(values_netcdf.keys()):
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        # epochs
                        list_epochs = sorted([float(k) for k in da["year"].values if k in [1980, 2035, 2100]])
                        # ensemble mean
                        da = da.mean(dim="member")
                        # select reference epoch
                        da_r = da.sel(indexers={"year": list_epochs[0]})
                        # loop on epochs
                        for epo, yy in zip(["present", "near", "far"], list_epochs):
                            # select epoch
                            da_o = da.sel(indexers={"year": yy})
                            # change
                            if yy != list_epochs[0]:
                                da_o -= da_r
                            # save
                            put_in_dict(mean_change, float(da_o.values), dia, epo, dat)
    print("mean_change", list(mean_change.keys()))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_1"
    for dia_y in list(mean_change.keys()):
        list_epo = list(mean_change[dia_y].keys())
        for ii, epo in enumerate(list_epo):
            list_dia_x = list(reversed(sorted(list(mean_bias.keys()), key=lambda s: s.lower())))
            for jj, dia_x in enumerate(list_dia_x):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dia_y, grp, pan)
                print(dia_y, epo.ljust(10), dia_x, tuple_k)
                # panel size
                tk = tuple_k[:-2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, k2, *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], [dia_x, str(dia_y) + "--" + str(epo)]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dd, panel):
                        for k1, k2 in fig_ticks[dd][panel].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd not in ["dataset"]:
                                d2 = dd.split("--")[0]
                                tmp = tmp.replace("VALUE", metadata_json[d2]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[d2]["units"])
                                tmp = tmp.replace(" ()", "")
                                tmp = tmp.replace(" (W21)", "").replace("PLEW", "WEP")
                                tmp = tmp.replace("PLEE", "EEP")
                            if (k3 == "x_nam" and epo != list_epo[-1]) or (k3 == "y_nam" and jj != 0):
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
                    put_in_dict(plot_data, fig_titles[dia_x], *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, fig_titles[epo], *tuple_k + ("title_row",))
                # -- marker
                lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                arr_x, arr_y = [], []
                if check_dict_levels(mean_bias, dia_x) and check_dict_levels(mean_change, dia_y, epo):
                    # list datasets
                    list_dat = list(set(list(mean_bias[dia_x].keys())) & set(list(mean_change[dia_y][epo].keys())))
                    for dat in list_dat:
                        lc += ["r"]
                        lm += ["o"]
                        ls += [80]
                        lw += [0.5]
                        arr_x += [mean_bias[dia_x][dat]]
                        arr_y += [mean_change[dia_y][epo][dat]]
                        lz += [numpy__random__randint(2, 9)]
                lx, ly = copy__deepcopy(arr_x), copy__deepcopy(arr_y)
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
                dict_ln = {}
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    dict_ln = linear_regression(numpy__array(arr_x), numpy__array(arr_y))
                    xx = [x1 + 5 * dx, x2 - 5 * dx]
                    lc += ["k"]
                    ls += ["-"]
                    lw += [2]
                    lx += [xx]
                    ly += [[dict_ln["s"] * k + dict_ln["i"] for k in xx]]
                    lz += [11]
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
                    if "r" in list(dict_ln.keys()) and "p" in list(dict_ln.keys()):
                        l1 = ["r", "p"] if epo == "far" else ["p", "r"]
                        for kk, vv in enumerate(l1):
                            lc += ["k"] * 2
                            lf += [15] * 2
                            lh += ["left", "right"] if epo == "near" else ["right"] * 2
                            lr += [0] * 2
                            lt += [str(vv) + "=", "{0:.3f}".format(round(dict_ln[vv], 3))]
                            lv += ["center"] * 2
                            lw += ["regular"] * 2
                            lx += [x1 + 2 * dx, x1 + 35 * dx] if epo == "near" else [x2 - 25 * dx, x2 - 2 * dx]
                            ly += [y2 - (kk * 8 + 5) * dy if epo == "far" else y1 + (kk * 8 + 5) * dy] * 2
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
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a4) + "rem"])
    # figure loop
    for dia in list(plot_data.keys()):
        # list plot names
        list_names = list(plot_data[dia].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(["", dia]))
        # plot
        fig_basic(plot_data[dia], list_names, 1, figure_axes[dia], fig_format, figure_name_and_path,
                  figure_psize[dia], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
z01_alessandra_gradx_plot()
