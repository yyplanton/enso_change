# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Relationship between ensemble mean and ensemble standard deviation with climate change
#    * ensemble size (x) vs. slope (y)
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
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.statistics import linear_regression
from enso_change.wrapper.xarray_base import numpy__array, concatenate
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["ave_ts_val_n34e", "std_ts_ano_n34e", "ske_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 5,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
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
        "panel_1": {"x_delt": 4, "x_size": 8, "y_delt": 2, "y_size": 8},
    },
    # colors
    "fig_colors": {
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "ACCESS-CM2": "s",
        "CanESM5-1": "^",
        "CNRM-ESM2-1": "D",
        "GISS-E2-1-G_p1f2": "*",
        "MIROC-ES2L": "X",
        "MPI-ESM1-2-HR": "p",
        "CanESM5_p1": "<",
        "CanESM5_p2": ">",
        "MIROC6": "P",
        "MPI-ESM1-2-LR": "H",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [-0.04, 0.04],
                    "nam": "slope (UNITS)",
                    "tic": [round(k / 100, 2) for k in range(-4, 5, 2)],
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-0.8, 0.8],
                    "nam": "slope (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-8, 9, 4)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-0.4, 0.4],
                    "nam": "slope (%)",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            "ensemble": {
                "panel_1": {
                    "lim": [4, 52],
                    "nam": "ensemble size",
                    "tic": list(range(4, 53, 12)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "VALUE (UNITS)",
        "ske_ts_ano_n34e": "VALUE",
        "std_ts_ano_n34e": "VALUE (%)",
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
def e11_em_vs_es_vs_cc_vs_mem_plot(
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
    # delete datasets if piControl or historical is not available
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                d1 = copy__deepcopy(values_json[dia][dur][pro])
                # list datasets
                list_dat = list(set([dat for exp in list(d1.keys()) for dat in list(d1[exp].keys())]))
                for dat in list_dat:
                    # list experiments
                    list_exp = [exp for exp in list(d1.keys()) if dat in list(d1[exp].keys())]
                    list_ssp = [exp for exp in list_exp if exp[:3] == "ssp"]
                    # delete dataset if piControl or historical is not available
                    if "historical" not in list_exp or len(list_ssp) == 0:
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
                        if len(lm) < smile_minimum_size:
                            continue
                        # select members
                        da_h, da_s = da_h.sel(indexers={"member": lm}), da_s.sel(indexers={"member": lm})
                        # concatenate along year axis
                        da_o = concatenate((da_h, da_s), "year")
                        # resave
                        put_in_dict(values_json, copy__deepcopy(da_o), dia, dur, pro, exp, dat)
                d1 = copy__deepcopy(values_json[dia][dur][pro])
                # list datasets
                list_dat = list(set([dat for exp in list(d1.keys()) for dat in list(d1[exp].keys())]))
                for dat in list_dat:
                    # list experiments
                    list_exp = [exp for exp in list(d1.keys()) if dat in list(d1[exp].keys())]
                    list_ssp = [exp for exp in list_exp if exp[:3] == "ssp"]
                    # delete dataset if piControl or historical is not available
                    if "historical" not in list_exp or len(list_ssp) == 0:
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                    # delete historical
                    del values_json[dia][dur][pro]["historical"][dat]
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
        if exp not in list(n_ensemble.keys()):
            continue
        list_experiments = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(12), "number of ensembles: " + str(len(list_experiments)).rjust(2))
    #
    # -- Processing
    #
    # change
    change = {}
    dur2 = data_epoch_lengths[-1]
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp][dat])
                        # select reference epoch
                        yy = 1850 + int(dur2.split("_")[0])
                        da2 = float(da2.sel(indexers={"year": yy}).mean(dim="member").values)
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(change, copy__deepcopy(da1), dia, dur1, pro, exp, dat)
    # mean and std
    mean_change, _, _ = compute_1var(change, "statistic",
                                   function_arguments={"dim": "member", "stat": "mea"})
    sdev_change, _, _ = compute_1var(change, "statistic",
                                   function_arguments={"dim": "member", "stat": "std"})
    # regression
    regression = {}
    for dia in list(change.keys()):
        for dur in list(change[dia].keys()):
            for pro in list(change[dia][dur].keys()):
                for exp in list(change[dia][dur][pro].keys()):
                    for dat in list(change[dia][dur][pro][exp].keys()):
                        if not check_dict_levels(mean_change, dia, dur, pro, exp, dat) or \
                                not check_dict_levels(sdev_change, dia, dur, pro, exp, dat) or \
                                not check_dict_levels(values_json, dia, dur, pro, exp, dat):
                            continue
                        # arrays
                        da_x = copy__deepcopy(mean_change[dia][dur][pro][exp][dat])
                        da_y = copy__deepcopy(sdev_change[dia][dur][pro][exp][dat])
                        da_n = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        # regression
                        ln = linear_regression(da_x.values, da_y.values)
                        # save
                        if da_n["member"].size >= 10 and dia == data_diagnostics1[0]:
                            print(dur, dia, exp, dat.ljust(17), da_n["member"].size)
                        put_in_dict(regression, float(da_n["member"].size), dur, dia, exp, dat, "x")
                        put_in_dict(regression, float(ln["s"]), dur, dia, exp, dat, "y")
    # min-max
    print("min-max")
    for dur in list(regression.keys()):
        for dia in list(regression[dur].keys()):
            for k in ["x", "y"]:
                tmp = []
                for exp in list(regression[dur][dia].keys()):
                    for dat in list(regression[dur][dia][exp].keys()):
                        if k in list(regression[dur][dia][exp][dat].keys()):
                            tmp.append(regression[dur][dia][exp][dat][k])
                print(dur, dia, k, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: scatterplot
    print("scatterplot")
    group, panel = "group_1", "panel_1"
    for dur in data_epoch_lengths:
        for jj, dia in enumerate(data_diagnostics1):
            if not check_dict_levels(regression, dur, dia):
                continue
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
            grp = copy__deepcopy(group)
            tuple_k = (dur, grp, pan)
            print(dur, dia, tuple_k)
            # panel size
            tk = tuple_k[:1]
            if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                        put_in_dict(figure_panel_size, k2, *tk + (k1,))
            # panel axes definition
            l1, l2 = ["x_", "y_"], ["ensemble", dia]
            for vv, dd in zip(l1, l2):
                if check_dict_levels(fig_ticks, dur, dd, panel):
                    for k1, k2 in fig_ticks[dur][dd][panel].items():
                        k3 = str(vv) + str(k1)
                        if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "nam" and dd not in ["correlation", "ensemble"]:
                            tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                            tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                            tmp = tmp.replace(" ()", "")
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
            val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
            val = val.replace("UNITS",  metadata_json[dia]["units"])
            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                # x = 0 and y = 0
                for k1, k2 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                    if k2[0] < 0 < k2[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if k1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if k1 == "x" else [0] * 2]
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
            lc, llc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], [], []
            for exp in data_experiments:
                if not check_dict_levels(regression, dur, dia, exp):
                    continue
                xx, yy = [], []
                for dat in list(regression[dur][dia][exp].keys()):
                    if check_dict_levels(regression, dur, dia, exp, dat, "x") and \
                            check_dict_levels(regression, dur, dia, exp, dat, "y"):
                        tx = read_from_dict(regression, dur, dia, exp, dat, "x")
                        ty = read_from_dict(regression, dur, dia, exp, dat, "y")
                        xx.append(tx)
                        yy.append(ty)
                        lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "grey"]
                        llc += ["k"]
                        lm += [fig_markers[dat] if dat in list(fig_markers.keys()) else "o"]
                        ls += [40]
                        lw += [0.5]
                        lx += [tx]
                        ly += [ty]
                        lz += [numpy__random__randint(2, 20)]
                # datasets mean for given experiment
                if len(xx) > 0 and len(xx) == len(yy):
                    lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "grey"]
                    llc += ["k"]
                    lm += ["$" + str(exp[3:4]) + "$"]
                    ls += [100]
                    lw += [0.5]
                    lx += [float(numpy__array(xx).mean())]
                    ly += [float(numpy__array(yy).mean())]
                    lz += [21]
            plot_type = "mar"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(llc), *tuple_k + (str(plot_type) + "_lc",))
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
        # list plot names
        list_names = list(plot_data[dur].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur)))
        # plot
        fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
                  figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e11_em_vs_es_vs_cc_vs_mem_plot()
