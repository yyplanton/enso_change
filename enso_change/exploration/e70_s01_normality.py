# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure S1: evolution of normality tests
#    * time (x) vs. SSTA std relative change (y) averaged across models
#        - departure from the end of the historical experiment (1965-2014)
#        - MMM-ME for the signal
#        - 1 random member for the spread as in IPCC
#        - normality tests: D'Agostino-Pearson & Shapiro-Wilk
#    * time (x) vs. PRA std relative change (y) averaged across models
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
from enso_change.wrapper.statistics import normality_test
from enso_change.wrapper.xarray_base import concatenate
from enso_change.wrapper.xarray_tools import statistic_resampled
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["std_ts_ano_n34e", "std_pr_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["piControl", "historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {
        "std_pr_ano_n34e": {"hs09--2": [1.07, 1, 5], "ld20": [1, 2, 5]},
        "std_ts_ano_n34e": {"hs09--2": [1.10, 1, 5], "ld20": [1, 2, 5]}
    },
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 1000,
    # threshold for classification
    "proc_th": {"ave": 4, "std": 2, "ske": 2},
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
        "panel_1": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "1m": "peru",
        "em": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "agostino": "--",
        "anderson": ":",
        "shapiro": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "std_pr_ano_n34e": {
                "panel_1": {
                    "lim": [0, 0.8],
                    "nam": "p-value",
                    "tic": [round(k / 10, 1) for k in range(0, 9, 2)],
                },
                "panel_2": {
                    "lim": [0, 0.8],
                    "nam": "p-value",
                    "tic": [round(k / 10, 1) for k in range(0, 9, 2)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [0, 0.8],
                    "nam": "p-value",
                    "tic": [round(k / 10, 1) for k in range(0, 9, 2)],
                },
                "panel_2": {
                    "lim": [0, 0.8],
                    "nam": "p-value",
                    "tic": [round(k / 10, 1) for k in range(0, 9, 2)],
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1860, 2020],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(1850, 2021, 50)),
                },
                "panel_2": {
                    "lim": [1990, 2090],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(2000, 2101, 50)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_pr_val_n34e": "change of VALUE (UNITS)",
        "ske_pr_ano_n34e": "change of VALUE",
        "std_pr_ano_n34e": "relative change of\nVALUE (%)",
        "ave_ts_val_n34e": "change of VALUE (UNITS)",
        "ske_ts_ano_n34e": "change of VALUE",
        "std_ts_ano_n34e": "relative change of\nVALUE (%)",
        "1m": "1 member",
        "em": "member mean",
        "agostino": "D'Agostino-Pearson",
        "anderson": "Anderson-Darling",
        "shapiro": "Shapiro-Wilk",
        "figure_title": "Normality of multimodel distribution",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e70_s01_normality_plot(
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
    proc_co = none_to_default(proc_co, default["proc_co"])
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
                if pro == "observations":
                    continue
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_ensembles = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_ensembles = sorted(list(set(list_ensembles)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_ensembles), exp)
        print(exp.ljust(10), "number of ensembles: " + str(len(list_ensembles)).rjust(2))
    #
    # -- Processing
    #
    # -- compute change compared to initial conditions
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # ensemble mean
    change_em = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if exp1 == "piControl" or not check_dict_levels(values_json, dia, dur2, pro, exp1) or \
                                not check_dict_levels(values_json, dia, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 2014}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # save
                        put_in_dict(change_em, da1, dur1, dia, exp1, pro, dat)
    # one random member
    change_1m = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                if pro == "observations":
                    continue
                for exp in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if exp == "piControl" or not check_dict_levels(values_json, dia, dur2, pro, exp):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp][dat])
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 2014})
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # random selection
                        da1 = statistic_resampled(da1, dim="member", pseudo_sample_size=proc_ps, sample_size=1)
                        # save
                        put_in_dict(change_1m, da1, dur1, dia, exp, pro, dat)
    # -- compute normality tests
    normality = {}
    for k1, d1 in zip(["1m", "em"], [change_1m, change_em]):
        for dur in list(d1.keys()):
            if dur == dur2:
                continue
            for dia in list(d1[dur].keys()):
                for exp in list(d1[dur][dia].keys()):
                    for pro in list(d1[dur][dia][exp].keys()):
                        if exp == "piControl" or pro == "observations":
                            continue
                        # concatenate datasets
                        da = ()
                        for dat in list(d1[dur][dia][exp][pro].keys()):
                            da_t = copy__deepcopy(d1[dur][dia][exp][pro][dat])
                            if "member" in list(da_t.dims):
                                da_t = da_t.mean(dim="member")
                            da += (da_t,)
                        da = concatenate(da, "dataset")
                        # normality
                        for k2 in ["agostino", "shapiro"]:
                            da_o = normality_test(da, dim="dataset", method=k2)
                            if "resample" in list(da_o.dims):
                                da_o = da_o.mean(dim="resample")
                            # save
                            put_in_dict(normality, da_o, dur, dia, exp, pro, k1, k2)
    # min-max
    print("normality min-max")
    for dur in list(normality.keys()):
        for dia in list(normality[dur].keys()):
            tmp = []
            for exp in list(normality[dur][dia].keys()):
                for pro in list(normality[dur][dia][exp].keys()):
                    for k1 in list(normality[dur][dia][exp][pro].keys()):
                        for da in normality[dur][dia][exp][pro][k1].values():
                            tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in list(normality.keys()):
        if not check_dict_levels(normality, dur):
            continue
        list_dia = list(set(list(normality[dur].keys())) & set(list(normality[dur].keys())))
        list_dia = [k for k in data_diagnostics1 if k in list_dia]
        list_exp1 = list(set([k2 for k1 in list_dia for k2 in list(normality[dur][k1].keys())]))
        list_exp2 = list(set([k2 for k1 in list_dia for k2 in list(normality[dur][k1].keys())]))
        list_exp = [k for k in data_experiments if k in list_exp1 + list_exp2]
        for ii, dia in enumerate(list_dia):
            for jj, exp in enumerate(list_exp):
                tp = "panel_1" if exp == "historical" else "panel_2"
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, grp, pan)
                print(dur, dia, exp.ljust(10), tuple_k)
                # panel size
                tk = tuple_k[:1]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    for k1, d1 in fig_panel_size.items():
                        if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k2, d2 in d1.items():
                                if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                    continue
                                tmp = copy__deepcopy(d2)
                                tt = (dur, "time", tp, "lim")
                                if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
                                    k3 = read_from_dict(fig_ticks, *tt)
                                    tmp *= (k3[-1] - k3[0])
                                put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                        elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time", dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, dd, tp):
                        for k1, k2 in fig_ticks[dur][dd][tp].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if k3 == "x_nam" and dia != list_dia[-1]:
                                tmp = ""
                            elif k3 == "y_nam" and jj > 0:
                                tmp = ""
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
                    val = str(fig_titles[exp]) + " (" + str(len(n_ensemble[exp])) + ")"
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                    val = val.replace("UNITS", metadata_json[dia]["units"])
                    put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                if check_dict_levels(normality, dur, dia, exp):
                    for pro in list(normality[dur][dia][exp].keys()):
                        for k1 in list(normality[dur][dia][exp][pro].keys()):
                            for k2, da in normality[dur][dia][exp][pro][k1].items():
                                lc += [fig_colors[k1] if k1 in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[k2] if k2 in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [4]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
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
                    # threshold
                    lc += ["g"]
                    ls += ["-"]
                    lw += [1]
                    lx += [[x1, x2]]
                    ly += [[1 - proc_ci / 100] * 2]
                    lz += [10]
                    # legend
                    if jj == 0 and check_dict_levels(normality, dur, dia, exp):
                        list_test = []
                        for pro in list(normality[dur][dia][exp].keys()):
                            for k1 in list(normality[dur][dia][exp][pro].keys()):
                                list_test += list(normality[dur][dia][exp][pro][k1].keys())
                        list_test = sorted(list(set(list_test)), key=lambda s: s.lower())
                        for k1, tt in enumerate(list_test):
                            lc += ["k"]
                            ls += [fig_linestyles[tt] if tt in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [[x2 - 10 * dx, x2 - 2 * dx]]
                            ly += [[y2 - (k1 * 8 + 5) * dy] * 2]
                            lz += [1]
                plot_type = "cur"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- text
                lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    if jj == 0 and check_dict_levels(normality, dur, dia, exp):
                        list1 = []
                        for pro in list(normality[dur][dia][exp].keys()):
                            for k1 in list(normality[dur][dia][exp][pro].keys()):
                                list1 += list(normality[dur][dia][exp][pro][k1].keys())
                        list1 = sorted(list(set(list1)), key=lambda s: s.lower())
                        for k1, tt in enumerate(list1):
                            lc += ["k"]
                            lf += [15]
                            lh += ["right"]
                            lt += [fig_titles[tt]]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x2 - 12 * dx]
                            ly += [y2 - (k1 * 8 + 5) * dy]
                        list1 = []
                        for pro in list(normality[dur][dia][exp].keys()):
                            list1 += list(normality[dur][dia][exp][pro].keys())
                        list1 = sorted(list(set(list1)), key=lambda s: s.lower())
                        for k1, tt in enumerate(list1):
                            lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["left"]
                            lt += [fig_titles[tt]]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x1 + 2 * dx]
                            ly += [y2 - (k1 * 8 + 15) * dy]
                    # plot title
                    if ii == 0 and jj == 0:
                        # plot size
                        l1 = len([k1 for k1 in list_exp if k1[:3] != "ssp"])
                        x1_d = fig_panel_size[tp]["x_delt"]
                        x1_m = fig_panel_size[tp]["x_size"]
                        x1_s = read_from_dict(fig_ticks, dur, "time", tp, "lim")
                        x1_s = max(x1_s) - min(x1_s)
                        l2 = len([k1 for k1 in list_exp if k1[:3] == "ssp"])
                        x2_d = fig_panel_size["panel_2"]["x_delt"]
                        x2_m = fig_panel_size["panel_2"]["x_size"]
                        x2_s = read_from_dict(fig_ticks, dur, "time", "panel_2", "lim")
                        x2_s = max(x2_s) - min(x2_s)
                        # add to text
                        lc += ["k"]
                        lf += [20]
                        lh += ["center"]
                        lt += [fig_titles["figure_title"]]
                        lv += ["center"]
                        lw += ["bold"]
                        a1 = x1_m * x1_s
                        a2 = (x1_m * x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
                        lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
                        ly += [y2 + 20 * dy]
                plot_type = "text"
                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
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
e70_s01_normality_plot()
