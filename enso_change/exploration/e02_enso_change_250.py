# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Time series of ENSO change
# IPCC ch4: Lee et al. (2021) figure 4.10
#     - standard deviation of the detrended Niño 3.4 over 30-year running windows
#     - every single model is normalized by its present-day standard deviation averaged from 1995 to 2014
#     - shading is the 5–95% range
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
from string import ascii_lowercase as string__ascii_lowercase
import sys
from typing import Literal

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import concatenate
from enso_change.wrapper.xarray_tools import statistic, statistic_resampled
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["std_ts_ano_n30e", "std_ts_ano_n34e", "var_ts_ano_n30e", "var_ts_ano_n34e"],
    # "data_diagnostics1": ["var_ts_ano_n30e", "var_ts_ano_n34e"],
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
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
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
        "panel_1": {"x_delt": 4, "x_size": 16, "y_delt": 2, "y_size": 6},
    },
    # colors
    "fig_colors": {
        "all": "k",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "decrease": ":",
        "increase": "--",
        "no_change": "-",
        "relative_change": "-",
        **default_arguments["fig_linestyles_dat"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "std_ts_ano_n30e": {
                "panel_1": {
                    "lim": [-90, 90],
                    "nam": "relative change\nVALUE (%)",
                    "tic": list(range(-90, 91, 45)),
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-90, 90],
                    "nam": "relative change\nVALUE (%)",
                    "tic": list(range(-90, 91, 45)),
                },
            },
            "var_ts_ano_n30e": {
                "panel_1": {
                    "lim": [-80, 240],
                    "nam": "relative change\nVALUE (%)",
                    "tic": list(range(-80, 241, 80)),
                },
            },
            "var_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-90, 270],
                    "nam": "relative change\nVALUE (%)",
                    "tic": list(range(-90, 271, 90)),
                },
            },
            "time_his": {
                "panel_1": {
                    "lim": [1850, 2015],
                    "nam": "Year",
                    "tic": list(range(1860, 2016, 20)),
                },
            },
            "time_ssp": {
                "panel_1": {
                    "lim": [1935, 2100],
                    "nam": "Year",
                    "tic": list(range(1940, 2101, 20)),
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
        "title_col_h": "left",
        "title_col_x": 0,
        "title_col_y": 3,
        "title_row_x": 45,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e02_enso_change_250_plot(
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
    #
    # -- Processing
    #
    # mean and std
    mean_json, _, _ = compute_1var(values_json, "statistic",
                                   function_arguments={"dim": "member", "stat": "mea"})
    sdev_json, _, _ = compute_1var(values_json, "statistic",
                                   function_arguments={"dim": "member", "stat": "std"})
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
                        da2 = da2.sel(indexers={"year": yy})
                        for typ in ["val", "mea"]:
                            da_i = copy__deepcopy(da1) if typ == "val" else da1.mean(dim="member")
                            da_r = copy__deepcopy(da2) if typ == "val" else da2.mean(dim="member")
                            # change
                            da_i -= da_r
                            if dia.split("_")[0] in ["std", "var"]:
                                da_i *= 100 / da_r
                            # centered time
                            da_i = da_i.assign_coords({"year": da_i["year"] - int(dur1.split("_")[0]) / 2})
                            # save
                            put_in_dict(change, copy__deepcopy(da_i), dur1, dia, typ, exp, dat)
    # dataset mean
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    to_plot_ci, to_plot_me, to_plot_r1, to_plot_r2 = {}, {}, {}, {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for typ in list(change[dur][dia].keys()):
                for exp in list(change[dur][dia][typ].keys()):
                    if typ == "mea":
                        # -- based on ensemble mean
                        # array (concatenate datasets)
                        da = concatenate(tuple(da for da in change[dur][dia][typ][exp].values()), "dataset")
                        # mean
                        da_m = da.mean(dim="dataset")
                        # range
                        da_r = statistic(da, dim="dataset", stat="per", per=percentiles)
                        # estimate the confidence interval on the datasets mean
                        da_c = statistic_resampled(
                            da, dim="dataset", pseudo_sample_size=proc_ps, sample_size=da["dataset"].size, stat="mea")
                        da_c = statistic(da_c, dim="resample", stat="per", per=percentiles)
                        # save
                        put_in_dict(to_plot_ci, copy__deepcopy(da_c), dur, dia, exp)
                        put_in_dict(to_plot_me, copy__deepcopy(da_m), dur, dia, exp)
                        put_in_dict(to_plot_r1, copy__deepcopy(da_r), dur, dia, exp)
                    else:
                        # -- based on all members
                        # array (concatenate datasets)
                        da = ()
                        for dd in change[dur][dia][typ][exp].values():
                            # minimum value
                            da += (dd.min(dim="member"),)
                            # maximum value
                            da += (dd.max(dim="member"),)
                        da = concatenate(da, "dataset")
                        # range
                        da_r = statistic(da, dim="dataset", stat="per", per=percentiles)
                        # estimate the confidence interval on the datasets mean
                        # save
                        put_in_dict(to_plot_r2, copy__deepcopy(da_r), dur, dia, exp)
    # min-max
    for dur in list(to_plot_r2.keys()):
        for dia in list(to_plot_r2[dur].keys()):
            tmp = []
            for exp in list(to_plot_r2[dur][dia].keys()):
                tmp += copy__deepcopy(to_plot_r2[dur][dia][exp]).values.flatten().tolist()
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
    for dur in list(to_plot_me.keys()):
        for dia in list(to_plot_me[dur].keys()):
            counter = 0
            for exp in data_experiments:
                if not check_dict_levels(to_plot_me, dur, dia, exp):
                    continue
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + counter % 3).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + counter // 3).zfill(2)
                tuple_k = (dur, dia, grp, pan)
                print(dur, dia, exp.ljust(12), tuple_k)
                # panel size
                tk = tuple_k[:2]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, k2, *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time_" + str(exp[:3]), dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, dd, panel):
                        for k1, k2 in fig_ticks[dur][dd][panel].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and "time" not in dd:
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
                val = str(string__ascii_lowercase[counter]) + ") " + str(fig_titles[exp])
                val += " (" + str(len(list(change[dur][dia]["mea"][exp].keys()))) + " models)"
                put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                lc += ["r"]
                ls += ["--"]
                lw += [1]
                lx += [to_plot_me[dur][dia][exp]["year"].values]
                ly += [to_plot_me[dur][dia][exp].values]
                lz += [5]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # y = 0
                    if y1 < 0 < y2:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [dt["x_lim"]]
                        ly += [[0] * 2]
                        lz += [10]
                    # legend
                    lc += ["r"]
                    ls += ["--"]
                    lw += [1]
                    lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                    ly += [[y1 + (0 * 8 + 5) * dy] * 2]
                    lz += [1]
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
                for i1, (c1, d1) in enumerate(zip(["orange", "blue", "limegreen"], [to_plot_r2, to_plot_r1, to_plot_ci])):
                    if not check_dict_levels(d1, dur, dia, exp):
                        continue
                    la += [1]
                    lc += [c1]
                    lx += [d1[dur][dia][exp]["year"].values]
                    ly1 += [d1[dur][dia][exp].min(dim="quantile").values]
                    ly2 += [d1[dur][dia][exp].max(dim="quantile").values]
                    lz += [i1 + 1]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    for i1, c1 in enumerate(["orange", "blue", "limegreen"]):
                        la += [1]
                        lc += [c1]
                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                        ly1 += [[y2 - (i1 * 8 + 5 + 3) * dy] * 2]
                        ly2 += [[y2 - (i1 * 8 + 5 - 3) * dy] * 2]
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
                lc, lf, lh, lt, lv, lx, ly = [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    lc += ["r"]
                    lf += [15]
                    lh += ["left"]
                    lt += ["MMM"]
                    lv += ["center"]
                    lx += [x1 + 12 * dx]
                    ly += [y1 + (0 * 8 + 5) * dy]
                    l1 = ["orange", "blue", "limegreen"]
                    l2 = ["range 1 member", "range MM", "CI on MMM"]
                    for i1, (c1, k1) in enumerate(zip(l1, l2)):
                        lc += [c1]
                        lf += [15]
                        lh += ["left"]
                        if proc_ci == 90:
                            lt += ["-".join([str(int(round(k2, 0))) for k2 in percentiles]) + "% " + str(k1)]
                        else:
                            lt += [str(int(proc_ci)) + "% " + str(k1)]
                        lv += ["center"]
                        lx += [x1 + 12 * dx]
                        ly += [y2 - (i1 * 8 + 5) * dy]
                    # reference
                    lc += ["k"]
                    lf += [15]
                    lh += ["right"]
                    lt += ["reference epoch: 1850-" + str(1850 + int(dur2.split("_")[0]) - 1)]
                    lv += ["center"]
                    lx += [x2 - 2 * dx]
                    ly += [y1 + (0 * 8 + 5) * dy]
                    # change at the end
                    val = to_plot_me[dur][dia][exp]
                    yy = val["year"].values[-1]
                    val = float(val.sel(indexers={"year": yy}).values)
                    lc += ["r"]
                    lf += [13]
                    lh += ["left"]
                    lt += ["{0:+.0f}".format(val, 0) + "%"]
                    lv += ["center"]
                    lx += [yy + 1 * dx]
                    ly += [val]
                plot_type = "text"
                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                    put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                counter += 1
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(smile_minimum_size).zfill(2) + "mem", str(int(proc_ci)) + "ci"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dia].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia)))
            # plot
            fig_basic(plot_data[dur][dia], list_names, 1, figure_axes[dur][dia], fig_format,
                      figure_name_and_path, figure_panel_size[dur][dia], fig_title_bool=False, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e02_enso_change_250_plot()
