# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Evolution of change and percentage of members beyond 2 piControl's internal variability
#    * time (x) vs. change [AVE & SKE] or relative change [STD] (y)
#    * time (x) vs. % of members with value below or above 2 piControl's internal variability (y)
# Both are per experiment averaged across models
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
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import concatenate
from enso_change.wrapper.xarray_tools import statistic
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
    "data_experiments": ["piControl", "historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 10,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": True,
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
        "frac": {"x": 0.005, "y": 0.005},
        "panel_1": {"x_delt": 400, "x_size": 5, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "mean": "-",
        "greater": "-",
        "less": "--",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-0.4, 0.4],
                    "nam": "change (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-15, 45],
                    "nam": "relative change (%)",
                    "tic": list(range(-15, 46, 15)),
                },
            },
            "percentage": {
                "panel_1": {
                    "lim": [0, 100],
                    "nam": "percentage of members (%)",
                    "tic": list(range(0, 101, 25)),
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1850, 2100],
                    "nam": "Year",
                    "tic": list(range(1850, 2101, 50)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "VALUE",
        "ske_ts_ano_n34e": "VALUE",
        "std_ts_ano_n34e": "VALUE",
        "change": "departure from 1850-1899",
        "likelihood": "significant change",
        "iqr": "models' IQR",
        "mean": "MMM",
        "greater": "mem > $+$2 E$_{\sigma}$(PI)",
        "less": "mem < $-$2 E$_{\sigma}$(PI)",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 3,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e12_change_250_pi_iv_plot(
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
        print(exp.ljust(10), "number of ensembles: " + str(len(list_experiments)).rjust(2))
    #
    # -- Processing
    #
    # -- change or relative change
    # compute change compared to initial conditions
    change = {}
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        yy = 1850 + int(dur2.split("_")[0])
                        da2 = float(da2.sel(indexers={"year": yy}).mean(dim="member").values)
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(change, copy__deepcopy(da1), dur1, dia, exp1, dat)
    # mean and std
    mean_change, _, _ = compute_1var(change, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    sdev_change, _, _ = compute_1var(change, "statistic",
                                     function_arguments={"dim": "member", "stat": "std"})
    # multimodel mean and iqr
    change_m, change_r = {}, {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            for exp in list(mean_change[dur][dia].keys()):
                if exp == "piControl":
                    continue
                list_dat = list(mean_change[dur][dia][exp].keys())
                print(mean_change[dur][dia][exp][list_dat[0]].shape, mean_change[dur][dia][exp][list_dat[0]].dims)
                da = concatenate(tuple(da for da in mean_change[dur][dia][exp].values()), "dataset")
                print(da.shape, da.dims)
                stop
                # array (concatenate datasets)
                # da = concatenate(tuple(da for da in mean_change[dur][dia][exp].values()), "dataset")
                # mean
                da_m = statistic(da, dim="dataset", stat="mea")
                # range
                da_r = statistic(da, dim="dataset", stat="per", per=[25, 75])
                # save
                put_in_dict(change_m, copy__deepcopy(da_m), dur, dia, exp, "mean")
                put_in_dict(change_r, copy__deepcopy(da_r), dur, dia, exp, "mean")
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia in list(change_r[dur].keys()):
            tmp = []
            for exp in list(change_r[dur][dia].keys()):
                for da in change_r[dur][dia][exp].values():
                    tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    # -- likelihood of change
    likelihood = {}
    exp2 = "piControl"
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for exp1 in list(change[dur][dia].keys()):
                for dat in list(change[dur][dia][exp1].keys()):
                    if exp1 == "piControl" or not check_dict_levels(sdev_change, dur, dia, exp2, dat):
                        continue
                    # arrays
                    da1 = copy__deepcopy(change[dur][dia][exp1][dat])
                    da2 = copy__deepcopy(sdev_change[dur][dia][exp2][dat])
                    # reference internal variability
                    if "year" in list(da2.dims):
                        da2 = da2.isel(indexers={"year": 0})
                    da2 = float(da2.values)
                    # change greater / less than 2 piControl's IV
                    for k in ["greater", "less"]:
                        # threshold
                        threshold = copy__deepcopy(da2) * (2 if k == "greater" else -2)
                        # test
                        da_t = copy__deepcopy(da1) / threshold
                        da_t = da_t.where(da_t > 1, other=0).where(da_t <= 1, other=100)
                        da_t = da_t.mean(dim="member")
                        # save
                        put_in_dict(likelihood, copy__deepcopy(da_t), dur, dia, exp1, k, dat)
    # multimodel mean and iqr
    likelihood_m, likelihood_r = {}, {}
    for dur in list(likelihood.keys()):
        for dia in list(likelihood[dur].keys()):
            for exp in list(likelihood[dur][dia].keys()):
                if exp == "piControl":
                    continue
                for typ in list(likelihood[dur][dia][exp].keys()):
                    # array (concatenate datasets)
                    da = concatenate(tuple(da for da in likelihood[dur][dia][exp][typ].values()), "dataset")
                    # mean
                    da_m = statistic(da, dim="dataset", stat="mea")
                    # range
                    da_r = statistic(da, dim="dataset", stat="per", per=[25, 75])
                    # save
                    put_in_dict(likelihood_m, copy__deepcopy(da_m), dur, dia, exp, typ)
                    put_in_dict(likelihood_r, copy__deepcopy(da_r), dur, dia, exp, typ)
    # min-max
    print("change min-max")
    for dur in list(likelihood_r.keys()):
        for dia in list(likelihood_r[dur].keys()):
            tmp = []
            for exp in list(likelihood_r[dur][dia].keys()):
                for da in likelihood_r[dur][dia][exp].values():
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
    for dur in data_epoch_lengths:
        ny = int(dur.split("_")[0])
        if dur == dur2:
            continue
        list_n = ["change", "likelihood"]
        list_m, list_r = [change_m, likelihood_m], [change_r, likelihood_r]
        for ii, (name, dict_m, dict_r) in enumerate(zip(list_n, list_m, list_r)):
            for jj, dia in enumerate(data_diagnostics1):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, grp, pan)
                print(dur, name.ljust(15), dia, tuple_k)
                # panel size
                tk = tuple_k[:1]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    for k1, d1 in fig_panel_size.items():
                        if k1[:6] == "panel_" and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k2, d2 in d1.items():
                                if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                    continue
                                tmp = copy__deepcopy(d2)
                                if k2 == "x_size" and check_dict_levels(fig_ticks, dur, "time", panel, "lim"):
                                    k3 = read_from_dict(fig_ticks, dur, "time", panel, "lim")
                                    tmp *= (k3[-1] - k3[0])
                                put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                        elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time", dia if name == "change" else "percentage"]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, dd, panel):
                        for k1, k2 in fig_ticks[dur][dd][panel].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd not in ["percentage", "time"]:
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if k3 == "x_nam" and name != list_n[-1]:
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
                    val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                    val = val.replace("UNITS", metadata_json[dia]["units"])
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, fig_titles[name], *tuple_k + ("title_row",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                if check_dict_levels(dict_m, dur, dia):
                    for exp in list(dict_m[dur][dia].keys()):
                        for typ, da in dict_m[dur][dia][exp].items():
                            if exp[:3] == "ssp":
                                da = da.sel(indexers={"year": [k for k in da["year"].values if k > 2015 - ny / 2]})
                            lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [3]
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
                    # legend
                    l1 = ["mean"] if name == "change" else ["greater", "less"]
                    for k1, typ in enumerate(l1):
                        lc += ["k"]
                        ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else ":"]
                        lw += [2]
                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                        ly += [[y2 - (k1 * 8 + 15) * dy] * 2]
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
                if ii == 0:
                    if check_dict_levels(dict_r, dur, dia):
                        for exp in list(dict_r[dur][dia].keys()):
                            for typ, da in dict_r[dur][dia][exp].items():
                                if exp[:3] == "ssp":
                                    da = da.sel(indexers={"year": [k for k in da["year"].values if k > 2015 - ny / 2]})
                                la += [0.3 if ii < 2 else 1]
                                lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                                lx += [da["year"].values]
                                ly1 += [da.min(dim="quantile").values]
                                ly2 += [da.max(dim="quantile").values]
                                lz += [2]
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    if ii == 0:
                        la += [0.3]
                        lc += ["k"]
                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                        ly1 += [[y2 - (1 * 8 + 15 + 3) * dy] * 2]
                        ly2 += [[y2 - (1 * 8 + 15 - 3) * dy] * 2]
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
                    # curves & shading
                    l1 = ["mean", "iqr"] if name == "change" else ["greater", "less"]
                    for k1, typ in enumerate(l1):
                        lc += ["k"]
                        lf += [15]
                        lh += ["left"]
                        lt += [fig_titles[typ]]
                        lv += ["center"]
                        lx += [x1 + 12 * dx]
                        ly += [y2 - (k1 * 8 + 15) * dy]
                    # colors
                    if ii == 0:
                        if check_dict_levels(dict_m, dur, dia):
                            for ll, (kk, da) in enumerate(dict_m[dur][dia].items()):
                                lc += [fig_colors[kk]]
                                lf += [15]
                                lh += ["left"]
                                lt += [str(fig_titles[kk]) + " (" + str(len(n_ensemble[kk])) + ")"] if ii == 0 else \
                                    [str(kk) + " unc."]
                                lv += ["center"]
                                if ii == 0 and dia.split("_")[0] == "ske":
                                    lx += [x1 + (2 * dx if ll == 0 else (((ll - 1) // 2) * 45 + 2) * dx)]
                                    ly += [y1 + ((2 * 8 + 5) * dy if ll == 0 else ((1 - ((ll - 1) % 2)) * 8 + 5) * dy)]
                                elif ii == 0:
                                    lx += [x1 + (2 * dx if ll == 0 else (((ll - 1) // 2) * 45 + 62) * dx)]
                                    ly += [y1 + (5 * dy if ll == 0 else ((1 - ((ll - 1) % 2)) * 8 + 5) * dy)]
                                else:
                                    lx += [x1 + (ll * 50 + 10) * dx]
                                    ly += [y2 - (0 * 8 + 5) * dy]
                plot_type = "text"
                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
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
e12_change_250_pi_iv_plot()
