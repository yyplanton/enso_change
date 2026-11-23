# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Evolution of change with ensembles (at least 5 members) or with 1 member
#    * time (x) vs. change [AVE & SKE] or relative change [STD] (y) per experiment averaged across models
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
    # method: {"combination", "random"}
    "proc_md": "combination",
    # number of pseudo samples: int
    "proc_ps": 10000,
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
        "panel_1": {"x_delt": 45, "x_size": 3, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 600, "x_size": 3, "y_delt": 100, "y_size": 800},
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
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from 1850-YEAR",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR-2014",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-0.4, 0.4],
                    "nam": "departure from 1850-YEAR",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-0.4, 0.4],
                    "nam": "departure from YEAR-2014",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-15, 45],
                    "nam": "departure from 1850-YEAR",
                    "tic": list(range(-15, 46, 15)),
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-15, 45],
                    "nam": "departure from YEAR-2014",
                    "tic": list(range(-15, 46, 15)),
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1850, 2100],
                    "nam": "Year",
                    "nbr_minor": 9,
                    "tic": list(range(1800, 2100, 100)),
                },
                "panel_2": {
                    "lim": [2000, 2100],
                    "nam": "Year",
                    "nbr_minor": 9,
                    "tic": list(range(2000, 2101, 100)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "change of VALUE (UNITS)",
        "ske_ts_ano_n34e": "change of VALUE",
        "std_ts_ano_n34e": "relative change of VALUE (%)",
        "1m": "1 member change",
        "1m_iqr": "mean models' IQR",
        "1m_mean": "mean",
        "em": "ensemble mean change",
        "em_iqr": "models' IQR",
        "em_mean": "mean",
        "figure_title": "MMM's IQR with or without ensembles",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50 * 1.46,
        "title_col_y": 5,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e15_change_with_iqr_plot(
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
    dur2 = data_epoch_lengths[-1]
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if exp == "piControl" or not check_dict_levels(values_json, dia, dur2, pro, exp):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp][dat])
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        for typ in ["1m", "em"]:
                            for ref in ["250", "100"]:
                                if exp == "piControl" or (exp == "historical" and ref != "250"):
                                    continue
                                # select reference epoch
                                yy = 1850 + int(dur2.split("_")[0]) if ref == "250" else 2015
                                da_r = copy__deepcopy(da2.sel(indexers={"year": yy}))
                                if typ == "em":
                                    da_r = float(da_r.mean(dim="member").values)
                                # change
                                da_o = copy__deepcopy(da1.mean(dim="member") if typ == "em" else da1)
                                da_o -= da_r
                                if dia.split("_")[0] in ["std", "var"]:
                                    da_o *= 100 / da_r
                                # save
                                put_in_dict(change, da_o, dur1, dia, typ, ref, exp, dat)
    # multimodel mean and iqr
    percentiles_to_compute = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_m, change_r = {}, {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for typ in list(change[dur][dia].keys()):
                for ref in list(change[dur][dia][typ].keys()):
                    for exp in list(change[dur][dia][typ][ref].keys()):
                        if exp == "piControl" or (exp == "historical" and ref != "250"):
                            continue
                        # concatenate datasets
                        da = ()
                        for dat in list(change[dur][dia][typ][ref][exp].keys()):
                            # array
                            da_t = copy__deepcopy(change[dur][dia][typ][ref][exp][dat])
                            if typ == "1m":
                                # randomly select one member to generate many possible 1 member multimodel mean
                                da_t = statistic_resampled(
                                    da_t, dim="member", pseudo_sample_size=proc_ps, sample_size=1)
                            # store array
                            da += (da_t,)
                        da = concatenate(da, "dataset")
                        # mean
                        da_m = statistic(da, dim="dataset", stat="mea")
                        if "resample" in list(da_m.dims):
                            da_m = da_m.mean(dim="resample")
                        # range
                        da_r = statistic(da, dim="dataset", stat="per", per=[25, 75])
                        if "resample" in list(da_r.dims):
                            da_r = statistic(da_r, dim="resample", stat="mea")
                            # da1 = statistic(da_r.isel(indexers={"quantile": 0}), dim="resample", stat="per",
                            #                 per=percentiles_to_compute[0])
                            # da2 = statistic(da_r.isel(indexers={"quantile": 1}), dim="resample", stat="per",
                            #                 per=percentiles_to_compute[-1])
                            # da_r = concatenate((da1, da2), "quantile")
                        # save
                        put_in_dict(change_m, copy__deepcopy(da_m), dur, dia, typ, ref, exp)
                        put_in_dict(change_r, copy__deepcopy(da_r), dur, dia, typ, ref, exp)
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia in list(change_r[dur].keys()):
            tmp = []
            for typ in list(change_r[dur][dia].keys()):
                for ref in list(change_r[dur][dia][typ].keys()):
                    for da in change_r[dur][dia][typ][ref].values():
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
    for dur in list(change_m.keys()):
        ny1 = int(dur.split("_")[0])
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        # list variables ('pr', 'ts', ...)
        list_var = sorted(list(set([k.split("_")[1] for k in list(change_m[dur].keys())])), key=lambda s: s.lower())
        for var in list_var:
            # list diagnostics
            list_dia = [k for k in list(change_m[dur].keys()) if k.split("_")[1] == var]
            # list figure types ('1m', 'em')
            list_typ = sorted(list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())])),
                              key=lambda s: s.lower())
            for ii, typ in enumerate(list_typ):
                cc = 0
                for jj, dia in enumerate(list_dia):
                    for kk, ref in enumerate(list(change_m[dur][dia][typ].keys())):
                        tp = "panel_1" if ref == "250" else "panel_2"
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + cc).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = (dur, var, grp, pan)
                        print(dur, var, typ, dia, ref, tuple_k)
                        # panel size
                        tk = tuple_k[:2]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k1, d1 in fig_panel_size.items():
                                if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                                    for k2, d2 in d1.items():
                                        if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                            continue
                                        tmp = copy__deepcopy(d2)
                                        if k2 == "x_size" and \
                                                check_dict_levels(fig_ticks, dur, "time", tp, "lim"):
                                            k3 = read_from_dict(fig_ticks, dur, "time", tp, "lim")
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
                                        if ref == "250":
                                            tmp = tmp.replace("YEAR", str(1850 + ny2 - 1))
                                        else:
                                            tmp = tmp.replace("YEAR", str(2015 - ny2))
                                        tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                        tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                        tmp = tmp.replace(" ()", "")
                                    if k3 == "x_nam" and typ != list_typ[-1]:
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
                        if ii == 0 and kk == 0:
                            val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                            val = val.replace("UNITS", metadata_json[dia]["units"])
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0 and kk == 0:
                            put_in_dict(plot_data, fig_titles[typ], *tuple_k + ("title_row",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        if check_dict_levels(change_m, dur, dia, typ, ref):
                            for exp, da in change_m[dur][dia][typ][ref].items():
                                if exp[:3] == "ssp":
                                    da = da.sel(
                                        indexers={"year": [k for k in da["year"].values if k >= 2015 - ny1 / 2]})
                                lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[exp] if exp in list(fig_linestyles.keys()) else ":"]
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
                            if kk == 0:
                                lc += ["k"]
                                ls += ["-"]
                                lw += [2]
                                lx += [[x1 + 10 * dx, x1 + 18 * dx]]
                                ly += [[y2 - 5 * dy] * 2]
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
                        if check_dict_levels(change_r, dur, dia, typ, ref):
                            for exp, da in change_r[dur][dia][typ][ref].items():
                                if exp[:3] == "ssp":
                                    da = da.sel(
                                        indexers={"year": [k for k in da["year"].values if k >= 2015 - ny1 / 2]})
                                la += [0.3]
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
                            if kk == 0:
                                la += [0.3]
                                lc += ["k"]
                                lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                                ly1 += [[y2 - (15 + 3) * dy] * 2]
                                ly2 += [[y2 - (15 - 3) * dy] * 2]
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
                        lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            if kk == 0:
                                # curve and shading
                                for k1, v1 in enumerate(["mean", "iqr"]):
                                    lc += ["k"]
                                    lf += [15]
                                    lh += ["left"]
                                    lt += [fig_titles[str(typ) + "_" + str(v1)]]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x1 + (20 if k1 == 0 else 12) * dx]
                                    ly += [y2 - (5 if k1 == 0 else 15) * dy]
                                # colors
                                if check_dict_levels(change_m, dur, dia, typ, ref):
                                    for ll, exp in enumerate(change_m[dur][dia][typ][ref].keys()):
                                        lc += [fig_colors[exp]]
                                        lf += [15]
                                        lh += ["left"]
                                        lt += [str(fig_titles[exp]) + " (" + str(len(n_ensemble[exp])) + ")"]
                                        lv += ["center"]
                                        lw += ["regular"]
                                        if dia.split("_")[0] == "ske":
                                            lx += [x1 + (2 if ll == 0 else (((ll - 1) // 2) * 45 + 2)) * dx]
                                            ly += [y1 + (
                                                (2 * 8 + 5) * dy if ll == 0 else ((1 - ((ll - 1) % 2)) * 8 + 5) * dy)]
                                        else:
                                            lx += [x1 + (2 if ll < 3 else ((ll - 3) * 50 + 2)) * dx]
                                            ly += [y2 - ((ll + 1) * 8 + 15) * dy if ll < 3 else y1 + 5 * dy]
                            # plot title
                            if ii == 0 and jj == 0 and kk == 0:
                                # plot size
                                x1_d = fig_panel_size[tp]["x_delt"]
                                x1_m = fig_panel_size[tp]["x_size"]
                                x1_s = read_from_dict(fig_ticks, dur, "time", tp, "lim")
                                x1_s = max(x1_s) - min(x1_s)
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
                                a2 = (a1 + x1_d + x2_m * x2_s) * len(list_dia) + x2_d * (len(list_dia) - 1)
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
                        cc += 1
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(smile_minimum_size).zfill(2) + "mem"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for var in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][var].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, var)))
            # plot
            fig_basic(plot_data[dur][var], list_names, 1, figure_axes[dur][var], fig_format,
                      figure_name_and_path, figure_panel_size[dur][var], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e15_change_with_iqr_plot()
