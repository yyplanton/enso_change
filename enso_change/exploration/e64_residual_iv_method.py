# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Illustration of the computation of RIV
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
from enso_change.wrapper.xarray_base import concatenate, polyval
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical", "ssp370", "ssp585"],
    # "data_experiments": ["piControl", "historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 50,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95,
    # method: {"combination", "random"}
    "proc_md": "random",
    # number of pseudo samples: int
    "proc_ps": 4000,
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
        "panel_01": {"x_delt": 200, "x_size": 6, "y_delt": 100, "y_size": 500},
        "panel_11": {"x_delt": 200, "x_size": 6, "y_delt": 300, "y_size": 500},
        "panel_21": {"x_delt": 200, "x_size": 6, "y_delt": 100, "y_size": 500},
    },
    # colors
    "fig_colors": {
        "member": "silver",
        "mean": "limegreen",
        "polynomial": "blue",
        "residual": "limegreen",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "member": "-",
        "mean": "-",
        "polynomial": "--",
        "residual": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "MPI-ESM1-2-LR": {
                "std_ts_ano_n34e": {
                    "panel_01": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
                    },
                    "panel_11": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
                    },
                    "panel_21": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
                    },
                },
                "time": {
                    "panel_01": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                    "panel_11": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                    "panel_21": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                },
            },
            "MIROC6": {
                "std_ts_ano_n34e": {
                    "panel_01": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 30)),
                    },
                    "panel_11": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 30)),
                    },
                    "panel_21": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 30)),
                    },
                },
                "time": {
                    "panel_01": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                    "panel_11": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                    "panel_21": {
                        "lim": [1860, 2090],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2101, 50)),
                    },
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "member": "member(s)",
        "mean": "multi-member mean",
        "polynomial": "polynomial fit",
        "residual": "residual (departure from polynomial fit)",
        "std_ts_ano_n34e": "relative change of VALUE (%)",
        "figure_title": "Estimating the Remaining Internal Variability (RIV)",
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
def e64_residual_iv_method_plot(
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
                    # select historical large enough
                    if "historical" in list_exp:
                        da_h = copy__deepcopy(values_json[dia][dur][pro]["historical"][dat])
                        del values_json[dia][dur][pro]["historical"][dat]
                        if da_h["member"].size < smile_minimum_size:
                            list_exp.remove("historical")
                        else:
                            put_in_dict(values_json, copy__deepcopy(da_h), dia, dur, pro, "historical", dat)
                    # delete all experiments if piControl or historical are not available
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
                        if len(lm) < smile_minimum_size:
                            continue
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
    n_ens = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ens, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        if exp not in list(n_ens.keys()):
            continue
        list_experiments = copy__deepcopy(n_ens[exp])
        del n_ens[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ens, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(10), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
    #
    # -- Processing
    #
    # compute change compared to initial conditions
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    change = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical"]:
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        d1 = values_json[dia][dur1][pro]
                        list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                        list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2) or (
                                exp1 == "historical" and len(list_ssp) > 0):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        # da2 = da2.sel(indexers={"year": 1850 + int(dur2.split("_")[0])}).mean(dim="member")
                        da2 = da2.sel(indexers={"year": 2006}).mean(dim="member")
                        # departure from the reference value
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # save
                        put_in_dict(change, da1, dur1, dia, dat, exp1)
    change_sample = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical"]:
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        d1 = values_json[dia][dur1][pro]
                        list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                        list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2) or (
                                exp1 == "historical" and len(list_ssp) > 0):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        # da2 = da2.sel(indexers={"year": 1850 + int(dur2.split("_")[0])})
                        da2 = da2.sel(indexers={"year": 2006})
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # ensemble size
                        members = da1["member"].values
                        es = da1["member"].size
                        # compute pseudo samples
                        if exp1 == "piControl":
                            list_pss = [es]
                        else:
                            list_pss = list(range(1, min(3, es + 1))) + [es]
                            list_pss = [k for k in sorted(list(set(list_pss))) if k <= es]
                        for kk in list_pss:
                            # samples
                            da_o = da1.sel(indexers={"member": members[:kk]})
                            da_r = da2.sel(indexers={"member": members[:kk]})
                            # ensemble mean for the reference (initial) value
                            da_r = da_r.mean(dim="member")
                            # departure from the reference value
                            da_o -= da_r
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o *= 100 / da_r
                            # save
                            put_in_dict(change_sample, da_o, dur1, dia, dat, kk, exp1)
    # -- compute residual internal variability
    polynomial_fit = {}
    internal = {}
    for dur in list(change_sample.keys()):
        if dur == dur2:
            continue
        for dia in list(change_sample[dur].keys()):
            for dat in list(change_sample[dur][dia].keys()):
                for met in ["hs09"]:
                    # l1 = list(range(1, 4)) if met == "hs09" else [20]
                    l1 = [2] if met == "hs09" else [20]
                    for k1 in l1:
                        for k2 in list(change_sample[dur][dia][dat].keys()):
                            da_o = ()
                            for exp in list(change_sample[dur][dia][dat][k2].keys()):
                                if exp == "piControl":
                                    continue
                                # array
                                da = copy__deepcopy(change_sample[dur][dia][dat][k2][exp])
                                size = da["member"].size
                                if met == "theo" and size < 3:
                                    continue
                                # compute residual internal variability (residual from polynomial)
                                # 'forced change' is a Nth order polynomial fit
                                if "member" in list(da.dims):
                                    da = da.mean(dim="member")
                                coef = da.polyfit("year", k1)
                                em = polyval(da["year"], coef["polyfit_coefficients"], degree_dim="degree")
                                # 'residual internal variability' is the temporal var of the residual (em - po)
                                es = (da - em).var(dim="year", ddof=1)
                                # save
                                put_in_dict(polynomial_fit, em, dur, dia, dat, exp, k2)
                                put_in_dict(internal, es, dur, dia, dat, exp, k2)
    # min-max
    print("change min-max")
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for dat in list(change[dur][dia].keys()):
                t1 = []
                for da in change[dur][dia][dat].values():
                    t1 += da.values.flatten().tolist()
                if len(t1) > 0:
                    print(dur, dia, dat.ljust(17), "{0:.1f}".format(round(min(t1), 1)).rjust(8),
                          "{0:.1f}".format(round(max(t1), 1)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group = "group_1"
    for dur in list(polynomial_fit.keys()):
        for dia in list(polynomial_fit[dur].keys()):
            for dat in list(polynomial_fit[dur][dia].keys()):
                for exp in list(polynomial_fit[dur][dia][dat].keys()):
                    list_size = list(polynomial_fit[dur][dia][dat][exp].keys())
                    for ii, siz in enumerate(list_size):
                        list_typ = ["change", "internal variability"]
                        for jj, typ in enumerate(list_typ):
                            panel = "panel_" + str(ii * 10 + 1).zfill(2)
                            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                            tuple_k = (dur, dia, dat, exp, grp, pan)
                            print(dur, dia, dat.ljust(15), exp, siz, typ, tuple_k)
                            # panel size
                            tk = tuple_k[:-2]
                            if not check_dict_levels(figure_ps, *tk + (pan,)):
                                for k1, d1 in fig_panel_size.items():
                                    if k1 == panel and not check_dict_levels(figure_ps, *tk + (pan,)):
                                        for k2, d2 in d1.items():
                                            if check_dict_levels(figure_ps, *tk + (pan, k2,)):
                                                continue
                                            tmp = copy__deepcopy(d2)
                                            tt = (dur, dat, "time", panel, "lim")
                                            if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
                                                k3 = read_from_dict(fig_ticks, *tt)
                                                tmp *= (k3[-1] - k3[0])
                                            put_in_dict(figure_ps, tmp, *tk + (pan, k2,))
                                    elif k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                        put_in_dict(figure_ps, copy__deepcopy(d1), *tk + (k1,))
                            # panel axes definition
                            l1, l2 = ["x_", "y_"], ["time", dia]
                            for vv, dd in zip(l1, l2):
                                if check_dict_levels(fig_ticks, dur, dat, dd, panel):
                                    for k1, k2 in fig_ticks[dur][dat][dd][panel].items():
                                        k3 = str(vv) + str(k1)
                                        if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                            continue
                                        tmp = copy__deepcopy(k2)
                                        if k3 == "x_nam" and siz != list_size[-1]:
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
                                put_in_dict(plot_data, typ, *tuple_k + ("title_col",))
                            # plot tests
                            pt1 = typ == "change" and check_dict_levels(change_sample, dur, dia, dat, siz, exp)
                            pt2 = typ == "internal variability" and \
                                  check_dict_levels(change_sample, dur, dia, dat, siz, exp)
                            list_leg = ["member", "mean", "polynomial"] if pt1 else (["residual"] if pt2 else [])
                            # -- curve
                            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                            if pt1:
                                # individual members
                                tmp = "member"
                                for da in change_sample[dur][dia][dat][siz][exp]:
                                    lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                    ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                    lw += [0.5]
                                    lx += [da["year"].values]
                                    ly += [da.values]
                                    lz += [2]
                                # member mean
                                tmp = "mean"
                                da = copy__deepcopy(change_sample[dur][dia][dat][siz][exp]).mean(dim="member")
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [3]
                                # fit
                                tmp = "polynomial"
                                da = copy__deepcopy(polynomial_fit[dur][dia][dat][exp][siz])
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [4]
                            elif pt2:
                                # residual
                                tmp = "residual"
                                da = copy__deepcopy(change_sample[dur][dia][dat][siz][exp]).mean(dim="member") - \
                                     copy__deepcopy(polynomial_fit[dur][dia][dat][exp][siz])
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [4]
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_ps, *tk)
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
                                for kk, tmp in enumerate(list_leg):
                                    lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                    ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                    lw += [0.5 if tmp == "member" else 2]
                                    lx += [[x1 + 10 * dx, x1 + 18 * dx] if kk == 0 else [x1 + 2 * dx, x1 + 10 * dx]]
                                    ly += [[y2 - (kk * 8 + 5) * dy] * 2]
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
                            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_ps, *tk)
                                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                                # riv
                                if pt2 and check_dict_levels(internal, dur, dia, dat, exp, siz):
                                    da = float(copy__deepcopy(internal[dur][dia][dat][exp][siz]).values)
                                    lc += ["k"] * 2
                                    lf += [15] * 2
                                    lh += ["right"] * 2
                                    lr += [0] * 2
                                    lt += ["RIV = ", str(int(round(da, 0))) + "%$^{2}$"]
                                    lv += ["center"] * 2
                                    lw += ["regular"] * 2
                                    lx += [x2 - 20 * dx, x2 - 2 * dx]
                                    ly += [y1 + 5 * dy] * 2
                                # legend
                                for kk, tmp in enumerate(list_leg):
                                    lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                    lf += [15]
                                    lh += ["left"]
                                    lr += [0]
                                    lt += [fig_titles[tmp]]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x1 + 20 * dx if kk == 0 else x1 + 12 * dx]
                                    ly += [y2 - (kk * 8 + 5) * dy]
                                # number of members
                                if jj == 0:
                                    lc += ["k"]
                                    lf += [15]
                                    lh += ["center"]
                                    lr += [90]
                                    lt += [str(siz) + " member" + ("" if siz < 2 else "s")]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x1 - 20 * dx]
                                    ly += [(y2 + y1) / 2]
                                # variable
                                if jj == 0 and siz == list_size[0]:
                                    # plot size
                                    y1_s = fig_panel_size[panel]["y_size"]
                                    y1_tot = 0
                                    list_panels = [k for k in list(fig_panel_size.keys()) if k[:6] == "panel_"]
                                    for kk in list_panels:
                                        y1_tot += fig_panel_size[kk]["y_size"] + \
                                                  (0 if kk == list_panels[-1] else fig_panel_size[kk]["y_delt"])
                                    # add to text
                                    lc += ["k"]
                                    lf += [15]
                                    lh += ["center"]
                                    lr += [90]
                                    val = copy__deepcopy(fig_titles[dia])
                                    val = val.replace("VALUE", metadata_json[dia]["name_short"])
                                    lt += [val.replace("UNITS", metadata_json[dia]["units"])]
                                    lv += ["center"]
                                    lw += ["bold"]
                                    lx += [x1 - 30 * dx]
                                    ly += [y2 - y1_tot * (y2 - y1) / (2 * y1_s)]
                                # plot title
                                if ii == 0 and jj == 0:
                                    # plot size
                                    x1_d = read_from_dict(figure_ps, *tk + (pan, "x_delt",))
                                    x1_m = read_from_dict(figure_ps, *tk + (pan, "x_size",))
                                    # add to text
                                    lc += ["k"]
                                    lf += [20]
                                    lh += ["center"]
                                    lr += [0]
                                    lt += [fig_titles["figure_title"]]
                                    lv += ["center"]
                                    lw += ["bold"]
                                    a2 = x1_m * len(list_typ) + x1_d * (len(list_typ) - 1)
                                    lx += [x1 + (a2 / (2 * x1_m)) * (x2 - x1)]
                                    ly += [y2 + 20 * dy]
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
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(len(data_experiments)).zfill(2) + "exp", str(smile_minimum_size).zfill(2) + "exp"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for dat in sorted(list(plot_data[dur][dia].keys()), key=lambda s: s.lower()):
                for exp in sorted(list(plot_data[dur][dia][dat].keys()), key=lambda s: s.lower()):
                    # list plot names
                    list_names = list(plot_data[dur][dia][dat][exp].keys())
                    # path to output figure
                    figure_name_and_path = os__path__join(figure_directory,
                                                          figure_name + "_".join(("", dur, dia, dat, exp)))
                    # plot
                    fig_basic(plot_data[dur][dia][dat][exp], list_names, 1,
                              figure_axes[dur][dia][dat][exp], fig_format, figure_name_and_path,
                              figure_ps[dur][dia][dat][exp], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e64_residual_iv_method_plot()
