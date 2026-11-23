# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Ratio of remaining internal variability estimated / internal variability a function of ensemble size
# Test the influence of the pseudo sample size on the remaining internal variability
#    * pseudo sample size (x) vs. ratio rIV_ps / IV_ref (y)
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
from enso_change.wrapper.xarray_base import numpy__array, align, concatenate, create_array, polyval
from enso_change.wrapper.xarray_tools import resample, statistic
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["std_ts_ano_n34e", "std_pr_ano_n34e"],
    "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585", "ssp370"],
    # "data_experiments": ["piControl", "historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6", "cmip5"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 20,
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
        "panel_1": {"x_delt": 200, "x_size": 800, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "1/x": "k",
        "hs09-1--cor": "limegreen",
        "hs09-1--raw": "limegreen",
        "hs09-2--cor": "r",
        "hs09-2--raw": "limegreen",
        "hs09-3--cor": "darkorchid",
        "hs09-3--raw": "darkorchid",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "1/x": ":",
        "hs09-1--cor": "--",
        "hs09-1--raw": "--",
        "hs09-2--cor": "-",
        "hs09-2--raw": "-.",
        "hs09-3--cor": "-.",
        "hs09-3--raw": "-.",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "std_pr_ano_n34e--mea": {
                "panel_1": {
                    "lim": [0, 1],
                    "nam": "RIV / IV",
                    "tic": [round(k / 100, 2) for k in range(0, 101, 25)],
                },
            },
            "std_ts_ano_n34e--mea": {
                "panel_1": {
                    "lim": [0, 1],
                    "nam": "RIV / IV",
                    "tic": [round(k / 100, 2) for k in range(0, 101, 25)],
                },
            },
            "sample": {
                "panel_1": {
                    # "lim": [20, 0],
                    # "nam": "pseudo sample size",
                    # "nbr_minor": 4,
                    # "tic": list(range(20, -1, -5)),
                    "lim": None,
                    "nam": "pseudo sample size",
                    "nbr_minor": None,
                    "tic": None,
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "hs09-1": "1$^{st}$ order",
        "hs09-2": "2$^{nd}$ order",
        "hs09-3": "3$^{rd}$ order",
        "cor": "correction",
        "raw": "estimation",
        # "figure_title": "Remaining internal variability estimation",
        "figure_title": "",
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
def e56b_residual_iv_plot(
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
                    if exp1 == "piControl":
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
                        # save
                        put_in_dict(change, da1.compute(), dur1, dia, dat, exp1)
    print("change", list(change.keys()))
    change_resample = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if exp1 == "piControl":
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
                        # ensemble size
                        es = da1["member"].size
                        # compute pseudo samples
                        if exp1 == "piControl":
                            list_pss = [es]
                        else:
                            list_pss = list(range(1, 11)) + [es]
                            if smile_minimum_size == 20:
                                list_pss += list(range(12, 21, 2))
                            elif smile_minimum_size > 20:
                                list_pss += list(range(25, 51, 5))
                            list_pss = [k for k in sorted(list(set(list_pss))) if k <= es]
                        for kk in list_pss:
                            # pseudo samples
                            da_o, da_r = resample(da1, arr_i2=da2, data_var_i1=dia, data_var_i2=dia, dim="member",
                                                  pseudo_sample_size=proc_ps, sample_method=proc_md, sample_size=kk)
                            # ensemble mean for the reference (initial) value
                            da_r = da_r.mean(dim="member")
                            # departure from the reference value
                            da_o -= da_r
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o *= 100 / da_r
                            # save
                            put_in_dict(change_resample, da_o.compute(), dur1, dia, dat, kk, exp1)
    print("change_resample", list(change_resample.keys()))
    # variance (reference internal variability)
    sdev_change, _, _ = compute_1var(change, "statistic",
                                     function_arguments={"dim": "member", "stat": "var", "ddof": 1})
    # -- compute residual internal variability
    internal = {}
    for dur in list(change_resample.keys()):
        if dur == dur2:
            continue
        for dia in list(change_resample[dur].keys()):
            for dat in list(change_resample[dur][dia].keys()):
                for met in ["hs09"]:
                    # l1 = list(range(1, 4)) if met == "hs09" else [20]
                    l1 = [2] if met == "hs09" else [20]
                    for k1 in l1:
                        for k2 in list(change_resample[dur][dia][dat].keys()):
                            da_o = ()
                            for exp in list(change_resample[dur][dia][dat][k2].keys()):
                                if exp == "piControl":
                                    continue
                                # array
                                da = copy__deepcopy(change_resample[dur][dia][dat][k2][exp])
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
                                # store residual internal variability
                                da_o += (es,)
                            if len(da_o) == 0:
                                continue
                            # multi-experiment mean
                            da_o = concatenate(da_o, "experiment").mean(dim="experiment")
                            # save
                            met_o = copy__deepcopy(met) if met == "theo" else str(met) + "-" + str(k1)
                            put_in_dict(internal, da_o.compute(), dur, dia, dat, met_o, "all", k2)
    print("internal", list(internal.keys()))
    internal_error = {}
    for dur in list(internal.keys()):
        for dia in list(internal[dur].keys()):
            for dat in list(internal[dur][dia].keys()):
                for met in list(internal[dur][dia][dat].keys()):
                    list_exp1 = list(internal[dur][dia][dat][met].keys())
                    list_exp2 = copy__deepcopy(list_exp1)
                    if met == "ld20" and "piControl" not in list_exp1 and "piControl" in data_experiments:
                        list_exp2.append("piControl")
                    for exp in list_exp2:
                        if exp == "piControl" and exp not in list_exp1:
                            da1 = copy__deepcopy(sdev_change[dur][dia][dat][exp])
                            if "year" in list(da1.dims):
                                da1 = da1.mean(dim="year")
                            da1 = concatenate((da1,) * 2, "ensemble_size")
                            list_k = [0, 100]
                        else:
                            d1 = copy__deepcopy(internal[dur][dia][dat][met][exp])
                            list_k = sorted(d1.keys())
                            # concatenate array (ensemble sizes)
                            if exp == "piControl":
                                da1 = concatenate((d1[list_k[0]],) * 2, "ensemble_size")
                                list_k = [0, 100]
                            else:
                                da1 = concatenate(tuple(d1[kk] for kk in list_k), "ensemble_size")
                        da1 = da1.assign_coords({"ensemble_size": numpy__array(list_k)})
                        # reference
                        d2 = copy__deepcopy(sdev_change[dur][dia][dat])
                        da2 = tuple(da.mean(dim="year") for k, da in d2.items() if k != "piControl")
                        da2 = concatenate(da2, "experiment")
                        da2 = float(da2.mean(dim="experiment").values)
                        # relative error
                        da1 = da1 / da2
                        # compute stats
                        # list_stat = ["mea"] if exp == "piControl" else ["5p", "25p", "mea", "75p", "95p"]
                        list_stat = ["mea"]
                        for sta in list_stat:
                            sta2 = "per" if sta[-1] == "p" else copy__deepcopy(sta)
                            per2 = float(sta[:-1]) if sta[-1] == "p" else None
                            da_o = copy__deepcopy(da1)
                            if "resample" in list(da_o.dims):
                                da_o = statistic(da1, dim="resample", stat=sta2, per=per2)
                            # save
                            put_in_dict(internal_error, da_o, dur, dia, sta, met, exp, dat)
    internal_plot = {}
    for dur in list(internal_error.keys()):
        for dia in list(internal_error[dur].keys()):
            for sta in list(internal_error[dur][dia].keys()):
                for met in list(internal_error[dur][dia][sta].keys()):
                    for exp in list(internal_error[dur][dia][sta][met].keys()):
                        d1 = copy__deepcopy(internal_error[dur][dia][sta][met][exp])
                        # -- compute multi-model mean
                        # find largest ensemble
                        da_o = {da["ensemble_size"].size: dat for dat, da in d1.items()}
                        dat = da_o[sorted(list(da_o.keys()))[-1]]
                        # create a dummy array with all pseudo ensemble sizes
                        da_t = create_array(d1[dat], data_var="mmm")
                        # concatenate all ensembles and average
                        da_o = tuple(align(da_t, da, join="outer")[1] for dat, da in d1.items())
                        da_o = concatenate(da_o, "dataset").mean(dim="dataset", skipna=True)
                        if exp == "all":
                            da_o = da_o.sel(indexers={
                                "ensemble_size": [k for k in da_o["ensemble_size"].values if k <= smile_minimum_size]})
                        # save
                        tk = (dur, dia, sta, "raw", met, exp, "mmm")
                        print(*tk, da_o[0])
                        put_in_dict(internal_plot, copy__deepcopy(da_o), *tk)
                        # -- correction
                        size = da_o["ensemble_size"].values
                        da_o = da_o * (1 / da_o[0]) * size ** -0.2
                        # da_o = da_o * (100 / da_o[0]) * size**0.02
                        # save
                        tk = (dur, dia, sta, "cor", met, exp, "mmm")
                        put_in_dict(internal_plot, copy__deepcopy(da_o), *tk)
    # min-max
    print("change min-max")
    for dur in list(internal_plot.keys()):
        for dia in list(internal_plot[dur].keys()):
            for sta in list(internal_plot[dur][dia].keys()):
                t1 = []
                for typ in list(internal_plot[dur][dia][sta].keys()):
                    t2 = []
                    for met in list(internal_plot[dur][dia][sta][typ].keys()):
                        for exp in list(internal_plot[dur][dia][sta][typ][met].keys()):
                            for dat, da in internal_plot[dur][dia][sta][typ][met][exp].items():
                                t2 += da.values.flatten().tolist()
                                print(dur, dia, sta.rjust(6), typ.rjust(25), met, exp, dat,
                                      "{0:.0f}".format(round(min(da.values.flatten().tolist()), 0)).rjust(4),
                                      "{0:.0f}".format(round(max(da.values.flatten().tolist()), 0)).rjust(4))
                    t1 += t2
                if len(t1) > 0:
                    print(dur, dia, sta.rjust(6), "{0:.1f}".format(round(min(t1), 1)).rjust(8),
                          "{0:.1f}".format(round(max(t1), 1)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in data_epoch_lengths:
        if not check_dict_levels(internal_plot, dur):
            continue
        sta = "mea"
        ny = int(dur.split("_")[0])
        if dur == dur2:
            continue
        list_dia = [k for k in data_diagnostics1 if k in list(internal_plot[dur].keys())]
        for jj, dia in enumerate(list_dia):
            if not check_dict_levels(internal_plot, dur, dia, sta):
                continue
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
            grp = copy__deepcopy(group)
            tuple_k = (dur, grp, pan)
            print(dur, dia, tuple_k)
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                        put_in_dict(figure_panel_size, k2, *tk + (k1,))
            # panel axes definition
            l1, l2 = ["x_", "y_"], ["sample", str(dia) + "--" + str(sta)]
            for vv, dd in zip(l1, l2):
                if check_dict_levels(fig_ticks, dur, dd, panel):
                    for k1, k2 in fig_ticks[dur][dd][panel].items():
                        k3 = str(vv) + str(k1)
                        if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if k3 == "x_lim" and tmp is None:
                            tmp = [smile_minimum_size, 0]
                        elif k3 == "x_nbr_minor" and tmp is None:
                            tmp = 1 if smile_minimum_size == 10 else (4 if smile_minimum_size == 20 else 1)
                        elif k3 == "x_tic" and tmp is None:
                            aa = 2 if smile_minimum_size == 10 else (5 if smile_minimum_size == 20 else 10)
                            tmp = list(range(smile_minimum_size, -1, -aa))
                        if k3 == "y_nam" and jj != 0:
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
            put_in_dict(plot_data, metadata_json[dia]["name_short"], *tuple_k + ("title_col",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            for typ in list(internal_plot[dur][dia][sta].keys()):
                for met in list(internal_plot[dur][dia][sta][typ].keys()):
                    for exp in list(internal_plot[dur][dia][sta][typ][met].keys()):
                        for da in internal_plot[dur][dia][sta][typ][met][exp].values():
                            tt = str(met) + "--" + str(typ)
                            lc += [fig_colors[exp] if exp == "piControl" else fig_colors[tt]]
                            ls += [fig_linestyles[exp] if exp == "piControl" else fig_linestyles[tt]]
                            lw += [2]
                            lx += [da["ensemble_size"].values]
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
                # 1/x
                val = list(range(max(1, min(x1, x2)), max(x1, x2) + 1, 1))
                lc += ["k"]
                ls += [":"]
                lw += [2]
                lx += [val]
                ly += [[1 / k for k in val]]
                lz += [4]
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
                list_m = []
                for k1 in list(internal_plot[dur][dia][sta].keys()):
                    for k2 in list(internal_plot[dur][dia][sta][k1].keys()):
                        list_m.append(str(k2) + "--" + str(k1))
                list_m += ["1/x"]
                for cc, met in enumerate(list_m):
                    k1, k2 = None, None
                    if "--" in met:
                        k2, k1 = met.split("--")
                    lc += [fig_colors[met] if met in list(fig_colors.keys()) else "k"]
                    lf += [15]
                    lh += ["left"]
                    if k1 is not None and len(list(internal_plot[dur][dia][sta][k1].keys())) == 1 and k1 == "raw":
                        # val = str(fig_titles[k1]) + " using " + str(fig_titles[k2])
                        val = str(fig_titles[k1])
                    elif k1 is not None and len(list(internal_plot[dur][dia][sta][k1].keys())) == 1 and k1 == "cor":
                        val = str(fig_titles[k1]) + ": y'=" + r"$\frac{1}{y0}x^{-0.2}$" + "y"
                    else:
                        val = copy__deepcopy(met)
                    lt += [val]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x1 + 2 * dx]
                    ly += [y2 - (cc * 10 + 15) * dy]
                # plot title
                if jj == 0:
                    # plot size
                    x1_d = fig_panel_size[panel]["x_delt"]
                    x1_m = fig_panel_size[panel]["x_size"]
                    # add to text
                    lc += ["k"]
                    lf += [20]
                    lh += ["center"]
                    lt += [fig_titles["figure_title"]]
                    lv += ["center"]
                    lw += ["bold"]
                    a2 = x1_m * len(list_dia) + x1_d * (len(list_dia) - 1)
                    lx += [x1 + (a2 / (2 * x1_m)) * (x2 - x1)]
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
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(len(data_experiments)).zfill(2) + "exp", str(smile_minimum_size).zfill(2) + "exp", proc_md[:2]])
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
e56b_residual_iv_plot()
