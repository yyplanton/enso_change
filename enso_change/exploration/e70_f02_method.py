# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure 2: method (ensemble sizes, normal distribution, STD error fct of ensemble size in SMILEs)
#    * ensemble size and cumulative % of members from largest
#    * snapshot of model distribution at several epochs (or flattened)
#    * STD error fct of ensemble size in SMILEs
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
from numpy.random import normal as numpy__random__normal
# scipy
from scipy.stats import normaltest as scipy__stats__normaltest
from scipy.stats import shapiro as scipy__stats__shapiro

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.numpy_tools import numpy__zeros
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import array_wrapper, numpy__array, align, concatenate, create_array, polyval
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
    "proc_ci": 90.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {
        "std_pr_ano_n34e": {"hs09--2": [1.07, 1, 5], "ld20": [1, 2, 5]},
        "std_ts_ano_n34e": {"hs09--2": [1.10, 1, 5], "ld20": [1, 2, 5]}
    },
    # method: {"combination", "random"}
    "proc_md": "random",
    # number of pseudo samples: int
    "proc_ps": 6000,
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
        "panel_01": {"x_delt": 300, "x_size": 750, "y_delt": 200, "y_size": 800},
        "panel_11": {"x_delt": 300, "x_size": 750, "y_delt": 100, "y_size": 800},
        "panel_16": {"x_delt": 300, "x_size": 750, "y_delt": 200, "y_size": 800},
        "panel_21": {"x_delt": 440, "x_size": 1205, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "historical": "b",
        "ssp126": "b",
        "ssp245": "r",
        "ssp370": "b",
        "ssp585": "r",
        "hs09-1": "limegreen",
        "hs09-2": "y",
        "hs09-3": "darkorchid",
        "ld20": "r",
        "05p": "k",
        "25p": "k",
        "mea": "k",
        "75p": "k",
        "95p": "k",
    },
    # line styles
    "fig_linestyles": {
        "historical": "-",
        "ssp126": "-",
        "ssp245": "--",
        "ssp370": "-",
        "ssp585": "--",
        "05p": "--",
        "25p": "-.",
        "mea": "-",
        "75p": "-.",
        "95p": "--",
    },
    # markers
    "fig_markers": {
        "historical": "v",
        "ssp126": "v",
        "ssp245": "o",
        "ssp370": "v",
        "ssp585": "o",
    },
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "std_pr_ano_n34e": {
                "panel_11": {
                    "_lim": [-4, 4],
                    "_nam": "normalized values",
                    "_tic": list(range(-4, 5, 2)),
                },
                "panel_16": {
                    "_lim": [-4, 4],
                    "_nam": "normalized values",
                    "_tic": list(range(-4, 5, 2)),
                },
                "panel_21": {
                    "_lim": [-50, 50],
                    "_nam": "relative error (%)",
                    "_tic": list(range(-50, 51, 25)),
                },
            },
            "std_ts_ano_n34e": {
                "panel_11": {
                    "_lim": [-4, 4],
                    "_nam": "normalized values",
                    "_tic": list(range(-4, 5, 2)),
                },
                "panel_16": {
                    "_lim": [-4, 4],
                    "_nam": "normalized values",
                    "_tic": list(range(-4, 5, 2)),
                },
                "panel_21": {
                    "_lim": [-50, 50],
                    "_nam": "relative error (%)",
                    "_tic": list(range(-50, 51, 25)),
                },
            },
            "percentage": {
                "panel_01": {
                    "_axe": "left",
                    "_col": "r",
                    "_lim": [0, 40],
                    "_nam": "members per model (%)",
                    "_tic": list(range(0, 41, 10)),
                    "2_axe": "right",
                    "2_lim": [0, 100],
                    "2_nam": "cumulative members (%)",
                    "2_tic": list(range(0, 101, 25)),
                },
            },
            "probability": {
                "panel_11": {
                    "_lim": [0, 0.4],
                    "_nam": "probability",
                    "_tic": [round(k / 10, 1) for k in range(0, 5, 1)],
                },
                "panel_16": {
                    "_lim": [0, 0.4],
                    "_nam": "probability",
                    "_tic": [round(k / 10, 1) for k in range(0, 5, 1)],
                },
            },
            "sample": {
                "panel_01": {
                    "_lim": [55, 0],
                    "_nam": "ensemble size",
                    "_tic": list(range(60, -1, -10)),
                },
                "panel_21": {
                    "_lim": [0, 55],
                    "_nam": "pseudo sample size",
                    "_tic": list(range(0, 61, 10)),
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
        "distribution": "distribution of\nVALUE",
        "ensemble_size": "ensemble size",
        "error": "internal variability of\nVALUE",
        "hs09-1": "Hawkins & Sutton (2009)",
        "hs09-2": "Hawkins & Sutton (2009)",
        "hs09-3": "Hawkins & Sutton (2009)",
        "ld20": "Lehner et al. (2020)",
        "05p": "5-95%",
        "25p": "25-75%",
        "75p": "25-75%",
        "95p": "5-95%",
        "5.0p": "5-95%",
        "25.0p": "25-75%",
        "75.0p": "25-75%",
        "95.0p": "5-95%",
        "mea": "mean",
        "figure_title": "Multimodel spread with or without internal variability",
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
def e70_f02_method_plot(
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
                        if not check_dict_levels(values_json, dia, dur2, pro, exp1) or \
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
                        put_in_dict(change_em, da1, dur1, pro, dia, exp1, dat)
    # -- ensemble size and cumulative
    ensemble_size, member_cum = {}, {}
    for dur in list(change_em.keys()):
        if dur == dur2:
            continue
        for pro in list(change_em[dur].keys()):
            for dia in list(change_em[dur][pro].keys()):
                for exp in list(change_em[dur][pro][dia].keys()):
                    if pro == "observations" or dia != data_diagnostics1[0] or exp == "piControl":
                        continue
                    # models of a given size
                    sizes = numpy__array(list(range(60)))
                    da_e = numpy__zeros(sizes.shape)
                    for da in change_em[dur][pro][dia][exp].values():
                        da_e[int(da["member"].size)] += 1
                    # cumulative sum
                    da_m = da_e * sizes
                    da_m = da_m[::-1].cumsum()[::-1]
                    # to xarray
                    da_e = array_wrapper(da_e * 100 / da_e.sum(), coords={"size": sizes}, dims=["size"], name="models")
                    da_m = array_wrapper(da_m * 100 / da_m[0], coords={"size": sizes}, dims=["size"], name="members")
                    # save
                    put_in_dict(ensemble_size, da_e, dur, pro, dia, exp)
                    put_in_dict(member_cum, da_m, dur, pro, dia, exp)
    # min-max
    print("ensemble_size min-max")
    for dur in list(ensemble_size.keys()):
        for pro in list(ensemble_size[dur].keys()):
            for dia in list(ensemble_size[dur][pro].keys()):
                tmp = []
                for da in ensemble_size[dur][pro][dia].values():
                    tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, pro, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    # -- compute normalized distributions
    # ensemble mean
    mean_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # mmm & mms
    change_m, change_s, change_t = {}, {}, {}
    for dur in list(mean_change.keys()):
        if dur == dur2:
            continue
        for pro in list(mean_change[dur].keys()):
            for dia in list(mean_change[dur][pro].keys()):
                for exp in list(mean_change[dur][pro][dia].keys()):
                    if pro == "observations" or exp == "piControl":
                        continue
                    # arrays
                    da = concatenate(tuple(da for da in mean_change[dur][pro][dia][exp].values()), "dataset")
                    # multimodel mean
                    da_m = da.mean(dim="dataset")
                    # multimodel spread
                    da_s = da.std(dim="dataset")
                    # normalize
                    da = (da - da_m) / da_s
                    # select given epoch
                    yy = da["year"].values
                    yy = [k for k in yy if k <= 1960] if exp == "historical" else [k for k in yy if k >= 2020]
                    da = da.sel(indexers={"year": yy[len(yy) // 2]})
                    # save
                    put_in_dict(change_t, da, dur, pro, dia, exp)
    # min-max
    print("change min-max")
    for dur in list(change_t.keys()):
        for pro in list(change_t[dur].keys()):
            for dia in list(change_t[dur][pro].keys()):
                tmp = []
                for da in change_t[dur][pro][dia].values():
                    tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, pro, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    # -- compute internal variability
    # std (reference internal variability)
    for dur in list(change_em.keys()):
        for pro in list(change_em[dur].keys()):
            for dia in list(change_em[dur][pro].keys()):
                for exp in list(change_em[dur][pro][dia].keys()):
                    for dat in list(change_em[dur][pro][dia][exp].keys()):
                        da = copy__deepcopy(change_em[dur][pro][dia][exp][dat])
                        del change_em[dur][pro][dia][exp][dat]
                        if da["member"].size >= 20:
                            put_in_dict(change_em, da, dur, pro, dia, exp, dat)
                    if len(list(change_em[dur][pro][dia][exp].keys())) == 0:
                        del change_em[dur][pro][dia][exp]
                if len(list(change_em[dur][pro][dia].keys())) == 0:
                    del change_em[dur][pro][dia]
            if len(list(change_em[dur][pro].keys())) == 0:
                del change_em[dur][pro]
        if len(list(change_em[dur].keys())) == 0:
            del change_em[dur]
    sdev_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "std", "ddof": 1})
    # compute change compared to initial conditions using pseudo samples
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                # dictionary
                d1 = values_json[dia][dur][pro]
                # list datasets
                list_dat = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                for dat in list_dat:
                    # list experiments
                    list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                    # check if ensembles are large enough
                    for exp in list_exp:
                        da = values_json[dia][dur][pro][exp][dat]
                        if da["member"].size < 20:
                            del values_json[dia][dur][pro][exp][dat]
                    # list experiments
                    d2 = copy__deepcopy(values_json[dia][dur][pro])
                    list_exp = [k1 for k1 in list(d2.keys()) if dat in list(d2[k1].keys())]
                    # delete datasets
                    if "historical" not in list_exp:
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                for exp in list(values_json[dia][dur][pro].keys()):
                    if len(list(values_json[dia][dur][pro][exp].keys())) == 0:
                        del values_json[dia][dur][pro][exp]
                if len(list(values_json[dia][dur][pro].keys())) == 0:
                    del values_json[dia][dur][pro]
            if len(list(values_json[dia][dur].keys())) == 0:
                del values_json[dia][dur]
        if len(list(values_json[dia].keys())) == 0:
            del values_json[dia]
    change_resample = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if pro == "observations":
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
                        da2 = da2.sel(indexers={"year": 2014})
                        # ensemble size
                        es1, es2 = da1["member"].size, da2["member"].size
                        if es1 < 20 or es2 < 20:
                            continue
                        # compute pseudo samples
                        if exp1 == "piControl":
                            list_pss = [es1]
                        else:
                            # list_pss = list(range(1, 11)) + list(range(12, 21, 2)) + list(range(25, es1 + 1, 5))
                            list_pss = list(range(1, 6)) + list(range(6, 21, 2)) + list(range(25, es1 + 1, 5))
                        if es1 not in list_pss:
                            list_pss.append(es1)
                        list_pss = sorted(list(set(list_pss)))
                        for kk in list_pss:
                            # pseudo samples
                            if exp1 == "piControl":
                                da_o, da_r = copy__deepcopy(da1), copy__deepcopy(da2)
                            else:
                                da_o, da_r = resample(da1, arr_i2=da2, data_var_i1=dia, data_var_i2=dia, dim="member",
                                                      pseudo_sample_size=proc_ps, sample_method=proc_md, sample_size=kk)
                            # ensemble mean for the reference (initial) value
                            da_r = da_r.mean(dim="member")
                            # departure from the reference value
                            da_o -= da_r
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o *= 100 / da_r
                            # save
                            put_in_dict(change_resample, da_o, dur1, pro, dia, dat, kk, exp1)
    # internal variability in pseudo samples
    internal = {}
    for dur in list(change_resample.keys()):
        for pro in list(change_resample[dur].keys()):
            for dia in list(change_resample[dur][pro].keys()):
                for dat in list(change_resample[dur][pro][dia].keys()):
                    typ = "corrected"  # "raw"  #
                    deg = proc_th[dia.split("_")[0]]
                    for met in ["hs09", "ld20"]:
                        for k2 in list(change_resample[dur][pro][dia][dat].keys()):
                            da_o = ()
                            for exp in list(change_resample[dur][pro][dia][dat][k2].keys()):
                                if (exp == "piControl" and met == "hs09") or (met == "ld20" and k2 < 2):
                                    continue
                                da = copy__deepcopy(change_resample[dur][pro][dia][dat][k2][exp])
                                size = da["member"].size if "member" in list(da.dims) else 1
                                if met == "hs09":
                                    # -- hs09
                                    # 'forced change' is a 2nd order polynomial fit: compute coefficient then fit
                                    em = copy__deepcopy(da)
                                    if "member" in list(em.dims):
                                        em = em.mean(dim="member")
                                    coef = em.polyfit("year", deg)
                                    em = polyval(da["year"], coef["polyfit_coefficients"], degree_dim="degree")
                                    # 'internal variability' is the temporal std of the residual
                                    es = (da - em).std(dim="year", ddof=1)
                                    # average across members
                                    if "member" in list(es.dims):
                                        es = es.mean(dim="member")
                                    # correction ?
                                    if typ == "corrected":
                                        if set(data_experiments) == {"piControl", "historical"}:
                                            es = es * (1.15 + 0.15 * deg / size)
                                        else:
                                            if "_".join(dia.split("_")[:3]) == "std_pr_ano":
                                                es = es * (1.07 + 6 * deg / (60 * size))
                                            else:
                                                es = es * (1.1 + 1 * deg / (10 * size))
                                else:
                                    # -- ld20
                                    # 'ensemble spread' is the multi-member std
                                    es = da.std(dim="member", ddof=1)
                                    # correction ?
                                    if typ == "corrected":
                                        es = es * (1. + 2 / (5 * size))
                                if "year" in list(es.dims):
                                    es = es.mean(dim="year")
                                # save
                                if exp == "piControl":
                                    # save
                                    put_in_dict(internal, es, dur, pro, dia, dat, met, exp, k2)
                                else:
                                    da_o += (es,)
                            if len(da_o) == 0:
                                continue
                            # multi-experiment mean
                            da_o = concatenate(da_o, "experiment").mean(dim="experiment")
                            # save
                            met_o = str(met) + "-" + str(deg) if met == "hs09" else copy__deepcopy(met)
                            put_in_dict(internal, da_o, dur, pro, dia, dat, met_o, "all", k2)
    # error relative to reference
    internal_error = {}
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    for dur in list(internal.keys()):
        for pro in list(internal[dur].keys()):
            for dia in list(internal[dur][pro].keys()):
                for dat in list(internal[dur][pro][dia].keys()):
                    for met in list(internal[dur][pro][dia][dat].keys()):
                        list_exp1 = list(internal[dur][pro][dia][dat][met].keys())
                        list_exp2 = copy__deepcopy(list_exp1)
                        if met == "ld20" and "piControl" not in list_exp1 and "piControl" in data_experiments:
                            list_exp2.append("piControl")
                        for exp in list_exp2:
                            if exp == "piControl" and exp not in list_exp1:
                                da1 = copy__deepcopy(sdev_change[dur][pro][dia][exp][dat])
                                if "year" in list(da1.dims):
                                    da1 = da1.mean(dim="year")
                                da1 = concatenate((da1,) * 2, "size")
                                list_k = [0, 100]
                            else:
                                d1 = copy__deepcopy(internal[dur][pro][dia][dat][met][exp])
                                list_k = sorted(d1.keys())
                                # concatenate array (ensemble sizes)
                                if exp == "piControl":
                                    da1 = concatenate((d1[list_k[0]],) * 2, "size")
                                    list_k = [0, 100]
                                else:
                                    da1 = concatenate(tuple(d1[kk] for kk in list_k), "size")
                            da1 = da1.assign_coords({"size": numpy__array(list_k)})
                            # reference
                            d2 = copy__deepcopy(sdev_change[dur][pro][dia])
                            d2 = tuple(d2[ex2][dat].mean(dim="year") for ex2 in list(d2.keys())
                                       if ex2 != "piControl" and dat in list(d2[ex2].keys()))
                            da2 = concatenate(d2, "experiment")
                            da2 = float(da2.mean(dim="experiment").values)
                            # relative error
                            da1 = (da1 - da2) * 100 / da2
                            # compute stats
                            list_stat = ["mea"]
                            if exp != "piControl":
                                list_stat += [
                                    (str(int(k)).zfill(2) if k == int(k) else str(k)) + "p" for k in percentiles]
                            for sta in list_stat:
                                sta2 = "per" if sta[-1] == "p" else copy__deepcopy(sta)
                                per2 = float(sta[:-1]) if sta[-1] == "p" else None
                                da_o = copy__deepcopy(da1)
                                if "resample" in list(da_o.dims):
                                    da_o = statistic(da1, dim="resample", stat=sta2, per=per2)
                                # save
                                put_in_dict(internal_error, da_o, dur, pro, dia, sta, met, exp, dat)
    # multimodel mean internal variability
    internal_cur = {}
    for dur in list(internal_error.keys()):
        for pro in list(internal_error[dur].keys()):
            for dia in list(internal_error[dur][pro].keys()):
                for sta in list(internal_error[dur][pro][dia].keys()):
                    for met in list(internal_error[dur][pro][dia][sta].keys()):
                        for exp in list(internal_error[dur][pro][dia][sta][met].keys()):
                            if exp == "piControl":
                                continue
                            d1 = copy__deepcopy(internal_error[dur][pro][dia][sta][met][exp])
                            # -- compute multi-model mean
                            # find largest ensemble
                            da_o = {da["size"].size: dat for dat, da in d1.items()}
                            dat = da_o[sorted(list(da_o.keys()))[-1]]
                            # create a dummy array with all pseudo ensemble sizes
                            da_t = create_array(d1[dat], data_var="mmm")
                            # concatenate all ensembles and average
                            da_o = ()
                            for dat, da in d1.items():
                                da_o += (align(da_t, da, join="outer")[1],)
                            da_o = concatenate(da_o, "dataset").mean(dim="dataset", skipna=True)
                            # save
                            put_in_dict(internal_cur, da_o, dur, pro, dia, str(met) + "--" + str(sta))
    # multimodel mean internal variability
    internal_sha = {}
    for dur in list(internal_cur.keys()):
        for pro in list(internal_cur[dur].keys()):
            for dia in list(internal_cur[dur][pro].keys()):
                list_met = [k.split("--")[0] for k in list(internal_cur[dur][pro][dia].keys())]
                for met in list_met:
                    list_sta = sorted(list(set([k.split("--")[1] for k in list(internal_cur[dur][pro][dia].keys())
                                                if met in k and "mea" not in k])), key=lambda s: s.lower())
                    # range
                    da1 = copy__deepcopy(internal_cur[dur][pro][dia][str(met) + "--" + str(list_sta[0])])
                    da2 = copy__deepcopy(internal_cur[dur][pro][dia][str(met) + "--" + str(list_sta[-1])])
                    # select
                    if met == "ld20":
                        da1 = da1.sel(indexers={"size": [k for k in da1["size"] if k > 4]})
                        da2 = da2.sel(indexers={"size": [k for k in da2["size"] if k > 4]})
                    else:
                        da1 = da1.sel(indexers={"size": [k for k in da1["size"] if k <= 4]})
                        da2 = da2.sel(indexers={"size": [k for k in da2["size"] if k <= 4]})
                    # concatenate
                    da_o = concatenate((da1, da2), "quantile")
                    # save
                    put_in_dict(internal_sha, da_o, dur, pro, dia, met)
    # min-max
    print("error min-max")
    for dur in list(internal_cur.keys()):
        for pro in list(internal_cur[dur].keys()):
            for dia in list(internal_cur[dur][pro].keys()):
                tmp = []
                for da in internal_cur[dur][pro][dia].values():
                    tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, pro.ljust(12), dia, "{0:.1f}".format(round(min(tmp), 1)).rjust(8),
                          "{0:.1f}".format(round(max(tmp), 1)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: methods
    print("methods")
    group, panel = "group_1", "panel_1"
    for dur in data_epoch_lengths:
        for pro in data_projects:
            list_c = [member_cum, None, internal_cur]
            list_h = [None, change_t, None]
            list_m = [ensemble_size, None, None]
            list_s = [None, None, internal_sha]
            ln = ["ensemble_size", "distribution", "error"]
            for c1, (dc, dh, dm, ds, nn) in enumerate(zip(list_c, list_h, list_m, list_s, ln)):
                list_dict = [
                    k for k in [dc, dh, dm, ds] if isinstance(k, dict) and check_dict_levels(k, dur, pro)]
                if len(list_dict) == 0:
                    continue
                list_dia = []
                for dd in list_dict:
                    if len(list_dia) == 0:
                        list_dia = list(dd[dur][pro].keys())
                    else:
                        list_dia = list(set(list_dia) & set(list(dd[dur][pro].keys())))
                list_dia = [k for k in data_diagnostics1 if k in list_dia]
                list_exp = []
                for dia in list_dia:
                    for dd in list_dict:
                        if check_dict_levels(dd, dur, pro, dia):
                            list_exp += list(dd[dur][pro][dia].keys())
                list_exp = sorted(list(set(list_exp)), key=lambda s: s.lower())
                list_exp_o = [k for k in data_experiments if k in list_exp]
                for ii, dia in enumerate(list_dia):
                    dict_e = {exp: [exp] for exp in list_exp_o}
                    if c1 in [0, 1]:
                        dict_e = {"historical": ["historical"], "ssp126": ["ssp126", "ssp245"],
                                  "ssp370": ["ssp370", "ssp585"]}
                    elif c1 == 2:
                        dict_e = {"method": list_exp}
                    list_e = list(dict_e.keys())
                    for jj, exp in enumerate(list_e):
                        if c1 == 2:
                            aa = c1 * 10
                            tp = "panel_" + str(aa + 1).zfill(2)
                            pan = str(tp.split("_")[0]) + "_" + str(int(tp.split("_")[-1]) + ii).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + aa).zfill(2)
                        else:
                            aa = c1 * 10 + ii * 5
                            tp = "panel_" + str(aa + 1).zfill(2)
                            pan = str(tp.split("_")[0]) + "_" + str(int(tp.split("_")[-1]) + jj).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + aa).zfill(2)
                        tuple_k = (dur, pro, grp, pan)
                        print(dur, pro, nn.ljust(15), dia, exp.ljust(10), tuple_k)
                        # panel size
                        tk = tuple_k[:2]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[tp]), *tk + (pan,))
                            for k1, k2 in fig_panel_size.items():
                                if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                    put_in_dict(figure_panel_size, k2, *tk + (k1,))
                        # panel axes definition
                        l1, l2 = ["x", "y"], []
                        if c1 == 0:
                            l2 = ["sample", "percentage"]
                        elif c1 == 1:
                            l2 = [dia, "probability"]
                        elif c1 == 2:
                            l2 = ["sample", dia]
                        for vv, dd in zip(l1, l2):
                            if check_dict_levels(fig_ticks, dur, dd, tp):
                                for k1, k2 in fig_ticks[dur][dd][tp].items():
                                    k3 = str(vv) + str(k1)
                                    if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                        continue
                                    tmp = copy__deepcopy(k2)
                                    if isinstance(tmp, str) and k1 == "_nam" and dd in [dia]:
                                        tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                        tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                        tmp = tmp.replace(" ()", "")
                                    if c1 != 2 and ((k3 == "x_nam" and dia != list_dia[-1]) or
                                                    (k3 == "y_nam" and jj != 0) or
                                                    (k3 == "y2_nam" and exp != list_e[-1])):
                                        tmp = ""
                                    put_in_dict(figure_axes, copy__deepcopy(tmp), *tuple_k + (k3,))
                        dt = {}
                        if check_dict_levels(figure_axes, *tuple_k):
                            dt = read_from_dict(figure_axes, *tuple_k)
                        # test possibility to plot text / legend
                        do_legend = False
                        if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None \
                                and dt["y_lim"] is not None:
                            do_legend = True
                        # -- title
                        if c1 == 0 and ii == 0:
                            val = ""
                            for ex2 in dict_e[exp]:
                                if val != "":
                                    val += " / "
                                val += str(fig_titles[ex2]) + " (" + str(len(n_ensemble[ex2])) + ")"
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0:
                            val = copy__deepcopy(fig_titles[nn]).replace("VALUE", metadata_json[dia]["name_short"])
                            put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        if isinstance(dc, dict) and check_dict_levels(dc, dur, pro, dia):
                            for ex2 in dict_e[exp]:
                                if check_dict_levels(dc, dur, pro, dia, ex2):
                                    v1, v2 = ex2.split("--") if c1 == 2 else (ex2, ex2)
                                    da = dc[dur][pro][dia][ex2]
                                    lc += [fig_colors[v1] if c1 == 2 and v1 in list(fig_colors.keys()) else "k"]
                                    ls += [fig_linestyles[v2] if v2 in list(fig_linestyles.keys()) else ":"]
                                    lw += [2]
                                    lx += [da["size"].values]
                                    ly += [da.values]
                                    lz += [1]
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + (dt["y2_lim"] if c1 == 0 else dt["y_lim"])
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # x = 0 and y = 0
                            if c1 == 2:
                                for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                                    if k1[0] < 0 < k1[1]:
                                        lc += ["k"]
                                        ls += [":"]
                                        lw += [1]
                                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                        lz += [10]
                            # legend
                            if isinstance(dc, dict):
                                list_ls = [k for k in dict_e[exp] if check_dict_levels(dc, dur, pro, dia, k)]
                                if c1 == 2:
                                    list_ls = list(set([k.split("--")[-1] for k in list_ls]))
                                    list_ls = list(reversed(sorted(
                                        [k for k in list_ls if k == "mea" or float(k[:-1]) < 50],
                                        key=lambda s: s.lower())))
                                for k1, k2 in enumerate(list_ls):
                                    lc += ["k"]
                                    ls += [fig_linestyles[k2] if k2 in list(fig_linestyles.keys()) else ":"]
                                    lw += [2]
                                    if c1 == 0:
                                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                                        ly += [[y2 - (k1 * 8 + 15) * dy] * 2]
                                    else:
                                        lx += [[x2 - 10 * dx, x2 - 2 * dx]]
                                        ly += [[y2 - (k1 * 8 + 5) * dy] * 2]
                                    lz += [1]
                        plot_type = "cur2" if c1 == 0 else "cur"  # "cur"
                        if len(lx) > 0 and len(lx) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                            put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                            put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                        # -- histogram
                        lc, lbr, lbw, ld, le, lf, lk, ls = [], [], [], [], [], [], [], []
                        if isinstance(dh, dict) and check_dict_levels(dh, dur, pro, dia):
                            for ex2 in dict_e[exp]:
                                if check_dict_levels(dh, dur, pro, dia, ex2):
                                    da = dh[dur][pro][dia][ex2]
                                    lc += [fig_colors[ex2] if ex2 in list(fig_colors.keys()) else "cyan"]
                                    lbr += [[-5, 5]]
                                    lbw += [0.5]
                                    ld += [da.values.flatten()]
                                    le += ["step"]
                                    lf += [False]
                                    lk += [True]
                                    ls += ["probability"]
                        if isinstance(do_legend, bool) and do_legend and isinstance(dh, dict):
                            # normal distribution
                            lc += ["k"]
                            lbr += [[-5, 5]]
                            lbw += [0.5]
                            ld += [numpy__random__normal(size=10000)]
                            le += ["step"]
                            lf += [False]
                            lk += [True]
                            ls += ["probability"]
                            # legend
                        plot_type = "his"  # "his2" if c1 == 0  else "his"
                        if len(ld) > 0:
                            put_in_dict(plot_data, copy__deepcopy(lbr), *tuple_k + (str(plot_type) + "_binrange",))
                            put_in_dict(plot_data, copy__deepcopy(lbw), *tuple_k + (str(plot_type) + "_binwidth",))
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_color",))
                            put_in_dict(plot_data, copy__deepcopy(ld), *tuple_k + (str(plot_type) + "_data",))
                            put_in_dict(plot_data, copy__deepcopy(le), *tuple_k + (str(plot_type) + "_element",))
                            put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fill",))
                            put_in_dict(plot_data, copy__deepcopy(lk), *tuple_k + (str(plot_type) + "_kde",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_stat",))
                        # -- marker
                        lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                        if isinstance(dm, dict) and check_dict_levels(dm, dur, pro, dia):
                            for ex2 in dict_e[exp]:
                                if check_dict_levels(dm, dur, pro, dia, ex2):
                                    da = dm[dur][pro][dia][ex2]
                                    for k1, k2 in zip(da["size"].values, da.values):
                                        if k2 == 0:
                                            continue
                                        lc += ["r"]
                                        lm += [fig_markers[ex2] if ex2 in list(fig_markers.keys()) else "o"]
                                        ls += [100]
                                        lw += [2]
                                        lx += [k1]
                                        ly += [k2]
                                        lz += [2]
                        if isinstance(do_legend, bool) and do_legend and isinstance(dm, dict):
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            l1 = [k for k in dict_e[exp] if check_dict_levels(dm, dur, pro, dia, k)]
                            for k1, ex2 in enumerate(l1):
                                lc += ["r"]
                                lm += [fig_markers[ex2] if ex2 in list(fig_markers.keys()) else "o"]
                                ls += [100]
                                lw += [2]
                                lx += [x1 + 14 * dx]
                                ly += [y2 - (k1 * 8 + 15) * dy]
                                lz += [1]
                        plot_type = "mar"  # "mar2" if c1 == 0  else "mar"
                        if len(lx) > 0 and len(lx) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                            put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                            put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                        # -- shading
                        la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
                        if isinstance(ds, dict) and check_dict_levels(ds, dur, pro, dia):
                            for ex2 in dict_e[exp]:
                                if check_dict_levels(ds, dur, pro, dia, ex2):
                                    da = ds[dur][pro][dia][ex2]
                                    la += [0.3]
                                    lc += [fig_colors[ex2] if ex2 in list(fig_colors.keys()) else "cyan"]
                                    lx += [da["size"].values]
                                    ly1 += [da.min(dim="quantile").values]
                                    ly2 += [da.max(dim="quantile").values]
                                    lz += [2]
                        if isinstance(do_legend, bool) and do_legend and isinstance(ds, dict):
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            la += [0.3]
                            lc += ["k"]
                            lx += [[x2 - k * dx for k in [10, 2]]]
                            ly1 += [[y1 + (5 + 3) * dy] * 2]
                            ly2 += [[y1 + (5 - 3) * dy] * 2]
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
                            if isinstance(dc, dict) or isinstance(dh, dict):
                                list_ls = [ex2 for k1, ex2 in enumerate(dict_e[exp])
                                           if (c1 in [0, 2] and check_dict_levels(dc, dur, pro, dia, ex2)) or
                                           (c1 == 1 and check_dict_levels(dh, dur, pro, dia, ex2))]
                                if c1 == 2:
                                    list_ls1 = list(set([k.split("--")[-1] for k in list_ls]))
                                    list_ls2 = list(set([k.split("--")[0] for k in list_ls]))
                                    list_ls = list(reversed(sorted(
                                        [k for k in list_ls1 if k == "mea" or float(k[:-1]) < 50],
                                        key=lambda s: s.lower())))
                                    list_ls += sorted(list(set([k.split("--")[0] for k in list_ls2])),
                                                      key=lambda s: s.lower())
                                for k1, ex2 in enumerate(list_ls):
                                    if c1 == 0:
                                        lc += ["k"]
                                        lf += [15]
                                        lh += ["left"]
                                        lt += [fig_titles[ex2]]
                                        lv += ["center"]
                                        lw += ["regular"]
                                        lx += [x1 + 18 * dx]
                                        ly += [y2 - (k1 * 8 + 15) * dy]
                                    elif c1 == 1:
                                        da = dh[dur][pro][dia][ex2].values.flatten()
                                        lc += [fig_colors[ex2] if ex2 in list(fig_colors.keys()) else "cyan"]
                                        lf += [12.5]
                                        lh += ["left" if len(dict_e[exp]) < 2 else "right"]
                                        val = fig_titles[ex2]
                                        s1, p1 = scipy__stats__normaltest(da)
                                        s2, p2 = scipy__stats__shapiro(da)
                                        val += " (AP=" + "{0:.2f}".format(round(p1, 2)) + " ; SW=" + \
                                               "{0:.2f}".format(round(p2, 2)) + ")"
                                        lt += [val]
                                        lv += ["center"]
                                        lw += ["regular"]
                                        if len(dict_e[exp]) < 2:
                                            lx += [x1 + 2 * dx]
                                            ly += [y2 - (k1 * 8 + 15) * dy]
                                        else:
                                            lx += [x2 - 2 * dx]
                                            ly += [y2 - (k1 * 8 + 5) * dy]
                                    elif c1 == 2:
                                        lc += [fig_colors[ex2] if ex2 in list(fig_colors.keys()) else "cyan"]
                                        lf += [15]
                                        lh += ["right"]
                                        lt += [fig_titles[ex2]]
                                        lv += ["center"]
                                        lw += ["regular"]
                                        lx += [x2 - 12 * dx]
                                        ly += [y2 - (k1 * 8 + 5) * dy if k1 < 2 else
                                               y1 + ((len(list_ls) - k1) * 8 + 5) * dy]
                                if isinstance(ds, dict):
                                    lc += ["k"]
                                    lf += [15]
                                    lh += ["right"]
                                    lt += ["range of application of each method"]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x2 - 12 * dx]
                                    ly += [y1 + 5 * dy]
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
    a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a2) + "po"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for pro in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][pro].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, pro)))
            # plot
            fig_basic(plot_data[dur][pro], list_names, 1, figure_axes[dur][pro], fig_format,
                      figure_name_and_path, figure_panel_size[dur][pro], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e70_f02_method_plot()
