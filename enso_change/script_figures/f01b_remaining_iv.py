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
from typing import Literal

# numpy
from numpy.random import default_rng as numpy__random__default_rng

# scipy
from scipy.stats import scoreatpercentile as scipy__stats__scoreatpercentile

# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_tools import remaining_internal_variability_hs09, resample
from enso_change.wrapper.xarray_base import numpy__array, numpy__ndarray, concatenate
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
    "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585", "ssp370"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
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
    "proc_ps": 8000,
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
        "panel_01": {"x_delt": 200, "x_size": 2320, "y_delt": 100, "y_size": 400},
        "panel_11": {"x_delt": 200, "x_size": 2320, "y_delt": 200, "y_size": 400},
        "panel_21": {"x_delt": 200, "x_size": 2320, "y_delt": 250, "y_size": 400},
        "panel_31": {"x_delt": 320, "x_size": 1000, "y_delt": 200, "y_size": 800},
        "panel_41": {"x_delt": 200, "x_size": 1000, "y_delt": 200, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "1/x": "k",
        "val": "limegreen",
        "cor": "r",
        "member": "silver",
        "mean": "limegreen",
        "polynomial": "blue",
        "residual": "limegreen",
        "cmip5": "green",
        "cmip6": "r",
        "obse": "k",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "1/x": ":",
        "val": "-",
        "cor": "--",
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
            "MPI-ESM": {
                "std_ts_ano_n34e": {
                    "panel_01": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_11": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_21": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_31": {
                        "lim": [0, 800],
                        "nam": "variance (%$^2$)",
                        "tic": list(range(0, 801, 200)),
                    },
                    "panel_41": {
                        "lim": [0, 1],
                        "nam": "RIV$_{estimate}$ / IV$_{ref}$",
                        "tic": [round(k / 100, 2) for k in range(0, 101, 25)],
                    },
                },
                "sample": {
                    "panel_31": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
                    },
                    "panel_41": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
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
            "MPI-ESM1-2-LR": {
                "std_ts_ano_n34e": {
                    "panel_01": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_11": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_21": {
                        "lim": [-60, 60],
                        "nam": "",
                        "tic": list(range(-60, 61, 60)),
                    },
                    "panel_31": {
                        "lim": [0, 320],
                        "nam": "variance (%$^2$)",
                        "tic": list(range(0, 321, 80)),
                    },
                    "panel_41": {
                        "lim": [0, 1],
                        "nam": "RIV$_{estimate}$ / IV$_{ref}$",
                        "tic": [round(k / 100, 2) for k in range(0, 101, 25)],
                    },
                },
                "sample": {
                    "panel_31": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
                    },
                    "panel_41": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
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
                    "panel_31": {
                        "lim": [0, 280],
                        "nam": "variance (%$^2$)",
                        "tic": list(range(0, 281, 70)),
                    },
                    "panel_41": {
                        "lim": [0, 1],
                        "nam": "RIV$_{estimate}$ / IV$_{ref}$",
                        "tic": [round(k / 100, 2) for k in range(0, 101, 25)],
                    },
                },
                "sample": {
                    "panel_31": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
                    },
                    "panel_41": {
                        "lim": None,
                        "nam": "pseudo sample size",
                        "nbr_minor": None,
                        "tic": None,
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
        # curves on plot
        "member": "member(s)",
        "mean": "multi-member mean",
        "polynomial": "polynomial fit",
        "residual": "residual used to estimate RIV",
        "riv": "RIV estimates in MIROC6",
        "riv_iv": "RIV estimates across 15 SMILEs",
        "1/x": r"$\frac{1}{x}$",
        "cor": "corrected: y'=" + r"$\frac{1}{y0}x^{-0.2}$" + "y",
        "val": "estimate",
        # variables
        "std_pr_ano_n34e": "relative change of VALUE (%)",
        # "std_ts_ano_n34e": "relative change of VALUE (%)",
        "std_ts_ano_n34e": "$\\Delta$E relative to 1956-2005 (%)",
        # figure title
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
def f01b_remaining_iv_plot(
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
                    if exp1 in ["piControl"]:
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
                        put_in_dict(change, da1.compute(), dur1, dia, dat, exp1)
    print("change", list(change.keys()))
    change_resample, change_resample_max, change_sample = {}, {}, {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in sorted(list(values_json[dia][dur1][pro].keys()), key=lambda s: s.lower()):
                    if exp1 in ["piControl"]:
                        continue
                    for dat in sorted(list(values_json[dia][dur1][pro][exp1].keys()), key=lambda s: s.lower()):
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
                            list_pss = list(range(1, 11)) + [es]
                            if smile_minimum_size == 20:
                                list_pss += list(range(12, 21, 2))
                            elif smile_minimum_size > 20:
                                list_pss += list(range(25, 51, 5))
                            list_pss = [k for k in sorted(list(set(list_pss))) if k <= es]
                        for kk in list_pss:
                            # pseudo samples
                            da_o1, da_r1 = resample(da1, arr_i2=da2, data_var_i1=dia, data_var_i2=dia, dim="member",
                                                    pseudo_sample_size=proc_ps, sample_method=proc_md, sample_size=kk)
                            # samples
                            da_o2 = da1.sel(indexers={"member": members[:kk]})
                            da_r2 = da2.sel(indexers={"member": members[:kk]})
                            # ensemble mean for the reference (initial) value
                            da_r1, da_r2 = da_r1.mean(dim="member"), da_r2.mean(dim="member")
                            # departure from the reference value
                            da_o1, da_o2 = da_o1 - da_r1, da_o2 - da_r2
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o1, da_o2 = da_o1 * 100 / da_r1, da_o2 * 100 / da_r2
                            da_o1, da_o2 = da_o1.compute(), da_o2.compute()
                            # save
                            put_in_dict(change_resample, da_o1, dur1, dia, dat, exp1, kk)
                            if kk == max(list_pss):
                                put_in_dict(change_resample_max, da_o1, dur1, dia, dat, exp1)
                            if es >= 50 and exp1 in ["rcp85", "ssp585"]:
                                put_in_dict(change_sample, da_o2, dur1, dia, dat, exp1, kk)
    print("change_resample", list(change_resample.keys()))
    # min-max
    print("change sample min-max")
    for dur in list(change_sample.keys()):
        for dia in list(change_sample[dur].keys()):
            for dat in list(change_sample[dur][dia].keys()):
                tmp = []
                for siz in list(change_sample[dur][dia][dat].keys()):
                    for da in change_sample[dur][dia][dat][siz].values():
                        tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, dia, dat.ljust(17), "{0:.1f}".format(round(min(tmp), 1)).rjust(8),
                          "{0:.1f}".format(round(max(tmp), 1)).rjust(8))
    # variance (reference internal variability)
    sdev_change1, _, _ = compute_1var(change, "statistic",
                                      function_arguments={"dim": "member", "stat": "var", "ddof": 1})
    sdev_change2, _, _ = compute_1var(change_resample_max, "statistic",
                                      function_arguments={"dim": "member", "stat": "var", "ddof": 1})
    # -- compute residual internal variability
    example_polynomial_fit, example_residual, example_riv = {}, {}, {}
    for dur in list(change_sample.keys()):
        if dur == dur2:
            continue
        for dia in list(change_sample[dur].keys()):
            for dat in list(change_sample[dur][dia].keys()):
                for exp in list(change_sample[dur][dia][dat].keys()):
                    if exp == "piControl":
                        continue
                    for siz in list(change_sample[dur][dia][dat][exp].keys()):
                        # array
                        da = copy__deepcopy(change_sample[dur][dia][dat][exp][siz])
                        # degree
                        degree = proc_th[dia.split("_")[0]]
                        # compute remaining internal variability
                        em, re, ri = remaining_internal_variability_hs09(da, degree=degree, dim="year", ddof=1)
                        # save (there should be only 1 experiment)
                        put_in_dict(example_polynomial_fit, em.compute(), dur, dia, dat, exp, siz)
                        put_in_dict(example_residual, re.compute(), dur, dia, dat, exp, siz)
                        put_in_dict(example_riv, ri.compute(), dur, dia, dat, exp, siz)
    print("example_polynomial_fit", list(example_polynomial_fit.keys()))
    resample_riv = {}
    for dur1 in list(change_resample.keys()):
        if dur1 == dur2:
            continue
        for dia in list(change_resample[dur1].keys()):
            for dat in list(change_resample[dur1][dia].keys()):
                for exp in list(change_resample[dur1][dia][dat].keys()):
                    if exp == "piControl":
                        continue
                    for siz in list(change_resample[dur1][dia][dat][exp].keys()):
                        # array
                        da = copy__deepcopy(change_resample[dur1][dia][dat][exp][siz])
                        # degree
                        degree = proc_th[dia.split("_")[0]]
                        # compute remaining internal variability
                        _, _, ri = remaining_internal_variability_hs09(da, degree=degree, dim="year", ddof=1)
                        # save
                        put_in_dict(resample_riv, ri.compute(), dur1, dia, dat, exp, siz)
    print("resample_riv", list(resample_riv.keys()))
    # min-max
    print("resample_riv min-max")
    for dur in list(resample_riv.keys()):
        for dia in list(resample_riv[dur].keys()):
            for dat in list(resample_riv[dur][dia].keys()):
                if not check_dict_levels(change_sample, dur, dia, dat):
                    continue
                tmp = []
                for exp in list(resample_riv[dur][dia][dat].keys()):
                    if not check_dict_levels(change_sample, dur, dia, dat, exp):
                        continue
                    for da in resample_riv[dur][dia][dat][exp].values():
                        tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, dia, dat.ljust(17), "{0:.1f}".format(round(min(tmp), 1)).rjust(8),
                          "{0:.1f}".format(round(max(tmp), 1)).rjust(8))
    # -- compute riv / iv ratio
    resample_riv_iv = {}
    for dur in list(resample_riv.keys()):
        for dia in list(resample_riv[dur].keys()):
            # -- multi-model multi-experiment riv / iv
            da_o = ()
            for dat in list(resample_riv[dur][dia].keys()):
                da1, da2 = (), ()
                for exp in list(resample_riv[dur][dia][dat].keys()):
                    if exp in ["piControl"]:
                        continue
                    # -- remaining internal variability
                    # concatenate array (ensemble sizes)
                    di = copy__deepcopy(resample_riv[dur][dia][dat][exp])
                    list_k = sorted(di.keys())
                    da = concatenate(tuple(di[kk] for kk in list_k), "ensemble_size")
                    da = da.assign_coords({"ensemble_size": numpy__array(list_k)})
                    da1 += (da,)
                    # da1 += (da.mean(dim="resample"),)
                    # -- internal variability
                    da = copy__deepcopy(sdev_change1[dur][dia][dat][exp]).mean(dim="year")
                    # da = copy__deepcopy(sdev_change2[dur][dia][dat][exp]).mean(dim="year")
                    # if "resample" in list(da.dims):
                    #     da = da.mean(dim="resample")
                    da2 += (da,)
                # concatenate (experiments)
                da1, da2 = concatenate(da1, "experiment"), concatenate(da2, "experiment")
                # means
                da1 = da1.mean(dim="experiment", skipna=True).mean(dim="resample")
                da2 = da2.mean(dim="experiment", skipna=True)
                if "resample" in list(da2.dims):
                    da2 = da2.mean(dim="resample")
                # ratio
                da_o += (da1 / da2,)
            # concatenate (dataset)
            da_o = concatenate(da_o, "dataset").mean(dim="dataset", skipna=True)
            # select pseudo ensemble size
            da_o = da_o.sel(indexers={
                "ensemble_size": [k for k in da_o["ensemble_size"].values if k <= smile_minimum_size]})
            # save
            tk = (dur, dia, "val")
            print(dur, dia, "{0:.3f}".format(round(float(da_o.sel(indexers={"ensemble_size": 1}).values), 3)).rjust(8))
            put_in_dict(resample_riv_iv, copy__deepcopy(da_o).compute(), *tk)
            # -- riv / iv correction
            size = da_o["ensemble_size"].values
            da_o = da_o * (1 / da_o[0]) * size ** -0.2
            # save
            tk = (dur, dia, "cor")
            put_in_dict(resample_riv_iv, copy__deepcopy(da_o).compute(), *tk)
    print("resample_riv_iv", list(resample_riv_iv.keys()))
    # min-max
    print("resample_riv_iv min-max")
    for dur in list(resample_riv_iv.keys()):
        for dia in list(resample_riv_iv[dur].keys()):
            tmp = []
            for da in resample_riv_iv[dur][dia].values():
                tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.3f}".format(round(min(tmp), 3)).rjust(8),
                      "{0:.3f}".format(round(max(tmp), 3)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group = "group_01"
    for dur in list(example_riv.keys()):
        for dia in list(example_riv[dur].keys()):
            for dat in list(example_riv[dur][dia].keys()):
                for exp in list(example_riv[dur][dia][dat].keys()):
                    list_size = sorted(list(example_riv[dur][dia][dat][exp].keys()))
                    list_size = list_size[:2] + list_size[-1:]
                    for ii, siz in enumerate(list_size):
                        list_typ = ["change"]
                        for jj, typ in enumerate(list_typ):
                            panel = "panel_" + str(ii * 10 + 1).zfill(2)
                            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                            tuple_k = ("_".join([dur, dia, dat, exp]), grp, pan)
                            print(dur, dia, dat.ljust(15), exp.ljust(6), str(siz).rjust(3), typ, tuple_k)
                            # panel size
                            tk = tuple_k[:-2]
                            if not check_dict_levels(figure_ps, *tk + (pan,)):
                                put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                                for k1, k2 in fig_panel_size.items():
                                    if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                        put_in_dict(figure_ps, k2, *tk + (k1,))
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
                            # plot tests
                            pt1 = typ == "change" and check_dict_levels(change_sample, dur, dia, dat, exp, siz) \
                                  and check_dict_levels(example_polynomial_fit, dur, dia, dat, exp, siz)
                            pt2 = typ == "internal variability" and \
                                  check_dict_levels(example_residual, dur, dia, dat, exp, siz)
                            list_leg = []
                            # -- curve
                            lc, lo, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                            if pt1:
                                # individual members
                                tmp = "member"
                                for da in change_sample[dur][dia][dat][exp][siz]:
                                    lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                    lo += [True]
                                    ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                    lw += [0.5]
                                    lx += [da["year"].values]
                                    ly += [da.values]
                                    lz += [2]
                                list_leg += [tmp]
                                # member mean
                                tmp = "mean"
                                da = copy__deepcopy(change_sample[dur][dia][dat][exp][siz]).mean(dim="member")
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                lo += [True]
                                ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [4]
                                list_leg += [tmp]
                                # fit
                                tmp = "polynomial"
                                da = copy__deepcopy(example_polynomial_fit[dur][dia][dat][exp][siz])
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                lo += [True]
                                ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                lw += [2]
                                lx += [da["year"].values]
                                ly += [da.values]
                                lz += [5]
                                list_leg += [tmp]
                            elif pt2:
                                # residual
                                tmp = "residual"
                                da = copy__deepcopy(example_residual[dur][dia][dat][exp][siz])
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                lo += [True]
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
                                        lo += [True]
                                        ls += [":"]
                                        lw += [1]
                                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                        lz += [10]
                                # legend
                                if ii == 0:
                                    for kk, tmp in enumerate(list_leg):
                                        lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                        lo += [True]
                                        ls += [fig_linestyles[tmp] if tmp in list(fig_linestyles.keys()) else ":"]
                                        lw += [0.5 if tmp == "member" else 2]
                                        lx += [[x1 + (((kk % 2) * 70) + ((kk // 2) * 170) + ll) * dx
                                                for ll in [10, 18]]]
                                        ly += [[y2 - 5 * dy] * 2]
                                        lz += [1]
                                # dash between panels
                                if siz == list_size[-2]:
                                    yd, ys = df[pan]["y_delt"], df[pan]["y_size"]
                                    lc += ["k"]
                                    lo += [False]
                                    ls += [":"]
                                    lw += [5]
                                    lx += [[(x1 + x2) / 2] * 2]
                                    ly += [[y1 - 5 * dy, y1 + 5 * dy - (yd / ys) * (y2 - y1)]]
                                    lz += [1]
                                    if jj == 0:
                                        lc += ["k"]
                                        lo += [False]
                                        ls += [":"]
                                        lw += [5]
                                        lx += [[x1 - 20 * dx] * 2]
                                        ly += [[y1 + 2 * dy, y1 - 2 * dy - (yd / ys) * (y2 - y1)]]
                                        lz += [1]
                            plot_type = "cur"
                            if len(lx) > 0 and len(lx) == len(ly):
                                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                                put_in_dict(plot_data, copy__deepcopy(lo), *tuple_k + (str(plot_type) + "_on",))
                                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                            # -- shading
                            la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
                            if pt1:
                                tmp = "residual"
                                # member mean
                                da1 = copy__deepcopy(change_sample[dur][dia][dat][exp][siz]).mean(dim="member")
                                # fit
                                da2 = copy__deepcopy(example_polynomial_fit[dur][dia][dat][exp][siz])
                                # fill between
                                la += [0.6]
                                lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                lx += [da1["year"].values]
                                ly1 += [da1.values]
                                ly2 += [da2.values]
                                lz += [3]
                                list_leg += [tmp]
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_ps, *tk)
                                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                                if ii == 0:
                                    for kk, tmp in enumerate(list_leg[-1:]):
                                        la += [0.6]
                                        lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                        lx += [[x2 - ll * dx for ll in [2, 10]]]
                                        ly1 += [[y1 + (15 - 3) * dy] * 2]
                                        ly2 += [[y1 + (15 + 3) * dy] * 2]
                                        lz += [1]
                                # reference epoch
                                la += [0.4]
                                lc += ["grey"]
                                lx += [[1956, 2005]]
                                ly1 += [[y1] * 2]
                                ly2 += [[y1 + 10 * dy] * 2]
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
                            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_ps, *tk)
                                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                                # riv
                                if pt2 and check_dict_levels(example_riv, dur, dia, dat, exp, siz):
                                    da = float(copy__deepcopy(example_riv[dur][dia][dat][exp][siz]).values)
                                    lc += ["k"] * 2
                                    lf += [15] * 2
                                    lh += ["right"] * 2
                                    lr += [0] * 2
                                    lt += ["RIV = ", str(int(round(da, 0))) + "%$^{2}$"]
                                    lv += ["center"] * 2
                                    lw += ["regular"] * 2
                                    lx += [x2 - 25 * dx, x2 - 2 * dx]
                                    ly += [y1 + 5 * dy] * 2
                                # reference epoch
                                lc += ["k"]
                                lf += [15]
                                lh += ["center"]
                                lr += [0]
                                lt += ["reference epoch"]
                                lv += ["center"]
                                lw += ["regular"]
                                lx += [(1956 + 2005) / 2]
                                ly += [y1 + 4 * dy]
                                # legend
                                if ii == 0:
                                    for kk, tmp in enumerate(list_leg):
                                        lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                                        lf += [15]
                                        lh += ["left" if tmp != list_leg[-1] else "right"]
                                        lr += [0]
                                        lt += [fig_titles[tmp]]
                                        lv += ["center"]
                                        lw += ["regular"]
                                        lx += [x1 + ((kk * 70 if kk == 1 else 0) + ((kk // 2) * 170) + 20) * dx
                                               if tmp != list_leg[-1] else x2 - 12 * dx]
                                        ly += [y2 - 5 * dy if tmp != list_leg[-1] else y1 + 15 * dy]
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
                                    list_panels = [k for k in list(fig_panel_size.keys())
                                                   if k[:6] == "panel_" and k[-2:] in ["01", "11", "21"]]
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
                                    lf += [19]
                                    lh += ["center"]
                                    lr += [0]
                                    lt += [fig_titles["figure_title"]]
                                    lv += ["center"]
                                    lw += ["bold"]
                                    a2 = x1_m * len(list_typ) + x1_d * (len(list_typ) - 1)
                                    lx += [x1 + (a2 / (2 * x1_m)) * (x2 - x1)]
                                    ly += [y2 + 10 * dy]
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
    # panel 31: riv distributions
    print("riv distributions")
    group, panel = "group_01", "panel_31"
    for dur in list(resample_riv.keys()):
        for dia in list(resample_riv[dur].keys()):
            for dat in list(resample_riv[dur][dia].keys()):
                for exp in list(resample_riv[dur][dia][dat].keys()):
                    if not check_dict_levels(change_sample, dur, dia, dat, exp):
                        continue
                    pan = copy__deepcopy(panel)
                    grp = list(read_from_dict(plot_data, "_".join([dur, dia, dat, exp])).keys())[-1]
                    grp = str(grp.split("_")[0]) + "_" + str(int(grp.split("_")[-1]) + 1).zfill(2)
                    tuple_k = ("_".join([dur, dia, dat, exp]), grp, pan)
                    print(dur, dia, dat.ljust(15), exp.ljust(6), tuple_k)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_ps, *tk + (pan,)):
                        put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                put_in_dict(figure_ps, k2, *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], ["sample", dia]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, dat, dd, panel):
                            for k1, k2 in fig_ticks[dur][dat][dd][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k3 == "y_nam":
                                    tmp = tmp.replace("UNITS", metadata_json[dia]["units"])
                                elif k3 == "x_lim" and tmp is None:
                                    tmp = [smile_minimum_size + 1, 0]
                                elif k3 == "x_nbr_minor" and tmp is None:
                                    tmp = 1 if smile_minimum_size == 10 else (4 if smile_minimum_size == 20 else 1)
                                elif k3 == "x_tic" and tmp is None:
                                    aa = 2 if smile_minimum_size == 10 else (5 if smile_minimum_size == 20 else 10)
                                    tmp = list(range(smile_minimum_size, -1, -aa))
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
                    put_in_dict(plot_data, fig_titles["riv"], *tuple_k + ("title_col",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    m1, m2 = None, None
                    if isinstance(do_legend, bool) and do_legend and \
                            check_dict_levels(sdev_change1, dur, dia, dat, exp):
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        # mean internal variability
                        da = copy__deepcopy(sdev_change1[dur][dia][dat][exp]).mean(dim="year")
                        if "resample" in list(da.dims):
                            da = da.mean(dim="resample")
                        lc += ["purple"]
                        ls += ["-"]
                        lw += [4]
                        lx += [[x1, x2]]
                        ly += [[float(da.values)] * 2]
                        lz += [11]
                        m1 = float(da.values)
                    lc += [fig_colors["residual"] if "residual" in list(fig_colors.keys()) else "cyan"]
                    ls += ["-"]
                    lw += [4]
                    lx += [[siz for siz, da in resample_riv[dur][dia][dat][exp].items()]]
                    ly += [[float(da.mean()) for siz, da in resample_riv[dur][dia][dat][exp].items()]]
                    lz += [11]
                    m2 = float(resample_riv[dur][dia][dat][exp][1].mean())
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        kk, typ = 0, "mean"
                        lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                        ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else ":"]
                        lw += [4]
                        lx += [[x1 + 2 * dx, x1 + 14 * dx]]
                        ly += [[y2 - (10 * kk + 35) * dy] * 2]
                        lz += [1]
                    plot_type = "cur"
                    if len(lx) > 0 and len(lx) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                    # -- violin
                    la, lbod_fc, lc, lf, lls, llw, lmea, lmea_lw = [], [], [], [], [], [], [], []
                    lmed, lmed_lw, lo, lq, lq_ls, lq_lw, lw, lx, ly, lz = [], [], [], [], [], [], [], [], [], []
                    for siz, da in resample_riv[dur][dia][dat][exp].items():
                        arr = [float(k) for k in da.values]
                        tmin = scipy__stats__scoreatpercentile(arr, 1)
                        tmax = scipy__stats__scoreatpercentile(arr, 99)
                        arr = [k for k in arr if tmin <= k <= tmax]
                        la += [0.35]
                        lbod_fc += ["k"]
                        lc += ["k"]
                        lf += [True]
                        lls += ["-"]
                        llw += [1]
                        lmea += [False]
                        lmea_lw += [0]
                        lmed += [True]
                        lmed_lw += [5]
                        lo += ["vertical"]
                        lq += [[0.25, 0.75]]
                        lq_ls += ["-"]
                        lq_lw += [1]
                        lw += [0.98]
                        lx += [siz]
                        ly += [arr]
                        lz += [2]
                    normal = None
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        normal = numpy__random__default_rng().normal(loc=x1 + 60 * dx, scale=abs(15 * dx), size=10000)
                        tmin = scipy__stats__scoreatpercentile(normal, 1)
                        tmax = scipy__stats__scoreatpercentile(normal, 99)
                        normal = numpy__array([k for k in normal if tmin <= k <= tmax])
                        la += [0.35]
                        lbod_fc += ["grey"]
                        lc += ["grey"]
                        lf += [True]
                        lls += ["-"]
                        llw += [1]
                        lmea += [False]
                        lmea_lw += [0]
                        lmed += [True]
                        lmed_lw += [4]
                        lo += ["horizontal"]
                        lq += [[0.25, 0.75]]
                        lq_ls += ["-"]
                        lq_lw += [1]
                        lw += [abs(8 * dy)]
                        lx += [y2 - 0.05 * (y2 - y1)]
                        ly += [normal]
                        lz += [2]
                    plot_type = "vio"
                    if len(lx) > 0 and len(lx) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
                        put_in_dict(plot_data, copy__deepcopy(lbod_fc), *tuple_k + (str(plot_type) + "_bod_fc",))
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_flie",))
                        put_in_dict(plot_data, copy__deepcopy(lls), *tuple_k + (str(plot_type) + "_ls",))
                        put_in_dict(plot_data, copy__deepcopy(llw), *tuple_k + (str(plot_type) + "_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lmea), *tuple_k + (str(plot_type) + "_mea",))
                        put_in_dict(plot_data, copy__deepcopy(lmea_lw), *tuple_k + (str(plot_type) + "_mea_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lmed), *tuple_k + (str(plot_type) + "_med",))
                        put_in_dict(plot_data, copy__deepcopy(lmed_lw), *tuple_k + (str(plot_type) + "_med_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lo), *tuple_k + (str(plot_type) + "_o",))
                        put_in_dict(plot_data, copy__deepcopy(lq), *tuple_k + (str(plot_type) + "_q",))
                        put_in_dict(plot_data, copy__deepcopy(lq_ls), *tuple_k + (str(plot_type) + "_qua_ls",))
                        put_in_dict(plot_data, copy__deepcopy(lq_lw), *tuple_k + (str(plot_type) + "_qua_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                    # -- marker
                    lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                    if isinstance(m1, float) and isinstance(m2, float):
                        lc += ["purple"] + [fig_colors["residual"] if "residual" in list(fig_colors.keys()) else "cyan"]
                        lm += ["D"] * 2
                        ls += [200] * 2
                        lw += [0] * 2
                        lx += [1] * 2
                        ly += [m1, m2]
                        lz += [1] * 2
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
                    lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend curve
                        kk, typ = 0, "mean"
                        lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                        lf += [15]
                        lh += ["left"]
                        lr += [0]
                        lt += ["mean estimated RIV"]
                        lv += ["center"]
                        lw += ["regular"]
                        lx += [x1 + 16 * dx]
                        ly += [y2 - (10 * kk + 35) * dy]
                        if isinstance(m1, float) and isinstance(m2, float):
                            lc += ["purple"] + \
                                  [fig_colors["residual"] if "residual" in list(fig_colors.keys()) else "cyan"]
                            lf += [14, 15]
                            lh += ["left", "right"]
                            lr += [0] * 2
                            lt += ["reference IV (from all 50 members)", "estimated IV"]
                            lv += ["bottom", "center"]
                            lw += ["regular"] * 2
                            lx += [x1 + 1 * dx, 3 - 2 * dx]
                            ly += [m1 + 2 * dy, m2]
                        # legend distribution
                        if isinstance(normal, numpy__ndarray):
                            lc += ["grey"] * 6
                            lf += [15] * 6
                            lh += ["center"] + ["left"] + ["center"] + ["right"] + ["center"] * 2
                            lr += [0] * 6
                            lt += [str(k) + ("" if 1 < k < 99 else "%") for k in [99, 75, 50, 25, 1]] + \
                                  ["distribution of estimated RIV"]
                            lv += ["top"] * 6
                            lw += ["regular"] * 6
                            lx += [scipy__stats__scoreatpercentile(normal, k) if 1 < k < 99 else
                                   (float(normal.min()) if k == 1 else float(normal.max()))
                                   for k in [1, 25, 50, 75, 99]] + [scipy__stats__scoreatpercentile(normal, 50)]
                            ly += [y2 - 0.05 * (y2 - y1) - abs(4.5 * dy)] * 5 + [y2 - 0.05 * (y2 - y1) - abs(12 * dy)]
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
    # panel 41: multi-model multi-ensemble riv / iv
    print("multi-model multi-ensemble riv / iv")
    group, panel = "group_01", "panel_41"
    for dur in list(resample_riv.keys()):
        for dia in list(resample_riv[dur].keys()):
            for dat in list(resample_riv[dur][dia].keys()):
                for exp in list(resample_riv[dur][dia][dat].keys()):
                    if not check_dict_levels(change_sample, dur, dia, dat, exp):
                        continue
                    pan = copy__deepcopy(panel)
                    grp = list(read_from_dict(plot_data, "_".join([dur, dia, dat, exp])).keys())[-1]
                    grp = str(grp.split("_")[0]) + "_" + str(int(grp.split("_")[-1])).zfill(2)
                    tuple_k = ("_".join([dur, dia, dat, exp]), grp, pan)
                    print(dur, dia, dat.ljust(15), exp.ljust(6), tuple_k)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_ps, *tk + (pan,)):
                        put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                put_in_dict(figure_ps, k2, *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], ["sample", dia]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, dat, dd, panel):
                            for k1, k2 in fig_ticks[dur][dat][dd][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k3 == "y_nam":
                                    tmp = tmp.replace("UNITS", metadata_json[dia]["units"])
                                elif k3 == "x_lim" and tmp is None:
                                    tmp = [smile_minimum_size + 1, 0]
                                elif k3 == "x_nbr_minor" and tmp is None:
                                    tmp = 1 if smile_minimum_size == 10 else (4 if smile_minimum_size == 20 else 1)
                                elif k3 == "x_tic" and tmp is None:
                                    aa = 2 if smile_minimum_size == 10 else (5 if smile_minimum_size == 20 else 10)
                                    tmp = list(range(smile_minimum_size, -1, -aa))
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
                    put_in_dict(plot_data, fig_titles["riv_iv"], *tuple_k + ("title_col",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    list_leg = []
                    if check_dict_levels(resample_riv_iv, dur, dia):
                        for typ, da in resample_riv_iv[dur][dia].items():
                            list_leg += [typ]
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else ":"]
                            lw += [3]
                            lx += [da["ensemble_size"].values]
                            ly += [da.values]
                            lz += [3]
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # 1/x
                        typ = "1/x"
                        list_leg += [typ]
                        val = list(range(max(1, min(x1, x2)), max(x1, x2) + 1, 1))
                        lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                        ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else "cyan"]
                        lw += [3]
                        lx += [val]
                        ly += [[1 / k for k in val]]
                        lz += [2]
                        # legend
                        for kk, typ in enumerate(list_leg):
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[typ] if typ in list(fig_linestyles.keys()) else "cyan"]
                            lw += [3]
                            lx += [[x1 + 2 * dx, x1 + 14 * dx]]
                            ly += [[y2 - (10 * kk + 15) * dy] * 2]
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
                    lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                    if check_dict_levels(resample_riv_iv, dur, dia):
                        for typ, da in resample_riv_iv[dur][dia].items():
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            lm += ["D"]
                            ls += [100]
                            lw += [0.5]
                            lx += [1]
                            ly += [float(da.sel(indexers={"ensemble_size": 1}).values)]
                            lz += [1]
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
                    lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        for kk, typ in enumerate(list_leg):
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["left"]
                            lt += [fig_titles[typ] if typ in list(fig_titles.keys()) else typ]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x1 + 17 * dx]
                            ly += [y2 - (kk * 10 + 15) * dy]
                        if check_dict_levels(resample_riv_iv, dur, dia):
                            for typ, da in resample_riv_iv[dur][dia].items():
                                val = float(da.sel(indexers={"ensemble_size": 1}).values)
                                lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["right"]
                                lt += ["{0:.2f}".format(round(val, 2))]
                                lv += ["top" if val >= y2 - 2 * dy else ("bottom" if val <= y1 + 2 * dy else "center")]
                                lw += ["regular"]
                                lx += [1 - 7 * dx]
                                ly += [val + (-2 if val >= y2 - 2 * dy else (2 if val <= y1 + 2 * dy else 0)) * dy]
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
        ["", str(len(data_experiments)).zfill(2) + "exp", str(smile_minimum_size).zfill(2) + "exp"])
    # figure loop
    for k in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[k].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", k)))
        # plot
        fig_basic(plot_data[k], list_names, 1, figure_axes[k], fig_format, figure_name_and_path,
                  figure_ps[k], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#
