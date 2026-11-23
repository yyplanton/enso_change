# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure 2: evolution of variable change in different experiments and regions
#    * time (x) vs. change averaged across models (y)
#        - departure from the end of the historical experiment (1956-2005)
#        - MMM for the signal
#        - Model range minus corrected remaining internal variability
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
from enso_change.wrapper.statistics import scipy__stats__t
from enso_change.wrapper.xarray_base import concatenate, numpy__array
from enso_change.wrapper.xarray_tools import remaining_internal_variability_hs09, statistic, statistic_resampled
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e",
    #                       "std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e"],
    # "data_diagnostics1": ["std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e"],
    # "data_diagnostics1": ["std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e"],
    # "data_diagnostics1": ["ske_ts_ano_n30e", "ske_ts_ano_n34e", "ske_ts_ano_n40e"],
    "data_diagnostics1": ["std_pr_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85",
    #                      "ssp126", "ssp245", "ssp585", "ssp370"],
    "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6", "cmip5"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    "smile_rejected": ["EC-EARTH", "MPI-ESM-LR"],  # default_arguments["smile_rejected"],
    # "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 60.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {
        # no resample 8000
        "ske_pr_ano_n30e": {"hs09--1": 0.77, "hs09--2": 0.63, "hs09--3": 54},
        "ske_pr_ano_n34e": {"hs09--1": 0.78, "hs09--2": 0.64, "hs09--3": 55},
        "ske_pr_ano_n40e": {"hs09--1": 0.76, "hs09--2": 0.63, "hs09--3": 53},
        "ske_ts_ano_n30e": {"hs09--1": 0.74, "hs09--2": 0.62, "hs09--3": 53},
        "ske_ts_ano_n34e": {"hs09--1": 0.74, "hs09--2": 0.62, "hs09--3": 53},
        "ske_ts_ano_n40e": {"hs09--1": 0.75, "hs09--2": 0.63, "hs09--3": 54},
        "std_pr_ano_n30e": {"hs09--1": 1.01, "hs09--2": 0.64, "hs09--3": 51},
        "std_pr_ano_n34e": {"hs09--1": 1.07, "hs09--2": 0.66, "hs09--3": 53},
        "std_pr_ano_n40e": {"hs09--1": 0.93, "hs09--2": 0.63, "hs09--3": 51},
        "std_ts_ano_n30e": {"hs09--1": 75, "hs09--2": 0.61, "hs09--3": 49},
        "std_ts_ano_n34e": {"hs09--1": 76, "hs09--2": 0.61, "hs09--3": 49},
        "std_ts_ano_n40e": {"hs09--1": 76, "hs09--2": 0.61, "hs09--3": 50},
        # resample 8000
        # "std_pr_ano_n34e": {"hs09--1": 1.11, "hs09--2": 0.68, "hs09--3": 54},
        # "std_ts_ano_n34e": {"hs09--1": 78, "hs09--2": 0.63, "hs09--3": 51},
    },
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
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
        "panel_1": {"x_delt": 180, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 180, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "cmip5": "forestgreen",
        "cmip5--mean": "forestgreen",
        "cmip5--pi_range": "limegreen",
        "cmip5--iv_on_mean": "darkgreen",
        "cmip6": "r",
        "cmip6--mean": "r",
        "cmip6--pi_range": "r",
        "cmip6--iv_on_mean": "firebrick",
        "observations": "k",
        "observations--mean": "k",
        "full--observations": "k",
        "full_m_iv--cmip5": "forestgreen",
        "full_m_iv--cmip6": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5--mean": "-",
        "cmip5--pi_range": "-.",
        "cmip5--iv_on_mean": "--",
        "cmip6--mean": "-",
        "cmip6--pi_range": "-.",
        "cmip6--iv_on_mean": "--",
        "observations--mean": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            60: {
                "std_pr_ano_n30e": {
                    "panel_1": {
                        "lim": [-50, 150],
                        "nam": "VALUE (%)",
                        "tic": list(range(-50, 151, 50)),
                    },
                    "panel_2": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                },
                "std_pr_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 150],
                        "nam": "VALUE (%)",
                        "tic": list(range(-50, 151, 50)),
                    },
                    "panel_2": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                },
                "std_pr_ano_n40e": {
                    "panel_1": {
                        "lim": [-25, 75],
                        "nam": "VALUE (%)",
                        "tic": list(range(-25, 76, 25)),
                    },
                    "panel_2": {
                        "lim": [-25, 75],
                        "nam": "",
                        "tic": list(range(-25, 76, 25)),
                    },
                },
                "std_ts_ano_n30e": {
                    "panel_1": {
                        "lim": [-30, 30],
                        "nam": "VALUE (%)",
                        "tic": list(range(-30, 31, 15)),
                    },
                    "panel_2": {
                        "lim": [-30, 30],
                        "nam": "",
                        "tic": list(range(-30, 31, 15)),
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [-30, 30],
                        "nam": "VALUE (%)",
                        "tic": list(range(-30, 31, 15)),
                    },
                    "panel_2": {
                        "lim": [-30, 30],
                        "nam": "",
                        "tic": list(range(-30, 31, 15)),
                    },
                },
                "std_ts_ano_n40e": {
                    "panel_1": {
                        "lim": [-30, 30],
                        "nam": "VALUE (%)",
                        "tic": list(range(-30, 31, 15)),
                    },
                    "panel_2": {
                        "lim": [-30, 30],
                        "nam": "",
                        "tic": list(range(-30, 31, 15)),
                    },
                },
                "time": {
                    "panel_1": {
                        "lim": [1860, 2005],
                        "nam": "Year",
                        "nbr_minor": 4,
                        "tic": list(range(1850, 2006, 50)),
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
    },
    # titles
    "fig_titles": {
        # variable
        "ave_pr_val": "change of PR " + r"$\bar{x}$" + " (UNITS)",
        "ave_pr_val_n30e": "change of VALUE (UNITS)",
        "ave_pr_val_n34e": "change of VALUE (UNITS)",
        "ave_pr_val_n40e": "change of VALUE (UNITS)",
        "ave_ts_val": "change of SST $\\bar{x}$ (UNITS)",
        "ave_ts_val_n30e": "change of VALUE (UNITS)",
        "ave_ts_val_n34e": "change of VALUE (UNITS)",
        "ave_ts_val_n40e": "change of VALUE (UNITS)",
        "ske_pr_ano": "change of PRA g$_1$",
        "ske_pr_ano_n30e": "change of VALUE",
        "ske_pr_ano_n34e": "change of VALUE",
        "ske_pr_ano_n40e": "change of VALUE",
        "ske_ts_ano": "change of SSTA g$_1$",
        "ske_ts_ano_n30e": "change of VALUE",
        "ske_ts_ano_n34e": "change of VALUE",
        "ske_ts_ano_n40e": "change of VALUE",
        "std_pr_ano": "relative change of PRA $\\sigma$ (%)",
        "std_pr_ano_n30e": "relative change of VALUE (%)",
        "std_pr_ano_n34e": "relative change of VALUE (%)",
        "std_pr_ano_n40e": "relative change of VALUE (%)",
        "std_ts_ano": "relative change of SSTA $\\sigma$ (%)",
        "std_ts_ano_n30e": "relative change of VALUE (%)",
        "std_ts_ano_n34e": "relative change of VALUE (%)",
        "std_ts_ano_n40e": "relative change of VALUE (%)",
        # statistics
        "full": "THR% full",
        "full_m_iv": "Forced-response\nspread",
        "iv_on_mean": "5-95% CI\non mean(IV)",
        "mean": "mean",
        "pi_range": "piControl range",
        # projects
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
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
def e70_s01_change_multiple_regions_plot(
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
                    for dat, da in values_json[dia][dur][pro][exp].items():
                        put_in_dict(n_ensemble, da["member"].size, pro, exp, dat)
    ens_size = {}
    for pro in list(n_ensemble.keys()):
        for exp in list(n_ensemble[pro].keys()):
            for dat, val in n_ensemble[pro][exp].items():
                exp_o = copy__deepcopy(exp) if exp in ["piControl", "historical", "ssp370"] else \
                    "scenario" + str(exp[-2:])
                put_in_dict(ens_size, val, pro, exp_o, dat)
    for pro in data_projects:
        for exp in data_experiments:
            if check_dict_levels(n_ensemble, pro, exp):
                n_ens = len(list(n_ensemble[pro][exp].keys()))
                n_mem = sum([k for k in n_ensemble[pro][exp].values()])
                val = str(n_ens) + " (" + str(n_mem) + ")"
                del n_ensemble[pro][exp]
                put_in_dict(n_ensemble, val, pro, exp)
                print(pro, exp.ljust(10), "number of ensembles: " + str(n_ens).rjust(2) + \
                      str("(" + str(n_mem) + ")").rjust(8))
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
                        ny = int(dur1.split("_")[0])
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if pro == "observations":
                            da1 = da1.sel(indexers={"year": [k for k in da1["year"].values if k >= 1900 + ny]})
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - ny / 2})
                        # save
                        exp_o = copy__deepcopy(exp1) if exp1 in ["piControl", "historical", "ssp370"] else \
                            "scenario" + str(exp1[-2:])
                        put_in_dict(change_em, da1.compute(), dur1, dia, pro, exp_o, dat)
    print("change_em", list(change_em.keys()))
    # -- compute remaining internal variability
    # remaining variability adjusted with corrected IV for each model
    remaining_iv = {}
    for dur1 in list(change_em.keys()):
        if dur1 == dur2:
            continue
        for dia in list(change_em[dur1].keys()):
            for pro in list(change_em[dur1][dia].keys()):
                for exp in list(change_em[dur1][dia][pro].keys()):
                    if exp == "piControl":
                        continue
                    for dat in list(change_em[dur1][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(change_em[dur1][dia][pro][exp][dat])
                        # ensemble size
                        size = da["member"].size
                        # degree
                        degree = proc_th[dia.split("_")[0]]
                        # compute remaining internal variability
                        _, _, ri = remaining_internal_variability_hs09(da, degree=degree, dim="year", ddof=1)
                        # correction: riv is underestimated as IV is underestimated for small ensembles
                        ri *= (1 / proc_co[dia]["hs09--" + str(degree)]) * size ** -0.2
                        # save
                        put_in_dict(remaining_iv, ri.compute(), dur1, dia, pro, exp, dat)
    print("remaining_iv", list(remaining_iv.keys()))
    # -- compute multimodel mean and std
    # mean
    mean_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # mmm and range
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_m, change_r = {}, {}
    for dur1 in list(mean_change.keys()):
        if dur1 == dur2:
            continue
        for dia in list(mean_change[dur1].keys()):
            for pro in list(mean_change[dur1][dia].keys()):
                for exp in list(mean_change[dur1][dia][pro].keys()):
                    if pro == "observations" or exp == "piControl" or \
                            not check_dict_levels(remaining_iv, dur1, dia, pro, exp):
                        continue
                    # dictionaries
                    d1 = mean_change[dur1][dia][pro][exp]
                    d3 = remaining_iv[dur1][dia][pro][exp]
                    # list datasets
                    list_datasets = list(set(list(d1.keys())) & set(list(d3.keys())))
                    # arrays
                    da1 = concatenate(tuple(d1[dat] for dat in list_datasets), "dataset")
                    da3 = concatenate(tuple(d3[dat] for dat in list_datasets), "dataset")
                    # -- scaled distribution
                    # initial multimodel mean using ensemble mean -> will not change
                    mi = da1.mean(dim="dataset", skipna=True)
                    # initial multimodel spread using ensemble mean -> that what will be scaled
                    si = da1.var(dim="dataset", skipna=True)
                    # initial multimodel mean residual internal variability -> will be used for scaling
                    ri = copy__deepcopy(da3).mean(dim="dataset", skipna=True)
                    # final spread (initial spread minus residual internal variability) -> must be set to 0 if negative
                    sf = (si - ri) / si
                    sf = sf.where(sf > 0, other=0)
                    # scaled distribution
                    da4 = mi + (da1 - mi) * sf
                    # -- influence of RIV on ensemble mean
                    # residual internal variability to compute SE -> sqrt of sum of squares / number of dataset
                    nbr = da3["dataset"].size
                    se = copy__deepcopy(da3).sum(dim="dataset", skipna=True) ** 0.5 / nbr
                    zscore = scipy__stats__t.ppf(0.5 + 90 / 200, nbr - 1)
                    da5 = concatenate((mi - zscore * se, mi + zscore * se), dim="quantile")
                    # save
                    lk = (dur1, dia, exp, pro)
                    put_in_dict(change_m, mi.compute(), *lk)
                    # compute range
                    da_o = statistic(da4, dim="dataset", stat="per", per=percentiles, skipna=True)
                    if "resample" in list(da_o.dims):
                        da_o = da_o.mean(dim="resample", skipna=True)
                    put_in_dict(change_r, da_o.compute(), *lk + ("full_m_iv",))
                    put_in_dict(change_r, da5.compute(), *lk + ("error",))
    print("change_m", list(change_m.keys()))
    # -- add piControl
    ex2 = "piControl"
    for dur in list(change_em.keys()):
        for dia in list(change_em[dur].keys()):
            for pro in list(change_em[dur][dia].keys()):
                if pro == "observations":
                    continue
                for ex1 in list(change_em[dur][dia][pro].keys()):
                    if ex1 == ex2 or not check_dict_levels(change_em, dur, dia, pro, ex2):
                        continue
                    da_o = ()
                    for dat in list(change_em[dur][dia][pro][ex1].keys()):
                        if not check_dict_levels(change_em, dur, dia, pro, ex2, dat):
                            continue
                        # array
                        da = copy__deepcopy(change_em[dur][dia][pro][ex2][dat]).isel(indexers={"year": 0})
                        # ensemble size in ex1
                        n_mem = ens_size[pro][ex1][dat]
                        # space ex2 data
                        spacing = da["member"].values
                        spacing = abs(int(dur.split("_")[0]) / (spacing[1] - spacing[0]))
                        # da = statistic_resampled(
                        #     da, dim="member", pseudo_sample_size=proc_ps, stat="mea", sample_size=n_mem)
                        da = statistic_resampled(da, dim="member", pseudo_sample_size=proc_ps, stat="mea",
                                                 sample_method="spaced", sample_size=n_mem, spacing=spacing)
                        # range
                        da_o += (statistic(da, dim="resample", stat="per", per=[5, 95], skipna=True),)
                    # compute interval in each model
                    da_o = concatenate(da_o, dim="dataset")
                    # average across models
                    da_o = da_o.mean(dim="dataset", skipna=True)
                    da_o = concatenate((da_o, da_o), dim="year")
                    da_o = da_o.assign_coords({"year": numpy__array([1850, 2100])})
                    # save
                    lk = (dur, dia, "--".join([ex1, ex2]), pro, "full")
                    put_in_dict(change_r, copy__deepcopy(da_o).compute(), *lk)
    # -- add observations
    for dur in list(change_em.keys()):
        for dia in list(change_em[dur].keys()):
            for pro in list(change_em[dur][dia].keys()):
                if pro != "observations":
                    continue
                for ex1 in list(change_em[dur][dia][pro].keys()):
                    # concatenate dataset
                    da = concatenate(tuple(da for da in change_em[dur][dia][pro][ex1].values()), dim="dataset")
                    da = da.mean(dim="member", skipna=True)
                    # dataset mean
                    da_m = da.mean(dim="dataset", skipna=True)
                    # dataset range
                    da_r = statistic(da, dim="dataset", stat="per", per=[0, 100], skipna=True)
                    # save
                    list_exp = []
                    if check_dict_levels(change_m, dur, dia):
                        list_exp += list(change_m[dur][dia].keys())
                    list_exp = sorted(list(set(list_exp)), key=lambda s: s.lower())
                    for ex2 in list_exp:
                        put_in_dict(change_m, da_m.compute(), dur, dia, ex2, pro)
                        put_in_dict(change_r, da_r.compute(), dur, dia, ex2, pro, "full")
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia in list(change_r[dur].keys()):
            tmp = []
            for exp in list(change_r[dur][dia].keys()):
                for pro in list(change_r[dur][dia][exp].keys()):
                    for da in change_r[dur][dia][exp][pro].values():
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
        if not check_dict_levels(change_r, dur):
            continue
        list_dia1 = []
        for k in data_diagnostics1:
            if "_".join(k.split("_")[:-1]) not in list_dia1:
                list_dia1.append("_".join(k.split("_")[:-1]))
        ii = -1
        for dia1 in list_dia1:
            list_dia2 = [k for k in data_diagnostics1 if dia1 in k]
            for dia2 in list_dia2:
                ii += 1
                list_exp = sorted(list(set([k for k in list(change_r[dur][dia2].keys()) if "piControl" not in k])),
                                  key=lambda s: s.lower())
                for jj, exp1 in enumerate(list_exp):
                    if not check_dict_levels(change_m, dur, dia2, exp1):
                        continue
                    tp = "panel_1" if exp1 == "historical" else "panel_2"
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, grp, pan)
                    print(dur, dia2, exp1.ljust(10), tuple_k)
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                        for k1, d1 in fig_panel_size.items():
                            if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                                for k2, d2 in d1.items():
                                    if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                        continue
                                    tmp = copy__deepcopy(d2)
                                    tt = (dur, proc_ci, "time", tp, "lim")
                                    if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
                                        k3 = read_from_dict(fig_ticks, *tt)
                                        tmp *= (k3[-1] - k3[0])
                                    put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                            elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], ["time", dia2]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
                            for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if (k3 == "y_nam" and jj > 0) or (
                                        k3 == "x_nam" and dia1 != list_dia1[-1] and dia2 != list_dia2[-1]):
                                    tmp = ""
                                elif k3 == "y_nam" and jj == 0:
                                    tmp = tmp.replace("VALUE", metadata_json[dia2]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dia2]["units"])
                                    tmp = tmp.replace(" ()", "")
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
                        put_in_dict(plot_data, fig_titles[exp1], *tuple_k + ("title_col",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    list_leg_ls = []
                    if check_dict_levels(change_m, dur, dia2, exp1):
                        for pro in list(change_m[dur][dia2][exp1].keys()):
                            tt = "--".join([pro, "mean"])
                            da = copy__deepcopy(change_m[dur][dia2][exp1][pro])
                            lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[tt] if tt in list(fig_linestyles.keys()) else ":"]
                            lw += [2 if pro[:4] == "cmip" else 1]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [8 if pro == "cmip6" else (7 if pro == "cmip5" else 9)]
                            if pro == list(change_m[dur][dia2][exp1].keys())[0] and tt in list(fig_linestyles.keys()):
                                list_leg_ls.append(tt)
                    exp2 = "piControl"
                    exp_t = "--".join([exp1, exp2])
                    if check_dict_levels(change_r, dur, dia2, exp_t):
                        for pro in list(change_r[dur][dia2][exp_t].keys()):
                            if check_dict_levels(change_r, dur, dia2, exp_t, pro, "full"):
                                tt = "--".join([pro, "pi_range"])
                                da = copy__deepcopy(change_r[dur][dia2][exp_t][pro]["full"])
                                lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"] * 2
                                ls += [fig_linestyles[tt] if tt in list(fig_linestyles.keys()) else ":"] * 2
                                lw += [1] * 2
                                lx += [da["year"].values] * 2
                                ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
                                lz += [1] * 2
                                if pro == list(change_r[dur][dia2][exp_t].keys())[0] and \
                                        tt in list(fig_linestyles.keys()):
                                    list_leg_ls.append(tt)
                    # if check_dict_levels(change_r, dur, dia2, exp1):
                    #     for pro in list(change_r[dur][dia2][exp1].keys()):
                    #         if check_dict_levels(change_r, dur, dia2, exp1, pro, "error"):
                    #             tt = "--".join([pro, "iv_on_mean"])
                    #             da = copy__deepcopy(change_r[dur][dia2][exp1][pro]["error"])
                    #             lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"] * 2
                    #             ls += [fig_linestyles[tt] if tt in list(fig_linestyles.keys()) else ":"] * 2
                    #             lw += [1] * 2
                    #             lx += [da["year"].values] * 2
                    #             ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
                    #             lz += [4] * 2
                    #             if pro == list(change_r[dur][dia2][exp1].keys())[0] and \
                    #                     tt in list(fig_linestyles.keys()):
                    #                 list_leg_ls.append(tt)
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
                        if ii == 0 and jj == 0:
                            for kk, mm in enumerate(list_leg_ls):
                                ln = [0] if mm.split("--")[-1] == "mean" else [-1, 1]
                                for nn in ln:
                                    lc += ["k"]
                                    ls += [fig_linestyles[mm]]
                                    lw += [1]
                                    lx += [[x1 + (mm + (8 if kk == 0 and dia2 not in ["ske_pr_ano_n40e"] else 0)) * dx
                                            for mm in [2, 9]]]
                                    ly += [[y1 + ((len(list_leg_ls) - 1 + kk) * 10 + nn * 1.5 + 4) * dy
                                            if dia2 in ["ske_pr_ano_n40e"] else y2 - (kk * 10 + nn * 1.5 + 4) * dy] * 2]
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
                    list_leg_sh = []
                    if check_dict_levels(change_r, dur, dia2, exp1):
                        for pro in list(change_r[dur][dia2][exp1].keys()):
                            for k2, (typ, da) in enumerate(change_r[dur][dia2][exp1][pro].items()):
                                if typ == "error":
                                    continue
                                tt = str(typ) + "--" + str(pro)
                                la += [0.3 + k2 * 0.2]
                                lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"]
                                lx += [da["year"].values]
                                ly1 += [da.min(dim="quantile").values]
                                ly2 += [da.max(dim="quantile").values]
                                lz += [3 if pro == "cmip6" else (2 if pro == "cmip5" else 4)]
                                if pro == list(change_r[dur][dia2][exp1].keys())[0] and \
                                        tt in list(fig_colors.keys()):
                                    list_leg_sh.append(tt)
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if jj == 0:
                            add = len(list_leg_ls)
                            for kk, mm in enumerate(list_leg_sh):
                                la += [0.3 + kk * 0.2]
                                lc += ["k"]
                                lx += [[x1 + (mm + (8 if kk + add == 0 else 0)) * dx for mm in [2, 9]]]
                                ly1 += [[y1 + (kk * 10 + 4 + 3) * dy if dia2 in ["ske_pr_ano_n40e"] else
                                         y2 - ((kk + add) * 10 + 4 + 3) * dy] * 2]
                                ly2 += [[y1 + (kk * 10 + 4 - 3) * dy if dia2 in ["ske_pr_ano_n40e"] else
                                         y2 - ((kk + add) * 10 + 4 - 3) * dy] * 2]
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
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if ii ==0 and jj == 0:
                            add = len(list_leg_sh)
                            for kk, mm in enumerate(list_leg_ls):
                                lc += ["k"]
                                lf += [15]
                                lh += ["left"]
                                lr += [0]
                                lt += [fig_titles[mm.split("--")[-1]]]
                                lv += ["center"]
                                lw += ["regular"]
                                lx += [x1 + (19 if (kk == 0 and dia2 not in ["ske_pr_ano_n40e"]) else 11) * dx]
                                ly += [y1 + ((len(list_leg_ls) - 2 + kk + add) * 7 + 5) * dy
                                       if dia2 in ["ske_pr_ano_n40e"] else y2 - (kk * 8 + 5) * dy]
                            add = len(list_leg_ls)
                            for kk, mm in enumerate(list_leg_sh):
                                tmp = fig_titles[mm.split("--")[0]].split("\n")
                                lc += ["k"] * len(tmp)
                                lf += [14] * len(tmp)
                                lh += ["left"] * len(tmp)
                                lr += [0] * len(tmp)
                                lt += copy__deepcopy(tmp)
                                lv += ["center"] * len(tmp) if len(tmp) == 1 else ["bottom", "top"]
                                lw += ["regular"] * len(tmp)
                                lx += [x1 + (19 if kk + add == 0 else 11) * dx] * len(tmp)
                                ly += [y2 - ((kk + add) * 10 + 4 + 0.4 * (1 - 2 * k4 if len(tmp) == 2 else 0)) * dy
                                       for k4 in range(len(tmp))]
                            for kk, mm in enumerate(list(reversed(data_projects))):
                                lc += [fig_colors[mm] if mm in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left"]
                                lr += [0]
                                lt += [fig_titles[mm] if mm in list(fig_titles.keys()) else mm]
                                lv += ["center"]
                                lw += ["regular"]
                                lx += [x1 + 2 * dx]
                                ly += [y2 - (kk * 8 + 46) * dy]
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
    a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    a5 = list(set(["_".join(k.split("_")[:-1]) for k in data_diagnostics1]))
    a5 = "" if len(a5) > 1 else "_" + str(a5[0])
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"]) + str(a5)
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
e70_s01_change_multiple_regions_plot()
