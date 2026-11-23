# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure 1: evolution of variable change in different experiments
#    * time (x) vs. SSTA std relative change (y) averaged across models
#        - departure from the end of the historical experiment (1965-2014)
#        - MMM-ME for the signal
#        - 1 random member for the spread as in IPCC
#        - ME minus corrected residual internal variability
#    * time (x) vs. PRA std relative change (y) averaged across models
# The distributions don't look normal. In the first code it was the assumption to reconstruct the distribution and
# estimate the range (e.g., 10-90% of models). In this code I will rescale the distribution:
# Df = Mi + (Di - Mi) * Sf / Si
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
from typing import Literal

# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.statistics import scipy__stats__t
from enso_change.wrapper.xarray_base import concatenate, numpy__array, polyval
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
    "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85",
    #                      "ssp126", "ssp245", "ssp585", "ssp370"],
    "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # "data_experiments": ["piControl", "historical", "rcp85", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6", "cmip5", "observations"],  # default_arguments["data_projects"],
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
        # "std_ts_ano_n34e": {"hs09--1": 78, "hs09--2": 0.63, "hs09--3": 51}
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
        "cmip6": "r",
        "observations": "k",
        "full--cmip5": "brown",
        "full--cmip6": "peru",
        "full--observations": "k",
        "full_m_iv--cmip5": "forestgreen",
        "full_m_iv--cmip6": "r",
        "iv_on_mean--cmip5": "darkgreen",
        "iv_on_mean--cmip6": "firebrick",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        "observations": "--",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            60: {
                "std_pr_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                    "panel_2": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                },
                "std_ts_ano_n30e": {
                    "panel_1": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
                    },
                    "panel_2": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [-34, 34],
                        "nam": "",
                        "tic": list(range(-34, 35, 17)),
                    },
                    "panel_2": {
                        "lim": [-34, 34],
                        "nam": "",
                        "tic": list(range(-34, 35, 17)),
                    },
                },
                "std_ts_ano_n40e": {
                    "panel_1": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
                    },
                    "panel_2": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
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
            70: {
                "std_pr_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                    "panel_2": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
                    },
                    "panel_2": {
                        "lim": [-40, 40],
                        "nam": "",
                        "tic": list(range(-40, 41, 20)),
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
            80: {
                "std_pr_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                    "panel_2": {
                        "lim": [-50, 150],
                        "nam": "",
                        "tic": list(range(-50, 151, 50)),
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
                    },
                    "panel_2": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
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
            90: {
                "std_pr_ano_n34e": {
                    "panel_1": {
                        "lim": [-70, 210],
                        "nam": "",
                        "tic": list(range(-70, 211, 70)),
                    },
                    "panel_2": {
                        "lim": [-70, 210],
                        "nam": "",
                        "tic": list(range(-70, 211, 70)),
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
                    },
                    "panel_2": {
                        "lim": [-50, 50],
                        "nam": "",
                        "tic": list(range(-50, 51, 25)),
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
        "ave_pr_val_n34e": "change of VALUE (UNITS)",
        "ske_pr_ano_n34e": "change of VALUE",
        "std_pr_ano_n34e": "relative change of VALUE (%)",
        "ave_ts_val_n34e": "change of VALUE (UNITS)",
        "ske_ts_ano_n34e": "change of VALUE",
        "std_ts_ano_n30e": "relative change of VALUE (%)",
        # "std_ts_ano_n34e": "relative change of VALUE (%)",
        "std_ts_ano_n34e": "$\\Delta$E relative to 1956-2005 (%)",
        "std_ts_ano_n40e": "relative change of VALUE (%)",
        "full--cmip5": "Raw model spread",  # "THR% full",
        "full--cmip6": "Raw model spread",  # "THR% full",
        "full--observations": "obs. range",
        "full_m_iv--cmip5": "Forced-response spread",  # "THR% no IV",
        "full_m_iv--cmip6": "Forced-response spread",  # "THR% no IV",
        "mean--cmip5": "MMM",  # "CMIP5 MMM",
        "mean--cmip6": "MMM",  # "CMIP6 MMM",
        "mean--observations": "obs. mean",
        # "figure_title": "Multimodel spread with or without internal variability",
        "figure_title": "ENSO amplitude change with and without internal variability",
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
        "title_row_x": 20,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def f02a_change_plot(
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
                            da1 = da1.sel(indexers={"year": [k for k in da1["year"].values if k > 1950 + ny]})
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        yy = 2006
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # save
                        exp_o = copy__deepcopy(exp1) if exp1 in ["piControl", "historical", "ssp370"] else \
                            "scenario" + str(exp1[-2:])
                        put_in_dict(change_em, da1.compute(), dur1, dia, pro, exp_o, dat)
    print("change_em", list(change_em.keys()))
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
                        yy = 2006
                        da2 = da2.sel(indexers={"year": yy})
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # random selection
                        da1 = statistic_resampled(da1, dim="member", pseudo_sample_size=proc_ps, sample_size=1)
                        # save
                        exp_o = copy__deepcopy(exp) if exp in ["piControl", "historical", "ssp370"] else \
                            "scenario" + str(exp[-2:])
                        put_in_dict(change_1m, da1.compute(), dur1, dia, pro, exp_o, dat)
    print("change_1m", list(change_1m.keys()))
    # -- compute residual internal variability
    # residual variability adjusted with corrected IV for each model
    residual_iv = {}
    for dur1 in list(change_em.keys()):
        if dur1 == dur2:
            continue
        for dia in list(change_em[dur1].keys()):
            for pro in list(change_em[dur1][dia].keys()):
                if pro == "observations":
                    continue
                for ex1 in list(change_em[dur1][dia][pro].keys()):
                    if ex1 == "piControl":
                        continue
                    for dat in list(change_em[dur1][dia][pro][ex1].keys()):
                        # array
                        da_em = copy__deepcopy(change_em[dur1][dia][pro][ex1][dat])
                        size = da_em["member"].size
                        if "member" in list(da_em.dims):
                            da_em = da_em.mean(dim="member")
                        # compute residual internal variability (residual from polynomial)
                        # order of the polynomial
                        nn = proc_th[dia.split("_")[0]]
                        # 'forced change' is a Nth order polynomial fit computed on the ensemble mean
                        coef = da_em.polyfit("year", nn)
                        da_p = polyval(da_em["year"], coef["polyfit_coefficients"], degree_dim="degree")
                        # 'residual internal variability' is the temporal std of the residual (da - da_p)
                        riv = (da_em - da_p).var(dim="year", ddof=1)
                        # correction: riv is underestimated as IV is underestimated for small ensembles
                        riv *= (1 / proc_co[dia]["hs09--" + str(nn)]) * size**-0.2
                        # save
                        put_in_dict(residual_iv, copy__deepcopy(riv), dur1, dia, pro, ex1, dat)
    print("residual_iv", list(residual_iv.keys()))
    # -- compute multimodel mean and std
    # mean
    mean_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # mmm and range
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_m, change_r = {}, {}
    for dur in list(mean_change.keys()):
        if dur == dur2:
            continue
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    if pro == "observations" or exp == "piControl" or \
                            not check_dict_levels(change_1m, dur, dia, pro, exp) or \
                            not check_dict_levels(change_em, dur, dia, pro, exp) or \
                            not check_dict_levels(residual_iv, dur, dia, pro, exp):
                        continue
                    # dictionaries
                    d1 = mean_change[dur][dia][pro][exp]
                    d2 = change_1m[dur][dia][pro][exp]
                    d3 = residual_iv[dur][dia][pro][exp]
                    # list datasets
                    list_datasets = list(set(list(d1.keys())) & set(list(d2.keys())) & set(list(d3.keys())))
                    # arrays
                    da1 = concatenate(tuple(d1[dat] for dat in list_datasets), "dataset")
                    da2 = concatenate(tuple(d2[dat] for dat in list_datasets), "dataset")
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
                    put_in_dict(change_m, mi, dur, dia, pro, exp)
                    # compute range
                    lk = (dur, dia, pro, exp)
                    for tt, dd in zip(["full", "full_m_iv"], [da2, da4]):
                        da_o = statistic(dd, dim="dataset", stat="per", per=percentiles, skipna=True)
                        if "resample" in list(da_o.dims):
                            da_o = da_o.mean(dim="resample", skipna=True)
                        put_in_dict(change_r, da_o.compute(), *lk + (tt,))
                    put_in_dict(change_r, da5.compute(), *lk + ("error",))
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
                        da_o += (statistic(da, dim="resample", stat="per", per=percentiles, skipna=True),)
                    # compute interval in each model
                    da_o = concatenate(da_o, dim="dataset")
                    # average across models
                    da_o = da_o.mean(dim="dataset", skipna=True)
                    da_o = concatenate((da_o, da_o), dim="year")
                    da_o = da_o.assign_coords({"year": numpy__array([1850, 2100])})
                    # save
                    lk = (dur, dia, pro, "--".join([ex1, ex2]), "full")
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
                        for k1 in list(change_m[dur][dia].keys()):
                            list_exp += list(change_m[dur][dia][k1].keys())
                    list_exp = sorted(list(set(list_exp)), key=lambda s: s.lower())
                    for ex2 in list_exp:
                        put_in_dict(change_m, da_m.compute(), dur, dia, pro, ex2)
                        put_in_dict(change_r, da_r.compute(), dur, dia, pro, ex2, "full")
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia in list(change_r[dur].keys()):
            tmp = []
            for pro in list(change_r[dur][dia].keys()):
                for exp in list(change_r[dur][dia][pro].keys()):
                    for da in change_r[dur][dia][pro][exp].values():
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
        ny1 = int(dur.split("_")[0])
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
        list_dia = [k for k in data_diagnostics1 if k in list_dia]
        list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
        list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
        ii = -1
        for dia in list_dia:
            for pro1 in list_pro:
                ii += 1
                list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
                                  key=lambda s: s.lower())
                for jj, exp in enumerate(list_exp):
                    if not check_dict_levels(change_m, dur, dia, pro1, exp):
                        continue
                    tp = "panel_1" if exp == "historical" else "panel_2"
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, grp, pan)
                    print(dur, dia, pro1, exp.ljust(10), tuple_k)
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
                    l1, l2 = ["x_", "y_"], ["time", dia]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
                            for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                    yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
                                    tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
                                    tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                    tmp = tmp.replace(" ()", "")
                                if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
                                                                  not (dia == list_dia[-1] and
                                                                       (pro1 == list_pro[-1] or exp == list_exp[-1]))):
                                    tmp = ""
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
                        put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    for pro2 in [pro1, "observations"]:
                        if check_dict_levels(change_m, dur, dia, pro2, exp):
                            da = copy__deepcopy(change_m[dur][dia][pro2][exp])
                            lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
                            lw += [2 if pro2[:4] == "cmip" else 1]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [4 if pro2[:4] == "cmip" else 5]
                    exp_t = "--".join([exp, "piControl"])
                    if check_dict_levels(change_r, dur, dia, pro1, exp_t, "full"):
                        da = copy__deepcopy(change_r[dur][dia][pro1][exp_t]["full"])
                        lc += [fig_colors["piControl"] if "piControl" in list(fig_colors.keys()) else "cyan"] * 2
                        ls += [fig_linestyles["piControl"] if "piControl" in list(fig_linestyles.keys()) else ":"] * 2
                        lw += [2] * 2
                        lx += [da["year"].values] * 2
                        ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
                        lz += [1] * 2
                    if check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
                        tmp = "iv_on_mean--" + str(pro1)
                        da = copy__deepcopy(change_r[dur][dia][pro1][exp]["error"])
                        lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
                        ls += ["--"] * 2
                        lw += [1] * 2
                        lx += [da["year"].values] * 2
                        ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
                        lz += [4] * 2
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
                            k1 = 0
                            for pro2 in [pro1, "observations"]:
                                if check_dict_levels(change_m, dur, dia, pro2, exp):
                                    t1 = pro2[:4] == "cmip"
                                    t2 = pro2[:4] != "cmip" and dia == "std_pr_ano_n34e"
                                    lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                                    ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
                                    lw += [2 if t1 else 1]
                                    lx += [[x1 + 10 * dx, x1 + 18 * dx] if t1 else [x2 - 10 * dx, x2 - 2 * dx]]
                                    ly += [[y2 - (k1 * 8 + 5) * dy if t1 else (
                                        y2 - ((k1 + 2) * 8 + 5) * dy if t2 else y1 + (k1 * 7 + 4) * dy)] * 2]
                                    lz += [1]
                                    k1 += 1
                        # if jj == 1 and check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
                        if ii == 0 and jj == 0 and check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
                            tmp = "iv_on_mean--" + str(pro1)
                            lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
                            ls += ["--"] * 2
                            lw += [1] * 2
                            lx += [[x1 + 10 * dx, x1 + 16 * dx]] * 2
                            ly += [[y2 - (5 + k * 1.5) * dy] * 2 for k in [-1, 1]]
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
                    for pro2 in [pro1, "observations"]:
                        if check_dict_levels(change_r, dur, dia, pro2, exp):
                            for k2, (typ, da) in enumerate(change_r[dur][dia][pro2][exp].items()):
                                if typ == "error":
                                    continue
                                vv = str(typ) + "--" + str(pro2)
                                la += [0.3 + k2 * 0.2]
                                lc += [fig_colors[vv] if vv in list(fig_colors.keys()) else "cyan"]
                                lx += [da["year"].values]
                                ly1 += [da.min(dim="quantile").values]
                                ly2 += [da.max(dim="quantile").values]
                                lz += [2 if pro2[:4] == "cmip" else 3]
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if ii == 0 and jj == 0:
                            k1 = 0
                            for pro2 in [pro1, "observations"]:
                                if check_dict_levels(change_r, dur, dia, pro2, exp):
                                    for k2, typ in enumerate(list(change_r[dur][dia][pro2][exp].keys())):
                                        if typ == "error":
                                            continue
                                        t1 = pro2[:4] == "cmip"
                                        t2 = pro2[:4] != "cmip" and dia == "std_pr_ano_n34e"
                                        vv = str(typ) + "--" + str(pro2)
                                        la += [0.3 + k2 * 0.2]
                                        lc += [fig_colors[vv] if vv in list(fig_colors.keys()) else "cyan"]
                                        lx += [[x1 + 2 * dx, x1 + 10 * dx] if t1 else [x2 - 10 * dx, x2 - 2 * dx]]
                                        ly1 += [[y2 - ((k2 + k1 * 3 + 1) * 8 + 5 + 3) * dy if t1 or t2 else
                                                 y1 + (5 + 3) * dy] * 2]
                                        ly2 += [[y2 - ((k2 + k1 * 3 + 1) * 8 + 5 - 3) * dy if t1 or t2 else
                                                 y1 + (4 - 3) * dy] * 2]
                                        lz += [1]
                                    k1 += 1
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
                        if ii == 0 and jj == 0:
                            for k1, pro2 in enumerate([pro1, "observations"]):
                                l1 = ["mean--" + str(pro2)]
                                if check_dict_levels(change_r, dur, dia, pro2, exp):
                                    l1 += [str(k2) + "--" + str(pro2) for k2 in
                                           list(change_r[dur][dia][pro2][exp].keys()) if k2 != "error"]
                                for k2, typ in enumerate(l1):
                                    t1 = pro2[:4] == "cmip"
                                    t2 = pro2[:4] != "cmip" and dia == "std_pr_ano_n34e"
                                    lc += ["k"]
                                    lf += [15]
                                    lh += ["left" if t1 else "right"]
                                    lr += [0]
                                    val = [str(int(k3)) if k3 == int(k3) else str(k3) for k3 in percentiles]
                                    lt += [fig_titles[typ].replace("THR", "-".join(val))]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x1 + (20 if k2 == 0 else 12) * dx if t1 else x2 - 12 * dx]
                                    ly += [y2 - ((k2 + k1 * 3) * 8 + 5) * dy if t1 or t2 else
                                           y1 + (7 * (len(l1) - k2 - 1) + 4) * dy]
                        if ii == 0 and jj == 0 and check_dict_levels(change_r, dur, dia, pro1, exp_t, "full"):
                            lc += [fig_colors["piControl"]]
                            lf += [15]
                            lh += ["center"]
                            lr += [0]
                            lt += ["piControl range%"]  # ["piControl\nrange 5-95%"]
                            lv += ["bottom"]
                            lw += ["regular"]
                            lx += [(x1 + x2) / 2]
                            ly += [change_r[dur][dia][pro1][exp_t]["full"].min().values + dy]
                        if ii == 0 and jj == 0 and check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
                            tmp = "iv_on_mean--" + str(pro1)
                            lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            lt += ["5–95% on MMM"]  # ["5-95% CI\non MMM(IV)"]
                            lv += ["top"]
                            lw += ["regular"]
                            lx += [x2 - 2 * dx]
                            ly += [y2 - 2 * dy]
                        # project
                        # if jj == 0:
                        #     lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
                        #     lf += [15]
                        #     lh += ["center"]
                        #     lr += [90]
                        #     lt += [fig_titles[pro1]]
                        #     lv += ["center"]
                        #     lw += ["regular"]
                        #     lx += [x1 - 20 * dx]
                        #     ly += [y1 + 50 * dy]
                        # variable
                        if jj == 0 and pro1 == list_pro[0]:
                            # plot size
                            y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
                            # add to text
                            lc += ["k"]
                            lf += [15]
                            lh += ["center"]
                            lr += [90]
                            val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                            lt += [val.replace("UNITS", metadata_json[dia]["units"])]
                            lv += ["center"]
                            lw += ["bold"]
                            tmp = panel_param["title_row_x"] if "title_row_x" in list(panel_param.keys()) else 30
                            lx += [x1 - tmp * dx]
                            ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
                        # ensemble size
                        # if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
                        #     k1 = 0
                        #     for pro2 in list(n_ensemble.keys()):
                        #         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                        #         lf += [15]
                        #         lh += ["left"]
                        #         lr += [0]
                        #         lt += [fig_titles[pro2]]
                        #         lv += ["top"]
                        #         lw += ["regular"]
                        #         lx += [x2 + 5 * dx]
                        #         ly += [y2 - (k1 * 8 + 10) * dy]
                        #         k1 += 1
                        #         for exp2 in list(n_ensemble[pro2].keys()):
                        #             if exp2 == "piControl":
                        #                 continue
                        #             n_ens = n_ensemble[pro2][exp2].split(" ")[0]
                        #             n_mem = n_ensemble[pro2][exp2].split(" ")[1]
                        #             lc += ["k"] * 3
                        #             lf += [15] * 3
                        #             lh += ["left", "right", "right"]
                        #             lr += [0] * 3
                        #             lt += [fig_titles[exp2], n_ens, n_mem]
                        #             lv += ["top"] * 3
                        #             lw += ["regular"] * 3
                        #             lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
                        #             ly += [y2 - (k1 * 8 + 10) * dy] * 3
                        #             k1 += 1
                        # plot title
                        if ii == 0 and jj == 0:
                            # plot size
                            l1 = len([k1 for k1 in list_exp if k1 == "historical"])
                            x1_d = fig_panel_size[tp]["x_delt"]
                            x1_m = fig_panel_size[tp]["x_size"]
                            x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
                            x1_s = max(x1_s) - min(x1_s)
                            l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
                            x2_d = fig_panel_size["panel_2"]["x_delt"]
                            x2_m = fig_panel_size["panel_2"]["x_size"]
                            x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
                            x2_s = max(x2_s) - min(x2_s)
                            # add to text
                            lc += ["k"]
                            lf += [19]
                            lh += ["center"]
                            lr += [0]
                            lt += [fig_titles["figure_title"]]
                            lv += ["center"]
                            lw += ["bold"]
                            a1 = x1_m * x1_s
                            a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
                            lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
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
    a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[dur].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "t1")))
        # plot
        fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
                  figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # # panel 1: time series
    # print("time series 1")
    # group, panel = "group_1", "panel_1"
    # for dur in list(change_m.keys()):
    #     if not check_dict_levels(change_r, dur):
    #         continue
    #     ny1 = int(dur.split("_")[0])
    #     ny2 = int(data_epoch_lengths[-1].split("_")[0])
    #     list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
    #     list_dia = [k for k in data_diagnostics1 if k in list_dia]
    #     list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
    #     list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
    #     ii = -1
    #     for dia in list_dia:
    #         for pro1 in list_pro:
    #             ii += 1
    #             list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
    #                               key=lambda s: s.lower())
    #             for jj, exp in enumerate(list_exp):
    #                 if not check_dict_levels(change_m, dur, dia, pro1, exp):
    #                     continue
    #                 tp = "panel_1" if exp == "historical" else "panel_2"
    #                 pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
    #                 grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
    #                 tuple_k = (dur, grp, pan)
    #                 print(dur, dia, pro1, exp.ljust(10), tuple_k)
    #                 # panel size
    #                 tk = tuple_k[:-2]
    #                 if not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                     for k1, d1 in fig_panel_size.items():
    #                         if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                             for k2, d2 in d1.items():
    #                                 if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
    #                                     continue
    #                                 tmp = copy__deepcopy(d2)
    #                                 tt = (dur, proc_ci, "time", tp, "lim")
    #                                 if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
    #                                     k3 = read_from_dict(fig_ticks, *tt)
    #                                     tmp *= (k3[-1] - k3[0])
    #                                 put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
    #                         elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
    #                             put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
    #                 # panel axes definition
    #                 l1, l2 = ["x_", "y_"], ["time", dia]
    #                 for vv, dd in zip(l1, l2):
    #                     if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
    #                         for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
    #                             k3 = str(vv) + str(k1)
    #                             if check_dict_levels(figure_axes, *tuple_k + (k3,)):
    #                                 continue
    #                             tmp = copy__deepcopy(k2)
    #                             if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
    #                                 yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
    #                                 tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
    #                                 tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
    #                                 tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
    #                                 tmp = tmp.replace(" ()", "")
    #                             if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
    #                                                               not (dia == list_dia[-1] and
    #                                                                    (pro1 == list_pro[-1] or exp == list_exp[-1]))):
    #                                 tmp = ""
    #                             put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
    #                 dt = {}
    #                 if check_dict_levels(figure_axes, *tuple_k):
    #                     dt = read_from_dict(figure_axes, *tuple_k)
    #                 # test possibility to plot text / legend
    #                 do_legend = False
    #                 if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
    #                         dt["y_lim"] is not None:
    #                     do_legend = True
    #                 # -- title
    #                 if ii == 0:
    #                     put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
    #                 # -- curve
    #                 lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_m, dur, dia, pro2, exp):
    #                         da = copy__deepcopy(change_m[dur][dia][pro2][exp])
    #                         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                         ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
    #                         lw += [2 if pro2[:4] == "cmip" else 1]
    #                         lx += [da["year"].values]
    #                         ly += [da.values]
    #                         lz += [4 if pro2[:4] == "cmip" else 5]
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # x = 0 and y = 0
    #                     for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
    #                         if k1[0] < 0 < k1[1]:
    #                             lc += ["k"]
    #                             ls += [":"]
    #                             lw += [1]
    #                             lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
    #                             ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
    #                             lz += [10]
    #                 plot_type = "cur"
    #                 if len(lx) > 0 and len(lx) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- text
    #                 lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # legend
    #                     # project
    #                     if jj == 0:
    #                         lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         lt += [fig_titles[pro1]]
    #                         lv += ["center"]
    #                         lw += ["regular"]
    #                         lx += [x1 - 20 * dx]
    #                         ly += [y1 + 50 * dy]
    #                     # variable
    #                     if jj == 0 and pro1 == list_pro[0]:
    #                         # plot size
    #                         y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
    #                         lt += [val.replace("UNITS", metadata_json[dia]["units"])]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         lx += [x1 - 30 * dx]
    #                         ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
    #                     # ensemble size
    #                     if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
    #                         k1 = 0
    #                         for pro2 in list(n_ensemble.keys()):
    #                             lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                             lf += [15]
    #                             lh += ["left"]
    #                             lr += [0]
    #                             lt += [fig_titles[pro2]]
    #                             lv += ["top"]
    #                             lw += ["regular"]
    #                             lx += [x2 + 5 * dx]
    #                             ly += [y2 - (k1 * 8 + 10) * dy]
    #                             k1 += 1
    #                             for exp2 in list(n_ensemble[pro2].keys()):
    #                                 if exp2 == "piControl":
    #                                     continue
    #                                 n_ens = n_ensemble[pro2][exp2].split(" ")[0]
    #                                 n_mem = n_ensemble[pro2][exp2].split(" ")[1]
    #                                 lc += ["k"] * 3
    #                                 lf += [15] * 3
    #                                 lh += ["left", "right", "right"]
    #                                 lr += [0] * 3
    #                                 lt += [fig_titles[exp2], n_ens, n_mem]
    #                                 lv += ["top"] * 3
    #                                 lw += ["regular"] * 3
    #                                 lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
    #                                 ly += [y2 - (k1 * 8 + 10) * dy] * 3
    #                                 k1 += 1
    #                     # plot title
    #                     if ii == 0 and jj == 0:
    #                         # plot size
    #                         l1 = len([k1 for k1 in list_exp if k1 == "historical"])
    #                         x1_d = fig_panel_size[tp]["x_delt"]
    #                         x1_m = fig_panel_size[tp]["x_size"]
    #                         x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
    #                         x1_s = max(x1_s) - min(x1_s)
    #                         l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
    #                         x2_d = fig_panel_size["panel_2"]["x_delt"]
    #                         x2_m = fig_panel_size["panel_2"]["x_size"]
    #                         x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
    #                         x2_s = max(x2_s) - min(x2_s)
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [20]
    #                         lh += ["center"]
    #                         lr += [0]
    #                         lt += [fig_titles["figure_title"]]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         a1 = x1_m * x1_s
    #                         a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
    #                         lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
    #                         ly += [y2 + 20 * dy]
    #                 plot_type = "text"
    #                 if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
    #                     put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
    #                     put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
    #                     put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # #
    # # -- Plot
    # #
    # # output plot directory (relative to current file directory)
    # figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # # figure name
    # a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    # a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # # figure loop
    # for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
    #     # list plot names
    #     list_names = list(plot_data[dur].keys())
    #     # path to output figure
    #     figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "01")))
    #     # plot
    #     fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
    #               figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
    # figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # # panel 1: time series
    # print("time series 2")
    # group, panel = "group_1", "panel_1"
    # for dur in list(change_m.keys()):
    #     if not check_dict_levels(change_r, dur):
    #         continue
    #     ny1 = int(dur.split("_")[0])
    #     ny2 = int(data_epoch_lengths[-1].split("_")[0])
    #     list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
    #     list_dia = [k for k in data_diagnostics1 if k in list_dia]
    #     list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
    #     list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
    #     ii = -1
    #     for dia in list_dia:
    #         for pro1 in list_pro:
    #             ii += 1
    #             list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
    #                               key=lambda s: s.lower())
    #             for jj, exp in enumerate(list_exp):
    #                 if not check_dict_levels(change_m, dur, dia, pro1, exp):
    #                     continue
    #                 tp = "panel_1" if exp == "historical" else "panel_2"
    #                 pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
    #                 grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
    #                 tuple_k = (dur, grp, pan)
    #                 print(dur, dia, pro1, exp.ljust(10), tuple_k)
    #                 # panel size
    #                 tk = tuple_k[:-2]
    #                 if not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                     for k1, d1 in fig_panel_size.items():
    #                         if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                             for k2, d2 in d1.items():
    #                                 if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
    #                                     continue
    #                                 tmp = copy__deepcopy(d2)
    #                                 tt = (dur, proc_ci, "time", tp, "lim")
    #                                 if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
    #                                     k3 = read_from_dict(fig_ticks, *tt)
    #                                     tmp *= (k3[-1] - k3[0])
    #                                 put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
    #                         elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
    #                             put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
    #                 # panel axes definition
    #                 l1, l2 = ["x_", "y_"], ["time", dia]
    #                 for vv, dd in zip(l1, l2):
    #                     if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
    #                         for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
    #                             k3 = str(vv) + str(k1)
    #                             if check_dict_levels(figure_axes, *tuple_k + (k3,)):
    #                                 continue
    #                             tmp = copy__deepcopy(k2)
    #                             if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
    #                                 yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
    #                                 tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
    #                                 tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
    #                                 tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
    #                                 tmp = tmp.replace(" ()", "")
    #                             if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
    #                                                               not (dia == list_dia[-1] and
    #                                                                    (pro1 == list_pro[-1] or exp == list_exp[-1]))):
    #                                 tmp = ""
    #                             put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
    #                 dt = {}
    #                 if check_dict_levels(figure_axes, *tuple_k):
    #                     dt = read_from_dict(figure_axes, *tuple_k)
    #                 # test possibility to plot text / legend
    #                 do_legend = False
    #                 if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
    #                         dt["y_lim"] is not None:
    #                     do_legend = True
    #                 # -- title
    #                 if ii == 0:
    #                     put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
    #                 # -- curve
    #                 lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_m, dur, dia, pro2, exp):
    #                         da = copy__deepcopy(change_m[dur][dia][pro2][exp])
    #                         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                         ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
    #                         lw += [2 if pro2[:4] == "cmip" else 1]
    #                         lx += [da["year"].values]
    #                         ly += [da.values]
    #                         lz += [4 if pro2[:4] == "cmip" else 5]
    #                 if check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
    #                     tmp = "iv_on_mean--" + str(pro1)
    #                     da = copy__deepcopy(change_r[dur][dia][pro1][exp]["error"])
    #                     lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
    #                     ls += ["--"] * 2
    #                     lw += [1] * 2
    #                     lx += [da["year"].values] * 2
    #                     ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
    #                     lz += [4] * 2
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # x = 0 and y = 0
    #                     for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
    #                         if k1[0] < 0 < k1[1]:
    #                             lc += ["k"]
    #                             ls += [":"]
    #                             lw += [1]
    #                             lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
    #                             ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
    #                             lz += [10]
    #                 plot_type = "cur"
    #                 if len(lx) > 0 and len(lx) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- text
    #                 lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # legend
    #                     # project
    #                     if jj == 0:
    #                         lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         lt += [fig_titles[pro1]]
    #                         lv += ["center"]
    #                         lw += ["regular"]
    #                         lx += [x1 - 20 * dx]
    #                         ly += [y1 + 50 * dy]
    #                     # variable
    #                     if jj == 0 and pro1 == list_pro[0]:
    #                         # plot size
    #                         y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
    #                         lt += [val.replace("UNITS", metadata_json[dia]["units"])]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         lx += [x1 - 30 * dx]
    #                         ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
    #                     # ensemble size
    #                     if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
    #                         k1 = 0
    #                         for pro2 in list(n_ensemble.keys()):
    #                             lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                             lf += [15]
    #                             lh += ["left"]
    #                             lr += [0]
    #                             lt += [fig_titles[pro2]]
    #                             lv += ["top"]
    #                             lw += ["regular"]
    #                             lx += [x2 + 5 * dx]
    #                             ly += [y2 - (k1 * 8 + 10) * dy]
    #                             k1 += 1
    #                             for exp2 in list(n_ensemble[pro2].keys()):
    #                                 if exp2 == "piControl":
    #                                     continue
    #                                 n_ens = n_ensemble[pro2][exp2].split(" ")[0]
    #                                 n_mem = n_ensemble[pro2][exp2].split(" ")[1]
    #                                 lc += ["k"] * 3
    #                                 lf += [15] * 3
    #                                 lh += ["left", "right", "right"]
    #                                 lr += [0] * 3
    #                                 lt += [fig_titles[exp2], n_ens, n_mem]
    #                                 lv += ["top"] * 3
    #                                 lw += ["regular"] * 3
    #                                 lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
    #                                 ly += [y2 - (k1 * 8 + 10) * dy] * 3
    #                                 k1 += 1
    #                     # plot title
    #                     if ii == 0 and jj == 0:
    #                         # plot size
    #                         l1 = len([k1 for k1 in list_exp if k1 == "historical"])
    #                         x1_d = fig_panel_size[tp]["x_delt"]
    #                         x1_m = fig_panel_size[tp]["x_size"]
    #                         x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
    #                         x1_s = max(x1_s) - min(x1_s)
    #                         l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
    #                         x2_d = fig_panel_size["panel_2"]["x_delt"]
    #                         x2_m = fig_panel_size["panel_2"]["x_size"]
    #                         x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
    #                         x2_s = max(x2_s) - min(x2_s)
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [20]
    #                         lh += ["center"]
    #                         lr += [0]
    #                         lt += [fig_titles["figure_title"]]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         a1 = x1_m * x1_s
    #                         a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
    #                         lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
    #                         ly += [y2 + 20 * dy]
    #                 plot_type = "text"
    #                 if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
    #                     put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
    #                     put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
    #                     put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # #
    # # -- Plot
    # #
    # # output plot directory (relative to current file directory)
    # figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # # figure name
    # a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    # a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # # figure loop
    # for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
    #     # list plot names
    #     list_names = list(plot_data[dur].keys())
    #     # path to output figure
    #     figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "02")))
    #     # plot
    #     fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
    #               figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
    # figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # # panel 1: time series
    # print("time series 3")
    # group, panel = "group_1", "panel_1"
    # for dur in list(change_m.keys()):
    #     if not check_dict_levels(change_r, dur):
    #         continue
    #     ny1 = int(dur.split("_")[0])
    #     ny2 = int(data_epoch_lengths[-1].split("_")[0])
    #     list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
    #     list_dia = [k for k in data_diagnostics1 if k in list_dia]
    #     list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
    #     list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
    #     ii = -1
    #     for dia in list_dia:
    #         for pro1 in list_pro:
    #             ii += 1
    #             list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
    #                               key=lambda s: s.lower())
    #             for jj, exp in enumerate(list_exp):
    #                 if not check_dict_levels(change_m, dur, dia, pro1, exp):
    #                     continue
    #                 tp = "panel_1" if exp == "historical" else "panel_2"
    #                 pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
    #                 grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
    #                 tuple_k = (dur, grp, pan)
    #                 print(dur, dia, pro1, exp.ljust(10), tuple_k)
    #                 # panel size
    #                 tk = tuple_k[:-2]
    #                 if not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                     for k1, d1 in fig_panel_size.items():
    #                         if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                             for k2, d2 in d1.items():
    #                                 if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
    #                                     continue
    #                                 tmp = copy__deepcopy(d2)
    #                                 tt = (dur, proc_ci, "time", tp, "lim")
    #                                 if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
    #                                     k3 = read_from_dict(fig_ticks, *tt)
    #                                     tmp *= (k3[-1] - k3[0])
    #                                 put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
    #                         elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
    #                             put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
    #                 # panel axes definition
    #                 l1, l2 = ["x_", "y_"], ["time", dia]
    #                 for vv, dd in zip(l1, l2):
    #                     if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
    #                         for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
    #                             k3 = str(vv) + str(k1)
    #                             if check_dict_levels(figure_axes, *tuple_k + (k3,)):
    #                                 continue
    #                             tmp = copy__deepcopy(k2)
    #                             if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
    #                                 yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
    #                                 tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
    #                                 tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
    #                                 tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
    #                                 tmp = tmp.replace(" ()", "")
    #                             if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
    #                                                               not (dia == list_dia[-1] and
    #                                                                    (pro1 == list_pro[-1] or exp == list_exp[-1]))):
    #                                 tmp = ""
    #                             put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
    #                 dt = {}
    #                 if check_dict_levels(figure_axes, *tuple_k):
    #                     dt = read_from_dict(figure_axes, *tuple_k)
    #                 # test possibility to plot text / legend
    #                 do_legend = False
    #                 if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
    #                         dt["y_lim"] is not None:
    #                     do_legend = True
    #                 # -- title
    #                 if ii == 0:
    #                     put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
    #                 # -- curve
    #                 lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_m, dur, dia, pro2, exp):
    #                         da = copy__deepcopy(change_m[dur][dia][pro2][exp])
    #                         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                         ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
    #                         lw += [2 if pro2[:4] == "cmip" else 1]
    #                         lx += [da["year"].values]
    #                         ly += [da.values]
    #                         lz += [4 if pro2[:4] == "cmip" else 5]
    #                 if check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
    #                     tmp = "iv_on_mean--" + str(pro1)
    #                     da = copy__deepcopy(change_r[dur][dia][pro1][exp]["error"])
    #                     lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
    #                     ls += ["--"] * 2
    #                     lw += [1] * 2
    #                     lx += [da["year"].values] * 2
    #                     ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
    #                     lz += [4] * 2
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # x = 0 and y = 0
    #                     for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
    #                         if k1[0] < 0 < k1[1]:
    #                             lc += ["k"]
    #                             ls += [":"]
    #                             lw += [1]
    #                             lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
    #                             ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
    #                             lz += [10]
    #                 plot_type = "cur"
    #                 if len(lx) > 0 and len(lx) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- shading
    #                 la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_r, dur, dia, pro2, exp):
    #                         for k2, (typ, da) in enumerate(change_r[dur][dia][pro2][exp].items()):
    #                             if typ != "full":
    #                                 continue
    #                             vv = str(typ) + "--" + str(pro2)
    #                             la += [0.3 + k2 * 0.2]
    #                             lc += [fig_colors[vv] if vv in list(fig_colors.keys()) else "cyan"]
    #                             lx += [da["year"].values]
    #                             ly1 += [da.min(dim="quantile").values]
    #                             ly2 += [da.max(dim="quantile").values]
    #                             lz += [2 if pro2[:4] == "cmip" else 3]
    #                 plot_type = "sha"
    #                 if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
    #                     put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(plot_type) + "_y1",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(plot_type) + "_y2",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- text
    #                 lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # legend
    #                     # project
    #                     if jj == 0:
    #                         lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         lt += [fig_titles[pro1]]
    #                         lv += ["center"]
    #                         lw += ["regular"]
    #                         lx += [x1 - 20 * dx]
    #                         ly += [y1 + 50 * dy]
    #                     # variable
    #                     if jj == 0 and pro1 == list_pro[0]:
    #                         # plot size
    #                         y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
    #                         lt += [val.replace("UNITS", metadata_json[dia]["units"])]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         lx += [x1 - 30 * dx]
    #                         ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
    #                     # ensemble size
    #                     if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
    #                         k1 = 0
    #                         for pro2 in list(n_ensemble.keys()):
    #                             lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                             lf += [15]
    #                             lh += ["left"]
    #                             lr += [0]
    #                             lt += [fig_titles[pro2]]
    #                             lv += ["top"]
    #                             lw += ["regular"]
    #                             lx += [x2 + 5 * dx]
    #                             ly += [y2 - (k1 * 8 + 10) * dy]
    #                             k1 += 1
    #                             for exp2 in list(n_ensemble[pro2].keys()):
    #                                 if exp2 == "piControl":
    #                                     continue
    #                                 n_ens = n_ensemble[pro2][exp2].split(" ")[0]
    #                                 n_mem = n_ensemble[pro2][exp2].split(" ")[1]
    #                                 lc += ["k"] * 3
    #                                 lf += [15] * 3
    #                                 lh += ["left", "right", "right"]
    #                                 lr += [0] * 3
    #                                 lt += [fig_titles[exp2], n_ens, n_mem]
    #                                 lv += ["top"] * 3
    #                                 lw += ["regular"] * 3
    #                                 lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
    #                                 ly += [y2 - (k1 * 8 + 10) * dy] * 3
    #                                 k1 += 1
    #                     # plot title
    #                     if ii == 0 and jj == 0:
    #                         # plot size
    #                         l1 = len([k1 for k1 in list_exp if k1 == "historical"])
    #                         x1_d = fig_panel_size[tp]["x_delt"]
    #                         x1_m = fig_panel_size[tp]["x_size"]
    #                         x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
    #                         x1_s = max(x1_s) - min(x1_s)
    #                         l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
    #                         x2_d = fig_panel_size["panel_2"]["x_delt"]
    #                         x2_m = fig_panel_size["panel_2"]["x_size"]
    #                         x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
    #                         x2_s = max(x2_s) - min(x2_s)
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [20]
    #                         lh += ["center"]
    #                         lr += [0]
    #                         lt += [fig_titles["figure_title"]]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         a1 = x1_m * x1_s
    #                         a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
    #                         lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
    #                         ly += [y2 + 20 * dy]
    #                 plot_type = "text"
    #                 if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
    #                     put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
    #                     put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
    #                     put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # #
    # # -- Plot
    # #
    # # output plot directory (relative to current file directory)
    # figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # # figure name
    # a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    # a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # # figure loop
    # for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
    #     # list plot names
    #     list_names = list(plot_data[dur].keys())
    #     # path to output figure
    #     figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "03")))
    #     # plot
    #     fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
    #               figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
    # figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # # panel 1: time series
    # print("time series 4")
    # group, panel = "group_1", "panel_1"
    # for dur in list(change_m.keys()):
    #     if not check_dict_levels(change_r, dur):
    #         continue
    #     ny1 = int(dur.split("_")[0])
    #     ny2 = int(data_epoch_lengths[-1].split("_")[0])
    #     list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
    #     list_dia = [k for k in data_diagnostics1 if k in list_dia]
    #     list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
    #     list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
    #     ii = -1
    #     for dia in list_dia:
    #         for pro1 in list_pro:
    #             ii += 1
    #             list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
    #                               key=lambda s: s.lower())
    #             for jj, exp in enumerate(list_exp):
    #                 if not check_dict_levels(change_m, dur, dia, pro1, exp):
    #                     continue
    #                 tp = "panel_1" if exp == "historical" else "panel_2"
    #                 pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
    #                 grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
    #                 tuple_k = (dur, grp, pan)
    #                 print(dur, dia, pro1, exp.ljust(10), tuple_k)
    #                 # panel size
    #                 tk = tuple_k[:-2]
    #                 if not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                     for k1, d1 in fig_panel_size.items():
    #                         if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                             for k2, d2 in d1.items():
    #                                 if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
    #                                     continue
    #                                 tmp = copy__deepcopy(d2)
    #                                 tt = (dur, proc_ci, "time", tp, "lim")
    #                                 if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
    #                                     k3 = read_from_dict(fig_ticks, *tt)
    #                                     tmp *= (k3[-1] - k3[0])
    #                                 put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
    #                         elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
    #                             put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
    #                 # panel axes definition
    #                 l1, l2 = ["x_", "y_"], ["time", dia]
    #                 for vv, dd in zip(l1, l2):
    #                     if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
    #                         for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
    #                             k3 = str(vv) + str(k1)
    #                             if check_dict_levels(figure_axes, *tuple_k + (k3,)):
    #                                 continue
    #                             tmp = copy__deepcopy(k2)
    #                             if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
    #                                 yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
    #                                 tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
    #                                 tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
    #                                 tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
    #                                 tmp = tmp.replace(" ()", "")
    #                             if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
    #                                                               not (dia == list_dia[-1] and
    #                                                                    (pro1 == list_pro[-1] or exp == list_exp[-1]))):
    #                                 tmp = ""
    #                             put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
    #                 dt = {}
    #                 if check_dict_levels(figure_axes, *tuple_k):
    #                     dt = read_from_dict(figure_axes, *tuple_k)
    #                 # test possibility to plot text / legend
    #                 do_legend = False
    #                 if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
    #                         dt["y_lim"] is not None:
    #                     do_legend = True
    #                 # -- title
    #                 if ii == 0:
    #                     put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
    #                 # -- curve
    #                 lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_m, dur, dia, pro2, exp):
    #                         da = copy__deepcopy(change_m[dur][dia][pro2][exp])
    #                         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                         ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
    #                         lw += [2 if pro2[:4] == "cmip" else 1]
    #                         lx += [da["year"].values]
    #                         ly += [da.values]
    #                         lz += [4 if pro2[:4] == "cmip" else 5]
    #                 if check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
    #                     tmp = "iv_on_mean--" + str(pro1)
    #                     da = copy__deepcopy(change_r[dur][dia][pro1][exp]["error"])
    #                     lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
    #                     ls += ["--"] * 2
    #                     lw += [1] * 2
    #                     lx += [da["year"].values] * 2
    #                     ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
    #                     lz += [4] * 2
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # x = 0 and y = 0
    #                     for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
    #                         if k1[0] < 0 < k1[1]:
    #                             lc += ["k"]
    #                             ls += [":"]
    #                             lw += [1]
    #                             lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
    #                             ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
    #                             lz += [10]
    #                 plot_type = "cur"
    #                 if len(lx) > 0 and len(lx) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- shading
    #                 la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_r, dur, dia, pro2, exp):
    #                         for k2, (typ, da) in enumerate(change_r[dur][dia][pro2][exp].items()):
    #                             if typ == "error":
    #                                 continue
    #                             vv = str(typ) + "--" + str(pro2)
    #                             la += [0.3 + k2 * 0.2]
    #                             lc += [fig_colors[vv] if vv in list(fig_colors.keys()) else "cyan"]
    #                             lx += [da["year"].values]
    #                             ly1 += [da.min(dim="quantile").values]
    #                             ly2 += [da.max(dim="quantile").values]
    #                             lz += [2 if pro2[:4] == "cmip" else 3]
    #                 plot_type = "sha"
    #                 if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
    #                     put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(plot_type) + "_y1",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(plot_type) + "_y2",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- text
    #                 lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # legend
    #                     # project
    #                     if jj == 0:
    #                         lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         lt += [fig_titles[pro1]]
    #                         lv += ["center"]
    #                         lw += ["regular"]
    #                         lx += [x1 - 20 * dx]
    #                         ly += [y1 + 50 * dy]
    #                     # variable
    #                     if jj == 0 and pro1 == list_pro[0]:
    #                         # plot size
    #                         y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
    #                         lt += [val.replace("UNITS", metadata_json[dia]["units"])]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         lx += [x1 - 30 * dx]
    #                         ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
    #                     # ensemble size
    #                     if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
    #                         k1 = 0
    #                         for pro2 in list(n_ensemble.keys()):
    #                             lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                             lf += [15]
    #                             lh += ["left"]
    #                             lr += [0]
    #                             lt += [fig_titles[pro2]]
    #                             lv += ["top"]
    #                             lw += ["regular"]
    #                             lx += [x2 + 5 * dx]
    #                             ly += [y2 - (k1 * 8 + 10) * dy]
    #                             k1 += 1
    #                             for exp2 in list(n_ensemble[pro2].keys()):
    #                                 if exp2 == "piControl":
    #                                     continue
    #                                 n_ens = n_ensemble[pro2][exp2].split(" ")[0]
    #                                 n_mem = n_ensemble[pro2][exp2].split(" ")[1]
    #                                 lc += ["k"] * 3
    #                                 lf += [15] * 3
    #                                 lh += ["left", "right", "right"]
    #                                 lr += [0] * 3
    #                                 lt += [fig_titles[exp2], n_ens, n_mem]
    #                                 lv += ["top"] * 3
    #                                 lw += ["regular"] * 3
    #                                 lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
    #                                 ly += [y2 - (k1 * 8 + 10) * dy] * 3
    #                                 k1 += 1
    #                     # plot title
    #                     if ii == 0 and jj == 0:
    #                         # plot size
    #                         l1 = len([k1 for k1 in list_exp if k1 == "historical"])
    #                         x1_d = fig_panel_size[tp]["x_delt"]
    #                         x1_m = fig_panel_size[tp]["x_size"]
    #                         x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
    #                         x1_s = max(x1_s) - min(x1_s)
    #                         l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
    #                         x2_d = fig_panel_size["panel_2"]["x_delt"]
    #                         x2_m = fig_panel_size["panel_2"]["x_size"]
    #                         x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
    #                         x2_s = max(x2_s) - min(x2_s)
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [20]
    #                         lh += ["center"]
    #                         lr += [0]
    #                         lt += [fig_titles["figure_title"]]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         a1 = x1_m * x1_s
    #                         a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
    #                         lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
    #                         ly += [y2 + 20 * dy]
    #                 plot_type = "text"
    #                 if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
    #                     put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
    #                     put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
    #                     put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # #
    # # -- Plot
    # #
    # # output plot directory (relative to current file directory)
    # figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # # figure name
    # a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    # a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # # figure loop
    # for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
    #     # list plot names
    #     list_names = list(plot_data[dur].keys())
    #     # path to output figure
    #     figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "04")))
    #     # plot
    #     fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
    #               figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
    # figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # # panel 1: time series
    # print("time series 5")
    # group, panel = "group_1", "panel_1"
    # for dur in list(change_m.keys()):
    #     if not check_dict_levels(change_r, dur):
    #         continue
    #     ny1 = int(dur.split("_")[0])
    #     ny2 = int(data_epoch_lengths[-1].split("_")[0])
    #     list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
    #     list_dia = [k for k in data_diagnostics1 if k in list_dia]
    #     list_pro = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
    #     list_pro = [k for k in data_projects if k in list_pro and k != "observations"]
    #     ii = -1
    #     for dia in list_dia:
    #         for pro1 in list_pro:
    #             ii += 1
    #             list_exp = sorted([k for k in list(change_r[dur][dia][pro1].keys()) if "piControl" not in k],
    #                               key=lambda s: s.lower())
    #             for jj, exp in enumerate(list_exp):
    #                 if not check_dict_levels(change_m, dur, dia, pro1, exp):
    #                     continue
    #                 tp = "panel_1" if exp == "historical" else "panel_2"
    #                 pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
    #                 grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
    #                 tuple_k = (dur, grp, pan)
    #                 print(dur, dia, pro1, exp.ljust(10), tuple_k)
    #                 # panel size
    #                 tk = tuple_k[:-2]
    #                 if not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                     for k1, d1 in fig_panel_size.items():
    #                         if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
    #                             for k2, d2 in d1.items():
    #                                 if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
    #                                     continue
    #                                 tmp = copy__deepcopy(d2)
    #                                 tt = (dur, proc_ci, "time", tp, "lim")
    #                                 if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
    #                                     k3 = read_from_dict(fig_ticks, *tt)
    #                                     tmp *= (k3[-1] - k3[0])
    #                                 put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
    #                         elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
    #                             put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
    #                 # panel axes definition
    #                 l1, l2 = ["x_", "y_"], ["time", dia]
    #                 for vv, dd in zip(l1, l2):
    #                     if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
    #                         for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
    #                             k3 = str(vv) + str(k1)
    #                             if check_dict_levels(figure_axes, *tuple_k + (k3,)):
    #                                 continue
    #                             tmp = copy__deepcopy(k2)
    #                             if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
    #                                 yy = 2006  # 2014 if data_projects[0] == "cmip6" else 2004
    #                                 tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
    #                                 tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
    #                                 tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
    #                                 tmp = tmp.replace(" ()", "")
    #                             if (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and
    #                                                               not (dia == list_dia[-1] and
    #                                                                    (pro1 == list_pro[-1] or exp == list_exp[-1]))):
    #                                 tmp = ""
    #                             put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
    #                 dt = {}
    #                 if check_dict_levels(figure_axes, *tuple_k):
    #                     dt = read_from_dict(figure_axes, *tuple_k)
    #                 # test possibility to plot text / legend
    #                 do_legend = False
    #                 if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
    #                         dt["y_lim"] is not None:
    #                     do_legend = True
    #                 # -- title
    #                 if ii == 0:
    #                     put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
    #                 # -- curve
    #                 lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_m, dur, dia, pro2, exp):
    #                         da = copy__deepcopy(change_m[dur][dia][pro2][exp])
    #                         lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                         ls += [fig_linestyles[pro2] if pro2 in list(fig_linestyles.keys()) else ":"]
    #                         lw += [2 if pro2[:4] == "cmip" else 1]
    #                         lx += [da["year"].values]
    #                         ly += [da.values]
    #                         lz += [4 if pro2[:4] == "cmip" else 5]
    #                 exp_t = "--".join([exp, "piControl"])
    #                 if check_dict_levels(change_r, dur, dia, pro1, exp_t, "full"):
    #                     da = copy__deepcopy(change_r[dur][dia][pro1][exp_t]["full"])
    #                     lc += [fig_colors["piControl"] if "piControl" in list(fig_colors.keys()) else "cyan"] * 2
    #                     ls += [fig_linestyles["piControl"] if "piControl" in list(fig_linestyles.keys()) else ":"] * 2
    #                     lw += [2] * 2
    #                     lx += [da["year"].values] * 2
    #                     ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
    #                     lz += [1] * 2
    #                 if check_dict_levels(change_r, dur, dia, pro1, exp, "error"):
    #                     tmp = "iv_on_mean--" + str(pro1)
    #                     da = copy__deepcopy(change_r[dur][dia][pro1][exp]["error"])
    #                     lc += [fig_colors[tmp] if tmp in list(fig_colors.keys()) else "cyan"] * 2
    #                     ls += ["--"] * 2
    #                     lw += [1] * 2
    #                     lx += [da["year"].values] * 2
    #                     ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
    #                     lz += [4] * 2
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # x = 0 and y = 0
    #                     for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
    #                         if k1[0] < 0 < k1[1]:
    #                             lc += ["k"]
    #                             ls += [":"]
    #                             lw += [1]
    #                             lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
    #                             ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
    #                             lz += [10]
    #                 plot_type = "cur"
    #                 if len(lx) > 0 and len(lx) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- shading
    #                 la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
    #                 for pro2 in [pro1, "observations"]:
    #                     if check_dict_levels(change_r, dur, dia, pro2, exp):
    #                         for k2, (typ, da) in enumerate(change_r[dur][dia][pro2][exp].items()):
    #                             if typ == "error":
    #                                 continue
    #                             vv = str(typ) + "--" + str(pro2)
    #                             la += [0.3 + k2 * 0.2]
    #                             lc += [fig_colors[vv] if vv in list(fig_colors.keys()) else "cyan"]
    #                             lx += [da["year"].values]
    #                             ly1 += [da.min(dim="quantile").values]
    #                             ly2 += [da.max(dim="quantile").values]
    #                             lz += [2 if pro2[:4] == "cmip" else 3]
    #                 plot_type = "sha"
    #                 if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
    #                     put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(plot_type) + "_y1",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(plot_type) + "_y2",))
    #                     put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
    #                 # -- text
    #                 lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
    #                 if isinstance(do_legend, bool) and do_legend:
    #                     df = read_from_dict(figure_panel_size, *tk)
    #                     x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
    #                     dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
    #                     dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
    #                     # legend
    #                     # project
    #                     if jj == 0:
    #                         lc += [fig_colors[pro1] if pro1 in list(fig_colors.keys()) else "cyan"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         lt += [fig_titles[pro1]]
    #                         lv += ["center"]
    #                         lw += ["regular"]
    #                         lx += [x1 - 20 * dx]
    #                         ly += [y1 + 50 * dy]
    #                     # variable
    #                     if jj == 0 and pro1 == list_pro[0]:
    #                         # plot size
    #                         y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [15]
    #                         lh += ["center"]
    #                         lr += [90]
    #                         val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
    #                         lt += [val.replace("UNITS", metadata_json[dia]["units"])]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         lx += [x1 - 30 * dx]
    #                         ly += [y2 - (y1_s * len(list_pro) + y1_d * (len(list_pro) - 1)) * (y2 - y1) / (2 * y1_s)]
    #                     # ensemble size
    #                     if dia == list_dia[0] and pro1 == "cmip5" and exp == list_exp[-1]:
    #                         k1 = 0
    #                         for pro2 in list(n_ensemble.keys()):
    #                             lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
    #                             lf += [15]
    #                             lh += ["left"]
    #                             lr += [0]
    #                             lt += [fig_titles[pro2]]
    #                             lv += ["top"]
    #                             lw += ["regular"]
    #                             lx += [x2 + 5 * dx]
    #                             ly += [y2 - (k1 * 8 + 10) * dy]
    #                             k1 += 1
    #                             for exp2 in list(n_ensemble[pro2].keys()):
    #                                 if exp2 == "piControl":
    #                                     continue
    #                                 n_ens = n_ensemble[pro2][exp2].split(" ")[0]
    #                                 n_mem = n_ensemble[pro2][exp2].split(" ")[1]
    #                                 lc += ["k"] * 3
    #                                 lf += [15] * 3
    #                                 lh += ["left", "right", "right"]
    #                                 lr += [0] * 3
    #                                 lt += [fig_titles[exp2], n_ens, n_mem]
    #                                 lv += ["top"] * 3
    #                                 lw += ["regular"] * 3
    #                                 lx += [x2 + 5 * dx, x2 + 50 * dx, x2 + 70 * dx]
    #                                 ly += [y2 - (k1 * 8 + 10) * dy] * 3
    #                                 k1 += 1
    #                     # plot title
    #                     if ii == 0 and jj == 0:
    #                         # plot size
    #                         l1 = len([k1 for k1 in list_exp if k1 == "historical"])
    #                         x1_d = fig_panel_size[tp]["x_delt"]
    #                         x1_m = fig_panel_size[tp]["x_size"]
    #                         x1_s = read_from_dict(fig_ticks, dur, proc_ci, "time", tp, "lim")
    #                         x1_s = max(x1_s) - min(x1_s)
    #                         l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "sce", "ssp"]])
    #                         x2_d = fig_panel_size["panel_2"]["x_delt"]
    #                         x2_m = fig_panel_size["panel_2"]["x_size"]
    #                         x2_s = read_from_dict(fig_ticks, dur, proc_ci, "time", "panel_2", "lim")
    #                         x2_s = max(x2_s) - min(x2_s)
    #                         # add to text
    #                         lc += ["k"]
    #                         lf += [20]
    #                         lh += ["center"]
    #                         lr += [0]
    #                         lt += [fig_titles["figure_title"]]
    #                         lv += ["center"]
    #                         lw += ["bold"]
    #                         a1 = x1_m * x1_s
    #                         a2 = (a1 + x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
    #                         lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
    #                         ly += [y2 + 20 * dy]
    #                 plot_type = "text"
    #                 if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
    #                     put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
    #                     put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
    #                     put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
    #                     put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
    #                     put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
    #                     put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
    #                     put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
    #                     put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
    #                     put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    # #
    # # -- Plot
    # #
    # # output plot directory (relative to current file directory)
    # figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # # figure name
    # a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    # a3, a4 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci), str(len(smile_rejected)).zfill(2)
    # figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
    #     ["", str(a1) + "mem", str(a2) + "po", str(a3) + "p", str(a4) + "rem"])
    # # figure loop
    # for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
    #     # list plot names
    #     list_names = list(plot_data[dur].keys())
    #     # path to output figure
    #     figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, "05")))
    #     # plot
    #     fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
    #               figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#
