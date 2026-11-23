# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure 1: evolution of variable change in different experiments
#    * time (x) vs. signal (y)
#    * time (x) vs. noise (y)
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
from enso_change.wrapper.xarray_base import concatenate, polyval
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
    "data_diagnostics1": ["ave_ts_val_glob"],
    # "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
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
    "proc_ci": 90.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {
        # no resample 8000
        "ave_ts_val_glob": {"hs09--1": 62.58, "hs09--2": 6.93, "hs09--3": 3.69,
                            "hs09--4": 1.80, "hs09--5": 0.92, "hs09--6": 0.56},
        "ske_pr_ano_n30e": {"hs09--1": 0.77, "hs09--2": 0.63, "hs09--3": 0.54},
        "ske_pr_ano_n34e": {"hs09--1": 0.78, "hs09--2": 0.64, "hs09--3": 0.55},
        "ske_pr_ano_n40e": {"hs09--1": 0.76, "hs09--2": 0.63, "hs09--3": 0.53},
        "ske_ts_ano_n30e": {"hs09--1": 0.73, "hs09--2": 0.62, "hs09--3": 0.53},
        "ske_ts_ano_n34e": {"hs09--1": 0.74, "hs09--2": 0.62, "hs09--3": 0.53},
        "ske_ts_ano_n40e": {"hs09--1": 0.75, "hs09--2": 0.63, "hs09--3": 0.54},
        "std_pr_ano_n30e": {"hs09--1": 1.01, "hs09--2": 0.64, "hs09--3": 0.51},
        "std_pr_ano_n34e": {"hs09--1": 1.07, "hs09--2": 0.66, "hs09--3": 0.53},
        "std_pr_ano_n40e": {"hs09--1": 0.93, "hs09--2": 0.63, "hs09--3": 0.51},
        "std_ts_ano_n30e": {"hs09--1": 0.75, "hs09--2": 0.61, "hs09--3": 0.49},
        "std_ts_ano_n34e": {"hs09--1": 0.76, "hs09--2": 0.61, "hs09--3": 0.49},
        "std_ts_ano_n40e": {"hs09--1": 0.76, "hs09--2": 0.61, "hs09--3": 0.50},
        # resample 8000
        # "std_pr_ano_n34e": {"hs09--1": 1.11, "hs09--2": 0.68, "hs09--3": 54},
        # "std_ts_ano_n34e": {"hs09--1": 78, "hs09--2": 0.63, "hs09--3": 51}
    },
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 1000,
    # threshold for classification
    "proc_th": {"ave": 5, "std": 2, "ske": 2},
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
        "020": "dodgerblue",
        "030": "green",
        "040": "orange",
        "050": "red",
        "max": "k",
        "sha": "k",
        "i": "k",
        "f": "k",
        "m": "k",
        "t": "k",
    },
    # line styles
    "fig_linestyles": {
        "020": "-",
        "030": "-",
        "040": "-",
        "050": "-",
        "max": "-",
        "sha": "-",
        "i": "--",
        "f": "-",
        "m": "-",
        "t": "-",
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            60: {
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 60],
                        "nam": "noise (%)",
                        "tic": list(range(0, 61, 15)),
                    },
                    "panel_2": {
                        "lim": [0, 60],
                        "nam": "",
                        "tic": list(range(0, 61, 15)),
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
                "variable": {
                    "panel_1": {
                        "lim": None,
                        "nam": "(UNITS)",
                        "tic": None,
                    },
                    "panel_2": {
                        "lim": None,
                        "nam": "",
                        "tic": None,
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
        # experiments
        "historical": "Historical",
        "scenario26": "SSP1-2.6",
        "scenario45": "SSP2-4.5",
        "scenario70": "SSP3-7.0",
        "scenario85": "SSP5-8.5",
        "ssp126": "SSP1-2.6",
        "ssp245": "SSP2-4.5",
        "ssp370": "SSP3-7.0",
        "ssp585": "SSP5-8.5",
        # statistic
        "i": "RIV",
        "f": "spread",
        "m": "mean",
        "t": "spread",
        "mea": "mean",
        "spr": "full spread",
        "spt": "full spread minus RIV",
        "sha": "95% CI",
        # variables
        "ave_ts_val": "$^\\circ$C",
        "ave_ts_val_glob": "change of GMST",
        "std_ts_ano": "%",
        "std_ts_ano_n34e": "relative change of N3.4 SSTA $\\sigma$",
        # figure title
        "figure_title": "Inter-model mean & spread (THR%) of VALUE",
        # "figure_title": "Inter-model spread (THR%) of relative change of SSTA $\\sigma$",
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 30,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e78_noise_vs_models_plot(
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
                        if dia == "ave_ts_val_glob" and dat in ["ACCESS-CM2", "GISS-E2-1-G_p3f1"]:
                            print(dia, dur, pro, exp.ljust(10), dat.ljust(17))
                            print("da_h", da_h["year"].values)
                            print("da_s", da_s["year"].values)
                        # concatenate along year axis
                        da_o = concatenate((da_h, da_s), "year")
                        if dia == "ave_ts_val_glob" and dat in ["ACCESS-CM2", "GISS-E2-1-G_p3f1"]:
                            print(dia, dur, pro, exp.ljust(10), dat.ljust(17))
                            print("da_o", da_o["year"].values)
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
    stop
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
                for exp1 in sorted(list(values_json[dia][dur1][pro].keys()), key=lambda s: s.lower()):
                    for dat in sorted(list(values_json[dia][dur1][pro][exp1].keys()), key=lambda s: s.lower()):
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
                        yy = 2000
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # save
                        exp_o = copy__deepcopy(exp1)
                        # exp_o = copy__deepcopy(exp1) if exp1 in ["piControl", "historical", "ssp370"] else \
                        #     "scenario" + str(exp1[-2:])
                        if dia == "ave_ts_val_glob" and dat in ["ACCESS-CM2", "GISS-E2-1-G_p3f1"]:
                            print(dur1, dia, pro, exp_o.ljust(10), dat.ljust(17), da1["year"].values)
                            # a1 = da1.sel(indexers={"year": slice(2075, 2080)}).mean(skipna=True)
                            # a2 = da1.sel(indexers={"year": slice(2081, 2100)}).mean(skipna=True)
                            # # if float(a2.values) < float(a1.values):
                            # #     print(dur1, dia, pro, exp_o.ljust(10), dat.ljust(17), float(a1.values), float(a2.values))
                            # print(dur1, dia, pro, exp_o.ljust(10), dat.ljust(17),
                            #       "{0:.2f}".format(round(float(a1.values), 2)).rjust(8),
                            #       "{0:.2f}".format(round(float(a2.values), 2)).rjust(8))
                            # print(da1.sel(indexers={"year": slice(2075, 2100)}).mean(dim="member").values)
                        put_in_dict(change_em, da1.compute(), dur1, dia, pro, exp_o, dat)
    print("change_em", list(change_em.keys()))
    stop
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
                        put_in_dict(residual_iv, copy__deepcopy(riv).compute(), dur1, dia, pro, ex1, dat)
    print("residual_iv", list(residual_iv.keys()))
    # -- compute multimodel mean and std
    # mean
    mean_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # mmm and range
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_signal, change_noise = {}, {}
    for dur in list(mean_change.keys()):
        if dur == dur2:
            continue
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    if pro == "observations" or exp == "piControl" or \
                            not check_dict_levels(change_em, dur, dia, pro, exp) or \
                            not check_dict_levels(residual_iv, dur, dia, pro, exp):
                        continue
                    # dictionaries
                    d1 = mean_change[dur][dia][pro][exp]
                    d3 = residual_iv[dur][dia][pro][exp]
                    # list datasets
                    list_datasets = list(set(list(d1.keys())) & set(list(d3.keys())))
                    # arrays
                    da1 = concatenate(tuple(d1[dat] for dat in list_datasets), "dataset")
                    da3 = concatenate(tuple(d3[dat] for dat in list_datasets), "dataset")
                    size = da1["dataset"].size
                    list_nbr = sorted(list(set(list(range(20, min(size, 50) + 1, 10)) + [size])))
                    # resample
                    for nbr in list_nbr:
                        if nbr == size:
                            da12, da32 = copy__deepcopy(da1), copy__deepcopy(da3)
                        else:
                            da12, da32 = resample(
                                da1, arr_i2=da3, dim="dataset", pseudo_sample_size=proc_ps, sample_method="combination",
                                sample_size=nbr)
                        # -- scaled distribution
                        # initial multimodel mean using ensemble mean -> will not change
                        mi = da12.mean(dim="dataset", skipna=True)
                        # initial multimodel spread using ensemble mean -> that what will be scaled
                        si = da12.var(dim="dataset", skipna=True)
                        # initial multimodel mean residual internal variability -> will be used for scaling
                        ri = da32.mean(dim="dataset", skipna=True)
                        # final spread (initial spread minus residual internal variability) -> must be set to 0 if
                        # negative
                        sf = (si - ri) / si
                        sf = sf.where(sf > 0, other=0)
                        # scaled distribution
                        da4 = mi + (da12 - mi) * sf
                        # mean
                        nbr_o = str(nbr).zfill(3) if nbr != size else "max"
                        mi = mi.rolling(**{"year": 5}, min_periods=1, center=True).mean()
                        put_in_dict(change_noise, mi.compute(),dur, dia, pro, "mea", exp, nbr_o, "m")
                        # full spread
                        nbr_o = str(nbr).zfill(3) if nbr != size else "max"
                        da12 = statistic(da12, dim="dataset", stat="per", per=percentiles, skipna=True)
                        da12 = da12.max(dim="quantile") - da12.min(dim="quantile")
                        da12 = da12.rolling(**{"year": 5}, min_periods=1, center=True).mean()
                        put_in_dict(change_noise, da12.compute(),dur, dia, pro, "spr", exp, nbr_o, "f")
                        # internal variability
                        put_in_dict(change_noise, ri.compute()**0.5,dur, dia, pro, "spr", exp, nbr_o, "i")
                        # true model spread
                        da4 = statistic(da4, dim="dataset", stat="per", per=percentiles, skipna=True)
                        da4 = da4.max(dim="quantile") - da4.min(dim="quantile")
                        da4 = da4.rolling(**{"year": 5}, min_periods=1, center=True).mean()
                        put_in_dict(change_noise, da4.compute(),dur, dia, pro, "spt", exp, nbr_o, "t")
    # min-max
    print("min-max")
    for dur in sorted(list(change_noise.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(change_noise[dur].keys()), key=lambda s: s.lower()):
            dt = {}
            for pro in list(change_noise[dur][dia].keys()):
                for typ in list(change_noise[dur][dia][pro].keys()):
                    for exp in list(change_noise[dur][dia][pro][typ].keys()):
                        for nbr in list(change_noise[dur][dia][pro][typ][exp].keys()):
                            for da in change_noise[dur][dia][pro][typ][exp][nbr].values():
                                put_in_dict(dt, da.values.flatten().tolist(), typ)
            for typ in sorted(list(dt.keys()), key=lambda s: s.lower()):
                if len(dt[typ]) > 0:
                    print(dur, dia, typ, "{0:.2f}".format(round(min(dt[typ]), 2)).rjust(8),
                          "{0:.2f}".format(round(max(dt[typ]), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    ranges = {
        "std_ts_ano_n34e": {
            "mea": list(range(-26, 27, 13)),
            "spr": list(range(0, 101, 25)),
            "spt": list(range(0, 81, 20)),
        },
        "ave_ts_val_glob": {
            "mea": [round(k / 10, 1) for k in range(-17, 52, 17)],
            "spr": [round(k / 10, 1) for k in range(0, 37, 9)],
            "spt": [round(k / 10, 1) for k in range(0, 37, 9)],
        },
    }
    stop
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in list(change_noise.keys()):
        for dia in list(change_noise[dur].keys()):
            for pro in list(change_noise[dur][dia].keys()):
                list_typ = list(change_noise[dur][dia][pro].keys())
                for ii, typ in enumerate(list_typ):
                    list_exp = sorted(list(change_noise[dur][dia][pro][typ].keys()), key=lambda s: s.lower())
                    for jj, exp in enumerate(list_exp):
                        tp = "panel_1" if exp == "historical" else "panel_2"
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = ("_".join([dur, dia, pro]), grp, pan)
                        print(dur, dia, pro, typ, exp.ljust(10), tuple_k)
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
                        l1, l2 = ["x_", "y_"], ["time", "variable"]
                        dia_t = "_".join(dia.split("_")[:-1])
                        for vv, dd in zip(l1, l2):
                            if check_dict_levels(fig_ticks, dur, proc_ci, dd, tp):
                                for k1, k2 in fig_ticks[dur][proc_ci][dd][tp].items():
                                    k3 = str(vv) + str(k1)
                                    if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                        continue
                                    tmp = copy__deepcopy(k2)
                                    if k3 in ["y_lim", "y_tic"] and tmp is None and \
                                            check_dict_levels(ranges, dia, typ):
                                        tmp = ranges[dia][typ]
                                        if k3 == "y_lim":
                                            tmp = [min(tmp), max(tmp)]
                                    elif (k3 == "y_nam" and jj > 0) or (k3 == "x_nam" and typ != list_typ[-1]):
                                        tmp = ""
                                    elif k3 == "y_nam" and "UNITS" in tmp and dia_t in list(fig_titles.keys()):
                                        tmp = tmp.replace("UNITS", fig_titles[dia_t])
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
                            val = str(fig_titles[exp]) + " (" + str(n_ensemble[pro][exp].split(" (")[0]) + ")"
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0:
                            put_in_dict(plot_data, fig_titles[typ], *tuple_k + ("title_row",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        leg_co, leg_ls = [], []
                        if check_dict_levels(change_noise, dur, dia, pro, typ, exp):
                            for nbr in list(change_noise[dur][dia][pro][typ][exp].keys()):
                                if nbr != "max":
                                    continue
                                for nn, da in change_noise[dur][dia][pro][typ][exp][nbr].items():
                                    da1 = da.mean(dim="resample") if "resample" in list(da.dims) else copy__deepcopy(da)
                                    lc += [fig_colors[nbr] if nbr in list(fig_colors.keys()) else "cyan"]
                                    ls += [fig_linestyles[nn] if nn in list(fig_linestyles.keys()) else ":"]
                                    lw += [2]
                                    lx += [[1850, 2100] if nn == "i" else da1["year"].values]
                                    ly += [[float(da1.values)] * 2 if nn == "i" else da1.values]
                                    lz += [11]
                                    leg_co += [nbr]
                                    leg_ls += [nn]
                        leg_co = sorted(list(set(leg_co)), key=lambda s: s.lower())
                        leg_ls = sorted(list(set(leg_ls)), key=lambda s: s.lower())
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
                                    lz += [12]
                            # legend
                            if jj == 0:
                                for k1, k2 in enumerate(leg_ls):
                                    lc += [fig_colors[k2] if k2 in list(fig_colors.keys()) else "cyan"]
                                    ls += [fig_linestyles[k2] if k2 in list(fig_linestyles.keys()) else ":"]
                                    lw += [2]
                                    lx += [[x1 + k3 * dx for k3 in [2, 15]]]
                                    ly += [[y2 - (k1 * 8 + 12) * dy] * 2]
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
                        if check_dict_levels(change_noise, dur, dia, pro, typ, exp):
                            for nbr in list(change_noise[dur][dia][pro][typ][exp].keys()):
                                for nn, da in change_noise[dur][dia][pro][typ][exp][nbr].items():
                                    if "resample" not in list(da.dims):
                                        continue
                                    da1 = statistic(da, dim="resample", stat="per", per=[2.5, 97.5], skipna=True)
                                    la += [0.8]
                                    lc += [fig_colors[nbr] if nbr in list(fig_colors.keys()) else "cyan"]
                                    lx += [[1850, 2100] if nn == "i" else da1["year"].values]
                                    ly1 += [da1.min(dim="quantile").values]
                                    ly2 += [da1.max(dim="quantile").values]
                                    lz += [2]
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            if jj == 0:
                                la += [0.5]
                                lc += ["k"]
                                lx += [[x1 + k3 * dx for k3 in [2, 15]]]
                                ly1 += [[y2 - (len(leg_ls) * 8 + 12 - 3) * dy] * 2]
                                ly2 += [[y2 - (len(leg_ls) * 8 + 12 + 3) * dy] * 2]
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
                            if jj == 0:
                                for k1, k2 in enumerate(leg_co):
                                    lc += [fig_colors[k2] if k2 in list(fig_colors.keys()) else "cyan"]
                                    lf += [15]
                                    lh += ["right"]
                                    lr += [0]
                                    lt += [fig_titles[k2] if k2 in list(fig_titles.keys()) else copy__deepcopy(k2)]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x2 - 2 * dx]
                                    ly += [y2 - (k1 * 8 + 5) * dy]
                                for k1, k2 in enumerate(leg_ls + ["sha"]):
                                    lc += [fig_colors[k2] if k2 in list(fig_colors.keys()) else "cyan"]
                                    lf += [15]
                                    lh += ["left"]
                                    lr += [0]
                                    lt += [fig_titles[k2] if k2 in list(fig_titles.keys()) else copy__deepcopy(k2)]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x1 + 17 * dx]
                                    ly += [y2 - (k1 * 8 + 12) * dy]
                            # variable
                            # if jj == 0 and typ == list_typ[0]:
                            #     # plot size
                            #     y1_d, y1_s = fig_panel_size[tp]["y_delt"], fig_panel_size[tp]["y_size"]
                            #     # add to text
                            #     lc += ["k"]
                            #     lf += [15]
                            #     lh += ["center"]
                            #     lr += [90]
                            #     lt += [fig_titles[dia]]
                            #     lv += ["center"]
                            #     lw += ["bold"]
                            #     tmp = panel_param["title_row_x"] if "title_row_x" in list(panel_param.keys()) else 30
                            #     lx += [x1 - tmp * dx]
                            #     ly += [y2 - (y1_s * len(list_typ) +
                            #                  y1_d * (len(list_typ) - 1)) * (y2 - y1) / (2 * y1_s)]
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
                                lt += [fig_titles["figure_title"].replace("THR", "-".join(
                                    [str(int(k4)) if k4 == int(k4) else str(k4) for k4 in percentiles])).replace(
                                    "VALUE", fig_titles[dia])]
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
    a3 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a3) + "p"])
    # figure loop
    for kk in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[kk].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", kk)))
        # plot
        fig_basic(plot_data[kk], list_names, 1, figure_axes[kk], fig_format, figure_name_and_path,
                  figure_panel_size[kk], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e78_noise_vs_models_plot()
