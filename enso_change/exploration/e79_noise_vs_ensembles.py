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

# numpy
from numpy import concatenate as numpy__concatenate
from numpy import full as numpy__full
from numpy import meshgrid as numpy__meshgrid
from numpy import nan as numpy__nan

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import array_wrapper, concatenate, polyval
from enso_change.wrapper.xarray_tools import numpy__array, resample, statistic
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
    # "data_experiments": ["historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical", "ssp245", "ssp370", "ssp585"],
    # "data_experiments": ["historical", "ssp245"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 5,
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
    # confidence interval as percentile: float
    "proc_cp": 90,
    # method: {"absolute", "relative"}
    "proc_md": "absolute",
    # number of pseudo samples: int
    "proc_ps": 500,
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
        "frac": {"x": 0.1, "y": 0.1},
        "panel_1": {"x_delt": 10, "x_size": 5, "y_delt": 6, "y_size": 5},
    },
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": {},
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            90: {
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "x_lim": None,
                        "x_nam": "nbr of datasets",
                        "x_nbr_minor": 0,
                        "x_tic": None,
                        "y_lim": None,
                        "y_nam": "members per dataset",
                        "y_nbr_minor": 0,
                        "y_tic": None,
                        "z_lim": None,
                        "z_nam": "",
                        "z_tic": None,
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
        "ssp": "Scenarios",
        "ssp126": "SSP1-2.6",
        "ssp245": "SSP2-4.5",
        "ssp370": "SSP3-7.0",
        "ssp585": "SSP5-8.5",
        # statistic
        "i": "riv",
        "f": "spread",
        "t": "spread",
        # "mea": "multi-model mean",
        # "spr": "multi-model spread",
        # "spt": "multi-model spread - RIV",
        "mea": "mean",
        "spr": "full spread",
        "spt": "full spread minus RIV",
        "true": "full spread - RIV",
        "raw": "full spread",
        "sha": "95% CI",
        # variables
        "std_ts_ano": "relative change of SSTA $\\sigma$",
        # figure title
        "absolute": "Absolute error on Inter-model mean & spread (THR%)",
        "relative": "Relative error on Inter-model mean & spread (THR%)",
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
def e79_noise_vs_ensembles_plot(
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
        proc_cp: float = None,
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
    proc_cp = none_to_default(proc_cp, default["proc_cp"])
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
    print("e79_noise_vs_ensembles_plot input")
    print("SMILE size   ", smile_minimum_size)
    print("ens spread CI", proc_ci)
    print("error CI     ", proc_cp)
    print("error method ", proc_md)
    print("error samples", proc_ps)
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
    change = {}
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
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 2000})
                        # resample members
                        size = da1["member"].size
                        # list_nm = sorted(list(range(1, 6, 1)) + list(range(6, 10, 2)) + [10, 15, 20])
                        # list_nm = sorted(list(range(1, 6, 1)) + [10, 15, 20])
                        list_nm = list(range(1, 10, 1))
                        if smile_minimum_size > 1:
                            list_nm = [k for k in list_nm if k <= smile_minimum_size] + ["max"]
                        else:
                            list_nm = [k for k in list_nm if k <= size] + ["max"]
                        for nm in list_nm:
                            nbr = copy__deepcopy(size) if nm == "max" else copy__deepcopy(nm)
                            da12, da22 = resample(da1, arr_i2=da2, dim="member", pseudo_sample_size=proc_ps,
                                                  sample_method="combination", sample_size=nbr)
                            da22 = da22.mean(dim="member")
                            # change
                            da12 -= da22
                            if dia.split("_")[0] in ["std", "var"]:
                                da12 *= 100 / da22
                            # rename resample dimension
                            da12 = da12.rename(new_name_or_name_dict={"resample": "re_mem"})
                            dims = list(set(list(da12.coords)) - set(list(da12.dims)))
                            if len(dims) > 0:
                                da12 = da12.drop_vars(dims)
                            # save
                            nm_o = "m" + (nm if nm == "max" else str(nm).zfill(3))
                            put_in_dict(change, copy__deepcopy(da12.compute()),dur1, dia, pro, exp1, dat, nm_o)
    print("change", list(change.keys()))
    # -- compute residual internal variability
    # residual variability adjusted with corrected IV for each model
    residual_iv = {}
    for dur1 in list(change.keys()):
        if dur1 == dur2:
            continue
        for dia in list(change[dur1].keys()):
            for pro in list(change[dur1][dia].keys()):
                if pro == "observations":
                    continue
                for ex1 in list(change[dur1][dia][pro].keys()):
                    if ex1 == "piControl":
                        continue
                    for dat in list(change[dur1][dia][pro][ex1].keys()):
                        for nme in list(change[dur1][dia][pro][ex1][dat].keys()):
                            # array
                            da_em = copy__deepcopy(change[dur1][dia][pro][ex1][dat][nme])
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
                            riv *= (1 / proc_co[dia]["hs09--" + str(nn)]) * size ** -0.2
                            # save
                            put_in_dict(residual_iv, copy__deepcopy(riv.compute()),dur1, dia, pro, ex1, dat, nme)
    print("residual_iv", list(residual_iv.keys()))
    # -- compute multimodel mean and std
    # mean
    mean_change, _, _ = compute_1var(change, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # mmm and range
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_sta = {}
    for dur in list(mean_change.keys()):
        ny = int(dur.split("_")[0])
        if dur == dur2:
            continue
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    dict_t = mean_change[dur][dia][pro][exp]
                    list_nme = sorted(list(set([nme for dat in list(dict_t.keys())
                                                for nme in list(dict_t[dat].keys())])), key=lambda s: s.lower())
                    dict_nme = {}
                    for nme in list_nme:
                        for dat in list(dict_t.keys()):
                            tmp = sorted(list(dict_t[dat].keys()), key=lambda s: s.lower())
                            tmp = copy__deepcopy(nme) if nme in tmp else tmp[-1]
                            put_in_dict(dict_nme, tmp, nme, dat)
                    for nme in list_nme:
                        if pro == "observations" or exp == "piControl" or \
                                not check_dict_levels(change, dur, dia, pro, exp) or \
                                not check_dict_levels(residual_iv, dur, dia, pro, exp):
                            continue
                        # dictionaries
                        d1, d3 = mean_change[dur][dia][pro][exp], residual_iv[dur][dia][pro][exp]
                        # list datasets
                        list_datasets = list(set(list(d1.keys())) & set(list(d3.keys())))
                        # arrays: concatenate datasets using the given number of members or the maximum number of
                        # members
                        da1 = concatenate(tuple(d1[dat][dict_nme[nme][dat]] for dat in list_datasets), "dataset")
                        da3 = concatenate(tuple(d3[dat][dict_nme[nme][dat]] for dat in list_datasets), "dataset")
                        y1, y2 = (1850 if exp == "historical" else 2000) + ny, 1950 if exp == "historical" else 2100
                        da1 = da1.sel(indexers={"year": slice(y1, y2)})
                        # da1, da3 = da1.mean(dim="re_mem"), da3.mean(dim="re_mem")
                        # print(dur, dia, pro, exp, nme, da1.shape, da1.dims)
                        size = da1["dataset"].size
                        if smile_minimum_size > 1:
                            list_nda = list(range(15, 20, 1)) + list(range(20, size, 5)) + ["max"]
                            # list_nda = list(range(10, min(size, 21), 5)) + ["max"]
                        else:
                            # list_nda = list(range(20, size, 5)) + ["max"]
                            list_nda = list(range(20, min(size, 51), 5)) + ["max"]
                        print(dur, dia, pro, exp.ljust(10), nme, da1.shape, list_nda)
                        # resample
                        for nda in list_nda:
                            nbr = copy__deepcopy(size) if nda == "max" else copy__deepcopy(nda)
                            da12, da32 = resample(da1, arr_i2=da3, dim="dataset", pseudo_sample_size=proc_ps,
                                                  sample_method="combination", sample_size=nbr)
                            # rename resample dimension
                            da12 = da12.rename(new_name_or_name_dict={"resample": "re_dat"})
                            dims = list(set(list(da12.coords)) - set(list(da12.dims)))
                            if len(dims) > 0:
                                da12 = da12.drop_vars(dims)
                            da32 = da32.rename(new_name_or_name_dict={"resample": "re_dat"})
                            dims = list(set(list(da32.coords)) - set(list(da32.dims)))
                            if len(dims) > 0:
                                da32 = da32.drop_vars(dims)
                            # stack datasets and members resamples
                            da12 = da12.stack(dim={"resample": ("re_dat", "re_mem")})
                            dims = list(set(list(da12.coords)) - set(list(da12.dims)))
                            if len(dims) > 0:
                                da12 = da12.drop_vars(dims)
                            da32 = da32.stack(dim={"resample": ("re_dat", "re_mem")})
                            dims = list(set(list(da32.coords)) - set(list(da32.dims)))
                            if len(dims) > 0:
                                da32 = da32.drop_vars(dims)
                            # -- scaled distribution
                            # initial multimodel mean using ensemble mean -> will not change
                            mi = da12.mean(dim="dataset", skipna=True)
                            # initial multimodel spread using ensemble mean -> that what will be scaled
                            si = da12.var(dim="dataset", skipna=True)
                            # initial multimodel mean residual internal variability -> will be used for scaling
                            ri = da32.mean(dim="dataset", skipna=True)
                            # final spread (initial spread minus residual internal variability) -> must be set to 0
                            # if negative
                            sf = (si - ri) / si
                            sf = sf.where(sf > 0, other=0)
                            # scaled distribution
                            da_o = mi + (da12 - mi) * sf
                            # -- mean
                            nm_o = "d" + (str(nda) + str(size).zfill(3) if nda == "max" else str(nda).zfill(3))
                            mi = mi.mean(dim="year")
                            lk = (dur, dia, pro, "mea", exp, nme, nm_o)
                            put_in_dict(change_sta, copy__deepcopy(mi.compute()), *lk)
                            # -- true model spread
                            da12 = statistic(da12, dim="dataset", stat="per", per=percentiles, skipna=True)
                            da12 = da12.max(dim="quantile") - da12.min(dim="quantile")
                            da12 = da12.mean(dim="year")
                            lk = (dur, dia, pro, "spr", exp, nme, nm_o)
                            put_in_dict(change_sta, copy__deepcopy(da12.compute()), *lk)
                            # -- true model spread
                            da_o = statistic(da_o, dim="dataset", stat="per", per=percentiles, skipna=True)
                            da_o = da_o.max(dim="quantile") - da_o.min(dim="quantile")
                            da_o = da_o.mean(dim="year")
                            lk = (dur, dia, pro, "spt", exp, nme, nm_o)
                            put_in_dict(change_sta, copy__deepcopy(da_o.compute()), *lk)
    print("change_sta", list(change_sta.keys()))
    # -- error
    to_plot = {}
    for dur in list(change_sta.keys()):
        for dia in list(change_sta[dur].keys()):
            for pro in list(change_sta[dur][dia].keys()):
                for sta in list(change_sta[dur][dia][pro].keys()):
                    dt = {}
                    for exp in list(change_sta[dur][dia][pro][sta].keys()):
                        dict_t = change_sta[dur][dia][pro][sta][exp]
                        # matrix of resampled data
                        da = ()
                        list_nme = sorted(list(dict_t.keys()), key=lambda s: s.lower())
                        for nme in list_nme:
                            list_nda = sorted(list(dict_t[nme].keys()), key=lambda s: s.lower())
                            da_t = concatenate(tuple(dict_t[nme][nda] for nda in list_nda), "dataset")
                            list_nda = numpy__array([int(k[1:]) if k[1:].isdigit() else int(1e6) for k in list_nda])
                            da += (da_t.assign_coords({"dataset": list_nda}),)
                        da = concatenate(da, "member")
                        list_nme = numpy__array([int(k[1:]) if k[1:].isdigit() else int(1e6) for k in list_nme])
                        da = da.assign_coords({"member": list_nme})
                        # matrix using maximum available data
                        nme = sorted(list(dict_t.keys()), key=lambda s: s.lower())[-1]
                        nda = sorted(list(dict_t[nme].keys()), key=lambda s: s.lower())[-1]
                        da_r = float(dict_t[nme][nda].mean().values)
                        # relative error
                        da = abs(da - da_r)
                        if proc_md == "relative":
                            da *= 100 / abs(da_r)
                        # percentiles of error
                        da = statistic(da, dim="resample", stat="per", per=proc_cp, skipna=True)
                        # da = da.compute()
                        # # reshape array
                        # list_mem, list_dat = da["member"].values.tolist(), da["dataset"].values.tolist()
                        # list_tot = list(range(20, 50, 5)) + list(range(50, 100, 10)) + list(range(100, 131, 15)) + \
                        #            list(range(150, 201, 25)) + [int(1e6)]
                        # da_o = array_wrapper(
                        #     numpy__full((len(list_tot), len(list_dat)), numpy__nan),
                        #     dims=["member", "dataset"],
                        #     coords={"member": numpy__array(list_tot), "dataset": numpy__array(list_dat)})
                        # for k1 in list_mem:
                        #     for k2 in list_dat:
                        #         nn = k1 * k2
                        #         if 20 <= nn <= 250 or k1 > 1e3 or k2 > 1e3:
                        #             vv = float(da.sel(indexers={"member": k1, "dataset": k2}).values)
                        #             closest = copy__deepcopy(k1) if k1 > 1e3 else (
                        #                 min(list_tot, key=lambda x: abs(x - nn / 1e3)) if k2 > 1e3 else
                        #                 min(list_tot, key=lambda x: abs(x - nn)))
                        #             da_o.loc[{"member": closest, "dataset": k2}] = vv
                        # temporary save
                        if exp == "historical":
                            put_in_dict(to_plot, copy__deepcopy(da.compute()), dur, dia, pro, sta, exp)
                        else:
                            put_in_dict(dt, copy__deepcopy(da.compute()), exp)
                    # average across projections
                    da_o = tuple(da for exp, da in dt.items() if exp[:3] in ["rcp", "sce", "ssp"])
                    da_o = concatenate(da_o, "experiment").mean(dim="experiment", skipna=True)
                    # save
                    put_in_dict(to_plot, copy__deepcopy(da_o.compute()), dur, dia, pro, sta, "ssp")
    # min-max
    print("min-max")
    for dur in sorted(list(to_plot.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(to_plot[dur].keys()), key=lambda s: s.lower()):
            dt = {}
            for pro in list(to_plot[dur][dia].keys()):
                for sta in list(to_plot[dur][dia][pro].keys()):
                    for da in to_plot[dur][dia][pro][sta].values():
                        put_in_dict(dt, da.values.flatten().tolist(), sta)
            for sta in sorted(list(dt.keys()), key=lambda s: s.lower()):
                if len(dt[sta]) > 0:
                    print(dur, dia, sta, "{0:.2f}".format(round(min(dt[sta]), 2)).rjust(8),
                          "{0:.2f}".format(round(max(dt[sta]), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    colorbar = {
        "absolute": {
            "mea": [round(k / 10, 1) for k in range(0, 61, 15)],
            "spr": list(range(0, 33, 8)),
            "spt": list(range(0, 17, 4)),
        },
        "relative": {
            "mea": list(range(0, 51, 10)),
            "spr": list(range(0, 121, 30)),
            "spt": list(range(0, 81, 20)),
        },
    }
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in sorted(list(to_plot.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(to_plot[dur].keys()), key=lambda s: s.lower()):
            for pro in sorted(list(to_plot[dur][dia].keys()), key=lambda s: s.lower()):
                list_sta = sorted(list(to_plot[dur][dia][pro].keys()), key=lambda s: s.lower())
                for ii, sta in enumerate(list_sta):
                    list_exp = sorted(list(to_plot[dur][dia][pro][sta].keys()), key=lambda s: s.lower())
                    for jj, exp in enumerate(list_exp):
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = ("_".join([dur, dia, pro]), grp, pan)
                        print(dur, dia, pro, sta, exp.ljust(10), tuple_k)
                        # data array
                        da = to_plot[dur][dia][pro][sta][exp]
                        x_nbr, y_nbr = da["dataset"].size, da["member"].size
                        x_tic, y_tic = da["dataset"].values.tolist(), da["member"].values.tolist()
                        cbar = colorbar[proc_md][sta]
                        # panel size
                        tk = tuple_k[:-2]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k1, d1 in fig_panel_size.items():
                                if k1 == panel and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                                    for k2, d2 in d1.items():
                                        if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                            continue
                                        tmp = copy__deepcopy(d2)
                                        if k2 in ["x_size", "y_size"]:
                                            tmp *= (x_nbr if k2 == "x_size" else y_nbr) + 1
                                        put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                                elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                    put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                        # panel axes definition
                        if check_dict_levels(fig_ticks, dur, proc_ci, dia, panel):
                            for k1, k2 in fig_ticks[dur][proc_ci][dia][panel].items():
                                if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if k1 in ["x_lim", "y_lim"] and tmp is None:
                                    tmp = [0, x_nbr] if k1 == "x_lim" else [0, y_nbr]
                                elif (k1 == "y_nam" and jj > 0) or (k1 == "x_nam" and sta != list_sta[-1]):
                                    tmp = ""
                                elif k1 in ["x_tic", "y_tic"] and tmp is None:
                                    tmp = list(range(x_nbr)) if k1 == "x_tic" else list(range(y_nbr))
                                    tmp = [k + 0.5 for k in tmp]
                                elif k1 in ["z_lim"] and tmp is None:
                                    tmp = [min(cbar), max(cbar)]
                                elif k1 in ["z_tic"] and tmp is None:
                                    tmp = copy__deepcopy(cbar)
                                put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                                if k1 in ["x_tic", "y_tic"]:
                                    k3 = k1.replace("tic", "lab")
                                    tmp = copy__deepcopy(x_tic) if k1 == "x_tic" else copy__deepcopy(y_tic)
                                    tmp = [str(k) if k < 1000 else "max" for k in tmp]
                                    put_in_dict(figure_axes, tmp,*tuple_k + (k3,))
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
                            val = str(fig_titles[exp])
                            if exp in list(n_ensemble[pro].keys()):
                                val += " (" + str(n_ensemble[pro][exp].split(" (")[0]) + ")"
                            elif exp == "ssp":
                                vv = max([int(n_ensemble[pro][ex2].split(" (")[0])
                                          for ex2 in list(n_ensemble[pro].keys()) if ex2[:3] == exp])
                                val += " (" + str(vv) + ")"
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0:
                            val = fig_titles[sta] if sta in list(fig_titles.keys()) else copy__deepcopy(sta)
                            put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                        # -- heatmap
                        xx = numpy__array([k + 0.5 for k in range(x_nbr)])
                        yy = numpy__array([k + 0.5 for k in range(y_nbr)])
                        xx = numpy__concatenate(([xx[0] - (xx[1] - xx[0]) / 2], (xx[:-1] + xx[1:]) / 2,
                                                 [xx[-1] + (xx[-1] - xx[-2]) / 2]))
                        yy = numpy__concatenate(([yy[0] - (yy[1] - yy[0]) / 2], (yy[:-1] + yy[1:]) / 2,
                                                 [yy[-1] + (yy[-1] - yy[-2]) / 2]))
                        xx, yy = numpy__meshgrid(xx, yy)
                        da_v = da.where(da.notnull(), other=-9999)
                        plot_type = "hea"
                        da_z = da.where(da > min(cbar), other=min(cbar)).where(da < max(cbar), other=max(cbar))
                        put_in_dict(plot_data, True, *tuple_k + (str(plot_type) + "_a_bool",))
                        put_in_dict(plot_data, "dodgerblue", *tuple_k + (str(plot_type) + "_a_color",))
                        put_in_dict(plot_data, 13, *tuple_k + (str(plot_type) + "_a_font",))
                        put_in_dict(plot_data, "{x:.0f}", *tuple_k + (str(plot_type) + "_a_form",))
                        put_in_dict(plot_data, "cmo.amp", *tuple_k + (str(plot_type) + "_c",))
                        val = True if exp == list_exp[-1] else False
                        put_in_dict(plot_data, val, *tuple_k + (str(plot_type) + "_cb",))
                        put_in_dict(plot_data, "right", *tuple_k + (str(plot_type) + "_cl",))
                        put_in_dict(plot_data, copy__deepcopy(xx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(yy), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(da.values), *tuple_k + (str(plot_type) + "_z",))
                        put_in_dict(plot_data, copy__deepcopy(da_v.values), *tuple_k + (str(plot_type) + "_zv",))
                        # -- text
                        lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # plot title
                            if ii == 0 and jj == 0:
                                # plot size
                                l1 = len(list_exp)
                                x1_d = fig_panel_size[panel]["x_delt"]
                                x1_m = fig_panel_size[panel]["x_size"]
                                x1_s = x2 - x1
                                # add to text
                                lc += ["k"]
                                lf += [19]
                                lh += ["center"]
                                lr += [0]
                                lt += [fig_titles[proc_md].replace("THR", "-".join(
                                    [str(int(k4)) if k4 == int(k4) else str(k4) for k4 in percentiles]))]
                                lv += ["center"]
                                lw += ["bold"]
                                a1 = x1_m * x1_s
                                a2 = (a1 + x1_d) * l1 - x1_d
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
    a1 = str(smile_minimum_size).zfill(2)
    a3 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci)
    a4 = str(int(proc_cp)) if proc_cp == int(proc_cp) else str(proc_cp)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(a1) + "m", str(a3) + "p", str(proc_md), str(a4) + "ci"])
    # figure loop
    for kk in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[kk].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", kk)))
        # plot
        fig_basic(plot_data[kk], list_names, 1, figure_axes[kk], fig_format, figure_name_and_path,
                  figure_panel_size[kk], fig_title_bool=False, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e79_noise_vs_ensembles_plot()
