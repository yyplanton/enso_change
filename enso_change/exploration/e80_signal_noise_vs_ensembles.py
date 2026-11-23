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
from math import ceil as math__ceil
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
import sys
from typing import Literal

# numpy
from numpy import concatenate as numpy__concatenate
from numpy import meshgrid as numpy__meshgrid

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper import numpy_tools
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import numpy__nan, array_wrapper, concatenate, polyval
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
    # confidence interval as percentile: float
    "proc_cp": 90,
    # method: {"absolute", "relative"}
    "proc_md": "relative",
    # number of pseudo samples: int
    "proc_ps": 300,
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
                        "y_nam": "total nbr of members",
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
        "mea": "multi-model mean",
        "spr": "multi-model spread",
        "spt": "multi-model spread - RIV",
        "true": "full spread - RIV",
        "raw": "full spread",
        "sha": "95% CI",
        # variables
        "std_ts_ano": "relative change of SSTA $\\sigma$",
        # figure title
        "figure_title": "Relative error on Inter-model mean & spread (THR%)",
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
def e80_signal_noise_vs_ensembles_plot(
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
    print("e80_signal_noise_vs_ensembles_plot input")
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
                        # list_nm = list(range(1, min(size + 1, 11), 1)) + ["max"]
                        list_nm = list(range(1, min(size + 1, 11), 1)) + ["max"]
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
    def stratified_sample(data1, data2, ensemble_sizes, n_samples, seed=None):
        """
        Select sample indices from ensembles, distributing
        samples as evenly as possible.
        ensemble_sizes: list of available samples per ensemble
        n_samples: total number of samples to select
        """
        if n_samples == "max":
            out1 = tuple(v.sel(indexers={"member": int(1e6)}) for k, v in data1.items())
            out2 = tuple(v.sel(indexers={"member": int(1e6)}) for k, v in data2.items())
        else:
            if n_samples > sum(ensemble_sizes.values()):
                raise ValueError("Not enough available samples")
            allocation = {key: 0 for key in ensemble_sizes}
            remaining = n_samples
            active = dict(ensemble_sizes)
            while remaining > 0:
                # equal quota among remaining ensembles
                quota = remaining // len(active)
                # handle cases where remaining < number of ensembles
                if quota == 0:
                    quota = 1
                for key in list(active):
                    available = active[key] - allocation[key]
                    take = min(quota, available, remaining)
                    allocation[key] += take
                    remaining -= take
                    if allocation[key] == active[key]:
                        del active[key]
                    if remaining == 0:
                        break
            out1 = tuple(data1[k].sel(indexers={"member": v}) for k, v in allocation.items())
            out2 = tuple(data2[k].sel(indexers={"member": v}) for k, v in allocation.items())
        out1, out2 = concatenate(out1, "dataset"), concatenate(out2, "dataset")
        return out1, out2
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_sta = {}
    for dur in list(mean_change.keys()):
        ny = int(dur.split("_")[0])
        if dur == dur2:
            continue
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    if pro == "observations" or exp == "piControl":
                        continue
                    # concatenate member resamples
                    dict_t1, dict_t2, dict_t3 = {}, {}, {}
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        if not check_dict_levels(change, dur, dia, pro, exp, dat) or \
                                not check_dict_levels(residual_iv, dur, dia, pro, exp, dat):
                            continue
                        list_nme = sorted(list(mean_change[dur][dia][pro][exp][dat].keys()), key=lambda s: s.lower())
                        da1 = tuple(mean_change[dur][dia][pro][exp][dat][k] for k in list_nme)
                        da2 = tuple(residual_iv[dur][dia][pro][exp][dat][k] for k in list_nme)
                        da1, da2 = concatenate(da1, "member"), concatenate(da2, "member")
                        list_nme = numpy__array([int(k[1:]) if k[1:].isdigit() else int(1e6) for k in list_nme])
                        da1, da2 = da1.assign_coords({"member": list_nme}), da2.assign_coords({"member": list_nme})
                        y1, y2 = (1850 if exp == "historical" else 2000) + ny, 1950 if exp == "historical" else 2100
                        da1 = da1.sel(indexers={"year": slice(y1, y2)})
                        dict_t1[dat] = da1
                        dict_t2[dat] = da2
                        dict_t3[dat] = max([k for k in list_nme if k < 1e5])
                    # number of datasets
                    dat_as_array = numpy__array(list(dict_t1.keys()))
                    nbr = dat_as_array.size
                    list_nda = list(range(20, min(nbr, 51), 5)) + ["max"]
                    # list_nda = list(range(20, 26, 5)) + ["max"]
                    for nda in list_nda:
                        if nda == "max":
                            resamples = ["max"]
                        else:
                            idx = numpy_tools.combination_indices(nbr, int(proc_ps * 1.5), nda)
                            resamples = dat_as_array[idx]
                        dt1, dt2 = {}, {}
                        for idx in resamples:
                            test1 = isinstance(idx, str) and idx == "max"
                            n_datasets = copy__deepcopy(nbr) if test1 or nda == "max" else copy__deepcopy(nda)
                            if test1:
                                d1, d2, d3 = dict_t1, dict_t2, dict_t3
                                mini = math__ceil(n_datasets / 5) * 5
                            else:
                                d1, d2 = {k: dict_t1[k] for k in idx}, {k: dict_t2[k] for k in idx}
                                d3 = {k: dict_t3[k] for k in idx}
                                mini = copy__deepcopy(n_datasets)
                            maxi = min(max(min(n_datasets * 5, sum(d3.values())), n_datasets), 150)
                            list_n = [k for k in range(max(30, mini), maxi + 1, 5) if (k <= 30 and k % 5 == 0) or
                                      (30 < k <= 50 and k % 10 == 0) or (50 < k <= 100 and k % 15 == 0) or
                                      (110 <= k <= 150 and (k - 10) % 20 == 0) or (150 < k <= 200 and k % 25 == 0) or
                                      (200 < k <= 300 and k % 50 == 0)]
                            for nme in list_n + ["max"]:
                                da1, da2 = stratified_sample(d1, d2, d3, nme)
                                nm_o = "m" + (nme if nme == "max" else str(nme).zfill(3))
                                dt1[nm_o] = dt1[nm_o] + (da1,) if nm_o in list(dt1.keys()) else (da1,)
                                dt2[nm_o] = dt2[nm_o] + (da2,) if nm_o in list(dt2.keys()) else (da2,)
                        for nme in sorted(list(dt1.keys()), key=lambda s: s.lower()):
                            if len(dt1[nme]) < proc_ps and nda != "max" and nme != "max":
                                continue
                            da1, da2 = concatenate(dt1[nme], "re_dat"), concatenate(dt2[nme], "re_dat")
                            nd_o = "d" + (nda if nda == "max" else str(nda).zfill(3))
                            print(dur, dia, pro, exp.ljust(10), nd_o, nme, da1.shape, da1.dims)
                            # stack datasets and members resamples
                            da1 = da1.stack(dim={"resample": ("re_dat", "re_mem")})
                            dims = list(set(list(da1.coords)) - set(list(da1.dims)))
                            if len(dims) > 0:
                                da1 = da1.drop_vars(dims)
                            da2 = da2.stack(dim={"resample": ("re_dat", "re_mem")})
                            dims = list(set(list(da2.coords)) - set(list(da2.dims)))
                            if len(dims) > 0:
                                da2 = da2.drop_vars(dims)
                            # -- scaled distribution
                            # initial multimodel mean using ensemble mean -> will not change
                            mi = da1.mean(dim="dataset", skipna=True)
                            # initial multimodel spread using ensemble mean -> that what will be scaled
                            si = da1.var(dim="dataset", skipna=True)
                            # initial multimodel mean residual internal variability -> will be used for scaling
                            ri = da2.mean(dim="dataset", skipna=True)
                            # final spread (initial spread minus residual internal variability) -> must be set to 0
                            # if negative
                            sf = (si - ri) / si
                            sf = sf.where(sf > 0, other=0)
                            # scaled distribution
                            da_o = mi + (da1 - mi) * sf
                            # -- mean
                            mi = mi.mean(dim="year")
                            lk = (dur, dia, pro, "mea", exp, nme, nd_o)
                            put_in_dict(change_sta, copy__deepcopy(mi.compute()), *lk)
                            # -- true model spread
                            da1 = statistic(da1, dim="dataset", stat="per", per=percentiles, skipna=True)
                            da1 = da1.max(dim="quantile") - da1.min(dim="quantile")
                            da1 = da1.mean(dim="year")
                            dims = list(set(list(da1.coords)) - set(list(da1.dims)))
                            if len(dims) > 0:
                                da1 = da1.drop_vars(dims)
                            lk = (dur, dia, pro, "spr", exp, nme, nd_o)
                            put_in_dict(change_sta, copy__deepcopy(da1.compute()), *lk)
                            # -- true model spread
                            da_o = statistic(da_o, dim="dataset", stat="per", per=percentiles, skipna=True)
                            da_o = da_o.max(dim="quantile") - da_o.min(dim="quantile")
                            da_o = da_o.mean(dim="year")
                            dims = list(set(list(da_o.coords)) - set(list(da_o.dims)))
                            if len(dims) > 0:
                                da_o = da_o.drop_vars(dims)
                            lk = (dur, dia, pro, "spt", exp, nme, nd_o)
                            put_in_dict(change_sta, copy__deepcopy(da_o.compute()), *lk)
    print("change_sta", list(change_sta.keys()))
    # -- error
    to_plot = {}
    for dur in list(change_sta.keys()):
        for dia in list(change_sta[dur].keys()):
            for pro in list(change_sta[dur][dia].keys()):
                for sta in list(change_sta[dur][dia][pro].keys()):
                    # ensures that the list of total number of members is the same for all experiments
                    list_nme_max = []
                    for exp in list(change_sta[dur][dia][pro][sta].keys()):
                        list_nme_max += list(change_sta[dur][dia][pro][sta][exp].keys())
                    list_nme_max = sorted(list(set(list_nme_max)), key=lambda s: s.lower())
                    # compute error on each experiments
                    dt1 = {}
                    for exp in list(change_sta[dur][dia][pro][sta].keys()):
                        dt2 = change_sta[dur][dia][pro][sta][exp]
                        # the number of dataset resamples is not the same for all total number of members so the error
                        # percentile must be computed before concatenation
                        nme = sorted(list(dt2.keys()), key=lambda s: s.lower())[-1]
                        nda = sorted(list(dt2[nme].keys()), key=lambda s: s.lower())[-1]
                        da_r = float(dt2[nme][nda].mean().values)
                        # matrix of resampled data
                        da_o = ()
                        list_nda = sorted(list(set([k2 for k1 in list(dt2.keys()) for k2 in list(dt2[nme].keys())])),
                                          key=lambda s: s.lower())
                        for nme in list_nme_max:
                            da_t = ()
                            for nda in list_nda:
                                # nan if given nme does not exist for this model
                                da = array_wrapper(numpy__nan)
                                if check_dict_levels(dt2, nme, nda):
                                    # error distribution
                                    da = abs(dt2[nme][nda] - da_r)
                                    if proc_md == "relative":
                                        da *= 100 / abs(da_r)
                                    # percentiles of error
                                    da = statistic(da, dim="resample", stat="per", per=proc_cp, skipna=True)
                                    dims = list(set(list(da.coords)) - set(list(da.dims)))
                                    if len(dims) > 0:
                                        da = da.drop_vars(dims)
                                da_t += (da,)
                            da_t = concatenate(da_t, "dataset")
                            coordinates = numpy__array([int(k[1:]) if k[1:].isdigit() else 1000 for k in list_nda])
                            da_o += (da_t.assign_coords({"dataset": coordinates}),)
                        da_o = concatenate(da_o, "member")
                        coordinates = numpy__array([int(k[1:]) if k[1:].isdigit() else 1000 for k in list_nme_max])
                        da_o = da_o.assign_coords({"member": coordinates})
                        # temporary save
                        if exp == "historical":
                            put_in_dict(to_plot, copy__deepcopy(da_o.compute()), dur, dia, pro, sta, exp)
                        else:
                            put_in_dict(dt1, copy__deepcopy(da_o.compute()), exp)
                    # average across projections
                    da_o = tuple(da for exp, da in dt1.items() if exp[:3] in ["rcp", "sce", "ssp"])
                    da_o = concatenate(da_o, "experiment").mean(dim="experiment", skipna=True)
                    dims = list(set(list(da_o.coords)) - set(list(da_o.dims)))
                    if len(dims) > 0:
                        da_o = da_o.drop_vars(dims)
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
            "mea": list(range(0, 5, 1)),  # [round(k / 10, 1) for k in range(0, 61, 15)],
            "spr": list(range(0, 25, 6)),
            "spt": list(range(0, 13, 3)),
        },
        "relative": {
            "mea": list(range(0, 51, 10)),
            "spr": list(range(0, 101, 25)),
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
                        plot_type = "hea"
                        put_in_dict(plot_data, True, *tuple_k + (str(plot_type) + "_a_bool",))
                        put_in_dict(plot_data, "dodgerblue", *tuple_k + (str(plot_type) + "_a_color",))
                        put_in_dict(plot_data, "{x:.0f}", *tuple_k + (str(plot_type) + "_a_form",))
                        put_in_dict(plot_data, "cmo.amp", *tuple_k + (str(plot_type) + "_c",))
                        val = True if exp == list_exp[-1] else False
                        put_in_dict(plot_data, val, *tuple_k + (str(plot_type) + "_cb",))
                        put_in_dict(plot_data, "right", *tuple_k + (str(plot_type) + "_cl",))
                        put_in_dict(plot_data, copy__deepcopy(xx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(yy), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(da.values), *tuple_k + (str(plot_type) + "_z",))
                        put_in_dict(plot_data, copy__deepcopy(da.values), *tuple_k + (str(plot_type) + "_zv",))
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
                                lt += [fig_titles["figure_title"].replace("THR", "-".join(
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
e80_signal_noise_vs_ensembles_plot()
