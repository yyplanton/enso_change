# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: scatterplots of bias (value at end of the 20th century) vs. bias, and bias vs. changes (departure from the end
# of the 20th century) in CMIP5 & CMIP6, raw & pcc
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from math import ceil as math__ceil
from math import floor as math__floor
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
import sys
from typing import Literal

# numpy
from numpy.random import randint as numpy__random__randint

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.numpy_tools import optimize_current_and_lagged_coefficients, optimize_lagged_coefficients
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.statistics import correlation_pearson, linear_regression, numpy__array, numpy__mean
from enso_change.wrapper.xarray_base import concatenate
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["ave_ts_val_glob", "ave_ts_val_gloo", "ave_ts_val_trop", "ave_ts_val_tron",
    #                       "ave_ts_val_n30e", "ave_ts_val_n34e", "ave_ts_val_n40e",
    #                       "ave_pr_val_n30e", "ave_pr_val_n34e", "ave_pr_val_n40e",
    #                       "std_ts_ano_n34e"],
    "data_diagnostics1": ["ave_ts_val_glob", "std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585", "ssp370"],
    # "data_experiments": ["historical", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},
    # list of projects: list[str]
    # "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    "data_projects": ["cmip5", "cmip6"],
    # "data_projects": ["cmip6"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    "smile_rejected": ["EC-EARTH", "MPI-ESM-LR"],
    # "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 25,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.2, "y": 0.2},
        "panel_01": {"x_delt": 7, "x_size": 20, "y_delt": 3, "y_size": 20},
    },
    # colors
    "fig_colors": {
        "cmip5": "green",
        "cmip6": "r",
        "observations": "k",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "s",
        "cmip6": "D",
        "observations": "o",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {},
    # titles
    "fig_titles": {
        "bia": "bias",
        "pcc": "change per GMST",
        "raw": "change",
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
        "title_row_x": 40,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e72_trend_vs_enso_plot(
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
    print("reader_json", str(len(list(values_json.keys()))).rjust(3))
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
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_experiments = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(12), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
    #
    # -- Processing
    #
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # -- change
    change_em = {}
    for dia1 in list(values_json.keys()):
        for dur1 in list(values_json[dia1].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia1][dur1].keys()):
                for exp1 in list(values_json[dia1][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical", exp2]:
                        continue
                    for dat in list(values_json[dia1][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia1, dur2, pro, exp1, dat) or \
                                not check_dict_levels(values_json, dia1, dur2, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia1][dur1][pro][exp1][dat])
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia1][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia1][dur2][pro][exp2][dat])
                        # select experiment epoch
                        da1 = da1.mean(dim="member")
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 2006}).mean(dim="member")
                        # change
                        da1 = da1 - da2
                        if dia1.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        exp_o = copy__deepcopy(exp1) if exp1 in ["piControl", "historical", "ssp370"] else \
                            "scenario" + str(exp1[-2:])
                        put_in_dict(change_em, da1, dur1, dia1, pro, exp_o, dat)
    print("change_em", list(change_em.keys()))
    # min-max
    print("change_em min-max")
    change_minmax = {}
    for dur in list(change_em.keys()):
        for dia in list(change_em[dur].keys()):
            tt = []
            for pro in list(change_em[dur][dia].keys()):
                for exp in list(change_em[dur][dia][pro].keys()):
                    for da in change_em[dur][dia][pro][exp].values():
                        tt += [float(da.min().values), float(da.max().values)]
            if len(tt) > 0:
                print(dur, dia.ljust(22), "{0:.2f}".format(round(min(tt), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tt), 2)).rjust(8))
                put_in_dict(change_minmax, copy__deepcopy(tt), dur, dia)
    # -- explain with mean change
    dia2 = "std_ts_ano_n34e"
    mc_coefficients, mc_r, mc_mse = {}, {}, {}
    for dur1 in list(change_em.keys()):
        if dur1 == dur2:
            continue
        for dia1 in list(change_em[dur1].keys()):
            if dia1 == dia2:
                continue
            for pro in list(change_em[dur1][dia1].keys()):
                for exp in list(change_em[dur1][dia1][pro].keys()):
                    if exp in ["piControl", "historical"]:
                        continue
                    for dat in list(change_em[dur1][dia1][pro][exp].keys()):
                        if not check_dict_levels(change_em, dur1, dia2, pro, exp, dat):
                            continue
                        # array
                        da_x = copy__deepcopy(change_em[dur1][dia1][pro][exp][dat])
                        da_y = copy__deepcopy(change_em[dur1][dia2][pro][exp][dat])
                        # select years with uniform delta between years
                        years_x, years_y = da_x["year"].values.tolist(), da_y["year"].values.tolist()
                        years = sorted(list(set(years_x) & set(years_y)))
                        dy = max(years_x[1] - years_x[0], years_y[1] - years_y[0])
                        years = [k for k in range(min(years), max(years) + 1, dy) if k in years]
                        da_x, da_y = da_x.sel(indexers={"year": years}), da_y.sel(indexers={"year": years})
                        # optimize lagged regression
                        cc, error = optimize_lagged_coefficients(da_x.values, da_y.values, max_lag=int(60 / dy))
                        a, b, lag = cc
                        if lag == 0:
                            y_pred = a * da_x + b
                        else:
                            y_pred = a * da_x[:-lag] + b
                        y_val = da_y[lag:]
                        r = correlation_pearson(y_pred, y_val)
                        # save
                        put_in_dict(mc_coefficients, cc, dur1, dia1, pro, exp, dat)
                        put_in_dict(mc_r, float(r.values), dur1, dia1, pro, exp, dat)
                        put_in_dict(mc_mse, error, dur1, dia1, pro, exp, dat)
    print("correlations: explain with mean change")
    for dur in list(mc_r.keys()):
        for dia in list(mc_r[dur].keys()):
            print(dur, dia)
            dt1, dt2 = {}, {}
            for pro in sorted(list(mc_r[dur][dia].keys()), key=lambda s: s.lower()):
                for exp in sorted(list(mc_r[dur][dia][pro].keys()), key=lambda s: s.lower()):
                    put_in_dict(dt1, [v for v in mc_r[dur][dia][pro][exp].values()], pro, exp)
                    put_in_dict(dt2, [v for v in mc_r[dur][dia][pro][exp].values()], exp)
            print(dur, dia, "per project per experiments")
            for pro in list(dt1.keys()):
                for exp, vv in dt1[pro].items():
                    tt = numpy__array(vv)
                    print(pro.rjust(20), exp.ljust(10), "r=" + "{0:.3f}".format(round(float(tt.mean()), 3)).rjust(6),
                          "(" + "{0:.3f}".format(round(float(tt.min()), 3)).rjust(6) + " to " +
                          "{0:.3f}".format(round(float(tt.max()), 3)).rjust(6) + ")")
            print(dur, dia, "per experiments")
            for exp, vv in dt2.items():
                tt = numpy__array(vv)
                print("".rjust(20), exp.ljust(10), "r=" + "{0:.3f}".format(round(float(tt.mean()), 3)).rjust(6),
                      "(" + "{0:.3f}".format(round(float(tt.min()), 3)).rjust(6) + " to " +
                      "{0:.3f}".format(round(float(tt.max()), 3)).rjust(6) + ")")
    # -- explain with mean change and lag
    mcl_coefficients, mcl_r, mcl_mse = {}, {}, {}
    for dur1 in list(change_em.keys()):
        if dur1 == dur2:
            continue
        for dia1 in list(change_em[dur1].keys()):
            if dia1 == dia2:
                continue
            for pro in list(change_em[dur1][dia1].keys()):
                for exp in list(change_em[dur1][dia1][pro].keys()):
                    for dat in list(change_em[dur1][dia1][pro][exp].keys()):
                        if not check_dict_levels(change_em, dur1, dia2, pro, exp, dat):
                            continue
                        # array
                        da_x = copy__deepcopy(change_em[dur1][dia1][pro][exp][dat])
                        da_y = copy__deepcopy(change_em[dur1][dia2][pro][exp][dat])
                        # select years with uniform delta between years
                        years_x, years_y = da_x["year"].values.tolist(), da_y["year"].values.tolist()
                        years = sorted(list(set(years_x) & set(years_y)))
                        dy = max(years_x[1] - years_x[0], years_y[1] - years_y[0])
                        years = [k for k in range(min(years), max(years) + 1, dy) if k in years]
                        da_x, da_y = da_x.sel(indexers={"year": years}), da_y.sel(indexers={"year": years})
                        # coefficients
                        cc, error = optimize_current_and_lagged_coefficients(
                            da_x.values, da_y.values, max_lag=int(60 / dy))
                        y_pred = cc[0] * da_x[cc[3]:] + cc[1] * (da_x[cc[3]:] - da_x[: -cc[3]]) + cc[2]
                        y_val = da_y[cc[3]:]
                        r = correlation_pearson(y_pred, y_val)
                        if float(r.values):
                            print(dur1, dia1, pro, exp, dat, float(r.values))
                            stop
                        # save
                        put_in_dict(mcl_coefficients, cc, dur1, dia1, pro, exp, dat)
                        put_in_dict(mcl_r, float(r.values), dur1, dia1, pro, exp, dat)
                        put_in_dict(mcl_mse, error, dur1, dia1, pro, exp, dat)
    print("correlations: explain with mean change and lag")
    for dur in list(mcl_r.keys()):
        for dia in list(mcl_r[dur].keys()):
            print(dur, dia)
            dt1, dt2 = {}, {}
            for pro in sorted(list(mcl_r[dur][dia].keys()), key=lambda s: s.lower()):
                for exp in sorted(list(mcl_r[dur][dia][pro].keys()), key=lambda s: s.lower()):
                    put_in_dict(dt1, [v for v in mcl_r[dur][dia][pro][exp].values()], pro, exp)
                    put_in_dict(dt2, [v for v in mcl_r[dur][dia][pro][exp].values()], exp)
            print(dur, dia, "per project per experiments")
            for pro in list(dt1.keys()):
                for exp, vv in dt1[pro].items():
                    tt = numpy__array(vv)
                    print(pro.rjust(20), exp.ljust(10), "r=" + "{0:.3f}".format(round(float(tt.mean()), 3)).rjust(6),
                          "(" + "{0:.3f}".format(round(float(tt.min()), 3)).rjust(6) + " to " +
                          "{0:.3f}".format(round(float(tt.max()), 3)).rjust(6) + ")")
            print(dur, dia, "per experiments")
            for exp, vv in dt2.items():
                tt = numpy__array(vv)
                print("".rjust(20), exp.ljust(10), "r=" + "{0:.3f}".format(round(float(tt.mean()), 3)).rjust(6),
                      "(" + "{0:.3f}".format(round(float(tt.min()), 3)).rjust(6) + " to " +
                          "{0:.3f}".format(round(float(tt.max()), 3)).rjust(6) + ")")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: lag vs. current_lag
    print("lag vs. current_lag")
    group, panel = "group_1", "panel_01"
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e72_trend_vs_enso_plot()
