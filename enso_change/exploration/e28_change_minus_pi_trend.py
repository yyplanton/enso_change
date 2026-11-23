# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Evolution of change members by members and ensemble mean when piControl is removed in SMILEs
#    * time (x) vs. change [AVE & SKE] or relative change [STD] (y) per experiment
# -> fail:
#    - MIROC6: 28 over 50 members branch late from piControl (according to metadata in netCDF files) so there is not
#              enough piControl left to cover the 250 years of historical+SSP, and 7 members branch too late to even
#              cover the 165 years of historical.
#    - MPI-ESM1-2-LR: 21 over 50 members branch at the same year (1850-01-01) without explanation so I think that there
#                     is an issue with the metadata
# I don't think that it is possible to detrend members using piControl
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
from enso_change.script_figures.branch_time import branch_time
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.xarray_base import array_wrapper, numpy__array, concatenate
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
    "data_experiments": ["piControl", "historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 50,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": True,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 90.,
    # method: {"combination", "random"}
    "proc_md": "random",
    # number of pseudo samples: int
    "proc_ps": 100000,
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
        "panel_1": {"x_delt": 250, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 150, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from 1850-YEAR",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-0.4, 0.4],
                    "nam": "departure from 1850-YEAR",
                    "tic": [round(k / 10, 1) for k in range(-4, 5, 2)],
                },
            },
            # "std_ts_ano_n34e": {
            #     "panel_1": {
            #         "lim": [-40, 40],
            #         "nam": "departure from 1850-YEAR",
            #         "tic": list(range(-40, 41, 20)),
            #     },
            #     "panel_2": {
            #         "lim": [-40, 40],
            #         "nam": "departure from YEAR-2014",
            #         "tic": list(range(-40, 41, 20)),
            #     },
            # },
            # "std_ts_ano_n34e": {
            #     "panel_1": {
            #         "lim": [-50, 50],
            #         "nam": "departure from 1850-YEAR",
            #         "tic": list(range(-50, 51, 25)),
            #     },
            #     "panel_2": {
            #         "lim": [-50, 50],
            #         "nam": "departure from YEAR-2014",
            #         "tic": list(range(-50, 51, 25)),
            #     },
            # },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-60, 60],
                    "nam": "departure from 1850-YEAR",
                    "tic": list(range(-60, 61, 30)),
                },
                "panel_2": {
                    "lim": [-60, 60],
                    "nam": "departure from YEAR-2014",
                    "tic": list(range(-60, 61, 30)),
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1850, 2015],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(1850, 2016, 50)),
                },
                "panel_2": {
                    "lim": [1970, 2100],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(2000, 2101, 50)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "change of VALUE (UNITS)",
        "ske_ts_ano_n34e": "change of VALUE",
        "std_ts_ano_n34e": "relative change of VALUE (%)",
        "full": "full variability",
        "ivc": "without IV (corrected)",
        "ivp": "without IV (piControl)",
        "ivr": "without IV",
        "mean": "MMM",
        "range": "THR%",
        "figure_title": "MMM's spread with or without internal variability",
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
def e28_change_minus_pi_trend_plot(
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
        data_observational_datasets=data_observational_datasets, data_projects=data_projects, smile_array=False,
        smile_minimum_size=smile_minimum_size, smile_rejected=smile_rejected,
        smile_require_picontrol=smile_require_picontrol)
    # concatenate historical and SSPs
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
                    # delete datasets
                    if "piControl" not in list_exp or "historical" not in list_exp or \
                            dat not in list(branch_time.keys()):
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                    # concatenate historical and SSPs
                    for exp in list_ssp:
                        # dictionaries
                        d2_h = copy__deepcopy(values_json[dia][dur][pro]["historical"][dat])
                        d2_s = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        # list members in both historical and ssp
                        list_members = sorted(list(set(list(d2_h.keys())) & set(list(d2_s.keys()))))
                        for mem in list_members:
                            da_h, da_s = d2_h[mem], d2_s[mem]
                            del values_json[dia][dur][pro][exp][dat][mem]
                            # concatenate along year axis
                            da_o = concatenate((da_h, da_s), "year")
                            # resave
                            put_in_dict(values_json, copy__deepcopy(da_o), dia, dur, pro, exp, dat, mem)
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
        print(exp.ljust(10), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
    #
    # -- Processing
    #
    # -- remove piControl time series corresponding to each member
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                # -- remove piControl time series
                # dictionary
                d1 = values_json[dia][dur][pro]
                # list datasets
                list_dataset = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                for dat in list_dataset:
                    if dat == "MPI-ESM1-2-LR":
                        l1 = [k for k in branch_time[dat].values()]
                        l2 = [k for k in l1 if k == 1850]
                        print(dia, dur, pro, dat, len(l1), len(list(set(l1))), len(l2))
                        stop
                #     # list experiments
                #     list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys()) if k1 != "piControl"]
                #     for exp in list_exp:
                #         if not check_dict_levels(values_json, dia, dur, pro, "piControl", dat):
                #             continue
                #         # dictionaries
                #         d2_e = values_json[dia][dur][pro][exp][dat]
                #         d2_p = values_json[dia][dur][pro]["piControl"][dat]
                #         for me1 in list(d2_e.keys()):
                #             if not check_dict_levels(branch_time, dat, me1):
                #                 continue
                #             bt = copy__deepcopy(branch_time[dat][me1]) + int(dur.split("_")[0])
                #             # array
                #             da_e = copy__deepcopy(d2_e[me1])
                #             # find piControl member corresponding to given branch time
                #             da_p = None
                #             for me2 in list(d2_p.keys()):
                #                 if bt in d2_p[me2]["year"].values:
                #                     da_p = copy__deepcopy(d2_p[me2])
                #             if not isinstance(da_p, array_wrapper):
                #                 continue
                #             # select piControl years corresponding to given exp years
                #             ye_e = da_e["year"].values
                #             print(dia, dur, pro, exp, dat, me1, ye_e + bt - ye_e[0])
                #             da_p = da_p.sel(indexers={"year": ye_e + bt - ye_e[0]})
                #             da_p = da_p.assign_coords({"year": ye_e})
                #             # remove piControl from given exp
                #             da_e -= da_p
                #             # save
                #             exp_o = str(exp) + "-pi"
                #             put_in_dict(values_json, copy__deepcopy(da_e), dia, dur, pro, exp_o, dat, me1)
                # # -- create array A[mem, yea]
                # for exp in list(values_json[dia][dur][pro].keys()):
                #     for dat in list(values_json[dia][dur][pro][exp].keys()):
                #         # concatenate members
                #         list_members = list(values_json[dia][dur][pro][exp][dat].keys())
                #         da = tuple(values_json[dia][dur][pro][exp][dat][mem] for mem in list_members)
                #         da = concatenate(da, "member")
                #         del values_json[dia][dur][pro][exp][dat]
                #         # assign member coordinates
                #         list_members = [int(k[1:].split("i")[0]) for k in list_members]
                #         da = da.assign_coords({"year": numpy__array(list_members)})
                #         # save
                #         put_in_dict(values_json, copy__deepcopy(da), dia, dur, pro, exp, dat)
    # -- compute change compared to initial conditions
    dur2, exp2 = data_epoch_lengths[-1], "piControl"
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
                        if not check_dict_levels(values_json, dia, dur2, pro, exp1) or \
                                not check_dict_levels(values_json, dia, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        da3 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        yy = 1850 + int(dur2.split("_")[0]) if "historical"in exp1 else 2014
                        da2, da3 = da2.sel(indexers={"year": yy}), da3.sel(indexers={"year": yy})
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        for typ in ["1m", "em"]:
                            da_o = copy__deepcopy(da1)
                            da22 = copy__deepcopy(da2) if typ == "1m" else da2.mean(dim="member")
                            da33 = copy__deepcopy(da3) if typ == "1m" else da3.mean(dim="member")
                            # change
                            da_o -= da22
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o *= 100 / da33
                            # save
                            put_in_dict(change, da1, dur1, dia, dat, typ, exp1)
    # min-max
    print("change min-max")
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            tmp = []
            for dat in list(change[dur][dia].keys()):
                for typ in list(change[dur][dia][dat].keys()):
                    for da in change[dur][dia][dat][typ].values():
                        tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e28_change_minus_pi_trend_plot()
