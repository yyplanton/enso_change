# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# NetCDF: takes json data and create ensemble mean time series
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from inspect import stack as inspect__stack
from os.path import dirname as os__path__dirname
from os.path import join as os__path__join
import sys
from typing import Literal

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, print_fail, sort_dict
from enso_change.wrapper.processors import reader_json, compute_1var
from enso_change.wrapper.xarray_base import numpy__array, concatenate, to_netcdf
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["ave_pr_val_glob", "ave_pr_val_tron", "ave_ts_val_glob", "ave_ts_val_tron",
                          "std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e",
                          "std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e"],
    "data_diagnostics2": [],
    # "data_diagnostics2": ["ave_pr_val", "ave_sst_val"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # "data_epoch_lengths": ["030_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["observations"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    "smile_rejected": ["EC-EARTH"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 100.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 40,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {},
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {},
    # titles
    "fig_titles": {},
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {},
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def c03_create_ensemble_mean_time_series_1d_obs_plot(
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
    print("reader_json", len(list(values_json.keys())), list(values_json.keys()))
    n_ensemble = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                if pro != "observations":
                    continue
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), pro, exp)
    for pro in data_projects:
        for exp in data_experiments:
            if not check_dict_levels(n_ensemble, pro, exp):
                continue
            list_ensembles = copy__deepcopy(n_ensemble[pro][exp])
            del n_ensemble[pro][exp]
            list_ensembles = sorted(list(set(list_ensembles)), key=lambda s: s.lower())
            put_in_dict(n_ensemble, copy__deepcopy(list_ensembles), pro, exp)
            print(pro, exp.ljust(10), "number of ensembles: " + str(len(list_ensembles)).rjust(2))
            if exp == "historical":
                print(list_ensembles)
    dict_datasets = {
        "20CRv2c": 1,
        "20CRv3": 2,
        "CMAP": 3,
        "COBE2": 4,
        "ERA5": 5,
        "ERSSTv5": 6,
        "GODAS": 7,
        "GPCPv2.3": 8,
        "HadISST": 9,
        "NCEP1": 10,
        "ORAS5": 11,
        "SODA3.15.2": 12,
    }
    #
    # -- Processing
    #
    # ensemble mean
    mean_json, _, _ = compute_1var(values_json, "statistic",
                             function_arguments={"dim": "member", "stat": "mea"})
    # regroup experiments
    print("regroup experiments")
    mean_json_exp1 = {}
    for dia in list(mean_json.keys()):
        for dur in list(mean_json[dia].keys()):
            for pro in list(mean_json[dia][dur].keys()):
                for exp in list(mean_json[dia][dur][pro].keys()):
                    for dat in list(mean_json[dia][dur][pro][exp].keys()):
                        if pro == "observations" and (
                                (dia[:10] == "ave_pr_val" and (
                                        (dur in ["030_year_epoch", "050_year_epoch"] and
                                         dat not in ["20CRv2c", "20CRv3", "ERA5", "NCEP1"]) or
                                        (dur in ["045_year_epoch"] and dat not in ["CMAP", "GPCPv2.3"]))) or
                                (dia[:10] == "ave_ts_val" and dat not in ["COBE2", "ERSSTv5", "HadISST"])):
                            continue
                        # array
                        da1 = values_json[dia][dur][pro][exp][dat]
                        da2 = mean_json[dia][dur][pro][exp][dat]
                        ensemble_size = da1["member"].size
                        da2.attrs.update({"ensemble_size": str(dat) + "=" + str(ensemble_size)})
                        # save
                        pro_o = copy__deepcopy(pro[:4])
                        exp_o = "scenario" + str(exp[-2:]) if exp[:3] in ["rcp", "ssp"] else copy__deepcopy(exp)
                        # for size in [1, 2, 3, 4, 5, 10]:
                        # for size in [1, 2, 3, 4, 5]:
                        for size in [1]:
                            if ensemble_size >= size:
                                size_o = str(size).zfill(3)
                                put_in_dict(mean_json_exp1, da2, dia, dur, size_o, pro_o, exp_o, dat)
    print("regroup experiments: done", len(list(mean_json_exp1.keys())))
    # select datasets for piControl and another experiment
    # output plot directory (relative to current file directory)
    data_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/data"
    print("select datasets for piControl and another experiment")
    for dia in list(mean_json_exp1.keys()):
        for dur in list(mean_json_exp1[dia].keys()):
            for siz in list(mean_json_exp1[dia][dur].keys()):
                for pro in list(mean_json_exp1[dia][dur][siz].keys()):
                    print(dia, dur, siz, pro)
                    # -- create multiple piControl experiments corresponding to models in each experiment
                    # list experiments
                    list_experiments = [k for k in list(mean_json_exp1[dia][dur][siz][pro].keys()) if k != "piControl"]
                    dt = {}
                    for exp in list_experiments:
                        for dat in list(mean_json_exp1[dia][dur][siz][pro][exp].keys()):
                            if not check_dict_levels(mean_json_exp1, dia, dur, siz, pro, "piControl", dat) and \
                                    pro == "cmip":
                                print(dia, dur, siz, exp.ljust(10), dat.ljust(17), "missing piControl")
                                continue
                            # arrays
                            da1 = mean_json_exp1[dia][dur][siz][pro][exp][dat]
                            da2 = None
                            if check_dict_levels(mean_json_exp1, dia, dur, siz, pro, "piControl", dat):
                                da2 = mean_json_exp1[dia][dur][siz][pro]["piControl"][dat]
                            # save
                            put_in_dict(dt, da1, exp, dat)
                            if da2 is not None:
                                put_in_dict(dt, da2, "piControl-" + str(exp), dat)
                    # -- regroup datasets in a single array
                    for exp in list(dt.keys()):
                        # concatenate
                        list_dat = [dat for dat in list(dict_datasets.keys()) if dat in list(dt[exp].keys())]
                        da = concatenate(tuple(dt[exp][dat] for dat in list_dat), "dataset")
                        da = da.assign_coords({"dataset": numpy__array([dict_datasets[dat] for dat in list_dat])})
                        # attributes
                        dat = list_dat[0]
                        att = str(dict_datasets[dat]).zfill(3) + " " + str(dt[exp][dat].attrs["ensemble_size"])
                        for dat in list_dat[1:]:
                            att += "; " + str(dict_datasets[dat]).zfill(3) + " " + \
                                   str(dt[exp][dat].attrs["ensemble_size"])
                        att, _, _ = sort_dict({"dataset": att, **metadata_json[dia]})
                        if "ensemble_size" in list(att.keys()):
                            del att["ensemble_size"]
                        # from json import dumps as json__dumps
                        # print(json__dumps(att, indent=4))
                        # stop
                        da = da.drop_attrs()
                        da.attrs.update(att)
                        # save
                        dat_o = "multi" + str("model" if pro == "cmip" else "obs") + str(siz)
                        print(pro, dat_o, str(exp).ljust(22), dur, dia, str(len(list_dat)).rjust(3))
                        file_o = "_".join([pro, dat_o, exp, dur, dia]) + ".nc"
                        to_netcdf(da, os__path__join(data_directory, file_o))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
c03_create_ensemble_mean_time_series_1d_obs_plot()
