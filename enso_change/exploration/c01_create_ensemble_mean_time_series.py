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
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_pr_val_glob",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw",
    #                       "ave_pr_val_peso", "ave_pr_val_pnee", "ave_pr_val_pnew", "ave_pr_val_pnqe",
    #                       "ave_pr_val_poee", "ave_pr_val_poeq",
    #                       "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw", "ave_pr_val_ptse",
    #                       "ave_pr_val_tron", "ave_pr_val_trop"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_ts_val_glob", "ave_ts_val_gloo",
    #                       "ave_ts_val_aequ", "ave_ts_val_atno", "ave_ts_val_atro",
    #                       "ave_ts_val_iequ", "ave_ts_val_isee", "ave_ts_val_itro", "ave_ts_val_iweq",
    #                       "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl",
    #                       "ave_ts_val_gxpn", "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
    #                       "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi",
    #                       "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
    #                       "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
    #                       "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s",
    #                       "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw",
    #                       "ave_ts_val_peso", "ave_ts_val_plee", "ave_ts_val_plew",
    #                       "ave_ts_val_pnee", "ave_ts_val_pnew", "ave_ts_val_pnqe", "ave_ts_val_poee",
    #                       "ave_ts_val_poeq", "ave_ts_val_ptre", "ave_ts_val_ptro", "ave_ts_val_ptrw",
    #                       "ave_ts_val_ptse", "ave_ts_val_tron", "ave_ts_val_trop"],
    # "data_diagnostics1": ["ave_pr_val_glob",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw",
    #                       "ave_pr_val_peso", "ave_pr_val_pnee", "ave_pr_val_pnew", "ave_pr_val_pnqe",
    #                       "ave_pr_val_poee", "ave_pr_val_poeq",
    #                       "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw", "ave_pr_val_ptse",
    #                       "ave_pr_val_tron", "ave_pr_val_trop",
    #                       "ave_ts_val_glob", "ave_ts_val_gloo",
    #                       "ave_ts_val_aequ", "ave_ts_val_atno", "ave_ts_val_atro",
    #                       "ave_ts_val_iequ", "ave_ts_val_isee", "ave_ts_val_itro", "ave_ts_val_iweq",
    #                       "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl",
    #                       "ave_ts_val_gxpn", "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
    #                       "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi",
    #                       "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
    #                       "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
    #                       "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s",
    #                       "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw",
    #                       "ave_ts_val_peso", "ave_ts_val_plee", "ave_ts_val_plew",
    #                       "ave_ts_val_pnee", "ave_ts_val_pnew", "ave_ts_val_pnqe", "ave_ts_val_poee",
    #                       "ave_ts_val_poeq", "ave_ts_val_ptre", "ave_ts_val_ptro", "ave_ts_val_ptrw",
    #                       "ave_ts_val_ptse", "ave_ts_val_tron", "ave_ts_val_trop",
    #                       "std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e",
    #                       "std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e"],
    # "data_diagnostics1": ["std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e",
    #                       "std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e"],
    # "data_diagnostics1": ["ave_pr_val_glob", "ave_pr_val_tron", "ave_ts_val_glob", "ave_ts_val_tron"],
    # "data_diagnostics1": ["ave_pr_val_glob", "ave_pr_val_tron",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw",
    #                       "ave_pr_val_peso",
    #                       "ave_ts_val_glob", "ave_ts_val_tron",
    #                       "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
    #                       "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
    #                       "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s",
    #                       "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw",
    #                       "ave_ts_val_peso",
    #                       "std_pr_ano_n30e", "std_pr_ano_n34e", "std_pr_ano_n40e",
    #                       "std_ts_ano_n30e", "std_ts_ano_n34e", "std_ts_ano_n40e"],
    "data_diagnostics1": ["ave_ts_val_glob", "ave_ts_val_tron", "std_ts_ano_n30e", "std_ts_ano_n34e"],
    # "data_diagnostics1": ["fbk_ts_n30e_to_hf_n30e", "fbk_ts_n34e_to_hf_n34e", "fbk_ts_n34e_to_tx_n40e",
    #                       "ske_ts_ano_n30e", "ske_ts_ano_n34e"],
    "data_diagnostics2": [],
    # "data_diagnostics2": ["ave_pr_val", "ave_sst_val"],
    # list of epoch lengths: list[str]
    # "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    "data_epoch_lengths": ["030_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    # "data_projects": ["observations"],
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
def c01_create_ensemble_mean_time_series_plot(
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
                    if "piControl" not in list_exp or "historical" not in list_exp:
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
            # if exp == "historical":
            #     print(list_ensembles)
    list_datasets, list_numbers = [], []
    for pro in data_projects:
        if not check_dict_levels(n_ensemble, pro):
            continue
        # list datasets for given project
        list_dat = sorted(list(set([dat for exp in list(n_ensemble[pro].keys()) for dat in n_ensemble[pro][exp]])),
                          key=lambda s: s.lower())
        list_datasets += copy__deepcopy(list_dat)
        list_numbers += [k + 1 + (100 if pro == "cmip6" else 0) for k in range(len(list_dat))]
    dict_datasets = {k: v for k, v in zip(list_datasets, list_numbers)}
    # check ensemble size of datasets
    # print("check ensemble size of datasets")
    # for dur in data_epoch_lengths:
    #     for pro in data_projects:
    #         for exp in data_experiments:
    #             # list datasets
    #             list_datasets = []
    #             for dia in data_diagnostics1:
    #                 if check_dict_levels(values_json, dia, dur, pro, exp):
    #                     list_datasets += list(values_json[dia][dur][pro][exp].keys())
    #             list_datasets = sorted(list(set(list_datasets)), key=lambda s: s.lower())
    #             for dat in list_datasets:
    #                 # list diagnostics available for given dataset
    #                 list_diagnostics = []
    #                 for dia in data_diagnostics1:
    #                     if check_dict_levels(values_json, dia, dur, pro, exp, dat):
    #                         list_diagnostics.append(dia)
    #                 missing = sorted(list(set(data_diagnostics1) - set(list_diagnostics)), key=lambda s: s.lower())
    #                 if len(missing) != 0:
    #                     error = "missing diagnostics\n" + ", ".join([dur, pro, exp.ljust(10), dat]) + "\nmissing: " + \
    #                             ", ".join(missing)
    #                     print_fail(inspect__stack(), error)
    #                 # ensemble size
    #                 dt = {dia: values_json[dia][dur][pro][exp][dat]["member"].size for dia in data_diagnostics1}
    #                 maximum_ensemble_size = max(list(dt.values()))
    #                 missing = [str(dia) + " (" + str(dt[dia]) + ")" for dia in data_diagnostics1
    #                            if (exp != "piControl" and dt[dia] != maximum_ensemble_size) or
    #                            (exp == "piControl" and dt[dia] != maximum_ensemble_size and dt[dia] < 40)]
    #                 if len(missing) != 0:
    #                     error = "missing members\n" + ", ".join([dur, pro, exp.ljust(10), dat]) + "\nmissing: " + \
    #                             ", ".join(missing) + "\nmaximum size: " + str(maximum_ensemble_size)
    #                     print_fail(inspect__stack(), error)
    # stop
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
                        # array
                        da1 = values_json[dia][dur][pro][exp][dat]
                        da2 = mean_json[dia][dur][pro][exp][dat]
                        ensemble_size = da1["member"].size
                        da2.attrs.update({"ensemble_size": str(dat) + "=" + str(ensemble_size)})
                        # save
                        exp_o = "scenario" + str(exp[-2:]) if exp[:3] in ["rcp", "ssp"] else copy__deepcopy(exp)
                        # for size in [1, 2, 3, 4, 5, 10]:
                        # for size in [1, 2, 3, 4, 5]:
                        for size in [1]:
                            if ensemble_size >= size:
                                size_o = str(size).zfill(3)
                                put_in_dict(mean_json_exp1, da2, dia, dur, size_o, exp_o, dat)
    print("regroup experiments: done", len(list(mean_json_exp1.keys())))
    # select datasets for piControl and another experiment
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/data"
    print("select datasets for piControl and another experiment")
    for dia in list(mean_json_exp1.keys()):
        for dur in list(mean_json_exp1[dia].keys()):
            for siz in list(mean_json_exp1[dia][dur].keys()):
                print(dia, dur, siz)
                # -- create multiple piControl experiments corresponding to models in each experiment
                # list experiments
                list_experiments = [k for k in list(mean_json_exp1[dia][dur][siz].keys()) if k != "piControl"]
                dt = {}
                for exp in list_experiments:
                    for dat in list(mean_json_exp1[dia][dur][siz][exp].keys()):
                        if not check_dict_levels(mean_json_exp1, dia, dur, siz, "piControl", dat):
                            print(dia, dur, siz, exp.ljust(10), dat.ljust(17), "missing piControl")
                            continue
                        # arrays
                        da1 = mean_json_exp1[dia][dur][siz][exp][dat]
                        da2 = mean_json_exp1[dia][dur][siz]["piControl"][dat]
                        # save
                        put_in_dict(dt, da1, exp, dat)
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
                        att += "; " + str(dict_datasets[dat]).zfill(3) + " " + str(dt[exp][dat].attrs["ensemble_size"])
                    att, _, _ = sort_dict({"dataset": att, **metadata_json[dia]})
                    if "ensemble_size" in list(att.keys()):
                        del att["ensemble_size"]
                    da = da.drop_attrs()
                    da.attrs.update(att)
                    # save
                    print("cmip", "multimodel" + str(siz), str(exp).ljust(22), dur, dia, str(len(list_dat)).rjust(3))
                    file_o = "_".join(["cmip", "multimodel" + str(siz), exp, dur, dia]) + ".nc"
                    to_netcdf(da, os__path__join(figure_directory, file_o))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
c01_create_ensemble_mean_time_series_plot()
