# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: Correlations between gradient change and ENSO
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
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import correlation_pearson
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, mean_weighted
from enso_change.wrapper.xarray_tools import select_region, create_land_sea_mask

# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": [],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val", "ave_ts_val_glob"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "scenario26", "scenario45", "scenario85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    # "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    "data_projects": ["cmip"],
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
    "proc_ci": 60.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co":{
        # no resample
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
    # threshold for selection
    "proc_se": 50,  # 100 / 3,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.25, "y": 0.25},
        "panel_01": {"x_delt": 8, "x_size": 5, "y_delt": 8, "y_size": 16},
    },
    # colors
    "fig_colors": {
        "cmip": "peru",
        "cmip5": "green",
        "cmip6": "r",
    },
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": {
        "scenario26": "H",
        "scenario45": "D",
        "scenario85": "s",
    },
    # ticks
    "fig_ticks": {
        "correlation": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "",
                "x_nbr_minor": 0,
                "x_tic": None,
                "y_lim": [-0.6, 0.6],
                "y_nam": "correlation",
                # "y_tic": [round(k / 100, 2) for k in range(-70, 71, 35)],
                "y_tic": [round(k / 100, 2) for k in range(-60, 61, 30)],
                "y_nbr_minor": 2,
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_gengts": "Gen24",
        "ave_ts_val_heedts": "HF21",
        "ave_ts_val_pla1ts": "Pla27$_{8-10}$",
        # "ave_ts_val_pla2ts": "Pla27$_{10-10}$",
        "ave_ts_val_pla2ts": "Pla27",
        "ave_ts_val_pla3ts": "Pla27$_{cmip}$",
        "ave_ts_val_seagts": "Sea22",
        "ave_ts_val_wangts": "Wan24",
        "ave_ts_val_watats": "Wat21",
        "cmip": "CMIP 5&6",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        # figure title
        "figure_title": "Gradient change and ENSO amplitude change",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
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
def e76_gradient_vs_change_plot(
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
        proc_se: float = None,
        proc_th: dict[str: int] = None,
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
    proc_se = none_to_default(proc_se, default["proc_se"])
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
    # values_netcdf = {"diagnostic": {"epoch_length": {"project": {"experiment": {"dataset": DataArray}}}}}
    values_netcdf, metadata_netcdf = reader_netcdf(
        data_diagnostics=data_diagnostics2, data_epoch_lengths=data_epoch_lengths, data_experiments=data_experiments,
        data_observational_datasets=data_observational_datasets, data_projects=data_projects,
        data_region={"X": (20, 380), "Y": (-90, 90)}, smile_minimum_size=smile_minimum_size,
        smile_rejected=smile_rejected, smile_require_picontrol=smile_require_picontrol)
    for dia in list(metadata_netcdf.keys()):
        if "_ts_" in dia and "name_short" in list(metadata_netcdf[dia].keys()):
            metadata_netcdf[dia]["name_short"] = metadata_netcdf[dia]["name_short"].replace("TS", "SST")
    print("reader_netcdf", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # remove datasets
    dia1, dur1, pro1 = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat1, exp1 = "multimodel001", "historical"
    da = values_netcdf[dia1][dur1][pro1][exp1][dat1]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    to_remove = [dict_datasets[k] for k in smile_rejected if k in list(dict_datasets.keys())]
    if len(to_remove) > 0:
        for dia in list(values_netcdf.keys()):
            for dur in list(values_netcdf[dia].keys()):
                for pro in list(values_netcdf[dia][dur].keys()):
                    if pro != "cmip":
                        continue
                    for exp in list(values_netcdf[dia][dur][pro].keys()):
                        for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                            # array
                            da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                            del values_netcdf[dia][dur][pro][exp][dat]
                            # list models by number
                            list_mod = da["dataset"].values.tolist()
                            # remove model numbers
                            list_mod = sorted(list(set(list_mod) - set(to_remove)))
                            # select models
                            da = da.sel(indexers={"dataset": list_mod})
                            # save
                            put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # number of datasets
    n_ensemble = {}
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 4:
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                if pro != "cmip":
                    continue
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat, da in values_netcdf[dia][dur][pro][exp].items():
                        n = da.isel(year=0).notnull().sum().values if exp == "historical" else da["dataset"].size
                        put_in_dict(n_ensemble, n, exp, dat)
    for exp in data_experiments:
        if not check_dict_levels(n_ensemble, exp):
            continue
        for dat in list(n_ensemble[exp].keys()):
            print(exp.ljust(10), dat, "number of datasets: " + str(n_ensemble[exp][dat]).rjust(3))
    #
    # -- Processing
    #
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical"]:
                        continue
                    for dat in list(values_netcdf[dia1][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_netcdf, dia1, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia1, dur2, pro, exp2) or \
                                not check_dict_levels(values_netcdf, dia2, dur1, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_netcdf[dia1][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp2][dat])
                        da3 = copy__deepcopy(values_netcdf[dia2][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp1][dat])
                        else:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp2][dat])
                        # select experiment epoch
                        yy = 1850 + int(dur1.split("_")[0]) if exp1 == "historical" else 2100
                        da1, da3 = da1.sel(indexers={"year": yy}), da3.sel(indexers={"year": yy})
                        # select reference epoch
                        yy = 2006
                        da2, da4 = da2.sel(indexers={"year": yy}), da4.sel(indexers={"year": yy})
                        # change
                        da1 = da1 - da2 if exp1 in ["scenario26", "scenario45", "scenario85"] else da2 - da1
                        if dia1.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # degree of GMST
                        da3 = da3 - da4 if exp1 in ["scenario26", "scenario45", "scenario85"] else da4 - da3
                        # change per CC
                        da1 /= da3
                        # save
                        put_in_dict(mean_change, copy__deepcopy(da1).compute(), dur1, dia1, pro, exp1, dat)
    print("mean_change", list(mean_change.keys()))
    default_region = {
        "east_gen": {"latitude": [-2.5, 2.5], "longitude": [210., 270.], "short_name": "EqEqP$_{Gen22}$"},
        "east_pla1": {"latitude": [-8., 10.], "longitude": [230., 270.], "short_name": "EaEqP$_{Pla27}$"},
        "east_pla2": {"latitude": [-10., 10.], "longitude": [230., 270.], "short_name": "EaEqP$_{Pla27}$"},
        "east_pla3": {"latitude": [-10., 10.], "longitude": [240., 280.], "short_name": "EaEqP$_{Pla27}$"},
        "east_sea": {"latitude": [-3., 3.], "longitude": [170., 270.], "short_name": "EqEqP$_{Sea22}$"},
        "east_wan": {"latitude": [-5., 5.], "longitude": [220., 280.], "short_name": "EaEqP$_{Wan24}$"},
        "east_wat": {"latitude": [-5., 5.], "longitude": [180., 280.], "short_name": "EaEqP$_{Wat21}$"},
        "nort_gen": {"latitude": [5., 10.], "longitude": [150., 270.], "short_name": "NoEqP$_{Gen24}$"},
        "sout_gen": {"latitude": [-10., -5.], "longitude": [160., 220.], "short_name": "SoEqP$_{Gen24}$"},
        "west_hee": {"latitude": [-5., 5.], "longitude": [80., 150.], "short_name": "WeEqP$_{HF21}$"},
        "west_gen": {"latitude": [-2.5, 2.5], "longitude": [120., 180.], "short_name": "WeEqP$_{Gen24}$"},
        "west_pla1": {"latitude": [-4., 4.], "longitude": [140., 230.], "short_name": "EaEqP$_{Pla27}$"},
        "west_pla2": {"latitude": [-4., 4.], "longitude": [135., 230.], "short_name": "EaEqP$_{Pla27}$"},
        "west_pla3": {"latitude": [-4., 4.], "longitude": [180., 230.], "short_name": "EaEqP$_{Pla27}$"},
        "west_sea": {"latitude": [-3., 3.], "longitude": [140., 170.], "short_name": "WeEqP$_{Sea22}$"},
        "west_wan": {"latitude": [-5., 5.], "longitude": [140., 170.], "short_name": "WeEqP$_{Wan21}$"},
        "west_wat": {"latitude": [-5., 5.], "longitude": [110., 180.], "short_name": "WeEqP$_{Wat21}$"},
        "gengts": {"regions": ["west_gen", "nort_gen", "sout_gen", "east_gen"], "short_name": "Gen24"},
        "heedts": {"regions": ["west_hee", "east_wat"], "short_name": "HF21"},
        "pla3pr": {"latitude": [-14., 7.], "longitude": [120., 175.], "short_name": "Pla27$_{west}$"},
        "pla1ts": {"regions": ["west_pla1", "east_pla1"], "short_name": "Pla27$_{8-10}$"},
        "pla2ts": {"regions": ["west_pla2", "east_pla2"], "short_name": "Pla27$_{10-10}$"},
        "pla3ts": {"regions": ["west_pla3", "east_pla3"], "short_name": "Pla27$_{cmip}$"},
        "seagts": {"regions": ["west_sea", "east_sea"], "short_name": "Sea22"},
        "wangts": {"regions": ["west_wan", "east_wan"], "short_name": "Wan24"},
        "watats": {"regions": ["west_wat", "east_wat"], "short_name": "Wat21"},
    }
    list_regions = ["heedts", "watats", "seagts", "wangts", "gengts", "pla2ts"]
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if dia not in ["ave_pr_val", "ave_ts_val"]:
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat])
                        mask = create_land_sea_mask(da)
                        da = da.where(mask == 0)
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            if not (reg[-2:] == "pr" and "_pr_" in dia) and not (reg[-2:] == "ts" and "_ts_" in dia):
                                continue
                            list_r2 = copy__deepcopy(default_region[reg]["regions"]) \
                                if reg in list(default_region.keys()) and \
                                   "regions" in list(default_region[reg].keys()) else [reg]
                            da_o, da_t = None, None
                            for r2 in list_r2:
                                # mask outside region
                                da_t = select_region(da, default_region[r2]["latitude"],
                                                     default_region[r2]["longitude"])
                                # horizontal average
                                da_t = mean_weighted(da_t, dim=[dim_y, dim_x], skipna=True)
                                # sum
                                if da_o is None:
                                    da_o = copy__deepcopy(da_t)
                                elif reg == "geng_ts" and r2 in ["west_gen", "nort_gen", "sout_gen"]:
                                    da_o += da_t
                                    if r2 == "sout_gen":
                                        da_o /= 3
                                else:
                                    da_o -= da_t
                            if isinstance(da_o, array_wrapper) and isinstance(da_t, array_wrapper):
                                da_o.attrs.update(da_t.attrs)
                            # save
                            put_in_dict(mean_change_ave, da_o.compute(),
                                        dur, str(dia) + "_" + str(reg), pro, exp, dat)
            if check_dict_levels(mean_change_ave, dur):
                for dia_o in list(mean_change_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1]
                        val = str(default_region[reg]["short_name"]) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    print("mean_change_ave", list(mean_change_ave.keys()))
    # -- correlations
    correlations = {}
    dia_y = "std_ts_ano_n34e"  # "ave_ts_val_pla2ts"  #
    for dur in list(mean_change_ave.keys()):
        for dia_x in list(mean_change_ave[dur].keys()):
            for pro1 in list(mean_change_ave[dur][dia_x].keys()):
                for exp in list(mean_change_ave[dur][dia_x][pro1].keys()):
                    for dat in list(mean_change_ave[dur][dia_x][pro1][exp].keys()):
                        if not check_dict_levels(mean_change, dur, dia_y, pro1, exp1, dat):
                            continue
                        # array
                        da_x = copy__deepcopy(mean_change_ave[dur][dia_x][pro1][exp][dat])
                        da_y = copy__deepcopy(mean_change[dur][dia_y][pro1][exp][dat])
                        # if not check_dict_levels(mean_change_ave, dur, dia_y, pro1, exp1, dat):
                        #     continue
                        # # array
                        # da_x = copy__deepcopy(mean_change_ave[dur][dia_x][pro1][exp][dat])
                        # da_y = copy__deepcopy(mean_change_ave[dur][dia_y][pro1][exp][dat])
                        # correlation CMIP5 & CMIP6
                        corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
                        # save
                        put_in_dict(correlations, corr,dur, dia_y, dia_x, pro1, exp)
                        # separate CMIP5 and CMIP6 datasets
                        list_mod = da_x["dataset"].values.tolist()
                        for pro2 in ["cmip5", "cmip6"]:
                            list_t = sorted([k for k in list_mod
                                             if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)])
                            # select models
                            da_xt = da_x.sel(indexers={"dataset": list_t})
                            da_yt = da_y.sel(indexers={"dataset": list_t})
                            # correlation CMIP5 or CMIP6
                            corr = float(correlation_pearson(da_xt, da_yt, dim="dataset").values)
                            if pro2 == "cmip6" and dia_x == "ave_ts_val_pla2ts":
                                print(dur, dia_y, dia_x, pro2, exp, da_xt.shape, corr)
                            # save
                            put_in_dict(correlations, corr, dur, dia_y, dia_x, pro2, exp)
    print("box ave min-max")
    for dur in list(correlations.keys()):
        t1 = []
        for dia_y in list(correlations[dur].keys()):
            for dia_x in list(correlations[dur][dia_y].keys()):
                for pro in list(correlations[dur][dia_y][dia_x].keys()):
                    t2 = []
                    for exp in list(correlations[dur][dia_y][dia_x][pro].keys()):
                        t2 += [correlations[dur][dia_y][dia_x][pro][exp]]
                    if len(t2) > 0:
                    # if pro == "cmip" and len(t2) > 0:
                        print(dur, dia_y, dia_x, pro.rjust(5), "{0:.2f}".format(round(sum(t2) / len(t2), 2)).rjust(8))
                    t1 += t2
        if len(t1) > 0:
            print(dur, "{0:.2f}".format(round(min(t1), 2)).rjust(8), "{0:.2f}".format(round(max(t1), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: scatter
    group, panel = "goup_01", "panel_01"
    print("scatter", group, panel)
    for dur in list(correlations.keys()):
        list_dia_y = sorted(list(correlations[dur].keys()), key=lambda s: s.lower())
        for ii, dia_y in enumerate(list_dia_y):
            pan = copy__deepcopy(panel)
            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
            tuple_k = (dur, grp, pan)
            print(dur, dia_y, tuple_k)
            # dictionary
            d1 = correlations[dur][dia_y]
            list_d = list(d1.keys())
            nn = len(list_d)
            # panel size
            tk = tuple_k[:-2]
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                for k1, d1 in fig_panel_size.items():
                    if k1 == panel and not check_dict_levels(figure_psize, *tk + (pan,)):
                        for k2, d2 in d1.items():
                            if check_dict_levels(figure_psize, *tk + (pan, k2,)):
                                continue
                            tmp = copy__deepcopy(d2)
                            if k2 == "x_size":
                                tmp *= nn
                            put_in_dict(figure_psize, tmp, *tk + (pan, k2,))
                    elif k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, copy__deepcopy(d1), *tk + (k1,))
            if not check_dict_levels(figure_psize, *tk + (pan,)):
                put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                        put_in_dict(figure_psize, k2, *tk + (k1,))
            # panel axes definition
            if check_dict_levels(fig_ticks, "correlation", panel):
                for k1, k2 in fig_ticks["correlation"][panel].items():
                    if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                        continue
                    tmp = copy__deepcopy(k2)
                    if tmp is None and k1 == "x_lim":
                        tmp = [-0.5, nn - 0.5]
                    elif tmp is None and k1 == "x_tic":
                        tmp = list(range(nn))
                    put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                    if isinstance(tmp, list) and k1 in ["x_tic"]:
                        put_in_dict(figure_axes, [fig_titles[k] for k in list_d], *tuple_k + ("x_lab",))
            dt = {}
            if check_dict_levels(figure_axes, *tuple_k):
                dt = read_from_dict(figure_axes, *tuple_k)
            # test possibility to plot text / legend
            do_legend = False
            if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                    dt["y_lim"] is not None:
                do_legend = True
            # -- title
            if len(list_dia_y) > 1:
                put_in_dict(plot_data, fig_titles[dia_y], *tuple_k + ("title_col",))
            # -- marker
            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
            leg_c, leg_m = [], []
            for kk, dia_x in enumerate(list(correlations[dur][dia_y].keys())):
                for pro in list(correlations[dur][dia_y][dia_x].keys()):
                    for exp in list(correlations[dur][dia_y][dia_x][pro].keys()):
                        mm = fig_markers[exp] if exp in list(fig_markers.keys()) else "o"
                        lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                        lm += [mm]
                        ls += [105 if mm == "H" else 100]
                        lw += [0.5]
                        lx += [kk + (0 if pro == "cmip" else (-0.2 if pro == "cmip5" else 0.2))]
                        ly += [correlations[dur][dia_y][dia_x][pro][exp]]
                        lz += [4 if mm == "H" else (3 if mm == "D" else 2)]
                        leg_c += [pro]
                        leg_m += [exp]
            leg_c = sorted(list(set(leg_c)), key=lambda s: s.lower())
            leg_m = sorted(list(set(leg_m)), key=lambda s: s.lower())
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # legend
                for kk, tt in enumerate(leg_m):
                    mm = fig_markers[tt] if tt in list(fig_markers.keys()) else "o"
                    lc += ["k"]
                    lm += [mm]
                    ls += [84 if mm == "H" else 80]
                    lw += [0.5]
                    lx += [x1 + (kk * 50 + 12) * dx]
                    ly += [y2 - 5 * dy]
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
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                # x = 0 and y = 0
                for n1, k1 in zip(["y"], [dt["y_lim"]]):
                    if k1[0] < 0 < k1[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                        lz += [1]
            plot_type = "cur"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- text
            lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_psize, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # legend
                for kk, tt in enumerate(leg_m):
                    lc += ["k"]
                    lf += [15]
                    lh += ["left"]
                    lr += [0]
                    lt += [fig_titles[tt] if tt in list(fig_titles.keys()) else tt]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x1 + (kk * 50 + 15) * dx]
                    ly += [y2 - 5 * dy]
                for kk, tt in enumerate(leg_c):
                    lc += [fig_colors[tt] if tt in list(fig_colors.keys()) else "cyan"]
                    lf += [15]
                    lh += ["left"]
                    lr += [0]
                    lt += [fig_titles[tt] if tt in list(fig_titles.keys()) else tt]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x1 + (kk * 50 + 15) * dx]
                    ly += [y2 - (5 + 8) * dy]
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
    a1, a4 = str(smile_minimum_size).zfill(2), str(len(smile_rejected)).zfill(2)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(a1) + "mem", str(a4) + "rem"])
    # figure loop
    for k in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[k].keys())
        numb = False if len(list_names) < 2 else True
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", k)))
        # plot
        fig_basic(plot_data[k], list_names, 1, figure_axes[k], fig_format, figure_name_and_path,
                  figure_psize[k], fig_title_bool=numb, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e76_gradient_vs_change_plot()
