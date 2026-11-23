# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: scatterplots of PR, rSST, SST, rTS, TS biases (values at the end of the 20th century) or changes  (departure
# from the end of the 20th century) vs. ENSO changes (departure from the end of the 20th century)
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from datetime import datetime, timedelta
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
from enso_change.plots.panels import default_plot, default_region
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, correlation_pearson, linear_regression, percentile
from enso_change.wrapper.xarray_base import cf_dim_to_dim, mean_weighted
from enso_change.wrapper.xarray_tools import create_land_sea_mask, select_region
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
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_glob", "ave_ts_val_tron",
    #                       "ave_sst_val", "ave_ts_val"],
    # "data_diagnostics2": ["std_pr_ano_n34e", "std_ts_ano_n34e", "ave_ts_val_glob", "ave_ts_val_tron",
    #                       "ave_pr_val", "ave_ts_val"],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_ts_val_glob", "ave_pr_val"],
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_glob"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["050_year_epoch"],
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
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 90.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 30,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.5, "y": 0.5},
        "panel_01": {"x_delt": 3, "x_size": 8, "y_delt": 3, "y_size": 8},
    },
    # colors
    "fig_colors": {
        "cmip5": "green",
        "cmip6": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "P",
        "cmip6": "D",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "bias--pcc": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "VALUE (UNITS)",
                "x_tic": None,
            },
        },
        "bias--raw": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "VALUE (UNITS)",
                "x_tic": None,
            },
        },
        "change--pcc": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "change VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                "x_tic": None,
            },
        },
        "change--raw": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "change VALUE (UNITS)",
                "x_tic": None,
            },
        },
        "std_pr_ano_n34e--pcc": {
            "panel_01": {
                # "y_lim": [-100, 100],
                # "y_nam": "VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                # "y_tic": list(range(-100, 101, 50)),
                "y_lim": None,
                "y_nam": "change VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                "y_tic": None,
            },
        },
        "std_pr_ano_n34e--raw": {
            "panel_01": {
                # "y_lim": [-100, 100],
                # "y_nam": "VALUE (UNITS)",
                # "y_tic": list(range(-100, 101, 50)),
                "y_lim": None,
                "y_nam": "change VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "std_ts_ano_n34e--pcc": {
            "panel_01": {
                # "y_lim": [-100, 100],
                # "y_nam": "VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                # "y_tic": list(range(-100, 101, 50)),
                "y_lim": None,
                "y_nam": "change VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                "y_tic": None,
            },
        },
        "std_ts_ano_n34e--raw": {
            "panel_01": {
                # "y_lim": [-100, 100],
                # "y_nam": "VALUE (UNITS)",
                # "y_tic": list(range(-100, 101, 50)),
                "y_lim": None,
                "y_nam": "change VALUE (UNITS)",
                "y_tic": None,
            },
        },
    },
    # titles
    "fig_titles": {
        "pcc": "change per GMST",
        "raw": "change",
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
def e50_change_explained_bias_change_scatter_plot(
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
    # values_netcdf = {"diagnostic": {"epoch_length": {"project": {"experiment": {"dataset": DataArray}}}}}
    values_netcdf, metadata_netcdf = reader_netcdf(
        data_diagnostics=data_diagnostics2, data_epoch_lengths=data_epoch_lengths, data_experiments=data_experiments,
        data_observational_datasets=data_observational_datasets, data_projects=data_projects,
        data_region={"X": (20, 380), "Y": (-90, 90)}, smile_minimum_size=smile_minimum_size,
        smile_rejected=smile_rejected, smile_require_picontrol=smile_require_picontrol)
    print("reader_netcdf", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # number of datasets
    n_ensemble = {}
    for dia in list(values_netcdf.keys()):
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                if pro == "observations":
                    continue
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        if dat[-3:] in ["003", "004", "005"]:
                            del values_netcdf[dia][dur][pro][exp][dat]
                            continue
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
    dia2, dur2, pro2 = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat2, exp2 = "multimodel001", "historical"
    da = values_netcdf[dia2][dur2][pro2][exp2][dat2]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    # # -- compute rSST
    dia2 = "ave_ts_val_tron"
    for dia1 in ["ave_sst_val", "ave_ts_val"]:
        if dia1 not in list(values_netcdf.keys()):
            continue
        dia_o = dia1.replace("_sst_", "_rsst_").replace("_ts_", "_rts_")
        for dur in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur].keys()):
                for exp in list(values_netcdf[dia1][dur][pro].keys()):
                    for dat in list(values_netcdf[dia1][dur][pro][exp].keys()):
                        if not check_dict_levels(values_netcdf, dia2, dur, pro, exp, dat):
                            continue
                        # arrays
                        da1 = copy__deepcopy(values_netcdf[dia1][dur][pro][exp][dat])
                        da2 = copy__deepcopy(values_netcdf[dia2][dur][pro][exp][dat])
                        # select common years
                        list_y = sorted(list(set(da1["year"].values.tolist()) & set(da2["year"].values.tolist())))
                        da1, da2 = da1.sel(indexers={"year": list_y}), da2.sel(indexers={"year": list_y})
                        # remove tropical mean SST
                        da1 -= da2
                        # save
                        put_in_dict(values_netcdf, da1, dia_o, dur, pro, exp, dat)
        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia1])
        val = metadata_netcdf[dia1]["name_short"].replace("SST", "rSST").replace("ts", "rTS")
        metadata_netcdf[dia_o]["name_short"] = val
        # del values_netcdf[dia1]
    for dia1 in [dia2]:  # [dia2, dia3]:  #
        if dia1 in list(values_netcdf.keys()):
            del values_netcdf[dia1]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- select bias year
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_bias = {}
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 3:
            continue
        for dur1 in list(values_netcdf[dia].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_netcdf[dia][dur1].keys()):
                for exp in list(values_netcdf[dia][dur1][pro].keys()):
                    if exp in ["historical", "piControl"]:
                        continue
                    for dat in list(values_netcdf[dia][dur1][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur1][pro][exp][dat])
                        # select year
                        yy = int(dur1.split("_")[0]) if exp == "piControl" else 2006
                        da = da.sel(indexers={"year": yy})
                        # save
                        put_in_dict(mean_bias, da, dur1, dia, pro, exp, dat, "pcc")
                        put_in_dict(mean_bias, da, dur1, dia, pro, exp, dat, "raw")
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 in ["historical", "piControl"]:
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
                        # change per CC or not
                        for war in ["raw", "pcc"]:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if war == "pcc" and exp1 in ["scenario26", "scenario45", "scenario85"]:
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- create land-sea mask
    dia2, pro2 = [k for k in data_diagnostics2 if len(k.split("_")) == 3][0], data_projects[0]
    exp2 = list(mean_bias[dur2][dia2][pro2].keys())[0]
    dat2 = list(mean_bias[dur2][dia2][pro2][exp2].keys())[0]
    war2 = list(mean_bias[dur2][dia2][pro2][exp2][dat2].keys())[0]
    mask = create_land_sea_mask(mean_bias[dur2][dia2][pro2][exp2][dat2][war2])
    # -- find the best Pacific PR bias region mean
    dia2, dur2, pro2, typ2 = "ave_pr_val", data_epoch_lengths[-1], data_projects[0], "bias"
    dt = {}
    t1 = datetime.now()
    print(t1.strftime("%m/%d %H-%M"))
    for exp in list(mean_bias[dur2][dia2][pro2].keys()):
        if exp == "historical":
            continue
        for dat in list(mean_bias[dur2][dia2][pro2][exp].keys()):
            for war in list(mean_bias[dur2][dia2][pro2][exp][dat].keys()):
                print(exp, dat, war)
                # array
                da1 = copy__deepcopy(mean_bias[dur2][dia2][pro2][exp][dat][war])
                da2 = copy__deepcopy(mean_change[dur2]["std_ts_ano_n34e"][pro2][exp][dat][war])
                # get lat-lon dimensions
                dim_x, dim_y = cf_dim_to_dim(da1, "X"), cf_dim_to_dim(da1, "Y")
                # loop on regions
                t2 = datetime.now()
                for la1 in range(-15, 0, 1):
                    for la2 in range(la1 + 5, 6, 1):
                        for lo1 in range(180, 231, 5):
                            for lo2 in range(lo1 + 50, 281, 5):
                                da_o = copy__deepcopy(da1)
                                # mask outside region
                                da_o = select_region(da_o, [la1, la2], [lo1, lo2])
                                # mask land
                                da_o = da_o.where(mask == 0)
                                # horizontal average
                                da_o = mean_weighted(da_o, dim=[dim_y, dim_x], skipna=True)
                                # correlation with ENSO change
                                corr = float(correlation_pearson(da_o, da2, dim="dataset").values)
                                # save correlation
                                reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2]])
                                put_in_dict(dt, [corr], reg)
                delta_t = datetime.now() - t2
                hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
                mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
                ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
                print(exp, dat, war, str(hh) + ":" + str(mm) + ":" + str(ss))
    delta_t = datetime.now() - t1
    hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    print(str(hh) + ":" + str(mm) + ":" + str(ss))
    # rank correlations
    dt = {k:  abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    # dt = dict(sorted(dt.items(), key=lambda item: item[1], reverse=True))
    dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    print("best regions")
    for k in list(dt.keys()):
        la1, la2, lo1, lo2 = int(k[:4]), int(k[4:8]), int(k[8:12]), int(k[12:])
        print(str(la1).rjust(10) + "N", str(la2).rjust(3) + "N", str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E",
              "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    stop
    # -- compute regional mean
    default_region = {
        # "hsou3050": {"latitude": [-30., -50.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou3060": {"latitude": [-30., -60.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou3070": {"latitude": [-30., -70.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou3080": {"latitude": [-30., -80.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou3090": {"latitude": [-30., -90.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou4050": {"latitude": [-40., -50.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou4060": {"latitude": [-40., -60.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou4070": {"latitude": [-40., -70.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou4080": {"latitude": [-40., -80.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hsou4090": {"latitude": [-40., -90.], "longitude": [20., 380.], "short_name": "HSOU"},
        # "hnor3050": {"latitude": [30., 50.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor3060": {"latitude": [30., 60.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor3070": {"latitude": [30., 70.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor3080": {"latitude": [30., 80.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor3090": {"latitude": [30., 90.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor4050": {"latitude": [40., 50.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor4060": {"latitude": [40., 60.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor4070": {"latitude": [40., 70.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor4080": {"latitude": [40., 80.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "hnor4090": {"latitude": [40., 90.], "longitude": [20., 380.], "short_name": "HNOR"},
        # "equa2020": {"latitude": [-20., 20.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1515": {"latitude": [-15., 15.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1414": {"latitude": [-14., 14.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1313": {"latitude": [-13., 13.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1212": {"latitude": [-12., 12.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1111": {"latitude": [-11., 11.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa1010": {"latitude": [-10., 10.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa0909": {"latitude": [-9., 9.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa0808": {"latitude": [-8., 8.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa0707": {"latitude": [-7., 7.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa0606": {"latitude": [-6., 6.], "longitude": [20., 380.], "short_name": "EQUA"},
        # "equa": {"latitude": [-5., 5.], "longitude": [20., 380.], "short_name": "N3&4"},
        "n20e1515": {"latitude": [-15., 15.], "longitude": [240., 280.], "short_name": "N2"},
        "n20e1010": {"latitude": [-10., 10.], "longitude": [240., 280.], "short_name": "N2"},
        "n20e": {"latitude": [-5., 5.], "longitude": [240., 280.], "short_name": "N2"},
        "n30e0703": {"latitude": [-7., 3.], "longitude": [210., 270.], "short_name": "N3"},
        "n30e1515": {"latitude": [-15., 15.], "longitude": [210., 270.], "short_name": "N3"},
        "n30e1010": {"latitude": [-10., 10.], "longitude": [210., 270.], "short_name": "N3"},
        "n30e": {"latitude": [-5., 5.], "longitude": [210., 270.], "short_name": "N3"},
        "n34e1515": {"latitude": [-15., 15.], "longitude": [190., 240.], "short_name": "N3.4"},
        "n34e1010": {"latitude": [-10., 10.], "longitude": [190., 240.], "short_name": "N3.4"},
        "n34e": {"latitude": [-5., 5.], "longitude": [190., 240.], "short_name": "N3.4"},
        "n36e0703": {"latitude": [-7., 3.], "longitude": [190., 280.], "short_name": "N3.6"},
        "n43e0703": {"latitude": [-7., 3.], "longitude": [160., 270.], "short_name": "N3&4"},
        "n43e1515": {"latitude": [-15., 15.], "longitude": [160., 270.], "short_name": "N3&4"},
        "n43e1010": {"latitude": [-10., 10.], "longitude": [160., 270.], "short_name": "N3&4"},
        "n43e": {"latitude": [-5., 5.], "longitude": [160., 270.], "short_name": "N3&4"},
        "n50e1515": {"latitude": [-15., 15.], "longitude": [120., 160.], "short_name": "N5"},
        "n50e1010": {"latitude": [-10., 10.], "longitude": [120., 160.], "short_name": "N5"},
        "n50e": {"latitude": [-5., 5.], "longitude": [120., 160.], "short_name": "N5"},
        "pceq": {"latitude": [-5., 5.], "longitude": [160., 230.], "short_name": "PCEQ"},
        # "pbhw": {"latitude": [0., -25., -25., -10., -5., 5., 5., 15., 15., 0.],
        #          "longitude": [130., 160., 220., 220., 170., 170., 210., 210., 130., 130.], "short_name": "PBHW"},
        # "prae": {"latitude": [-10., -10., 0., 0., 12., 12., -10.],
        #          "longitude": [325., 355., 355., 335., 330., 290., 325.], "short_name": "AWpr"},
        # "prie": {"latitude": [-15., 0.], "longitude": [50., 100.], "short_name": "ISpr"},
        # "prpe": {"latitude": [-8., 8.], "longitude": [160., 270.], "short_name": "N3&4pr"},
        # "prps": {"latitude": [-30., -15.], "longitude": [170., 220.], "short_name": "N3&4pr"},
    }
    list_regions = list(default_region.keys())
    l1, l2, l3 = ["n30e", "n20e"], ["n50e", "pceq"], ["pgx1", "pgx2"]
    bias_change_ave = {}
    for typ, d1 in zip(["bias", "change"], [mean_bias, mean_change]):
        for dur in list(d1.keys()):
            for dia in list(d1[dur].keys()):
                if len(dia.split("_")) != 3:
                    continue
                for pro in list(d1[dur][dia].keys()):
                    for exp in list(d1[dur][dia][pro].keys()):
                        for dat in list(d1[dur][dia][pro][exp].keys()):
                            for war in list(d1[dur][dia][pro][exp][dat].keys()):
                                # array
                                da = copy__deepcopy(d1[dur][dia][pro][exp][dat][war])
                                # get lat-lon dimensions
                                dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                                # loop on regions
                                for reg in list_regions:
                                    da_o = copy__deepcopy(da)
                                    # mask outside region
                                    da_o = select_region(da_o, default_region[reg]["latitude"],
                                                         default_region[reg]["longitude"])
                                    # mask land
                                    if dia == "ave_pr_val":
                                        da_o = da_o.where(mask == 0)
                                    # horizontal average
                                    da_o = mean_weighted(da_o, dim=[dim_y, dim_x], skipna=True)
                                    # save
                                    lk = (dur, typ, pro, dat, str(dia) + "_" + str(reg), war, exp)
                                    put_in_dict(bias_change_ave, da_o, *lk)
                                # regional differences
                                for r1, r2, nn in zip(l1, l2, l3):
                                    if r1 not in list_regions or r2 not in list_regions or dia == "ave_rsst_val" or \
                                            dia == "ave_rts_val":
                                        continue
                                    da1 = read_from_dict(bias_change_ave,
                                                         dur, typ, pro, dat, str(dia) + "_" + str(r1), war, exp)
                                    da2 = read_from_dict(bias_change_ave,
                                                         dur, typ, pro, dat, str(dia) + "_" + str(r2), war, exp)
                                    lk = (dur, typ, pro, dat, str(dia) + "_" + str(nn), war, exp)
                                    put_in_dict(bias_change_ave, da2 - da1, *lk)
                            if check_dict_levels(bias_change_ave, dur, typ, pro, dat):
                                for dia_o in list(bias_change_ave[dur][typ][pro][dat].keys()):
                                    if not check_dict_levels(metadata_netcdf, dia_o):
                                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                                        r1 = dia_o.split("_")[-1]
                                        reg = r1[:4].replace("asou", "ATL$_S$").replace("hsou", "SH")
                                        reg = reg.replace("anor", "ATL$_N$").replace("hnor", "NH")
                                        reg = reg.replace("n20e", "N2")
                                        reg = reg.replace("n30e", "N3").replace("n34e", "N3.4")
                                        reg = reg.replace("n43e", "N3&4").replace("n36e", "N3&3.4")
                                        reg = reg.replace("n50e", "N5").replace("pgx1", "$\Delta_{W-E}$")
                                        reg = reg.replace("pgx2", "$\Delta_{C-E}$")
                                        reg = reg.replace("pceq", "CEP")
                                        reg = reg.replace("pbhw", "PWHS").replace("prae", "ATL$_{WE}$")
                                        reg = reg.replace("prie", "IND$_E$").replace("prpe", "EP$_B$")
                                        reg = reg.replace("prps", "PAC$_{WCS}$")
                                        if len(r1) == 8:
                                            reg += "[" + str(int(r1[4:6]))
                                            if r1[4:6] != r1[6:8] and r1[:2] in ["n2", "n3", "n4", "n5"]:
                                                reg += "S-" + str(int(r1[6:8])) + "N]"
                                            elif r1[4:6] != r1[6:8]:
                                                reg += "-" + str(int(r1[6:8]))
                                                reg += "S]" if "sou" in r1 else "N]"
                                            else:
                                                reg += "S/N]"
                                        val = str(reg) + " " + metadata_netcdf[dia_o]["name_short"]
                                        metadata_netcdf[dia_o]["name_short"] = val
    print("regional mean", list(bias_change_ave.keys()))
    # # print worst bias
    # dia2, dat2, dur2, exp2 = "ave_pr_val_n36e0703", "multimodel001", data_epoch_lengths[-1], data_experiments[0]
    # pro2, typ2, war2 = data_projects[0], "bias", "raw"
    # dt2 = {v: k for k, v in dict_datasets.items()}
    # for thr in [5, 10, 15]:
    #     da = copy__deepcopy(bias_change_ave[dur2][typ2][pro2][dat2][dia2][war2][exp2]).compute()
    #     t1 = percentile(da, per=100 - thr)
    #     l1 = da.where(da >= t1, drop=True)["dataset"].values.tolist()
    #     l1, l2 = [k for k in l1 if k < 100], [k for k in l1 if k > 100]
    #     print(thr, exp2, da["dataset"].size)
    #     print("cmip5".rjust(10), ", ".join([dt2[k] for k in l1]))
    #     print("cmip6".rjust(10), ", ".join([dt2[k] for k in l2]))
    # stop
    # -- order regions
    dia2 = data_diagnostics2[0]
    # bias_change_ave_sort = {}
    # for dur in list(bias_change_ave.keys()):
    #     for typ in list(bias_change_ave[dur].keys()):
    #         for pro1 in list(bias_change_ave[dur][typ].keys()):
    #             # list diagnostics
    #             d1 = bias_change_ave[dur][typ][pro1]
    #             list_dia = sorted(list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())])),
    #                               key=lambda s: s.lower())
    #             dt1 = {}
    #             for dia in list_dia:
    #                 for dat in list(d1.keys()):
    #                     if not check_dict_levels(bias_change_ave, dur, typ, pro1, dat, dia):
    #                         continue
    #                     for war in list(d1[dat][dia].keys()):
    #                         for exp in list(d1[dat][dia][war].keys()):
    #                             if not check_dict_levels(mean_change, dur, dia2, pro1, exp, dat, war):
    #                                 continue
    #                             # array
    #                             da1 = copy__deepcopy(bias_change_ave[dur][typ][pro1][dat][dia][war][exp])
    #                             da2 = copy__deepcopy(mean_change[dur][dia2][pro1][exp][dat][war]).compute()
    #                             l1 = sorted(list(set(da1["dataset"].values.tolist()) &
    #                                              set(da2["dataset"].values.tolist())))
    #                             # separate cmip5/6
    #                             dt2 = {}
    #                             for pro2 in ["cmip5", "cmip6"]:
    #                                 # model numbers corresponding to cmip5/6
    #                                 l2 = [k for k in l1
    #                                       if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
    #                                 da3 = copy__deepcopy(da2).sel(indexers={"dataset": l2})
    #                                 # find top / bottom of ENSO SSTA change
    #                                 thr = "bottom"
    #                                 t1 = (copy__deepcopy(proc_th) if thr == "bottom" else 100 - proc_th) / 100
    #                                 t2 = da3.quantile(t1)
    #                                 cond = da3 <= t2 if thr == "bottom" else da3 >= t2
    #                                 l2 = da3.where(cond, drop=True)["dataset"].values.tolist()
    #                                 # select models
    #                                 dt2[pro2] = copy__deepcopy(da1).sel(indexers={"dataset": l2})
    #                             # compare cmip5/6
    #                             da11, da22 = dt2["cmip5"], dt2["cmip6"]
    #                             # difference of average
    #                             if da11["dataset"].size >= 4 and da22["dataset"].size >= 4:
    #                                 # previous value
    #                                 da_o = copy__deepcopy(dt1[dia]) if dia in list(dt1.keys()) else 0
    #                                 # normalized difference
    #                                 da_t = da22.mean(dim="dataset") - da11.mean(dim="dataset")
    #                                 da_t = float((da_t / da1.mean(dim="dataset")).values)
    #                                 # save
    #                                 dt1[dia] = da_o + da_t
    #             # sort regions
    #             sorted_regions = sorted(dt1, key=dt1.get)
    #             for dat in list(d1.keys()):
    #                 for dia in sorted_regions:
    #                     for war in list(d1[dat][dia].keys()):
    #                         for exp, da in d1[dat][dia][war].items():
    #                             # save
    #                             put_in_dict(bias_change_ave_sort, da, dur, typ, pro1, dat, dia, war, exp)
    bias_change_ave_sort = {}
    for dur in list(bias_change_ave.keys()):
        for typ in list(bias_change_ave[dur].keys()):
            for pro in list(bias_change_ave[dur][typ].keys()):
                for dia_y in list(mean_change[dur].keys()):
                    if not check_dict_levels(mean_change, dur, dia_y, pro) or \
                            not (len(dia_y.split("_")) == 4 and dia_y.split("_")[0] in ["std", "var"]):
                        continue
                    # list diagnostics
                    d1 = bias_change_ave[dur][typ][pro]
                    list_dia = sorted(list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())])),
                                      key=lambda s: s.lower())
                    dt1 = {}
                    for dia_x in list_dia:
                        for dat in list(d1.keys()):
                            if not check_dict_levels(bias_change_ave, dur, typ, pro, dat, dia_x):
                                continue
                            for war in list(d1[dat][dia_x].keys()):
                                for exp in list(d1[dat][dia_x][war].keys()):
                                    if not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat, war):
                                        continue
                                    # array
                                    da1 = copy__deepcopy(bias_change_ave[dur][typ][pro][dat][dia_x][war][exp])
                                    da2 = copy__deepcopy(mean_change[dur][dia_y][pro][exp][dat][war])
                                    # linearregression
                                    ln = linear_regression(da1.values, da2.values)
                                    # save
                                    put_in_dict(dt1, [float(ln["r"])], dia_x)
                    # average correlation
                    dt1 = {k: abs(float(numpy__array(v).mean())) for k, v in dt1.items()}
                    # sort regions
                    sorted_regions = sorted(dt1, key=dt1.get, reverse=True)
                    for dat in list(d1.keys()):
                        for dia_x in sorted_regions:
                            for war in list(d1[dat][dia_x].keys()):
                                for exp, da in d1[dat][dia_x][war].items():
                                    # save
                                    put_in_dict(bias_change_ave_sort, da,
                                                dur, typ, pro, dat, dia_y, dia_x, war, exp)
    print("order regional mean", list(bias_change_ave_sort.keys()))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_01"
    for dur in list(bias_change_ave_sort.keys()):
        if not check_dict_levels(mean_change, dur):
            continue
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for typ in list(bias_change_ave_sort[dur].keys()):
            for pro in list(bias_change_ave_sort[dur][typ].keys()):
                for dat in list(bias_change_ave_sort[dur][typ][pro].keys()):
                    for dia_y in list(bias_change_ave_sort[dur][typ][pro][dat].keys()):
                        if not check_dict_levels(mean_change, dur, dia_y, pro) or \
                                not (len(dia_y.split("_")) == 4 and dia_y.split("_")[0] in ["std", "var"]):
                            continue
                        for ii, dia_x in enumerate(list(bias_change_ave_sort[dur][typ][pro][dat][dia_y].keys())):
                            list_war  = sorted(list(bias_change_ave_sort[dur][typ][pro][dat][dia_y][dia_x].keys()),
                                               key=lambda s: s.lower())
                            for jj, war in enumerate(list_war):
                                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                                tuple_k = (dur, typ, dat, dia_y, grp, pan)
                                print(dur, typ, pro, dat, dia_y, dia_x.ljust(21), war, tuple_k)
                                # dictionaries
                                dix = copy__deepcopy(bias_change_ave_sort[dur][typ][pro][dat][dia_y][dia_x][war])
                                diy = copy__deepcopy(mean_change[dur][dia_y][pro])
                                print("diy", list(diy.keys()))
                                x_mm = [float(da.min().values) for da in dix.values()]
                                x_mm += [float(da.max().values) for da in dix.values()]
                                x1, x2 = min(x_mm), max(x_mm)
                                x_mm = [math__floor(x1), math__ceil(x2)]
                                if x2 - x1 < 5:
                                    n1 = 3 if x2 - x1 < 0.05 else (2 if x2 - x1 < 0.5 else 1)
                                    n2 = 10**n1
                                    x_mm = [round(math__floor(x1 * n2) / n2, n1), round(math__ceil(x2 * n2) / n2, n1)]
                                y_mm = [float(diy[exp][dat][war].min().values) for exp in list(dix.keys())]
                                y_mm += [float(diy[exp][dat][war].max().values) for exp in list(dix.keys())]
                                y1, y2 = min(y_mm), max(y_mm)
                                y_mm = [math__floor(y1), math__ceil(y2)]
                                if y2 - y1 < 5:
                                    n1 = 3 if y2 - y1 < 0.05 else (2 if y2 - y1 < 0.5 else 1)
                                    n2 = 10 ** n1
                                    y_mm = [round(math__floor(y1 * n2) / n2, n1), round(math__ceil(y2 * n2) / n2, n1)]
                                # panel size
                                tk = tuple_k[:-2]
                                if not check_dict_levels(figure_ps, *tk + (pan,)):
                                    put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                                    for k1, k2 in fig_panel_size.items():
                                        if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                            put_in_dict(figure_ps, k2, *tk + (k1,))
                                # panel axes definition
                                for dd1 in [dia_x, dia_y]:
                                    dd2 = str(dd1 if dd1 == dia_y else typ) + "--" + str(war)
                                    if check_dict_levels(fig_ticks, dd2, panel):
                                        for k1, k2 in fig_ticks[dd2][panel].items():
                                            if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                                                continue
                                            tmp = copy__deepcopy(k2)
                                            if isinstance(tmp, str) and k1 in ["x_nam", "y_nam"]:
                                                tmp = tmp.replace("VALUE", metadata_netcdf[dd1]["name_short"])
                                                tmp = tmp.replace("UNITS", metadata_netcdf[dd1]["units"])
                                                tmp = tmp.replace(" ()", "")
                                            elif tmp is None and k1 == "x_lim":
                                                tmp = copy__deepcopy(x_mm)
                                            elif tmp is None and k1 == "x_tic":
                                                tmp = copy__deepcopy(x_mm)
                                            elif tmp is None and k1 == "y_lim":
                                                tmp = copy__deepcopy(y_mm)
                                            elif tmp is None and k1 == "y_tic":
                                                tmp = copy__deepcopy(y_mm)
                                            put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                                dt = {}
                                if check_dict_levels(figure_axes, *tuple_k):
                                    dt = read_from_dict(figure_axes, *tuple_k)
                                # test possibility to plot text / legend
                                do_legend = False
                                if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and \
                                        dt["x_lim"] is not None and dt["y_lim"] is not None:
                                    do_legend = True
                                # -- title
                                if ii == 0:
                                    put_in_dict(plot_data, fig_titles[war], *tuple_k + ("title_col",))
                                if jj == 0:
                                    val = copy__deepcopy(metadata_netcdf[dia_x]["name_short"])
                                    put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                                # -- marker
                                lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                                for exp in list(dix.keys()):
                                    for pro2 in ["cmip5", "cmip6"]:
                                        da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp][dat][war])
                                        # model numbers corresponding to cmip5/6
                                        list_mod = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                               set(da_y["dataset"].values.tolist())))
                                        list_mod = [k for k in list_mod if
                                                    (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                                        da_x = da_x.sel(indexers={"dataset": list_mod}).values.tolist()
                                        da_y = da_y.sel(indexers={"dataset": list_mod}).values.tolist()
                                        nn = len(da_x)
                                        lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"] * nn
                                        lm += [fig_markers[pro2] if pro2 in list(fig_markers.keys()) else "o"] * nn
                                        ls += [80] * nn
                                        lw += [0.5] * nn
                                        lx += copy__deepcopy(da_x)
                                        ly += copy__deepcopy(da_y)
                                        lz += [numpy__random__randint(3, 9) for _ in range(nn)]
                                pt = "mar"
                                if len(lx) > 0 and len(lx) == len(ly):
                                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(pt) + "_c",))
                                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(pt) + "_lw",))
                                    put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(pt) + "_m",))
                                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(pt) + "_s",))
                                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(pt) + "_y",))
                                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(pt) + "_z",))
                                # -- curve
                                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                                # regression
                                for exp in list(dix.keys()):
                                    da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp][dat][war])
                                    ln = linear_regression(da_x.values, da_y.values)
                                    lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                                    ls += ["-"]
                                    lw += [1]
                                    lx += [x_mm]
                                    ly += [[ln["s"] * k + ln["i"] for k in x_mm]]
                                    lz += [20]
                                if isinstance(do_legend, bool) and do_legend:
                                    # x = 0 and y = 0
                                    for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                                        if k1[0] < 0 < k1[1]:
                                            lc += ["k"]
                                            ls += [":"]
                                            lw += [1]
                                            lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                            ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                            lz += [10]
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
                                    df = read_from_dict(figure_ps, *tk)
                                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                    dx = (x2 - x1) * default_plot["size_x"] / (
                                            100 * df[pan]["x_size"] * df["frac"]["x"])
                                    dy = (y2 - y1) * default_plot["size_y"] / (
                                            100 * df[pan]["y_size"] * df["frac"]["y"])
                                    # regression
                                    for kk, exp in enumerate(list(dix.keys())):
                                        da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp][dat][war])
                                        ln = linear_regression(da_x.values, da_y.values)
                                        lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                                        lf += [12]
                                        lh += ["right"]
                                        lr += [0]
                                        val = "<10$^{-3}$" if ln["p"] < 1e-3 else (
                                            "<10$^{-2}$" if ln["p"] < 1e-2 else ("<0.05" if ln["p"] < 0.05 else (
                                                "<0.10" if ln["p"] <= 0.1 else ">0.10")))
                                        lt += [fig_titles[exp] + " r=" + "{0:.2f}".format(round(ln["r"], 2)) + " p" +
                                               str(val)]
                                        lv += ["top" if dia_y == "std_pr_ano_n34e" else "bottom"]
                                        lw += ["regular"]
                                        lx += [x2 - 2 * dx]
                                        ly += [y2 - (kk * 6 + 2) * dy if dia_y == "std_pr_ano_n34e" else
                                               y1 + (kk * 6 + 2) * dy]
                                ptype = "text"
                                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (ptype,))
                                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(ptype) + "_c",))
                                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(ptype) + "_fs",))
                                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(ptype) + "_ha",))
                                    put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(ptype) + "_r",))
                                    put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(ptype) + "_va",))
                                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(ptype) + "_w",))
                                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(ptype) + "_x",))
                                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(ptype) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1 = int(proc_ci) if proc_ci == int(proc_ci) else copy__deepcopy(proc_ci)
    a2 = int(proc_th) if proc_ci == int(proc_th) else copy__deepcopy(proc_th)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_" + str(a1) + "ci_" + str(a2) + "th"
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for typ in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for dat in sorted(list(plot_data[dur][typ].keys()), key=lambda s: s.lower()):
                for dia_y in sorted(list(plot_data[dur][typ][dat].keys()), key=lambda s: s.lower()):
                    # list plot names
                    list_names = list(plot_data[dur][typ][dat][dia_y].keys())
                    # path to output figure
                    figure_name_and_path = os__path__join(figure_directory,
                                                          figure_name + "_".join(("", dur, typ, dat, dia_y)))
                    # plot
                    fig_basic(plot_data[dur][typ][dat][dia_y], list_names, 1,
                              figure_axes[dur][typ][dat][dia_y], fig_format, figure_name_and_path,
                              figure_ps[dur][typ][dat][dia_y], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e50_change_explained_bias_change_scatter_plot()
