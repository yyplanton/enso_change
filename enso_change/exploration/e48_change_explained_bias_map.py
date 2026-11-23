# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Maps of rSST biases (values at the end of the 20th century)
#  - Distributions of ENSO changes (departure from the end of the 20th century) as a function of selected bias
#
# Layout
#
# bias    CMIP5 rSST   CMIP6 rSST   difference
#        ------------ ------------ ------------
#        |          | |          | |          |
#        ------------ ------------ ------------
# bias     CMIP5 PR     CMIP6 PR    difference
#        ------------ ------------ ------------
#        |          | |          | |          |
#        ------------ ------------ ------------
# ENSO   ------------ ------------ ------------
# change |          | |          | |          |
#        ------------ ------------ ------------
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

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import student_confidence_interval
from enso_change.wrapper.xarray_base import numpy__array, cf_dim_to_dim, mean_weighted
from enso_change.wrapper.xarray_tools import select_region
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
    # "data_diagnostics2": ["std_pr_ano_n34e", "std_ts_ano_n34e", "ave_pr_val", "ave_sst_val", "ave_ts_val_glob",
    #                       "ave_ts_val_tron"],
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val", "ave_ts_val_tron"],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_tron", "ave_sst_val"],
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
        "frac": {"x": 0.02, "y": 0.02},
        "panel_01": {"x_delt": 140, "x_size": 360, "y_delt": 30, "y_size": 160},
        "panel_11": {"x_delt": 70, "x_size": 160, "y_delt": 20, "y_size": 160},
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
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "ave_pr_val": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "VALUE (UNITS)",
                "map_c_tic": list(range(-2, 3, 1)),
                "region": None,
            },
        },
        "ave_rsst_val": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "VALUE (UNITS)",
                "map_c_tic": list(range(-2, 3, 1)),
                "region": None,
            },
        },
        "ave_rsst_val_anor": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rsst_val_asou": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rsst_val_hnor": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rsst_val_hsou": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rts_val": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "VALUE (UNITS)",
                "map_c_tic": list(range(-2, 3, 1)),
                "region": None,
            },
        },
        "ave_rts_val_anor": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rts_val_asou": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rts_val_hnor": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_rts_val_hsou": {
            "panel_11": {
                "x_lim": None,
                "x_nam": "",
                "x_tic": None,
                "y_lim": None,
                "y_nam": "VALUE (UNITS)",
                "y_tic": None,
            },
        },
        "ave_sst_val": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "VALUE (UNITS)",
                "map_c_tic": list(range(-2, 3, 1)),
                "region": None,
            },
        },
        "ave_ts_val": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "VALUE (UNITS)",
                "map_c_tic": list(range(-2, 3, 1)),
                "region": None,
            },
        },
    },
    # titles
    "fig_titles": {
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        "diffe": "CMIP6 minus CMIP5",
        "scenario": "RCP-SSP",
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
def e48_change_explained_bias_map_plot(
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
    dia, dur, pro = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat, exp = "multimodel001", "historical"
    da = values_netcdf[dia][dur][pro][exp][dat]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    # -- compute rSST
    dia2, dia3 = "ave_ts_val_tron", data_diagnostics2[-1]
    for dia1 in list(values_netcdf.keys()):
        if dia1 != dia3:
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
    for dia in [dia2, dia3]:
        if dia in list(values_netcdf.keys()):
            del values_netcdf[dia]
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
                    for dat in list(values_netcdf[dia][dur1][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur1][pro][exp][dat])
                        # select year
                        yy = int(dur1.split("_")[0]) if exp == "piControl" else 2006
                        da = da.sel(indexers={"year": yy})
                        # save
                        put_in_dict(mean_bias, da, dur1, dia, pro, exp, dat)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 == "piControl":
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
                        l1 = ["raw", "pcc"] if exp1 in ["scenario26", "scenario45", "scenario85"] else ["raw"]
                        for k1 in l1:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if k1 == "pcc":
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, k1)
    print("mean_change", list(mean_change.keys()))
    # -- compute regional mean
    default_region = {
        "hsou3050": {"latitude": [-30., -50.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou3060": {"latitude": [-30., -60.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou3070": {"latitude": [-30., -70.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou3080": {"latitude": [-30., -80.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou3090": {"latitude": [-30., -90.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou4050": {"latitude": [-40., -50.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou4060": {"latitude": [-40., -60.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou4070": {"latitude": [-40., -70.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou4080": {"latitude": [-40., -80.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hsou4090": {"latitude": [-40., -90.], "longitude": [20., 380.], "short_name": "HSOU"},
        "hnor3050": {"latitude": [30., 50.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor3060": {"latitude": [30., 60.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor3070": {"latitude": [30., 70.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor3080": {"latitude": [30., 80.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor3090": {"latitude": [30., 90.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor4050": {"latitude": [40., 50.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor4060": {"latitude": [40., 60.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor4070": {"latitude": [40., 70.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor4080": {"latitude": [40., 80.], "longitude": [20., 380.], "short_name": "HNOR"},
        "hnor4090": {"latitude": [40., 90.], "longitude": [20., 380.], "short_name": "HNOR"},
    }
    list_regions = list(default_region.keys())
    mean_bias_ave = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            for pro in list(mean_bias[dur][dia].keys()):
                for exp in list(mean_bias[dur][dia][pro].keys()):
                    for dat in list(mean_bias[dur][dia][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(mean_bias[dur][dia][pro][exp][dat])
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # loop on regions
                        for reg in list_regions:
                            da_o = copy__deepcopy(da)
                            # mask outside region
                            da_o = select_region(da_o, default_region[reg]["latitude"],
                                                 default_region[reg]["longitude"])
                            # horizontal average
                            da_o = mean_weighted(da_o, dim=[dim_y, dim_x], skipna=True)
                            # save
                            put_in_dict(mean_bias_ave, da_o, dur, str(dia) + "_" + str(reg), pro, exp, dat)
            if check_dict_levels(mean_bias_ave, dur):
                for dia_o in list(mean_bias_ave[dur].keys()):
                    if not check_dict_levels(metadata_netcdf, dia_o):
                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                        reg = dia_o.split("_")[-1][:4].replace("asou", "ATL$_S$").replace("hsou", "SH")
                        reg = reg.replace("anor", "ATL$_N$").replace("hnor", "NH") + "["
                        reg += str(int(dia_o.split("_")[-1][4:6])) + "-" + str(int(dia_o.split("_")[-1][6:8])) + "S]"
                        val = str(reg) + " " + metadata_netcdf[dia_o]["name_short"]
                        metadata_netcdf[dia_o]["name_short"] = val
    # -- compute mean maps
    remove_giss = False  # True  #
    mean_map = {}
    for dur in list(mean_bias.keys()):
        for dia in list(mean_bias[dur].keys()):
            for pro in list(mean_bias[dur][dia].keys()):
                # dictionary
                d1 = mean_bias[dur][dia][pro]
                # list datasets (i.e., minimum ensemble size)
                list_dat = sorted(list(set([dat for exp in list(d1.keys()) for dat in list(d1[exp].keys())])),
                                  key=lambda s: s.lower())
                for dat in list_dat:
                    # select models with projections
                    list_mod = []
                    for exp in list(d1.keys()):
                        if exp[:-2] == "scenario" and dat in list(d1[exp].keys()):
                            # models by number
                            list_mod += d1[exp][dat]["dataset"].values.tolist()
                    list_mod = sorted(list(set(list_mod)))
                    # remove datasets
                    if isinstance(remove_giss, bool) and remove_giss:
                        l1 = [v for k, v in dict_datasets.items() if "GISS-E2-H" in k or "GISS-E2-R" in k]
                        list_mod = [k for k in list_mod if k not in l1]
                    # separate cmip5/6
                    dt = {}
                    for pro2 in ["cmip5", "cmip6"]:
                        # model numbers corresponding to cmip5/6
                        l2 = [k for k in list_mod if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                        # select historical bias of models with projections
                        exp = "historical"
                        da = copy__deepcopy(d1[exp][dat]).sel(indexers={"dataset": l2})
                        # save
                        put_in_dict(dt, copy__deepcopy(da), pro2)
                    # compare cmip5/6
                    da1, da2 = dt["cmip5"], dt["cmip6"]
                    # average
                    da1_a, da2_a = da1.mean(dim="dataset"), da2.mean(dim="dataset")
                    # confidence interval
                    da1_c = student_confidence_interval(da1, dim="dataset", confidence_level=proc_ci)
                    da2_c = student_confidence_interval(da2, dim="dataset", confidence_level=proc_ci)
                    # overlap (0 = no overlap)
                    da_o = ((da1_c.sel(interval=0) <= da2_c.sel(interval=1)) &
                            (da2_c.sel(interval=0) <= da1_c.sel(interval=1))).astype(int)
                    # save
                    put_in_dict(mean_map, copy__deepcopy(da2_a - da1_a), dur, dat, dia, "a")
                    put_in_dict(mean_map, copy__deepcopy(da_o), dur, dat, dia, "c")
    # -- region distributions
    mean_box = {}
    for dur in list(mean_bias_ave.keys()):
        for dia in list(mean_bias_ave[dur].keys()):
            for pro in list(mean_bias_ave[dur][dia].keys()):
                # dictionary
                d1 = mean_bias_ave[dur][dia][pro]
                # list datasets (i.e., minimum ensemble size)
                list_dat = sorted(list(set([dat for exp in list(d1.keys()) for dat in list(d1[exp].keys())])),
                                  key=lambda s: s.lower())
                for dat in list_dat:
                    # select models with projections
                    list_mod = []
                    for exp in list(d1.keys()):
                        if exp[:-2] == "scenario" and dat in list(d1[exp].keys()):
                            # models by number
                            list_mod += d1[exp][dat]["dataset"].values.tolist()
                    list_mod = sorted(list(set(list_mod)))
                    # remove datasets
                    if isinstance(remove_giss, bool) and remove_giss:
                        l1 = [v for k, v in dict_datasets.items() if "GISS-E2-H" in k or "GISS-E2-R" in k]
                        list_mod = [k for k in list_mod if k not in l1]
                    # separate cmip5/6
                    for pro2 in ["cmip5", "cmip6"]:
                        # model numbers corresponding to cmip5/6
                        l2 = [k for k in list_mod if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                        # select historical bias of models with projections
                        exp = "historical"
                        da = copy__deepcopy(d1[exp][dat]).sel(indexers={"dataset": l2})
                        # save
                        put_in_dict(mean_box, copy__deepcopy(da), dur, dat, dia, pro2)
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group, panel = "group_1", "panel_01"
    for dur in list(mean_map.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for dat in list(mean_map[dur].keys()):
            ii = 0
            for jj, dia in enumerate(list(mean_map[dur][dat].keys())):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, dat, grp, pan)
                print(dur, dat, dia.ljust(12), tuple_k)
                # panel size
                tk = tuple_k[:2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, k2, *tk + (k1,))
                # panel axes definition
                if check_dict_levels(fig_ticks, dia, panel):
                    for k1, k2 in fig_ticks[dia][panel].items():
                        if check_dict_levels(figure_axes, *tuple_k + (k1,)) or \
                                (k1 == "map_c_land" and dia in ["ave_pr_val", "ave_rts_val", "ave_ts_val"]):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "map_c_nam":
                            tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                            tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                            tmp = tmp.replace(" ()", "")
                        put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                dt = {}
                if check_dict_levels(figure_axes, *tuple_k):
                    dt = read_from_dict(figure_axes, *tuple_k)
                # test possibility to plot text / legend
                do_legend = False
                if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None \
                        and dt["y_lim"] is not None:
                    do_legend = True
                # -- title
                if ii == 0:
                    val = copy__deepcopy(metadata_netcdf[dia]["name_short"]) + ": historical mean"
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, fig_titles["diffe"], *tuple_k + ("title_row",))
                # -- map
                d1 = mean_map[dur][dat][dia]
                da1, da2 = d1["a"], d1["c"] if "c" in list(d1.keys()) else None
                plot_type = "map"
                # shading
                put_in_dict(plot_data, copy__deepcopy(da1), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, "right", *tuple_k + ("legend_position",))
                # hatching
                if da2 is not None:
                    put_in_dict(plot_data, copy__deepcopy(da2), *tuple_k + (str(plot_type) + "_h",))
                    density = 3
                    kwargs_h = {"hatches": [density * "/"], "levels": [0.5, 1.5]}
                    put_in_dict(plot_data, kwargs_h, *tuple_k + ("kwargs_h",))
                # -- text
                lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
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
    # panel 11: distribution
    print("distribution")
    group, panel = "group_2", "panel_11"
    for dur in list(mean_box.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for dat in list(mean_box[dur].keys()):
            ii = 0
            list_dia = [k for k in list(mean_box[dur][dat].keys()) if "_rts_" in k or "_rsst_" in k]
            for jj, dia in enumerate(list_dia):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + (jj % 5)).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + (jj // 5)).zfill(2)
                tuple_k = (dur, dat, grp, pan)
                print(dur, dat, dia.ljust(12), tuple_k)
                # dictionary
                d1 = mean_box[dur][dat][dia]
                list_d = list(d1.keys())
                nn = len(list_d)
                minmax = [float(da.min().values) for da in d1.values()] + [float(da.max().values) for da in d1.values()]
                minmax = [math__floor(min(minmax)), math__ceil(max(minmax))]
                # panel size
                tk = tuple_k[:2]
                if not check_dict_levels(figure_psize, *tk + (pan,)):
                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_psize, *tk + (k1,)):
                            put_in_dict(figure_psize, k2, *tk + (k1,))
                # panel axes definition
                dia2 = dia[:-4]
                if check_dict_levels(fig_ticks, dia2, panel):
                    for k1, k2 in fig_ticks[dia2][panel].items():
                        if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "y_nam":
                            tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                            tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                            tmp = tmp.replace(" ()", "")
                        elif tmp is None and k1 == "x_lim":
                            tmp = [-0.5, nn - 0.5]
                        elif tmp is None and k1 == "x_tic":
                            tmp = list(range(nn))
                        elif tmp is None and k1 == "y_lim":
                            tmp = copy__deepcopy(minmax)
                        elif tmp is None and k1 == "y_tic":
                            tmp = copy__deepcopy(minmax)
                        put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                        if isinstance(tmp, list) and k1 == "x_tic":
                            put_in_dict(figure_axes, [fig_titles[k3] for k3 in list_d], *tuple_k + ("x_lab",))
                dt = {}
                if check_dict_levels(figure_axes, *tuple_k):
                    dt = read_from_dict(figure_axes, *tuple_k)
                # test possibility to plot text / legend
                do_legend = False
                if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None \
                        and dt["y_lim"] is not None:
                    do_legend = True
                # -- title
                # if ii == 0:
                #     put_in_dict(plot_data, fig_titles[pro], *tuple_k + ("title_col",))
                # if jj == 0:
                #     val = copy__deepcopy(metadata_netcdf[dia]["name_short"]) + "\nhistorical mean"
                #     put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                # -- violin
                la, lbod_fc, lc, lf, lls, llw, lmea, lmea_lw = [], [], [], [], [], [], [], []
                lmed, lmed_lw, lq, lq_ls, lq_lw, lw, lx, ly, lz = [], [], [], [], [], [], [], [], []
                for kk, pro in enumerate(list_d):
                    la += [0.4]
                    lbod_fc += ["k"]
                    lc += [fig_colors[pro] if pro in list(fig_colors.keys()) else "cyan"]
                    lf += [True]
                    lls += ["-"]
                    llw += [2]
                    lmea += [False]
                    lmea_lw += [0]
                    lmed += [True]
                    lmed_lw += [4]
                    lq += [[0.25, 0.75]]
                    lq_ls += ["--"]
                    lq_lw += [1]
                    lw += [0.9]
                    ly += [d1[pro].values]
                    lx += [kk]
                    lz += [3]
                plot_type = "vio"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
                    put_in_dict(plot_data, copy__deepcopy(lbod_fc), *tuple_k + (str(plot_type) + "_bod_fc",))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_flie",))
                    put_in_dict(plot_data, copy__deepcopy(lls), *tuple_k + (str(plot_type) + "_ls",))
                    put_in_dict(plot_data, copy__deepcopy(llw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lmea), *tuple_k + (str(plot_type) + "_mea",))
                    put_in_dict(plot_data, copy__deepcopy(lmea_lw), *tuple_k + (str(plot_type) + "_mea_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lmed), *tuple_k + (str(plot_type) + "_med",))
                    put_in_dict(plot_data, copy__deepcopy(lmed_lw), *tuple_k + (str(plot_type) + "_med_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lq), *tuple_k + (str(plot_type) + "_q",))
                    put_in_dict(plot_data, copy__deepcopy(lq_ls), *tuple_k + (str(plot_type) + "_qua_ls",))
                    put_in_dict(plot_data, copy__deepcopy(lq_lw), *tuple_k + (str(plot_type) + "_qua_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- text
                lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_psize, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
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
    a1 = int(proc_ci) if proc_ci == int(proc_ci) else copy__deepcopy(proc_ci)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_" + str(a1) + "ci"
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dat in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dat].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dat)))
            if isinstance(remove_giss, bool) and remove_giss:
                figure_name_and_path += "_noGISS-E2"
            # plot
            fig_basic(plot_data[dur][dat], list_names, 1, figure_axes[dur][dat], fig_format,
                      figure_name_and_path, figure_psize[dur][dat], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e48_change_explained_bias_map_plot()
