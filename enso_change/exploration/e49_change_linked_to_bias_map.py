# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: maps of PR, rTS, TS biases (values at the end of the 20th century) difference between top CMIP6 and bottom
# CMIP5 in terms of ENSO change
#
# Layout
#
# bias     SSP1-2.6     SSP2-4.5     SSP5-8.5
#        ------------ ------------ ------------
#        |          | |          | |          |
#        ------------ ------------ ------------
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
from numpy import nan as numpy__nan
from numpy import sign as numpy__sign

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import student_confidence_interval
from enso_change.wrapper.xarray_base import concatenate, create_array, xarray
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
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_ts_val_glob", "ave_ts_val_tron", "ave_ts_val"],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_ts_val_glob", "ave_ts_val_tron", "ave_ts_val"],
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
        "frac": {"x": 0.02, "y": 0.02},
        "panel_01": {"x_delt": 70, "x_size": 360, "y_delt": 20, "y_size": 160},
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
        "ave_pr_val--bias": {
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
        "ave_pr_val--change": {
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
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "region": None,
            },
        },
        "ave_rsst_val--bias": {
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
        "ave_rsst_val--change": {
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
        "ave_rts_val--bias": {
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
        "ave_rts_val--change": {
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
        "ave_sst_val--bias": {
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
        "ave_sst_val--change": {
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
        "ave_ts_val--bias": {
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
        "ave_ts_val--change": {
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
        "diffe": "CMIP6 bottom THR% minus CMIP5 bottom THR%",
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
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        if dat[-3:] in ["004", "005"]:
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
    for dia in [dia2]:  # [dia2, dia3]:
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
    # -- select top / bottom,  CMIP6 / CMIP5
    remove_giss = False  # True  #
    mean_map = {}
    for typ, d1 in zip(["bias", "change"], [mean_bias, mean_change]):
        for dur in list(d1.keys()):
            if not check_dict_levels(mean_change, dur):
                continue
            for dia_x in list(d1[dur].keys()):
                if len(dia_x.split("_")) != 3:
                    continue
                for pro in list(d1[dur][dia].keys()):
                    for exp in list(d1[dur][dia][pro].keys()):
                        if exp == "historical":
                            continue
                        for dat in list(d1[dur][dia][pro][exp].keys()):
                            for dia_y in list(mean_change[dur].keys()):
                                if not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat) or \
                                        not (len(dia_y.split("_")) == 4 and dia_y.split("_")[0] in ["std", "var"]):
                                    continue
                                for war in list(mean_change[dur][dia_y][pro][exp][dat].keys()):
                                    # array
                                    da_x = copy__deepcopy(d1[dur][dia_x][pro][exp][dat])
                                    if typ == "change":
                                        da_x = da_x[war]
                                    da_y = copy__deepcopy(mean_change[dur][dia_y][pro][exp][dat][war])
                                    # list models
                                    list_mod = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                           set(da_y["dataset"].values.tolist())))
                                    # remove datasets
                                    if isinstance(remove_giss, bool) and remove_giss:
                                        l1 = [v for k, v in dict_datasets.items()
                                              if "GISS-E2-H" in k or "GISS-E2-R" in k]
                                        list_mod = [k for k in list_mod if k not in l1]
                                    # separate cmip5/6
                                    dt = {}
                                    for pro2 in ["cmip5", "cmip6"]:
                                        # model numbers corresponding to cmip5/6
                                        l2 = [k for k in list_mod
                                              if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                                        # select models
                                        da_xt = da_x.sel(indexers={"dataset": l2})
                                        da_yt = da_y.sel(indexers={"dataset": l2}).compute()
                                        # for thr in ["bottom", "top"]:
                                        thr = "bottom"
                                        # find top / bottom of ENSO SSTA change
                                        t1 = (copy__deepcopy(proc_th) if thr == "bottom" else 100 - proc_th) / 100
                                        t2 = da_yt.quantile(t1)
                                        cond = da_yt <= t2 if thr == "bottom" else da_yt >= t2
                                        l2 = da_yt.where(cond, drop=True)["dataset"].values.tolist()
                                        # select top / bottom of ENSO SSTA change
                                        da_o = copy__deepcopy(da_xt).sel(indexers={"dataset": l2})
                                        put_in_dict(dt, da_o, thr, pro2)
                                    # compare cmip5/6
                                    for thr in list(dt.keys()):
                                        da1, da2 = dt[thr]["cmip5"], dt[thr]["cmip6"]
                                        # difference of average
                                        da_a = da2.mean(dim="dataset") - da1.mean(dim="dataset")
                                        if da1["dataset"].size < 4 or da2["dataset"].size < 4:
                                            da_a, da_c = create_array(da_a, value=None), create_array(da_a, value=0)
                                        else:
                                            # confidence interval
                                            dk = {"dim": "dataset", "confidence_level": proc_ci}
                                            da1_c = student_confidence_interval(da1, **dk)
                                            da1_mi, da1_ma = da1.min(dim="dataset"), da1.max(dim="dataset")
                                            da1_c = da1_c.where(da1_c >= da1_mi, other=da1_mi)
                                            da1_c = da1_c.where(da1_c <= da1_ma, other=da1_ma)
                                            da2_c = student_confidence_interval(da2, **dk)
                                            da2_mi, da2_ma = da2.min(dim="dataset"), da2.max(dim="dataset")
                                            da2_c = da2_c.where(da2_c >= da2_mi, other=da2_mi)
                                            da2_c = da2_c.where(da2_c <= da2_ma, other=da2_ma)
                                            # overlap (0 = no overlap)
                                            da_c = ((da1_c.sel(interval=0) <= da2_c.sel(interval=1)) &
                                                    (da2_c.sel(interval=0) <= da1_c.sel(interval=1))).astype(int)
                                        # save
                                        tk = (dur, typ, dia_y, dia_x, war, dat, exp)
                                        put_in_dict(mean_map, copy__deepcopy(da_a), *tk + ("a",))
                                        put_in_dict(mean_map, copy__deepcopy(da_c), *tk + ("c",))
    # -- agreement
    for dur in list(mean_map.keys()):
        for typ in list(mean_map[dur].keys()):
            for dia_y in list(mean_map[dur][typ].keys()):
                for dia_x in list(mean_map[dur][typ][dia_y].keys()):
                    for war in list(mean_map[dur][typ][dia_y][dia_x].keys()):
                        for ex1 in ["scenario"]:  # ["historical", "scenario"]:
                            # regroup correlation values
                            da = ()
                            for dat in list(mean_map[dur][typ][dia_y][dia_x][war].keys()):
                                for ex2, d2 in mean_map[dur][typ][dia_y][dia_x][war][dat].items():
                                    if ex1 in ex2 and not d2["a"].isnull().all().values:
                                        da += (d2["a"],)
                            if len(da) == 0:
                                continue
                            da = concatenate(da, dim="map")
                            # median correlation across maps
                            median_va = da.median(dim="map", skipna=True)
                            # sign of the median
                            median_sign = xarray.where(median_va > 0, 1, xarray.where(median_va < 0, -1, numpy__nan))
                            # does each map have the same sign as the median?
                            same_sign = xarray.where(numpy__sign(da) == median_sign, 1.0, 0.0)
                            # fraction of maps agreeing with the median sign
                            agreement = same_sign.where(da.notnull()).mean(dim="map")
                            # save
                            tk = (dur, typ, dia_y, dia_x, war, "median & agreement", ex1)
                            put_in_dict(mean_map, median_va, *tk + ("a",))
                            put_in_dict(mean_map, agreement, *tk + ("c",))
    print("agreement", list(mean_map.keys()))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group, panel = "group_1", "panel_01"
    for dur in list(mean_map.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for typ in list(mean_map[dur].keys()):
            for dia_y in list(mean_map[dur][typ].keys()):
                for dia_x in list(mean_map[dur][typ][dia_y].keys()):
                    for war in list(mean_map[dur][typ][dia_y][dia_x].keys()):
                        list_dat = list(mean_map[dur][typ][dia_y][dia_x][war].keys())
                        for ii, dat in enumerate(list_dat):
                            list_exp = list(mean_map[dur][typ][dia_y][dia_x][war][dat].keys())
                            for jj, exp in enumerate(list_exp):
                                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                                tuple_k = (dur, typ, dia_y, dia_x, war, grp, pan)
                                print(dur, typ, dia_y, dia_x.ljust(12), war, tuple_k)
                                # panel size
                                tk = tuple_k[:-2]
                                if not check_dict_levels(figure_psize, *tk + (pan,)):
                                    put_in_dict(figure_psize, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                                    for k1, k2 in fig_panel_size.items():
                                        if k1[:6] != "panel_" and \
                                                not check_dict_levels(figure_psize, *tk + (k1,)):
                                            put_in_dict(figure_psize, k2, *tk + (k1,))
                                # panel axes definition
                                dia_x2 = copy__deepcopy(dia_x) + "--" + str(typ)
                                if check_dict_levels(fig_ticks, dia_x2, panel):
                                    for k1, k2 in fig_ticks[dia_x2][panel].items():
                                        if check_dict_levels(figure_axes, *tuple_k + (k1,)) or \
                                                (k1 == "map_c_land" and
                                                 dia_x in ["ave_pr_val", "ave_rts_val", "ave_ts_val"]):
                                            continue
                                        tmp = copy__deepcopy(k2)
                                        k3 = str(typ) + "--" + str(dia_x) + "--" + str(dia_y)
                                        if isinstance(tmp, str) and k1 == "map_c_nam":
                                            tmp = tmp.replace("VALUE", metadata_netcdf[dia_x]["name_short"])
                                            tmp = tmp.replace("UNITS", metadata_netcdf[dia_x]["units"])
                                            tmp = tmp.replace(" ()", "")
                                        elif k1 == "region" and k3 == "bias--ave_pr_val--std_ts_ano_n34e":
                                            tmp = ["prae", "prie", "prpe", "prps"]
                                        elif k1 == "region" and k3 == "bias--ave_rts_val--std_ts_ano_n34e":
                                            tmp = ["pbhw", "prps", "n50e", "n30e"]
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
                                    val = copy__deepcopy(fig_titles[exp])
                                    val += " " + copy__deepcopy(fig_titles["diffe"]).replace("THR", str(int(proc_th)))
                                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                                if jj == 0:
                                    val = str(int(dat[-3:])) + " members" if "multimodel" in dat else str(dat)
                                    put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                                # -- map
                                d1 = mean_map[dur][typ][dia_y][dia_x][war][dat][exp]
                                da1, da2 = d1["a"], d1["c"] if "c" in list(d1.keys()) else None
                                ptype = "map"
                                # shading
                                put_in_dict(plot_data, copy__deepcopy(da1), *tuple_k + (str(ptype) + "_c",))
                                val = "right" if exp == list_exp[-1] else None
                                put_in_dict(plot_data, val, *tuple_k + ("legend_position",))
                                # hatching
                                if da2 is not None:
                                    put_in_dict(plot_data, copy__deepcopy(da2), *tuple_k + (str(ptype) + "_h",))
                                    kwargs_h = {"hatches": [3 * "/"],
                                                "levels": [0.5, 1.5] if "multimodel" in dat else [0.0, 0.8]}
                                    put_in_dict(plot_data, kwargs_h, *tuple_k + ("kwargs_h",))
                                # -- text
                                lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                                if isinstance(do_legend, bool) and do_legend and ii == 0 and jj == 0:
                                    df = read_from_dict(figure_psize, *tk)
                                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                    tx = df[pan]["x_size"] * len(list_exp) + df[pan]["x_delt"] * (len(list_exp) - 1)
                                    dy = (y2 - y1) * default_plot["size_y"] / (
                                            100 * df[pan]["y_size"] * df["frac"]["y"])
                                    # super title
                                    lc += ["k"]
                                    lf += [30]
                                    lh += ["center"]
                                    lr += [0]
                                    lt += [typ]
                                    lv += ["bottom"]
                                    lw += ["bold"]
                                    lx += [x1 - (x2 - x1) / 2 + tx / 2]
                                    ly += [y2 + 30 * dy]
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
            for dia_y in sorted(list(plot_data[dur][typ].keys()), key=lambda s: s.lower()):
                for dia_x in sorted(list(plot_data[dur][typ][dia_y].keys()), key=lambda s: s.lower()):
                    for war in sorted(list(plot_data[dur][typ][dia_y][dia_x].keys()), key=lambda s: s.lower()):
                        # list plot names
                        list_names = list(plot_data[dur][typ][dia_y][dia_x][war].keys())
                        # path to output figure
                        figure_name_and_path = os__path__join(
                            figure_directory, figure_name + "_".join(("", dur, typ, dia_y, dia_x, war)))
                        if isinstance(remove_giss, bool) and remove_giss:
                            figure_name_and_path += "_noGISS-E2"
                        # plot
                        fig_basic(plot_data[dur][typ][dia_y][dia_x][war], list_names, 1,
                                  figure_axes[dur][typ][dia_y][dia_x][war], fig_format, figure_name_and_path,
                                  figure_psize[dur][typ][dia_y][dia_x][war], fig_title_bool=True,
                                  panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e48_change_explained_bias_map_plot()
