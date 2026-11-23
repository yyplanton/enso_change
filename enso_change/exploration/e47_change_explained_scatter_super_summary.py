# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: Scatterplots of biases (values at the end of the 20th century) or mean climate changes (departure from the end
# of the 20th century) and ENSO changes (departure from the end of the 20th century) in selected regions
#
# Layout
#
# bias        region1      region2      region N
#           ------------ ------------ ------------
#       his |          | |          | |          |
#           ------------ ------------ ------------
#           ------------ ------------ ------------
#       r/s |          | |          | |          |
#           ------------ ------------ ------------
# change      region1      region2      region N
#           ------------ ------------ ------------
#       his |          | |          | |          |
#           ------------ ------------ ------------
#           ------------ ------------ ------------
#       r/s |          | |          | |          |
#           ------------ ------------ ------------
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
from numpy.random import randint as numpy__random__randint

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.plots.panels import default_plot, default_region
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import correlation_pearson
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
    "data_diagnostics2": ["std_pr_ano_n34e", "std_ts_ano_n34e", "ave_pr_val", "ave_sst_val", "ave_ts_val_glob",
                          "ave_ts_val_tron"],
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
    "proc_ci": 95.,
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
    "fig_panel_size": {
        "frac": {"x": 0.5, "y": 0.5},
        "panel_01": {"x_delt": 3, "x_size": 8, "y_delt": 3, "y_size": 8},
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
    "fig_markers": {
        "multimodel001": "o",
        "multimodel002": "v",
        "multimodel003": "X",
        "multimodel004": "P",
        "multimodel005": ">",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "correlation": {
            "panel_01": {
                "x_lim": [-1, 1],
                "x_nam": "corr(VALUE1, VALUE2)",
                "x_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "y_lim": [-1, 1],
                "y_nam": "per GMST",
                "y_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
            },
        },
    },
    # titles
    "fig_titles": {
        "scenario": "RCP-SSP",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x":40,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e47_change_explained_map_super_summary_plot(
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
                        n = values_netcdf[dia][dur][pro][exp][dat]["dataset"].size
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
    dia2 = "ave_ts_val_tron"
    for dia1 in list(values_netcdf.keys()):
        if dia1 != "ave_sst_val":
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
        metadata_netcdf[dia_o]["name_short"] = metadata_netcdf[dia1]["name_short"].replace("SST", "rSST")
    if dia2 in list(values_netcdf.keys()):
        del values_netcdf[dia2]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- select bias year
    dur2, exp2 = data_epoch_lengths[-1], "historical"
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
    dia2 = "ave_ts_val_glob"
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
                        for war in l1:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if war == "pcc":
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- compute region averages
    list_regions = ["anor", "apme", "asou", "hnor", "hsou", "itse", "itsw", "lono", "mari", "n20b", "n20e", "n20n",
                    "n20s", "n30b", "n30e", "n30n", "n30s", "n34b", "n34e", "n34n", "n34s", "n40b", "n40e", "n40n",
                    "n40s", "n43b", "n43e", "n43n", "n43s", "n50b", "n50e", "n50n", "n50s", "pbes", "pbhw", "pbnw",
                    "pbse", "pbsw", "pece", "pecw", "pnee", "pnww", "pone", "ptsc", "ptse", "soip", "sooc"]
    l1, l2 = ["n43e", "n20e", "asou", "hsou", "n34s", "n43s"], ["n50e", "pecw", "anor", "hnor", "n34n", "n43n"]
    l3 = ["pgx1", "pgx2", "agy1", "hgy1", "pgy1", "pgy2"]
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
                        # regional differences
                        for r1, r2, nn in zip(l1, l2, l3):
                            if r1 not in list_regions or r2 not in list_regions or dia == "ave_rsst_val":
                                continue
                            da1 = read_from_dict(mean_bias_ave, dur, str(dia) + "_" + str(r1), pro, exp, dat)
                            da2 = read_from_dict(mean_bias_ave, dur, str(dia) + "_" + str(r2), pro, exp, dat)
                            put_in_dict(mean_bias_ave, da2 - da1, dur, str(dia) + "_" + str(nn), pro, exp, dat)
                        if check_dict_levels(mean_bias_ave, dur):
                            for dia_o in list(mean_bias_ave[dur].keys()):
                                if not check_dict_levels(metadata_netcdf, dia_o):
                                    metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                                    reg = dia_o.split("_")[-1]
                                    val = str(reg).upper() + " " + metadata_netcdf[dia_o]["name_short"]
                                    metadata_netcdf[dia_o]["name_short"] = val
    mean_change_ave = {}
    for dur in list(mean_change.keys()):
        for dia in list(mean_change[dur].keys()):
            if len(dia.split("_")) != 3:
                continue
            for pro in list(mean_change[dur][dia].keys()):
                for exp in list(mean_change[dur][dia][pro].keys()):
                    for dat in list(mean_change[dur][dia][pro][exp].keys()):
                        for war in list(mean_change[dur][dia][pro][exp][dat].keys()):
                            # array
                            da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat][war])
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
                                lk = (dur, str(dia) + "_" + str(reg), pro, exp, dat, war)
                                put_in_dict(mean_change_ave, da_o,*lk)
                            # regional differences
                            for r1, r2, nn in zip(l1, l2, l3):
                                if r1 not in list_regions or r2 not in list_regions or dia == "ave_rsst_val":
                                    continue
                                lk = (pro, exp, dat, war)
                                da1 = read_from_dict(mean_change_ave, dur, str(dia) + "_" + str(r1), *lk)
                                da2 = read_from_dict(mean_change_ave, dur, str(dia) + "_" + str(r2), *lk)
                                lk = (dur, str(dia) + "_" + str(nn), pro, exp, dat, war)
                                put_in_dict(mean_change_ave, da2 - da1, *lk)
                            if check_dict_levels(mean_change_ave, dur):
                                for dia_o in list(mean_change_ave[dur].keys()):
                                    if not check_dict_levels(metadata_netcdf, dia_o):
                                        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
                                        reg = dia_o.split("_")[-1]
                                        val = str(reg).upper() + " " + metadata_netcdf[dia_o]["name_short"]
                                        metadata_netcdf[dia_o]["name_short"] = val
    # -- compute correlations
    dict_reg = {
        "std_pr_ano_n34e": {
            "bias": {
                "historical": {
                    "ave_pr_val": ["asou", "n34n", "n34s", "pgy1"],
                    "ave_rsst_val": ["lono"],
                    "ave_sst_val": ["lono"],
                },
                "scenario": {
                    "ave_pr_val": ["anor", "n34e", "pbsw", "pece"],
                    "ave_rsst_val": ["asou", "n50b", "pbhw", "pece", "pgx1"],
                    "ave_sst_val": ["asou", "n50b", "pbhw", "pgx1"],
                },
            },
            "change": {
                "historical": {
                    "ave_pr_val": ["pbes", "pbsw", "pone"],
                    "ave_rsst_val": ["asou", "hgy1", "n50b", "pbes", "pbsw", "pgx1", "sooc"],
                    "ave_sst_val": ["pbes"],
                },
                "scenario": {
                    "ave_pr_val": ["lono", "pnee", "pnww", "sooc"],
                    "ave_rsst_val": ["agy1", "anor", "asou", "hgy1", "n20b", "ptse", "ptsc"],
                    "ave_sst_val": ["anor"],
                },
            },
        },
        "std_ts_ano_n34e": {
            "bias": {
                "historical": {
                    "ave_pr_val": ["pone", "pbes", "pbsw"],
                    "ave_rsst_val": ["agy1", "anor", "pbnw", "pgx1", "sooc"],
                    "ave_sst_val": ["agy1", "anor", "pbnw", "sooc"],
                },
                "scenario": {
                    "ave_pr_val": ["asou", "n30e", "pbes", "pbsw"],
                    "ave_rsst_val": ["agy1", "asou", "mari", "pbnw", "sooc"],
                    "ave_sst_val": ["agy1", "asou", "mari", "pbnw", "sooc"],
                },
            },
            "change": {
                "historical": {
                    "ave_pr_val": ["n40e", "pbes", "pbsw"],
                    "ave_rsst_val": ["n50b", "pbse", "pgx1"],
                    "ave_sst_val": ["pbse"],
                },
                "scenario": {
                    "ave_pr_val": ["apme", "itsw", "pgx1", "pone"],
                    "ave_rsst_val": ["agy1", "anor", "asou", "hgy1", "itsw", "mari", "n20b", "pbnw", "pecw", "pgx2",
                                     "ptsc", "soip"],
                    "ave_sst_val": ["anor", "", "soip"],
                },
            },
        },
    }
    remove_giss = False  # True  #
    corr_ave = {}
    for typ, d1 in zip(["bias", "change"], [mean_bias_ave, mean_change_ave]):
        for dur in list(d1.keys()):
            if not check_dict_levels(mean_change, dur):
                continue
            for dia_x in list(d1[dur].keys()):
                if len(dia_x.split("_")) != 4:
                    continue
                reg, var = dia_x.split("_")[-1], "_".join(dia_x.split("_")[:3])
                for dia_y in list(mean_change[dur].keys()):
                    if not (len(dia_y.split("_")) == 4 and dia_y.split("_")[0] in ["std", "var"]):
                        continue
                    for pro in list(d1[dur][dia_x].keys()):
                        for exp in list(d1[dur][dia_x][pro].keys()):
                            ex2 = copy__deepcopy(exp) if exp == "historical" else exp[:-2]
                            # if not check_dict_levels(dict_reg, dia_y, typ, ex2, var) or \
                            #         reg not in dict_reg[dia_y][typ][ex2][var]:
                            #     continue
                            for dat in list(d1[dur][dia_x][pro][exp].keys()):
                                if not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat):
                                    continue
                                for war in list(mean_change[dur][dia_y][pro][exp][dat].keys()):
                                    if not check_dict_levels(mean_change_ave, dur, dia_x, pro, exp, dat, war) \
                                            and typ == "change":
                                        continue
                                    # arrays
                                    da_x = d1[dur][dia_x][pro][exp][dat]
                                    if typ == "change":
                                        da_x = da_x[war]
                                    da_y = mean_change[dur][dia_y][pro][exp][dat][war]
                                    # remove datasets
                                    if isinstance(remove_giss, bool) and remove_giss:
                                        l1 = [v for k, v in dict_datasets.items()
                                              if "GISS-E2-H" in k or "GISS-E2-R" in k]
                                        list_d = [k for k in da_x["dataset"].values if k not in l1]
                                        da_x = da_x.sel(indexers={"dataset": list_d})
                                        da_y = da_y.sel(indexers={"dataset": list_d})
                                    # correlation and p-value
                                    da_r = correlation_pearson(da_x, da_y, dim="dataset", skipna=True)
                                    da_r = float(da_r.values)
                                    # save
                                    put_in_dict(corr_ave, da_r, dur, pro, dia_y, typ, dia_x, dat, exp, war)
    print("corr_ave", list(corr_ave.keys()))
    # -- correlation agreement
    agreement_ave = {}
    for dur in list(corr_ave.keys()):
        for pro in list(corr_ave[dur].keys()):
            for dia_y in list(corr_ave[dur][pro].keys()):
                for typ in list(corr_ave[dur][pro][dia_y].keys()):
                    for ex1 in ["historical", "scenario"]:
                        dt = {}
                        for dia_x in list(corr_ave[dur][pro][dia_y][typ].keys()):
                            # dictionary
                            d1 = corr_ave[dur][pro][dia_y][typ][dia_x]
                            # regroup correlation values
                            da = numpy__array([v for dat in list(d1.keys()) for ex2 in list(d1[dat].keys())
                                               for v in d1[dat][ex2].values() if ex1 in ex2])
                            if da.size == 0:
                                continue
                            # save if good correlation
                            a1, a2 = max(abs(da)), float(da.mean())
                            if a1 > 0.3 or a2 > 0.3:
                                dt[dia_x] = copy__deepcopy(a2)
                                # if a1 > 0.5 and a2 > 0.5:
                                #     dt[dia_x] = max([a1, a2])
                                # elif a1 > 0.5 >= a2:
                                #     dt[dia_x] = copy__deepcopy(a1)
                                # else:
                                #     dt[dia_x] = copy__deepcopy(a2)
                        print(dur, pro, dia_y, typ, ex1)
                        if len(list(dt.keys())) > 0:
                            dt = dict(sorted(dt.items(), key=lambda item: abs(item[1]), reverse=True))
                            for thr in [0.1, 0.2, 0.3, 0.4]:
                                if len(dt.keys()) > 10:
                                    dt2 = {k: v for k, v in dt.items() if abs(v) >= thr}
                                    if len(list(dt2.keys())) > 4:
                                        dt = copy__deepcopy(dt2)
                                    else:
                                        break
                            aa = ", ".join([str(k) + " (" + "{0:.2f}".format(round(v, 2)).rjust(5) + ")"
                                            for k, v in dt.items()])
                            print(dur, pro, dia_y, typ, ex1, aa)
                            put_in_dict(agreement_ave, list(dt.keys()), dur, pro, dia_y, typ, ex1)
    print("agreement_ave", list(agreement_ave.keys()))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_psize, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group, panel = "group_1", "panel_01"
    for dur in list(agreement_ave.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for pro in list(agreement_ave[dur].keys()):
            for dia_y in list(agreement_ave[dur][pro].keys()):
                ii = 0
                for typ in list(agreement_ave[dur][pro][dia_y].keys()):
                    list_exp = list(agreement_ave[dur][pro][dia_y][typ].keys())
                    for exp in list_exp:
                        list_dia_x = agreement_ave[dur][pro][dia_y][typ][exp]
                        for jj, dia_x in enumerate(list_dia_x):
                            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                            tuple_k = (dur, pro, dia_y, grp, pan)
                            print(dur, pro, dia_y, typ.ljust(6), exp.ljust(10), dia_x.ljust(17), tuple_k)
                            # dictionary
                            dict_t = corr_ave[dur][pro][dia_y][typ][dia_x]
                            # panel size
                            tk = tuple_k[:3]
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
                                    if isinstance(tmp, str) and k1 in ["x_nam", "y_nam"]:
                                        tmp = tmp.replace("VALUE1", metadata_netcdf[dia_x]["name_short"])
                                        tmp = tmp.replace("VALUE2", metadata_netcdf[dia_y]["name_short"])
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
                            if jj == 0:
                                val = str(fig_titles[exp]) + ": "
                                val += "1850-" + str(1850 + ny2 - 1) if exp == "historical" else \
                                    str(2099 - ny2 + 1) + "-2099"
                                val += "\nvs. " + str(2005 - ny2 + 1) + "-2005"
                                put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                            # -- marker
                            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                            for dat in list(dict_t.keys()):
                                for ex2 in list(dict_t[dat].keys()):
                                    if exp not in ex2 or "raw" not in list(dict_t[dat][ex2].keys()):
                                        continue
                                    lc += [fig_colors[ex2] if ex2 in list(fig_colors.keys()) else "cyan"]
                                    lm += [fig_markers[dat] if dat in list(fig_markers.keys()) else "o"]
                                    ls += [50]
                                    lw += [0.5]
                                    lx += [dict_t[dat][ex2]["raw"]]
                                    ly += [dict_t[dat][ex2]["pcc"] if "pcc" in list(dict_t[dat][ex2].keys()) else 0]
                                    lz += [numpy__random__randint(3, 9)]
                            # mean
                            mean = 0
                            if len(lx) > 0 and len(lx) == len(ly):
                                lc += ["royalblue"]
                                lm += ["D"]
                                ls += [80]
                                lw += [0.5]
                                lx += [float(numpy__array(lx).mean())]
                                ly += [float(numpy__array(ly).mean())]
                                lz += [2]
                                if len(list(set(ly))) == 1 and ly[0] == 0:
                                    mean = float(numpy__array(lx).mean())
                                else:
                                    mean = float(numpy__array(lx + ly).mean())
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
                            lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_psize, *tk)
                                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                                # variable
                                lc += ["k"]
                                lf += [15]
                                lh += ["right"]
                                lt += [dia_x.split("_")[1].upper().replace("RSST", "rSST")]
                                lv += ["top"]
                                lw += ["regular"]
                                lx += [x2 - 2 * dx]
                                ly += [y2 - 2 * dy]
                                # average correlation
                                lc += ["k"]
                                lf += [15]
                                lh += ["right"]
                                lt += ["r=" + "{0:.2f}".format(round(mean, 2))]
                                lv += ["top"]
                                lw += ["regular"]
                                lx += [x2 - 2 * dx]
                                ly += [y2 - 7 * dy]
                            plot_type = "text"
                            if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                                put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                                put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                                put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_w",))
                                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                        ii += 1
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
        for pro in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for dia in sorted(list(plot_data[dur][pro].keys()), key=lambda s: s.lower()):
                # list plot names
                list_names = list(plot_data[dur][pro][dia].keys())
                # path to output figure
                figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia, pro)))
                if isinstance(remove_giss, bool) and remove_giss:
                    figure_name_and_path += "_noGISS-E2"
                # plot
                fig_basic(plot_data[dur][pro][dia], list_names, 1, figure_axes[dur][pro][dia], fig_format,
                          figure_name_and_path, figure_psize[dur][pro][dia], fig_title_bool=True,
                          panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e47_change_explained_map_super_summary_plot()
