# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: Maps of correlations between biases (values at the end of the 20th century) and changes (beginning or end of
# experiments minus end of 20th century)
# Here biases (historical) are selected
# Here changes (departure from the end of the 20th century) are computed per degree of GMST
# Maps of correlations are plotted
#
#     ------------ ------------ ------------ ------------
# 001 |historical| |scenario24| |scenario45| |scenario85| correlation maps using members with at least 1 member
#     ------------ ------------ ------------ ------------
#     ------------ ------------ ------------ ------------
# 002 |historical| |scenario24| |scenario45| |scenario85| correlation maps using members with at least 2 members
#     ------------ ------------ ------------ ------------
#     ------------ ------------ ------------ ------------
# n   |historical| |scenario24| |scenario45| |scenario85| correlation maps using members with at least n members
#     ------------ ------------ ------------ ------------
#
# Notes:
# - southern ocean correlations could be caused by GISS-E2-H / GISS-E2-R
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
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import correlation_pearson, correlation_pearson_pvalue
from enso_change.wrapper.xarray_base import cf_dim_to_dim, mean_weighted
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics1": [],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_sst_val", "ave_ts_val_tron"],
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_sst_val"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    "data_experiments": ["historical", "scenario26", "scenario45", "scenario85"],
    # "data_experiments": ["historical", "scenario26"],
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
        "frac": {"x": 0.02, "y": 0.02},
        "panel_01": {"x_delt": 50, "x_size": 360, "y_delt": 20, "y_size": 160},
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
        "correlation": {
            "panel_01": {
                "x_lim": [20, 380],
                "x_nam": "",
                "x_tic": [90, 180, 270, 360],  # list(range(60, 381, 60)),
                "y_lim": [-80, 80],
                "y_nam": "",
                "y_tic": list(range(-60, 61, 30)),
                "map_c_cbar": "cmo.balance",
                "map_c_land": "dimgrey",
                "map_c_nam": "correlation",
                "map_c_tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                "region": None,
            },
        },
    },
    # titles
    "fig_titles": {
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
def e41_change_explained_map_bias_cc_plot(
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
    # for dd in [values_json, values_netcdf]:
    #     for dia in list(dd.keys()):
    #         for dur in list(dd[dia].keys()):
    #             for pro in list(dd[dia][dur].keys()):
    #                 for exp in list(dd[dia][dur][pro].keys()):
    #                     for dat in list(dd[dia][dur][pro][exp].keys()):
    #                         if dat != "multimodel" + str(smile_minimum_size).zfill(3):
    #                             del dd[dia][dur][pro][exp][dat]
    #                     if len(list(dd[dia][dur][pro][exp].keys())) == 0:
    #                         del dd[dia][dur][pro][exp]
    #                 if len(list(dd[dia][dur][pro].keys())) == 0:
    #                     del dd[dia][dur][pro]
    #             if len(list(dd[dia][dur].keys())) == 0:
    #                 del dd[dia][dur]
    #         if len(list(dd[dia].keys())) == 0:
    #             del dd[dia]
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
    # dat = "multimodel005"
    # dd = {}
    # for exp in list(values_netcdf[dia][dur][pro].keys()):
    #     if exp not in ["scenario26", "scenario45"]:
    #         continue
    #     da = values_netcdf[dia][dur][pro][exp][dat]
    #     ld = [kk.split(" ")[1].split("=")[0] for kk in da.attrs["dataset"].split("; ")]
    #     dd[exp] = ld
    # for exp1 in list(dd.keys()):
    #     for exp2 in list(dd.keys()):
    #         if exp1 != exp2:
    #             print("in " + str(exp1) + " but not in " + str(exp2),
    #                   ", ".join(sorted(list(set(dd[exp1]) - set(dd[exp2])), key=lambda s: s.lower())))
    # stop
    # -- compute rSST
    dia2 = "ave_ts_val_tron"
    # for dia1 in list(values_netcdf.keys()):
    #     if dia1 != "ave_sst_val":
    #         continue
    #     for dur in list(values_netcdf[dia1].keys()):
    #         for pro in list(values_netcdf[dia1][dur].keys()):
    #             for exp in list(values_netcdf[dia1][dur][pro].keys()):
    #                 for dat in list(values_netcdf[dia1][dur][pro][exp].keys()):
    #                     if not check_dict_levels(values_netcdf, dia2, dur, pro, exp, dat):
    #                         continue
    #                     print(dia1, dur, pro, exp, dat, dia2)
    #                     # arrays
    #                     da1, da2 = values_netcdf[dia1][dur][pro][exp][dat], values_netcdf[dia2][dur][pro][exp][dat]
    #                     print("da1", da1.shape, da1.dims)
    #                     print("da2", da2.shape, da2.dims)
    #                     # remove tropical mean SST
    #                     da1 -= da2
    #                     # save
    #                     dia_o = dia1.replace("_sst_", "_rsst_").replace("_ts_", "_rts_")
    #                     put_in_dict(values_netcdf, da1, dia_o, dur, pro, exp, dat)
    for dia in list(values_netcdf.keys()):
        if dia != "ave_sst_val":
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # arrays
                        da = values_netcdf[dia][dur][pro][exp][dat]
                        # remove horizontal mean
                        dim_lon, dim_lat = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        da_r = da.sel(indexers={dim_lat: slice(*(-20, 20))})
                        da = da - mean_weighted(da_r, dim=[dim_lat, dim_lon], skipna=True)
                        # save
                        dia_o = dia.replace("_sst_", "_rsst_").replace("_ts_", "_rts_")
                        put_in_dict(values_netcdf, da, dia_o, dur, pro, exp, dat)
    if dia2 in list(values_netcdf.keys()):
        del values_netcdf[dia2]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # -- select bias year
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
                        # -- check outliers
                        # if ' ':
                        # #     da_t = da.sel(indexers={"latitude": [k for k in da["latitude"].values if abs(k) < 60]})
                        #     da_t = da.sel(indexers={"latitude": [k for k in da["latitude"].values if -60 <= k <= -40]})
                        #     mask = da_t.isnull().any(dim="dataset")
                        #     da_t = da_t.where(~mask)
                        #     da_t = (da_t - da_t.mean(dim="dataset", skipna=True)) / da_t.std(dim="dataset", skipna=True)
                        #     if abs(da_t).max() > 5:
                        #         print(dia, dur1, pro, exp, dat, "outliers")
                        #         bad_datasets = da_t.where(abs(da_t) > 5).notnull().any(dim=("latitude", "longitude")).compute()
                        #         bad_datasets_coords = da_t["dataset"].where(bad_datasets, drop=True)
                        #         dd = {}
                        #         for kk in da.attrs["dataset"].split("; "):
                        #             dd[int(kk.split(" ")[0])] = kk.split(" ")[1]
                        #         aa = ", ".join([dd[int(kk)] for kk in bad_datasets_coords.values])
                        #         print("outliers: " + aa)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 4:
            continue
        for dur1 in list(values_netcdf[dia].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_netcdf[dia][dur1].keys()):
                for exp1 in list(values_netcdf[dia][dur1][pro].keys()):
                    if exp1 == "piControl":
                        continue
                    for dat in list(values_netcdf[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_netcdf, dia, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_netcdf[dia][dur1][pro][exp1][dat])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da2 = copy__deepcopy(values_netcdf[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_netcdf[dia][dur2][pro][exp2][dat])
                        # select experiment epoch
                        yy = da1["year"].values
                        yy = yy[0] if exp1 == "historical" else 2100
                        da1 = da1.sel(indexers={"year": yy})
                        # select reference epoch
                        yy = 2006
                        da2 = da2.sel(indexers={"year": yy})
                        # change
                        da1 = da1 - da2 if exp1 in ["scenario26", "scenario45", "scenario85"] else da2 - da1
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(mean_change, da1, dur1, dia, pro, exp1, dat)
    print("mean_change", list(mean_change.keys()))
    # -- compute correlations
    corr_map = {}
    for dur in list(mean_bias.keys()):
        if not check_dict_levels(mean_change, dur):
            continue
        for dia_x in list(mean_bias[dur].keys()):
            for dia_y in list(mean_change[dur].keys()):
                for pro in list(mean_bias[dur][dia_x].keys()):
                    if not check_dict_levels(mean_change, dur, dia_y, pro):
                        continue
                    for exp in list(mean_bias[dur][dia_x][pro].keys()):
                        if not check_dict_levels(mean_change, dur, dia_y, pro, exp):
                            continue
                        for dat in list(mean_bias[dur][dia_x][pro][exp].keys()):
                            if not check_dict_levels(mean_change, dur, dia_y, pro, exp, dat):
                                continue
                            # arrays
                            da_x = mean_bias[dur][dia_x][pro][exp][dat]
                            da_y = mean_change[dur][dia_y][pro][exp][dat]
                            # remove datasets
                            l1 = [v for k, v in dict_datasets.items() if "GISS-E2-H" in k or "GISS-E2-R" in k]
                            list_d = [k for k in da_x["dataset"].values if k not in l1]
                            da_x, da_y = da_x.sel(indexers={"dataset": list_d}), da_y.sel(indexers={"dataset": list_d})
                            # correlation and p-value
                            da_r = correlation_pearson(da_x, da_y, dim="dataset", skipna=True)
                            da_p = correlation_pearson_pvalue(da_x, da_y, dim="dataset", skipna=True)
                            # save
                            dia_o = str(dia_x) + "--" + str(dia_y)
                            put_in_dict(corr_map, da_p, dur, dia_o, pro, dat, exp, "p")
                            put_in_dict(corr_map, da_r, dur, dia_o, pro, dat, exp, "r")
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: map
    print("map")
    group, panel = "group_1", "panel_01"
    for dur in list(corr_map.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for dia in list(corr_map[dur].keys()):
            for pro in list(corr_map[dur][dia].keys()):
                list_dat = list(corr_map[dur][dia][pro].keys())
                for ii, dat in enumerate(list_dat):
                    list_exp = sorted(list(corr_map[dur][dia][pro][dat].keys()), key=lambda s: s.lower())
                    for jj, exp in enumerate(list_exp):
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = (dur, dia, pro, grp, pan)
                        print(dur, dia, pro, dat, exp, tuple_k)
                        # panel size
                        tk = tuple_k[:3]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                            for k1, k2 in fig_panel_size.items():
                                if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                    put_in_dict(figure_panel_size, k2, *tk + (k1,))
                        # panel axes definition
                        if check_dict_levels(fig_ticks, "correlation", panel):
                            for k1, k2 in fig_ticks["correlation"][panel].items():
                                if check_dict_levels(figure_axes, *tuple_k + (k1,)) or \
                                        (k1 == "map_c_land" and "ave_pr_val" in dia):
                                    continue
                                tmp = copy__deepcopy(k2)
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
                            val = str(fig_titles[exp]) + ": "
                            val += "1850-" + str(1850 + ny2 - 1) if exp == "historical" else \
                                str(2099 - ny2 + 1) + "-2099"
                            val += " vs. " + str(2005 - ny2 + 1) + "-2005"
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        if jj == 0:
                            put_in_dict(plot_data, str(int(dat[-3:])) + " members", *tuple_k + ("title_row",))
                        # -- map
                        da_p, da_r = corr_map[dur][dia][pro][dat][exp]["p"], corr_map[dur][dia][pro][dat][exp]["r"]
                        plot_type = "map"
                        # shading
                        put_in_dict(plot_data, copy__deepcopy(da_r), *tuple_k + (str(plot_type) + "_c",))
                        val = "right" if exp == list_exp[-1] else None
                        put_in_dict(plot_data, val, *tuple_k + ("legend_position",))
                        # hatching
                        put_in_dict(plot_data, copy__deepcopy(da_p), *tuple_k + (str(plot_type) + "_h",))
                        density = 3
                        kwargs_h = {"hatches": [density * "/"], "levels": [1 - proc_ci / 100, 1]}
                        put_in_dict(plot_data, kwargs_h, *tuple_k + ("kwargs_h",))
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
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for pro in sorted(list(plot_data[dur][dia].keys()), key=lambda s: s.lower()):
                # list plot names
                list_names = list(plot_data[dur][dia][pro].keys())
                # path to output figure
                figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia, pro)))
                figure_name_and_path += "_noGISS-E2"
                # plot
                fig_basic(plot_data[dur][dia][pro], list_names, 1, figure_axes[dur][dia][pro], fig_format,
                          figure_name_and_path, figure_panel_size[dur][dia][pro], fig_title_bool=True,
                          panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e41_change_explained_map_bias_cc_plot()
