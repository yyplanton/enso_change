# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: scatterplots of non-convective area (Geng et al. 2024; https://doi.org/10.1038/s41467-024-50156-9) biases
# (values at the end of the 20th century) or changes  (departure from the end of the 20th century) vs. ENSO changes
# (departure from the end of the 20th century)
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
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import linear_regression
from enso_change.wrapper.xarray_base import cf_dim_to_dim, mean_weighted
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_ts_val_glob", "ave_ts_val_tron", "ave_sst_val"],
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
        "end--pcc": {
            "panel_01": {
                "x_lim": None,
                "x_nam": "VALUE (UNITS)",
                "x_tic": None,
            },
        },
        "end--raw": {
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
        "bias": "reference [YEAR1-YEAR2]",
        "end": "end [YEAR1-YEAR2]",
        "change": "change [end minus reference]",
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
def e54_change_convective_scatter_plot(
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
    # -- compute rSST
    dia2 = "ave_ts_val_tron"
    for dia3 in data_diagnostics2[-1:]:
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
        del values_netcdf[dia3]
    for dia in [dia2]:
        if dia in list(values_netcdf.keys()):
            del values_netcdf[dia]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- compute non-convective area
    for dia in list(values_netcdf.keys()):
        if ("ave_rsst_val" in list(values_netcdf.keys()) and dia != "ave_rsst_val") or \
                ("ave_rsst_val" not in list(values_netcdf.keys()) and dia != "ave_rts_val"):
            continue
        dia_o = "ave_nc_val_pequ"
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        # get lat-lon dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # mask outside equatorial Pacific [5S/N, 120-280E]
                        da = select_region(da, [-5, 5], [120, 280])
                        # select non-convective area: 1 where rSST < 0 else 1 (masked outside equatorial Pacific)
                        da_n = da.where(da < 0, other=0).where(da >= 0, other=1).where(da.notnull())
                        # spatial average of equatorial Pacific and non-convective area within equatorial Pacific
                        da_nc = mean_weighted(da_n, dim=[dim_y, dim_x], skipna=True) * 100
                        # save
                        put_in_dict(values_netcdf, da_nc, dia_o, dur, pro, exp, dat)
        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia])
        val = metadata_netcdf[dia]["name_short"].replace("SST", "EP NC area").replace("ts", "EP NC area")
        metadata_netcdf[dia_o]["name_short"] = val
        metadata_netcdf[dia_o]["units"] = "%"
        del values_netcdf[dia]
    print("non-convective area", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- select beginning bias year
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_bias_b = {}
    for dia in list(values_netcdf.keys()):
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
                        put_in_dict(mean_bias_b, da, dur1, pro, dia, dat, "pcc", exp)
                        put_in_dict(mean_bias_b, da, dur1, pro, dia, dat, "raw", exp)
    print("mean_bias beginning", list(mean_bias_b.keys()))
    # -- select end bias year
    mean_bias_e = {}
    for dia in list(values_netcdf.keys()):
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
                        yy = 1850 + int(dur1.split("_")[0]) if exp == "historical" else 2100
                        da = da.sel(indexers={"year": yy})
                        # save
                        put_in_dict(mean_bias_e, da, dur1, pro, dia, dat, "pcc", exp)
                        put_in_dict(mean_bias_e, da, dur1, pro, dia, dat, "raw", exp)
    print("mean_bias end", list(mean_bias_e.keys()))
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
                            if war not in ["historical", "piControl"]:
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, pro, dia1, dat, war, exp1)
    print("mean_change", list(mean_change.keys()))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_01"
    dia_x = "ave_nc_val_pequ"
    for dur in list(mean_bias_b.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for pro in list(mean_bias_b[dur].keys()):
            if not check_dict_levels(mean_bias_e, dur, pro, dia_x) or \
                    not check_dict_levels(mean_change, dur, pro, dia_x):
                continue
            for dia_y in list(mean_change[dur][pro].keys()):
                if not (len(dia_y.split("_")) == 4 and dia_y.split("_")[0] in ["std", "var"]):
                    continue
                ii = -1
                list_typ = ["bias", "end", "change"]
                list_dic = [mean_bias_b[dur][pro][dia_x], mean_bias_e[dur][pro][dia_x], mean_change[dur][pro][dia_x]]
                for typ, d1 in zip(list_typ, list_dic):
                    list_dat = sorted(list(d1.keys()), key=lambda s: s.lower())
                    for dat in list_dat:
                        ii += 1
                        list_war  = sorted(list(d1[dat].keys()), key=lambda s: s.lower())
                        for jj, war in enumerate(list_war):
                            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                            grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                            tuple_k = (dur, pro, dia_y, grp, pan)
                            print(dur, pro, dia_y, typ, dat, war, tuple_k)
                            # dictionaries
                            dix = copy__deepcopy(d1[dat][war])
                            diy = copy__deepcopy(mean_change[dur][pro][dia_y][dat][war])
                            x_mm = [float(da.min().values) for da in dix.values()]
                            x_mm += [float(da.max().values) for da in dix.values()]
                            x1, x2 = min(x_mm), max(x_mm)
                            x_mm = [math__floor(x1), math__ceil(x2)]
                            if x2 - x1 < 5:
                                n1 = 3 if x2 - x1 < 0.05 else (2 if x2 - x1 < 0.5 else 1)
                                n2 = 10**n1
                                x_mm = [round(math__floor(x1 * n2) / n2, n1), round(math__ceil(x2 * n2) / n2, n1)]
                            y_mm = [float(da.min().values) for da in diy.values()]
                            y_mm += [float(da.max().values) for da in diy.values()]
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
                                        elif tmp is None and k1 in ["x_lim", "x_tic"]:
                                            tmp = copy__deepcopy(x_mm)
                                        elif tmp is None and k1 in ["y_lim", "y_tic"]:
                                            tmp = copy__deepcopy(y_mm)
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
                                put_in_dict(plot_data, fig_titles[war], *tuple_k + ("title_col",))
                            if jj == 0:
                                val = str(int(dat[-3:])) + " members" if "multimodel" in dat else str(dat)
                                put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                            # -- marker
                            lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                            for exp in list(dix.keys()):
                                for pro2 in ["cmip5", "cmip6"]:
                                    da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp])
                                    # model numbers corresponding to cmip5/6
                                    list_mod = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                           set(da_y["dataset"].values.tolist())))
                                    list_mod = [k for k in list_mod if
                                                (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                                    da_x = da_x.sel(indexers={"dataset": list_mod}).values
                                    da_y = da_y.sel(indexers={"dataset": list_mod}).values
                                    # nn = len(da_x)
                                    # lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"] * nn
                                    # lm += [fig_markers[pro2] if pro2 in list(fig_markers.keys()) else "o"] * nn
                                    # ls += [80] * nn
                                    # lw += [0.5] * nn
                                    # lx += copy__deepcopy(da_x.tolist())
                                    # ly += copy__deepcopy(da_y.tolist())
                                    # lz += [numpy__random__randint(3, 9) for _ in range(nn)]
                                    # mean
                                    lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                                    lm += [fig_markers[pro2] if pro2 in list(fig_markers.keys()) else "o"]
                                    ls += [150]
                                    lw += [1]
                                    lx += [float(da_x.mean())]
                                    ly += [float(da_y.mean())]
                                    lz += [numpy__random__randint(10, 15)]
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
                                da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp])
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
                                    da_x, da_y = copy__deepcopy(dix[exp]), copy__deepcopy(diy[exp])
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
                                    lv += ["bottom"]
                                    lw += ["regular"]
                                    lx += [x2 - 2 * dx]
                                    ly += [y1 + (kk * 6 + 2) * dy]
                                # type (bias, end, change)
                                if dat == list_dat[0] and jj == 0:
                                    # plot size
                                    y1_d = fig_panel_size[panel]["y_delt"]
                                    y1_s = fig_panel_size[panel]["y_size"]
                                    # add to text
                                    lc += ["k"]
                                    lf += [20]
                                    lh += ["center"]
                                    lr += [90]
                                    val = copy__deepcopy(fig_titles[typ])
                                    val = val.replace("YEAR1", str((2006 if typ == "bias" else 2100) - ny2))
                                    lt += [val.replace("YEAR2", str(2005 if typ == "bias" else 2099))]
                                    lv += ["center"]
                                    lw += ["bold"]
                                    a2 = (y1_s * len(list_dat) + y1_d * (len(list_dat) - 1)) / y1_s
                                    lx += [x1 - 50 * dx]
                                    ly += [y2 - 100 * a2 * dy / 2]
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
    figure_name = __file__.split("/")[-1].split(".")[0]
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for pro in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            for dia_y in sorted(list(plot_data[dur][pro].keys()), key=lambda s: s.lower()):
                # list plot names
                list_names = list(plot_data[dur][pro][dia_y].keys())
                # path to output figure
                figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, pro, dia_y)))
                # plot
                fig_basic(plot_data[dur][pro][dia_y], list_names, 1, figure_axes[dur][pro][dia_y],
                          fig_format, figure_name_and_path, figure_ps[dur][pro][dia_y], fig_title_bool=True,
                          panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e54_change_convective_scatter_plot()
