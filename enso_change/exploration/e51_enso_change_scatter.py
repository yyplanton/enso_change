# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: scatterplots of ENSO changes (departure from the end of the 20th century) in CMIP5 vs. CMIP6, raw vs. pcc
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from math import ceil as math__ceil
from math import floor as math__floor
from matplotlib.markers import MarkerStyle
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
from enso_change.wrapper.statistics import linear_regression, numpy__array, percentile, student_confidence_interval
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_ts_val_glob"],
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
        "frac": {"x": 0.5, "y": 0.5},
        "panel_01": {"x_delt": 2, "x_size": 8, "y_delt": 2, "y_size": 8},
    },
    # colors
    "fig_colors": {
        "cmip5": "limegreen",
        "cmip6": "r",
        "cmip52": "darkgreen",
        "cmip62": "firebrick",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "o",
        "cmip6": "o",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "std_ts_ano_n34e--pcc": {
            "panel_01": {
                "y_lim": [-50, 50],
                "y_nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                "y_tic": list(range(-50, 51, 25)),
                # "y_lim": None,
                # "y_nam": "change VALUE per CC (UNITS.$^\circ$C$^{-1}$)",
                # "y_tic": None,
            },
        },
        "std_ts_ano_n34e--raw": {
            "panel_01": {
                "x_lim": [-70, 70],
                "x_nam": "VALUE (UNITS)",
                "x_tic": list(range(-70, 71, 35)),
                # "x_lim": None,
                # "x_nam": "change VALUE (UNITS)",
                # "x_tic": None,
            },
        },
    },
    # titles
    "fig_titles": {
        "pcc": "change per GMST",
        "raw": "change",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 40,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e51_enso_change_scatter_plot(
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
    dia2, dur2, pro2 = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat2, exp2 = "multimodel001", "historical"
    da = values_netcdf[dia2][dur2][pro2][exp2][dat2]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    # -- compute change compared to initial conditions
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            if dur1 != dur2:
                continue
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
                        # change per CC or not
                        for war in ["raw", "pcc"]:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if war == "pcc":
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- remove datasets?
    remove_dat = 10
    dict_remove = {
        # 5: ["FGOALS-g2", "GISS-E2-R-CC", "GISS-E2-R_p1", "GISS-E2-R_p2", "GISS-E2-R_p3",
        #     "CAS-ESM2-0"],
        # 10: ["FGOALS-g2", "FIO-ESM", "GISS-E2-H-CC", "GISS-E2-R-CC", "GISS-E2-R_p1", "GISS-E2-R_p2", "GISS-E2-R_p3",
        #      "MIROC-ESM", "MIROC-ESM-CHEM",
        #      "CAS-ESM2-0", "GISS-E2-1-G_p1f1", "MIROC-ES2L"],
        5: ["FIO-ESM", "GISS-E2-H-CC", "GISS-E2-R-CC", "GISS-E2-R_p1", "GISS-E2-R_p2", "GISS-E2-R_p3"],
        10: ["FGOALS-g2", "FIO-ESM", "GISS-E2-H-CC", "GISS-E2-H_p1", "GISS-E2-H_p2", "GISS-E2-H_p3", "GISS-E2-R-CC",
             "GISS-E2-R_p1", "GISS-E2-R_p2", "GISS-E2-R_p3",
             "CAS-ESM2-0", "GISS-E2-1-G_p1f1"],
    }
    if remove_dat in list(dict_remove.keys()):
        for dur in list(mean_change.keys()):
            for dia in list(mean_change[dur].keys()):
                for pro in list(mean_change[dur][dia].keys()):
                    for exp in list(mean_change[dur][dia][pro].keys()):
                        for dat in list(mean_change[dur][dia][pro][exp].keys()):
                            for war in list(mean_change[dur][dia][pro][exp][dat].keys()):
                                # array
                                da = copy__deepcopy(mean_change[dur][dia][pro][exp][dat][war])
                                del mean_change[dur][dia][pro][exp][dat][war]
                                # list datasets
                                list_mod = da["dataset"].values.tolist()
                                # select datasets
                                to_remove = [dict_datasets[k] for k in dict_remove[remove_dat]]
                                list_mod = [k for k in list_mod if k not in to_remove]
                                da = da.sel(indexers={"dataset": list_mod})
                                # save
                                put_in_dict(mean_change, copy__deepcopy(da), dur, dia, pro, exp, dat, war)
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_01"
    for dur in list(mean_change.keys()):
        if not check_dict_levels(mean_change, dur):
            continue
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for dia in list(mean_change[dur].keys()):
            for pro in list(mean_change[dur][dia].keys()):
                list_exp = [k for k in list(mean_change[dur][dia][pro].keys()) if k not in ["piControl", "historical"]]
                for jj, exp in enumerate(list_exp):
                    list_dat = [k for k in list(mean_change[dur][dia][pro][exp].keys()) if k[-3:] not in ["004", "005"]]
                    for ii, dat in enumerate(list_dat):
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = (dur, dia, grp, pan)
                        print(dur, dia, pro, exp, dat, tuple_k)
                        # dictionaries
                        da_x = copy__deepcopy(mean_change[dur][dia][pro][exp][dat]["raw"])
                        da_y = copy__deepcopy(mean_change[dur][dia][pro][exp][dat]["pcc"])
                        x_mm = [math__floor(float(da_x.min().values))] + [math__ceil(float(da_x.max().values))]
                        y_mm = [math__floor(float(da_y.min().values))] + [math__ceil(float(da_y.max().values))]
                        # panel size
                        tk = tuple_k[:-2]
                        if not check_dict_levels(figure_ps, *tk + (pan,)):
                            put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                            for k1, k2 in fig_panel_size.items():
                                if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                    put_in_dict(figure_ps, k2, *tk + (k1,))
                        # panel axes definition
                        for war in list(mean_change[dur][dia][pro][exp][dat].keys()):
                            dd = str(dia) + "--" + str(war)
                            if check_dict_levels(fig_ticks, dd, panel):
                                for k1, k2 in fig_ticks[dd][panel].items():
                                    if check_dict_levels(figure_axes, *tuple_k + (k1,)):
                                        continue
                                    tmp = copy__deepcopy(k2)
                                    if isinstance(tmp, str) and k1 in ["x_nam", "y_nam"]:
                                        tmp = tmp.replace("VALUE", metadata_netcdf[dia]["name_short"])
                                        tmp = tmp.replace("UNITS", metadata_netcdf[dia]["units"])
                                        tmp = tmp.replace(" ()", "")
                                    elif tmp is None and k1 in ["x_lim", "x_tic"]:
                                        tmp = copy__deepcopy(x_mm)
                                    elif tmp is None and k1 in ["y_lim", "y_tic"]:
                                        tmp = copy__deepcopy(y_mm)
                                    if (k1 == "x_nam" and dat != list_dat[-1]) or (k1 == "y_nam" and jj != 0):
                                        tmp = ""
                                    put_in_dict(figure_axes, tmp, *tuple_k + (k1,))
                        dt = {}
                        if check_dict_levels(figure_axes, *tuple_k):
                            dt = read_from_dict(figure_axes, *tuple_k)
                        # test possibility to plot text / legend
                        do_legend = False
                        if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                                dt["y_lim"] is not None:
                            do_legend = True
                        # -- title
                        if ii == 0:
                            put_in_dict(plot_data, fig_titles[exp], *tuple_k + ("title_col",))
                        if jj == 0:
                            val = str(int(dat[-3:])) + " members" if "multimodel" in dat else str(dat)
                            put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                        # -- marker
                        lcf, lcl, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], [], []
                        dict_t = {}
                        for pro2 in ["cmip5", "cmip6"]:
                            # model numbers corresponding to cmip5/6
                            list_mod = sorted(list(set(da_x["dataset"].values.tolist()) &
                                                   set(da_y["dataset"].values.tolist())))
                            list_mod = [k for k in list_mod
                                        if (pro2 == "cmip5" and k < 100) or (pro2 == "cmip6" and k > 100)]
                            da1 = da_x.sel(indexers={"dataset": list_mod}).values
                            da2 = da_y.sel(indexers={"dataset": list_mod}).values
                            dict_t[pro2] = {"x": copy__deepcopy(da1), "y": copy__deepcopy(da2)}
                            da1, da2 = da1.tolist(), da2.tolist()
                            nn = len(da1)
                            # values
                            lcf += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"] * nn
                            lcl += ["k"] * nn
                            lm += [fig_markers[pro2] if pro2 in list(fig_markers.keys()) else "o"] * nn
                            ls += [10] * nn
                            lw += [0.5] * nn
                            lx += copy__deepcopy(da1)
                            ly += copy__deepcopy(da2)
                            lz += [numpy__random__randint(3, 9) for _ in range(nn)]
                        # mean & confidence interval
                        dict_ci = {k1: {k2: student_confidence_interval(v2, axis=0, confidence_level=proc_ci).tolist()
                                        for k2, v2 in dict_t[k1].items()} for k1 in list(dict_t.keys())}
                        overlap_x = 0 if min(dict_ci["cmip6"]["x"]) > max(dict_ci["cmip5"]["x"]) or \
                                         min(dict_ci["cmip5"]["x"]) > max(dict_ci["cmip6"]["x"]) else 1
                        overlap_y = 0 if min(dict_ci["cmip6"]["y"]) > max(dict_ci["cmip5"]["y"]) or \
                                         min(dict_ci["cmip5"]["y"]) > max(dict_ci["cmip6"]["y"]) else 1
                        for pro2 in list(dict_t.keys()):
                            pro3 = str(pro2) + "2"
                            val = fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"
                            lcf += ["none" if k1 == 1 else val for k1 in [overlap_x, overlap_y]]
                            lcl += ["k" for _ in [overlap_x, overlap_y]]
                            lm += [MarkerStyle("D", fillstyle=k1).rotated(deg=45) for k1 in ["bottom", "top"]]
                            ls += [80] * 2
                            lw += [1] * 2
                            lx += [float(dict_t[pro2]["x"].mean())] * 2
                            ly += [float(dict_t[pro2]["y"].mean())] * 2
                            lz += [10] * 2
                        # top / bottom & confidence interval
                        dict_a, dict_aci = {}, {}
                        for thr in ["bot", "top"]:
                            t1 = copy__deepcopy(proc_th) if thr == "bot" else 100 - proc_th
                            for k1 in list(dict_t.keys()):
                                for k2, v2 in dict_t[k1].items():
                                    t2 = percentile(v2, per=t1)
                                    v3 = numpy__array(
                                        [k3 for k3 in v2 if (thr == "bot" and k3 <= t2) or (thr == "top" and k3 >= t2)])
                                    v4 = student_confidence_interval(v3, axis=0, confidence_level=proc_ci)
                                    put_in_dict(dict_a, v3, k1, thr, k2)
                                    put_in_dict(dict_aci, v4, k1, thr, k2)
                            overlap_x = 0 if min(dict_aci["cmip6"][thr]["x"]) > max(dict_aci["cmip5"][thr]["x"]) or \
                                             min(dict_aci["cmip5"][thr]["x"]) > max(dict_aci["cmip6"][thr]["x"]) else 1
                            overlap_y = 0 if min(dict_aci["cmip6"][thr]["y"]) > max(dict_aci["cmip5"][thr]["y"]) or \
                                             min(dict_aci["cmip5"][thr]["y"]) > max(dict_aci["cmip6"][thr]["y"]) else 1
                            for pro2 in list(dict_a.keys()):
                                pro3 = str(pro2) + "2"
                                val = fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"
                                lcf += ["none" if k1 == 1 else val for k1 in [overlap_x, overlap_y]]
                                lcl += ["k" for _ in [overlap_x, overlap_y]]
                                lm += [MarkerStyle("o", fillstyle=k1).rotated(deg=45)
                                       for k1 in ["bottom", "top"]]
                                ls += [80] * 2
                                lw += [1] * 2
                                lx += [float(dict_a[pro2][thr]["x"].mean())] * 2
                                ly += [float(dict_a[pro2][thr]["y"].mean())] * 2
                                lz += [10] * 2
                        pt = "mar"
                        if len(lx) > 0 and len(lx) == len(ly):
                            put_in_dict(plot_data, copy__deepcopy(lcf), *tuple_k + (str(pt) + "_c",))
                            put_in_dict(plot_data, copy__deepcopy(lcl), *tuple_k + (str(pt) + "_lc",))
                            put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(pt) + "_lw",))
                            put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(pt) + "_m",))
                            put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(pt) + "_s",))
                            put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                            put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(pt) + "_y",))
                            put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(pt) + "_z",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        for pro2 in ["cmip5", "cmip6"]:
                            # regression
                            ln = linear_regression(dict_t[pro2]["x"], dict_t[pro2]["y"])
                            lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                            ls += ["-"]
                            lw += [1]
                            lx += [x_mm]
                            ly += [[ln["s"] * k + ln["i"] for k in x_mm]]
                            lz += [20]
                            # confidence interval on the mean
                            if isinstance(do_legend, bool) and do_legend:
                                df = read_from_dict(figure_ps, *tk)
                                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                                da1_c, da2_c = dict_ci[pro2]["x"], dict_ci[pro2]["y"]
                                xxm, yym = float(dict_t[pro2]["x"].mean()), float(dict_t[pro2]["y"].mean())
                                pro3 = str(pro2) + "2"
                                lc += [fig_colors[pro3] if pro3 in list(fig_colors.keys()) else "cyan"] * 6
                                ls += ["-"] * 6
                                lw += [1] * 6
                                lx += [da1_c, [xxm] * 2, [min(da1_c)] * 2, [max(da1_c)] * 2] + \
                                       [[xxm + k * 3 * dx for k in [-1, 1]]] * 2
                                ly += [[yym] * 2, da2_c] + [[yym + k * 3 * dy for k in [-1, 1]]] * 2 + \
                                      [[min(da2_c)] * 2, [max(da2_c)] * 2]
                                lz += [11] * 6
                        if isinstance(do_legend, bool) and do_legend:
                            # x = 0 and y = 0
                            for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                                if k1[0] < 0 < k1[1]:
                                    lc += ["k"]
                                    ls += [":"]
                                    lw += [1]
                                    lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                    ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                    lz += [20]
                            # diagonal
                            lc += ["k"]
                            ls += [":"]
                            lw += [1]
                            lx += [[-max(abs(numpy__array(x_mm + y_mm))), max(abs(numpy__array(x_mm + y_mm)))]]
                            ly += [[-max(abs(numpy__array(x_mm + y_mm))), max(abs(numpy__array(x_mm + y_mm)))]]
                            lz += [20]
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
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            l1 = list(dict_t.keys())
                            for kk, pro2 in enumerate(l1):
                                lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left"]
                                lr += [0]
                                val = str(fig_titles[pro2]) + " " + str(dict_t[pro2]["x"].size)
                                for ll in list(dict_a[pro2].keys())[:1]:
                                    val += "/" + str(dict_a[pro2][ll]["x"].size)
                                lt += [val]
                                lv += ["bottom"]
                                lw += ["regular"]
                                lx += [x1 + 2 * dx]
                                ly += [y1 + ((len(l1) - 1 - kk) * 8 + 2) * dy]
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
    a3 = "_" + str(remove_dat).zfill(2) if remove_dat in list(dict_remove.keys()) else ""
    figure_name = __file__.split("/")[-1].split(".")[0] + "_" + str(a1) + "ci_" + str(a2) + "th" + str(a3)
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dia].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia)))
            # plot
            fig_basic(plot_data[dur][dia], list_names, 1, figure_axes[dur][dia], fig_format,
                      figure_name_and_path, figure_ps[dur][dia], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e51_enso_change_scatter_plot()
