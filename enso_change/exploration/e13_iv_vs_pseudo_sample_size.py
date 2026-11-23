# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Test the influence of the pseudo sample size on the error on the magnitude of the internal variability (y)
#    * pseudo sample size (x) vs. relative error [abs(IV_ps - IV_ref) * 100 / IV_ref] (y)
# Pseudo samples are generated as combination from a kernel representing the input distribution
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
from enso_change.wrapper.processors import compute_1var, reader_json
from enso_change.wrapper.xarray_base import apply__kernel, concatenate
from enso_change.wrapper.xarray_tools import statistic, statistic_resampled
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["ave_ts_val_n34e", "std_ts_ano_n34e", "ske_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 30,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 97.5,
    # method: {"combination", "random"}
    "proc_md": "random",
    # number of pseudo samples: int
    "proc_ps": 100000,
    # threshold for classification
    "proc_th": 10,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.005, "y": 0.005},
        "panel_1": {"x_delt": 400, "x_size": 1250, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        # "ACCESS-ESM1-5": "b",
        # "CanESM5-1": "r",
        # "CanESM5_p2": "y",
        # "HadGEM3-GC31-LL": "mediumorchid",
        # "MIROC6": "orange",
        # "MPI-ESM1-2-LR": "limegreen",
        "ACCESS-ESM1-5": "b",
        "CanESM5-1": "r",
        "CanESM5_p2": "y",
        "HadGEM3-GC31-LL": "mediumorchid",
        "IPSL-CM6A-LR": "orange",
        "MIROC-ES2L": "limegreen",
        "MIROC6": "deepskyblue",
        "MPI-ESM1-2-LR": "peru",
        "NorCPM1": "forestgreen",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        # "ACCESS-ESM1-5": "-",
        # "CanESM5-1": "-",
        # "CanESM5_p2": "-",
        # "HadGEM3-GC31-LL": "--",
        # "MIROC6": "--",
        # "MPI-ESM1-2-LR": "--",
        "ACCESS-ESM1-5": "-",
        "CanESM5-1": "-",
        "CanESM5_p2": "-",
        "HadGEM3-GC31-LL": "--",
        "IPSL-CM6A-LR": "--",
        "MIROC-ES2L": "--",
        "MIROC6": "-.",
        "MPI-ESM1-2-LR": "-.",
        "NorCPM1": "-.",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [0, 80],
                    "nam": "absolute relative error (%)",
                    "tic": list(range(0, 81, 20)),
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [0, 80],
                    "nam": "absolute relative error (%)",
                    "tic": list(range(0, 81, 20)),
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [0, 80],
                    "nam": "absolute relative error (%)",
                    "tic": list(range(0, 81, 20)),
                },
            },
            "sample": {
                "panel_1": {
                    "lim": [0, 50],
                    "nam": "pseudo sample size",
                    "tic": list(range(0, 51, 10)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "VALUE",
        "ske_ts_ano_n34e": "VALUE",
        "std_ts_ano_n34e": "VALUE",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 3,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e13_iv_vs_pseudo_sample_size_plot(
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
    n_ensemble = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_experiments = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(10), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print(str().ljust(10), ", ".join(list_experiments))
    #
    # -- Processing
    #
    # compute change compared to initial conditions
    change = {}
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        yy = 1850 + int(dur2.split("_")[0])
                        da2 = float(da2.sel(indexers={"year": yy}).mean(dim="member").values)
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(change, da1, dur1, dia, exp1, dat)
    # -- 'real' internal variability according to existing SMILEs
    sdev_change, _, _ = compute_1var(change, "statistic",
                                     function_arguments={"dim": "member", "stat": "std"})
    # -- estimate internal variability with pseudo samples generated using a kernel
    # iv_error = {}
    # for dur in list(change.keys()):
    #     for dia in list(change[dur].keys()):
    #         for exp in list(change[dur][dia].keys()):
    #             for dat in list(change[dur][dia][exp].keys()):
    #                 if not check_dict_levels(sdev_change, dur, dia, exp, dat):
    #                     continue
    #                 # get array
    #                 da_i = copy__deepcopy(change[dur][dia][exp][dat])
    #                 da_r = copy__deepcopy(sdev_change[dur][dia][exp][dat])
    #                 if da_i["member"].size < 5:
    #                     continue
    #                 # estimate kernel and create data fitting this kernel
    #                 da_k = apply__kernel(da_i, "member", 100)
    #                 # compute internal variability of pseudo samples
    #                 da_o = ()
    #                 list_pss = list(range(5, 51, 5))
    #                 for kk in list_pss:
    #                     # compute std on pseudo samples
    #                     da_t = statistic_resampled(da_k, dim="member", pseudo_sample_size=proc_ps,
    #                                                sample_method=proc_md, sample_size=kk, stat="std")
    #                     # absolute relative error
    #                     da_t = abs(da_t - da_r) * 100 / abs(da_r)
    #                     # confidence interval (1 sided as the absolute error is computed)
    #                     da_t = statistic(da_t, per=proc_ci, stat="per")
    #                     # store value
    #                     da_o += (da_t,)
    #                 # concatenate
    #                 da_o = concatenate(da_o, "ensemble_size")
    #                 da_o = da_o.assign_coords({"ensemble_size": list_pss})
    #                 # save
    #                 put_in_dict(iv_error, da_o, dur, dia, exp, dat)
    iv_error = {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for exp in list(change[dur][dia].keys()):
                for dat in list(change[dur][dia][exp].keys()):
                    if not check_dict_levels(sdev_change, dur, dia, exp, dat):
                        continue
                    # get array
                    da_i = copy__deepcopy(change[dur][dia][exp][dat])
                    da_r = copy__deepcopy(sdev_change[dur][dia][exp][dat])
                    if da_i["member"].size < 5:
                        continue
                    # compute internal variability of pseudo samples
                    da_o = ()
                    list_pss = list(range(4, da_i["member"].size + 1, 2))
                    for kk in list_pss:
                        # compute std on pseudo samples
                        da_t = statistic_resampled(da_i, dim="member", pseudo_sample_size=proc_ps,
                                                   sample_method=proc_md, sample_size=kk, stat="std")
                        # absolute relative error
                        da_t = abs(da_t - da_r) * 100 / abs(da_r)
                        # confidence interval (1 sided as the absolute error is computed)
                        da_t = statistic(da_t, per=proc_ci, stat="per")
                        # store value
                        da_o += (da_t,)
                    # concatenate
                    da_o = concatenate(da_o, "ensemble_size")
                    da_o = da_o.assign_coords({"ensemble_size": list_pss})
                    # save
                    put_in_dict(iv_error, da_o, dur, dia, exp, dat)
    # min-max
    print("change min-max")
    for dur in list(iv_error.keys()):
        for dia in list(iv_error[dur].keys()):
            tmp = []
            for exp in list(iv_error[dur][dia].keys()):
                for da in iv_error[dur][dia][exp].values():
                    tmp += da.values.flatten().tolist()
            if len(tmp) > 0:
                print(dur, dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in data_epoch_lengths:
        ny = int(dur.split("_")[0])
        if dur == dur2:
            continue
        for jj, dia in enumerate(data_diagnostics1):
            pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
            grp = copy__deepcopy(group)
            tuple_k = (dur, grp, pan)
            print(dur, dia, tuple_k)
            # panel size
            tk = tuple_k[:1]
            if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                        put_in_dict(figure_panel_size, k2, *tk + (k1,))
            # panel axes definition
            l1, l2 = ["x_", "y_"], ["sample", dia]
            for vv, dd in zip(l1, l2):
                if check_dict_levels(fig_ticks, dur, dd, panel):
                    for k1, k2 in fig_ticks[dur][dd][panel].items():
                        k3 = str(vv) + str(k1)
                        if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "nam" and dd not in ["percentage", "sample", "time"]:
                            tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                            tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                            tmp = tmp.replace(" ()", "")
                        put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
            dt = {}
            if check_dict_levels(figure_axes, *tuple_k):
                dt = read_from_dict(figure_axes, *tuple_k)
            # test possibility to plot text / legend
            do_legend = False
            if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None \
                    and dt["y_lim"] is not None:
                do_legend = True
            # -- title
            val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
            val = val.replace("UNITS", metadata_json[dia]["units"])
            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            for exp in list(iv_error[dur][dia].keys()):
                for dat, da in iv_error[dur][dia][exp].items():
                    lc += [fig_colors[dat] if dat in list(fig_colors.keys()) else "cyan"]
                    ls += [fig_linestyles[dat] if dat in list(fig_linestyles.keys()) else ":"]
                    lw += [2]
                    lx += [da["ensemble_size"].values]
                    ly += [da.values]
                    lz += [3]
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_panel_size, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # x = 0 and y = 0
                for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                    if k1[0] < 0 < k1[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                        lz += [10]
                # legend
                list_dat = sorted(list(set([dat for exp in list(iv_error[dur][dia].keys())
                                            for dat in list(iv_error[dur][dia][exp].keys())])), key=lambda s: s.lower())
                for cc, dat in enumerate(list_dat):
                    lc += [fig_colors[dat] if dat in list(fig_colors.keys()) else "cyan"]
                    ls += [fig_linestyles[dat] if dat in list(fig_linestyles.keys()) else ":"]
                    lw += [2]
                    lx += [[x2 - 10 * dx, x2 - 2 * dx] if cc < 6 else [x1 + 2 * dx, x1 + 10 * dx]]
                    ly += [[y2 - (cc * 8 + 5) * dy if cc < 6 else y1 + ((len(list_dat) - cc - 1) * 8 + 5) * dy] * 2]
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
            lc, lf, lh, lt, lv, lx, ly = [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_panel_size, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                # legend
                list_dat = sorted(list(set([dat for exp in list(iv_error[dur][dia].keys())
                                            for dat in list(iv_error[dur][dia][exp].keys())])), key=lambda s: s.lower())
                for cc, dat in enumerate(list_dat):
                    lc += [fig_colors[dat] if dat in list(fig_colors.keys()) else "cyan"]
                    lf += [15]
                    lh += ["right" if cc < 6 else "left"]
                    lt += [dat]
                    lv += ["center"]
                    lx += [x2 - 12 * dx if cc < 6 else x1 + 12 * dx]
                    ly += [y2 - (cc * 8 + 5) * dy if cc < 6 else y1 + ((len(list_dat) - cc - 1) * 8 + 5) * dy]
            plot_type = "text"
            if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1 = 100 - (100 - proc_ci) * 2
    if int(a1) == a1:
        a1 = int(a1)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(smile_minimum_size).zfill(2) + "mem", str(a1) + str(proc_md[:2])])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        # list plot names
        list_names = list(plot_data[dur].keys())
        # path to output figure
        figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur)))
        # plot
        fig_basic(plot_data[dur], list_names, 1, figure_axes[dur], fig_format, figure_name_and_path,
                  figure_panel_size[dur], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e13_iv_vs_pseudo_sample_size_plot()
