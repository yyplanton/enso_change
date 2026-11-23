# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Test if piControl's values must be detrended to define piControl's internal variability
#    * ensemble's standard deviation without (x) vs. with (y) correction in piControl
#    * ensemble's standard deviation without (x) and with (x) correction in piControl vs. historical's (y)
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
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.statistics import linear_regression
from enso_change.wrapper.xarray_tools import detrend
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
    "data_epoch_lengths": ["030_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["piControl", "historical"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 10,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": True,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
    # method: {"members", "percent"}
    "proc_md": "percent",
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
        "frac": {"x": 0.5, "y": 0.5},
        "panel_1": {"x_delt": 4, "x_size": 8, "y_delt": 2, "y_size": 8},
    },
    # colors
    "fig_colors": {
        "cor": "r",
        "raw": "b",
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
        "cor": "P",
        "raw": "o",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [0.02, 0.18],
                    "nam": "EXP E$_{\sigma}$ (UNITS)",
                    "tic": [round(k / 100, 2) for k in range(2, 19, 4)],
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [0.10, 0.54],
                    "nam": "EXP E$_{\sigma}$ (UNITS)",
                    "tic": [round(k / 100, 2) for k in range(10, 55, 11)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [3, 23],
                    "nam": "EXP E$_{\sigma}$ (%)",
                    "tic": list(range(3, 24, 5)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "VALUE (UNITS)",
        "ske_ts_ano_n34e": "VALUE",
        "std_ts_ano_n34e": "VALUE (%)",
        "cor": "detrended",
        "raw": "normal",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 4,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e07_iv_correction_plot(
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
        print(exp.ljust(12), "number of ensembles: " + str(len(list_experiments)).rjust(2))
    #
    # -- Processing
    #
    # the first time step(s)
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    for dat in list(values_json[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        del values_json[dia][dur][pro][exp][dat]
                        # select time
                        if exp == "piControl":
                            # piControl does not have 'year', it has a size of one and is the epoch length
                            da = da.isel(indexers={"year": 0})
                            # for the moment 'member' are non overlapping epochs, change member coordinates to have dur
                            # interval
                            da = da.assign_coords({"member": da["member"] * int(dur.split("_")[0])})
                        else:
                            # select the first few time steps to increase the sample size
                            # da = da.sel(indexers={"year": [k for k in da["year"].values if k <= 1900]})
                            da = da.isel(indexers={"year": list(range(5))})
                        # save
                        put_in_dict(values_json, copy__deepcopy(da), dia, dur, pro, exp, dat)
    # -- compute the internal variability
    internal = {}
    exp2 = "piControl"
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp1 in list(values_json[dia][dur][pro].keys()):
                    for dat in list(values_json[dia][dur][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur][pro][exp1][dat])
                        da2 = copy__deepcopy(values_json[dia][dur][pro][exp2][dat])
                        da2 = float(da2.mean(dim="member").values)
                        # correct trend?
                        for typ in ["cor", "raw"]:
                            da_t = copy__deepcopy(da1)
                            # detrend data
                            if typ == "cor" and exp1 == "piControl":
                                da_t = detrend(da_t, deg=1, dim="member")
                            # standard deviation
                            da_t = da_t.std(dim="member")
                            if "year" in list(da_t.dims):
                                da_t = da_t.mean(dim="year")
                            da_t = float(da_t.values)
                            if dia.split("_")[0] in ["std", "var"]:
                                da_t *= 100 / da2
                            # save
                            grp = "detrend or not" if exp1 == "piControl" else "piControl vs. historical"
                            k1 = copy__deepcopy(exp1) if exp1 == "piControl" else copy__deepcopy(typ)
                            k2 = "x" if typ == "raw" and exp1 == "piControl" else "y"
                            put_in_dict(internal, copy__deepcopy(da_t), dur, grp, dia, k1, dat, k2)
                            if exp1 == "piControl":
                                grp = "piControl vs. historical"
                                put_in_dict(internal, copy__deepcopy(da_t), dur, grp, dia, typ, dat, "x")
    # min-max
    print("internal variability min-max")
    for dur in list(internal.keys()):
        dt = {}
        for grp in list(internal[dur].keys()):
            for dia in list(internal[dur][grp].keys()):
                lt = []
                for typ in list(internal[dur][grp][dia].keys()):
                    for dat in list(internal[dur][grp][dia][typ].keys()):
                        lt += [da for da in internal[dur][grp][dia][typ][dat].values()]
                if len(lt) > 0:
                    dt[dia] = copy__deepcopy(lt)
        for dia, lt in dt.items():
            print(dur, dia, "{0:.2f}".format(round(min(lt), 2)).rjust(8), "{0:.2f}".format(round(max(lt), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: scatterplot
    print("scatterplot")
    group, panel = "group_1", "panel_1"
    for dur in list(internal.keys()):
        for ii, lin in enumerate(list(internal[dur].keys())):
            for jj, dia in enumerate(list(internal[dur][lin].keys())):
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, grp, pan)
                print(dur, grp.ljust(8), dia, tuple_k)
                # panel size
                tk = tuple_k[:1]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, k2, *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], [dia, dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, dd, panel):
                        for k1, k2 in fig_ticks[dur][dd][panel].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1.split("_")[-1] == "nam":
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if k3 == "x_nam" and lin == "detrend or not":
                                tmp = tmp.replace("EXP", "normal PI")
                            elif k3 == "x_nam" and lin == "piControl vs. historical":
                                tmp = tmp.replace("EXP", "normal & detrended PI")
                            elif k3 == "y_nam" and lin == "detrend or not":
                                tmp = tmp.replace("EXP", "detrended PI")
                            elif k3 == "y_nam" and lin == "piControl vs. historical":
                                tmp = tmp.replace("EXP", "historical")
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
                if ii == 0:
                    val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                    val = val.replace("UNITS",  metadata_json[dia]["units"])
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, lin, *tuple_k + ("title_row",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                for typ in list(internal[dur][lin][dia].keys()):
                    xx, yy = [], []
                    for dat in list(internal[dur][lin][dia][typ].keys()):
                        if "x" not in list(internal[dur][lin][dia][typ][dat].keys()) or \
                                "y" not in list(internal[dur][lin][dia][typ][dat].keys()):
                            continue
                        xx.append(internal[dur][lin][dia][typ][dat]["x"])
                        yy.append(internal[dur][lin][dia][typ][dat]["y"])
                    ln = linear_regression(xx, yy)
                    lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "k"]
                    ls += ["-"]
                    lw += [1]
                    lx += [[min(xx), max(xx)]]
                    ly += [[ln["s"] * k + ln["i"] for k in [min(xx), max(xx)]]]
                    lz += [20]
                if isinstance(do_legend, bool) and do_legend:
                    # diagonal
                    lc += ["k"]
                    ls += [":"]
                    lw += [1]
                    lx += [[min(dt["x_lim"] + dt["y_lim"]), max(dt["x_lim"] + dt["y_lim"])]]
                    ly += [[min(dt["x_lim"] + dt["y_lim"]), max(dt["x_lim"] + dt["y_lim"])]]
                    lz += [10]
                plot_type = "cur"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- marker
                lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                for typ in list(internal[dur][lin][dia].keys()):
                    for dat in list(internal[dur][lin][dia][typ].keys()):
                        if "x" not in list(internal[dur][lin][dia][typ][dat].keys()) or \
                                "y" not in list(internal[dur][lin][dia][typ][dat].keys()):
                            continue
                        lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "grey"]
                        lm += [fig_markers[typ] if typ in list(fig_markers.keys()) else "o"]
                        ls += [40]
                        lw += [0.5]
                        lx += [internal[dur][lin][dia][typ][dat]["x"]]
                        ly += [internal[dur][lin][dia][typ][dat]["y"]]
                        lz += [2]
                plot_type = "mar"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- text
                lc, lf, lh, lr, lt, lv, lx, ly = [], [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    for kk, typ in enumerate(list(internal[dur][lin][dia].keys())):
                        lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "grey"]
                        lf += [15]
                        lh += ["right"]
                        lr += [0]
                        lt += [fig_titles[typ]]
                        lv += ["center"]
                        lx += [x2 - 2 * dx]
                        ly += [y1 + (kk * 8 + 5) * dy]
                plot_type = "text"
                if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (plot_type,))
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(plot_type) + "_fs",))
                    put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(plot_type) + "_ha",))
                    put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(plot_type) + "_r",))
                    put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(plot_type) + "_va",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", str(smile_minimum_size).zfill(2) + "mem"])
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
e07_iv_correction_plot()
