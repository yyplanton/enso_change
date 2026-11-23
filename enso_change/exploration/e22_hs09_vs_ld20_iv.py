# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Ratio of internal variability estimated using HS09 vs. LD20 as a function of ensemble size
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
from enso_change.plots.panels import default_plot
from enso_change.plots.templates import fig_basic
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict, read_from_dict
from enso_change.wrapper.processors import reader_json
from enso_change.wrapper.xarray_base import numpy__array, concatenate, polyval
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["std_ts_ano_n34e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["piControl", "historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": {},  # default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 5,
    # list of rejected SMILEs: list[str]
    "smile_rejected": default_arguments["smile_rejected"],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
    # method: {"combination", "random"}
    "proc_md": "combination",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # polynomial order to represent the long-term trend
    "proc_th": {"ave": 4, "std": 2, "ske": 2},
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
        "panel_1": {"x_delt": 4, "x_size": 12, "y_delt": 2, "y_size": 8},
    },
    # colors
    "fig_colors": {
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "hs09": "-",
        "ld20": "--",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_ts_val_n34e": {
            },
            "ske_ts_ano_n34e": {
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-50, 50],
                    "nam": "HS09 / LD20 (%)",
                    "tic": list(range(-50, 51, 25)),
                },
            },
            "ensemble": {
                "panel_1": {
                    "lim": [0, 60],
                    "nam": "ensemble size",
                    "tic": list(range(0, 61, 10)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_ts_val_n34e": "VALUE: NBR polynomial vs. E$_{\sigma}$",
        "ske_ts_ano_n34e": "VALUE: NBR polynomial vs. E$_{\sigma}$",
        "std_ts_ano_n34e": "VALUE: NBR polynomial vs. E$_{\sigma}$",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e22_hs09_vs_ld20_iv_plot(
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
    # concatenate historical and SSPs
    print("concatenate historical and SSPs")
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                # dictionary
                d1 = values_json[dia][dur][pro]
                # list datasets
                list_dataset = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                for dat in list_dataset:
                    # list experiments
                    list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                    list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                    # delete datasets
                    if "piControl" not in list_exp or "historical" not in list_exp:
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                    # concatenate historical and SSPs
                    for exp in list_ssp:
                        # arrays
                        da_h = copy__deepcopy(values_json[dia][dur][pro]["historical"][dat])
                        da_s = copy__deepcopy(values_json[dia][dur][pro][exp][dat])
                        del values_json[dia][dur][pro][exp][dat]
                        # list members in both historical and ssp
                        lm = sorted(list(set(da_h["member"].values.tolist()) & set(da_s["member"].values.tolist())))
                        # select members
                        da_h, da_s = da_h.sel(indexers={"member": lm}), da_s.sel(indexers={"member": lm})
                        # concatenate along year axis
                        da_o = concatenate((da_h, da_s), "year")
                        # resave
                        put_in_dict(values_json, copy__deepcopy(da_o), dia, dur, pro, exp, dat)
                # delete experiment if no dataset available
                for exp in list(values_json[dia][dur][pro].keys()):
                    if len(list(values_json[dia][dur][pro][exp].keys())) == 0:
                        del values_json[dia][dur][pro][exp]
    n_ens = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ens, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        if exp not in list(n_ens.keys()):
            continue
        list_experiments = copy__deepcopy(n_ens[exp])
        del n_ens[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ens, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(10), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
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
                        if not check_dict_levels(values_json, dia, dur2, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] == "ssp":
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select reference epoch
                        da2 = da2.sel(indexers={"year": 1850 + int(dur2.split("_")[0])}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(change, da1, dur1, dia, pro, exp1, dat)
    # internal variability
    internal = {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for pro in list(change[dur][dia].keys()):
                for exp1 in list(change[dur][dia][pro].keys()):
                    for dat in list(change[dur][dia][pro][exp1].keys()):
                        if not check_dict_levels(change, dur, dia, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(change[dur][dia][pro][exp1][dat])
                        da2 = copy__deepcopy(change[dur][dia][pro][exp2][dat]) if exp1 == "piControl" else \
                            copy__deepcopy(da1)
                        l1 = ["ld20", "ld20"] if exp1 == "piControl" else ["hs09", "ld20"]
                        dd = {}
                        for kk, (met, da) in enumerate(zip(l1, [da1, da2])):
                            if met == "hs09":
                                # -- hs09
                                # 'forced change' is a 2nd order polynomial fit: compute coefficient then fit
                                em = copy__deepcopy(da)
                                if "member" in list(em.dims):
                                    em = em.mean(dim="member")
                                coef = em.polyfit("year", proc_th[dia.split("_")[0]])
                                em = polyval(da["year"], coef["polyfit_coefficients"], degree_dim="degree")
                                # 'internal variability' is the temporal std of the residual
                                es = (da - em).std(dim="year")
                                # average across members
                                if "member" in list(es.dims):
                                    es = es.mean(dim="member")
                            elif met == "ld20":
                                # -- ld20
                                # 'ensemble spread' is the multi-member std
                                es = da.std(dim="member")
                            else:
                                # -- pv26
                                # 'internal variability' is the residual internal variability
                                es = da.std(dim="member") / (da["member"].size - 1) ** 0.5
                            if "year" in list(es.dims):
                                es = es.mean(dim="year")
                            dd[kk] = copy__deepcopy(es)
                        # -- ratio
                        da1, da2 = dd[0], dd[1]
                        da_o = (da1 - da2) * 100 / da2
                        # -- save
                        put_in_dict(internal, da_o, dur, dia, pro, dat, exp1)
    # min-max
    print("change min-max")
    for dur in list(internal.keys()):
        for dia in list(internal[dur].keys()):
            tmp = []
            for pro in list(internal[dur][dia].keys()):
                for dat in list(internal[dur][dia][pro].keys()):
                    for da in internal[dur][dia][pro][dat].values():
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
    for dur in list(internal.keys()):
        for dia in list(internal[dur].keys()):
            pan = copy__deepcopy(panel)
            grp = copy__deepcopy(group)
            tuple_k = (dur, dia, grp, pan)
            print(dur, dia, tuple_k)
            # dictionary
            di1 = {}
            if check_dict_levels(internal, dur, dia):
                di1 = read_from_dict(internal, dur, dia)
            # panel size
            tk = tuple_k[:2]
            if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                for k1, k2 in fig_panel_size.items():
                    if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                        put_in_dict(figure_panel_size, k2, *tk + (k1,))
            # panel axes definition
            l1, l2 = ["x_", "y_"], ["ensemble", dia]
            for vv, dd in zip(l1, l2):
                if check_dict_levels(fig_ticks, dur, dd, panel):
                    for k1, k2 in fig_ticks[dur][dd][panel].items():
                        k3 = str(vv) + str(k1)
                        if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                            continue
                        tmp = copy__deepcopy(k2)
                        if isinstance(tmp, str) and k1 == "nam" and dd not in ["ensemble"]:
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
            nn = proc_th[dia.split("_")[0]]
            nn = str(nn) + "$^{" + ("st" if nn == 1 else ("nd" if nn == 2 else ("rd" if nn == 3 else "th"))) +"}$"
            put_in_dict(plot_data, val.replace("NBR", nn), *tuple_k + ("title_col",))
            # -- curve
            lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                # x = 0 and y = 0
                for k1, k2 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                    if k2[0] < 0 < k2[1]:
                        lc += ["k"]
                        ls += [":"]
                        lw += [1]
                        lx += [[0] * 2 if k1 == "x" else dt["x_lim"]]
                        ly += [dt["y_lim"] if k1 == "x" else [0] * 2]
                        lz += [1]
            plot_type = "cur"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- marker
            lc, llc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], [], []
            dict_t = {}
            for pro in list(di1.keys()):
                for dat in list(di1[pro].keys()):
                    for exp1, da in di1[pro][dat].items():
                        if check_dict_levels(change, dur, dia, pro, exp2, dat):
                            ensemble_size = change[dur][dia][pro][exp2][dat]["member"].size
                            lc += [fig_colors[exp1] if exp1 in list(fig_colors.keys()) else "grey"]
                            llc += ["k"]
                            lm += ["o"]
                            ls += [40]
                            lw += [0.5]
                            lx += [ensemble_size]
                            ly += [float(da.values)]
                            lz += [numpy__random__randint(2, 20)]
                            put_in_dict(dict_t, [float(da.values)], ensemble_size)
            for kk, vv in dict_t.items():
                if isinstance(vv, list) and len(vv) > 1:
                    lc += ["k"]
                    llc += ["k"]
                    lm += ["D"]
                    ls += [50]
                    lw += [0.5]
                    lx += [kk]
                    ly += [numpy__array(vv).mean()]
                    lz += [21]
            plot_type = "mar"
            if len(lx) > 0 and len(lx) == len(ly):
                put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                put_in_dict(plot_data, copy__deepcopy(llc), *tuple_k + (str(plot_type) + "_lc",))
                put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
            # -- text
            lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
            if isinstance(do_legend, bool) and do_legend:
                df = read_from_dict(figure_panel_size, *tk)
                x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                for kk, exp in enumerate(data_experiments):
                    lc += [fig_colors[exp]]
                    lf += [15]
                    lh += ["left"]
                    lt += [str(fig_titles[exp]) +
                           (" (" + str(len(n_ens[exp])) + ")" if exp in list(n_ens.keys()) else "")]
                    lv += ["center"]
                    lw += ["regular"]
                    lx += [x1 + ((kk % 3) * 54 + 2) * dx]
                    ly += [y2 - 5 * dy if kk < 3 else y1 + 5 * dy]
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
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", proc_md, str(a1) + "mem", str(a2) + "po"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for var in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][var].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, var)))
            # plot
            fig_basic(plot_data[dur][var], list_names, 1, figure_axes[dur][var], fig_format,
                      figure_name_and_path, figure_panel_size[dur][var], fig_title_bool=False, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e22_hs09_vs_ld20_iv_plot()
