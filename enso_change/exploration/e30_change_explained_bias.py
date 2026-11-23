# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: evolution of variable change in different experiments for models with low / large ENSO change
# Here biases (piControl or historical) are selected
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
from enso_change.wrapper.xarray_base import concatenate
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_pr_val_glob",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw", "ave_pr_val_peso",
    #                       "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw",
    #                       "ave_ts_val_glob",
    #                       "ave_ts_val_aequ", "ave_ts_val_atno", "ave_ts_val_atro",
    #                       "ave_ts_val_iequ", "ave_ts_val_isee", "ave_ts_val_itro", "ave_ts_val_iweq",
    #                       "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl", "ave_ts_val_gxpn",
    #                       "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
    #                       "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi",
    #                       "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
    #                       "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
    #                       "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s",
    #                       "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw", "ave_ts_val_peso",
    #                       "ave_ts_val_pnee", "ave_ts_val_pnew", "ave_ts_val_pnqe",
    #                       "ave_ts_val_ptre", "ave_ts_val_ptro", "ave_ts_val_ptrw", "ave_ts_val_trop"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_pr_val_glob",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw", "ave_pr_val_peso",
    #                       "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw"],
    "data_diagnostics1": ["std_ts_ano_n34e",
                          "ave_ts_val_glob",
                          "ave_ts_val_aequ", "ave_ts_val_atno", "ave_ts_val_atro",
                          "ave_ts_val_iequ", "ave_ts_val_isee", "ave_ts_val_itro", "ave_ts_val_iweq",
                          "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl", "ave_ts_val_gxpn",
                          "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
                          "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi",
                          "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
                          "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
                          "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s",
                          "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw", "ave_ts_val_peso",
                          "ave_ts_val_pnee", "ave_ts_val_pnew", "ave_ts_val_pnqe",
                          "ave_ts_val_ptre", "ave_ts_val_ptro", "ave_ts_val_ptrw", "ave_ts_val_trop"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    # "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["piControl", "historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # "data_projects": ["cmip5", "observations"],
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
    "proc_ci": 100.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 50,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.05, "y": 0.05},
        "panel_1": {"x_delt": 20, "x_size": 80, "y_delt": 10, "y_size": 80},
    },
    # colors
    "fig_colors": {
        "all": "k",
        "piControl--bot": "b",
        "piControl--top": "r",
        "historical--bot": "g",
        "historical--top": "sienna",
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
        "050_year_epoch": {
            "std_ts_ano_n34e--change": {
                "panel_1": {
                    "lim": [-25, 75],
                    "nam": "relative departure\nfrom YEAR1-YEAR2\nVALUE (%)",
                    "tic": list(range(-25, 76, 25)),
                },
            },
            "dataset": {
                "panel_1": {
                    "lab": ["" for _ in range(0, 5, 1)],
                    "lim": [-0.5, 4.5],
                    "nam": "",
                    "tic": list(range(0, 5, 1)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "all": "all models",
        "piControl--bot": "bottom piControl",
        "piControl--top": "top piControl",
        "historical--bot": "bottom historical",
        "historical--top": "top historical",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50,
        "title_col_y": 5,
        "title_row_x": 50,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e30_change_explained_bias_plot(
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
    # concatenate historical and SSPs
    print("concatenate historical and SSPs")
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                if pro == "observations":
                    continue
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
    n_ensemble = {}
    for dia in list(values_json.keys()):
        for dur in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur].keys()):
                if pro == "observations":
                    continue
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_ensembles = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_ensembles = sorted(list(set(list_ensembles)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_ensembles), exp)
        print(exp.ljust(10), "number of ensembles: " + str(len(list_ensembles)).rjust(2))
    #
    # -- Processing
    #
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # -- compute bias
    mean_bias = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp in list(values_json[dia][dur1][pro].keys()):
                    if exp not in ["piControl", "historical"]:
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if check_dict_levels(values_json, dia, dur1, pro) and \
                                len(list({"piControl", "historical"} -
                                         set(list(values_json[dia][dur1][pro].keys())))) != 0:
                            continue
                        # array
                        da = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        # select year
                        yy = int(dur1.split("_")[0]) if exp == "piControl" else (2005 if pro == "cmip5" else 2015)
                        da = da.sel(indexers={"year": yy})
                        # ensemble mean
                        da = float(da.mean(dim="member").values)
                        # save
                        put_in_dict(mean_bias, da, dur1, dia, exp, pro, dat)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia in list(values_json.keys()):
        if dia != data_diagnostics1[0]:
            continue
        for dur1 in list(values_json[dia].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
                    if exp1 == "piControl":
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp1) or \
                                not check_dict_levels(values_json, dia, dur2, pro, exp2) or \
                                (dat == "JRA-55" and data_projects[0] == "cmip5"):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp1][dat])
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia][dur2][pro][exp2][dat])
                        # select experiment epoch
                        yy = da1["year"].values
                        yy = yy[0] if exp1 == "historical" else yy[-1]
                        da1 = da1.sel(indexers={"year": yy}).mean(dim="member")
                        # select reference epoch
                        # yy = 2014 if data_projects[0] == "cmip6" else 2004
                        yy = 2015 if data_projects[0] == "cmip6" else 2005
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        put_in_dict(mean_change, float(da1.values), dur1, dia, exp1, pro, dat)
    print("mean_change", list(mean_change.keys()))
    # -- select top / bottom models
    change_r = {}
    print("select models for all variables")
    for dur1 in list(mean_bias.keys()):
        if dur1 != dur2 or not check_dict_levels(mean_change, dur1):
            continue
        for dia1 in list(mean_bias[dur1].keys()):
            for exp1 in list(mean_bias[dur1][dia1].keys()):
                for pro in list(mean_bias[dur1][dia1][exp1].keys()):
                    for dia2 in list(mean_change[dur1].keys()):
                        for exp2 in list(mean_change[dur1][dia2].keys()):
                            if not check_dict_levels(mean_change, dur1, dia2, exp2, pro):
                                continue
                            # dictionaries
                            d1, d2 = mean_bias[dur1][dia1][exp1][pro], mean_change[dur1][dia2][exp2][pro]
                            # datasets
                            list_datasets = list(set(list(d1.keys())) & set(list(d2.keys())))
                            # regroup biases
                            dict_t = {dat: d1[dat] for dat in list_datasets}
                            # select top / bottom datasets according to bias
                            dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
                            list_datasets = list(dict_t.keys())
                            nn = int(round(len(list_datasets) * proc_th / 100, 0))
                            list_d_bot, list_d_top = list_datasets[:nn], list_datasets[-nn:]
                            # select change values
                            change_all = [d2[dat] for dat in list_datasets]
                            change_bot, change_top = [d2[dat] for dat in list_d_bot], [d2[dat] for dat in list_d_top]
                            # save
                            dia_o = str(dia2) + "--change"
                            if exp1 == list(mean_bias[dur1][dia1].keys())[0]:
                                put_in_dict(change_r, change_all, dur1, dia_o, dia1, exp2, pro, "all")
                            put_in_dict(change_r, change_bot, dur1, dia_o, dia1, exp2, pro, str(exp1) + "--bot")
                            put_in_dict(change_r, change_top, dur1, dia_o, dia1, exp2, pro, str(exp1) + "--top")
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia1 in list(change_r[dur].keys()):
            tmp = []
            for dia2 in list(change_r[dur][dia1].keys()):
                for exp in list(change_r[dur][dia1][dia2].keys()):
                    for pro in list(change_r[dur][dia1][dia2][exp].keys()):
                        for da in change_r[dur][dia1][dia2][exp][pro].values():
                            tmp += copy__deepcopy(da)
            if len(tmp) > 0:
                print(dur, dia1, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                      "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: time series
    print("time series")
    group, panel = "group_1", "panel_1"
    for dur in list(change_r.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for dia1 in list(change_r[dur].keys()):
            list_dia = list(set(list(change_r[dur][dia1].keys())) & set(list(change_r[dur][dia1].keys())))
            list_dia = [k for k in data_diagnostics1 if k in list_dia]
            list_exp = list(set([k2 for k1 in list_dia for k2 in list(change_r[dur][dia1][k1].keys())]))
            list_exp = [k for k in data_experiments if k in list_exp]
            for ii, dia2 in enumerate(list_dia):
                for jj, exp in enumerate(list_exp):
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, dia1, grp, pan)
                    print(dur, dia1, dia2, exp.ljust(10), tuple_k)
                    # panel size
                    tk = tuple_k[:2]
                    if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                        put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                put_in_dict(figure_panel_size, k2, *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], ["dataset", dia1]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, dd, panel):
                            for k1, k2 in fig_ticks[dur][dd][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k1 == "nam" and dd not in ["dataset"]:
                                    d2 = dd.split("--")[0]
                                    # yy = 2014 if data_projects[0] == "cmip6" else 2004
                                    yy = 2015 if data_projects[0] == "cmip6" else 2005
                                    tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
                                    tmp = tmp.replace("VALUE", metadata_json[d2]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[d2]["units"])
                                    tmp = tmp.replace(" ()", "")
                                if k3 == "x_nam" and dia2 != list_dia[-1]:
                                    tmp = ""
                                elif k3 == "y_nam" and jj > 0:
                                    tmp = ""
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
                        val = str(fig_titles[exp]) + " (" + str(len(n_ensemble[exp])) + ")\n"
                        val += "1850-" + str(1850 + ny2 - 1) if exp == "historical" else "2050-" + str(2050 + ny2 - 1)
                        put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, metadata_json[dia2]["name_short"], *tuple_k + ("title_row",))
                    # -- violin
                    la, lbod_fc, lc, lf, lls, llw, lmea, lmea_lw = [], [], [], [], [], [], [], []
                    lmed, lmed_lw, lq, lq_ls, lq_lw, lw, lx, ly, lz = [], [], [], [], [], [], [], [], []
                    if check_dict_levels(change_r, dur, dia1, dia2, exp):
                        for pro in list(change_r[dur][dia1][dia2][exp].keys()):
                            for k2, (sel, da) in enumerate(change_r[dur][dia1][dia2][exp][pro].items()):
                                la += [0.3]
                                lbod_fc += ["k"]
                                lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
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
                                ly += [da]
                                lx += [k2]
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
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # x = 0 and y = 0
                        for n1, k1 in zip(["y"], [dt["y_lim"]]):
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
                    if isinstance(do_legend, bool) and do_legend and jj == 0:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if check_dict_levels(change_r, dur, dia1, dia2, exp):
                            for pro in list(change_r[dur][dia1][dia2][exp].keys()):
                                for k2, (sel, da) in enumerate(change_r[dur][dia1][dia2][exp][pro].items()):
                                    lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
                                    lf += [15]
                                    lh += ["right"]
                                    lt += [fig_titles[sel]]
                                    lv += ["center"]
                                    lw += ["regular"]
                                    lx += [x2 - 2 * dx]
                                    ly += [y2 - (k2 * 8 + 5) * dy]
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
    a1, a2 = str(smile_minimum_size).zfill(2), str(proc_th).zfill(2)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(a1) + "mem", str(a2) + "grp", str(data_projects[0]), str(data_diagnostics1[1].split("_")[1]) + "1"])
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dia].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia)))
            # plot
            fig_basic(plot_data[dur][dia], list_names, 1, figure_axes[dur][dia], fig_format,
                      figure_name_and_path, figure_panel_size[dur][dia], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e30_change_explained_bias_plot()
