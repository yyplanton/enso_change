# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Time series of std, internal variability and forced response according to HS09 (1 member / model, 2nd order
# polynomial) and LD20 (LE >= 40 member)
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
from enso_change.wrapper.xarray_base import concatenate, polyval
from enso_change.wrapper.xarray_tools import resample, statistic
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
    "data_experiments": ["historical", "ssp245"],
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
    "proc_ci": 95.,
    # method: {"combination", "random"}
    "proc_md": "combination",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # polynomial order to represent the long-term trend
    "proc_th": {"ave": 4, "std": 3, "ske": 2},
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
        # "panel_1": {"x_delt": 45, "x_size": 3, "y_delt": 100, "y_size": 800},
        "panel_1": {"x_delt": 150, "x_size": 6, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 600, "x_size": 3, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "1m": "turquoise",
        "01m": "turquoise",
        "5m": "y",
        "05m": "y",
        "10m": "darkorchid",
        "30m": "r",
        "35m": "r",
        "40m": "r",
        "50m": "r",
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
                "panel_1": {
                    "lim": [-2, 6],
                    "nam": "departure from 1850-YEAR",
                    "tic": list(range(-2, 7, 2)),
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-2, 6],
                    "nam": "departure from YEAR-2014",
                    "tic": list(range(-2, 7, 2)),
                },
            },
            "ske_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-1, 1],
                    "nam": "departure from 1850-YEAR",
                    "tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-1, 1],
                    "nam": "departure from YEAR-2014",
                    "tic": [round(k / 10, 1) for k in range(-10, 11, 5)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-20, 60],
                    "nam": "departure from 1850-YEAR",
                    "tic": list(range(-20, 61, 20)),
                },
                "panel_2": {
                    "axe": "right",
                    "lim": [-20, 60],
                    "nam": "departure from YEAR-2014",
                    "tic": list(range(-20, 61, 20)),
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1850, 2100],
                    "nam": "Year",
                    "nbr_minor": 9,
                    "tic": list(range(1800, 2100, 100)),
                },
                "panel_2": {
                    "lim": [2000, 2100],
                    "nam": "Year",
                    "nbr_minor": 9,
                    "tic": list(range(2000, 2101, 100)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        # ensemble size
        "1m": "1-mem",
        "3m": "3-mem",
        "5m": "5-mem",
        "01m": "1-mem",
        "03m": "3-mem",
        "05m": "5-mem",
        "10m": "10-mem",
        "30m": "30-mem",
        "35m": "35-mem",
        "40m": "40-mem",
        "50m": "50-mem",
        # uncertainty sources
        "internal": "internal variability",
        "forced": "forced change",
        # methods
        "hs09": "HS09",
        "ld20": "LD20",
        # experiments
        "figure_title": "VALUE (%)",
        **default_arguments["fig_names_exp"],
        **default_arguments["fig_names_reg"],
    },
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {
        "title_col_h": "center",
        "title_col_x": 50 * 1.46,
        "title_col_y": 5,
        "title_row_x": 35,
    },
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e21_hs09_vs_ld20_in_smiles_plot(
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
                    if "historical" not in list_exp or len(list_ssp) == 0:
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
                # delete historical is ssp not available
                d1 = copy__deepcopy(values_json[dia][dur][pro])
                list_dataset = list(set([k2 for k1 in list(d1.keys()) for k2 in list(d1[k1].keys())]))
                for dat in list_dataset:
                    list_exp = [k1 for k1 in list(d1.keys()) if dat in list(d1[k1].keys())]
                    list_ssp = [k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]]
                    if "historical" not in list_exp or len(list_ssp) == 0:
                        for exp in list_exp:
                            del values_json[dia][dur][pro][exp][dat]
                        continue
                    if "historical" in list_exp:
                        del values_json[dia][dur][pro]["historical"][dat]
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
    # -- change or relative change
    # compute change compared to initial conditions
    change = {}
    dur2 = data_epoch_lengths[-1]
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp in list(values_json[dia][dur1][pro].keys()):
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if exp == "piControl" or not check_dict_levels(values_json, dia, dur2, pro, exp):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp][dat])
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        for typ in ["01m", "05m", "10m", str(da1["member"].size) + "m"]:
                            for ref in ["250"]:
                                if exp == "piControl" or (exp == "historical" and ref != "250"):
                                    continue
                                da_o = copy__deepcopy(da1)
                                # select reference epoch
                                yy = 1850 + int(dur2.split("_")[0]) if ref == "250" else 2015
                                da_r = copy__deepcopy(da2.sel(indexers={"year": yy}))
                                # resample
                                if typ in ["1m", "01m"]:
                                    # randomly select 1 member to generate many possible 1 member multimodel mean
                                    da_o, da_r = resample(
                                        da_o, arr_i2=da_r, data_var_i1=dia, data_var_i2=dia, dim="member",
                                        pseudo_sample_size=proc_ps, sample_method="random", sample_size=1)
                                else:
                                    nn = int(typ[:-1])
                                    # randomly select nn members to generate many possible nn member multimodel mean
                                    da_o, da_r = resample(
                                        da_o, arr_i2=da_r, data_var_i1=dia, data_var_i2=dia, dim="member",
                                        pseudo_sample_size=proc_ps, sample_method=proc_md, sample_size=nn)
                                    da_r = da_r.mean(dim="member")
                                # change
                                da_o -= da_r
                                if dia.split("_")[0] in ["std", "var"]:
                                    da_o *= 100 / da_r
                                # save
                                put_in_dict(change, da_o, dur1, dia, typ, ref, exp, dat)
    # internal variability vs. forced change
    list_per = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    uncertainty_m, uncertainty_r = {}, {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for met in ["hs09", "ld20"]:
                for typ in list(change[dur][dia].keys()):
                    if typ in ["1m", "2m", "01m", "02m"] and met == "ld20":
                        continue
                    for ref in list(change[dur][dia][typ].keys()):
                        for exp in list(change[dur][dia][typ][ref].keys()):
                            for dat in list(change[dur][dia][typ][ref][exp].keys()):
                                # array
                                da = copy__deepcopy(change[dur][dia][typ][ref][exp][dat])
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
                                    es = es.expand_dims(dim={"year": em["year"]}, axis=1)
                                    # average across members
                                    if "member" in list(es.dims):
                                        es = es.mean(dim="member")
                                elif met == "ld20":
                                    # -- ld20
                                    # 'forced change' is the multi-member mean
                                    em = da.mean(dim="member")
                                    # 'ensemble spread' is the multi-member std
                                    es = da.std(dim="member")
                                else:
                                    # -- pv26
                                    # 'internal variability' is the residual internal variability
                                    es = da.std(dim="member") / (da["member"].size - 1) ** 0.5
                                    # 'forced change' is the multi-member mean
                                    em = da.mean(dim="member")
                                # -- mean & IQR
                                mea_es, mea_em = es.mean(dim="resample"), em.mean(dim="resample")
                                iqr_es = statistic(es, dim="resample", stat="per", per=[25, 75])
                                iqr_em = statistic(em, dim="resample", stat="per", per=[25, 75])
                                # -- save
                                put_in_dict(uncertainty_m, mea_es, dur, dia, dat, "internal", met, typ)
                                put_in_dict(uncertainty_m, mea_em, dur, dia, dat, "forced", met, typ)
                                put_in_dict(uncertainty_r, iqr_es, dur, dia, dat, "internal", met, typ)
                                put_in_dict(uncertainty_r, iqr_em, dur, dia, dat, "forced", met, typ)
    # min-max
    print("change min-max")
    for dur in list(uncertainty_r.keys()):
        for dia in list(uncertainty_r[dur].keys()):
            tmp = []
            for dat in list(uncertainty_r[dur][dia].keys()):
                for src in list(uncertainty_r[dur][dia][dat].keys()):
                    for met in list(uncertainty_r[dur][dia][dat][src].keys()):
                        for da in uncertainty_r[dur][dia][dat][src][met].values():
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
    for dur in list(uncertainty_m.keys()):
        ny1 = int(dur.split("_")[0])
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        # list diagnostics
        for dia in list(uncertainty_m[dur].keys()):
            # list dataset
            list_dat = sorted(list(uncertainty_m[dur][dia].keys()), key=lambda s: s.lower())
            for ii, dat in enumerate(list_dat):
                list_src = sorted(list(uncertainty_m[dur][dia][dat].keys()), key=lambda s: s.lower())
                for jj, src in enumerate(list_src):
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, dia, grp, pan)
                    print(dur, dia, dat.ljust(17), src.ljust(10), tuple_k)
                    # dictionaries
                    di1 = {}
                    if check_dict_levels(uncertainty_m, dur, dia, dat, src):
                        di1 = read_from_dict(uncertainty_m, dur, dia, dat, src)
                    di2 = {}
                    if check_dict_levels(uncertainty_r, dur, dia, dat, src):
                        di2 = read_from_dict(uncertainty_r, dur, dia, dat, src)
                    # panel size
                    tk = tuple_k[:2]
                    if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                        for k1, d1 in fig_panel_size.items():
                            if k1 == panel and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                                for k2, d2 in d1.items():
                                    if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                        continue
                                    tmp = copy__deepcopy(d2)
                                    if k2 == "x_size" and check_dict_levels(fig_ticks, dur, "time", panel, "lim"):
                                        k3 = read_from_dict(fig_ticks, dur, "time", panel, "lim")
                                        tmp *= (k3[-1] - k3[0])
                                    put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                            elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], ["time", dia]
                    for vv, dd in zip(l1, l2):
                        if check_dict_levels(fig_ticks, dur, dd, panel):
                            for k1, k2 in fig_ticks[dur][dd][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                    tmp = tmp.replace("YEAR", str(1850 + ny2 - 1))
                                    tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                    tmp = tmp.replace(" ()", "")
                                if k3 == "x_nam" and dat != list_dat[-1]:
                                    tmp = ""
                                elif k3 == "y_nam" and jj != 0:
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
                        put_in_dict(plot_data, fig_titles[src], *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, dat, *tuple_k + ("title_row",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    list_met = [k for k in list(di1.keys())]
                    for met in list_met:
                        for typ, da in di1[met].items():
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[met] if met in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [3]
                        for typ, da in di2[met].items():
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"] * 2
                            ls += [fig_linestyles[met] if met in list(fig_linestyles.keys()) else ":"] * 2
                            lw += [0.5]
                            lx += [da["year"].values] * 2
                            ly += [da.min(dim="quantile").values, da.max(dim="quantile").values]
                            lz += [4] * 2
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
                        for ll, met in enumerate(list_met):
                            lc += ["k"]
                            ls += [fig_linestyles[met] if met in list(fig_linestyles.keys()) else ":"]
                            lw += [1]
                            vv = x1 + ((ll % 2) * 33 + (10 if ll < 2 else 2)) * dx
                            lx += [[vv, vv + 8 * dx]]
                            ly += [[y2 - (5 if ll < 2 else 15) * dy] * 2]
                            lz += [1]
                    plot_type = "cur"
                    if len(lx) > 0 and len(lx) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_ls",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                    # -- shading
                    la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
                    for met in list(di2.keys()):
                        for typ, da in di2[met].items():
                            la += [0.3]
                            lc += [fig_colors[typ] if typ in list(fig_colors.keys()) else "cyan"]
                            lx += [da["year"].values]
                            ly1 += [da.min(dim="quantile").values]
                            ly2 += [da.max(dim="quantile").values]
                            lz += [2]
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if len(list(di2.keys())) > 0:
                            la += [0.3]
                            lc += ["k"]
                            vv = x1 + ((len(list_met) % 2) * 33 + (10 if len(list_met) < 2 else 2)) * dx
                            lx += [[vv, vv + 8 * dx]]
                            ly1 += [[y2 - (15 + 3) * dy] * 2]
                            ly2 += [[y2 - (15 - 3) * dy] * 2]
                            lz += [1]
                    plot_type = "sha"
                    if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
                        put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(plot_type) + "_a",))
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(plot_type) + "_y1",))
                        put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(plot_type) + "_y2",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                    # -- text
                    lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_panel_size, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # legend
                        if len(list(di2.keys())) > 0:
                            list_met += ["iqr"]
                        for ll, met in enumerate(list_met):
                            lc += ["k"]
                            lf += [15]
                            lh += ["left"]
                            a1 = int(proc_ci) if proc_ci == int(proc_ci) else copy__deepcopy(proc_ci)
                            lt += ["IQR" if met == "iqr" else fig_titles[met]]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x1 + ((ll % 2) * 33 + (20 if ll < 2 else 12)) * dx]
                            ly += [y2 - (5 if ll < 2 else 15) * dy]
                        # colors
                        list_typ = sorted(list(set([k2 for k1 in list(di1.keys()) for k2 in list(di1[k1].keys())])),
                                          key=lambda s: s.lower())
                        for ll, typ in enumerate(list_typ):
                            lc += [fig_colors[typ]]
                            lf += [15]
                            lh += ["left"]
                            lt += [fig_titles[typ]]
                            lv += ["center"]
                            lw += ["regular"]
                            lx += [x1 + 2 * dx]
                            ly += [y2 - ((ll + 1) * 8 + 15) * dy]
                        # plot title
                        if ii == 0 and jj == 0:
                            # plot size
                            x1_d = fig_panel_size[panel]["x_delt"]
                            x1_m = fig_panel_size[panel]["x_size"]
                            x1_s = read_from_dict(fig_ticks, dur, "time", panel, "lim")
                            x1_s = max(x1_s) - min(x1_s)
                            # add to text
                            lc += ["k"]
                            lf += [20]
                            lh += ["center"]
                            val = copy__deepcopy(fig_titles["figure_title"])
                            val = val.replace("VALUE", metadata_json[dia]["name_short"])
                            lt += [val.replace("UNITS", metadata_json[dia]["units"])]
                            lv += ["center"]
                            lw += ["bold"]
                            a1 = x1_m * x1_s
                            a2 = a1 * len(list_src) + x1_d * (len(list_src) - 1)
                            lx += [x1 + (a2 / (2 * a1)) * (x2 - x1)]
                            ly += [y2 + 20 * dy]
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
                      figure_name_and_path, figure_panel_size[dur][var], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e21_hs09_vs_ld20_in_smiles_plot()
