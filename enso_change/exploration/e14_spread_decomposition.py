# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Time series of mean, std, ske changes uncertainty source (as in HS09)
#    * internal variability
#    * model
#    * scenario
# 3 methods are compared:
#  1) HS09 (or BB21): 1 member / model, 2nd order polynomial
#  2) LD20: SMILEs (I chose SMILE >= 5 member)
#  3) PV2?: SMILEs (I chose SMILE >= 5 member) mean -> residual IV and model spread - residual IV
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
from enso_change.wrapper.xarray_base import align, concatenate, create_array, polyval
from enso_change.wrapper.xarray_tools import statistic
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
    "data_experiments": ["historical", "ssp126", "ssp245", "ssp370", "ssp585"],
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
    # method: {"250", "150", "100"}
    "proc_md": "250",
    # number of pseudo samples: int
    "proc_ps": 100000,
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
        "frac": {"x": 0.005, "y": 0.005},
        "panel_1": {"x_delt": 400, "x_size": 5, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "all": "darkviolet",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "hs09": "-",
        "ld20": "--",
        "pv26": ":",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "100": {
                "ave_ts_val_n34e": {
                    "panel_1": {
                        "lim": [0, 1.2],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 10, 1) for k in range(0, 13, 3)],
                    },
                },
                "ske_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 0.48],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 100, 2) for k in range(0, 49, 12)],
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 36],
                        "nam": "uncertainty (%)",
                        "tic": list(range(0, 37, 9)),
                    },
                },
                "time": {
                    "panel_1": {
                        "lim": [2000, 2100],
                        "nam": "Year",
                        "tic": list(range(2000, 2101, 50)),
                    },
                },
            },
            "150": {
                "ave_ts_val_n34e": {
                    "panel_1": {
                        "lim": [0, 1.2],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 10, 1) for k in range(0, 13, 3)],
                    },
                },
                "ske_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 0.52],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 100, 2) for k in range(0, 53, 13)],
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 32],
                        "nam": "uncertainty (%)",
                        "tic": list(range(0, 33, 8)),
                    },
                },
                "time": {
                    "panel_1": {
                        "lim": [1950, 2100],
                        "nam": "Year",
                        "tic": list(range(1950, 2101, 50)),
                    },
                },
            },
            "250": {
                "ave_ts_val_n34e": {
                    "panel_1": {
                        "lim": [0, 1.2],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 10, 1) for k in range(0, 13, 3)],
                    },
                },
                "ske_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 0.6],
                        "nam": "uncertainty (UNITS)",
                        "tic": [round(k / 100, 2) for k in range(0, 61, 15)],
                    },
                },
                "std_ts_ano_n34e": {
                    "panel_1": {
                        "lim": [0, 40],
                        "nam": "uncertainty (%)",
                        "tic": list(range(0, 41, 10)),
                    },
                },
                "time": {
                    "panel_1": {
                        "lim": [1850, 2100],
                        "nam": "Year",
                        "tic": list(range(1850, 2101, 50)),
                    },
                },
            },
        },
    },
    # titles
    "fig_titles": {
        # diagnostics
        "ave_ts_val_n34e": "VALUE (UNITS)",
        "ske_ts_ano_n34e": "VALUE",
        "std_ts_ano_n34e": "VALUE (%)",
        # uncertainty sources
        "internal": "internal variability",
        "model": "model",
        "scenario": "scenario",
        # methods
        "hs09": "HS09",
        "ld20": "LD20",
        "pv26": "ours",
        # experiments
        "all": "all exp.",
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
def e14_spread_decomposition_plot(
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
        proc_th: dict[str, int] = None,
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
                    if "historical" not in list_exp:
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
    # compute change compared to initial conditions
    change = {}
    dur2 = data_epoch_lengths[-1]
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dur1 == dur2:
                continue
            for pro in list(values_json[dia][dur1].keys()):
                for exp in list(values_json[dia][dur1][pro].keys()):
                    if exp == "piControl" or (exp == "historical" and proc_md != "250"):
                        continue
                    for dat in list(values_json[dia][dur1][pro][exp].keys()):
                        if not check_dict_levels(values_json, dia, dur2, pro, exp):
                            continue
                        ny1 = int(dur1.split("_")[0])
                        ny2 = int(dur2.split("_")[0])
                        # array
                        da1 = copy__deepcopy(values_json[dia][dur1][pro][exp][dat])
                        da2 = copy__deepcopy(values_json[dia][dur2][pro][exp][dat])
                        # select reference epoch
                        yy = 1850 + ny2 if proc_md == "250" else (2000 if proc_md == "150" else 2015)
                        da2 = da2.sel(indexers={"year": yy})
                        # select epochs
                        if proc_md in ["100", "150"]:
                            da1 = da1.sel(indexers={"year": [k for k in da1["year"] if k >= 2100 - int(proc_md)]})
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - ny1 / 2})
                        for typ in ["val", "mea"]:
                            da_o = copy__deepcopy(da1.isel(indexers={"member": 0}) if typ == "val" else da1)
                            da_r = float(da2.isel(indexers={"member": 0}).values if typ == "val" else
                                         da2.mean(dim="member").values)
                            # change
                            da_o -= da_r
                            if dia.split("_")[0] in ["std", "var"]:
                                da_o *= 100 / da_r
                            # save
                            put_in_dict(change, da_o, dur1, dia, typ, exp, dat)
    # uncertainty source (internal variability, model, scenario)
    percentiles_to_compute = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    uncertainty_m, uncertainty_r = {}, {}
    for dur in list(change.keys()):
        for dia in list(change[dur].keys()):
            for met in ["hs09", "ld20", "pv26"]:
                typ = "val" if met == "hs09" else "mea"
                if typ not in list(change[dur][dia].keys()):
                    continue
                ny = int(dur.split("_")[0])
                d1 = copy__deepcopy(change[dur][dia][typ])
                dict_em = {}
                for exp in list(d1.keys()):
                    # array (concatenate datasets)
                    da = concatenate(tuple(da for da in d1[exp].values()), "dataset")
                    if met == "hs09" and typ == "val":
                        # -- hs09
                        # 'ensemble mean' is a 2nd order polynomial fit: compute coefficient then fit
                        coef = da.polyfit("year", proc_th[dia.split("_")[0]])
                        em = polyval(da["year"], coef["polyfit_coefficients"], degree_dim="degree")
                        # 'ensemble spread' is the temporal std of the residual from the 'ensemble mean'
                        es = (da - em).std(dim="year")
                        es = es.expand_dims(dim={"year": em["year"]}, axis=1)
                    elif met == "ld20" and typ == "mea":
                        # -- ld20
                        # 'ensemble mean' is the multi-member mean
                        em = da.mean(dim="member")
                        # 'ensemble spread' is the multi-member std
                        es = da.std(dim="member")
                    else:
                        # -- pv26
                        # 'ensemble spread' is the uncertainty on the multi-member mean
                        es = da.std(dim="member") / (da["member"].size - 1) ** 0.5
                        # 'ensemble mean' is the multi-member mean
                        em = da.mean(dim="member")
                    # -- uncertainty sources
                    # 'model uncertainty' is the multi-model-mean std
                    mod_m = em.std(dim="dataset")
                    if met == "pv26" and typ == "mea":
                        # spread still contains some internal variability residual as the ensemble size is not infinit
                        mod_m -= es.mean(dim="dataset")
                        mod_m = mod_m.where(mod_m >= 0, other=0)
                    # 'internal variability' is the multi-model-mean
                    iva_m = es.mean(dim="dataset")
                    iva_r = statistic(es, dim="dataset", stat="per", per=percentiles_to_compute)
                    # -- save
                    put_in_dict(uncertainty_m, mod_m, dur, dia, "model", met, exp)
                    put_in_dict(uncertainty_m, iva_m, dur, dia, "internal", met, exp)
                    put_in_dict(uncertainty_r, iva_r, dur, dia, "internal", met, exp)
                    if exp not in ["piControl", "historical"]:
                        dict_em[exp] = copy__deepcopy(em)
                # -- scenario spread
                if not check_dict_levels(uncertainty_r, dur, dia, "internal", met) or met == "pv26":
                    continue
                list_ssp = [exp for exp in list(uncertainty_r[dur][dia]["internal"][met].keys()) if exp[:3] == "ssp"]
                if len(list_ssp) == 0 or len(list(dict_em.keys())) < 2:
                    continue
                # create arrays of zeros
                sce_m = create_array(uncertainty_m[dur][dia]["internal"][met][list_ssp[0]], data_var="scenario_mean")
                sce_r = create_array(uncertainty_r[dur][dia]["internal"][met][list_ssp[0]], data_var="scenario_range")
                # concatenate experiments (datasets mean)
                da_t1 = concatenate(tuple(da.mean(dim="dataset") for da in dict_em.values()), "experiment")
                # concatenate datasets (experiments spread)
                list_dat = sorted(list(set([dat for exp in list(dict_em.keys()) for dat in list(d1[exp].keys())])),
                                  key=lambda s: s.lower())
                da_t2 = ()
                for dat in list_dat:
                    list_exp = [exp for exp in list(dict_em.keys()) if dat in list(d1[exp].keys())]
                    if len(list_exp) < 3:
                        continue
                    # concatenate experiments
                    tmp = tuple(d1[exp][dat].mean(dim="member") if "member" in list(d1[exp][dat].dims) else
                                d1[exp][dat] for exp in list_exp)
                    tmp = concatenate(tmp, "experiment")
                    # store experiments spread
                    da_t2 += (tmp.std(dim="experiment"),)
                da_t2 = concatenate(da_t2, "dataset")
                # select scenario epochs
                da_t1 = da_t1.sel(indexers={"year": [k for k in da_t1["year"] if k >= 2015 - ny / 2]})
                da_t2 = da_t2.sel(indexers={"year": [k for k in da_t2["year"] if k >= 2015 - ny / 2]})
                # 'scenario' is the multi-model-mean intermodel std
                sce_m += align(sce_m, da_t1.std(dim="experiment"), fill_value=0, join="outer")[1]
                da_t2 = statistic(da_t2, dim="dataset", stat="per", per=percentiles_to_compute)
                print(dur, dia, met, "sce_r", sce_r.shape, sce_r.dims, "da_t2", da_t2.shape, da_t2.dims)
                sce_r += align(sce_r, da_t2, fill_value=0, join="outer")[1]
                # -- save
                put_in_dict(uncertainty_m, sce_m, dur, dia, "scenario", met, "all")
                put_in_dict(uncertainty_r, sce_r, dur, dia, "scenario", met, "all")
    # min-max
    print("change min-max")
    for dur in list(uncertainty_m.keys()):
        for dia in list(uncertainty_m[dur].keys()):
            tmp = []
            for src in list(uncertainty_m[dur][dia].keys()):
                for met in list(uncertainty_m[dur][dia][src].keys()):
                    for exp, da in uncertainty_m[dur][dia][src][met].items():
                        tmp += da.values.flatten().tolist()
                        if check_dict_levels(uncertainty_r, dur, dia, src, met, exp):
                            tmp += uncertainty_r[dur][dia][src][met][exp].values.flatten().tolist()
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
        list_n = ["internal", "model", "scenario"]
        for ii, src in enumerate(list_n):
            for jj, dia in enumerate(data_diagnostics1):
                ny = int(dur.split("_")[0])
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, grp, pan)
                print(dur, src.ljust(8), dia, tuple_k)
                # panel size
                tk = tuple_k[:1]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    for k1, d1 in fig_panel_size.items():
                        if k1[:6] == "panel_" and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k2, d2 in d1.items():
                                if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                    continue
                                tmp = copy__deepcopy(d2)
                                if k2 == "x_size" and \
                                        check_dict_levels(fig_ticks, dur, proc_md, "time", panel, "lim"):
                                    k3 = read_from_dict(fig_ticks, dur, proc_md, "time", panel, "lim")
                                    tmp *= (k3[-1] - k3[0])
                                put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                        elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time", dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, proc_md, dd, panel):
                        for k1, k2 in fig_ticks[dur][proc_md][dd][panel].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd != "time":
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if k3 == "x_nam" and src != list_n[-1]:
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
                    val = str(metadata_json[dia]["name_short"])
                    val += "" if dia.split("_")[0] == "ske" else (" (%)" if dia.split("_")[0] in ["std", "var"] else
                                                                  " (" + str(metadata_json[dia]["units"]) + ")")
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    put_in_dict(plot_data, fig_titles[src], *tuple_k + ("title_row",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                if check_dict_levels(uncertainty_m, dur, dia, src):
                    for met in list(uncertainty_m[dur][dia][src].keys()):
                        for exp, da in uncertainty_m[dur][dia][src][met].items():
                            if proc_md == "250" and exp[:3] == "ssp":
                                da = da.sel(indexers={"year": [k for k in da["year"].values if k > 2015 - ny / 2]})
                            lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[met] if met in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [da["year"].values]
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
                    if check_dict_levels(uncertainty_m, dur, dia, src):
                        list_met = [k for k in list(uncertainty_m[dur][dia][src].keys())]
                        for kk, met in enumerate(list_met):
                            lc += ["k"]
                            ls += [fig_linestyles[met] if met in list(fig_linestyles.keys()) else ":"]
                            lw += [1]
                            lx += [[x1 + (kk * 33 + 10) * dx, x1 + (kk * 33 + 18) * dx]]
                            ly += [[y2 - 5 * dy] * 2]
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
                if check_dict_levels(uncertainty_r, dur, dia, src):
                    for met in list(uncertainty_r[dur][dia][src].keys()):
                        for exp, da in uncertainty_r[dur][dia][src][met].items():
                            if proc_md == "250" and exp[:3] == "ssp":
                                da = da.sel(indexers={"year": [k for k in da["year"].values if k > 2015 - ny / 2]})
                            la += [0.3]
                            lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
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
                    if check_dict_levels(uncertainty_r, dur, dia, src):
                        la += [0.3]
                        lc += ["k"]
                        lx += [[x1 + (3 * 33 + 10) * dx, x1 + (3 * 33 + 18) * dx]]
                        ly1 += [[y2 - (5 + 3) * dy] * 2]
                        ly2 += [[y2 - (5 - 3) * dy] * 2]
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
                lc, lf, lh, lt, lv, lx, ly = [], [], [], [], [], [], []
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    if check_dict_levels(uncertainty_m, dur, dia, src):
                        list_met = [k for k in list(uncertainty_m[dur][dia][src].keys())]
                        for kk, met in enumerate(list_met):
                            # line styles
                            lc += ["k"]
                            lf += [15]
                            lh += ["left"]
                            lt += [fig_titles[met]]
                            lv += ["center"]
                            lx += [x1 + (kk * 33 + 20) * dx]
                            ly += [y2 - 5 * dy]
                    if check_dict_levels(uncertainty_r, dur, dia, src):
                        # shading
                        lc += ["k"]
                        lf += [15]
                        lh += ["left"]
                        val = [str(int(k)) if k == int(5) else str(k) for k in percentiles_to_compute]
                        lt += ["-".join(val) + "%"]
                        lv += ["center"]
                        lx += [x1 + (3 * 33 + 20) * dx]
                        ly += [y2 - 5 * dy]
                    if check_dict_levels(uncertainty_m, dur, dia, src):
                        list_exp = sorted(list(set([k2 for k1 in list(uncertainty_m[dur][dia][src].keys())
                                                    for k2 in list(uncertainty_m[dur][dia][src][k1].keys())])),
                                          key=lambda s: s.lower())
                        for kk, exp in enumerate(list_exp):
                            # colors
                            lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["left"]
                            lt += [str(fig_titles[exp]) + (" (" + str(len(n_ensemble[exp])) + ")"
                                                           if exp in list(n_ensemble.keys()) else "")]
                            lv += ["center"]
                            lx += [x1 + 2 * dx]
                            ly += [y2 - (kk * 8 + 15) * dy]
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
    a1, a2 = str(smile_minimum_size).zfill(2), "".join([str(v) for v in proc_th.values()])
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(["", proc_md, str(a1) + "mem", str(a2) + "po"])
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
e14_spread_decomposition_plot()
