# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Look for best correlation between ENSO sstA change (departure from the end of the 20th century) and zonal gradient
#    changes. Here the zonal gradient change is computed on deltaT and W-E regions are defined in each model
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from datetime import datetime
from os.path import dirname as os__path__dirname
import sys
from typing import Literal

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, scipy__stats__pearsonr, correlation_pearson
from enso_change.wrapper.xarray_base import array_wrapper, cf_dim_to_dim, concatenate, mean_weighted
from enso_change.wrapper.xarray_tools import create_land_sea_mask, select_region
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
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_ts_val_glob", "ave_ts_val_tron", "ave_ts_val"],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    "data_experiments": ["historical", "scenario26", "scenario45", "scenario85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    # "smile_rejected": ["CMCC-CM", "EC-EARTH", "KACE-1-0-G"],
    # "smile_rejected": ["EC-EARTH", "KACE-1-0-G"],
    "smile_rejected": ["EC-EARTH", "MPI-ESM-LR"],
    # "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 60.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co":{},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": {"ave": 4, "std": 2, "ske": 2},
    # threshold for selection
    "proc_se": 100 / 3,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {},
    # colors
    "fig_colors": {},
    # line styles
    "fig_linestyles": {},
    # markers
    "fig_markers": {},
    # ticks
    "fig_ticks": {},
    # titles
    "fig_titles": {},
    # panel parameters (to modify default values in enso_evolution_with_cc.lib_plot.panels)
    "panel_param": {},
}
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Main
# ---------------------------------------------------------------------------------------------------------------------#
def e71_change_best_correlation_zonal_grad_plot(
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
        proc_se: float = None,
        proc_th: dict[str: int] = None,
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
    proc_se = none_to_default(proc_se, default["proc_se"])
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
    # remove datasets
    dia1, dur1, pro1 = data_diagnostics2[0], data_epoch_lengths[0], data_projects[0]
    dat1, exp1 = "multimodel001", "historical"
    da = values_netcdf[dia1][dur1][pro1][exp1][dat1]
    dict_datasets = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0]) for k in da.attrs["dataset"].split("; ")}
    to_remove = [dict_datasets[k] for k in smile_rejected if k in list(dict_datasets.keys())]
    if len(to_remove) > 0:
        for dia in list(values_netcdf.keys()):
            for dur in list(values_netcdf[dia].keys()):
                for pro in list(values_netcdf[dia][dur].keys()):
                    if pro != "cmip":
                        continue
                    for exp in list(values_netcdf[dia][dur][pro].keys()):
                        for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                            # array
                            da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                            del values_netcdf[dia][dur][pro][exp][dat]
                            # list models by number
                            list_mod = da["dataset"].values.tolist()
                            # remove model numbers
                            list_mod = sorted(list(set(list_mod) - set(to_remove)))
                            # select models
                            da = da.sel(indexers={"dataset": list_mod})
                            # save
                            put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # select cmip6 only
    # for dia in list(values_netcdf.keys()):
    #     for dur in list(values_netcdf[dia].keys()):
    #         for pro in list(values_netcdf[dia][dur].keys()):
    #             if pro != "cmip":
    #                 continue
    #             for exp in list(values_netcdf[dia][dur][pro].keys()):
    #                 for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
    #                     # array
    #                     da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
    #                     del values_netcdf[dia][dur][pro][exp][dat]
    #                     # list models by number
    #                     list_mod = da["dataset"].values.tolist()
    #                     # remove model numbers
    #                     list_mod = sorted([k for k in list_mod if k > 100])
    #                     # select models
    #                     da = da.sel(indexers={"dataset": list_mod})
    #                     # save
    #                     put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # number of datasets
    n_ensemble = {}
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 4:
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                if pro != "cmip":
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
    # -- compute rSST
    dia2, dia3 = "ave_ts_val_tron", data_diagnostics2[-1]
    if dia2 in list(values_netcdf.keys()) and dia3 in list(values_netcdf.keys()):
        dia_o = dia3.replace("_sst_", "_rsst_").replace("_ts_", "_rts_")
        for dur in list(values_netcdf[dia3].keys()):
            for pro in list(values_netcdf[dia3][dur].keys()):
                for exp in list(values_netcdf[dia3][dur][pro].keys()):
                    for dat in list(values_netcdf[dia3][dur][pro][exp].keys()):
                        if not check_dict_levels(values_netcdf, dia2, dur, pro, exp, dat):
                            continue
                        # arrays
                        da1 = copy__deepcopy(values_netcdf[dia3][dur][pro][exp][dat])
                        da2 = copy__deepcopy(values_netcdf[dia2][dur][pro][exp][dat])
                        # select common years
                        list_y = sorted(list(set(da1["year"].values.tolist()) & set(da2["year"].values.tolist())))
                        da1, da2 = da1.sel(indexers={"year": list_y}), da2.sel(indexers={"year": list_y})
                        # remove tropical mean SST
                        da1 -= da2
                        # save
                        put_in_dict(values_netcdf, da1, dia_o, dur, pro, exp, dat)
        metadata_netcdf[dia_o] = copy__deepcopy(metadata_netcdf[dia3])
        val = metadata_netcdf[dia3]["name_short"].replace("SST", "rSST").replace("ts", "rTS")
        metadata_netcdf[dia_o]["name_short"] = val
    for dia1 in [dia2]:  # [dia2, dia3]:
        if dia1 in list(values_netcdf.keys()):
            del values_netcdf[dia1]
    print("compute rSST", str(len(list(values_netcdf.keys()))).rjust(3), list(values_netcdf.keys()))
    # -- create land-sea mask
    dia2 = [k for k in data_diagnostics2 if len(k.split("_")) == 3][0]
    dur2 = data_epoch_lengths[0]
    pro2 = data_projects[0]
    exp2 = data_experiments[0]
    dat2 = list(values_netcdf[dia2][dur2][pro2][exp2].keys())[0]
    mask = create_land_sea_mask(values_netcdf[dia2][dur2][pro2][exp2][dat2])
    for dia in list(values_netcdf.keys()):
        if len(dia.split("_")) != 3:
            continue
        for dur in list(values_netcdf[dia].keys()):
            for pro in list(values_netcdf[dia][dur].keys()):
                for exp in list(values_netcdf[dia][dur][pro].keys()):
                    for dat in list(values_netcdf[dia][dur][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia][dur][pro][exp][dat])
                        del values_netcdf[dia][dur][pro][exp][dat]
                        # mask land
                        da = da.where(mask == 0)
                        # save
                        put_in_dict(values_netcdf, da, dia, dur, pro, exp, dat)
    # -- compute change compared to initial conditions
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
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
                        attributes = copy__deepcopy(da1.attrs)
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
                            da_o.attrs.update(**attributes)
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o), dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- compute sst change as a change of gradient
    # dia_y, war2 = "std_ts_ano_n34e", "pcc"
    # print(datetime.now().strftime("%m/%d %H-%M"))
    # for dur1 in list(mean_change.keys()):
    #     if dur1 != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_change[dur1].keys()):
    #         if dia_x != "ave_ts_val":
    #             continue
    #         for pro in list(mean_change[dur1][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             # -- get data per dataset per experiment
    #             list_exp = [k for k in list(mean_change[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat1 = sorted(list(mean_change[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
    #             t1 = datetime.now()
    #             print(dia_x, dia_y, dur1, pro, t1.strftime("%m/%d %H-%M"), "get data", list_exp)
    #             dict_dat_x, dict_dat_y = {}, {}
    #             for exp1 in list_exp:
    #                 if not check_dict_levels(mean_change, dur1, dia_x, pro, exp1, dat1, war2) or \
    #                         not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat1, war2):
    #                     continue
    #                 # array
    #                 da_x = copy__deepcopy(mean_change[dur1][dia_x][pro][exp1][dat1][war2])
    #                 da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat1][war2])
    #                 # list datasets
    #                 dx = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0])
    #                       for k in da_x.attrs["dataset"].split("; ")}
    #                 dx = {k: v for k, v in dx.items() if v in da_x["dataset"].values}
    #                 dy = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0])
    #                       for k in da_y.attrs["dataset"].split("; ")}
    #                 dy = {k: v for k, v in dy.items() if v in da_y["dataset"].values}
    #                 list_dat = sorted(list(set(list(dx.keys())) & set(list(dy.keys()))), key=lambda s: s.lower())
    #                 for dat2 in list_dat:
    #                     da_xt = da_x.sel(indexers={"dataset": dict_datasets[dat2]})
    #                     da_yt = da_y.sel(indexers={"dataset": dict_datasets[dat2]})
    #                     put_in_dict(dict_dat_x, da_xt, dat2, exp1)
    #                     put_in_dict(dict_dat_y, da_yt, exp1, dat2)
    #             t2 = datetime.now()
    #             delta_t = t2 - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dia_x, dia_y, dur1, pro, t2.strftime("%m/%d %H-%M"), "get data", ":".join([hh, mm, ss]))
    #             # -- select zonal gradient change in each dataset
    #             list_datasets = sorted(list(dict_dat_x.keys()), key=lambda s: s.lower())
    #             print(dia_x, dia_y, dur1, pro, datetime.now().strftime("%m/%d %H-%M"), "zonal gr", len(list_datasets))
    #             grad = {}
    #             for dat2 in list_datasets:
    #                 t2 = datetime.now()
    #                 # concatenate experiments of given dataset
    #                 da = tuple(da for da in dict_dat_x[dat2].values())
    #                 da = concatenate(da, "experiment").mean(dim="experiment")
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # maximize the gradient change
    #                 tmp = {}
    #                 la1, la2, lo3 = -5, 5, 280
    #                 for lo1 in range(80, 141, 5):
    #                     for lo2 in range(lo1 + 70, 211, 5):
    #                         # mask outside region
    #                         da_w = select_region(da, [la1, la2], [lo1, lo2])
    #                         da_e = select_region(da, [la1, la2], [lo2, lo3])
    #                         # horizontal average
    #                         da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
    #                         da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
    #                         # save correlation
    #                         reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2, lo3]])
    #                         put_in_dict(tmp, abs(float(da_w - da_e)), reg)
    #                 # rank gradient change
    #                 tmp = dict(sorted(tmp.items(), key=lambda item: item[1]))
    #                 # select the region maximizing the change
    #                 re = list(tmp.keys())[-1]
    #                 la1, la2, lo1, lo2, lo3 = int(re[:4]), int(re[4:8]), int(re[8:12]), int(re[12:16]), int(re[16:])
    #                 for exp1, da in dict_dat_x[dat2].items():
    #                     # mask outside region
    #                     da_w = select_region(da, [la1, la2], [lo1, lo2])
    #                     da_e = select_region(da, [la1, la2], [lo2, lo3])
    #                     # horizontal average
    #                     da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
    #                     da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
    #                     # gradient
    #                     put_in_dict(grad, float(da_w - da_e), exp1, dat2)
    #                 delta_t = datetime.now() - t2
    #                 hh = str().ljust(5) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #                 mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #                 ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #                 print(dia_x, dia_y, dur1, pro, dat2.rjust(20), str(hh) + ":" + str(mm) + ":" + str(ss))
    #             delta_t = datetime.now() - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dia_x, dia_y, dur1, pro, "maximum gradient change", str(hh) + ":" + str(mm) + ":" + str(ss))
    #             t2 = datetime.now()
    #             delta_t = t2 - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dia_x, dia_y, dur1, pro, t2.strftime("%m/%d %H-%M"), "zonal gr", ":".join([hh, mm, ss]))
    #             # -- correlation between maximum gradient change and ENSO change
    #             list_exp = sorted(list(grad.keys()), key=lambda s: s.lower())
    #             print(dia_x, dia_y, dur1, pro, datetime.now().strftime("%m/%d %H-%M"), "correlat", list_exp)
    #             dt = {}
    #             for exp1 in list_exp:
    #                 # array
    #                 arr_x = numpy__array([v for v in grad[exp1].values()])
    #                 arr_y = numpy__array([v for v in dict_dat_y[exp1].values()])
    #                 # correlation
    #                 r_value, p_value = scipy__stats__pearsonr(arr_x, arr_y)
    #                 put_in_dict(dt, float(r_value), exp1)
    #             t2 = datetime.now()
    #             delta_t = t2 - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dia_x, dia_y, dur1, pro, t2.strftime("%m/%d %H-%M"), "correlat", ":".join([hh, mm, ss]))
    #             for k, v in dt.items():
    #                 print(str(k).ljust(20), "{0:.3f}".format(round(v, 3)).rjust(6))
    # -- compute sst change of maximum gradient, tailored to each dataset
    dia_y, war2 = "std_ts_ano_n34e", "pcc"
    print("compute maximum gradient", datetime.now().strftime("%m/%d %H-%M"))
    grad_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 != "ave_ts_val":
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 in ["piControl", "historical"]:
                        continue
                    for dat1 in list(values_netcdf[dia1][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_netcdf, dia1, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia1, dur2, pro, exp2) or \
                                not check_dict_levels(values_netcdf, dia2, dur1, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp1) or \
                                not check_dict_levels(values_netcdf, dia2, dur2, pro, exp2):
                            continue
                        # array
                        t1 = datetime.now()
                        print(dia1, dur1, pro, exp1, dat1, t1.strftime("%m/%d %H-%M"), "get data")
                        da1 = copy__deepcopy(values_netcdf[dia1][dur1][pro][exp1][dat1])
                        attributes = copy__deepcopy(da1.attrs)
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp1][dat1])
                        else:
                            da2 = copy__deepcopy(values_netcdf[dia1][dur2][pro][exp2][dat1])
                        da3 = copy__deepcopy(values_netcdf[dia2][dur1][pro][exp1][dat1])
                        if exp1 in ["scenario26", "scenario45", "scenario85"]:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp1][dat1])
                        else:
                            da4 = copy__deepcopy(values_netcdf[dia2][dur2][pro][exp2][dat1])
                        # dimensions
                        dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                        # select experiment epoch
                        yy = 1850 + int(dur1.split("_")[0]) if exp1 == "historical" else 2100
                        da1, da3 = da1.sel(indexers={"year": yy}), da3.sel(indexers={"year": yy})
                        # select reference epoch
                        yy = 2006
                        da2, da4 = da2.sel(indexers={"year": yy}), da4.sel(indexers={"year": yy})
                        # list datasets
                        list_dat = []
                        for dd in [da1, da2, da3, da4]:
                            dict_d = {k.split(" ")[1].split("=")[0]: int(k.split(" ")[0])
                                      for k in dd.attrs["dataset"].split("; ")}
                            dict_d = {k: v for k, v in dict_d.items() if v in dd["dataset"].values}
                            if len(list_dat) == 0:
                                list_dat = list(dict_d.keys())
                            else:
                                list_dat = list(set(list_dat) & set(list(dict_d.keys())))
                        list_dat = [k for k in list(dict_datasets.keys()) if k in list_dat]
                        print(dia1, dur1, pro, exp1, dat1, datetime.now().strftime("%m/%d %H-%M"), "dat loop",
                              len(list_dat))
                        grad = []
                        for dat2 in list_dat:
                            t2 = datetime.now()
                            list_array = [da1, da2] if exp1[:3] in ["rcp", "sce", "ssp"] else [da2, da1]
                            change_max_grad = None
                            for da in list_array:
                                # select dataset
                                da_t = da.sel(indexers={"dataset": dict_datasets[dat2]})
                                # maximize the zonal gradient
                                tmp = {}
                                la1, la2, lo3 = -5, 5, 280
                                for lo1 in range(80, 141, 5):
                                    for lo2 in range(lo1 + 70, 211, 5):
                                        # mask outside region
                                        da_w = select_region(da_t, [la1, la2], [lo1, lo2])
                                        da_e = select_region(da_t, [la1, la2], [lo2, lo3])
                                        # horizontal average
                                        da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
                                        da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
                                        # save correlation
                                        reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2, lo3]])
                                        put_in_dict(tmp, float(da_w - da_e), reg)
                                # rank gradient change
                                tmp = dict(sorted(tmp.items(), key=lambda item: item[1]))
                                if change_max_grad is None:
                                    change_max_grad = copy__deepcopy(tmp[list(tmp.keys())[-1]])
                                else:
                                    change_max_grad -= copy__deepcopy(tmp[list(tmp.keys())[-1]])
                            t3 = datetime.now()
                            delta_t = t3 - t2
                            hh = str().ljust(2) + str(int(round(delta_t.seconds // 3600))).zfill(2)
                            mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
                            ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
                            print(dia1, dur1, pro, exp1, dat1, t3.strftime("%m/%d %H-%M"), dat2.rjust(20),
                                  ":".join([hh, mm, ss]))
                            grad.append(change_max_grad)
                        grad = array_wrapper(numpy__array(grad), attrs=attributes, name="zonal_gradient",
                                             coords={"dataset": numpy__array([dict_datasets[k] for k in list_dat])})
                        # degree of GMST
                        da3 = (da3 - da4 if exp1[:3] in ["rcp", "sce", "ssp"] else da4 - da3).compute()
                        # change per CC or not
                        l1 = ["raw", "pcc"] if exp1[:3] in ["rcp", "sce", "ssp"] else ["raw"]
                        t2 = datetime.now()
                        print(dia1, dur1, pro, exp1, dat1, t2.strftime("%m/%d %H-%M"), "war loop", l1)
                        for war in l1:
                            da_o = copy__deepcopy(grad)
                            # per degree of GMST
                            if war == "pcc":
                                da_o /= da3
                            da_o.attrs.update(**attributes)
                            # save
                            put_in_dict(grad_change, copy__deepcopy(da_o),
                                        dur1, dia1, pro, exp1, dat1, war)
    print("compute correlations", datetime.now().strftime("%m/%d %H-%M"))
    for dur1 in list(grad_change.keys()):
        if dur1 != "050_year_epoch":
            continue
        for dia_x in list(grad_change[dur1].keys()):
            if dia_x != "ave_ts_val":
                continue
            for pro in list(grad_change[dur1][dia_x].keys()):
                if pro != "cmip":
                    continue
                # -- get data per dataset per experiment
                list_exp = [k for k in list(grad_change[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
                dat1 = sorted(list(grad_change[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
                t1 = datetime.now()
                print(dia_x, dia_y, dur1, pro, t1.strftime("%m/%d %H-%M"), "correlat", list_exp)
                dt = {}
                for exp1 in list_exp:
                    if not check_dict_levels(grad_change, dur1, dia_x, pro, exp1, dat1, war2) or \
                            not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat1, war2):
                        continue
                    # array
                    da_x = copy__deepcopy(grad_change[dur1][dia_x][pro][exp1][dat1][war2])
                    da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat1][war2])
                    # correlation
                    dt[exp1] = float(correlation_pearson(da_x, da_y, dim="dataset").values)
                t2 = datetime.now()
                delta_t = t2 - t1
                hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
                mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
                ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
                print(dia_x, dia_y, dur1, pro, t2.strftime("%m/%d %H-%M"), "correlat", ":".join([hh, mm, ss]))
                for k, v in dt.items():
                    print(str(k).rjust(30), "{0:.3f}".format(round(v, 3)).rjust(8))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e71_change_best_correlation_zonal_grad_plot()
