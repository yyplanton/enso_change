# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure:
#  - Look for best correlation between ENSO sstA change (departure from the end of the 20th century) and PR, TS, rTS, Tx
#    bias (values at the end of the 20th century)
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from datetime import datetime
from math import ceil as math__ceil
from math import floor as math__floor
from os.path import dirname as os__path__dirname
import sys
from typing import Literal

# add local package directory (relative to current file directory) to sys
sys.path.insert(0, "/".join(os__path__dirname(__file__).split("/")[:-2]))
# local package
from enso_change.script_figures.default_arguments import default_arguments
from enso_change.tools import check_dict_levels, none_to_default, put_in_dict
from enso_change.wrapper.processors import reader_json, reader_netcdf
from enso_change.wrapper.statistics import numpy__array, correlation_pearson
from enso_change.wrapper.xarray_base import cf_dim_to_dim, mean_weighted
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
    # "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_taux_val",
    #                       "ave_ts_val_glob", "ave_ts_val_tron", "ave_ts_val"],
    "data_diagnostics2": ["std_ts_ano_n34e", "ave_pr_val", "ave_tx_val", "ave_ts_val_glob", "ave_ts_val_tron",
                          "ave_ts_val"],
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
def e62_bias_best_correlation_plot(
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
    # -- select bias year
    dia2, dur2, exp2 = "ave_ts_val_glob", data_epoch_lengths[-1], "historical"
    mean_bias = {}
    for dia1 in list(values_netcdf.keys()):
        if len(dia1.split("_")) != 3:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp in list(values_netcdf[dia1][dur1][pro].keys()):
                    for dat in list(values_netcdf[dia1][dur1][pro][exp].keys()):
                        # array
                        da = copy__deepcopy(values_netcdf[dia1][dur1][pro][exp][dat])
                        # select year
                        if pro == "cmip" or (pro == "obse" and dur1 == "050_year_epoch"):
                            yy = int(dur1.split("_")[0]) if exp == "piControl" else 2006
                        else:
                            yy = 2025
                        da = da.sel(indexers={"year": yy})
                        # save
                        put_in_dict(mean_bias, da.compute(), dur1, dia1, pro, exp, dat)
    print("mean_bias", list(mean_bias.keys()))
    # -- compute change compared to initial conditions
    mean_change = {}
    for dia1 in list(values_netcdf.keys()):
        if dia1 == dia2:
            continue
        for dur1 in list(values_netcdf[dia1].keys()):
            for pro in list(values_netcdf[dia1][dur1].keys()):
                for exp1 in list(values_netcdf[dia1][dur1][pro].keys()):
                    if exp1 == "piControl":
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
                        l1 = ["raw", "pcc"] if exp1 in ["scenario26", "scenario45", "scenario85"] else ["raw"]
                        for war in l1:
                            da_o = copy__deepcopy(da1)
                            # per degree of GMST
                            if war == "pcc":
                                da_o /= da3
                            # save
                            put_in_dict(mean_change, copy__deepcopy(da_o).compute(),
                                        dur1, dia1, pro, exp1, dat, war)
    print("mean_change", list(mean_change.keys()))
    # -- compute regional mean
    # dia_y, war2 = "std_ts_ano_n34e", "pcc"
    # print(datetime.now().strftime("%m/%d %H-%M"), "region")
    # # dlat = 6
    # # lat_lon_combinations = []
    # # for a in range(-15, 16 - dlat, 1):
    # #     for b in range(a + dlat, 16, 1):
    # #         dlon = max(math__ceil(500 / (5 * abs(b - a))) * 5, 40)
    # #         for c in range(120, 281 - dlon, 1):
    # #             for d in range(c + dlon, 281, 5):
    # #                 lat_lon_combinations.append((a, b, c, d))
    # dlat = 5
    # lat_lon_combinations = []
    # for a in range(-15, 1 - dlat, 1):
    #     for b in range(dlat, 16, 1):
    #         dlon = max(math__ceil(500 / (5 * abs(b - a))) * 5, 40)
    #         for c in range(120, 281 - dlon, 1):
    #             for d in range(c + dlon, 281, 5):
    #                 lat_lon_combinations.append((a, b, c, d))
    # print("nbr regions", str(len(lat_lon_combinations)).rjust(7))
    # for dur1 in list(mean_bias.keys()):
    #     if dur1 != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_bias[dur1].keys()):
    #         if dia_x != "ave_tx_val":
    #             continue
    #         for pro in list(mean_bias[dur1][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             list_exp = [k for k in list(mean_bias[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat = sorted(list(mean_bias[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
    #             t1 = datetime.now()
    #             print(dur1, dia_x, dia_y, pro, t1.strftime("%m/%d %H-%M"), list_exp)
    #             dt = {}
    #             for exp1 in list_exp:
    #                 if not check_dict_levels(mean_bias, dur1, dia_x, pro, exp1, dat) or \
    #                         not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat, war2):
    #                     continue
    #                 t2 = datetime.now()
    #                 print(dur1, dia_x, dia_y, pro, dat, exp1, t2.strftime("%m/%d %H-%M"))
    #                 # array
    #                 da = copy__deepcopy(mean_bias[dur1][dia_x][pro][exp1][dat]).compute()
    #                 da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat][war2]).compute()
    #                 # get lat-lon dimensions
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # region average
    #                 for la1, la2, lo1, lo2 in lat_lon_combinations:
    #                     # mask outside region
    #                     da_x = select_region(da, [la1, la2], [lo1, lo2])
    #                     if da_x.count().values / da_x["dataset"].size < 500:
    #                         continue
    #                     # horizontal average
    #                     da_x = mean_weighted(da_x, dim=[dim_y, dim_x], skipna=True)
    #                     # correlation with ENSO change
    #                     corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
    #                     # save correlation
    #                     reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2]])
    #                     put_in_dict(dt, [corr], reg)
    #                 delta_t = datetime.now() - t2
    #                 hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #                 mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #                 ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #                 print(dur1, dia_x, dia_y, pro, dat, exp1, str(hh) + ":" + str(mm) + ":" + str(ss))
    #             # rank correlations
    #             dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    #             dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    #             for thr in range(10, 81, 1):
    #                 tmp = {k: v for k, v in dt.items() if v >= round(thr / 100, 2)}
    #                 if len(list(tmp.keys())) > 500:
    #                     dt = copy__deepcopy(tmp)
    #                 else:
    #                     break
    #             delta_t = datetime.now() - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dur1, dia_x, dia_y, pro, str(hh) + ":" + str(mm) + ":" + str(ss))
    #             list_k = list(dt.keys())
    #             if dt[list_k[-1]] < 0.4:
    #                 list_k = list_k[-20:]
    #             for k in list_k:
    #                 la1, la2, lo1, lo2 = int(k[:4]), int(k[4:8]), int(k[8:12]), int(k[12:])
    #                 print(dur1, dia_x, dia_y, pro, str(la1).rjust(10) + "N", str(la2).rjust(3) + "N",
    #                       str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E",
    #                       "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    # -- compute regional mean zonal difference
    # dia_y, war2 = "std_ts_ano_n34e", "pcc"
    # print(datetime.now().strftime("%m/%d %H-%M"), "zonal")
    # lat_lon_combinations = []
    # for a in range(-10, -2, 1):
    #     for b in range(max(math__ceil(abs(a) * 2 / 3), 3), min(math__ceil(abs(a) * 3 / 2), 10) + 1, 1):
    #         dlon = max(math__ceil(500 / (5 * abs(b - a))) * 5, 60)
    #         for c in range(80, 281 - 2 * dlon, 1):
    #             for d in range(c + dlon, 281 - dlon, 5):
    #                 for e in range(d + dlon, 281, 5):
    #                     lat_lon_combinations.append((a, b, c, d, e))
    # for dur1 in list(mean_bias.keys()):
    #     if dur1 != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_bias[dur1].keys()):
    #         if dia_x != "ave_ts_val":
    #             continue
    #         for pro in list(mean_bias[dur1][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             list_exp = [k for k in list(mean_bias[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat = sorted(list(mean_bias[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
    #             t1 = datetime.now()
    #             print(dur1, "zonal", dia_x, dia_y, pro, t1.strftime("%m/%d %H-%M"), list_exp)
    #             dt = {}
    #             for exp1 in list_exp:
    #                 if not check_dict_levels(mean_bias, dur1, dia_x, pro, exp1, dat) or \
    #                         not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat, war2):
    #                     continue
    #                 t2 = datetime.now()
    #                 print(dur1, "zonal", dia_x, dia_y, pro, dat, exp1, t2.strftime("%m/%d %H-%M"))
    #                 # array
    #                 da = copy__deepcopy(mean_bias[dur1][dia_x][pro][exp1][dat]).compute()
    #                 da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat][war2]).compute()
    #                 # get lat-lon dimensions
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # region average
    #                 for la1, la2, lo1, lo2, lo3 in lat_lon_combinations:
    #                     da_w = select_region(da, [la1, la2], [lo1, lo2])
    #                     da_e = select_region(da, [la1, la2], [lo2, lo3])
    #                     if da_w.count().values / da_w["dataset"].size < 500 or \
    #                             da_e.count().values / da_e["dataset"].size < 500:
    #                         continue
    #                     # horizontal average
    #                     da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
    #                     da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
    #                     # change of zonal gradient
    #                     da_x = da_w - da_e
    #                     # correlation with ENSO change
    #                     corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
    #                     # save correlation
    #                     reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2, lo3]])
    #                     put_in_dict(dt, [corr], reg)
    #                 delta_t = datetime.now() - t2
    #                 hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #                 mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #                 ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #                 print(dur1, "zonal", dia_x, dia_y, pro, dat, exp1, str(hh) + ":" + str(mm) + ":" + str(ss))
    #                 # rank correlations
    #             dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    #             dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    #             for thr in range(10, 81, 1):
    #                 tmp = {k: v for k, v in dt.items() if v >= round(thr / 100, 2)}
    #                 if len(list(tmp.keys())) > 500:
    #                     dt = copy__deepcopy(tmp)
    #                 else:
    #                     break
    #             delta_t = datetime.now() - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dur1, "zonal", dia_x, dia_y, pro, str(hh) + ":" + str(mm) + ":" + str(ss))
    #             list_k = list(dt.keys())
    #             if dt[list_k[-1]] < 0.4:
    #                 list_k = list_k[-20:]
    #             for k in list_k:
    #                 la1, la2, lo1, lo2, lo3 = int(k[:4]), int(k[4:8]), int(k[8:12]), int(k[12:16]), int(k[16:])
    #                 print(dur1, "zonal", dia_x, dia_y, pro, str(la1).rjust(10) + "N", str(la2).rjust(3) + "N",
    #                       str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E", str(lo3).rjust(3) + "E",
    #                       "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    # dia_y, war2 = "std_ts_ano_n34e", "pcc"
    # print(datetime.now().strftime("%m/%d %H-%M"), "zonal symmetric")
    # lat_lon_combinations = []
    # for a in range(-10, -2, 1):
    #     b = abs(a)
    #     dlon = max(math__ceil(500 / (5 * abs(b - a))) * 5, 60)
    #     for c in range(80, 281 - 2 * dlon, 1):
    #         for d in range(c + dlon, 281 - dlon, 5):
    #             for e in range(d, 281 - dlon, 1):
    #                 for f in range(e + dlon, 281, 5):
    #                     # if not (c % 2 == 0 or c % 5 == 0) or not (e % 2 == 0 or e % 5 == 0):
    #                     #     continue
    #                     lat_lon_combinations.append((a, b, c, d, e, f))
    # for dur1 in list(mean_change.keys()):
    #     if dur1 != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_change[dur1].keys()):
    #         if dia_x != "ave_ts_val":
    #             continue
    #         for pro in list(mean_change[dur1][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             list_exp = [k for k in list(mean_bias[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat = sorted(list(mean_bias[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
    #             t1 = datetime.now()
    #             print(dur1, "zonal symmetric", dia_x, dia_y, pro, t1.strftime("%m/%d %H-%M"), list_exp)
    #             dt = {}
    #             for exp1 in list_exp:
    #                 if not check_dict_levels(mean_bias, dur1, dia_x, pro, exp1, dat) or \
    #                         not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat, war2):
    #                     continue
    #                 t2 = datetime.now()
    #                 print(dur1, "zonal symmetric", dia_x, dia_y, pro, dat, exp1, t2.strftime("%m/%d %H-%M"))
    #                 # array
    #                 da = copy__deepcopy(mean_bias[dur1][dia_x][pro][exp1][dat]).compute()
    #                 da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat][war2]).compute()
    #                 # get lat-lon dimensions
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # region average
    #                 for la1, la2, lo1, lo2, lo3, lo4 in lat_lon_combinations:
    #                     # mask outside region
    #                     da_w = select_region(da, [la1, la2], [lo1, lo2])
    #                     da_e = select_region(da, [la1, la2], [lo3, lo4])
    #                     if da_w.count().values / da_w["dataset"].size < 500 or \
    #                             da_e.count().values / da_e["dataset"].size < 500:
    #                         continue
    #                     # horizontal average
    #                     da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
    #                     da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
    #                     # change of zonal gradient
    #                     da_x = da_w - da_e
    #                     # correlation with ENSO change
    #                     corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
    #                     # save correlation
    #                     reg = "".join([str(k).zfill(4) for k in [la1, la2, lo1, lo2, lo3, lo4]])
    #                     put_in_dict(dt, [corr], reg)
    #                 delta_t = datetime.now() - t2
    #                 hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #                 mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #                 ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #                 print(dur1, "zonal symmetric", dia_x, dia_y, pro, dat, exp1, ":".join([hh, mm, ss]))
    #             # rank correlations
    #             dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    #             dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    #             for thr in range(10, 81, 1):
    #                 tmp = {k: v for k, v in dt.items() if v >= round(thr / 100, 2)}
    #                 if len(list(tmp.keys())) > 500:
    #                     dt = copy__deepcopy(tmp)
    #                 else:
    #                     break
    #             delta_t = datetime.now() - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dur1, "zonal symmetric", dia_x, dia_y, pro, ":".join([hh, mm, ss]))
    #             list_k = list(dt.keys())
    #             if dt[list_k[-1]] < 0.4:
    #                 list_k = list_k[-20:]
    #             for k in list_k:
    #                 la1, la2 = int(k[:4]), int(k[4:8])
    #                 lo1, lo2, lo3, lo4 = int(k[8:12]), int(k[12:16]), int(k[16:20]), int(k[20:])
    #                 print(dur1, "zonal symmetric", dia_x, dia_y, pro, str(la1).rjust(10) + "N", str(la2).rjust(3) + "N",
    #                       str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E", str(lo3).rjust(3) + "E",
    #                       str(lo4).rjust(3) + "E", "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    dia_y, war2 = "std_ts_ano_n34e", "pcc"
    print(datetime.now().strftime("%m/%d %H-%M"), "zonal boxes")
    dlat = 6
    lat_lon_combinations = []
    for la1_w in range(-10, 1 - math__ceil(dlat / 2), 2):
        la2_w = abs(la1_w)
        dlon_w = max(math__ceil(400 / (5 * abs(la2_w - la1_w))) * 5, 40)
        for lo1_w in range(80, 231 - dlon_w, 5):
            if lo1_w % 2 != 0 and lo1_w % 5 != 0:
                continue
            for lo2_w in range(lo1_w + dlon_w, 231, 5):
                for la1_e in range(-10, 1 - math__ceil(dlat / 2), 2):
                    la2_e = abs(la1_e)
                    dlon_e = max(math__ceil(400 / (5 * abs(la2_e - la1_e))) * 5, 40)
                    for lo1_e in range(max(lo2_w, 210), 281 - dlon_e, 5):
                        if lo1_e % 2 != 0 and lo1_e % 5 != 0:
                            continue
                        for lo2_e in range(lo1_e + dlon_e, 281, 5):
                            lat_lon_combinations.append((la1_w, la2_w, lo1_w, lo2_w, la1_e, la2_e, lo1_e, lo2_e))
    print("nbr regions", str(len(lat_lon_combinations)).rjust(10))
    for dur1 in list(mean_change.keys()):
        if dur1 != "050_year_epoch":
            continue
        for dia_x in list(mean_change[dur1].keys()):
            if dia_x != "ave_ts_val":
                continue
            for pro in list(mean_change[dur1][dia_x].keys()):
                if pro != "cmip":
                    continue
                list_exp = [k for k in list(mean_bias[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
                dat = sorted(list(mean_bias[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
                t1 = datetime.now()
                print(dur1, "zonal boxes".ljust(15), dia_x.ljust(16), dia_y, pro, t1.strftime("%m/%d %H-%M"), list_exp)
                dt = {}
                for exp1 in list_exp:
                    if not check_dict_levels(mean_bias, dur1, dia_x, pro, exp1, dat) or \
                            not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat, war2):
                        continue
                    t2 = datetime.now()
                    print(dur1, "zonal boxes".ljust(15), dia_x.ljust(16), dia_y, pro, exp1, dat,
                          t2.strftime("%m/%d %H-%M"))
                    # array
                    da = copy__deepcopy(mean_bias[dur1][dia_x][pro][exp1][dat]).compute()
                    da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat][war2]).compute()
                    # get lat-lon dimensions
                    dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
                    # region average
                    for la1_w, la2_w, lo1_w, lo2_w, la1_e, la2_e, lo1_e, lo2_e in lat_lon_combinations:
                        # mask outside region
                        da_w = select_region(da, [la1_w, la2_w], [lo1_w, lo2_w])
                        da_e = select_region(da, [la1_e, la2_e], [lo1_e, lo2_e])
                        if da_w.count().values / da_w["dataset"].size < 400 or \
                                da_e.count().values / da_e["dataset"].size < 400:
                            continue
                        # horizontal average
                        da_w = mean_weighted(da_w, dim=[dim_y, dim_x], skipna=True)
                        da_e = mean_weighted(da_e, dim=[dim_y, dim_x], skipna=True)
                        # change of zonal gradient
                        da_x = da_w - da_e
                        # correlation with ENSO change
                        corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
                        # save correlation
                        reg = "".join([str(k).zfill(4)
                                       for k in [la1_w, la2_w, lo1_w, lo2_w, la1_e, la2_e, lo1_e, lo2_e]])
                        put_in_dict(dt, [corr], reg)
                    delta_t = datetime.now() - t2
                    hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
                    mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
                    ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
                    print(dur1, "zonal boxes".ljust(15), dia_x.ljust(16), dia_y, pro, ":".join([hh, mm, ss]))
                # rank correlations
                dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
                dt = dict(sorted(dt.items(), key=lambda item: item[1]))
                for thr in range(10, 81, 1):
                    tmp = {k: v for k, v in dt.items() if v >= round(thr / 100, 2)}
                    if len(list(tmp.keys())) > 500:
                        dt = copy__deepcopy(tmp)
                    else:
                        break
                delta_t = datetime.now() - t1
                hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
                mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
                ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
                print(dur1, "zonal boxes".ljust(15), dia_x.ljust(16), dia_y, pro, ":".join([hh, mm, ss]))
                list_k = list(dt.keys())
                if dt[list_k[-1]] < 0.4:
                    list_k = list_k[-20:]
                for k in list_k:
                    la1_w, la2_w, lo1_w, lo2_w = int(k[: 4]), int(k[4: 8]), int(k[8: 12]), int(k[12: 16])
                    la1_e, la2_e, lo1_e, lo2_e = int(k[16: 20]), int(k[20: 24]), int(k[24: 28]), int(k[28:])
                    print(dia_x, "zonal boxes", dia_y, dur1, pro, "west".ljust(10), str(la1_w).rjust(3) + "N",
                          str(la2_w).rjust(3) + "N", str(lo1_w).rjust(3) + "E", str(lo2_w).rjust(3) + "E", "east",
                          str(la1_e).rjust(3) + "N", str(la2_e).rjust(3) + "N", str(lo1_e).rjust(3) + "E",
                          str(lo2_e).rjust(3) + "E", "{0:.3f}".format(round(dt[k], 3)).rjust(6))
    # -- compute regional mean meridional difference
    # dia_y, war2 = "std_ts_ano_n34e", "pcc"
    # print(datetime.now().strftime("%m/%d %H-%M"), "meridional")
    # dlat, dlon = 5, 50
    # lat_lon_combinations = []
    # for a in range(-15, 16 - 2 * dlat, 1):
    #     for b in range(a + 5, 16 - 2 * dlat, 1):
    #         c_min = max(b + math__ceil(abs(b - a) * 2 / 3), b + dlat)
    #         c_max = min(b + math__floor(abs(b - a) * 3 / 2), 15)
    #         for c in range(c_min, c_max + 1, 1):
    #             for d in range(180, 281 - dlon, 1):
    #                 for e in range(d + dlon, 281, 5):
    #                     lat_lon_combinations.append((a, b, c, d, e))
    # for dur1 in list(mean_bias.keys()):
    #     if dur1 != "050_year_epoch":
    #         continue
    #     for dia_x in list(mean_bias[dur1].keys()):
    #         if dia_x != "ave_ts_val":
    #             continue
    #         for pro in list(mean_bias[dur1][dia_x].keys()):
    #             if pro != "cmip":
    #                 continue
    #             list_exp = [k for k in list(mean_bias[dur1][dia_x][pro].keys()) if k[:3] in ["rcp", "sce", "ssp"]]
    #             dat = sorted(list(mean_bias[dur1][dia_x][pro][list_exp[0]].keys()), key=lambda s: s.lower())[0]
    #             t1 = datetime.now()
    #             print(dur1, "meridional", dia_x, dia_y, pro, t1.strftime("%m/%d %H-%M"), list_exp)
    #             dt = {}
    #             for exp1 in list_exp:
    #                 if not check_dict_levels(mean_bias, dur1, dia_x, pro, exp1, dat) or \
    #                         not check_dict_levels(mean_change, dur1, dia_y, pro, exp1, dat, war2):
    #                     continue
    #                 t2 = datetime.now()
    #                 print(dur1, "meridional", dia_x, dia_y, pro, dat, exp1, t2.strftime("%m/%d %H-%M"))
    #                 # array
    #                 da = copy__deepcopy(mean_bias[dur1][dia_x][pro][exp1][dat]).compute()
    #                 da_y = copy__deepcopy(mean_change[dur1][dia_y][pro][exp1][dat][war2]).compute()
    #                 # get lat-lon dimensions
    #                 dim_x, dim_y = cf_dim_to_dim(da, "X"), cf_dim_to_dim(da, "Y")
    #                 # region average
    #                 for la1, la2, la3, lo1, lo2 in lat_lon_combinations:
    #                     # mask outside region
    #                     da_s = select_region(da, [la1, la2], [lo1, lo2])
    #                     da_n = select_region(da, [la2, la3], [lo1, lo2])
    #                     ns = da_s.count().values / da_s["dataset"].size
    #                     nn = da_n.count().values / da_n["dataset"].size
    #                     if ns < 200 or nn < 200 or ns + nn < 450:
    #                         continue
    #                     # horizontal average
    #                     da_s = mean_weighted(da_s, dim=[dim_y, dim_x], skipna=True)
    #                     da_n = mean_weighted(da_n, dim=[dim_y, dim_x], skipna=True)
    #                     # change of meridional gradient
    #                     da_x = da_n - da_s
    #                     # correlation with ENSO change
    #                     corr = float(correlation_pearson(da_x, da_y, dim="dataset").values)
    #                     # save correlation
    #                     reg = "".join([str(k).zfill(4) for k in [la1, la2, la3, lo1, lo2]])
    #                     put_in_dict(dt, [corr], reg)
    #                 delta_t = datetime.now() - t2
    #                 hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #                 mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #                 ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #                 print(dur1, "meridional", dia_x, dia_y, pro, dat, exp1, str(hh) + ":" + str(mm) + ":" + str(ss))
    #                 # rank correlations
    #             dt = {k: abs(float(numpy__array(v).mean())) for k, v in dt.items()}
    #             dt = dict(sorted(dt.items(), key=lambda item: item[1]))
    #             for thr in range(10, 81, 1):
    #                 tmp = {k: v for k, v in dt.items() if v >= round(thr / 100, 2)}
    #                 if len(list(tmp.keys())) > 500:
    #                     dt = copy__deepcopy(tmp)
    #                 else:
    #                     break
    #             delta_t = datetime.now() - t1
    #             hh = str().ljust(20) + str(int(round(delta_t.seconds // 3600))).zfill(2)
    #             mm = str(int(round((delta_t.seconds % 3600) // 60))).zfill(2)
    #             ss = str(int(round((delta_t.seconds % 3600) % 60))).zfill(2)
    #             print(dur1, "meridional", dia_x, dia_y, pro, str(hh) + ":" + str(mm) + ":" + str(ss))
    #             list_k = list(dt.keys())
    #             if dt[list_k[-1]] < 0.4:
    #                 list_k = list_k[-20:]
    #             for k in list_k:
    #                 la1, la2, la3, lo1, lo2 = int(k[:4]), int(k[4:8]), int(k[8:12]), int(k[12:16]), int(k[16:])
    #                 print(dur1, "meridional", dia_x, dia_y, pro, str(la1).rjust(10) + "N", str(la2).rjust(3) + "N",
    #                       str(la3).rjust(3) + "N", str(lo1).rjust(3) + "E", str(lo2).rjust(3) + "E",
    #                       "{0:.3f}".format(round(dt[k], 3)).rjust(6))
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e62_bias_best_correlation_plot()
