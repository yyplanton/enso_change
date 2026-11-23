# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: evolution of variable change in different experiments for models with low / large ENSO change
# Here scatterplots of mean state change vs. ENSO change are used
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
from enso_change.wrapper.statistics import linear_regression, student_confidence_interval
from enso_change.wrapper.xarray_base import numpy__array, concatenate
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Default arguments
# ---------------------------------------------------------------------------------------------------------------------#
default = {
    #
    # -- Data
    #
    # list of diagnostics: list[str]
    "data_diagnostics1": ["std_ts_ano_n34e",
                          "ave_pr_val_glob",
                          "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
                          "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
                          "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
                          "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw", "ave_pr_val_peso",
                          "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw",
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
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_pr_val_glob",
    #                       "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
    #                       "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
    #                       "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
    #                       "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw", "ave_pr_val_peso",
    #                       "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
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
    # "data_diagnostics1": ["std_ts_ano_n34e", "ave_pr_val_n34s", "ave_pr_val_n40s", "ave_ts_val_aequ",
    #                       "ave_ts_val_atno", "ave_ts_val_isee"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
    # "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["piControl", "historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp585"],
    "data_experiments": ["historical", "ssp126", "ssp245", "ssp370", "ssp585"],
    # "data_experiments": ["piControl", "historical", "rcp26", "rcp45", "rcp85"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets": default_arguments["data_observational_datasets"],
    # list of projects: list[str]
    "data_projects": ["cmip6"],  # default_arguments["data_projects"],
    # "data_projects": ["cmip5", "observations"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 3,
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
    "proc_th": 40,
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
        "panel_01": {"x_delt": 400, "x_size": 800, "y_delt": 100, "y_size": 800},
        "panel_11": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_21": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_31": {"x_delt": 400, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "all": "k",
        "bot": "b",
        "top": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "all": "-",
        "bot": "-",
        "top": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_pr_val_n34s": {
                "panel_11": {
                    "lim": [-1.2, 1.2],
                    "nam": "departure from YEAR1-YEAR2\nVALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
                "panel_21": {
                    "lim": [-1.2, 1.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
                "panel_31": {
                    "lim": [-1.2, 1.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
            },
            "ave_pr_val_n40s": {
                "panel_11": {
                    "lim": [-1.2, 1.2],
                    "nam": "departure from YEAR1-YEAR2\nVALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
                "panel_21": {
                    "lim": [-1.2, 1.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
                "panel_31": {
                    "lim": [-1.2, 1.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 13, 6)],
                },
            },
            "ave_ts_val_aequ": {
                "panel_11": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2\nVALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_21": {
                    "lim": [-1.4, 4.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_31": {
                    "lim": [-1.4, 4.2],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_atno": {
                "panel_11": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2\nVALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_21": {
                    "lim": [-1.2, 3.6],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_31": {
                    "lim": [-1.2, 3.6],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
            },
            "ave_ts_val_isee": {
                "panel_11": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2\nVALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_21": {
                    "lim": [-1.2, 3.6],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_31": {
                    "lim": [-1.2, 3.6],
                    "nam": "",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_11": {
                    "lim": [-25, 75],
                    "nam": "relative departure\nfrom YEAR1-YEAR2\nVALUE (%)",
                    "tic": list(range(-25, 76, 25)),
                },
                "panel_21": {
                    "lim": [-25, 75],
                    "nam": "",
                    "tic": list(range(-25, 76, 25)),
                },
                "panel_31": {
                    "lim": [-25, 75],
                    "nam": "",
                    "tic": list(range(-25, 76, 25)),
                },
            },
            "time": {
                "panel_11": {
                    "lim": [1860, 2020],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(1850, 2021, 50)),
                },
                "panel_21": {
                    "lim": [1990, 2090],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(2000, 2101, 50)),
                },
                "panel_31": {
                    "lim": [1990, 2090],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(2000, 2101, 50)),
                },
            },
        },
        "050_year_epoch": {
            "ave_pr_val_glob": {
                "panel_01": {
                    "lim": [2.5, 3.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(25, 34, 2)],
                },
            },
            "ave_pr_val_n30e": {
                "panel_01": {
                    "lim": [0.5, 4.9],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(5, 50, 11)],
                },
            },
            "ave_pr_val_n30n": {
                "panel_01": {
                    "lim": [3.5, 8.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(35, 84, 12)],
                },
            },
            "ave_pr_val_n30s": {
                "panel_01": {
                    "lim": [0.8, 5.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(8, 57, 12)],
                },
            },
            "ave_pr_val_n34e": {
                "panel_01": {
                    "lim": [0.2, 4.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(2, 47, 11)],
                },
            },
            "ave_pr_val_n34n": {
                "panel_01": {
                    "lim": [3.3, 8.9],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(33, 90, 14)],
                },
            },
            "ave_pr_val_n34s": {
                "panel_01": {
                    "lim": [2.5, 7.7],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(25, 78, 13)],
                },
            },
            "ave_pr_val_n40e": {
                "panel_01": {
                    "lim": [0.9, 7.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(9, 74, 16)],
                },
            },
            "ave_pr_val_n40n": {
                "panel_01": {
                    "lim": [4.6, 9.0],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(46, 91, 11)],
                },
            },
            "ave_pr_val_n40s": {
                "panel_01": {
                    "lim": [5.6, 10.4],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(56, 105, 12)],
                },
            },
            "ave_pr_val_peno": {
                "panel_01": {
                    "lim": [4.3, 8.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(43, 84, 10)],
                },
            },
            "ave_pr_val_peqe": {
                "panel_01": {
                    "lim": [0.8, 5.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(8, 53, 11)],
                },
            },
            "ave_pr_val_pequ": {
                "panel_01": {
                    "lim": [0.8, 5.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(8, 53, 11)],
                },
            },
            "ave_pr_val_peqw": {
                "panel_01": {
                    "lim": [2.0, 7.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(20, 77, 14)],
                },
            },
            "ave_pr_val_peso": {
                "panel_01": {
                    "lim": [2.9, 6.5],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(29, 66, 9)],
                },
            },
            "ave_pr_val_ptre": {
                "panel_01": {
                    "lim": [1.4, 5.8],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(14, 59, 11)],
                },
            },
            "ave_pr_val_ptro": {
                "panel_01": {
                    "lim": [3.4, 4.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(34, 47, 3)],
                },
            },
            "ave_pr_val_ptrw": {
                "panel_01": {
                    "lim": [3.8, 8.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(38, 87, 12)],
                },
            },
            "ave_ts_val_aequ": {
                "panel_01": {
                    "lim": [25.8, 29.0],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(258, 291, 8)],
                },
            },
            "ave_ts_val_atno": {
                "panel_01": {
                    "lim": [24.0, 27.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(240, 272, 8)],
                },
            },
            "ave_ts_val_atro": {
                "panel_01": {
                    "lim": [23.4, 26.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(234, 263, 7)],
                },
            },
            "ave_ts_val_glob": {
                "panel_01": {
                    "lim": [12.5, 16.5],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(125, 166, 10)],
                },
            },
            "ave_ts_val_gxip": {
                "panel_01": {
                    "lim": [1.3, 3.7],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(13, 38, 6)],
                },
            },
            "ave_ts_val_gxpb": {
                "panel_01": {
                    "lim": [0.6, 3.8],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(6, 39, 9)],
                },
            },
            "ave_ts_val_gxpe": {
                "panel_01": {
                    "lim": [0.6, 3.8],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(6, 39, 8)],
                },
            },
            "ave_ts_val_gxpl": {
                "panel_01": {
                    "lim": [1.0, 4.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(10, 43, 8)],
                },
            },
            "ave_ts_val_gxpn": {
                "panel_01": {
                    "lim": [0.9, 4.5],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(9, 46, 9)],
                },
            },
            "ave_ts_val_gxpr": {
                "panel_01": {
                    "lim": [-1.7, 0.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-17, 4, 5)],
                },
            },
            "ave_ts_val_gxps": {
                "panel_01": {
                    "lim": [1.6, 4.8],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(16, 49, 8)],
                },
            },
            "ave_ts_val_gxpw": {
                "panel_01": {
                    "lim": [1.4, 4.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(14, 43, 7)],
                },
            },
            "ave_ts_val_gypb": {
                "panel_01": {
                    "lim": [0.1, 2.1],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(1, 22, 5)],
                },
            },
            "ave_ts_val_gype": {
                "panel_01": {
                    "lim": [1.8, 4.2],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(18, 43, 6)],
                },
            },
            "ave_ts_val_gypi": {
                "panel_01": {
                    "lim": [-0.6, 2.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(-6, 27, 8)],
                },
            },
            "ave_ts_val_iequ": {
                "panel_01": {
                    "lim": [27.3, 30.1],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(273, 302, 7)],
                },
            },
            "ave_ts_val_isee": {
                "panel_01": {
                    "lim": [26.7, 30.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(267, 304, 9)],
                },
            },
            "ave_ts_val_itro": {
                "panel_01": {
                    "lim": [24.5, 28.1],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(245, 282, 9)],
                },
            },
            "ave_ts_val_iweq": {
                "panel_01": {
                    "lim": [26.6, 29.4],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(266, 295, 7)],
                },
            },
            "ave_ts_val_n30e": {
                "panel_01": {
                    "lim": [23.8, 27.4],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(238, 275, 9)],
                },
            },
            "ave_ts_val_n30n": {
                "panel_01": {
                    "lim": [26.0, 28.8],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(260, 289, 7)],
                },
            },
            "ave_ts_val_n30s": {
                "panel_01": {
                    "lim": [24.5, 28.5],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(245, 286, 10)],
                },
            },
            "ave_ts_val_n34e": {
                "panel_01": {
                    "lim": [24.0, 28.4],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(240, 285, 11)],
                },
            },
            "ave_ts_val_n34n": {
                "panel_01": {
                    "lim": [25.8, 28.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(258, 287, 7)],
                },
            },
            "ave_ts_val_n34s": {
                "panel_01": {
                    "lim": [26.2, 29.4],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(262, 295, 8)],
                },
            },
            "ave_ts_val_n40e": {
                "panel_01": {
                    "lim": [25.3, 29.7],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(253, 298, 11)],
                },
            },
            "ave_ts_val_n40n": {
                "panel_01": {
                    "lim": [25.9, 29.1],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(259, 292, 8)],
                },
            },
            "ave_ts_val_n40s": {
                "panel_01": {
                    "lim": [27.2, 30.0],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(272, 301, 7)],
                },
            },
            "ave_ts_val_peno": {
                "panel_01": {
                    "lim": [26.2, 29.0],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(262, 291, 7)],
                },
            },
            "ave_ts_val_peqe": {
                "panel_01": {
                    "lim": [24.0, 27.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(240, 277, 9)],
                },
            },
            "ave_ts_val_pequ": {
                "panel_01": {
                    "lim": [24.7, 28.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(247, 284, 9)],
                },
            },
            "ave_ts_val_peqw": {
                "panel_01": {
                    "lim": [26.3, 29.9],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(263, 300, 9)],
                },
            },
            "ave_ts_val_peso": {
                "panel_01": {
                    "lim": [25.4, 28.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(254, 287, 8)],
                },
            },
            "ave_ts_val_ptre": {
                "panel_01": {
                    "lim": [24.4, 28.0],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(244, 281, 9)],
                },
            },
            "ave_ts_val_ptro": {
                "panel_01": {
                    "lim": [24.5, 27.5],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(245, 276, 7)],
                },
            },
            "ave_ts_val_ptrw": {
                "panel_01": {
                    "lim": [26.5, 30.1],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(265, 302, 9)],
                },
            },
            "ave_ts_val_trop": {
                "panel_01": {
                    "lim": [24.5, 27.3],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(245, 274, 7)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_01": {
                    "lim": [0.4, 1.6],
                    "nam": "VALUE (UNITS)",
                    "tic": [round(k / 10, 1) for k in range(4, 17, 3)],
                },
            },
            "std_ts_ano_n34e--change": {
                "panel_01": {
                    "lim": [-35, 105],
                    "nam": "relative departure\nfrom YEAR1-YEAR2\nVALUE (%)",
                    "tic": list(range(-35, 106, 35)),
                },
            },
            "dataset": {
                "panel_01": {
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
        "all": "all",
        "bot": "bottom THR%",
        "top": "top THR%",
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
def e34_change_explained_bias_scatter_plot(
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
                        put_in_dict(mean_bias, da, dur1, pro, dia, exp, dat)
    # -- compute change compared to initial conditions
    mean_change, mean_change_ts = {}, {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            if dia != data_diagnostics1[0] and dur1 == dur2:
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
                        # ensemble mean (and select epoch?)
                        da1 = da1.mean(dim="member")
                        if dur1 == dur2:
                            yy = da1["year"].values
                            yy = yy[0] if exp1 == "historical" else yy[-1]
                            da1 = da1.sel(indexers={"year": yy})
                        # select reference epoch
                        # yy = 2014 if data_projects[0] == "cmip6" else 2004
                        yy = 2015 if data_projects[0] == "cmip6" else 2005
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # save
                        if dur1 == dur2:
                            put_in_dict(mean_change, float(da1.values), dur1, pro, dia, exp1, dat)
                        else:
                            put_in_dict(mean_change_ts, da1, dur1, pro, dia, exp1, dat)
    # -- correlation
    corr = {}
    for dur in list(mean_bias.keys()):
        for pro in list(mean_bias[dur].keys()):
            if not check_dict_levels(mean_change, dur, pro):
                continue
            for dia_x in list(mean_bias[dur][pro].keys()):
                for dia_y in list(mean_change[dur][pro].keys()):
                    for exp1 in list(mean_bias[dur][pro][dia_x].keys()):
                        for exp2 in list(mean_change[dur][pro][dia_y].keys()):
                            list_dat = list(set(list(mean_bias[dur][pro][dia_x][exp1].keys())) &
                                            set(list(mean_change[dur][pro][dia_y][exp2].keys())))
                            # array
                            da_x = numpy__array([mean_bias[dur][pro][dia_x][exp1][dat] for dat in list_dat])
                            da_y = numpy__array([mean_change[dur][pro][dia_y][exp2][dat] for dat in list_dat])
                            if exp2 == "historical":
                                da_y *= -1
                            # correlation
                            ln = linear_regression(da_x, da_y)
                            # save
                            dia_o = str(dia_y) + "--change"
                            put_in_dict(corr, ln, dur, pro, dia_o, dia_x, str(exp1) + "--bias", exp2)
    # check correlations and p-values
    print("check correlations")
    for dur in list(corr.keys()):
        for pro in list(corr[dur].keys()):
            for dia_y in list(corr[dur][pro].keys()):
                for dia_x in list(corr[dur][pro][dia_y].keys()):
                    cc = []
                    for exp1 in list(corr[dur][pro][dia_y][dia_x].keys()):
                        for exp2 in list(corr[dur][pro][dia_y][dia_x][exp1].keys()):
                            cc.append(corr[dur][pro][dia_y][dia_x][exp1][exp2]["r"])
                    cc = numpy__array(cc)
                    if cc.size > 0:
                        c1, c2 = float(cc.mean()), float(abs(cc).max())
                        if abs(c1) >= 0.3 or c2 >= 0.4:
                            print(dur, pro, dia_y, dia_x, "{0:.2f}".format(round(c1, 2)).rjust(8),
                                  ["{0:.2f}".format(round(k, 2)).rjust(8) for k in cc])
    # -- min-max
    print("bias min-max")
    for dur in list(mean_bias.keys()):
        for pro in list(mean_bias[dur].keys()):
            for dia in list(mean_bias[dur][pro].keys()):
                tmp = []
                for exp in list(mean_bias[dur][pro][dia].keys()):
                    for da in mean_bias[dur][pro][dia][exp].values():
                        tmp += [da]
                if len(tmp) > 0:
                    print(dur, pro.ljust(12), dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    print("change min-max")
    for dur in list(mean_change.keys()):
        for pro in list(mean_change[dur].keys()):
            for dia in list(mean_change[dur][pro].keys()):
                tmp = []
                for exp in list(mean_change[dur][pro][dia].keys()):
                    for da in mean_change[dur][pro][dia][exp].values():
                        tmp += [da]
                if len(tmp) > 0:
                    print(dur, pro.ljust(12), dia, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    # -- select top / bottom models
    change_m, change_r = {}, {}
    print("select models for all variables")
    for dur1 in list(mean_change_ts.keys()):
        if dur1 == dur2 or not check_dict_levels(mean_bias, dur2):
            continue
        for pro in list(mean_change_ts[dur1].keys()):
            if not check_dict_levels(mean_bias, dur2, pro):
                continue
            for dia_x in list(mean_bias[dur2][pro].keys()):
                for exp_x in list(mean_bias[dur2][pro][dia_x].keys()):
                    for dia_y in list(mean_change_ts[dur1][pro].keys()):
                        if dia_y not in data_diagnostics1[:1] + [dia_x]:
                            continue
                        for exp_y in list(mean_change_ts[dur1][pro][dia_y].keys()):
                            # dictionaries
                            dx, dy = mean_bias[dur2][pro][dia_x][exp_x], mean_change_ts[dur1][pro][dia_y][exp_y]
                            # list datasets
                            list_datasets = list(set(list(dx.keys())) & set(list(dy.keys())))
                            # regroup biases
                            dict_t = {dat: dx[dat] for dat in list_datasets}
                            # select top / bottom datasets according to bias
                            dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
                            list_datasets = list(dict_t.keys())
                            nn = int(round(len(list_datasets) * proc_th / 100, 0))
                            list_d_bot, list_d_top = list_datasets[:nn], list_datasets[-nn:]
                            # select change values
                            # for k1, l1 in zip(["all", "top", "bot"], [list_datasets, list_d_top, list_d_bot]):
                            for k1, l1 in zip(["top", "bot"], [list_d_top, list_d_bot]):
                                # concatenate datasets
                                da = concatenate(tuple(dy[dat] for dat in l1), "dataset")
                                # multimodel mean
                                da_m = da.mean(dim="dataset")
                                # confidence interval on mean
                                da_r = student_confidence_interval(da, dim="dataset", confidence_level=95)
                                # save
                                dia_o = str(dia_x) + "--" + str(exp_x)
                                put_in_dict(change_m, da_m, dur1, pro, dia_o, dia_y, exp_y, k1)
                                put_in_dict(change_r, da_r, dur1, pro, dia_o, dia_y, exp_y, k1)
    print("change ts min-max")
    for dur in list(change_r.keys()):
        for pro in list(change_r[dur].keys()):
            list_y = [k2 for k1 in list(change_r[dur][pro].keys()) for k2 in list(change_r[dur][pro][k1].keys())]
            list_y = sorted(list(set(list_y)), key=lambda s: s.lower())
            for dia_y in list_y:
                tmp = []
                list_x = [k1 for k1 in list(change_r[dur][pro].keys()) if dia_y in list(change_r[dur][pro][k1].keys())]
                for dia_x in list_x:
                    for exp_y in list(change_r[dur][pro][dia_x][dia_y].keys()):
                        for da in change_r[dur][pro][dia_x][dia_y][exp_y].values():
                            tmp += da.values.flatten().tolist()
                if len(tmp) > 0:
                    print(dur, pro.ljust(12), dia_y, "{0:.2f}".format(round(min(tmp), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tmp), 2)).rjust(8))
    #
    # -- Organize data for plot
    #
    figure_axes, figure_panel_size, plot_data = {}, {}, {}
    # panel 1: scatterplot
    print("scatterplot")
    group, panel = "group_1", "panel_01"
    for dur in list(mean_bias.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for pro in list(mean_bias[dur].keys()):
            list_var_x = [str(k1) + "--" + str(k2) for k1 in list(mean_bias[dur][pro].keys())
                          for k2 in list(mean_bias[dur][pro][k1].keys())]
            for ii, var_x in enumerate(list_var_x):
                dia_x, exp_x, dia_y = var_x.split("--")[0], var_x.split("--")[1], data_diagnostics1[0]
                if not check_dict_levels(mean_change, dur, pro, dia_y):
                    continue
                pan = copy__deepcopy(panel)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, pro, grp, pan)
                print(dur, pro, var_x, tuple_k)
                # dictionaries
                dict_x, dict_y = mean_bias[dur][pro][dia_x][exp_x], mean_change[dur][pro][dia_y]
                # list experiments
                list_exp_y = [k for k in data_experiments if k in list(dict_y.keys())]
                # panel size
                tk = tuple_k[:2]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    put_in_dict(figure_panel_size, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                    for k1, k2 in fig_panel_size.items():
                        if k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, k2, *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], [dia_x, str(dia_y) + "--change"]
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
                            if k3 == "x_nam" and var_x != list_var_x[-1]:
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
                    val = "bias vs. Hist. (1850-" + str(1850 + ny2 - 1) + ")\n & SSPs (2050-" + str(2050 + ny2 - 1) + \
                          ")"
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                val = str(exp_x) + "\n" + str(metadata_json[dia_x]["name_short"])
                put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                # -- marker
                lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                for exp_y in list_exp_y:
                    list_datasets = list(set(list(dict_x.keys())) & set(list(dict_y[exp_y].keys())))
                    for dat in list_datasets:
                        lc += [fig_colors[exp_y] if exp_y in list(fig_colors.keys()) else "cyan"]
                        lm += [fig_markers[exp_y] if exp_y in list(fig_markers.keys()) else "o"]
                        ls += [40]
                        lw += [0.5]
                        lx += [dict_x[dat]]
                        ly += [dict_y[exp_y][dat]]
                        lz += [numpy__random__randint(2, 9)]
                plot_type = "mar"
                if len(lx) > 0 and len(lx) == len(ly):
                    put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(plot_type) + "_c",))
                    put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(plot_type) + "_lw",))
                    put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(plot_type) + "_m",))
                    put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(plot_type) + "_s",))
                    put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(plot_type) + "_x",))
                    put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(plot_type) + "_y",))
                    put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(plot_type) + "_z",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                for exp_y in list_exp_y:
                    # regression
                    list_datasets = list(set(list(dict_x.keys())) & set(list(dict_y[exp_y].keys())))
                    arr_x = numpy__array([dict_x[dat] for dat in list_datasets])
                    arr_y = numpy__array([dict_y[exp_y][dat] for dat in list_datasets])
                    ln = linear_regression(arr_x, arr_y)
                    # plot
                    lc += [fig_colors[exp_y] if exp_y in list(fig_colors.keys()) else "cyan"]
                    ls += ["-"]
                    lw += [2]
                    lx += [[arr_x.min(), arr_x.max()]]
                    ly += [[ln["s"] * k + ln["i"] for k in [arr_x.min(), arr_x.max()]]]
                    lz += [10]
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
                if isinstance(do_legend, bool) and do_legend:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    for k1, exp_y in enumerate(list_exp_y):
                        # regression
                        list_datasets = list(set(list(dict_x.keys())) & set(list(dict_y[exp_y].keys())))
                        arr_x = numpy__array([dict_x[dat] for dat in list_datasets])
                        arr_y = numpy__array([dict_y[exp_y][dat] for dat in list_datasets])
                        ln = linear_regression(arr_x, arr_y)
                        # plot
                        lc += [fig_colors[exp_y] if exp_y in list(fig_colors.keys()) else "cyan"]
                        lf += [12]
                        lh += ["right"]
                        val = copy__deepcopy(fig_titles[exp_y])
                        val += " (r= + " + "{0:.2f}".format(round(ln["r"], 2)) + " ; p=" + \
                               "{0:.2f}".format(round(ln["p"], 2)) + ")"
                        lt += [val]
                        lv += ["center"]
                        lw += ["regular"]
                        lx += [x2 - 2 * dx]
                        ly += [y2 - (k1 * 8 + 5) * dy]
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
    # panel 11: time series
    print("time series")
    group, panel = "group_1", "panel_11"
    for dur1 in list(change_m.keys()):
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        for pro in list(change_m[dur1].keys()):
            list_var_x = list(change_m[dur1][pro].keys())
            for ii, var_x in enumerate(list_var_x):
                for jj, dia_y in enumerate(list(change_m[dur1][pro][var_x].keys())):
                    list_e = [k for k in data_experiments if k in list(change_m[dur1][pro][var_x][dia_y].keys())]
                    for kk, exp_y in enumerate(list_e):
                        tp = "panel_11" if exp_y == list_e[0] else ("panel_31" if exp_y == list_e[-1] else "panel_21")
                        pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj * 30 + kk).zfill(2)
                        grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                        tuple_k = (data_epoch_lengths[-1], pro, grp, pan)
                        print(dur1, pro, var_x, tuple_k)
                        # panel size
                        tk = tuple_k[:2]
                        if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k1, d1 in fig_panel_size.items():
                                if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                                    for k2, d2 in d1.items():
                                        if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                            continue
                                        tmp = copy__deepcopy(d2)
                                        tt = (dur1, "time", tp, "lim")
                                        if k2 == "x_size" and check_dict_levels(fig_ticks, *tt):
                                            k3 = read_from_dict(fig_ticks, *tt)
                                            tmp *= (k3[-1] - k3[0])
                                        put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                                elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                                    put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                        # panel axes definition
                        l1, l2 = ["x_", "y_"], ["time", dia_y]
                        for vv, dd in zip(l1, l2):
                            if check_dict_levels(fig_ticks, dur1, dd, tp):
                                for k1, k2 in fig_ticks[dur1][dd][tp].items():
                                    k3 = str(vv) + str(k1)
                                    if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                        continue
                                    tmp = copy__deepcopy(k2)
                                    if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                        # yy = 2014 if data_projects[0] == "cmip6" else 2004
                                        yy = 2015 if data_projects[0] == "cmip6" else 2005
                                        tmp = tmp.replace("YEAR1",
                                                          str(yy - ny2)).replace("YEAR2", str(yy - 1))
                                        tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                        tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                        tmp = tmp.replace(" ()", "")
                                    if k3 == "x_nam" and var_x != list_var_x[-1]:
                                        tmp = ""
                                    elif k3 == "y_nam" and kk > 0:
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
                            val = str(fig_titles[exp_y]) + " (" + str(len(n_ensemble[exp_y])) + ")"
                            put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                        # -- curve
                        lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                        for sel, da in change_m[dur1][pro][var_x][dia_y][exp_y].items():
                            lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
                            ls += [fig_linestyles[sel] if sel in list(fig_linestyles.keys()) else ":"]
                            lw += [2]
                            lx += [da["year"].values]
                            ly += [da.values]
                            lz += [4]
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
                            if kk == 0:
                                lc += ["k"]
                                ls += ["-"]
                                lw += [2]
                                lx += [[x1 + 10 * dx, x1 + 18 * dx]]
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
                        for sel, da in change_r[dur1][pro][var_x][dia_y][exp_y].items():
                            la += [0.4]
                            lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
                            lx += [da["year"].values]
                            ly1 += [da.min(dim="interval").values]
                            ly2 += [da.max(dim="interval").values]
                            lz += [2]
                        if isinstance(do_legend, bool) and do_legend:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            if kk == 0:
                                la += [0.4]
                                lc += ["k"]
                                lx += [[x1 + 40 * dx, x1 + 48 * dx]]
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
                        lc, lf, lh, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], []
                        if isinstance(do_legend, bool) and do_legend and kk == 0:
                            df = read_from_dict(figure_panel_size, *tk)
                            x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                            dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                            dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                            # legend
                            list_m = ["MMM", "95% CI"]
                            for k1, v1 in enumerate(list_m):
                                lc += ["k"]
                                lf += [15]
                                lh += ["left"]
                                lt += [v1]
                                lv += ["center"]
                                lw += ["regular"]
                                lx += [x1 + (30 * k1 + 20) * dx]
                                ly += [y2 - 5 * dy]
                            for k1, sel in enumerate(list(change_r[dur1][pro][var_x][dia_y][exp_y].keys())):
                                lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
                                lf += [15]
                                lh += ["left"]
                                lt += [copy__deepcopy(fig_titles[sel]).replace("THR", str(proc_th))]
                                lv += ["center"]
                                lw += ["regular"]
                                lx += [x1 + 2 * dx]
                                ly += [y2 - (k1 * 8 + 15) * dy]
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
        for pro in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][pro].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, pro)))
            # plot
            fig_basic(plot_data[dur][pro], list_names, 1, figure_axes[dur][pro], fig_format,
                      figure_name_and_path, figure_panel_size[dur][pro], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e34_change_explained_bias_scatter_plot()
