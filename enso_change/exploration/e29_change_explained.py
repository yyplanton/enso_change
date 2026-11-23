# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: evolution of variable change in different experiments for models with low / large ENSO change
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
from enso_change.wrapper.xarray_base import concatenate
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
    "data_diagnostics1": ["std_ts_ano_n34e",
                          "ave_pr_val_glob",
                          "ave_pr_val_n30e", "ave_pr_val_n30n", "ave_pr_val_n30s",
                          "ave_pr_val_n34e", "ave_pr_val_n34n", "ave_pr_val_n34s",
                          "ave_pr_val_n40e", "ave_pr_val_n40n", "ave_pr_val_n40s",
                          "ave_pr_val_peno", "ave_pr_val_peqe", "ave_pr_val_pequ", "ave_pr_val_peqw", "ave_pr_val_peso",
                          "ave_pr_val_ptre", "ave_pr_val_ptro", "ave_pr_val_ptrw"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_ts_val_glob",
    #                       "ave_ts_val_aequ", "ave_ts_val_atno", "ave_ts_val_atro",
    #                       "ave_ts_val_iequ", "ave_ts_val_isee", "ave_ts_val_itro", "ave_ts_val_iweq"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl", "ave_ts_val_gxpn",
    #                       "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
    #                       "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_ts_val_n30e", "ave_ts_val_n30n", "ave_ts_val_n30s",
    #                       "ave_ts_val_n34e", "ave_ts_val_n34n", "ave_ts_val_n34s",
    #                       "ave_ts_val_n40e", "ave_ts_val_n40n", "ave_ts_val_n40s"],
    # "data_diagnostics1": ["std_ts_ano_n34e",
    #                       "ave_ts_val_peno", "ave_ts_val_peqe", "ave_ts_val_pequ", "ave_ts_val_peqw", "ave_ts_val_peso",
    #                       "ave_ts_val_ptre", "ave_ts_val_ptro", "ave_ts_val_ptrw", "ave_ts_val_trop"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["030_year_epoch", "050_year_epoch"],
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
    "proc_th": 30,
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
        "panel_1": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
        "panel_2": {"x_delt": 200, "x_size": 4, "y_delt": 100, "y_size": 800},
    },
    # colors
    "fig_colors": {
        "all": "k",
        "bottom": "b",
        "top": "r",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "all": "-",
        "bottom": "-",
        "top": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": default_arguments["fig_markers_dat"],
    # ticks
    "fig_ticks": {
        "030_year_epoch": {
            "ave_pr_val_glob": {
                "panel_1": {
                    "lim": [-0.12, 0.36],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 100, 2) for k in range(-12, 37, 12)],
                },
                "panel_2": {
                    "lim": [-0.12, 0.36],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 100, 2) for k in range(-12, 37, 12)],
                },
            },
            "ave_pr_val_n30e": {
                "panel_1": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_2": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
            },
            "ave_pr_val_n30n": {
                "panel_1": {
                    "lim": [-2, 2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-2, 3, 1)),
                },
                "panel_2": {
                    "lim": [-2, 2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-2, 3, 1)),
                },
            },
            "ave_pr_val_n30s": {
                "panel_1": {
                    "lim": [-1.8, 1.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-18, 19, 9)],
                },
                "panel_2": {
                    "lim": [-1.8, 1.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-18, 19, 9)],
                },
            },
            "ave_pr_val_n34e": {
                "panel_1": {
                    "lim": [-1.3, 3.9],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-13, 40, 13)],
                },
                "panel_2": {
                    "lim": [-1.3, 3.9],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-13, 40, 13)],
                },
            },
            "ave_pr_val_n34n": {
                "panel_1": {
                    "lim": [-0.8, 2.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-8, 25, 8)],
                },
                "panel_2": {
                    "lim": [-0.8, 2.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-8, 25, 8)],
                },
            },
            "ave_pr_val_n34s": {
                "panel_1": {
                    "lim": [-2.2, 2.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-22, 23, 11)],
                },
                "panel_2": {
                    "lim": [-2.2, 2.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-22, 23, 11)],
                },
            },
            "ave_pr_val_n40e": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_pr_val_n40n": {
                "panel_1": {
                    "lim": [-1.6, 1.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 17, 8)],
                },
                "panel_2": {
                    "lim": [-1.6, 1.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 17, 8)],
                },
            },
            "ave_pr_val_n40s": {
                "panel_1": {
                    "lim": [-1.8, 1.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-18, 19, 9)],
                },
                "panel_2": {
                    "lim": [-1.8, 1.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-18, 19, 9)],
                },
            },
            "ave_pr_val_peno": {
                "panel_1": {
                    "lim": [-1.6, 1.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 17, 8)],
                },
                "panel_2": {
                    "lim": [-1.6, 1.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 17, 8)],
                },
            },
            "ave_pr_val_peqe": {
                "panel_1": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_2": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
            },
            "ave_pr_val_pequ": {
                "panel_1": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
                "panel_2": {
                    "lim": [-1.2, 3.6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 37, 12)],
                },
            },
            "ave_pr_val_peqw": {
                "panel_1": {
                    "lim": [-1, 3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-1, 4, 1)),
                },
                "panel_2": {
                    "lim": [-1, 3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-1, 4, 1)),
                },
            },
            "ave_pr_val_peso": {
                "panel_1": {
                    "lim": [-1.4, 1.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 15, 7)],
                },
                "panel_2": {
                    "lim": [-1.4, 1.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 15, 7)],
                },
            },
            "ave_pr_val_ptre": {
                "panel_1": {
                    "lim": [-1.1, 3.3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-11, 34, 11)],
                },
                "panel_2": {
                    "lim": [-1.1, 3.3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-11, 34, 11)],
                },
            },
            "ave_pr_val_ptro": {
                "panel_1": {
                    "lim": [-0.3, 0.9],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-3, 10, 3)],
                },
                "panel_2": {
                    "lim": [-0.3, 0.9],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-3, 10, 3)],
                },
            },
            "ave_pr_val_ptrw": {
                "panel_1": {
                    "lim": [-1, 3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-1, 4, 1)),
                },
                "panel_2": {
                    "lim": [-1, 3],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-1, 4, 1)),
                },
            },
            "ave_ts_val_aequ": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
                "panel_2": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_atno": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_atro": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_glob": {
                "panel_1": {
                    "lim": [-2, 6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-2, 7, 2)),
                },
                "panel_2": {
                    "lim": [-2, 6],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-2, 7, 2)),
                },
            },
            "ave_ts_val_iequ": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_isee": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_itro": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_iweq": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
                "panel_2": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_gxip": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gxpb": {
                "panel_1": {
                    "lim": [-1.2, 0.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 5, 4)],
                },
                "panel_2": {
                    "lim": [-1.2, 0.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 5, 4)],
                },
            },
            "ave_ts_val_gxpe": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gxpl": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gxpn": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gxpr": {
                "panel_1": {
                    "lim": [-0.4, 1.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-4, 13, 4)],
                },
                "panel_2": {
                    "lim": [-0.4, 1.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-4, 13, 4)],
                },
            },
            "ave_ts_val_gxps": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },"ave_ts_val_gxpw": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gypb": {
                "panel_1": {
                    "lim": [-1.2, 0.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 5, 4)],
                },
                "panel_2": {
                    "lim": [-1.2, 0.4],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-12, 5, 4)],
                },
            },
            "ave_ts_val_gype": {
                "panel_1": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
                "panel_2": {
                    "lim": [-1.5, 0.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 6, 5)],
                },
            },
            "ave_ts_val_gypi": {
                "panel_1": {
                    "lim": [-0.7, 2.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-7, 22, 7)],
                },
                "panel_2": {
                    "lim": [-0.7, 2.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-7, 22, 7)],
                },
            },
            "ave_ts_val_n30e": {
                "panel_1": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
                "panel_2": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
            },
            "ave_ts_val_n30n": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_n30s": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_n34e": {
                "panel_1": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
                "panel_2": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
            },
            "ave_ts_val_n34n": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_n34s": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_n40e": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
                "panel_2": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_n40n": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_n40s": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_peno": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_peqe": {
                "panel_1": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
                "panel_2": {
                    "lim": [-1.7, 5.1],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-17, 52, 17)],
                },
            },
            "ave_ts_val_pequ": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
                "panel_2": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_peqw": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_peso": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_ptre": {
                "panel_1": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
                "panel_2": {
                    "lim": [-1.6, 4.8],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-16, 49, 16)],
                },
            },
            "ave_ts_val_ptro": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "ave_ts_val_ptrw": {
                "panel_1": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
                "panel_2": {
                    "lim": [-1.5, 4.5],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-15, 46, 15)],
                },
            },
            "ave_ts_val_trop": {
                "panel_1": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
                "panel_2": {
                    "lim": [-1.4, 4.2],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": [round(k / 10, 1) for k in range(-14, 43, 14)],
                },
            },
            "std_ts_ano_n34e": {
                "panel_1": {
                    "lim": [-80, 80],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-80, 81, 40)),
                },
                "panel_2": {
                    "lim": [-80, 80],
                    "nam": "departure from YEAR1-YEAR2",
                    "tic": list(range(-80, 81, 40)),
                },
            },
            "time": {
                "panel_1": {
                    "lim": [1860, 2020],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(1850, 2021, 50)),
                },
                "panel_2": {
                    "lim": [1990, 2090],
                    "nam": "Year",
                    "nbr_minor": 4,
                    "tic": list(range(2000, 2101, 50)),
                },
            },
        },
    },
    # titles
    "fig_titles": {
        "ave_pr_val_glob": "change of\nVALUE (UNITS)",
        "ave_pr_val_n30e": "change of\nVALUE (UNITS)",
        "ave_pr_val_n30n": "change of\nVALUE (UNITS)",
        "ave_pr_val_n30s": "change of\nVALUE (UNITS)",
        "ave_pr_val_n34e": "change of\nVALUE (UNITS)",
        "ave_pr_val_n34n": "change of\nVALUE (UNITS)",
        "ave_pr_val_n34s": "change of\nVALUE (UNITS)",
        "ave_pr_val_n40e": "change of\nVALUE (UNITS)",
        "ave_pr_val_n40n": "change of\nVALUE (UNITS)",
        "ave_pr_val_n40s": "change of\nVALUE (UNITS)",
        "ave_pr_val_peno": "change of\nVALUE (UNITS)",
        "ave_pr_val_peqe": "change of\nVALUE (UNITS)",
        "ave_pr_val_pequ": "change of\nVALUE (UNITS)",
        "ave_pr_val_peqw": "change of\nVALUE (UNITS)",
        "ave_pr_val_peso": "change of\nVALUE (UNITS)",
        "ave_pr_val_ptre": "change of\nVALUE (UNITS)",
        "ave_pr_val_ptro": "change of\nVALUE (UNITS)",
        "ave_pr_val_ptrw": "change of\nVALUE (UNITS)",
        "ave_ts_val_glob": "change of\nVALUE (UNITS)",
        "ave_ts_val_aequ": "change of\nVALUE (UNITS)",
        "ave_ts_val_atno": "change of\nVALUE (UNITS)",
        "ave_ts_val_atro": "change of\nVALUE (UNITS)",
        "ave_ts_val_iequ": "change of\nVALUE (UNITS)",
        "ave_ts_val_isee": "change of\nVALUE (UNITS)",
        "ave_ts_val_itro": "change of\nVALUE (UNITS)",
        "ave_ts_val_iweq": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxip": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpb": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpe": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpl": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpn": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpr": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxps": "change of\nVALUE (UNITS)",
        "ave_ts_val_gxpw": "change of\nVALUE (UNITS)",
        "ave_ts_val_gypb": "change of\nVALUE (UNITS)",
        "ave_ts_val_gype": "change of\nVALUE (UNITS)",
        "ave_ts_val_gypi": "change of\nVALUE (UNITS)",
        "ave_ts_val_n30e": "change of\nVALUE (UNITS)",
        "ave_ts_val_n30n": "change of\nVALUE (UNITS)",
        "ave_ts_val_n30s": "change of\nVALUE (UNITS)",
        "ave_ts_val_n34e": "change of\nVALUE (UNITS)",
        "ave_ts_val_n34n": "change of\nVALUE (UNITS)",
        "ave_ts_val_n34s": "change of\nVALUE (UNITS)",
        "ave_ts_val_n40e": "change of\nVALUE (UNITS)",
        "ave_ts_val_n40n": "change of\nVALUE (UNITS)",
        "ave_ts_val_n40s": "change of\nVALUE (UNITS)",
        "ave_ts_val_peno": "change of\nVALUE (UNITS)",
        "ave_ts_val_peqe": "change of\nVALUE (UNITS)",
        "ave_ts_val_pequ": "change of\nVALUE (UNITS)",
        "ave_ts_val_peqw": "change of\nVALUE (UNITS)",
        "ave_ts_val_peso": "change of\nVALUE (UNITS)",
        "ave_ts_val_pnee": "change of\nVALUE (UNITS)",
        "ave_ts_val_pnew": "change of\nVALUE (UNITS)",
        "ave_ts_val_pnqe": "change of\nVALUE (UNITS)",
        "ave_ts_val_ptre": "change of\nVALUE (UNITS)",
        "ave_ts_val_ptro": "change of\nVALUE (UNITS)",
        "ave_ts_val_ptrw": "change of\nVALUE (UNITS)",
        "ave_ts_val_trop": "change of\nVALUE (UNITS)",
        "std_pr_ano_n34e": "relative change of\nVALUE (%)",
        "std_ts_ano_n34e": "relative change of\nVALUE (%)",
        "full--cmip5": "THR% full",
        "full--cmip6": "THR% full",
        "mean--cmip5": "MMM",
        "mean--cmip6": "MMM",
        "figure_title": "Mean state influence on ENSO change",
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
def e29_change_explained_plot(
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
    # -- compute change compared to initial conditions
    dur2, exp2 = data_epoch_lengths[-1], "historical"
    # ensemble mean
    change_em = {}
    for dia in list(values_json.keys()):
        for dur1 in list(values_json[dia].keys()):
            for pro in list(values_json[dia][dur1].keys()):
                for exp1 in list(values_json[dia][dur1][pro].keys()):
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
                        # select reference epoch
                        # yy = 2014 if data_projects[0] == "cmip6" else 2004
                        yy = 2015 if data_projects[0] == "cmip6" else 2005
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 -= da2
                        if dia.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # centered time
                        da1 = da1.assign_coords({"year": da1["year"] - int(dur1.split("_")[0]) / 2})
                        # save
                        put_in_dict(change_em, da1, dur1, dia, exp1, pro, dat)
    # mean
    mean_change, _, _ = compute_1var(change_em, "statistic",
                                     function_arguments={"dim": "member", "stat": "mea"})
    # -- select top / bottom models
    # select models based on variability
    print("select models based on variability")
    selected = {}
    for dur1 in list(mean_change.keys()):
        if dur1 != dur2:
            continue
        for dia in list(mean_change[dur1].keys()):
            if dia.split("_")[0] not in ["std", "var"]:
                continue
            for exp in list(mean_change[dur1][dia].keys()):
                for pro in list(mean_change[dur1][dia][exp].keys()):
                    # regroup models
                    dict_t = {}
                    for dat in list(mean_change[dur1][dia][exp][pro].keys()):
                        # array
                        da = copy__deepcopy(mean_change[dur1][dia][exp][pro][dat])
                        # select time
                        yy = da["year"].values
                        yy = yy[0] if exp == "historical" else yy[-1]
                        da = da.sel(indexers={"year": yy})
                        # save
                        dict_t[dat] = float(da.values)
                    # select top / bottom models
                    dict_t = {k: v for k, v in sorted(dict_t.items(), key=lambda item: item[1])}
                    list_datasets = list(dict_t.keys())
                    nn = int(round(len(list_datasets) * proc_th / 100, 0))
                    # save
                    put_in_dict(selected, list_datasets[:nn], exp, pro, "bottom")
                    put_in_dict(selected, list_datasets[-nn:], exp, pro, "top")
    # select models for all variables
    percentiles = [50 + k * (proc_ci / 2) for k in [-1, 1]]
    change_m, change_r = {}, {}
    print("select models for all variables")
    for dur1 in list(mean_change.keys()):
        if dur1 == dur2:
            continue
        for dia in list(mean_change[dur1].keys()):
            for exp in list(mean_change[dur1][dia].keys()):
                for pro in list(mean_change[dur1][dia][exp].keys()):
                    if not check_dict_levels(selected, exp, pro):
                        continue
                    for grp in list(selected[exp][pro].keys()):
                        # concatenate
                        da = ()
                        for dat in list(mean_change[dur1][dia][exp][pro].keys()):
                            if dat in selected[exp][pro][grp]:
                                da += (mean_change[dur1][dia][exp][pro][dat],)
                        da = concatenate(da, "dataset")
                        # mean
                        da_m = da.mean(dim="dataset")
                        # range
                        da_r = statistic(da, dim="dataset", stat="per", per=percentiles, skipna=True)
                        # save
                        put_in_dict(change_m, da_m, dur1, dia, exp, pro, grp)
                        put_in_dict(change_r, da_r, dur1, dia, exp, pro, grp)
    # min-max
    print("change min-max")
    for dur in list(change_r.keys()):
        for dia in list(change_r[dur].keys()):
            tmp = []
            for exp in list(change_r[dur][dia].keys()):
                for pro in list(change_r[dur][dia][exp].keys()):
                    for da in change_r[dur][dia][exp][pro].values():
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
    for dur in list(change_m.keys()):
        if not check_dict_levels(change_r, dur):
            continue
        ny1 = int(dur.split("_")[0])
        ny2 = int(data_epoch_lengths[-1].split("_")[0])
        list_dia = list(set(list(change_m[dur].keys())) & set(list(change_r[dur].keys())))
        list_dia = [k for k in data_diagnostics1 if k in list_dia]
        list_exp1 = list(set([k2 for k1 in list_dia for k2 in list(change_m[dur][k1].keys())]))
        list_exp2 = list(set([k2 for k1 in list_dia for k2 in list(change_r[dur][k1].keys())]))
        list_exp = [k for k in data_experiments if k in list_exp1 + list_exp2]
        for ii, dia in enumerate(list_dia):
            for jj, exp in enumerate(list_exp):
                tp = "panel_1" if exp == "historical" else "panel_2"
                pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                tuple_k = (dur, grp, pan)
                print(dur, dia, exp.ljust(10), tuple_k)
                # panel size
                tk = tuple_k[:1]
                if not check_dict_levels(figure_panel_size, *tk + (pan,)):
                    for k1, d1 in fig_panel_size.items():
                        if k1 == tp and not check_dict_levels(figure_panel_size, *tk + (pan,)):
                            for k2, d2 in d1.items():
                                if check_dict_levels(figure_panel_size, *tk + (pan, k2,)):
                                    continue
                                tmp = copy__deepcopy(d2)
                                if k2 == "x_size" and check_dict_levels(fig_ticks, dur, "time", tp, "lim"):
                                    k3 = read_from_dict(fig_ticks, dur, "time", tp, "lim")
                                    tmp *= (k3[-1] - k3[0])
                                put_in_dict(figure_panel_size, tmp, *tk + (pan, k2,))
                        elif k1[:6] != "panel_" and not check_dict_levels(figure_panel_size, *tk + (k1,)):
                            put_in_dict(figure_panel_size, copy__deepcopy(d1), *tk + (k1,))
                # panel axes definition
                l1, l2 = ["x_", "y_"], ["time", dia]
                for vv, dd in zip(l1, l2):
                    if check_dict_levels(fig_ticks, dur, dd, tp):
                        for k1, k2 in fig_ticks[dur][dd][tp].items():
                            k3 = str(vv) + str(k1)
                            if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                continue
                            tmp = copy__deepcopy(k2)
                            if isinstance(tmp, str) and k1 == "nam" and dd not in ["time"]:
                                # yy = 2014 if data_projects[0] == "cmip6" else 2004
                                yy = 2015 if data_projects[0] == "cmip6" else 2005
                                tmp = tmp.replace("YEAR1", str(yy - ny2)).replace("YEAR2", str(yy - 1))
                                tmp = tmp.replace("VALUE", metadata_json[dd]["name_short"])
                                tmp = tmp.replace("UNITS", metadata_json[dd]["units"])
                                tmp = tmp.replace(" ()", "")
                            if k3 == "x_nam" and dia != list_dia[-1]:
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
                    val = str(fig_titles[exp]) + " (" + str(len(n_ensemble[exp])) + ")"
                    put_in_dict(plot_data, val, *tuple_k + ("title_col",))
                if jj == 0:
                    val = copy__deepcopy(fig_titles[dia]).replace("VALUE", metadata_json[dia]["name_short"])
                    val = val.replace("UNITS", metadata_json[dia]["units"])
                    put_in_dict(plot_data, val, *tuple_k + ("title_row",))
                # -- curve
                lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                if check_dict_levels(change_m, dur, dia, exp):
                    for pro in list(change_m[dur][dia][exp].keys()):
                        for sel, da in change_m[dur][dia][exp][pro].items():
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
                    if jj == 0 and check_dict_levels(change_m, dur, dia, exp):
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
                if check_dict_levels(change_r, dur, dia, exp):
                    for pro in list(change_r[dur][dia][exp].keys()):
                        for sel, da in change_r[dur][dia][exp][pro].items():
                            la += [0.6]
                            lc += [fig_colors[sel] if sel in list(fig_colors.keys()) else "cyan"]
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
                    if jj == 0 and check_dict_levels(change_r, dur, dia, exp):
                        la += [0.6]
                        lc += ["k"]
                        lx += [[x1 + 2 * dx, x1 + 10 * dx]]
                        ly1 += [[y2 - (1 * 8 + 5 + 3) * dy] * 2]
                        ly2 += [[y2 - (1 * 8 + 5 - 3) * dy] * 2]
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
                if isinstance(do_legend, bool) and do_legend and jj == 0:
                    df = read_from_dict(figure_panel_size, *tk)
                    x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                    dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                    dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                    # legend
                    list_m = ["MMM", "-".join([str(int(k3)) if k3 == int(k3) else str(k3) for k3 in percentiles]) + "%"]
                    for k1, v1 in enumerate(list_m):
                        lc += ["k"]
                        lf += [15]
                        lh += ["left"]
                        lt += [v1]
                        lv += ["center"]
                        lw += ["regular"]
                        lx += [x1 + (20 if k1 == 0 else 12) * dx]
                        ly += [y2 - (k1 * 8 + 5) * dy]
                    list_c = ["bottom", "top"]
                    for k1, v1 in enumerate(list_c):
                        lc += [fig_colors[v1] if v1 in list(fig_colors.keys()) else "cyan"]
                        lf += [15]
                        lh += ["right"]
                        lt += [str(v1) + " " + str(proc_th) + "%"]
                        lv += ["center"]
                        lw += ["regular"]
                        lx += [x2 - 2 * dx]
                        ly += [y1 + (k1 * 8 + 5) * dy]
                    # plot title
                    if ii == 0 and jj == 0:
                        # plot size
                        l1 = len([k1 for k1 in list_exp if k1[:3] not in ["rcp", "ssp"]])
                        x1_d = fig_panel_size[tp]["x_delt"]
                        x1_m = fig_panel_size[tp]["x_size"]
                        x1_s = read_from_dict(fig_ticks, dur, "time", tp, "lim")
                        x1_s = max(x1_s) - min(x1_s)
                        l2 = len([k1 for k1 in list_exp if k1[:3] in ["rcp", "ssp"]])
                        x2_d = fig_panel_size["panel_2"]["x_delt"]
                        x2_m = fig_panel_size["panel_2"]["x_size"]
                        x2_s = read_from_dict(fig_ticks, dur, "time", "panel_2", "lim")
                        x2_s = max(x2_s) - min(x2_s)
                        # add to text
                        lc += ["k"]
                        lf += [20]
                        lh += ["center"]
                        lt += [fig_titles["figure_title"]]
                        lv += ["center"]
                        lw += ["bold"]
                        a1 = x1_m * x1_s
                        a2 = (x1_m * x1_d) * l1 + x2_m * x2_s * l2 + x2_d * (l2 - 1)
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
    a1, a2 = str(smile_minimum_size).zfill(2), str(proc_th).zfill(2)
    a3 = str(int(proc_ci)) if proc_ci == int(proc_ci) else str(proc_ci)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_".join(
        ["", str(a1) + "mem", str(a2) + "grp", str(a3) + "p", str(data_projects[0]),
         str(data_diagnostics1[1].split("_")[1]) + "1"])
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
e29_change_explained_plot()
