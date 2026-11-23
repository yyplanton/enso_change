# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Figure: scatterplots of change vs. changes (departure from the end of the 20th century) in CMIP5 & CMIP6, raw & pcc
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy as copy__deepcopy
from math import ceil as math__ceil
from math import floor as math__floor
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
from enso_change.wrapper.statistics import linear_regression, numpy__array
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
    "data_diagnostics1": ["ave_ts_val_glob", "ave_pr_val_pequ", "ave_pr_val_n30e", "ave_pr_val_n34e",
                          "ave_ts_val_n30e", "ave_ts_val_n34e", "ave_ts_val_n40e",
                          "fbk_ts_n30e_to_hf_n30e", "fbk_ts_n34e_to_hf_n34e", "fbk_ts_n34e_to_tx_n40e",
                          "ske_pr_ano_n30e", "ske_ts_ano_n30e", "std_pr_ano_n34e", "std_ts_ano_n34e"],
    # "data_diagnostics1": ["ave_ts_val_glob", "ave_pr_val_n30e", "std_ts_ano_n34e"],
    # "data_diagnostics1": ["ave_ts_val_glob", "ave_pr_val_pequ", "ave_pr_val_n30e", "ave_pr_val_n34e",
    #                       "ave_ts_val_gxip", "ave_ts_val_gxpb", "ave_ts_val_gxpe", "ave_ts_val_gxpl",
    #                       "ave_ts_val_gxpn", "ave_ts_val_gxpr", "ave_ts_val_gxps", "ave_ts_val_gxpw",
    #                       "ave_ts_val_gypb", "ave_ts_val_gype", "ave_ts_val_gypi", "ave_ts_val_n30e",
    #                       "ave_ts_val_n34e", "fbk_ts_n34e_to_tx_n40e", "ske_pr_ano_n30e", "ske_ts_ano_n30e",
    #                       "std_pr_ano_n34e", "std_ts_ano_n34e"],
    # "data_diagnostics1": ["ave_ts_val_glob", "ave_pr_val_n30e"],
    "data_diagnostics2": [],
    # list of epoch lengths: list[str]
    "data_epoch_lengths": ["050_year_epoch"],
    # list of experiments: list[str]
    # "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585", "ssp370"],
    "data_experiments": ["historical", "rcp26", "rcp45", "rcp85", "ssp126", "ssp245", "ssp585"],
    # dictionary of observational datasets desired for each diagnostic: dict[str, list[str]]
    "data_observational_datasets":  {
        "ave_ts_val_glob": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_pequ": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_pr_val_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ave_ts_val_gxip": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpb": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpe": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpl": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpn": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpr": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxps": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gxpw": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gypb": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gype": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_gypi": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_n30e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_n34e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ave_ts_val_n40e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "fbk_ts_n30e_to_hf_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "fbk_ts_n34e_to_hf_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "fbk_ts_n34e_to_tx_n40e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_n30e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_pr_ano_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "ske_ts_ano_n30e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "ske_ts_ano_n34e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
        "std_pr_ano_n34e": ["20CRv2c", "20CRv3", "ERA5", "NCEP1"],
        "std_ts_ano_n34e": ["20CRv2c", "20CRv3", "COBE2", "ERA5", "ERSSTv5", "HadISST", "NCEP1"],
    },
    # list of projects: list[str]
    # "data_projects": ["cmip5", "cmip6"],  # default_arguments["data_projects"],
    "data_projects": ["cmip5", "cmip6"],
    # minimum SMILE size: int[1, 100]
    "smile_minimum_size": 1,
    # list of rejected SMILEs: list[str]
    # "smile_rejected": default_arguments["smile_rejected"],
    # "smile_rejected": ["CMCC-CM", "EC-EARTH", "KACE-1-0-G"],
    # "smile_rejected": ["EC-EARTH", "KACE-1-0-G"],
    "smile_rejected": ["EC-EARTH"],
    # "smile_rejected": [],
    # reject SMILE if piControl is not available (or too small): bool
    "smile_require_picontrol": False,
    #
    # -- Data processing
    #
    # confidence interval: float
    "proc_ci": 95.,
    # IV correction: dict[str: dict[str, list[float]]]
    "proc_co": {},
    # method: {"piControl", "experiment"}
    "proc_md": "experiment",
    # number of pseudo samples: int
    "proc_ps": 10000,
    # threshold for classification
    "proc_th": 25,
    #
    # -- Figure
    #
    # figure format: eps, pdf, png, svg
    "fig_format": "pdf",
    # selected datasets
    "fig_datasets": {},
    # size of each panel
    "fig_panel_size": {
        "frac": {"x": 0.2, "y": 0.2},
        "panel_01": {"x_delt": 7, "x_size": 20, "y_delt": 3, "y_size": 20},
    },
    # colors
    "fig_colors": {
        "cmip5": "green",
        "cmip6": "r",
        "observations": "k",
        **default_arguments["fig_colors_dat"],
        **default_arguments["fig_colors_exp"],
    },
    # line styles
    "fig_linestyles": {
        "cmip5": "-",
        "cmip6": "-",
        **default_arguments["fig_linestyles_dat"],
        **default_arguments["fig_linestyles_exp"],
    },
    # markers
    "fig_markers": {
        "cmip5": "s",
        "cmip6": "D",
        "observations": "o",
        **default_arguments["fig_markers_dat"],
    },
    # ticks
    "fig_ticks": {
        "ave_pr_val_pequ--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_pr_val_pequ--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_pr_val_n30e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_pr_val_n30e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_pr_val_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_pr_val_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_glob--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_glob--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxip--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxip--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpb--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpb--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpe--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpe--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpl--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpl--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpn--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpn--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpr--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpr--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxps--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxps--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpw--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gxpw--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gypb--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gypb--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gype--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gype--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_gypi--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_gypi--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_n30e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_n30e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ave_ts_val_n40e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ave_ts_val_n40e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "fbk_ts_n30e_to_hf_n30e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "fbk_ts_n30e_to_hf_n30e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "fbk_ts_n34e_to_hf_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "fbk_ts_n34e_to_hf_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "fbk_ts_n34e_to_tx_n40e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(UNITS.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (UNITS.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "fbk_ts_n34e_to_tx_n40e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (UNITS)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(UNITS)",
                "tic": None,
            },
        },
        "ske_pr_ano_n30e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ske_pr_ano_n30e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
        "ske_pr_ano_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ske_pr_ano_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
        "ske_ts_ano_n30e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ske_ts_ano_n30e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
        "ske_ts_ano_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "ske_ts_ano_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
        "std_pr_ano_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "std_pr_ano_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
        "std_ts_ano_n34e--pcc": {
            "panel_01": {
                # "lim": [-50, 50],
                # "nam": "VALUE per CC\n(%.$^\circ$C$^{-1}$)",
                # "tic": list(range(-50, 51, 25)),
                "lim": None,
                "nam": "change VALUE\nper CC (%.$^\circ$C$^{-1}$)",
                "tic": None,
            },
        },
        "std_ts_ano_n34e--raw": {
            "panel_01": {
                # "lim": [-70, 70],
                # "nam": "VALUE (%)",
                # "tic": list(range(-70, 71, 35)),
                "lim": None,
                "nam": "change VALUE\n(%)",
                "tic": None,
            },
        },
    },
    # titles
    "fig_titles": {
        "bia": "bias",
        "pcc": "change per GMST",
        "raw": "change",
        "cmip5": "CMIP5",
        "cmip6": "CMIP6",
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
def e60_change_vs_change_scatter_plot(
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
    print("reader_json", str(len(list(values_json.keys()))).rjust(3), list(values_json.keys()))
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
                for exp in list(values_json[dia][dur][pro].keys()):
                    put_in_dict(n_ensemble, list(values_json[dia][dur][pro][exp].keys()), exp)
    for exp in data_experiments:
        list_experiments = copy__deepcopy(n_ensemble[exp])
        del n_ensemble[exp]
        list_experiments = sorted(list(set(list_experiments)), key=lambda s: s.lower())
        put_in_dict(n_ensemble, copy__deepcopy(list_experiments), exp)
        print(exp.ljust(12), "number of ensembles: " + str(len(list_experiments)).rjust(2))
        if len(list_experiments) < 10:
            print("".ljust(20), ", ".join(list_experiments))
    #
    # -- Processing
    #
    dia2, dur2, exp2 = data_diagnostics1[0], data_epoch_lengths[-1], "historical"
    # -- change
    change_em = {}
    for dia1 in list(values_json.keys()):
        for dur1 in list(values_json[dia1].keys()):
            if dur1 != dur2:
                continue
            for pro in list(values_json[dia1][dur1].keys()):
                for exp1 in list(values_json[dia1][dur1][pro].keys()):
                    if exp1 == exp2:
                        continue
                    for dat in list(values_json[dia1][dur1][pro][exp1].keys()):
                        if not check_dict_levels(values_json, dia1, dur2, pro, exp1, dat) or \
                                not check_dict_levels(values_json, dia1, dur2, pro, exp2, dat) or \
                                not check_dict_levels(values_json, dia2, dur1, pro, exp1, dat) or \
                                not check_dict_levels(values_json, dia2, dur2, pro, exp1, dat) or \
                                not check_dict_levels(values_json, dia2, dur2, pro, exp2, dat):
                            continue
                        # array
                        da1 = copy__deepcopy(values_json[dia1][dur1][pro][exp1][dat])
                        if exp1[:3] in ["rcp", "ssp"]:
                            da2 = copy__deepcopy(values_json[dia1][dur2][pro][exp1][dat])
                        else:
                            da2 = copy__deepcopy(values_json[dia1][dur2][pro][exp2][dat])
                        da3 = copy__deepcopy(values_json[dia2][dur1][pro][exp1][dat])
                        if exp1[:3] in ["rcp", "ssp"]:
                            da4 = copy__deepcopy(values_json[dia2][dur2][pro][exp1][dat])
                        else:
                            da4 = copy__deepcopy(values_json[dia2][dur2][pro][exp2][dat])
                        # select experiment epoch
                        yy = 1850 + int(dur1.split("_")[0]) if exp1 == "historical" else 2100
                        da1 = da1.sel(indexers={"year": yy}).mean(dim="member")
                        da3 = da3.sel(indexers={"year": yy}).mean(dim="member")
                        # select reference epoch
                        yy = 2006
                        da2 = da2.sel(indexers={"year": yy}).mean(dim="member")
                        da4 = da4.sel(indexers={"year": yy}).mean(dim="member")
                        # change
                        da1 = da1 - da2 if exp1[:3] in ["rcp", "ssp"] else da2 - da1
                        if dia1.split("_")[0] in ["std", "var"]:
                            da1 *= 100 / da2
                        # degree of GMST
                        da_c = da3 - da4 if exp1[:3] in ["rcp", "ssp"] else da4 - da3
                        exp_o = copy__deepcopy(exp1) if exp1 in ["piControl", "historical", "ssp370"] else \
                            "scenario" + str(exp1[-2:])
                        for war in ["raw", "pcc"]:
                            da_o = copy__deepcopy(da1)
                            # raw change or change per degree of GMST
                            if war == "pcc" and exp1[:3] in ["rcp", "ssp"]:
                                da_o /= da_c
                            # save
                            put_in_dict(change_em, float(da_o.values), dur1, dia1, war, exp_o, pro, dat)
    print("change_em", list(change_em.keys()))
    # min-max
    print("change_em min-max")
    change_minmax = {}
    for dur in list(change_em.keys()):
        for dia in list(change_em[dur].keys()):
            for war in list(change_em[dur][dia].keys()):
                tt = [vv for k1 in list(change_em[dur][dia][war].keys())
                      for k2 in list(change_em[dur][dia][war][k1].keys())
                      for vv in change_em[dur][dia][war][k1][k2].values()]
                if len(tt) > 0:
                    print(dur, dia.ljust(22), war, "{0:.2f}".format(round(min(tt), 2)).rjust(8),
                          "{0:.2f}".format(round(max(tt), 2)).rjust(8))
                    t1, t2 = min(tt), max(tt)
                    tt = [math__floor(t1), math__ceil(t2)]
                    if t2 - t1 < 5:
                        n1 = 3 if t2 - t1 < 0.05 else (2 if t2 - t1 < 0.5 else 1)
                        n2 = 10 ** n1
                        tt = [round(math__floor(t1 * n2) / n2, n1), round(math__ceil(t2 * n2) / n2, n1)]
                    if min(tt) == max(tt):
                        tt = [math__floor(min(tt)) - 0.5, math__ceil(min(tt) + 0.5)]
                    put_in_dict(change_minmax, copy__deepcopy(tt), dur, dia, war)
    #
    # -- Organize data for plot
    #
    figure_axes, figure_ps, plot_data = {}, {}, {}
    # panel 1: scatter
    print("scatter")
    group, panel = "group_1", "panel_01"
    for dur in list(change_em.keys()):
        for dia_x in list(change_em[dur].keys()):
            list_y1, list_y2 = ["raw", "pcc"], [change_em[dur], change_em[dur]]
            list_y3 = [change_minmax[dur], change_minmax[dur]]
            for ii, (typ, d1, d2) in enumerate(zip(list_y1, list_y2, list_y3)):
                list_dia_y = list(d1.keys())
                for jj, dia_y in enumerate(list_dia_y):
                    pan = str(panel.split("_")[0]) + "_" + str(int(panel.split("_")[-1]) + jj).zfill(2)
                    grp = str(group.split("_")[0]) + "_" + str(int(group.split("_")[-1]) + ii).zfill(2)
                    tuple_k = (dur, dia_x, grp, pan)
                    print(dur, dia_x.ljust(22), dia_y.ljust(22), tuple_k)
                    # dictionaries
                    dix, diy = change_em[dur][dia_x][typ], d1[dia_y][typ]
                    # min-max
                    x_mm, y_mm = change_minmax[dur][dia_x][typ], d2[dia_y][typ]
                    # panel size
                    tk = tuple_k[:-2]
                    if not check_dict_levels(figure_ps, *tk + (pan,)):
                        put_in_dict(figure_ps, copy__deepcopy(fig_panel_size[panel]), *tk + (pan,))
                        for k1, k2 in fig_panel_size.items():
                            if k1[:6] != "panel_" and not check_dict_levels(figure_ps, *tk + (k1,)):
                                put_in_dict(figure_ps, k2, *tk + (k1,))
                    # panel axes definition
                    l1, l2 = ["x_", "y_"], [str(dia_x) + "--" + str(typ), str(dia_y) + "--" + str(typ)]
                    for vv, dd1 in zip(l1, l2):
                        dd2 = dd1.split("--")[0]
                        if check_dict_levels(fig_ticks, dd1, panel):
                            for k1, k2 in fig_ticks[dd1][panel].items():
                                k3 = str(vv) + str(k1)
                                if check_dict_levels(figure_axes, *tuple_k + (k3,)):
                                    continue
                                tmp = copy__deepcopy(k2)
                                if isinstance(tmp, str) and k1 == "nam":
                                    tmp = tmp.replace("VALUE", metadata_json[dd2]["name_short"])
                                    tmp = tmp.replace("UNITS", metadata_json[dd2]["units"])
                                    tmp = tmp.replace(" ()", "").replace("()", "")
                                elif tmp is None and k1 in ["lim", "tic"]:
                                    tmp = copy__deepcopy(x_mm) if vv == "x_" else copy__deepcopy(y_mm)
                                if k3 == "x_nam" and typ != list_y1[-1]:
                                    tmp = ""
                                put_in_dict(figure_axes, tmp, *tuple_k + (k3,))
                    dt = {}
                    if check_dict_levels(figure_axes, *tuple_k):
                        dt = read_from_dict(figure_axes, *tuple_k)
                    # test possibility to plot text / legend
                    do_legend = False
                    if "x_lim" in list(dt.keys()) and "y_lim" in list(dt.keys()) and dt["x_lim"] is not None and \
                            dt["y_lim"] is not None:
                        do_legend = True
                    # -- title
                    if ii == 0:
                        put_in_dict(plot_data, metadata_json[dia_y]["name_short"], *tuple_k + ("title_col",))
                    if jj == 0:
                        put_in_dict(plot_data, fig_titles[typ], *tuple_k + ("title_row",))
                    # -- marker
                    lc, lm, ls, lw, lx, ly, lz = [], [], [], [], [], [], []
                    for exp in list(dix.keys()):
                        for pro in list(dix[exp].keys()):
                            for dat in list(dix[exp][pro].keys()):
                                if not check_dict_levels(diy, exp, pro, dat):
                                    continue
                                da_x, da_y = copy__deepcopy(dix[exp][pro][dat]), copy__deepcopy(diy[exp][pro][dat])
                                val = "cyan"
                                if pro == "observations" and pro in list(fig_colors.keys()):
                                    val = fig_colors[pro]
                                elif pro != "observations" and exp in list(fig_colors.keys()):
                                    val = fig_colors[exp]
                                lc += [val]
                                lm += [fig_markers[pro] if pro in list(fig_markers.keys()) else "o"]
                                ls += [40]
                                lw += [0.5]
                                lx += [da_x]
                                ly += [da_y]
                                lz += [numpy__random__randint(3, 9) ]
                    pt = "mar"
                    if len(lx) > 0 and len(lx) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(pt) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(pt) + "_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lm), *tuple_k + (str(pt) + "_m",))
                        put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(pt) + "_s",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(pt) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(pt) + "_z",))
                    # -- curve
                    lc, ls, lw, lx, ly, lz = [], [], [], [], [], []
                    for exp in list(dix.keys()):
                        da_x, da_y = [], []
                        for pro in list(dix[exp].keys()):
                            if pro[:4] != "cmip":
                                continue
                            for dat in list(dix[exp][pro].keys()):
                                if not check_dict_levels(diy, exp, pro, dat):
                                    continue
                                da_x += [dix[exp][pro][dat]]
                                da_y += [diy[exp][pro][dat]]
                        if len(list(set(da_x))) < 2 or len(list(set(da_y))) < 2:
                            continue
                        ln = linear_regression(numpy__array(da_x), numpy__array(da_y))
                        lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                        ls += ["-"]
                        lw += [1]
                        lx += [x_mm]
                        ly += [[ln["s"] * k + ln["i"] for k in x_mm]]
                        lz += [20]
                    if isinstance(do_legend, bool) and do_legend:
                        # x = 0 and y = 0
                        for n1, k1 in zip(["x", "y"], [dt["x_lim"], dt["y_lim"]]):
                            if k1[0] < 0 < k1[1]:
                                lc += ["k"]
                                ls += [":"]
                                lw += [1]
                                lx += [[0] * 2 if n1 == "x" else dt["x_lim"]]
                                ly += [dt["y_lim"] if n1 == "x" else [0] * 2]
                                lz += [20]
                    pt = "cur"
                    if len(lx) > 0 and len(lx) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(pt) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(ls), *tuple_k + (str(pt) + "_ls",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(pt) + "_lw",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(pt) + "_y",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(pt) + "_z",))
                    # -- shading
                    la, lc, lx, ly1, ly2, lz = [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        for exp in list(dix.keys()):
                            if exp != "historical":
                                continue
                            da_x, da_y = [], []
                            for pro in list(dix[exp].keys()):
                                if pro != "observations":
                                    continue
                                for dat in list(dix[exp][pro].keys()):
                                    da_x += [dix[exp][pro][dat]]
                                    if check_dict_levels(diy, exp, pro, dat):
                                        da_y += [diy[exp][pro][dat]]
                            if len(da_x) == 0:
                                da_x = copy__deepcopy(dt["x_lim"])
                            if len(da_y) == 0:
                                da_y = copy__deepcopy(dt["y_lim"])
                            pro2 = "observations"
                            la += [0.5]
                            lc += [fig_colors[pro2] if pro2 in list(fig_colors.keys()) else "cyan"]
                            lx += [[min(da_x), max(da_x)]]
                            ly1 += [[min(da_y)] * 2]
                            ly2 += [[max(da_y)] * 2]
                            lz += [2]
                    pt = "sha"
                    if len(lx) > 0 and len(lx) == len(ly1) and len(lx) == len(ly2):
                        put_in_dict(plot_data, copy__deepcopy(la), *tuple_k + (str(pt) + "_a",))
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(pt) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly1), *tuple_k + (str(pt) + "_y1",))
                        put_in_dict(plot_data, copy__deepcopy(ly2), *tuple_k + (str(pt) + "_y2",))
                        put_in_dict(plot_data, copy__deepcopy(lz), *tuple_k + (str(pt) + "_z",))
                    # -- text
                    lc, lf, lh, lr, lt, lv, lw, lx, ly = [], [], [], [], [], [], [], [], []
                    if isinstance(do_legend, bool) and do_legend:
                        df = read_from_dict(figure_ps, *tk)
                        x1, x2, y1, y2 = dt["x_lim"] + dt["y_lim"]
                        dx = (x2 - x1) * default_plot["size_x"] / (100 * df[pan]["x_size"] * df["frac"]["x"])
                        dy = (y2 - y1) * default_plot["size_y"] / (100 * df[pan]["y_size"] * df["frac"]["y"])
                        # regression
                        dt = {}
                        for exp in list(dix.keys()):
                            da_x, da_y = [], []
                            for pro in list(dix[exp].keys()):
                                if pro[:4] != "cmip":
                                    continue
                                for dat in list(dix[exp][pro].keys()):
                                    if not check_dict_levels(diy, exp, pro, dat):
                                        continue
                                    da_x += [dix[exp][pro][dat]]
                                    da_y += [diy[exp][pro][dat]]
                            if len(list(set(da_x))) < 2 or len(list(set(da_y))) < 2:
                                continue
                            ln = linear_regression(numpy__array(da_x), numpy__array(da_y))
                            dt[exp] = ln["r"]
                        tt = 0 if len(list(dt.keys())) == 0 else float(numpy__array(list(dt.values())).mean())
                        l3 = sorted(list(dt.keys()), key=lambda s: s.lower())
                        if tt >= 0:
                            l3 = list(reversed(l3))
                        for kk, exp in enumerate(l3):
                            lc += [fig_colors[exp] if exp in list(fig_colors.keys()) else "cyan"]
                            lf += [15]
                            lh += ["right"]
                            lr += [0]
                            val = (fig_titles[exp][:2] if fig_titles[exp].lower() == "historical" else
                                   fig_titles[exp]) + " " + "{0:+.2f}".format(round(dt[exp], 2))
                            lt += [val]
                            lv += ["bottom" if tt >= 0 else "top"]
                            lw += ["regular"]
                            lx += [x2 - 2 * dx]
                            ly += [y1 + (7.5 * kk + 1) * dy if tt >= 0 else y2 - (7.5 * kk + 1) * dy]
                    pt = "text"
                    if len(lt) > 0 and len(lt) == len(lx) and len(lt) == len(ly):
                        put_in_dict(plot_data, copy__deepcopy(lt), *tuple_k + (pt,))
                        put_in_dict(plot_data, copy__deepcopy(lc), *tuple_k + (str(pt) + "_c",))
                        put_in_dict(plot_data, copy__deepcopy(lf), *tuple_k + (str(pt) + "_fs",))
                        put_in_dict(plot_data, copy__deepcopy(lh), *tuple_k + (str(pt) + "_ha",))
                        put_in_dict(plot_data, copy__deepcopy(lr), *tuple_k + (str(pt) + "_r",))
                        put_in_dict(plot_data, copy__deepcopy(lv), *tuple_k + (str(pt) + "_va",))
                        put_in_dict(plot_data, copy__deepcopy(lw), *tuple_k + (str(pt) + "_w",))
                        put_in_dict(plot_data, copy__deepcopy(lx), *tuple_k + (str(pt) + "_x",))
                        put_in_dict(plot_data, copy__deepcopy(ly), *tuple_k + (str(pt) + "_y",))
    #
    # -- Plot
    #
    # output plot directory (relative to current file directory)
    figure_directory = "/".join(os__path__dirname(__file__).split("/")[:-2]) + "/figures"
    # figure name
    a1, a2 = str(len(data_diagnostics1)).zfill(2), str(len(smile_rejected)).zfill(2)
    figure_name = __file__.split("/")[-1].split(".")[0] + "_" + str(a1) + "dia_" + str(a2) + "rem"
    # figure loop
    for dur in sorted(list(plot_data.keys()), key=lambda s: s.lower()):
        for dia in sorted(list(plot_data[dur].keys()), key=lambda s: s.lower()):
            # list plot names
            list_names = list(plot_data[dur][dia].keys())
            # path to output figure
            figure_name_and_path = os__path__join(figure_directory, figure_name + "_".join(("", dur, dia)))
            # plot
            fig_basic(plot_data[dur][dia], list_names, 1, figure_axes[dur][dia], fig_format,
                      figure_name_and_path, figure_ps[dur][dia], fig_title_bool=True, panel_param=panel_param)
# ---------------------------------------------------------------------------------------------------------------------#


# call function
e60_change_vs_change_scatter_plot()
